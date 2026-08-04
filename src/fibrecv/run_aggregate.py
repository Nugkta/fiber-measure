"""CLI: group per-image profiles by A_B and build registered averages.

Dependencies
------------
``numpy``, ``pandas``, ``matplotlib`` (Agg), the ``register`` + ``io_utils`` +
``config`` modules. Reads the CSV/meta artifacts produced by ``run_measure``.

Inputs
------
Command-line flags: ``--out`` (output root holding ``per_image/``), selector
(``--groups`` / ``--all``), and a few CONFIG overrides (``--ppu``, ``--max-shift``,
``--min-corr``, ``--min-coverage``, ``--rep-dev-frac``, ``--anomaly-exclude``).

Output
------
Per A_B sample: ``per_sample/csv/sample_<A_B>_registered.csv`` (mean+/-std curve),
``per_sample/plots/sample_<A_B>_registered.png`` (mean + +/-std band),
``per_sample/shifts/sample_<A_B>_shifts.json``, one row per sample in
``summary/master_summary.csv``, and one row per image (excluded ones included,
with anomaly flags + replicate_outlier deviations and the exclusion reason) in
``summary/per_image_summary.csv``.

Pos
---
Second entrypoint of the two-stage workflow. Cheap; runs after ``run_measure``.
Implements the locked decision that the 3 replicates are the same fibre segment
re-photographed -> register + pointwise mean/variance.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .anomaly import detect_replicate_outliers, exclusion_reason  # noqa: E402
from .config import CONFIG  # noqa: E402
from .io_utils import natural_key, parse_name  # noqa: E402
from .register import register_sample  # noqa: E402

DEFAULT_OUT = "./fibrecv_output"


def _load_profiles(out_root: Path, cfg: CONFIG) -> tuple[dict[str, list[dict]], list[dict]]:
    """Read all per-image CSVs (+meta); bucket registrable profiles by group.

    Returns ``(groups, image_rows)``. ``groups`` holds only the replicates that
    pass ``exclusion_reason`` (band_mismatch / coverage / optionally anomaly);
    ``image_rows`` keeps one summary row for EVERY parseable image -- excluded
    ones included -- for ``summary/per_image_summary.csv``. The rows'
    ``anomaly_flags`` are still lists here; main() appends ``replicate_outlier``
    and joins them.
    """
    csv_dir = out_root / "per_image" / "csv"
    meta_dir = out_root / "per_image" / "diagnostics"
    groups: dict[str, list[dict]] = defaultdict(list)
    image_rows: list[dict] = []
    for csv in sorted(csv_dir.glob("*_profile.csv")):
        base = csv.stem[:-len("_profile")]
        try:
            group, replicate = parse_name(base)
        except ValueError:
            continue
        df = pd.read_csv(csv)
        coverage = None
        band_mismatch = False
        anomaly: dict = {}
        median_um = None
        meta_path = meta_dir / f"{base}_meta.json"
        if meta_path.exists():
            with open(meta_path) as fh:
                meta = json.load(fh)
            coverage = meta.get("coverage")
            band_mismatch = bool(meta.get("band_mismatch"))
            anomaly = meta.get("anomaly") or {}  # absent in pre-anomaly trees
            median_um = meta.get("median_diameter_um")
        flags = list(anomaly.get("flags") or [])
        reason = exclusion_reason(band_mismatch, coverage, flags, cfg)
        image_rows.append(
            {
                "name": base,
                "group": group,
                "replicate": replicate,
                "median_diameter_um": median_um,
                "coverage": coverage,
                "anomaly_flags": flags,
                "max_jump_px": anomaly.get("max_jump_px"),
                "longest_gap_frac": anomaly.get("longest_gap_frac"),
                "step_frac": anomaly.get("step_frac"),
                "rep_dev_frac": None,
                "excluded": reason is not None,
                "excluded_reason": reason or "",
            }
        )
        if reason is not None:
            continue
        groups[group].append(
            {
                "replicate": replicate,
                "coverage": coverage,
                "x": df["x_px"].to_numpy(float),
                "diameter_px_raw": df["diameter_px_raw"].to_numpy(float),
                "diameter_px_smooth": df["diameter_px_smooth"].to_numpy(float),
                "valid": df["valid"].to_numpy(bool),
            }
        )
    return groups, image_rows


def _plot_sample(table: dict, group: str, out: Path) -> None:
    """Mean diameter curve with a +/- std shaded band."""
    x = table["x_aligned_px"]
    mean = table["mean_um"]
    std = table["std_um"]
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(x, mean, "-", lw=1.5, color="tab:blue", label="mean")
    band = np.where(np.isfinite(std), std, 0.0)
    ax.fill_between(x, mean - band, mean + band, alpha=0.25, color="tab:blue", label="±std")
    ax.set_xlabel("aligned x position (px)")
    ax.set_ylabel("diameter (µm)")
    ax.set_title(f"sample {group}  (n_replicates={int(np.nanmax(table['n']))})")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Aggregate per-image profiles into registered averages.")
    ap.add_argument("--out", default=DEFAULT_OUT, help="output root (holds per_image/)")
    sel = ap.add_mutually_exclusive_group()
    sel.add_argument("--groups", nargs="+", help="group labels, e.g. 3_1 10_5")
    sel.add_argument("--all", action="store_true", help="aggregate every group found")
    ap.add_argument("--ppu", type=float, default=None)
    ap.add_argument("--max-shift", dest="max_shift", type=int, default=None)
    ap.add_argument("--min-corr", dest="min_corr", type=float, default=None)
    ap.add_argument("--min-coverage", dest="min_coverage", type=float, default=None)
    ap.add_argument("--anomaly-exclude", dest="anomaly_exclude",
                    action=argparse.BooleanOptionalAction, default=None,
                    help="image-level anomalies also exclude a replicate from "
                         "registration (default: advisory only)")
    ap.add_argument("--rep-dev-frac", dest="rep_dev_frac", type=float, default=None)
    args = ap.parse_args(argv)

    cfg = CONFIG()
    ov = {k: v for k, v in {
        "ppu": args.ppu, "max_shift": args.max_shift,
        "min_corr": args.min_corr, "min_coverage": args.min_coverage,
        "anomaly_exclude": args.anomaly_exclude, "rep_dev_frac": args.rep_dev_frac,
    }.items() if v is not None}
    cfg = replace(cfg, **ov)

    out_root = Path(args.out)
    groups, image_rows = _load_profiles(out_root, cfg)
    if args.groups:
        wanted = set(args.groups)
        groups = {g: v for g, v in groups.items() if g in wanted}

    # per-image summary before the no-groups bailout: when exclusions drop
    # every replicate this file is the audit of WHY, so it must still land.
    # It covers every image in the tree, not just --groups selections.
    if image_rows:
        medians_by_group: dict[str, dict] = defaultdict(dict)
        for r in image_rows:
            medians_by_group[r["group"]][r["name"]] = r["median_diameter_um"]
        for g, medians in medians_by_group.items():
            devs, outliers = detect_replicate_outliers(medians, cfg)
            for r in image_rows:
                if r["group"] != g:
                    continue
                r["rep_dev_frac"] = devs.get(r["name"])
                if r["name"] in outliers:
                    r["anomaly_flags"] = [*r["anomaly_flags"], "replicate_outlier"]
        for r in image_rows:
            r["anomaly_flags"] = ";".join(r["anomaly_flags"])
        (out_root / "summary").mkdir(parents=True, exist_ok=True)
        per_image = pd.DataFrame(
            sorted(image_rows, key=lambda r: natural_key(r["name"]))
        )
        per_image_path = out_root / "summary" / "per_image_summary.csv"
        per_image.to_csv(per_image_path, index=False)
        print(f"Wrote {len(image_rows)} image rows -> {per_image_path}")

    if not groups:
        print("No groups matched / no per-image CSVs found.")
        return 1

    for d in ("per_sample/csv", "per_sample/plots", "per_sample/shifts", "summary"):
        (out_root / d).mkdir(parents=True, exist_ok=True)

    rows = []
    for group in sorted(groups, key=natural_key):
        profiles = groups[group]
        if not profiles:
            continue
        table, shifts, summary = register_sample(profiles, cfg)

        # per-sample CSV
        pd.DataFrame(table).to_csv(
            out_root / "per_sample" / "csv" / f"sample_{group}_registered.csv", index=False
        )
        # per-sample plot
        _plot_sample(table, group, out_root / "per_sample" / "plots" / f"sample_{group}_registered.png")
        # shifts json
        with open(out_root / "per_sample" / "shifts" / f"sample_{group}_shifts.json", "w") as fh:
            json.dump({"group": group, "shifts": shifts, "summary": summary}, fh, indent=2)

        covs = [p["coverage"] for p in profiles if p["coverage"] is not None]
        rows.append(
            {
                "group": group,
                "mean_um": summary["mean_um"],
                "std_um": summary["std_um"],
                "CV": summary["cv"],
                "n_points": summary["n_points"],
                "mean_coverage": float(np.mean(covs)) if covs else None,
                "overlap_px": summary["overlap_px"],
                "n_replicates_used": summary["n_replicates_used"],
                "low_confidence": summary["n_replicates_used"] < 2,
                "registration_uncertain": summary["registration_uncertain"],
            }
        )
        print(f"  {group}: mean={summary['mean_um']:.2f}um std={summary['std_um']:.2f} "
              f"CV={summary['cv']:.3f} reps={summary['n_replicates_used']} "
              f"overlap={summary['overlap_px']}px")

    master = pd.DataFrame(rows).sort_values("group", key=lambda s: s.map(natural_key))
    master_path = out_root / "summary" / "master_summary.csv"
    master.to_csv(master_path, index=False)
    print(f"Wrote {len(rows)} sample rows -> {master_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
