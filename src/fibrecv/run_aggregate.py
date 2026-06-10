"""CLI: group per-image profiles by A_B and build registered averages.

Dependencies
------------
``numpy``, ``pandas``, ``matplotlib`` (Agg), the ``register`` + ``io_utils`` +
``config`` modules. Reads the CSV/meta artifacts produced by ``run_measure``.

Inputs
------
Command-line flags: ``--out`` (output root holding ``per_image/``), selector
(``--groups`` / ``--all``), and a few CONFIG overrides (``--ppu``, ``--max-shift``,
``--min-corr``, ``--min-coverage``).

Output
------
Per A_B sample: ``per_sample/csv/sample_<A_B>_registered.csv`` (mean+/-std curve),
``per_sample/plots/sample_<A_B>_registered.png`` (mean + +/-std band),
``per_sample/shifts/sample_<A_B>_shifts.json``, and one row per sample in
``summary/master_summary.csv``.

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

from .config import CONFIG  # noqa: E402
from .io_utils import parse_name  # noqa: E402
from .register import register_sample  # noqa: E402

DEFAULT_OUT = "/net/scratch/j56806hx/spins-cv/output"


def _load_profiles(out_root: Path, cfg: CONFIG) -> dict[str, list[dict]]:
    """Read all per-image CSVs (+meta) and bucket valid profiles by A_B group."""
    csv_dir = out_root / "per_image" / "csv"
    meta_dir = out_root / "per_image" / "diagnostics"
    groups: dict[str, list[dict]] = defaultdict(list)
    for csv in sorted(csv_dir.glob("*_profile.csv")):
        base = csv.stem[:-len("_profile")]
        try:
            group, replicate = parse_name(base)
        except ValueError:
            continue
        df = pd.read_csv(csv)
        coverage = None
        meta_path = meta_dir / f"{base}_meta.json"
        if meta_path.exists():
            with open(meta_path) as fh:
                meta = json.load(fh)
            coverage = meta.get("coverage")
            if meta.get("band_mismatch"):
                continue  # drop replicate whose detector locked inside a blur band
        if coverage is not None and coverage < cfg.min_coverage:
            continue  # drop replicate that fails coverage QC
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
    return groups


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
    sel.add_argument("--groups", nargs="+", help="A_B groups, e.g. 3_1 10_5")
    sel.add_argument("--all", action="store_true", help="aggregate every group found")
    ap.add_argument("--ppu", type=float, default=None)
    ap.add_argument("--max-shift", dest="max_shift", type=int, default=None)
    ap.add_argument("--min-corr", dest="min_corr", type=float, default=None)
    ap.add_argument("--min-coverage", dest="min_coverage", type=float, default=None)
    args = ap.parse_args(argv)

    cfg = CONFIG()
    ov = {k: v for k, v in {
        "ppu": args.ppu, "max_shift": args.max_shift,
        "min_corr": args.min_corr, "min_coverage": args.min_coverage,
    }.items() if v is not None}
    cfg = replace(cfg, **ov)

    out_root = Path(args.out)
    groups = _load_profiles(out_root, cfg)
    if args.groups:
        wanted = set(args.groups)
        groups = {g: v for g, v in groups.items() if g in wanted}
    if not groups:
        print("No groups matched / no per-image CSVs found.")
        return 1

    for d in ("per_sample/csv", "per_sample/plots", "per_sample/shifts", "summary"):
        (out_root / d).mkdir(parents=True, exist_ok=True)

    rows = []
    for group in sorted(groups, key=lambda g: tuple(int(t) for t in g.split("_"))):
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

    master = pd.DataFrame(rows).sort_values("group", key=lambda s: s.map(
        lambda g: tuple(int(t) for t in g.split("_"))))
    master_path = out_root / "summary" / "master_summary.csv"
    master.to_csv(master_path, index=False)
    print(f"Wrote {len(rows)} sample rows -> {master_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
