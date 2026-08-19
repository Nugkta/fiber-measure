"""Three-way comparison of full-pipeline runs at different edge_frac.

Why this exists
---------------
`edge_frac=0.30` was shipped on the strength of a full 450-image run, but the
corrected calibration (`recalibrate_edge_frac.py`) found the pair-disagreement
criterion favours LARGER fractions, and the synthetic per-side bias falls from
+2.12 px (0.30) to +1.03 px (0.40). Neither of those criteria is the metric the
study actually reports, so this script runs the real comparison: the same
per-part ellipse residual and low-confidence counts used to justify the study,
at 0.30 vs 0.40, against the same old baseline.

Input
-----
Three pipeline output roots, each with `per_part/csv/xsec_C1_*_p*.csv` and
`summary/xsection_summary.csv`:
  OLD  legacy max(R,G,B) z-map, edge_frac=0.65, k_band=4.0
  B    median z-map, edge_frac=0.30, k_band=6.0   (the shipped run)
  C    median z-map, edge_frac=0.40, k_band=6.0   (the candidate)

Output
------
`scripts/compare_edge_frac_runs.json` and a printed table:
- per-part rms distribution (p50/p90/p95/max, n>6px) for each run;
- paired Wilcoxon + bootstrap CI for C vs B on the parts they share;
- low-confidence fiber counts and total positions;
- median measured width per run (to size the geometric shrink).

Pos
---
Study-03 reproducibility artifact. Decides edge_frac on the study's own
headline metric rather than on calibration proxies.

Usage
    uv run python scripts/compare_edge_frac_runs.py

Reminder: once I am updated, update my header comments and the folder's md.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

WORKTREES = Path("/Users/stan/Documents/UOM/spins/fiber-measure/.claude/worktrees")
RUNS = {
    "old_frac0.65": WORKTREES / "multiangle-xsection" / "fibrecv_output" / "multiangle_c1",
    "new_frac0.30": WORKTREES / "edge-criteria" / "fibrecv_output" / "multiangle_c1_s03b",
    "new_frac0.40": WORKTREES / "edge-criteria" / "fibrecv_output" / "multiangle_c1_s03c",
}
RMS_FLAG_PX = 6.0
PART_RE = re.compile(r"xsec_C1_(\d+)_p(\d+)\.csv$")
BOOT_N, BOOT_SEED = 20000, 20260819


def load_parts(root: Path) -> dict[tuple[int, int], dict]:
    out: dict[tuple[int, int], dict] = {}
    for p in sorted((root / "per_part" / "csv").glob("xsec_C1_*_p*.csv")):
        m = PART_RE.search(p.name)
        if not m:
            continue
        df = pd.read_csv(p)
        valid = df["valid"].astype(bool)
        rms = pd.to_numeric(df.loc[valid, "rms_resid_px"], errors="coerce").dropna()
        out[(int(m.group(1)), int(m.group(2)))] = {
            "n_valid_cols": int(valid.sum()),
            "rms_med_px": float(rms.median()) if len(rms) else float("nan"),
        }
    return out


def dist(vals: np.ndarray) -> dict:
    v = vals[np.isfinite(vals)]
    return {
        "n_parts": int(v.size),
        "p50": float(np.percentile(v, 50)),
        "p90": float(np.percentile(v, 90)),
        "p95": float(np.percentile(v, 95)),
        "max": float(v.max()),
        f"n_over_{RMS_FLAG_PX}px": int((v > RMS_FLAG_PX).sum()),
    }


def summary_stats(root: Path) -> dict:
    df = pd.read_csv(root / "summary" / "xsection_summary.csv")
    return {
        "n_fibers": int(len(df)),
        "n_low_confidence": int(df["low_confidence"].astype(bool).sum()),
        "low_conf_fibers": sorted(int(f) for f in df.loc[df["low_confidence"].astype(bool), "fiber"]),
        "total_positions": int(df["n_positions"].sum()),
        "axis_ratio_grand_median": float(df["axis_ratio_median"].median()),
        "A_min_um2_median": float(df["A_min_um2"].median()),
        "part_rms_med_max_px_p50": float(df["part_rms_med_max_px"].median()),
    }


def median_width_px(root: Path) -> float:
    """Median over parts of the median per-angle width, as a size reference."""
    meds = []
    for p in sorted((root / "per_part" / "csv").glob("xsec_C1_*_p*.csv")):
        df = pd.read_csv(p)
        cols = [c for c in df.columns if re.fullmatch(r"w_a\d_px", c)]
        v = pd.to_numeric(df[cols].stack(), errors="coerce").dropna()
        if len(v):
            meds.append(float(v.median()))
    return float(np.median(meds)) if meds else float("nan")


def bootstrap_median_ci(d: np.ndarray) -> tuple[float, float]:
    rng = np.random.default_rng(BOOT_SEED)
    idx = rng.integers(0, d.size, size=(BOOT_N, d.size))
    m = np.median(d[idx], axis=1)
    return tuple(float(x) for x in np.percentile(m, [2.5, 97.5]))


def paired(a: dict, b: dict, name_a: str, name_b: str) -> dict:
    """Paired change b - a over the parts present in both."""
    keys = sorted(set(a) & set(b))
    va = np.array([a[k]["rms_med_px"] for k in keys])
    vb = np.array([b[k]["rms_med_px"] for k in keys])
    ok = np.isfinite(va) & np.isfinite(vb)
    va, vb, keys = va[ok], vb[ok], [k for k, o in zip(keys, ok) if o]
    d = vb - va
    ca = np.array([a[k]["n_valid_cols"] for k in keys], dtype=float)
    cb = np.array([b[k]["n_valid_cols"] for k in keys], dtype=float)
    res = {
        "comparison": f"{name_b} minus {name_a}",
        "n_parts": int(d.size),
        "n_improved": int((d < 0).sum()),
        "n_worsened": int((d > 0).sum()),
        "median_change_px": float(np.median(d)),
        "mean_change_px": float(d.mean()),
        "column_count_change_pct_median": float(np.median((cb - ca) / np.maximum(ca, 1) * 100)),
    }
    if np.any(d != 0):
        st, p = wilcoxon(d)
        res["wilcoxon_statistic"], res["wilcoxon_p_value"] = float(st), float(p)
        res["bootstrap_median_change_ci95_px"] = bootstrap_median_ci(d)
    return res


def main() -> int:
    missing = [n for n, r in RUNS.items() if not (r / "summary" / "xsection_summary.csv").exists()]
    if missing:
        print(f"missing xsection output for: {missing}")
        return 1

    parts = {n: load_parts(r) for n, r in RUNS.items()}
    out = {"runs": {}, "paired": []}

    print(f"{'run':>14} {'p50':>7} {'p90':>7} {'p95':>7} {'max':>7} "
          f"{'>6px':>5} {'lowconf':>8} {'positions':>10} {'med_w':>8}")
    print("-" * 82)
    for name, root in RUNS.items():
        d = dist(np.array([v["rms_med_px"] for v in parts[name].values()]))
        s = summary_stats(root)
        w = median_width_px(root)
        out["runs"][name] = {"root": str(root), "rms": d, "summary": s, "median_width_px": w}
        print(f"{name:>14} {d['p50']:>7.2f} {d['p90']:>7.2f} {d['p95']:>7.2f} "
              f"{d['max']:>7.2f} {d[f'n_over_{RMS_FLAG_PX}px']:>5d} "
              f"{s['n_low_confidence']:>8d} {s['total_positions']:>10d} {w:>8.1f}")

    print()
    for a, b in [("old_frac0.65", "new_frac0.30"), ("old_frac0.65", "new_frac0.40"),
                 ("new_frac0.30", "new_frac0.40")]:
        r = paired(parts[a], parts[b], a, b)
        out["paired"].append(r)
        p = r.get("wilcoxon_p_value")
        ci = r.get("bootstrap_median_change_ci95_px")
        print(f"{r['comparison']}: median {r['median_change_px']:+.3f} px, "
              f"{r['n_improved']}/{r['n_parts']} improved"
              + (f", p={p:.3g}" if p is not None else "")
              + (f", CI95 [{ci[0]:+.3f}, {ci[1]:+.3f}]" if ci else ""))

    dst = Path("scripts/compare_edge_frac_runs.json")
    with open(dst, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n-> {dst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
