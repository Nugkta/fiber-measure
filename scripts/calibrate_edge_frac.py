"""Calibration sweep: edge_frac on synthetic + real C1 images.

Criterion 1 — synthetic circle delta (per-side edge bias): must stay in [-1, +4] px.
Criterion 2 — real C1 width stability across the sweep.

Usage:
    uv run python scripts/calibrate_edge_frac.py
    uv run python scripts/calibrate_edge_frac.py --c1-root /path/to/C1 --n-sample 30
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from fibrecv.compute import compute_measurement
from fibrecv.config import CONFIG

W_IMG, H_IMG = 1200, 500
CIRCLE_D = 190.0


def _render_circle(seed: int = 200, noise: float = 0.03) -> np.ndarray:
    """Horizontal circle fibre (same geometry as test_xsection_synthetic)."""
    h = CIRCLE_D / 2.0
    y = np.arange(H_IMG, dtype=float)[:, None]
    dist = np.abs(y - H_IMG / 2)
    t = np.clip((h + 1.0 - dist) / 2.0, 0.0, 1.0)[..., None]
    bg = np.array([0.21, 0.20, 0.19])
    fg = np.array([0.85, 0.84, 0.83])
    img = bg + (fg - bg) * t
    rng = np.random.default_rng(seed)
    img = img + rng.normal(0.0, noise, img.shape)
    return np.clip(img, 0.0, 1.0).astype(np.float32)


def synthetic_delta(cfg: CONFIG) -> float:
    """Per-side edge bias: (measured_median - true) / 2."""
    rgb = _render_circle()
    mr = compute_measurement(rgb, cfg, "synth_circle")
    med = float(np.nanmedian(mr.res.diameter_raw))
    return (med - CIRCLE_D) / 2.0


def measure_c1_sample(c1_root: Path, cfg: CONFIG, n_sample: int = 30,
                      seed: int = 0) -> dict:
    """Measure a deterministic random sample of C1 images."""
    from fibrecv.io_utils import load_rgb
    all_imgs = sorted(p for p in c1_root.glob("C1_*_a*_part*.tiff")
                      if "s.tiff" not in p.name and "_metadata" not in p.name)
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(all_imgs), size=min(n_sample, len(all_imgs)), replace=False)
    sample = [all_imgs[i] for i in sorted(idx)]

    meds, coverages = [], []
    low_conf = 0
    for p in sample:
        rgb = load_rgb(str(p))
        mr = compute_measurement(rgb, cfg, p.stem)
        med = float(np.nanmedian(mr.res.diameter_raw))
        meds.append(med)
        coverages.append(mr.meta["coverage"])
        if mr.meta.get("low_confidence"):
            low_conf += 1

    return {
        "n": len(sample),
        "median_width_px": float(np.median(meds)),
        "std_width_px": float(np.std(meds)),
        "mean_coverage": float(np.mean(coverages)),
        "low_confidence_count": low_conf,
        "width_iqr_px": float(np.subtract(*np.percentile(meds, [75, 25]))),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--c1-root", type=Path, default=None)
    ap.add_argument("--n-sample", type=int, default=30)
    args = ap.parse_args()

    fracs = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
    cap = 0.50

    hdr = f"{'frac':>6} {'delta':>8} {'band?':>6}"
    if args.c1_root:
        hdr += f" {'med_w':>8} {'std_w':>8} {'iqr_w':>8} {'cov':>6} {'lowc':>5}"
    print(hdr)
    print("-" * len(hdr))

    results = []
    for frac in fracs:
        cfg = CONFIG(feature_mode="bright", edge_frac=frac, edge_cap=cap)
        delta = synthetic_delta(cfg)
        ok = -1.0 <= delta <= 4.0

        row = {"edge_frac": frac, "delta_px": round(delta, 3), "in_band": ok}
        line = f"{frac:>6.2f} {delta:>8.2f} {'YES' if ok else 'NO':>6}"

        if args.c1_root:
            c1 = measure_c1_sample(args.c1_root, cfg, n_sample=args.n_sample)
            row.update(c1)
            line += (f" {c1['median_width_px']:>8.1f} {c1['std_width_px']:>8.1f} "
                     f"{c1['width_iqr_px']:>8.1f} {c1['mean_coverage']:>6.1%} "
                     f"{c1['low_confidence_count']:>5d}")

        print(line)
        sys.stdout.flush()
        results.append(row)

    # Selection logic
    valid = [r for r in results if r["in_band"]]
    if not valid:
        print("\nNo edge_frac passes the delta band. Check code.")
        return 1

    print(f"\n--- Valid candidates (delta in [-1, +4] px) ---")
    for r in valid:
        print(f"  frac={r['edge_frac']:.2f}  delta={r['delta_px']:.2f} px")

    out = Path("scripts/calibration_results.json")
    out.parent.mkdir(exist_ok=True)
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
