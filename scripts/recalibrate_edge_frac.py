"""Re-calibrate edge_frac, replacing the two broken criteria of the first pass.

Why this exists
---------------
The original `calibrate_edge_frac.py` selected edge_frac=0.30 on two criteria
that a code review showed were both invalid:

  * Criterion 2 (synthetic delta) was measured on `_render_circle`, which
    silently produced a 500x1 image -- W_IMG was declared but never broadcast.
    On one column the delta is seed-unstable (old: 4.36 +/- 0.33, range
    3.87-5.10); at full width sd drops to ~0.006.
  * Criterion 1 was `low_confidence_count`, which is 0 at every edge_frac --
    zero discriminating power. Its 30-image sample also excluded all six
    defocused images and fiber 07 entirely, i.e. it excluded the phenomenon
    being calibrated.

This script fixes both:
  * the synthetic is rendered at full width and averaged over seeds;
  * defocus stability is measured the way the study plan actually specified --
    as within-pair width disagreement. The six angles are 60deg-spaced
    (0/60/120/180/240/300), so mod 180 they form three pairs (a1,a4), (a2,a5),
    (a3,a6) that image the SAME direction and must therefore agree. Pair
    disagreement is a direct, ground-truth-free measure of defocus-induced
    width error, evaluated on parts that contain known-defocused angles.

Usage: uv run python scripts/recalibrate_edge_frac.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from fibrecv.compute import compute_measurement
from fibrecv.config import CONFIG
from fibrecv.io_utils import load_rgb

C1 = Path("/Users/stan/Documents/UOM/spins/multiangle/C1")

W_IMG, H_IMG = 1200, 500
CIRCLE_D = 190.0
PAIRS = [(1, 4), (2, 5), (3, 6)]

# Parts carrying known-defocused angles (fiber 01 p4 has a2/a3/a5 defocused;
# p2/p3 carry the other flagged images), plus sharp controls from other fibers.
DEFOCUSED_PARTS = [(1, 4), (1, 3), (1, 2)]
CONTROL_PARTS = [(3, 1), (14, 5), (10, 2)]

FRACS = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
SEEDS = list(range(200, 212))


def _render_circle(seed: int, noise: float = 0.03) -> np.ndarray:
    """Full-width circular fibre. NB the original broadcast only to (H,1,3)."""
    h = CIRCLE_D / 2.0
    y = np.arange(H_IMG, dtype=float)[:, None]
    dist = np.abs(y - H_IMG / 2.0)
    t = np.clip((h + 1.0 - dist) / 2.0, 0.0, 1.0)
    t = np.broadcast_to(t, (H_IMG, W_IMG))[..., None]
    bg = np.array([0.21, 0.20, 0.19])
    fg = np.array([0.85, 0.84, 0.83])
    img = bg + (fg - bg) * t
    rng = np.random.default_rng(seed)
    img = img + rng.normal(0.0, noise, img.shape)
    return np.clip(img, 0.0, 1.0).astype(np.float32)


def synthetic_delta(cfg: CONFIG) -> tuple[float, float]:
    """Per-side edge bias, mean +/- sd over seeds, at full width."""
    ds = []
    for s in SEEDS:
        mr = compute_measurement(_render_circle(s), cfg, f"synth{s}")
        ds.append((float(np.nanmedian(mr.res.diameter_raw)) - CIRCLE_D) / 2.0)
    return float(np.mean(ds)), float(np.std(ds))


def _profile(fiber: int, angle: int, part: int, cfg: CONFIG):
    p = C1 / f"C1_{fiber:02d}_a{angle}_part{part}.tiff"
    mr = compute_measurement(load_rgb(str(p)), cfg, p.stem)
    w = np.where(mr.res.valid, mr.res.diameter_smooth, np.nan)
    return w, mr.meta["coverage"]


def pair_disagreement(parts, cfg: CONFIG) -> dict:
    """Median |w_a - w_b| over the three same-direction angle pairs.

    Profiles are compared on the overlapping x-range without cross-angle
    registration: stage-1 output is in image coordinates and the same part is
    imaged at the same stage position, so a shared x index is the same place
    along the fibre up to the (small) inter-angle repositioning. This is a
    relative comparison across edge_frac, so any residual misregistration is
    common-mode and cancels in the sweep.
    """
    out = []
    for fib, prt in parts:
        for a, b in PAIRS:
            try:
                wa, _ = _profile(fib, a, prt, cfg)
                wb, _ = _profile(fib, b, prt, cfg)
            except Exception:
                continue
            n = min(wa.size, wb.size)
            d = np.abs(wa[:n] - wb[:n])
            d = d[np.isfinite(d)]
            if d.size:
                out.append(float(np.median(d)))
    return {"median": float(np.median(out)) if out else float("nan"),
            "n_pairs": len(out), "per_pair": out}


def main() -> int:
    results = []
    print(f"{'frac':>6} {'delta_px':>16} {'defocus_pair_dis':>17} "
          f"{'control_pair_dis':>17}")
    print("-" * 60)
    for frac in FRACS:
        cfg = CONFIG(feature_mode="bright", edge_frac=frac,
                     edge_cap=0.50, k_band=6.0)
        dm, ds = synthetic_delta(cfg)
        dfc = pair_disagreement(DEFOCUSED_PARTS, cfg)
        ctl = pair_disagreement(CONTROL_PARTS, cfg)
        results.append({"edge_frac": frac, "delta_mean": dm, "delta_sd": ds,
                        "defocus_pair_median": dfc["median"],
                        "defocus_n": dfc["n_pairs"],
                        "control_pair_median": ctl["median"],
                        "control_n": ctl["n_pairs"]})
        print(f"{frac:>6.2f} {dm:>9.3f}+/-{ds:<5.3f} "
              f"{dfc['median']:>17.3f} {ctl['median']:>17.3f}")

    out = Path("scripts/recalibrate_edge_frac.json")
    with open(out, "w") as f:
        json.dump({"seeds": SEEDS, "defocused_parts": DEFOCUSED_PARTS,
                   "control_parts": CONTROL_PARTS, "results": results}, f, indent=2)
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
