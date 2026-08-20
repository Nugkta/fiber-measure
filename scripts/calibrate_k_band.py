"""Recalibrate k_band for the median-RGB bright z-map.

The median z-map has a smaller background MAD than max(R,G,B), so z-scores
inflate and more of the defocus halo crosses the absolute k_band threshold.
On the most defocused C1 images this balloons the coarse band mask, collapsing
band_ratio = width / band_thickness below band_ratio_min and raising a false
band_mismatch — even though the measured width is correct.

Sweeps k_band over a problem set (the 6 images flagged in the s03 run) and a
control set (healthy images), reporting band_half, band_ratio and width.

Usage: uv run python scripts/calibrate_k_band.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from fibrecv.compute import compute_measurement
from fibrecv.config import CONFIG
from fibrecv.io_utils import load_rgb

C1 = Path("/Users/stan/Documents/UOM/spins/multiangle/C1")

PROBLEM = [
    "C1_01_a2_part4", "C1_01_a3_part4", "C1_01_a5_part2",
    "C1_01_a5_part3", "C1_01_a5_part4", "C1_01_a6_part3",
]
# healthy controls: the part-4 siblings that passed, plus a spread of fibers
CONTROL = [
    "C1_01_a1_part4", "C1_01_a4_part4", "C1_01_a6_part4",
    "C1_03_a1_part1", "C1_05_a3_part2", "C1_07_a2_part3",
    "C1_08_a4_part1", "C1_10_a1_part1", "C1_12_a5_part4",
    "C1_14_a6_part5", "C1_15_a2_part2",
]

K_VALUES = [4.0, 5.0, 6.0, 7.0, 8.0]


def measure(name: str, k_band: float) -> dict:
    rgb = load_rgb(str(C1 / f"{name}.tiff"))
    cfg = CONFIG(feature_mode="bright", edge_frac=0.30, edge_cap=0.50, k_band=k_band)
    mr = compute_measurement(rgb, cfg, name)
    med_px = float(np.nanmedian(mr.res.diameter_raw))
    bh = float(mr.bnd.band_half)
    return {
        "width_px": med_px,
        "band_half": bh,
        "band_ratio": med_px / (2 * bh) if bh else float("nan"),
        "coverage": float(mr.meta["coverage"]),
        "mismatch": bool(mr.meta.get("band_mismatch")),
    }


def main() -> int:
    results = {}
    for k in K_VALUES:
        prob = [measure(n, k) for n in PROBLEM]
        ctrl = [measure(n, k) for n in CONTROL]
        results[k] = {"problem": prob, "control": ctrl}

        pr = np.array([r["band_ratio"] for r in prob])
        cr = np.array([r["band_ratio"] for r in ctrl])
        pw = np.array([r["width_px"] for r in prob])
        cw = np.array([r["width_px"] for r in ctrl])
        pbh = np.array([r["band_half"] for r in prob])
        n_mm = sum(r["mismatch"] for r in prob) + sum(r["mismatch"] for r in ctrl)
        cov = np.array([r["coverage"] for r in prob + ctrl])

        print(f"k_band={k:.1f}  "
              f"problem: ratio min={pr.min():.3f} med={np.median(pr):.3f} "
              f"band_half max={pbh.max():.0f} w_med={np.median(pw):.1f}  |  "
              f"control: ratio min={cr.min():.3f} med={np.median(cr):.3f} "
              f"w_med={np.median(cw):.1f}  |  mismatch={n_mm}  cov_min={cov.min():.2f}")

    out = Path("scripts/k_band_results.json")
    with open(out, "w") as f:
        json.dump({str(k): v for k, v in results.items()}, f, indent=2)
    print(f"\nResults -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
