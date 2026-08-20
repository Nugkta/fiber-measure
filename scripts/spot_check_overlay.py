"""Visual spot-check: old vs new detected boundaries on representative C1 images.

Reads the y_top_px / y_bot_px columns from the study-02 (old) and study-03b
(new, k_band=6) per-image profile CSVs and draws both boundary pairs on a
zoomed crop of the source image, so the edge shift is directly visible.

Usage: uv run python scripts/spot_check_overlay.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fibrecv.io_utils import load_rgb

C1 = Path("/Users/stan/Documents/UOM/spins/multiangle/C1")
OLD = Path("/Users/stan/Documents/UOM/spins/fiber-measure/.claude/worktrees/"
           "multiangle-xsection/fibrecv_output/multiangle_c1/per_image/csv")
NEW = Path("fibrecv_output/multiangle_c1_s03b/per_image/csv")

# 1 formerly-bad (defocused, was falsely flagged), 2 good, 1 dim
CASES = [
    ("C1_01_a5_part4", "formerly-bad (defocused; false band_mismatch at k_band=4)"),
    ("C1_03_a1_part1", "good / sharp"),
    ("C1_14_a6_part5", "good / sharp"),
    ("C1_04_a1_part1", "dim (lowest coverage in set)"),
]

CROP_W = 420   # px of x-range shown
ZOOM_PAD = 55  # px above/below the fibre


def panel(ax, name: str, label: str) -> None:
    rgb = load_rgb(str(C1 / f"{name}.tiff"))
    o = pd.read_csv(OLD / f"{name}_profile.csv")
    n = pd.read_csv(NEW / f"{name}_profile.csv")

    # centre the crop on the densest measured stretch
    valid = n[n["valid"] == True]
    xc = int(valid["x_px"].median())
    x0, x1 = max(0, xc - CROP_W // 2), min(rgb.shape[1], xc + CROP_W // 2)

    yc = np.nanmean([valid["y_top_px"].median(), valid["y_bot_px"].median()])
    half = np.nanmedian(valid["y_bot_px"] - valid["y_top_px"]) / 2
    y0 = max(0, int(yc - half - ZOOM_PAD))
    y1 = min(rgb.shape[0], int(yc + half + ZOOM_PAD))

    ax.imshow(rgb[y0:y1, x0:x1])
    for df, colour, tag in ((o, "#E8501E", "old (max-V, ez4)"),
                            (n, "#00B4D8", "new (median-V, frac .30)")):
        d = df[(df["x_px"] >= x0) & (df["x_px"] < x1)]
        xs = d["x_px"].to_numpy() - x0
        for col, lbl in (("y_top_px", tag), ("y_bot_px", None)):
            ax.plot(xs, d[col].to_numpy() - y0, color=colour, lw=1.3,
                    label=lbl, alpha=0.95)

    ow = np.nanmedian(o[o["valid"] == True]["diameter_px_raw"])
    nw = np.nanmedian(valid["diameter_px_raw"])
    ax.set_title(f"{name}\n{label}\nwidth {ow:.1f} → {nw:.1f} px "
                 f"({nw - ow:+.1f}, {(nw - ow) / 2:+.1f}/side)", fontsize=8.5)
    ax.set_xticks([]); ax.set_yticks([])


def main() -> int:
    fig, axes = plt.subplots(len(CASES), 1, figsize=(9, 3.1 * len(CASES)))
    for ax, (name, label) in zip(np.atleast_1d(axes), CASES):
        panel(ax, name, label)
    np.atleast_1d(axes)[0].legend(loc="upper right", fontsize=7.5, framealpha=0.9)
    fig.suptitle("Study 03 — detected boundaries, old vs new", fontsize=11, y=0.997)
    fig.tight_layout()
    out = Path("../../../docs/report/03_edge_criteria/spot_check_overlay.png").resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=155, bbox_inches="tight")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
