"""Generate annotated overlay figures for the study 03 report.

Creates two figures:
1. spot_check_overlay.png — old vs new boundaries on 4 representative images
   (defocused, sharp, dim, regression)
2. rms_comparison.png — per-fibre rms waterfall showing old vs new

Usage: uv run python scripts/report_figures.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

from fibrecv.io_utils import load_rgb

C1 = Path("/Users/stan/Documents/UOM/spins/multiangle/C1")
OLD = Path("/Users/stan/Documents/UOM/spins/fiber-measure/.claude/worktrees/"
           "multiangle-xsection/fibrecv_output/multiangle_c1/per_image/csv")
NEW = Path("fibrecv_output/multiangle_c1_s03d/per_image/csv")
OLD_SUM = Path("/Users/stan/Documents/UOM/spins/fiber-measure/.claude/worktrees/"
               "multiangle-xsection/fibrecv_output/multiangle_c1/summary/xsection_summary.csv")
NEW_SUM = Path("fibrecv_output/multiangle_c1_s03d/summary/xsection_summary.csv")
OUT = Path("/Users/stan/Documents/UOM/spins/fiber-measure/docs/report/03_edge_criteria")

CASES = [
    ("C1_01_a5_part4", "Defocused\nFibre 01 a5 part 4 — was flagged band_mismatch at k=4"),
    ("C1_08_a3_part1", "Dim\nFibre 08 a3 part 1 — previously worst rms (11.84 px)"),
    ("C1_03_a1_part1", "Sharp\nFibre 03 a1 part 1 — always good, check for regression"),
    ("C1_10_a2_part1", "Registration lock\nFibre 10 a2 part 1 — the one new low-confidence"),
]

CROP_W = 420
ZOOM_PAD = 55
OLD_COLOR = "#E8501E"
NEW_COLOR = "#00B4D8"


def overlay_panel(ax, name: str, label: str) -> None:
    rgb = load_rgb(str(C1 / f"{name}.tiff"))
    o = pd.read_csv(OLD / f"{name}_profile.csv")
    n = pd.read_csv(NEW / f"{name}_profile.csv")

    valid_n = n[n["valid"] == True]
    valid_o = o[o["valid"] == True]

    xc = int(valid_n["x_px"].median()) if len(valid_n) else int(valid_o["x_px"].median())
    x0, x1 = max(0, xc - CROP_W // 2), min(rgb.shape[1], xc + CROP_W // 2)

    ref = valid_n if len(valid_n) else valid_o
    yc = np.nanmean([ref["y_top_px"].median(), ref["y_bot_px"].median()])
    half = np.nanmedian(ref["y_bot_px"] - ref["y_top_px"]) / 2
    y0 = max(0, int(yc - half - ZOOM_PAD))
    y1 = min(rgb.shape[0], int(yc + half + ZOOM_PAD))

    ax.imshow(rgb[y0:y1, x0:x1])
    for df, colour, tag in ((o, OLD_COLOR, "Old (max V, edge_z=4)"),
                            (n, NEW_COLOR, "New (median V, frac=0.30)")):
        d = df[(df["x_px"] >= x0) & (df["x_px"] < x1)]
        xs = d["x_px"].to_numpy() - x0
        for col, lbl in (("y_top_px", tag), ("y_bot_px", None)):
            ax.plot(xs, d[col].to_numpy() - y0, color=colour, lw=1.2,
                    label=lbl, alpha=0.92)

    ow = np.nanmedian(valid_o["diameter_px_raw"]) if len(valid_o) else float("nan")
    nw = np.nanmedian(valid_n["diameter_px_raw"]) if len(valid_n) else float("nan")
    shift = (nw - ow) / 2 if np.isfinite(ow) and np.isfinite(nw) else float("nan")
    ax.set_title(label, fontsize=8, fontweight="bold", loc="left", pad=4)
    info = f"width {ow:.1f} → {nw:.1f} px" if np.isfinite(ow) else ""
    if np.isfinite(shift):
        info += f"  ({shift:+.1f} px/side)"
    ax.text(0.99, 0.04, info, transform=ax.transAxes, fontsize=7,
            ha="right", va="bottom", color="white",
            bbox=dict(facecolor="black", alpha=0.55, pad=2, edgecolor="none"))
    ax.set_xticks([])
    ax.set_yticks([])


def make_overlay():
    fig, axes = plt.subplots(len(CASES), 1, figsize=(9, 3.0 * len(CASES)))
    for ax, (name, label) in zip(np.atleast_1d(axes), CASES):
        overlay_panel(ax, name, label)
    np.atleast_1d(axes)[0].legend(loc="upper right", fontsize=7, framealpha=0.85)
    fig.tight_layout(pad=0.8)
    out = OUT / "spot_check_overlay.png"
    fig.savefig(out, dpi=155, bbox_inches="tight")
    print(f"wrote {out}")
    plt.close(fig)


def make_rms_waterfall():
    old_df = pd.read_csv(OLD_SUM)
    new_df = pd.read_csv(NEW_SUM)
    merged = old_df.merge(new_df, on="fiber", suffixes=("_old", "_new"))
    merged = merged.sort_values("part_rms_med_max_px_old", ascending=False)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(merged))
    w = 0.35

    bars_old = ax.bar(x - w/2, merged["part_rms_med_max_px_old"], w,
                      color=OLD_COLOR, alpha=0.8, label="Old pipeline")
    bars_new = ax.bar(x + w/2, merged["part_rms_med_max_px_new"], w,
                      color=NEW_COLOR, alpha=0.8, label="New pipeline")

    ax.axhline(6.0, color="#888", ls="--", lw=0.8, alpha=0.7)
    ax.text(len(merged) - 0.5, 6.15, "rms flag = 6.0 px", fontsize=7,
            ha="right", color="#888")

    for i, row in merged.reset_index(drop=True).iterrows():
        if row["low_confidence_old"]:
            ax.annotate("was flagged", (i - w/2, row["part_rms_med_max_px_old"]),
                        textcoords="offset points", xytext=(0, 5),
                        fontsize=5.5, ha="center", color=OLD_COLOR, fontweight="bold")
        if row["low_confidence_new"]:
            ax.annotate("flagged", (i + w/2, row["part_rms_med_max_px_new"]),
                        textcoords="offset points", xytext=(0, 5),
                        fontsize=5.5, ha="center", color=NEW_COLOR, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([f"F{int(f):02d}" for f in merged["fiber"]], fontsize=8)
    ax.set_ylabel("Worst-part rms (px)", fontsize=9)
    ax.set_xlabel("Fibre (sorted by old rms, descending)", fontsize=9)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_ylim(0, 13)
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.grid(axis="y", alpha=0.3, which="both")
    fig.tight_layout()
    out = OUT / "rms_waterfall.png"
    fig.savefig(out, dpi=155, bbox_inches="tight")
    print(f"wrote {out}")
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    make_overlay()
    make_rms_waterfall()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
