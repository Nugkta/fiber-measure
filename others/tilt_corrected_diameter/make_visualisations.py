"""Generate the two explanatory figures for the tilt-corrected diameter algorithm.

Run from the repo root:  uv run python others/tilt_corrected_diameter/make_visualisations.py
Outputs (next to this script):
  1_perpendicular_chord.png    -- why (y_bot - y_top) over-reads and the cos(tilt) fix
  2_axis_following_average.png -- why straight-across column averaging smears the
                                  wall and how the axis-following average avoids it

Both figures are embedded by algorithm_visualisation.html -- the
self-contained explainer page (open it in a browser; it also has an
interactive tilt demo). Re-run this script after algorithm changes so the
page's figures stay pipeline-truth.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import uniform_filter1d

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests"))

from test_edges_tilt import TRUE_W, _inclined_fibre  # noqa: E402

from fibrecv.band import locate_band  # noqa: E402
from fibrecv.config import CONFIG  # noqa: E402
from fibrecv.edges import _axis_average, detect_edges  # noqa: E402
from fibrecv.features import rgb_to_desaturation  # noqa: E402
from fibrecv.qc import run_qc  # noqa: E402

# palette (validated): blue = corrected / axis-following, orange = old / straight
BLUE = "#2a78d6"
ORANGE = "#eb6834"
INK = "#0b0b0b"
INK2 = "#52514e"
SURF = "#fcfcfb"
GRID = "#e4e3df"

plt.rcParams.update({
    "figure.facecolor": SURF,
    "axes.facecolor": SURF,
    "savefig.facecolor": SURF,
    "text.color": INK,
    "axes.edgecolor": INK2,
    "axes.labelcolor": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 11,
})

OUT = Path(__file__).resolve().parent
ANGLE = 30.0  # demo tilt for the geometry panels
W_IMG, H_IMG = 800, 700


def _geometry(angle_deg: float):
    """True centerline / boundary geometry of the synthetic fibre."""
    m = math.tan(math.radians(angle_deg))
    cth = 1.0 / math.sqrt(1.0 + m * m)
    xc, yc = W_IMG / 2.0, H_IMG / 2.0
    half_vert = (TRUE_W / 2.0) / cth  # vertical offset of each boundary line
    return m, cth, xc, yc, half_vert


# --------------------------------------------------------------------------- #
# Figure 1 -- vertical chord vs perpendicular width                            #
# --------------------------------------------------------------------------- #
def figure_1() -> None:
    rgb = _inclined_fibre(ANGLE)
    m, cth, xc, yc, half_vert = _geometry(ANGLE)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    fig.suptitle("Fix 1 — measure perpendicular to the fibre axis, not vertically",
                 fontsize=14, fontweight="bold", color=INK)

    # -- panel A: the geometry on the actual synthetic image ----------------- #
    crop = 130
    ax_a.imshow(rgb, interpolation="nearest")
    ax_a.set_xlim(xc - crop, xc + crop)
    ax_a.set_ylim(yc + crop, yc - crop)  # y down, zoomed
    ax_a.grid(False)

    xs = np.array([xc - crop, xc + crop])
    for s in (-1.0, 1.0):  # true boundary lines
        ax_a.plot(xs, m * (xs - xc) + yc + s * half_vert,
                  color=INK2, lw=1.2, ls=":", alpha=0.9)

    # old: vertical chord through the centre point
    ax_a.plot([xc, xc], [yc - half_vert, yc + half_vert],
              color=ORANGE, lw=3, solid_capstyle="butt")
    ax_a.annotate(f"vertical chord\n(y_bot − y_top) = {TRUE_W / cth:.1f} px  "
                  f"(+{100 * (1 / cth - 1):.0f} %)",
                  xy=(xc, yc - half_vert), xytext=(xc - crop + 8, yc - crop + 14),
                  color=ORANGE, fontsize=10, fontweight="bold",
                  arrowprops=dict(arrowstyle="-", color=ORANGE, lw=1))

    # new: perpendicular width through the same point (normal = (-sin, cos))
    nx, ny = -math.sin(math.radians(ANGLE)), math.cos(math.radians(ANGLE))
    off = 38  # draw it a little along the axis so both segments stay readable
    px, py = xc + off * cth, yc + off * cth * m
    ax_a.plot([px - nx * TRUE_W / 2, px + nx * TRUE_W / 2],
              [py - ny * TRUE_W / 2, py + ny * TRUE_W / 2],
              color=BLUE, lw=3, solid_capstyle="butt")
    ax_a.annotate(f"perpendicular width\n= vertical × cos θ = {TRUE_W:.0f} px  ✓",
                  xy=(px + nx * TRUE_W / 2, py + ny * TRUE_W / 2),
                  xytext=(xc + 8, yc + crop - 18),
                  color=BLUE, fontsize=10, fontweight="bold",
                  arrowprops=dict(arrowstyle="-", color=BLUE, lw=1))

    ax_a.annotate("", xy=(xc + 105, yc + 105 * m), xytext=(xc + 55, yc + 55 * m),
                  arrowprops=dict(arrowstyle="->", color=INK2, lw=1.4))
    ax_a.text(xc + 78, yc + 78 * m - 14, f"fibre axis (θ = {ANGLE:.0f}°)",
              color=INK2, fontsize=9, rotation=-ANGLE, ha="center",
              rotation_mode="anchor")
    ax_a.set_title("A — same fibre, two ways to measure", color=INK)
    ax_a.set_xlabel("x (px)")
    ax_a.set_ylabel("y (px)")

    # -- panel B: measured width vs tilt, from the real pipeline ------------ #
    cfg = CONFIG()
    angles = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
    old_px, new_px = [], []
    for a in angles:
        D, *_ = rgb_to_desaturation(_inclined_fibre(float(a)), cfg)
        band = locate_band(D, cfg)
        edg = detect_edges(D, band, cfg)
        res = run_qc(edg, band, cfg)
        v = res.valid
        old_px.append(float(np.nanmedian((edg.y_bot - edg.y_top)[v])))
        new_px.append(float(np.nanmedian(edg.diameter[v])))

    base = new_px[0]
    old_pct = [100 * (o - base) / base for o in old_px]
    new_pct = [100 * (n - base) / base for n in new_px]

    ax_b.axhline(0, color=INK2, lw=1, ls="--", alpha=0.7)
    ax_b.plot(angles, old_pct, color=ORANGE, lw=2, marker="o", ms=5,
              label="vertical chord (old)")
    ax_b.plot(angles, new_pct, color=BLUE, lw=2, marker="o", ms=5,
              label="perpendicular (fixed)")
    ax_b.text(angles[-1] + 0.6, old_pct[-1], "old", color=ORANGE,
              fontsize=10, fontweight="bold", va="center")
    ax_b.text(angles[-1] + 0.6, new_pct[-1], "fixed", color=BLUE,
              fontsize=10, fontweight="bold", va="center")
    ax_b.set_xlim(-1, 52)
    ax_b.set_title("B — measured width vs tilt (same synthetic fibre)", color=INK)
    ax_b.set_xlabel("fibre tilt θ (degrees)")
    ax_b.set_ylabel("error vs horizontal measurement (%)")
    ax_b.legend(frameon=False, loc="upper left")

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "1_perpendicular_chord.png", dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Figure 2 -- straight-across vs axis-following column averaging               #
# --------------------------------------------------------------------------- #
def figure_2() -> None:
    cfg = CONFIG()
    rgb = _inclined_fibre(ANGLE)
    D, *_ = rgb_to_desaturation(rgb, cfg)
    m, cth, xc, yc, half_vert = _geometry(ANGLE)
    xc_i = int(xc)

    old_avg = uniform_filter1d(D, size=cfg.wcol, axis=1, mode="nearest")
    new_avg = _axis_average(D, m, cfg.wcol)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    fig.suptitle("Fix 2 — the 41-column average must follow the fibre axis",
                 fontsize=14, fontweight="bold", color=INK)

    # -- panel A: which pixels get averaged together ------------------------- #
    y_wall = yc - half_vert  # top boundary row at the centre column
    ax_a.imshow(rgb, interpolation="nearest")
    ax_a.set_xlim(xc - 30, xc + 30)
    ax_a.set_ylim(y_wall + 26, y_wall - 26)
    ax_a.grid(False)
    xs = np.array([xc - 30, xc + 30])
    ax_a.plot(xs, m * (xs - xc) + y_wall, color=INK2, lw=1.2, ls=":", alpha=0.9)

    ks = np.arange(-24, 25, 4)
    ax_a.scatter(xc + ks, np.full_like(ks, y_wall, dtype=float), s=42,
                 color=ORANGE, edgecolors="white", linewidths=1.2, zorder=3,
                 label="old: same row — crosses the wall")
    ax_a.scatter(xc + ks, y_wall + m * ks, s=42, color=BLUE,
                 edgecolors="white", linewidths=1.2, zorder=3,
                 label="fixed: slides with the axis — stays on the wall")
    ax_a.legend(frameon=False, loc="upper left", fontsize=9)
    ax_a.set_title("A — pixels averaged into column x (zoom on the top wall)",
                   color=INK)
    ax_a.set_xlabel("x (px)")
    ax_a.set_ylabel("y (px)")

    # -- panel B: the averaged column profile each one produces -------------- #
    rows = np.arange(H_IMG)
    lo, hi = int(yc - half_vert - 45), int(yc + half_vert + 45)
    ax_b.plot(rows[lo:hi], D[lo:hi, xc_i], color=INK2, lw=1, alpha=0.45,
              label="raw single column")
    ax_b.plot(rows[lo:hi], old_avg[lo:hi, xc_i], color=ORANGE, lw=2,
              label="old straight average — walls smeared outward")
    ax_b.plot(rows[lo:hi], new_avg[lo:hi, xc_i], color=BLUE, lw=2,
              label="axis-following average — walls stay sharp")
    for s in (-1.0, 1.0):
        ax_b.axvline(yc + s * half_vert, color=INK2, lw=1, ls="--", alpha=0.7)
    ax_b.text(yc, ax_b.get_ylim()[1] * 0.02, "true walls (dashed)",
              color=INK2, fontsize=9, ha="center")
    ax_b.set_title("B — resulting profile of one column (θ = 30°)", color=INK)
    ax_b.set_xlabel("image row y (px)")
    ax_b.set_ylabel("desaturation D (z units)")
    ax_b.legend(frameon=False, loc="upper right", fontsize=9)

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "2_axis_following_average.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    figure_1()
    figure_2()
    print("wrote", OUT / "1_perpendicular_chord.png")
    print("wrote", OUT / "2_axis_following_average.png")
