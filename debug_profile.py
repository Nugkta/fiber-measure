"""Scratch debug: windowed desaturation profiles for a named image.

Dependencies: numpy, matplotlib(Agg), fibrecv. Inputs: image stem via argv[1],
optional column list via argv[2] (comma-separated).
Output: PNG of windowed D[:,x] with the edge level, window extent and detected
band mask marked, plus printed window/bg diagnostics per column.
Pos: throwaway diagnostic for tuning edges.py; not part of the pipeline.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d, uniform_filter1d

from fibrecv import band as band_mod
from fibrecv import edges as edges_mod
from fibrecv import features as feat_mod
from fibrecv import io_utils
from fibrecv.config import CONFIG

cfg = CONFIG()
stem = sys.argv[1] if len(sys.argv) > 1 else "masp2 10_5_1"
cols = [int(c) for c in sys.argv[2].split(",")] if len(sys.argv) > 2 else [400, 900, 1300, 1800, 2200]

path = f"/net/scratch/j56806hx/spins-cv/Images MasP2/{stem}.jpg"
rgb = io_utils.load_rgb(path)
D, S, s_bg, mad = feat_mod.rgb_to_desaturation(rgb, cfg)
bnd = band_mod.locate_band(D, cfg)
hw = edges_mod.half_window_px(bnd, cfg, D.shape[0])
print(f"{stem}: s_bg={s_bg:.4f} mad={mad:.5f} band_half={bnd.band_half:.1f} "
      f"hw={hw} span=({bnd.x0},{bnd.x1}) slope={bnd.slope:.4f} lowconf={bnd.low_confidence}")

Davg = uniform_filter1d(D, size=cfg.wcol, axis=1, mode="nearest")
Dsm = gaussian_filter1d(Davg, sigma=cfg.sigma_y, axis=0, mode="nearest")

fig, axes = plt.subplots(len(cols), 1, figsize=(9, 2.6 * len(cols)))
for ax, x in zip(np.atleast_1d(axes), cols):
    c = bnd.c_fit[x]
    y0 = int(np.clip(round(c - hw), 0, D.shape[0] - 1))
    y1 = int(np.clip(round(c + hw), 0, D.shape[0] - 1))
    d = Dsm[y0:y1 + 1, x]
    yy = np.arange(y0, y1 + 1)
    g = max(1, min(cfg.guard, d.size // 2))
    bg_top = float(np.median(d[:g])); bg_bot = float(np.median(d[-g:]))
    local_bg = min(bg_top, bg_bot)
    A = float(d.max()) - local_bg
    level = local_bg + min(cfg.edge_z, cfg.edge_frac * A)
    ax.plot(yy, d, "-", color="tab:blue")
    ax.axhline(level, color="tab:red", ls=":", lw=1.0)
    ax.axhline(local_bg, color="k", ls=":", lw=0.6)
    ax.axvline(c, color="k", ls="--", lw=0.7)
    bm = bnd.mask[y0:y1 + 1, x]
    if bm.any():
        ax.plot(yy[bm], np.full(bm.sum(), -1.0), "s", color="tab:green", ms=2)
    ax.set_title(f"x={x} win=({y0},{y1}) bg_top={bg_top:.1f} bg_bot={bg_bot:.1f} "
                 f"A={A:.1f} level={level:.1f}", fontsize=9)
    ax.set_ylabel("D")
    print(f"  x={x}: win=({y0},{y1}) bg_top={bg_top:.2f} bg_bot={bg_bot:.2f} "
          f"peak={d.max():.2f} level={level:.2f}")
fig.tight_layout()
out = f"/net/scratch/j56806hx/spins-cv/output/_hardcases/_D_{stem.replace(' ', '_')}.png"
fig.savefig(out, dpi=105)
print("wrote", out)
