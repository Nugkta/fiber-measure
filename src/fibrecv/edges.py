"""Per-column tight-inner-edge detection (the shadow-critical recipe).

Dependencies
------------
``numpy``, ``scipy.ndimage`` (uniform_filter1d, gaussian_filter1d).

Inputs
------
- ``D``: desaturation z-map (H, W) from ``features``.
- ``BandResult`` (centerline + band thickness) from ``band``.
- ``CONFIG`` for ``wcol``, ``sigma_y``, ``edge_frac``, ``guard``, ``amin``,
  ``k_band``, ``wall_gap_frac`` and window sizing.

Output
------
``EdgeResult`` with float arrays ``y_top``, ``y_bot`` (sub-pixel boundary rows,
NaN where invalid), ``diameter`` = ``y_bot - y_top``, plus per-column ``amp``
(band amplitude A), ``y_core`` and integer ``flags`` for QC.

Signal model & strictness logic
-------------------------------
Each fibre column is a broad desaturation **hump**: D rises from ~0 (saturated
pink background) up to a peak of order 10 over the fibre body, with small
specular wiggles near the top. A shadow band beside the fibre shows up as an
intermediate plateau partway up one flank.

The two failure mechanisms this recipe defends against (both observed in the
``10_5`` replicate group):
  1. **Internal reflections** (bright specular rim/iridescence inside the tube)
     create steep desaturation flanks *inside* the fibre that can out-compete
     the true wall -- naive gradient-max or fraction-of-peak levels get dragged
     inward (false narrowing).
  2. **Shadow / vignette ramps** alongside the fibre desaturate the background
     gently over tens of px -- naive low absolute thresholds (or backgrounds
     estimated at far window extremes, which vignetting can push *negative*)
     place the boundary far outside the fibre (false widening).

The discriminator, measured empirically: true fibre walls rise steeply
(~0.2-0.5 z/px after smoothing) while shadow/vignette ramps are gentle
(~0.03 z/px). Hence, per column (a ``wcol`` (default 41, >=15) neighbourhood
average of ``D`` to suppress iridescence banding, vertically smoothed):
  * per side, find the **outermost steep rising wall**: the first run,
    scanning inward from the window extreme, whose slope toward the core
    exceeds ``max(slope_min, slope_rel * side_max_slope)`` and whose total rise
    is at least ``rise_min`` z-units (so noise blips never qualify). The run
    may bridge short flat plateaus (``wall_gap_frac`` of band thickness) so a
    defocus-fragmented soft wall is not rejected piecewise, but a falling
    stretch still ends it. Outermost beats reflection flanks (they are inner);
    the slope+rise gates exclude shadow ramps (too gentle).
  * estimate ``base`` = median of ``D`` over the ``guard`` px just *outside*
    that wall (so an adjacent shadow plateau becomes the reference, per the
    strict never-count-shadow rule), and the wall amplitude
    ``A_side = D[wall top] - base``,
  * place the boundary at the level ``base + min(edge_z, edge_frac*A_side)``
    as a sub-pixel crossing searched **on the wall run itself** -- it cannot
    land on a shadow shoulder or an interior reflection.
  * an outside-recovery guard flags columns where no true background
    (``D < k_band/2``) exists beyond the wall (fibre fills the window).

``edge_z`` (z-units above the wall-local base) is the strictness knob: higher
-> the crossing sits higher up the wall -> tighter. ``edge_frac`` caps the
level for weak walls so the crossing stays on them.

Pos
---
Fourth stage of the per-image pipeline. Consumes ``features`` + ``band``,
produces the raw diameter profile that ``qc`` cleans and ``register`` aggregates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter1d, uniform_filter1d

from .band import BandResult
from .config import CONFIG

# bit flags recorded per column (0 == clean)
FLAG_OK = 0
FLAG_NO_CORE = 1        # amplitude below amin / no band in window
FLAG_NO_BG = 2          # no true background (D<k_band/2) beyond a boundary -> shadow/fill risk
FLAG_TOP_CLIP = 4       # top edge hit the search-window boundary
FLAG_BOT_CLIP = 8       # bottom edge hit the search-window boundary
FLAG_BAD_GRAD = 16      # could not find a usable crossing


@dataclass
class EdgeResult:
    y_top: np.ndarray     # float (W,) sub-pixel top boundary, NaN if invalid
    y_bot: np.ndarray     # float (W,) sub-pixel bottom boundary, NaN if invalid
    diameter: np.ndarray  # float (W,) y_bot - y_top, NaN if invalid
    amp: np.ndarray       # float (W,) band amplitude A (z units)
    y_core: np.ndarray    # float (W,) core row (argmax of smoothed D)
    flags: np.ndarray     # int   (W,) per-column QC bit flags
    half_window: int      # vertical search half-window actually used (px)


def half_window_px(band: BandResult, cfg: CONFIG, H: int) -> int:
    """Vertical search half-window from band thickness (clamped to image)."""
    hw = int(cfg.window_thick_mult * 2.0 * band.band_half + cfg.window_pad + 3 * cfg.guard)
    return int(np.clip(hw, 8, H // 2))


def _outer_crossing(d: np.ndarray, y_core: int, y_outer: int, level: float) -> float | None:
    """Sub-pixel index of the outermost crossing of ``level``.

    Scans from ``y_outer`` (background side) toward ``y_core`` and returns the
    first index where ``d`` rises through ``level`` (background -> fibre),
    linearly interpolated. Returns None if there is no such crossing.
    """
    step = 1 if y_core > y_outer else -1
    prev = y_outer
    for i in range(y_outer + step, y_core + step, step):
        a, b = d[prev], d[i]
        if a < level <= b or a <= level < b:
            if b != a:
                frac = (level - a) / (b - a)
                return prev + step * frac
            return float(i)
        prev = i
    return None


def _find_wall(
    d: np.ndarray, g_in: np.ndarray, outer: int, core: int, step: int, cfg: CONFIG,
    gap: int = 0,
) -> tuple[int, int] | None:
    """Outermost steep rising wall on one side of the core.

    ``g_in`` is the profile slope *toward the core* (positive = rising into the
    fibre); we walk from the window extreme ``outer`` toward ``core`` and return
    the first run with slope >= max(slope_min, slope_rel*side_max) (sustained
    above half that) whose total rise is >= ``rise_min`` z-units -- noise blips
    fail the rise test, shadow/vignette ramps fail the slope test.

    ``gap`` px of flat (not falling) profile may be bridged inside a run: a
    defocus-softened wall is a long gentle ramp broken by small plateaus, and
    without bridging each fragment fails ``rise_min`` so the finder hops inward
    onto internal reflections (observed across the MasP2 set). A falling
    stretch (slope < -g_lo) still terminates the run, so distinct features are
    never merged. ``gap=0`` reproduces the strict contiguous behaviour.
    Returns ``(i_outer, i_inner)`` window indices, or None.
    """
    idxs = list(range(outer, core, step))
    if len(idxs) < 3:
        return None
    gvals = g_in[idxs]
    gmax = float(gvals.max())
    if gmax <= 0:
        return None
    # relative gate capped at slope_cap: a defocus-softened true wall still
    # qualifies even when a sharp internal reflection dominates the side max
    g_hi = max(cfg.slope_min, min(cfg.slope_rel * gmax, cfg.slope_cap))
    g_lo = 0.5 * g_hi

    k, n = 0, len(idxs)
    while k < n:
        if gvals[k] >= g_hi:
            j = k
            last_good = k     # last genuinely-rising index of the run
            pending = 0       # px spent inside the current sub-threshold gap
            while j + 1 < n:
                nxt = gvals[j + 1]
                if nxt >= g_lo:
                    j += 1
                    last_good = j
                    pending = 0
                elif nxt > -g_lo and pending < gap:
                    j += 1    # bridge a flat (not falling) stretch
                    pending += 1
                else:
                    break
            i_outer, i_inner = idxs[k], idxs[last_good]
            if d[i_inner] - d[i_outer] >= cfg.rise_min:
                return i_outer, i_inner
            k = j + 1
        else:
            k += 1
    return None


def _side_edge(
    d: np.ndarray, g_in: np.ndarray, outer: int, core_idx: int, step: int, cfg: CONFIG,
    gap: int = 0,
) -> tuple[float, int]:
    """One boundary (top side: step=+1, bottom side: step=-1) -> (edge, flag)."""
    wall = _find_wall(d, g_in, outer, core_idx, step, cfg, gap=gap)
    if wall is None:
        return np.nan, FLAG_BAD_GRAD
    i_outer, i_inner = wall
    flag = FLAG_OK

    # wall-local background: guard px just outside the wall (a shadow plateau
    # there becomes the reference level -> shadow never counted as fibre)
    if step > 0:
        lo, hi = max(0, i_outer - cfg.guard), i_outer
        outer_zone = d[: i_outer + 1]
    else:
        lo, hi = i_outer + 1, min(d.size, i_outer + 1 + cfg.guard)
        outer_zone = d[i_outer:]
    base_val = float(np.median(d[lo:hi])) if hi > lo else float(d[i_outer])

    A_side = float(d[i_inner]) - base_val
    if A_side <= cfg.eps:
        return np.nan, FLAG_BAD_GRAD

    level = base_val + min(cfg.edge_z, cfg.edge_frac * A_side)
    scan_outer = int(np.clip(i_outer - step * cfg.guard, 0, d.size - 1))
    edge = _outer_crossing(d, i_inner, scan_outer, level)
    if edge is None:
        edge = float(i_outer)  # conservative fallback: outer base of the wall

    # outside-recovery guard: genuine background must exist beyond the wall
    if outer_zone.size == 0 or float(np.min(outer_zone)) >= cfg.k_band / 2.0:
        flag |= FLAG_NO_BG
    return edge, flag


def _detect_column(d: np.ndarray, g: np.ndarray, cfg: CONFIG,
                   gap: int = 0) -> tuple[float, float, float, float, int]:
    """Detect (top, bottom, amp, core_idx, flag) for one windowed column profile.

    ``d``/``g`` are the smoothed desaturation profile and its y-gradient over
    the search window; indices are window-local (0 = top of window).
    """
    n = d.size
    gk = max(1, min(cfg.guard, n // 2))
    core_idx = int(np.argmax(d))
    core_val = float(d[core_idx])

    # lenient amplitude gate vs the window extremes (vignetting only ever
    # lowers the extremes' median, so a real fibre always passes)
    rough_bg = float(min(np.median(d[:gk]), np.median(d[-gk:])))
    A = core_val - rough_bg
    if A < cfg.amin:
        return np.nan, np.nan, A, float(core_idx), FLAG_NO_CORE

    top, ft = _side_edge(d, g, 0, core_idx, +1, cfg, gap=gap)       # slope toward core = +g
    bot, fb = _side_edge(d, -g, n - 1, core_idx, -1, cfg, gap=gap)  # slope toward core = -g
    flag = ft | fb
    if np.isnan(top) or np.isnan(bot):
        return np.nan, np.nan, A, float(core_idx), flag | FLAG_BAD_GRAD
    if top <= 1:
        flag |= FLAG_TOP_CLIP
    if bot >= n - 2:
        flag |= FLAG_BOT_CLIP
    return top, bot, A, float(core_idx), flag


def detect_edges(D: np.ndarray, band: BandResult, cfg: CONFIG) -> EdgeResult:
    """Run the tight-inner-edge recipe over every column in the band span."""
    H, W = D.shape
    hw = half_window_px(band, cfg, H)

    # plateau-bridging length for the wall finder, scaled to the fibre size so
    # thin fibres never bridge across their own thickness (bounds 4..16 px)
    bh = band.band_half if np.isfinite(band.band_half) else 0.0
    gap = int(np.clip(round(cfg.wall_gap_frac * 2.0 * bh), 4, 16))

    # column-neighbourhood average, then vertical smoothing + gradient (vectorised)
    Davg = uniform_filter1d(D, size=max(1, cfg.wcol), axis=1, mode="nearest")
    Dsm = gaussian_filter1d(Davg, sigma=cfg.sigma_y, axis=0, mode="nearest")
    G = np.gradient(Dsm, axis=0)

    y_top = np.full(W, np.nan, dtype=np.float32)
    y_bot = np.full(W, np.nan, dtype=np.float32)
    amp = np.full(W, np.nan, dtype=np.float32)
    y_core_arr = np.full(W, np.nan, dtype=np.float32)
    flags = np.zeros(W, dtype=np.int32)

    for x in range(band.x0, band.x1 + 1):
        c = band.c_fit[x]
        y_lo = int(np.clip(round(c - hw), 0, H - 1))
        y_hi = int(np.clip(round(c + hw), 0, H - 1))
        if y_hi - y_lo < 8:
            flags[x] = FLAG_NO_CORE
            continue

        d = Dsm[y_lo:y_hi + 1, x]
        g = G[y_lo:y_hi + 1, x]
        top, bot, A, core_idx, flag = _detect_column(d, g, cfg, gap=gap)
        flags[x] = flag
        amp[x] = A
        y_core_arr[x] = y_lo + core_idx
        if np.isnan(top) or np.isnan(bot):
            continue
        y_top[x] = y_lo + top
        y_bot[x] = y_lo + bot

    diameter = y_bot - y_top
    return EdgeResult(
        y_top=y_top,
        y_bot=y_bot,
        diameter=diameter,
        amp=amp,
        y_core=y_core_arr,
        flags=flags,
        half_window=hw,
    )
