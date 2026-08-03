"""Erf edge refinement: refit each detected wall as a blurred step, shift to its midpoint.

Dependencies
------------
``numpy``, ``scipy.optimize.curve_fit`` (the per-block least-squares fit) and
``scipy.special.erf`` (the model). Consumes ``band`` (``BandResult``,
``tilt_geometry``) and ``edges`` (``EdgeResult``, ``FLAG_OK``).

Inputs
------
- ``D``: desaturation z-map (H, W) from ``features``, used *unsmoothed* --
  pre-smoothing would inflate the fitted sigma; the block-mean + LSQ fit over
  ~86-127 samples is the noise handling instead.
- ``EdgeResult`` from ``edges.detect_edges``: the legacy per-column boundaries
  this stage refines, or falls back to unchanged.
- ``BandResult`` (centerline, tilt slope, band_half) from ``band``.
- ``CONFIG`` for ``refine_on`` and the ``refine_*`` block/window/gate
  parameters.

Output
------
``refine_edges(D, edg, bnd, cfg)`` -> ``(EdgeResult, RefineResult)``. When
``cfg.refine_on`` is False, the input ``EdgeResult`` is returned unchanged --
the SAME object, not a copy -- alongside an empty ``RefineResult``: bit-
identical off-path. ``RefineResult`` carries, per column: which columns were
actually refined (``refined_top``/``refined_bot``), the fitted blur widths
(``sigma_top``/``sigma_bot``), relative fit residuals
(``resid_top``/``resid_bot``), the applied perpendicular offsets
(``o_top``/``o_bot``) and the attempted/passing block counts.

Erf model & gates
-----------------
Each wall is modelled as a Gaussian-blurred step along the perpendicular
profile coordinate ``t`` (0 at the legacy edge, +t into the fibre):
``D(t) = a + (b - a) * 0.5 * (1 + erf((t - t0) / (sigma * sqrt(2))))``.
By the knife-edge principle the fitted midpoint ``t0`` is invariant to the
optical blur sigma, unlike the legacy fixed-level crossing (``edges.py``),
which drifts with focus. A block of ``refine_block`` columns is fit at once
(the mean profile over the block) rather than column-by-column, since ``D``
is never pre-smoothed here. A block's fit is accepted only when: ``b - a >
0`` (a genuine rising step); ``rms_residual / (b - a) < refine_relmax``;
``sigma`` within ``[refine_sigma_min, refine_sigma_max]``; and ``|t0| <=
refine_maxshift``. Accepted offset/sigma/residual values are interpolated
between passing block centres (chains break where the gap exceeds
``refine_gap_blocks`` blocks; no extrapolation beyond the outermost passing
centre). Offsets are applied only to anchor columns -- finite ``y_top`` &
``y_bot`` AND ``flags == FLAG_OK`` -- within interpolation coverage; flagged
columns never move.

Sampling & frames
-----------------
Blocks are the consecutive full ``refine_block``-column runs starting at
``bnd.x0`` (a trailing partial block is skipped) and a block is fitted only if
at least 70% of its columns are usable -- anchor columns whose whole sample
window stays inside the image (border-clipped columns are dropped *before*
that quorum test, so replicated border rows never enter a profile). The
profile runs from ``-refine_out`` (background side) to ``+inside`` px, where
``inside = clip(refine_in_frac * band_half * cth, 8, refine_in_max)``, in
0.5 px steps: ~86-127 samples per fit. Everything is measured in the frame
perpendicular to the fibre axis (``m, cth = band.tilt_geometry(bnd.slope)``,
the single tilt source): rows are sampled at ``y = anchor +- t/cth``, the
fitted ``t0``/``sigma`` are already perpendicular px (so the gates apply
unconverted) and the applied *vertical* shift is ``t0/cth``. Every conversion
is an exact identity for a horizontal fibre. ``diameter`` is recomputed as
``(y_bot - y_top) * cth`` after the shifts, exactly as ``edges.detect_edges``
does -- ``qc`` trusts ``EdgeResult.diameter`` and never recomputes it.
``flags``/``amp``/``y_core``/``half_window`` are passed through untouched.

Pos
---
Fifth stage of the per-image pipeline, inserted between ``edges`` and ``qc``:
consumes ``D`` + the raw ``EdgeResult`` + ``BandResult``, refines (or passes
through unchanged) the boundaries that ``qc`` then cleans. Manual edits
(``manual_edit.py``) run after this stage and override it with no re-refine.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit
from scipy.special import erf

from .band import BandResult, tilt_geometry
from .config import CONFIG
from .edges import FLAG_OK, EdgeResult

_T_STEP = 0.5     # profile sampling step along t (perpendicular px)
_QUORUM = 0.7     # min fraction of a block's columns usable for a fit
_N_END = 6        # samples averaged at each profile end for the a/b guesses
_SIGMA0 = 3.0     # initial sigma guess (perpendicular px)
_SIGMA_FLOOR = 0.05  # hard optimiser bound on sigma (the gate is much tighter)
_SQRT2 = float(np.sqrt(2.0))


@dataclass
class RefineResult:
    refined_top: np.ndarray  # bool (W,) column's top boundary was actually refined
    refined_bot: np.ndarray  # bool (W,) column's bottom boundary was actually refined
    sigma_top: np.ndarray    # float32 (W,) fitted blur width, perpendicular px, NaN unrefined
    sigma_bot: np.ndarray
    resid_top: np.ndarray    # float32 (W,) relative fit residual, NaN unrefined
    resid_bot: np.ndarray
    o_top: np.ndarray        # float32 (W,) applied offset (+into fibre), perpendicular px, NaN unrefined
    o_bot: np.ndarray
    n_blocks: int            # blocks attempted (quorum met on at least one side --
    #                          the same set for both unless a wall's sample window
    #                          runs off the image on one side only)
    n_pass_top: int          # blocks whose top fit passed all gates
    n_pass_bot: int          # blocks whose bottom fit passed all gates

    @classmethod
    def empty(cls, W: int) -> RefineResult:
        """All-unrefined result for image width ``W`` (off-path / no-op placeholder)."""
        nan = np.full(W, np.nan, dtype=np.float32)
        return cls(
            refined_top=np.zeros(W, dtype=bool),
            refined_bot=np.zeros(W, dtype=bool),
            sigma_top=nan.copy(),
            sigma_bot=nan.copy(),
            resid_top=nan.copy(),
            resid_bot=nan.copy(),
            o_top=nan.copy(),
            o_bot=nan.copy(),
            n_blocks=0,
            n_pass_top=0,
            n_pass_bot=0,
        )


def _erf_model(t: np.ndarray, a: float, b: float, t0: float, sigma: float) -> np.ndarray:
    """Gaussian-blurred step: level ``a`` outside, ``b`` inside, midpoint ``t0``."""
    return a + (b - a) * 0.5 * (1.0 + erf((t - t0) / (sigma * _SQRT2)))


def _block_profile(
    D: np.ndarray, y_anchor: np.ndarray, cols: np.ndarray, t_grid: np.ndarray,
    sign: float, cth: float,
) -> np.ndarray | None:
    """Mean aligned wall profile over one block of columns (bilinear in y).

    Each column is sampled at ``y = y_anchor[x] + sign * t / cth`` (``sign``
    +1 for the top wall, -1 for the bottom, so ``+t`` always points into the
    fibre) and the block's columns are averaged -- per-column alignment, so
    fibre tilt never smears the ramp. One vectorised gather, the
    ``edges._vshift`` pattern. Columns without an anchor, and columns whose
    sample window leaves the image, are dropped *before* the ``_QUORUM``
    check; returns None when too few columns survive it.
    """
    H = D.shape[0]
    ya = y_anchor[cols]
    y_lo = ya + (sign / cth) * float(t_grid[0])
    y_hi = ya + (sign / cth) * float(t_grid[-1])
    keep = (
        np.isfinite(ya)
        & (np.minimum(y_lo, y_hi) >= 0.0)
        & (np.maximum(y_lo, y_hi) <= H - 1)
    )
    if keep.sum() < _QUORUM * cols.size:
        return None

    cx = cols[keep][None, :]
    ys = ya[keep][None, :] + (sign / cth) * t_grid[:, None]
    # every ys is inside [0, H-1] by the window test above, so clipping the
    # lower row to H-2 only moves the last row's weight onto f (no clamping)
    y0 = np.clip(np.floor(ys), 0, H - 2).astype(np.int32)
    f = (ys - y0).astype(np.float32)
    prof = (1.0 - f) * D[y0, cx] + f * D[y0 + 1, cx]
    return prof.mean(axis=1).astype(np.float32)


def _fit_block(
    t: np.ndarray, prof: np.ndarray, cfg: CONFIG
) -> tuple[float, float, float] | None:
    """Least-squares erf fit of one block profile -> ``(t0, sigma, resid)``.

    ``resid`` is the rms residual relative to the step height ``b - a``.
    Returns None when the optimiser fails or any gate rejects the fit: the
    step must rise into the fibre (``b - a > 0``), be explained to within
    ``refine_relmax``, have ``sigma`` inside ``[refine_sigma_min,
    refine_sigma_max]`` and sit within ``refine_maxshift`` of the legacy edge.
    """
    a0 = float(prof[:_N_END].mean())
    b0 = float(prof[-_N_END:].mean())
    t0_0 = float(t[int(np.argmin(np.abs(prof - 0.5 * (a0 + b0))))])
    p0 = [a0, b0, t0_0, _SIGMA0]
    bounds = (
        [-np.inf, -np.inf, float(t[0]), _SIGMA_FLOOR],
        [np.inf, np.inf, float(t[-1]), np.inf],
    )
    try:
        with warnings.catch_warnings():
            # a rank-deficient covariance is not an error here: the gates below
            # judge the fit, and the warning would pollute every caller's output
            warnings.simplefilter("ignore", OptimizeWarning)
            popt, _ = curve_fit(_erf_model, t, prof, p0=p0, bounds=bounds, maxfev=2000)
    except (RuntimeError, ValueError, TypeError):
        return None

    a, b, t0, sigma = (float(v) for v in popt)
    step = b - a
    if step <= 0.0:
        return None
    resid = float(np.sqrt(np.mean((prof - _erf_model(t, *popt)) ** 2)) / step)
    if not resid < cfg.refine_relmax:
        return None
    if not cfg.refine_sigma_min <= sigma <= cfg.refine_sigma_max:
        return None
    if not abs(t0) <= cfg.refine_maxshift:
        return None
    return t0, sigma, resid


def _interp_side(
    W: int, blk_idx: np.ndarray, centres: np.ndarray, vals: np.ndarray, cfg: CONFIG
) -> np.ndarray:
    """Spread the passing blocks' values over the columns between their centres.

    ``vals`` is ``(K, N)`` -- K quantities over the N passing blocks, in block
    order. Consecutive passing centres separated by more than
    ``refine_gap_blocks`` failed blocks start a new chain, and nothing is
    extrapolated beyond a chain's outermost centres, so the returned ``(K, W)``
    array is NaN wherever no fit vouches for the column.
    """
    out = np.full((vals.shape[0], W), np.nan, dtype=np.float32)
    if centres.size == 0:
        return out
    breaks = np.where(np.diff(blk_idx) - 1 > cfg.refine_gap_blocks)[0] + 1
    for chain in np.split(np.arange(centres.size), breaks):
        xs = np.arange(int(centres[chain[0]]), int(centres[chain[-1]]) + 1)
        for k in range(vals.shape[0]):
            out[k, xs] = np.interp(xs, centres[chain], vals[k, chain])
    return out


def _masked(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """``values`` where ``mask``, NaN elsewhere (float32 diagnostics array)."""
    out = np.full(values.size, np.nan, dtype=np.float32)
    out[mask] = values[mask]
    return out


def refine_edges(
    D: np.ndarray, edg: EdgeResult, bnd: BandResult, cfg: CONFIG
) -> tuple[EdgeResult, RefineResult]:
    """Refit each wall as a blurred step and shift the boundary to its midpoint.

    Returns ``(edg, RefineResult.empty(W))`` with ``edg`` the SAME object
    (identity, not a copy) when ``cfg.refine_on`` is False -- bit-identical
    off-path -- and likewise when there is nothing to fit. Otherwise both walls
    are fitted block by block, the accepted midpoint offsets are interpolated
    between passing block centres and applied to the anchor columns they cover,
    and a fresh ``EdgeResult`` (shifted ``y_top``/``y_bot``, recomputed
    ``diameter``) is returned with the refinement diagnostics.
    """
    W = D.shape[1]
    if not cfg.refine_on:
        return edg, RefineResult.empty(W)

    _m, cth = tilt_geometry(bnd.slope)
    anchor = np.isfinite(edg.y_top) & np.isfinite(edg.y_bot) & (edg.flags == FLAG_OK)

    blk = max(1, int(cfg.refine_block))
    x0 = int(bnd.x0)
    n_blk = max(0, int(bnd.x1) - x0 + 1) // blk  # trailing partial block skipped
    if n_blk == 0 or not anchor.any():
        return edg, RefineResult.empty(W)

    # profile window, perpendicular px: fixed reach outward, band-scaled inward
    # (capped so the specular rim inside the fibre stays out of the fit)
    inside = float(np.clip(cfg.refine_in_frac * bnd.band_half * cth, 8.0, cfg.refine_in_max))
    n_t = int(round((float(cfg.refine_out) + inside) / _T_STEP)) + 1
    t_grid = (-float(cfg.refine_out) + _T_STEP * np.arange(n_t)).astype(np.float32)

    attempted = np.zeros(n_blk, dtype=bool)
    fields: dict[str, np.ndarray] = {}
    n_pass: dict[str, int] = {}
    for name, y_edge, sign in (("top", edg.y_top, 1.0), ("bot", edg.y_bot, -1.0)):
        ya = np.where(anchor, y_edge, np.nan).astype(np.float32)
        idx: list[int] = []
        centres: list[int] = []
        vals: list[tuple[float, float, float]] = []
        for k in range(n_blk):
            cols = np.arange(x0 + k * blk, x0 + (k + 1) * blk)
            prof = _block_profile(D, ya, cols, t_grid, sign, cth)
            if prof is None:
                continue
            attempted[k] = True
            fit = _fit_block(t_grid, prof, cfg)
            if fit is None:
                continue
            idx.append(k)
            centres.append(x0 + k * blk + blk // 2)
            vals.append(fit)
        n_pass[name] = len(idx)
        fields[name] = _interp_side(
            W,
            np.asarray(idx, dtype=np.int64),
            np.asarray(centres, dtype=np.int64),
            np.asarray(vals, dtype=np.float64).reshape(len(vals), 3).T,
            cfg,
        )

    o_top, sig_top, res_top = fields["top"]
    o_bot, sig_bot, res_bot = fields["bot"]
    ref_top = anchor & np.isfinite(o_top)
    ref_bot = anchor & np.isfinite(o_bot)

    # +t0 points into the fibre on both sides: the top edge moves DOWN, the
    # bottom edge UP. t0 is perpendicular, the shift applied is vertical.
    y_top = edg.y_top.copy()
    y_bot = edg.y_bot.copy()
    y_top[ref_top] += o_top[ref_top] / cth
    y_bot[ref_bot] -= o_bot[ref_bot] / cth
    # qc.run_qc reads EdgeResult.diameter as-is -- recompute it here or the
    # shifts never reach the measurement
    diameter = (y_bot - y_top) * cth

    ref = RefineResult(
        refined_top=ref_top,
        refined_bot=ref_bot,
        sigma_top=_masked(sig_top, ref_top),
        sigma_bot=_masked(sig_bot, ref_bot),
        resid_top=_masked(res_top, ref_top),
        resid_bot=_masked(res_bot, ref_bot),
        o_top=_masked(o_top, ref_top),
        o_bot=_masked(o_bot, ref_bot),
        n_blocks=int(attempted.sum()),
        n_pass_top=n_pass["top"],
        n_pass_bot=n_pass["bot"],
    )
    return replace(edg, y_top=y_top, y_bot=y_bot, diameter=diameter), ref
