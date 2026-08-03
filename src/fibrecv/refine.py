"""Erf edge refinement: refit each detected wall as a blurred step, shift to its midpoint.

Dependencies
------------
``numpy`` (``scipy.optimize.curve_fit`` joins from Task 2 onward, when the
fitting itself lands). Consumes ``band`` (``BandResult``, ``tilt_geometry``)
and ``edges`` (``EdgeResult``).

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

NOTE (M1 skeleton): this module currently ships the dataclass and wiring
only. ``refine_edges``'s on-path is a no-op placeholder -- it returns the
input ``EdgeResult`` unchanged and an empty ``RefineResult`` (``n_blocks=0``)
exactly like the off-path -- so nothing above is implemented yet. The fitting
algorithm itself lands in Task 2.

Pos
---
Fifth stage of the per-image pipeline, inserted between ``edges`` and ``qc``:
consumes ``D`` + the raw ``EdgeResult`` + ``BandResult``, refines (or passes
through unchanged) the boundaries that ``qc`` then cleans. Manual edits
(``manual_edit.py``) run after this stage and override it with no re-refine.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .band import BandResult
from .config import CONFIG
from .edges import EdgeResult


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
    n_blocks: int            # blocks attempted (quorum met; same count both sides)
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


def refine_edges(
    D: np.ndarray, edg: EdgeResult, bnd: BandResult, cfg: CONFIG
) -> tuple[EdgeResult, RefineResult]:
    """Refit each wall as a blurred step and shift the boundary to its midpoint.

    Returns ``(edg, RefineResult.empty(W))`` with ``edg`` the SAME object
    (identity, not a copy) when ``cfg.refine_on`` is False -- bit-identical
    off-path. In this M1 skeleton the on-path is also a no-op placeholder,
    behaving exactly like the off-path; the real batched erf-model fit over
    ``refine_block``-column blocks lands in Task 2.
    """
    W = D.shape[1]
    if not cfg.refine_on:
        return edg, RefineResult.empty(W)
    # TODO(Task 2): batched erf-model fit + gate + interpolate per refine_block
    # columns; for now, no-op (identical to the off-path).
    return edg, RefineResult.empty(W)
