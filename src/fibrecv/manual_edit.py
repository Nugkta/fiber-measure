"""Pure manual boundary correction: anchors/nudges -> corrected edges + re-QC.

Dependencies
------------
``numpy`` plus the pipeline modules ``qc`` (re-run after an edit) and ``edges``
(``FLAG_OK``). Deliberately imports no Streamlit so the logic is headlessly
testable; the GUI is a thin caller.

Inputs
------
- A (possibly cached) ``MeasureResult`` and a per-image edits dict
  ``{"top": [set, ...], "bot": [...], "nudge_top": float, "nudge_bot": float}``
  where each set is an independent list of ``(x, y)`` anchors in native
  full-image pixel coordinates (sub-pixel floats are fine). Sets correct
  disjoint stretches independently — the detected line between two sets is
  never redrawn; where sets overlap, the later one wins.
- ``CONFIG`` for the QC re-run and the µm conversion.

Output
------
- ``apply_manual_edits`` -> a NEW ``MeasureResult`` (the input is never mutated,
  so ``st.cache_data`` entries stay pristine) plus per-column edited masks.
- ``corrected_boundary`` / ``display_to_native`` / ``has_edits`` /
  ``empty_edits``: the small pure pieces the GUI composes.

Pos
---
Sits between the cached compute results and every downstream consumer in the
GUI (overlay, profile plot, registration, export): ``gui_app`` applies the
edits once, right after loading, so corrections flow everywhere with no further
plumbing. Batch/CLI never see these edits by construction.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import replace
from typing import TYPE_CHECKING, Sequence

import numpy as np

from . import qc as qc_mod
from .config import CONFIG
from .edges import FLAG_OK

if TYPE_CHECKING:  # avoid importing compute at runtime (only needed for types)
    from .compute import MeasureResult

RAMP_PX = 25  # blend length at each end of an edited range (px)


def display_to_native(
    x_disp: float,
    y_disp: float,
    disp_w: int,
    disp_h: int,
    crop_w: int,
    crop_h: int,
    y_offset: int,
    img_w: int,
    img_h: int,
) -> tuple[float, float]:
    """Map a click on the displayed crop to native full-image pixel coords.

    The crop may be stretched independently in x and y (the GUI stretches
    vertically for precision), so the two scale factors are separate. The
    result is clamped to the native image bounds.
    """
    sx = crop_w / disp_w
    sy = crop_h / disp_h
    x_nat = float(np.clip(x_disp * sx, 0.0, img_w - 1))
    y_nat = float(np.clip(y_disp * sy + y_offset, 0.0, img_h - 1))
    return x_nat, y_nat


def corrected_boundary(
    y: np.ndarray,
    anchors: Sequence[tuple[float, float]],
    ramp: int = RAMP_PX,
) -> tuple[np.ndarray, np.ndarray]:
    """Redraw a boundary through user anchors; return (y_corrected, edited).

    The corrected line is the linear interpolation through the anchors over
    ``[first_anchor_x, last_anchor_x]`` only (never extrapolated), blended into
    the detected line over ``ramp`` px at each end so no step is introduced.
    Where the detected line is NaN the correction replaces it outright. A
    single anchor edits exactly its own column. Pure: ``y`` is not modified.
    """
    W = y.size
    y_new = y.copy()
    edited = np.zeros(W, dtype=bool)
    if not anchors:
        return y_new, edited

    # clamp x into the image; collapse duplicate integer columns keeping the
    # last-clicked anchor (click order = user intent)
    by_col: dict[int, tuple[float, float]] = {}
    for ax, ay in anchors:
        axc = float(np.clip(float(ax), 0.0, W - 1.0))
        by_col[int(round(axc))] = (axc, float(ay))
    pts = sorted(by_col.values())
    xs = np.array([p[0] for p in pts])
    ys = np.array([p[1] for p in pts])

    xa, xb = int(round(xs[0])), int(round(xs[-1]))
    cols = np.arange(xa, xb + 1)
    target = np.interp(cols, xs, ys)

    # end-ramp blend, entirely inside the edited range
    half = min(ramp, (cols.size + 1) // 2)
    w = np.minimum(1.0, np.minimum((cols - xa + 1) / half, (xb - cols + 1) / half))
    det = y[cols]
    finite = np.isfinite(det)
    w = np.where(finite, w, 1.0)        # nothing to blend against under NaN
    base = np.where(finite, det, target)
    y_new[cols] = w * target + (1.0 - w) * base
    edited[cols] = True
    return y_new, edited


def corrected_segments(
    y: np.ndarray,
    segments: Sequence[Sequence[tuple[float, float]]],
    ramp: int = RAMP_PX,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply independent anchor sets to a boundary; return (y_new, edited).

    Each set redraws only its own ``[first_anchor_x, last_anchor_x]`` range via
    ``corrected_boundary``; the line between sets is untouched. Sets apply in
    list order, each blending against the result of the previous ones, so a
    later overlapping set wins. Empty sets are skipped. Pure.
    """
    y_new = y.copy()
    edited = np.zeros(y.size, dtype=bool)
    for anchors in segments:
        if not anchors:
            continue
        y_new, e = corrected_boundary(y_new, anchors, ramp)
        edited |= e
    return y_new, edited


def empty_edits() -> dict:
    """A fresh per-image edits dict (the GUI's session-state value).

    ``top``/``bot`` hold lists of independent anchor sets (each a list of
    ``(x, y)`` points).
    """
    return {"top": [], "bot": [], "nudge_top": 0.0, "nudge_bot": 0.0}


def has_edits(edits: dict | None) -> bool:
    """True if the edits dict contains any anchor set with points or a nudge."""
    if not edits:
        return False
    return bool(
        any(edits.get("top", []))
        or any(edits.get("bot", []))
        or edits.get("nudge_top", 0.0)
        or edits.get("nudge_bot", 0.0)
    )


def apply_manual_edits(
    mr: "MeasureResult",
    edits: dict,
    cfg: CONFIG,
    ramp: int = RAMP_PX,
) -> tuple["MeasureResult", np.ndarray, np.ndarray]:
    """Apply manual edits to a (possibly cached) ``MeasureResult``.

    Nudges move each whole line first, then each anchor set independently
    redraws its own range (``corrected_segments``). Edited columns where both
    boundaries are finite and the diameter is positive get their flags cleared
    to ``FLAG_OK`` (user override), then QC is re-run in full so smoothing,
    coverage and confidence flags refresh.

    Returns ``(new_mr, edited_top, edited_bot)``. ``mr`` is never mutated
    (cache-safe): new ``EdgeResult``/``QCResult``/``meta`` are built via
    ``dataclasses.replace`` and copies.
    """
    y_top = mr.edg.y_top.copy()
    y_bot = mr.edg.y_bot.copy()
    flags = mr.edg.flags.copy()
    W = y_top.size
    edited_top = np.zeros(W, dtype=bool)
    edited_bot = np.zeros(W, dtype=bool)

    nudge_top = float(edits.get("nudge_top", 0.0) or 0.0)
    nudge_bot = float(edits.get("nudge_bot", 0.0) or 0.0)
    if nudge_top != 0.0:
        m = np.isfinite(y_top)
        y_top[m] += nudge_top
        edited_top |= m
    if nudge_bot != 0.0:
        m = np.isfinite(y_bot)
        y_bot[m] += nudge_bot
        edited_bot |= m

    y_top, et = corrected_segments(y_top, edits.get("top", []), ramp)
    edited_top |= et
    y_bot, eb = corrected_segments(y_bot, edits.get("bot", []), ramp)
    edited_bot |= eb

    diameter = y_bot - y_top
    ok = (
        (edited_top | edited_bot)
        & np.isfinite(y_top)
        & np.isfinite(y_bot)
        & (diameter > 0)
    )
    flags[ok] = FLAG_OK  # user override; one-sided NaN columns keep their flag

    new_edg = replace(mr.edg, y_top=y_top, y_bot=y_bot, diameter=diameter,
                      flags=flags)
    new_res = qc_mod.run_qc(new_edg, mr.bnd, cfg)
    diameter_um = np.where(new_res.valid, new_res.diameter_raw / cfg.ppu, np.nan)

    # refresh the meta fields derived from edges/QC (same formulas as
    # compute.compute_measurement) + JSON-safe provenance of the edit
    span = slice(mr.bnd.x0, mr.bnd.x1 + 1)
    flag_counts = Counter(int(f) for f in new_res.reason[span])
    meta = dict(mr.meta)
    meta.update({
        "coverage": new_res.coverage,
        "n_valid": int(new_res.valid.sum()),
        "low_confidence": bool(new_res.low_confidence),
        "band_mismatch": bool(new_res.band_mismatch),
        "flag_counts": {str(k): int(v) for k, v in flag_counts.items()},
        "median_diameter_um": (float(np.nanmedian(diameter_um))
                               if new_res.valid.any() else None),
        "manual_edit": {
            "top_sets": [[[float(x), float(y)] for x, y in s]
                         for s in edits.get("top", []) if s],
            "bot_sets": [[[float(x), float(y)] for x, y in s]
                         for s in edits.get("bot", []) if s],
            "nudge_top": nudge_top,
            "nudge_bot": nudge_bot,
        },
    })

    new_mr = replace(mr, edg=new_edg, res=new_res, diameter_um=diameter_um,
                     meta=meta)
    return new_mr, edited_top, edited_bot
