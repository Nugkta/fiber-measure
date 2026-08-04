"""Anomaly flagging: image-level profile defects and group-level outliers.

Dependencies
------------
``numpy`` and ``.config`` only -- pure functions, no I/O, no pipeline imports
(``qc`` imports this module, never the reverse).

Inputs
------
- Per-image: edge rows, smoothed diameter, validity mask, span and centerline
  slope from ``EdgeResult``/``QCResult``/``BandResult`` (passed as plain arrays).
- Group-level: per-replicate median diameters.
- ``CONFIG`` thresholds (``jump_thresh_px``, ``gap_frac``, ``step_frac``,
  ``step_window_px``, ``rep_dev_frac``, ``anomaly_exclude``).

Output
------
``AnomalyResult`` (flags + evidence scalars, JSON-safe via ``as_dict``),
``detect_image_anomalies``, ``detect_replicate_outliers``, and the shared
``exclusion_reason`` policy used by both the CLI aggregator and the GUI.

Pos
---
Side-car of the QC stage: ``run_qc`` attaches an ``AnomalyResult`` to every
image; ``run_aggregate`` and the GUI add the group-level ``replicate_outlier``
flag and apply ``exclusion_reason``. Advisory by default -- flags only exclude
a replicate when ``anomaly_exclude`` is on, and ``replicate_outlier`` never
excludes (diameter genuinely varies along a fibre).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .config import CONFIG

# image-level flag names (replicate_outlier is added at aggregation time)
FLAG_EDGE_JUMP = "edge_jump"
FLAG_LARGE_GAP = "large_gap"
FLAG_DIAMETER_STEP = "diameter_step"
FLAG_REPLICATE_OUTLIER = "replicate_outlier"

_JUMP_COLS_CAP = 20  # as_dict keeps at most this many example columns


@dataclass
class AnomalyResult:
    """Per-image anomaly flags plus the evidence behind each one."""

    flags: list[str] = field(default_factory=list)
    jump_cols: list[int] = field(default_factory=list)  # columns after a jump
    max_jump_px: float = 0.0        # largest detrended edge displacement seen
    n_jumps: int = 0                # total flagged pairs across both edges
    longest_gap_frac: float = 0.0   # longest invalid run / span length
    gap_start_col: int | None = None
    step_frac: float = float("nan")  # max window-median shift / global median
    step_col: int | None = None      # centre of the maximal-shift plateau

    def as_dict(self) -> dict:
        """JSON-safe snapshot (NaN -> None, numpy scalars -> python)."""
        step = self.step_frac
        step = None if step is None or not np.isfinite(step) else float(step)
        return {
            "flags": [str(f) for f in self.flags],
            "jump_cols": [int(c) for c in self.jump_cols[:_JUMP_COLS_CAP]],
            "max_jump_px": float(self.max_jump_px),
            "n_jumps": int(self.n_jumps),
            "longest_gap_frac": float(self.longest_gap_frac),
            "gap_start_col": None if self.gap_start_col is None else int(self.gap_start_col),
            "step_frac": step,
            "step_col": None if self.step_col is None else int(self.step_col),
        }


def _edge_jumps(y: np.ndarray, valid: np.ndarray, slope: float,
                cfg: CONFIG) -> tuple[list[int], float, int]:
    """Detrended jumps between consecutive *valid* columns of one edge.

    Diffing valid neighbours (not adjacent columns) means a gap does not hide
    a jump; subtracting ``slope * dx`` keeps a tilted fibre from false-flagging
    across that gap.
    """
    idx = np.where(valid & np.isfinite(y))[0]
    if idx.size < 2:
        return [], 0.0, 0
    dy = np.diff(y[idx])
    dx = np.diff(idx).astype(float)
    resid = np.abs(dy - slope * dx)
    bad = np.where(resid > cfg.jump_thresh_px)[0]
    # flag the later column of each offending pair
    return [int(idx[i + 1]) for i in bad], float(resid.max()), int(bad.size)


def _longest_invalid_run(valid: np.ndarray, x0: int, x1: int) -> tuple[int, int | None]:
    """(length, start column) of the longest invalid run inside the span."""
    inv = ~valid[x0:x1 + 1]
    if not inv.any():
        return 0, None
    edges = np.diff(inv.astype(np.int8))
    starts = np.where(edges == 1)[0] + 1
    ends = np.where(edges == -1)[0] + 1
    if inv[0]:
        starts = np.r_[0, starts]
    if inv[-1]:
        ends = np.r_[ends, inv.size]
    lengths = ends - starts
    k = int(np.argmax(lengths))
    return int(lengths[k]), x0 + int(starts[k])


def _diameter_step(diameter_smooth: np.ndarray, x0: int, x1: int,
                   cfg: CONFIG) -> tuple[float, int | None]:
    """Max adjacent-window median shift as a fraction of the global median.

    Two ``step_window_px`` windows slide along the span; a genuine level shift
    produces a plateau of maximal shift, so the reported column is the plateau
    centre (the breakpoint). Returns (NaN, None) when the span is shorter than
    two windows or the profile has no finite global median.
    """
    w = int(cfg.step_window_px)
    span_len = x1 - x0 + 1
    if w < 1 or span_len < 2 * w:
        return float("nan"), None
    seg = np.asarray(diameter_smooth[x0:x1 + 1], dtype=float)
    fin = seg[np.isfinite(seg)]
    if fin.size == 0:
        return float("nan"), None
    g_med = float(np.median(fin))
    if g_med <= 0:
        return float("nan"), None
    stride = max(1, w // 25)
    best = -1.0
    plateau: list[int] = []
    for c in range(w, span_len - w + 1, stride):
        left = seg[c - w:c]
        right = seg[c:c + w]
        left = left[np.isfinite(left)]
        right = right[np.isfinite(right)]
        if left.size == 0 or right.size == 0:
            continue
        shift = abs(float(np.median(right)) - float(np.median(left))) / g_med
        if shift > best + 1e-12:
            best = shift
            plateau = [c]
        elif shift > best - 1e-12:
            plateau.append(c)
    if best < 0:
        return float("nan"), None
    return float(best), x0 + plateau[len(plateau) // 2]


def detect_image_anomalies(
    y_top: np.ndarray,
    y_bot: np.ndarray,
    diameter_smooth: np.ndarray,
    valid: np.ndarray,
    x0: int,
    x1: int,
    slope: float,
    cfg: CONFIG,
) -> AnomalyResult:
    """Run the three per-image detectors (edge_jump, large_gap, diameter_step)."""
    res = AnomalyResult()
    span_len = max(1, x1 - x0 + 1)

    cols: set[int] = set()
    for y in (y_top, y_bot):
        edge_cols, edge_max, edge_n = _edge_jumps(y, valid, slope, cfg)
        cols.update(edge_cols)
        res.max_jump_px = max(res.max_jump_px, edge_max)
        res.n_jumps += edge_n
    res.jump_cols = sorted(cols)
    if res.n_jumps > 0:
        res.flags.append(FLAG_EDGE_JUMP)

    gap_len, res.gap_start_col = _longest_invalid_run(valid, x0, x1)
    res.longest_gap_frac = gap_len / span_len
    if res.longest_gap_frac > cfg.gap_frac:
        res.flags.append(FLAG_LARGE_GAP)

    res.step_frac, res.step_col = _diameter_step(diameter_smooth, x0, x1, cfg)
    if np.isfinite(res.step_frac) and res.step_frac > cfg.step_frac:
        res.flags.append(FLAG_DIAMETER_STEP)

    return res


def detect_replicate_outliers(
    medians: dict, cfg: CONFIG
) -> tuple[dict, set]:
    """Replicates whose median diameter deviates > rep_dev_frac from the group.

    ``medians`` maps any hashable per-image key (image name in the CLI, rep
    index in the GUI) -> median diameter; None/NaN entries are dropped. Needs
    >= 3 finite medians -- with fewer, the group median is not meaningful and
    everything passes. Returns ({key: deviation_frac}, {outlier keys}).
    Advisory only: callers must never exclude on this flag.
    """
    finite = {k: float(v) for k, v in medians.items()
              if v is not None and np.isfinite(v)}
    if len(finite) < 3:
        return {}, set()
    g_med = float(np.median(list(finite.values())))
    if g_med <= 0:
        return {}, set()
    devs = {k: abs(v - g_med) / g_med for k, v in finite.items()}
    return devs, {k for k, d in devs.items() if d > cfg.rep_dev_frac}


def exclusion_reason(
    band_mismatch: bool,
    coverage: float | None,
    anomaly_flags: list[str] | None,
    cfg: CONFIG,
) -> str | None:
    """Why a replicate is dropped from registration, or None to keep it.

    Single policy shared by run_aggregate and the GUI. Priority:
    band_mismatch, then coverage, then image-level anomalies (only when
    ``anomaly_exclude`` is on; ``replicate_outlier`` never excludes).
    """
    if band_mismatch:
        return "band_mismatch"
    if coverage is not None and coverage < cfg.min_coverage:
        return f"coverage {coverage:.0%} < {cfg.min_coverage:.0%}"
    if cfg.anomaly_exclude:
        image_flags = [f for f in (anomaly_flags or [])
                       if f != FLAG_REPLICATE_OUTLIER]
        if image_flags:
            return "anomaly: " + ", ".join(image_flags)
    return None
