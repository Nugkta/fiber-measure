"""Quality control: NaN policy, outlier rejection, smoothing, coverage.

Dependencies
------------
``numpy``, ``scipy.ndimage.median_filter``, ``scipy.signal.savgol_filter``.

Inputs
------
- ``EdgeResult`` (raw per-column diameters, edge rows, flags).
- ``BandResult`` (span; its centerline is only a fallback -- the deviation check
  refits a robust line to the measured centers, since the coarse-mask centroid
  can be biased toward a bright reflection rim).
- ``CONFIG`` for the rejection thresholds, smoothing kernels and coverage gate.

Output
------
``QCResult`` with ``diameter_raw`` (NaN where rejected), ``diameter_smooth``
(median + Savitzky-Golay on a gap-filled series), ``valid`` and ``interpolated``
boolean masks, per-column ``reason`` flags, and scalar ``coverage`` /
``low_confidence``.

Pos
---
Fifth stage of the per-image pipeline. Turns noisy raw edges into a clean,
NaN-honest profile + a smoothed view, and decides per-image confidence. Feeds
the per-image CSV and the registration stage.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import savgol_filter
from scipy.stats import theilslopes

from .band import BandResult
from .config import CONFIG
from .edges import EdgeResult, FLAG_OK


@dataclass
class QCResult:
    diameter_raw: np.ndarray     # float (W,) accepted raw diameters, NaN elsewhere
    diameter_smooth: np.ndarray  # float (W,) denoised, NaN outside the band span
    valid: np.ndarray            # bool (W,) column passed all QC
    interpolated: np.ndarray     # bool (W,) value was gap-filled for smoothing
    reason: np.ndarray           # int (W,) accumulated rejection flags
    coverage: float              # fraction of span columns that are valid
    low_confidence: bool         # coverage/band/band-mismatch problem
    band_mismatch: bool          # diameter << coarse band thickness (defocus trap)


# extra QC reasons (kept distinct from edge flags' bit space, reused additively)
FLAG_CENTER_DEV = 32
FLAG_ROLL_OUTLIER = 64


def _rolling_mad_outliers(d: np.ndarray, valid: np.ndarray, cfg: CONFIG) -> np.ndarray:
    """Boolean mask of columns whose diameter is a local rolling-MAD outlier."""
    out = np.zeros_like(valid)
    n = d.size
    if valid.sum() < 5:
        return out
    half = max(1, cfg.roll_window // 2)
    idx = np.where(valid)[0]
    for i in idx:
        lo, hi = max(0, i - half), min(n, i + half + 1)
        window = d[lo:hi]
        window = window[np.isfinite(window)]
        if window.size < 5:
            continue
        med = np.median(window)
        mad = np.median(np.abs(window - med))
        if mad <= cfg.eps:
            continue
        if abs(d[i] - med) > cfg.roll_k * cfg.mad_scale * mad:
            out[i] = True
    return out


def run_qc(edges: EdgeResult, band: BandResult, cfg: CONFIG) -> QCResult:
    """Apply the NaN policy, reject outliers, smooth, and score coverage."""
    W = edges.diameter.size
    diam = edges.diameter.astype(np.float32).copy()
    reason = edges.flags.copy()

    span = np.zeros(W, dtype=bool)
    span[band.x0:band.x1 + 1] = True

    # start from edge flags: any non-OK flag invalidates the column
    valid = span & np.isfinite(diam) & (edges.flags == FLAG_OK)

    # centerline deviation: measured center vs a robust line fit of the
    # *measured* centers themselves (the coarse-mask centroid can be biased
    # toward a bright reflection rim, so band.c_fit is only used for windowing)
    center = (edges.y_top + edges.y_bot) / 2.0
    fin = np.where(np.isfinite(center) & valid)[0]
    if fin.size >= 10:
        slope, intercept, _, _ = theilslopes(center[fin], fin)
        line = slope * np.arange(W) + intercept
    else:
        line = band.c_fit
    dev = np.abs(center - line)
    center_bad = np.isfinite(dev) & (dev > cfg.reject_dev)
    reason[center_bad] |= FLAG_CENTER_DEV
    valid &= ~center_bad

    # non-positive / absurd diameters
    valid &= np.isfinite(diam) & (diam > 0)

    # rolling-MAD outlier rejection (only among currently-valid columns)
    roll_bad = _rolling_mad_outliers(diam, valid, cfg)
    reason[roll_bad] |= FLAG_ROLL_OUTLIER
    valid &= ~roll_bad

    diameter_raw = np.where(valid, diam, np.nan).astype(np.float32)

    # smoothed series: linear-interpolate gaps within the span, then median+savgol
    diameter_smooth = np.full(W, np.nan, dtype=np.float32)
    interpolated = np.zeros(W, dtype=bool)
    xs = np.where(valid)[0]
    if xs.size >= 2:
        x_span = np.arange(band.x0, band.x1 + 1)
        filled = np.interp(x_span, xs, diameter_raw[xs])
        interp_mask = ~np.isin(x_span, xs)
        interpolated[x_span[interp_mask]] = True
        # median filter (kernel clamped to odd <= series length)
        k = min(cfg.median_k, filled.size if filled.size % 2 == 1 else filled.size - 1)
        k = max(1, k if k % 2 == 1 else k - 1)
        med = median_filter(filled, size=k, mode="nearest")
        # Savitzky-Golay (window clamped to odd <= series length)
        win = min(cfg.savgol_window, med.size if med.size % 2 == 1 else med.size - 1)
        win = max(cfg.savgol_poly + 2 - (cfg.savgol_poly % 2), win if win % 2 == 1 else win - 1)
        if win >= cfg.savgol_poly + 2 and win <= med.size:
            sg = savgol_filter(med, window_length=win, polyorder=cfg.savgol_poly)
        else:
            sg = med
        diameter_smooth[x_span] = sg.astype(np.float32)

    n_span = int(span.sum())
    coverage = float(valid.sum()) / n_span if n_span > 0 else 0.0

    # band-consistency check: a measured diameter far thinner than the coarse
    # desaturation band means the detector likely locked onto a sharp internal
    # feature inside a defocused blur band; far WIDER means it likely grabbed
    # a shadow/halo outside the fibre -> confident garbage either way, flag it
    band_thickness = 2.0 * band.band_half
    med_diam = float(np.nanmedian(diameter_raw)) if valid.any() else np.nan
    band_mismatch = bool(
        np.isfinite(med_diam)
        and band_thickness > 0
        and (med_diam < cfg.band_ratio_min * band_thickness
             or med_diam > cfg.band_ratio_max * band_thickness)
    )

    low_conf = bool(
        band.low_confidence or coverage < cfg.min_coverage or band_mismatch
    )

    return QCResult(
        diameter_raw=diameter_raw,
        diameter_smooth=diameter_smooth,
        valid=valid,
        interpolated=interpolated,
        reason=reason,
        coverage=coverage,
        low_confidence=low_conf,
        band_mismatch=band_mismatch,
    )
