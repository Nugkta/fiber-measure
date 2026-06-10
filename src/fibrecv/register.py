"""Replicate registration + pointwise mean/variance statistics.

Dependencies
------------
``numpy``, ``scipy.signal.correlate``.

Inputs
------
- A list of replicate profiles for one ``A_B`` sample, each a dict with keys
  ``x`` (column index), ``diameter_px_raw``, ``diameter_px_smooth``, ``valid``,
  ``replicate``.
- ``CONFIG`` for ``ppu``, ``max_shift`` and ``min_corr``.

Output
------
- ``register_sample(profiles, cfg)`` -> ``(table, shifts, summary)`` where
  ``table`` is the aligned pointwise mean+/-std curve in microns, ``shifts`` are
  per-replicate lag/corr/flags, and ``summary`` is the scalar row for
  master_summary.csv.

Method
------
Each replicate's smoothed profile is cross-correlated against replicate ``_1``
(zero-meaned, NaN-filled), integer lag refined by 3-point parabolic interpolation
and bounded by ``max_shift``. Replicates are shifted onto a common integer grid;
pointwise statistics use the **raw** (NaN-honest) diameters so the variance is
real. Weak correlation peaks fall back to zero shift + ``registration_uncertain``.

Pos
---
Aggregation stage (run by ``run_aggregate.py``) -- the only stage that combines
the three same-fibre replicate photographs into one mean+/-variance curve.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import correlate

from .config import CONFIG


def _fill_for_corr(diam: np.ndarray) -> np.ndarray:
    """Zero-mean, NaN-filled copy of a profile suitable for cross-correlation."""
    d = diam.astype(np.float64)
    finite = np.isfinite(d)
    if finite.sum() < 2:
        return np.zeros_like(d)
    xs = np.where(finite)[0]
    d = np.interp(np.arange(d.size), xs, d[xs])
    d = d - np.mean(d)
    return d


def _parabolic_peak(corr: np.ndarray, k: int) -> float:
    """3-point parabolic sub-pixel offset around integer peak index ``k``."""
    if k <= 0 or k >= corr.size - 1:
        return 0.0
    ym1, y0, yp1 = corr[k - 1], corr[k], corr[k + 1]
    denom = ym1 - 2 * y0 + yp1
    if abs(denom) < 1e-12:
        return 0.0
    return 0.5 * (ym1 - yp1) / denom


def estimate_shift(ref: np.ndarray, other: np.ndarray, cfg: CONFIG) -> tuple[float, float, bool]:
    """Sub-pixel lag of ``other`` relative to ``ref``; returns (shift, corr, uncertain).

    ``shift`` > 0 means ``other`` must move right (+x) to align with ``ref``.
    """
    a = _fill_for_corr(ref)
    b = _fill_for_corr(other)
    if np.allclose(a, 0) or np.allclose(b, 0):
        return 0.0, 0.0, True
    corr = correlate(a, b, mode="full")
    lags = np.arange(-(b.size - 1), a.size)
    # bound the search to +/- max_shift
    keep = np.abs(lags) <= cfg.max_shift
    corr_b = corr[keep]
    lags_b = lags[keep]
    k = int(np.argmax(corr_b))
    lag = int(lags_b[k])
    sub = _parabolic_peak(corr_b, k)
    shift = lag + sub
    # normalised peak strength (Pearson-like at best lag)
    norm = np.sqrt(np.sum(a * a) * np.sum(b * b)) + 1e-12
    peak = float(corr_b[k] / norm)
    uncertain = peak < cfg.min_corr
    if uncertain:
        shift = 0.0
    return float(shift), peak, bool(uncertain)


def register_sample(profiles: list[dict], cfg: CONFIG) -> tuple[dict, list[dict], dict]:
    """Align replicates and compute pointwise + scalar statistics.

    Returns ``(table, shifts, summary)``. ``profiles`` must be non-empty; the
    first (lowest replicate number) is the reference.
    """
    profiles = sorted(profiles, key=lambda p: p["replicate"])
    ref = profiles[0]
    ref_x = ref["x"].astype(np.float64)

    shifts = []
    aligned = []  # (shifted_x, raw_diam_px) per replicate
    for p in profiles:
        if p is ref:
            shift, peak, uncertain = 0.0, 1.0, False
        else:
            shift, peak, uncertain = estimate_shift(
                ref["diameter_px_smooth"], p["diameter_px_smooth"], cfg
            )
        shifts.append(
            {
                "replicate": int(p["replicate"]),
                "shift_px": float(shift),
                "corr_peak": float(peak),
                "registration_uncertain": bool(uncertain),
            }
        )
        aligned.append((p["x"].astype(np.float64) + shift, p["diameter_px_raw"].astype(np.float64)))

    # common integer grid spanning the union of shifted, valid x
    lo = int(np.floor(min(ax[np.isfinite(d)].min() for ax, d in aligned if np.isfinite(d).any())))
    hi = int(np.ceil(max(ax[np.isfinite(d)].max() for ax, d in aligned if np.isfinite(d).any())))
    grid = np.arange(lo, hi + 1)

    # resample each replicate onto the grid (NaN outside its own valid range)
    stack = np.full((len(aligned), grid.size), np.nan)
    for i, (ax, d) in enumerate(aligned):
        finite = np.isfinite(d)
        if finite.sum() < 2:
            continue
        axf, df = ax[finite], d[finite]
        order = np.argsort(axf)
        axf, df = axf[order], df[order]
        inside = (grid >= axf[0]) & (grid <= axf[-1])
        stack[i, inside] = np.interp(grid[inside], axf, df)

    n = np.sum(np.isfinite(stack), axis=0)
    import warnings
    with warnings.catch_warnings(), np.errstate(invalid="ignore"):
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean_px = np.nanmean(np.where(n >= 1, stack, np.nan), axis=0)
        std_px = np.nanstd(np.where(n >= 2, stack, np.nan), axis=0, ddof=1)
    std_px[n < 2] = np.nan

    mean_um = mean_px / cfg.ppu
    std_um = std_px / cfg.ppu
    var_um2 = std_um ** 2

    keep = n >= 1
    table = {
        "x_aligned_px": grid[keep],
        "mean_um": mean_um[keep],
        "std_um": std_um[keep],
        "var_um2": var_um2[keep],
        "n": n[keep],
    }

    # scalar summary over the overlap span (all replicates present)
    overlap = n >= max(2, len([p for p in profiles]))
    overlap_px = int(overlap.sum())
    if overlap_px > 0:
        pooled_mean = float(np.nanmean(mean_um[overlap]))
        pooled_std = float(np.nanmean(std_um[overlap]))
    else:
        # fall back to wherever >=1 replicate exists
        pooled_mean = float(np.nanmean(mean_um[keep])) if keep.any() else float("nan")
        pooled_std = float(np.nanmean(std_um[n >= 2])) if (n >= 2).any() else float("nan")
    cv = pooled_std / pooled_mean if pooled_mean and np.isfinite(pooled_mean) else float("nan")

    summary = {
        "mean_um": pooled_mean,
        "std_um": pooled_std,
        "cv": cv,
        "n_points": int(keep.sum()),
        "n_replicates_used": len(profiles),
        "overlap_px": overlap_px,
        "registration_uncertain": any(s["registration_uncertain"] for s in shifts),
    }
    return table, shifts, summary
