"""Locate the single full-width fibre band and fit its centerline.

Dependencies
------------
``numpy``, ``scipy.ndimage`` (closing, labelling), ``skimage.morphology``
(remove_small_objects), ``scipy.stats.theilslopes`` (robust line fit).

Inputs
------
- ``D``: desaturation z-map from ``features.rgb_to_desaturation``.
- ``CONFIG`` for ``k_band``, ``close_width``, ``min_object``, ``min_width``.

Output
------
``BandResult`` with: ``mask`` (selected component), ``c_fit`` (per-column
centerline, length W), ``slope``/``intercept`` (tilt), ``band_half`` (median
half-thickness), ``x0``/``x1`` (valid column span), ``centroid`` (per-column,
NaN where absent) and ``low_confidence`` flag.

Pos
---
Third stage of the per-image pipeline. The full-width (>=85% of image width)
component rule is what rejects left-edge graticule digits, dust blobs and
specular flecks. Feeds the search window + outlier rejection in ``edges``/``qc``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi
from scipy.stats import theilslopes
from skimage.morphology import remove_small_objects

from .config import CONFIG


@dataclass
class BandResult:
    mask: np.ndarray          # bool (H, W) selected component
    c_fit: np.ndarray         # float (W,) fitted centerline row per column
    slope: float              # centerline slope (px/px) -> tilt magnitude
    intercept: float
    band_half: float          # median half-thickness of the band (px)
    x0: int                   # first column the band spans
    x1: int                   # last column the band spans (inclusive)
    centroid: np.ndarray      # float (W,) per-column band centroid, NaN if absent
    low_confidence: bool      # no component reached min_width
    n_components: int         # labelled components in the coarse mask (diagnostic)


def coarse_mask(D: np.ndarray, cfg: CONFIG) -> np.ndarray:
    """Threshold ``D > k_band``, bridge specular gaps, drop small objects."""
    raw = D > cfg.k_band
    # horizontal closing to bridge specular highlights along the fibre
    structure = np.ones((1, max(1, cfg.close_width)), dtype=bool)
    closed = ndi.binary_closing(raw, structure=structure)
    closed = remove_small_objects(closed, min_size=cfg.min_object)
    return closed


def select_band(mask: np.ndarray, cfg: CONFIG) -> tuple[np.ndarray, bool, int]:
    """Keep the component spanning >= ``min_width`` of the image width.

    Returns ``(band_mask, low_confidence, n_components)``. If no component
    reaches the width threshold, the widest component is kept and
    ``low_confidence`` is True. An empty mask yields an all-False band.
    """
    labels, n = ndi.label(mask)
    if n == 0:
        return np.zeros_like(mask, dtype=bool), True, 0

    W = mask.shape[1]
    # column coverage per label (number of distinct columns the component touches)
    coverage = np.zeros(n + 1, dtype=int)
    for lbl in range(1, n + 1):
        cols = np.any(labels == lbl, axis=0)
        coverage[lbl] = int(cols.sum())

    widest = int(np.argmax(coverage[1:]) + 1)
    if coverage[widest] >= cfg.min_width * W:
        return labels == widest, False, n
    # nothing full-width: fall back to widest, flag it
    return labels == widest, True, n


def centerline_fit(band_mask: np.ndarray, cfg: CONFIG) -> BandResult:
    """Per-column centroid -> robust (Theil-Sen) line fit of the centerline.

    Also records the median band half-thickness (for sizing the vertical search
    window) and the column span ``[x0, x1]`` over which the band exists.
    """
    H, W = band_mask.shape
    rows = np.arange(H, dtype=np.float32)[:, None]
    col_counts = band_mask.sum(axis=0)
    present = col_counts > 0

    centroid = np.full(W, np.nan, dtype=np.float32)
    with np.errstate(invalid="ignore"):
        centroid[present] = (rows * band_mask).sum(axis=0)[present] / col_counts[present]

    xs = np.where(present)[0]
    if xs.size >= 2:
        ys = centroid[xs]
        slope, intercept, _, _ = theilslopes(ys, xs)
        x0, x1 = int(xs.min()), int(xs.max())
    else:
        # degenerate: flat line through image middle, whole width
        slope, intercept = 0.0, H / 2.0
        x0, x1 = 0, W - 1

    c_fit = (slope * np.arange(W, dtype=np.float32) + intercept).astype(np.float32)
    band_half = float(np.median(col_counts[present]) / 2.0) if present.any() else H / 4.0

    return BandResult(
        mask=band_mask,
        c_fit=c_fit,
        slope=float(slope),
        intercept=float(intercept),
        band_half=band_half,
        x0=x0,
        x1=x1,
        centroid=centroid,
        low_confidence=False,
        n_components=0,
    )


def locate_band(D: np.ndarray, cfg: CONFIG) -> BandResult:
    """Full band stage: coarse mask -> component selection -> centerline fit."""
    mask = coarse_mask(D, cfg)
    band_mask, low_conf, n = select_band(mask, cfg)
    res = centerline_fit(band_mask, cfg)
    res.low_confidence = low_conf
    res.n_components = n
    return res
