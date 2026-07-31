"""Tilt-invariance tests: perpendicular diameter + axis-following averaging."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.ndimage import uniform_filter1d

from fibrecv.compute import compute_measurement
from fibrecv.config import CONFIG
from fibrecv.edges import _axis_average

TRUE_W = 60.0  # true perpendicular width of the synthetic fibre (px)


def _inclined_fibre(angle_deg: float, width: float = TRUE_W,
                    W: int = 800, H: int = 700, seed: int = 0) -> np.ndarray:
    """Pink background + pale band at ``angle_deg``, anti-aliased boundary.

    The edge blends linearly over ~2 px of *perpendicular* distance so the
    boundary stays smooth at any angle (a hard-edged mask aliases when tilted
    and adds ~2 % noise to the measured width).
    """
    m = math.tan(math.radians(angle_deg))
    y, x = np.mgrid[0:H, 0:W].astype(np.float64)
    dist = np.abs(y - (m * (x - W / 2) + H / 2)) / math.sqrt(1 + m * m)
    t = np.clip((width / 2 + 1.0 - dist) / 2.0, 0.0, 1.0)[..., None]
    bg = np.array([0.95, 0.45, 0.65])   # saturated pink background
    fg = np.array([0.92, 0.90, 0.91])   # pale desaturated fibre
    img = bg + (fg - bg) * t
    rng = np.random.default_rng(seed)
    img = img + rng.normal(0.0, 0.006, img.shape)
    return np.clip(img, 0.0, 1.0).astype(np.float32)


def _median_diameter_px(rgb: np.ndarray, cfg: CONFIG) -> float:
    mr = compute_measurement(rgb, cfg, "synthetic")
    assert mr.res.valid.any()
    return float(np.nanmedian(mr.res.diameter_raw))


def test_axis_average_slope0_matches_uniform_filter():
    rng = np.random.default_rng(1)
    D = rng.normal(size=(50, 80)).astype(np.float32)
    ref = uniform_filter1d(D, size=41, axis=1, mode="nearest")
    assert np.allclose(_axis_average(D, 0.0, 41), ref)
    # non-finite slope degrades to the plain average too
    assert np.allclose(_axis_average(D, float("nan"), 41), ref)


def test_axis_average_follows_the_axis():
    """A pattern constant along a 45-degree axis must survive axis-averaging.

    D[r, c] = 1 exactly on the line r = 10 + c. With slope=1 every sample in
    the sheared window lands back on the line, so interior columns are
    reproduced exactly; a straight column average would dilute the line to
    ~1/wcol.
    """
    H, W, wcol = 60, 41, 11
    D = np.zeros((H, W), dtype=np.float32)
    c = np.arange(W)
    D[10 + c, c] = 1.0
    A = _axis_average(D, 1.0, wcol)
    hw = wcol // 2
    interior = slice(hw, W - hw)
    assert np.allclose(A[:, interior], D[:, interior], atol=1e-6)
