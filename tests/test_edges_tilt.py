"""Tilt-invariance tests: perpendicular diameter + axis-following averaging."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.ndimage import uniform_filter1d

from fibrecv.band import tilt_geometry
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
    the sheared window lands back on the line, so ALL columns -- including
    the image borders, where the window must extend along the axis rather
    than smear vertically -- are reproduced exactly; a straight column
    average would dilute the line to ~1/wcol.
    """
    H, W, wcol = 60, 41, 11
    D = np.zeros((H, W), dtype=np.float32)
    c = np.arange(W)
    D[10 + c, c] = 1.0
    A = _axis_average(D, 1.0, wcol)
    assert np.allclose(A, D, atol=1e-6)


def test_axis_average_even_wcol_keeps_exact_width():
    """An even wcol must mean the same window in both branches (the sheared
    path must not silently widen a user-configured even width to odd)."""
    rng = np.random.default_rng(2)
    D = rng.normal(size=(50, 80)).astype(np.float32)
    ref = uniform_filter1d(D, size=40, axis=1, mode="nearest")
    # a vanishingly small tilt takes the sheared path but is numerically
    # indistinguishable from horizontal -- the result must match size=40
    assert np.allclose(_axis_average(D, 1e-12, 40), ref, atol=1e-4)


def test_tilt_geometry_clamps_and_handles_nan():
    m, cth = tilt_geometry(float("nan"))
    assert m == 0.0 and cth == 1.0
    m, cth = tilt_geometry(0.3)
    assert m == 0.3 and cth == pytest.approx(1.0 / math.sqrt(1.09))
    # pathological band fits are clamped to 60 degrees, bounding every
    # 1/cos(tilt) compensation at 2x
    m, cth = tilt_geometry(1e9)
    assert m == pytest.approx(math.tan(math.radians(60.0)))
    assert cth == pytest.approx(0.5)
    m, cth = tilt_geometry(-1e9)
    assert m < 0.0 and cth == pytest.approx(0.5)


def test_tilt_invariance_and_absolute_width():
    """Median measured width must not drift with tilt (validated to 40 deg).

    Tolerances: pairwise spread <= 2 % relative (widen to at most 3 % only if
    the 0-degree reference itself sits within the absolute band below).
    The absolute band is wide because the boundary level intentionally sits
    partway down the smoothed wall (edge_z above local background), which on
    this fixture reads systematically wide -- the property under test is
    *invariance*, not absolute accuracy.
    """
    cfg = CONFIG()
    meds = {a: _median_diameter_px(_inclined_fibre(float(a)), cfg)
            for a in (0, 10, 20, 30, 40, 45)}
    ref = meds[0]
    for a, v in meds.items():
        assert v == pytest.approx(ref, rel=0.02), f"angle {a}: {v:.2f} vs {ref:.2f}"
    assert TRUE_W * 0.95 <= ref <= TRUE_W * 1.35
