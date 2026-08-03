"""Unit tests for fibrecv.refine (erf edge refinement).

Group A (this file, M1): config + identity wiring -- the off-path is
bit-identical, ``MeasureResult`` still constructs without ``ref``, and the
meta ``"refine"`` sub-dict is present, correctly placed and JSON-safe.
Task 2 extends this file with groups B (fit/gates), C (offset field) and D
(full-pipeline integration), reusing ``_blurred_fibre`` below.
"""

from __future__ import annotations

import json
import math
from dataclasses import replace

import numpy as np
from scipy.special import erf

from fibrecv.band import BandResult, locate_band
from fibrecv.compute import MeasureResult, compute_measurement
from fibrecv.config import CONFIG
from fibrecv.edges import EdgeResult, detect_edges
from fibrecv.features import rgb_to_desaturation
from fibrecv.qc import run_qc
from fibrecv.refine import RefineResult, refine_edges


def _blurred_fibre(sigma_px: float, angle_deg: float = 0.0, width: float = 60.0,
                   W: int = 800, H: int = 700, seed: int = 0) -> np.ndarray:
    """Pink background + pale band at ``angle_deg``, a TRUE erf-blurred boundary.

    Adapted from ``_inclined_fibre`` (tests/test_edges_tilt.py:19-37), but here
    each wall is the exact convolution of a hard rectangular band with a
    Gaussian of std ``sigma_px`` (perpendicular px) -- i.e. the blend fraction
    ``t`` at signed perpendicular distance ``y`` from the centerline is
    ``0.5*(erf((a-y)/(sigma*sqrt2)) + erf((a+y)/(sigma*sqrt2)))`` with
    ``a = width/2`` -- matching the erf step model ``refine.py`` fits against,
    with a controllable, known ground-truth blur ``sigma_px``.

    ``fg = (0.95, 0.90, 0.91)`` keeps the R channel identical to ``bg``'s
    (0.95): R is the max channel at both t=0 and t=1 and G is the min at both
    (checked: G never overtakes B for width=60), so HSV saturation
    ``S = (R-G)/R`` is an EXACTLY AFFINE function of the blend fraction ``t``
    for every column, i.e. the erf shape in ``t`` becomes an erf shape in
    ``S``/``D`` unchanged -- no HSV nonlinearity to contaminate the ground truth.
    """
    m = math.tan(math.radians(angle_deg))
    y, x = np.mgrid[0:H, 0:W].astype(np.float64)
    # signed perpendicular distance from the centerline (positive = below/one side)
    dist = (y - (m * (x - W / 2) + H / 2)) / math.sqrt(1 + m * m)
    a = width / 2.0
    s = max(sigma_px, 1e-6) * math.sqrt(2.0)
    t = 0.5 * (erf((a - dist) / s) + erf((a + dist) / s))
    t = t[..., None]
    bg = np.array([0.95, 0.45, 0.65])   # saturated pink background
    fg = np.array([0.95, 0.90, 0.91])   # pale desaturated fibre core (R unchanged)
    img = bg + (fg - bg) * t
    rng = np.random.default_rng(seed)
    img = img + rng.normal(0.0, 0.006, img.shape)
    return np.clip(img, 0.0, 1.0).astype(np.float32)


def _pipeline_up_to_edges(rgb: np.ndarray, cfg: CONFIG) -> tuple[np.ndarray, BandResult, EdgeResult]:
    """features -> band -> edges only (no refine, no qc) -- for direct refine_edges tests."""
    D, _S, _s_bg, _mad = rgb_to_desaturation(rgb, cfg)
    bnd = locate_band(D, cfg)
    edg = detect_edges(D, bnd, cfg)
    return D, bnd, edg


# --------------------------------------------------------------------------- #
# off-path bit-identical + identity                                           #
# --------------------------------------------------------------------------- #
def test_refine_off_returns_same_edge_object_and_bit_identical_arrays():
    cfg = replace(CONFIG(), refine_on=False)
    rgb = _blurred_fibre(sigma_px=2.0)
    D, bnd, edg = _pipeline_up_to_edges(rgb, cfg)

    edg2, ref = refine_edges(D, edg, bnd, cfg)

    assert edg2 is edg  # identity, not a copy
    for name in ("y_top", "y_bot", "diameter", "amp", "y_core", "flags"):
        np.testing.assert_array_equal(getattr(edg2, name), getattr(edg, name))

    W = D.shape[1]
    empty = RefineResult.empty(W)
    assert ref.n_blocks == empty.n_blocks == 0
    assert ref.n_pass_top == empty.n_pass_top == 0
    assert ref.n_pass_bot == empty.n_pass_bot == 0
    np.testing.assert_array_equal(ref.refined_top, empty.refined_top)
    np.testing.assert_array_equal(ref.refined_bot, empty.refined_bot)
    assert np.isnan(ref.sigma_top).all() and np.isnan(ref.sigma_bot).all()
    assert np.isnan(ref.resid_top).all() and np.isnan(ref.resid_bot).all()
    assert np.isnan(ref.o_top).all() and np.isnan(ref.o_bot).all()


def test_refine_on_placeholder_is_also_a_noop_in_m1():
    """M1 resolution: with refine_on=True the on-path is a placeholder that
    behaves exactly like the off-path (same object, n_blocks=0) -- the real
    algorithm lands in Task 2."""
    cfg = replace(CONFIG(), refine_on=True)
    rgb = _blurred_fibre(sigma_px=2.0)
    D, bnd, edg = _pipeline_up_to_edges(rgb, cfg)

    edg2, ref = refine_edges(D, edg, bnd, cfg)

    assert edg2 is edg
    assert ref.n_blocks == 0


# --------------------------------------------------------------------------- #
# MeasureResult without ref                                                   #
# --------------------------------------------------------------------------- #
def test_measure_result_constructs_without_ref():
    """MeasureResult must still be constructable without `ref` (default None) --
    tests/test_manual_edit.py:80 relies on exactly this."""
    W, H = 40, 40
    bnd = BandResult(
        mask=np.zeros((H, W), dtype=bool),
        c_fit=np.full(W, 20.0),
        slope=0.0,
        intercept=20.0,
        band_half=5.0,
        x0=2,
        x1=37,
        centroid=np.full(W, 20.0),
        low_confidence=False,
        n_components=1,
    )
    edg = EdgeResult(
        y_top=np.full(W, 15.0),
        y_bot=np.full(W, 25.0),
        diameter=np.full(W, 10.0),
        amp=np.full(W, 5.0),
        y_core=np.full(W, 20.0),
        flags=np.zeros(W, dtype=np.int32),
        half_window=10,
    )
    cfg = CONFIG()
    res = run_qc(edg, bnd, cfg)
    mr = MeasureResult(
        rgb=None, D=None, bnd=bnd, edg=edg, res=res,
        diameter_um=np.full(W, np.nan), name=None, group=None, replicate=None,
    )
    assert mr.ref is None


# --------------------------------------------------------------------------- #
# meta "refine" sub-dict                                                      #
# --------------------------------------------------------------------------- #
def test_meta_refine_subdict_placement_and_json_safety():
    cfg = CONFIG()
    rgb = _blurred_fibre(sigma_px=2.0)
    mr = compute_measurement(rgb, cfg, name="synthetic")

    meta = mr.meta
    assert "refine" in meta
    keys = list(meta.keys())
    i_flags, i_refine, i_median = (
        keys.index("flag_counts"), keys.index("refine"), keys.index("median_diameter_um"),
    )
    assert i_flags < i_refine < i_median, keys

    refine_meta = meta["refine"]
    expected_keys = {
        "enabled", "n_blocks", "n_pass_top", "n_pass_bot",
        "coverage_top", "coverage_bot",
        "median_sigma_top", "median_sigma_bot",
        "median_abs_t0_top", "median_abs_t0_bot",
    }
    assert set(refine_meta.keys()) == expected_keys
    assert refine_meta["enabled"] is True

    # M1: nothing is ever refined yet, so the None-safe stats must all be None
    # and coverage must be 0.0 (never NaN, even though anchor columns exist).
    assert refine_meta["n_blocks"] == 0
    assert refine_meta["coverage_top"] == 0.0
    assert refine_meta["coverage_bot"] == 0.0
    assert refine_meta["median_sigma_top"] is None
    assert refine_meta["median_sigma_bot"] is None
    assert refine_meta["median_abs_t0_top"] is None
    assert refine_meta["median_abs_t0_bot"] is None

    json.dumps(meta)  # must not raise


def test_meta_refine_enabled_reflects_config_when_off():
    cfg = replace(CONFIG(), refine_on=False)
    rgb = _blurred_fibre(sigma_px=2.0)
    mr = compute_measurement(rgb, cfg, name="synthetic")
    assert mr.meta["refine"]["enabled"] is False
    json.dumps(mr.meta)
