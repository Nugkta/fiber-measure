"""Unit tests for fibrecv.refine (erf edge refinement).

Group A: config + identity wiring -- the off-path is bit-identical,
``MeasureResult`` still constructs without ``ref``, and the meta ``"refine"``
sub-dict is present, correctly placed and JSON-safe.
Group B: the 1-D fit and its gates, on synthetic profiles.
Group C: the per-column offset field (quorum, interpolation chains, fallback).
Group D: the full pipeline on ``_blurred_fibre`` images -- sigma-independence,
tilt invariance and sigma-map accuracy.
"""

from __future__ import annotations

import json
import math
from dataclasses import replace

import numpy as np
import pytest
from scipy.special import erf

from fibrecv.band import BandResult, locate_band, tilt_geometry
from fibrecv.compute import MeasureResult, compute_measurement
from fibrecv.config import CONFIG
from fibrecv.edges import EdgeResult, detect_edges
from fibrecv.features import rgb_to_desaturation
from fibrecv.qc import run_qc
from fibrecv.refine import (
    RefineResult,
    _block_profile,
    _erf_model,
    _fit_block,
    _interp_side,
    refine_edges,
)


def _blurred_fibre(sigma_px: float, angle_deg: float = 0.0, width: float = 60.0,
                   W: int = 800, H: int = 700, seed: int = 0,
                   noise: float = 0.006) -> np.ndarray:
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

    ``noise`` sets the RGB noise std, which fixes the *contrast* of the z-map
    ``D``: ``D`` is normalised by the background MAD, so its amplitude is
    roughly ``0.4/noise`` z-units. The 0.006 default is a near-noiseless
    fixture (amplitude ~68 z); group D uses 0.03 (amplitude ~13 z) to match the
    contrast of real MasP2 images, where the legacy fixed level of ``edge_z``
    sits partway up the wall rather than in its very foot.
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
    img = img + rng.normal(0.0, noise, img.shape)
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


def test_refine_on_moves_both_walls_inward_and_keeps_untouched_fields():
    """On-path: a fresh EdgeResult with shifted walls, the rest passed through."""
    cfg = replace(CONFIG(), refine_on=True)
    rgb = _blurred_fibre(sigma_px=2.0)
    D, bnd, edg = _pipeline_up_to_edges(rgb, cfg)

    edg2, ref = refine_edges(D, edg, bnd, cfg)

    assert edg2 is not edg
    assert ref.n_blocks > 0 and ref.n_pass_top > 0 and ref.n_pass_bot > 0
    # the legacy edges sit outside the true wall, so both move INTO the fibre
    assert np.nanmedian(edg2.y_top) > np.nanmedian(edg.y_top)
    assert np.nanmedian(edg2.y_bot) < np.nanmedian(edg.y_bot)
    # flags / amp / y_core / half_window are never touched by refinement
    np.testing.assert_array_equal(edg2.flags, edg.flags)
    np.testing.assert_array_equal(edg2.amp, edg.amp)
    np.testing.assert_array_equal(edg2.y_core, edg.y_core)
    assert edg2.half_window == edg.half_window


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

    assert refine_meta["n_blocks"] > 0
    assert 0.8 <= refine_meta["coverage_top"] <= 1.0
    assert 0.8 <= refine_meta["coverage_bot"] <= 1.0
    for key in ("median_sigma_top", "median_sigma_bot",
                "median_abs_t0_top", "median_abs_t0_bot"):
        assert isinstance(refine_meta[key], float), key
    # the fixture is rendered with a 2 px blur
    assert refine_meta["median_sigma_top"] == pytest.approx(2.0, abs=0.4)

    json.dumps(meta)  # must not raise


def test_meta_refine_stats_are_none_safe_when_nothing_is_refined():
    """Coverage must be 0.0 (never NaN) and the medians None when no fit passes."""
    cfg = replace(CONFIG(), refine_sigma_min=100.0)  # no fit can ever pass
    rgb = _blurred_fibre(sigma_px=2.0)
    mr = compute_measurement(rgb, cfg, name="synthetic")

    refine_meta = mr.meta["refine"]
    assert refine_meta["n_pass_top"] == 0 and refine_meta["n_pass_bot"] == 0
    assert refine_meta["coverage_top"] == 0.0
    assert refine_meta["coverage_bot"] == 0.0
    assert refine_meta["median_sigma_top"] is None
    assert refine_meta["median_abs_t0_bot"] is None
    json.dumps(mr.meta)


def test_meta_refine_enabled_reflects_config_when_off():
    cfg = replace(CONFIG(), refine_on=False)
    rgb = _blurred_fibre(sigma_px=2.0)
    mr = compute_measurement(rgb, cfg, name="synthetic")
    assert mr.meta["refine"]["enabled"] is False
    json.dumps(mr.meta)


# --------------------------------------------------------------------------- #
# group B -- the 1-D fit and its gates                                        #
# --------------------------------------------------------------------------- #
def _t_grid(out: float = 35.0, inside: float = 28.0) -> np.ndarray:
    """The production profile abscissa: -out .. +inside in 0.5 px steps."""
    n = int(round((out + inside) / 0.5)) + 1
    return (-out + 0.5 * np.arange(n)).astype(np.float32)


def _profile(t: np.ndarray, t0: float = 0.0, sigma: float = 4.0, a: float = 0.0,
             b: float = 12.0, seed: int = 0, noise_frac: float = 1 / 40,
             n_cols: int = CONFIG.refine_block) -> np.ndarray:
    """A block's mean profile: blurred step + per-column noise ``(b-a)*noise_frac``.

    ``_fit_block`` never sees a single column -- it sees the mean over a block
    of ``refine_block`` columns, so the per-column noise amplitude of
    ``(b-a)/40`` (roughly twice what real MasP2 blocks show) is averaged down
    by ``sqrt(n_cols)`` here exactly as it is in ``_block_profile``.
    """
    rng = np.random.default_rng(seed)
    noise = rng.normal(0.0, (b - a) * noise_frac, (t.size, n_cols)).mean(axis=1)
    return (_erf_model(t, a, b, t0, sigma) + noise).astype(np.float32)


def _legacy_crossing(t: np.ndarray, prof: np.ndarray, cfg: CONFIG) -> float:
    """The edges.py rule: first crossing of ``base + min(edge_z, edge_frac*A)``.

    ``base``/``A`` are read off the profile's two ends exactly as the refine
    fit's own initial guess does, so the two estimators see identical data.
    """
    base = float(prof[:6].mean())
    amp = float(prof[-6:].mean()) - base
    level = base + min(cfg.edge_z, cfg.edge_frac * amp)
    above = np.where(prof >= level)[0]
    i = int(above[0])
    if i == 0:
        return float(t[0])
    lo, hi = float(prof[i - 1]), float(prof[i])
    frac = (level - lo) / (hi - lo)
    return float(t[i - 1]) + frac * float(t[i] - t[i - 1])


@pytest.mark.parametrize("sigma", [1.0, 2.0, 4.0, 8.0, 12.0, 15.0])
@pytest.mark.parametrize("t0", [-4.0, 0.0, 5.0])
def test_fit_recovers_the_midpoint_across_sigma_and_offset(sigma, t0):
    """The fitted midpoint is the true step position to <= 0.3 px, any blur."""
    cfg = CONFIG()
    t = _t_grid()
    fit = _fit_block(t, _profile(t, t0=t0, sigma=sigma, seed=int(sigma * 10)), cfg)
    assert fit is not None, f"fit rejected for sigma={sigma}, t0={t0}"
    t0_fit, sigma_fit, resid = fit
    assert abs(t0_fit - t0) <= 0.3, f"sigma={sigma}, t0={t0}: fitted {t0_fit:.3f}"
    assert sigma_fit == pytest.approx(sigma, rel=0.15), f"sigma {sigma_fit:.3f} vs {sigma}"
    assert resid < cfg.refine_relmax


def test_erf_midpoint_is_blur_invariant_where_the_legacy_level_drifts():
    """The contrast that motivates the whole stage.

    Same true edge at t=0, blur swept 1 -> 15 px: the fitted erf midpoint must
    stay put while the fixed-level crossing walks outward with sigma. That
    drift IS the focus-dependent diameter bias refinement removes.
    """
    cfg = CONFIG()
    t = _t_grid()
    sigmas = [1.0, 2.0, 4.0, 8.0, 12.0, 15.0]
    fitted, legacy = [], []
    for s in sigmas:
        prof = _profile(t, t0=0.0, sigma=s, seed=int(s * 10))
        fit = _fit_block(t, prof, cfg)
        assert fit is not None, f"fit rejected for sigma={s}"
        fitted.append(fit[0])
        legacy.append(_legacy_crossing(t, prof, cfg))

    spread_fit = max(fitted) - min(fitted)
    spread_legacy = max(legacy) - min(legacy)
    assert spread_fit <= 0.3, f"erf midpoints drifted {spread_fit:.2f} px: {fitted}"
    assert max(abs(v) for v in fitted) <= 0.3, fitted
    assert spread_legacy > 2.0, f"fixed-level crossings only drifted {spread_legacy:.2f} px"
    # and the drift is monotone outward (more blur -> earlier crossing)
    assert legacy[-1] < legacy[0] - 2.0, legacy


def test_fit_rejects_a_specular_bump_on_the_wall():
    """A bright specular blip inside the fibre must fail the residual gate."""
    cfg = CONFIG()
    t = _t_grid()
    clean = _profile(t, t0=0.0, sigma=4.0)
    assert _fit_block(t, clean, cfg) is not None  # the same profile passes clean
    bump = (0.6 * 12.0 * np.exp(-0.5 * ((t - 12.0) / 5.0) ** 2)).astype(np.float32)
    assert _fit_block(t, (clean + bump).astype(np.float32), cfg) is None


def test_fit_rejects_a_shadow_ramp_outside_the_wall():
    """A vignette/shadow ramp on the background side is not a blurred step."""
    cfg = CONFIG()
    t = _t_grid()
    clean = _profile(t, t0=0.0, sigma=4.0)
    ramp = (0.4 * 12.0 * np.clip((-t - 10.0) / 25.0, 0.0, 1.0)).astype(np.float32)
    assert _fit_block(t, (clean + ramp).astype(np.float32), cfg) is None


def test_fit_rejects_an_inverted_step():
    """b - a <= 0: desaturation must RISE into the fibre."""
    cfg = CONFIG()
    t = _t_grid()
    assert _fit_block(t, _profile(t, t0=0.0, sigma=4.0)[::-1].copy(), cfg) is None


@pytest.mark.parametrize("sigma", [0.25, 25.0])
def test_fit_rejects_sigma_outside_the_accepted_range(sigma):
    """Blur below refine_sigma_min or above refine_sigma_max -> no refinement."""
    cfg = CONFIG()
    t = _t_grid()
    assert not cfg.refine_sigma_min <= sigma <= cfg.refine_sigma_max
    assert _fit_block(t, _profile(t, t0=0.0, sigma=sigma), cfg) is None
    # ... and the same profile is accepted once the range is widened to admit it
    wide = replace(cfg, refine_sigma_min=0.01, refine_sigma_max=40.0)
    fit = _fit_block(t, _profile(t, t0=0.0, sigma=sigma), wide)
    assert fit is not None and fit[1] == pytest.approx(sigma, rel=0.25)


def test_fit_rejects_a_shift_beyond_maxshift():
    """|t0| > refine_maxshift: the legacy edge is too far away to be trusted."""
    cfg = CONFIG()
    t = _t_grid()
    prof = _profile(t, t0=18.0, sigma=3.0)
    assert _fit_block(t, prof, cfg) is None
    wide = replace(cfg, refine_maxshift=25.0)
    fit = _fit_block(t, prof, wide)
    assert fit is not None and fit[0] == pytest.approx(18.0, abs=0.3)


# --------------------------------------------------------------------------- #
# group C -- the per-column offset field                                      #
# --------------------------------------------------------------------------- #
def _slab_case(H: int = 240, W: int = 160, top: float = 80.0, bot: float = 160.0,
               sigma: float = 3.0, amp: float = 12.0, shift: float = 3.0,
               seed: int = 0, noise: float = 0.3, slope: float = 0.0):
    """A slab z-map plus the (band, edges) a legacy pass would give.

    At ``slope = 0`` the TRUE walls sit at rows ``top``/``bot`` with an erf blur
    of ``sigma``. At ``slope != 0`` the whole slab is sheared about the image
    centre column, keeping its VERTICAL thickness ``bot - top`` (so
    ``band_half``, a vertical count in ``band.py``, is unchanged) while its
    perpendicular width becomes ``(bot - top) * cth`` -- the blur ``sigma`` and
    the offset ``shift`` stay PERPENDICULAR quantities, exactly the frame
    ``refine.py`` fits in.

    The supplied ``EdgeResult`` places each wall ``shift`` perpendicular px
    OUTSIDE the true one (i.e. ``shift / cth`` vertical px), so a correct
    refinement must recover ``t0 = +shift`` on both walls, move each wall by
    ``shift / cth`` VERTICAL px onto the true row, and pull the perpendicular
    diameter back to ``(bot - top) * cth``. ``W`` is a multiple of
    ``refine_block``, so there are exactly ``W // 16`` blocks with centres at
    8, 24, 40, ...
    """
    cth = 1.0 / math.sqrt(1.0 + slope * slope)
    mid = float(top + bot) / 2.0
    half = (bot - top) / 2.0                    # vertical half-thickness
    x = np.arange(W, dtype=np.float64)[None, :]
    y = np.arange(H, dtype=np.float64)[:, None]
    centre = mid + slope * (x - (W - 1) / 2.0)  # tilted centerline row
    dist = (y - centre) * cth                   # PERPENDICULAR distance from it
    s = sigma * math.sqrt(2.0)
    blend = 0.5 * (erf((half * cth + dist) / s) + erf((half * cth - dist) / s))
    rng = np.random.default_rng(seed)
    D = (amp * blend + rng.normal(0.0, noise, (H, W))).astype(np.float32)

    c_fit = centre[0].astype(np.float32)
    bnd = BandResult(
        mask=np.zeros((H, W), dtype=bool),
        c_fit=c_fit,
        slope=slope,
        intercept=float(c_fit[0]),
        band_half=half,
        x0=0,
        x1=W - 1,
        centroid=c_fit.copy(),
        low_confidence=False,
        n_components=1,
    )
    y_top = (c_fit - half - shift / cth).astype(np.float32)
    y_bot = (c_fit + half + shift / cth).astype(np.float32)
    edg = EdgeResult(
        y_top=y_top,
        y_bot=y_bot,
        diameter=((y_bot - y_top) * cth).astype(np.float32),
        amp=np.full(W, amp, dtype=np.float32),
        y_core=c_fit.copy(),
        flags=np.zeros(W, dtype=np.int32),
        half_window=60,
    )
    return D, bnd, edg


def test_interp_side_bridges_a_short_gap_breaks_a_long_one_and_never_extrapolates():
    cfg = CONFIG()  # refine_gap_blocks = 2
    assert cfg.refine_gap_blocks == 2
    # blocks 0, 3 and 7 passed: 2 failed blocks between 0 and 3 (bridged),
    # 3 failed blocks between 3 and 7 (chain broken)
    out = _interp_side(
        160,
        np.array([0, 3, 7]),
        np.array([8, 56, 120]),
        np.array([[0.0, 3.0, 10.0]]),
        cfg,
    )[0]

    assert np.isnan(out[:8]).all(), "extrapolated before the first passing centre"
    assert np.all(np.isfinite(out[8:57])), "the 2-block gap was not bridged"
    assert out[8] == pytest.approx(0.0)
    assert out[32] == pytest.approx(1.5)   # linear halfway between the centres
    assert out[56] == pytest.approx(3.0)
    assert np.isnan(out[57:120]).all(), "the 3-block gap was bridged but must not be"
    assert out[120] == pytest.approx(10.0)  # a lone chain covers only its own centre
    assert np.isnan(out[121:]).all(), "extrapolated past the last passing centre"


def test_refine_recovers_a_known_offset_and_recomputes_the_diameter():
    """The whole point: both walls move inward by t0 and diameter follows.

    Written as an explicit check on ``edg2.diameter`` because ``qc.run_qc``
    consumes that array directly -- a refinement that moves y_top/y_bot but
    leaves the stale diameter behind would change nothing downstream.
    """
    cfg = CONFIG()
    D, bnd, edg = _slab_case(shift=3.0, sigma=3.0)
    edg2, ref = refine_edges(D, edg, bnd, cfg)

    assert ref.n_blocks == 10
    assert ref.n_pass_top == 10 and ref.n_pass_bot == 10
    r = ref.refined_top
    assert r.sum() >= 140  # columns 8..152 inclusive
    assert np.allclose(ref.o_top[r], 3.0, atol=0.3), np.nanmedian(ref.o_top)
    assert np.allclose(ref.o_bot[r], 3.0, atol=0.3), np.nanmedian(ref.o_bot)
    assert np.allclose(ref.sigma_top[r], 3.0, rtol=0.15)
    assert np.all(ref.resid_top[r] < cfg.refine_relmax)

    # top edge DOWN, bottom edge UP -- both into the fibre
    assert np.allclose(edg2.y_top[r], 80.0, atol=0.3)
    assert np.allclose(edg2.y_bot[r], 160.0, atol=0.3)
    # diameter recomputed from the shifted walls (86 -> 80), not left stale
    np.testing.assert_allclose(edg2.diameter[r], edg2.y_bot[r] - edg2.y_top[r], rtol=1e-6)
    assert np.allclose(edg2.diameter[r], 80.0, atol=0.6)
    assert np.all(edg2.diameter[r] < edg.diameter[r] - 4.0)


def test_perpendicular_offsets_are_converted_to_vertical_shifts_exactly_once():
    """The tilt conversion of global constraint #1, on a slab that can see it.

    Every other group-C case is horizontal, where ``cth == 1`` and dropping a
    ``/cth`` is invisible. Here the slab is tilted hard enough (slope 0.6,
    ``cth = 0.857``) that the vertical shift and the perpendicular offset differ
    by 0.83 px -- far outside the fit's ~0.05 px noise -- so an offset applied
    unconverted, or converted twice, moves the wall off its true row.
    """
    cfg = CONFIG()
    slope, shift = 0.6, 5.0
    D, bnd, edg = _slab_case(H=320, top=120.0, bot=200.0, slope=slope, shift=shift)
    _m, cth = tilt_geometry(bnd.slope)
    assert cth == pytest.approx(1.0 / math.sqrt(1.0 + slope ** 2))
    assert shift / cth - shift > 0.8  # the margin this test lives on

    edg2, ref = refine_edges(D, edg, bnd, cfg)
    rt, rb = ref.refined_top, ref.refined_bot
    assert rt.sum() >= 140 and rb.sum() >= 140

    # 1. the FITTED offset is perpendicular: it recovers the rendered shift
    assert np.allclose(ref.o_top[rt], shift, atol=0.3), np.nanmedian(ref.o_top)
    assert np.allclose(ref.o_bot[rb], shift, atol=0.3), np.nanmedian(ref.o_bot)
    assert np.allclose(ref.sigma_top[rt], 3.0, rtol=0.15)

    # 2. the APPLIED shift is vertical -- o/cth, the conversion done once ...
    np.testing.assert_allclose(edg2.y_top[rt] - edg.y_top[rt], ref.o_top[rt] / cth,
                               rtol=1e-5, err_msg="top shift is not o_top/cth")
    np.testing.assert_allclose(edg.y_bot[rb] - edg2.y_bot[rb], ref.o_bot[rb] / cth,
                               rtol=1e-5, err_msg="bot shift is not o_bot/cth")
    # ... and it is measurably NOT the unconverted offset
    assert np.all(np.abs(np.abs(edg2.y_top[rt] - edg.y_top[rt]) - ref.o_top[rt]) > 0.5)

    # 3. ground truth, independent of the formulas above: each wall lands on its
    #    true row, which an unconverted (or doubly converted) shift misses by 0.83
    true_top = bnd.c_fit - bnd.band_half
    true_bot = bnd.c_fit + bnd.band_half
    assert np.allclose(edg2.y_top[rt], true_top[rt], atol=0.3), (
        f"top off by {np.max(np.abs(edg2.y_top[rt] - true_top[rt])):.3f} px"
    )
    assert np.allclose(edg2.y_bot[rb], true_bot[rb], atol=0.3)

    # 4. and the diameter is the PERPENDICULAR chord of the shifted walls
    both = rt & rb
    np.testing.assert_allclose(edg2.diameter[both],
                               (edg2.y_bot[both] - edg2.y_top[both]) * cth, rtol=1e-6)
    assert np.allclose(edg2.diameter[both], 2 * bnd.band_half * cth, atol=0.6)


def test_a_two_block_hole_is_bridged_and_a_three_block_hole_is_not():
    cfg = CONFIG()
    for n_bad, bridged in ((2, True), (3, False)):
        D, bnd, edg = _slab_case()
        bad = slice(4 * 16, (4 + n_bad) * 16)
        D[:, bad] = -D[:, bad]  # inverted step -> b - a <= 0 -> those blocks fail
        edg2, ref = refine_edges(D, edg, bnd, cfg)

        assert ref.n_pass_top == 10 - n_bad, f"n_bad={n_bad}: {ref.n_pass_top} passed"
        hole = np.arange(4 * 16, (4 + n_bad) * 16)
        got = bool(np.isfinite(ref.o_top[hole]).all())
        assert got == bridged, f"n_bad={n_bad}: hole finite={got}, expected {bridged}"
        if not bridged:
            # only the chain interiors survive; the hole itself is untouched
            np.testing.assert_array_equal(edg2.y_top[hole], edg.y_top[hole])
            np.testing.assert_array_equal(edg2.diameter[hole], edg.diameter[hole])


@pytest.mark.parametrize("n_flagged, attempted", [(4, 10), (5, 9)])
def test_a_block_below_the_column_quorum_is_skipped(n_flagged, attempted):
    """>= 70% of a block's 16 columns must be anchors (11 is not enough, 12 is)."""
    cfg = CONFIG()
    D, bnd, edg = _slab_case()
    edg.flags[32:32 + n_flagged] = 4  # block 2 loses n_flagged anchor columns
    edg2, ref = refine_edges(D, edg, bnd, cfg)
    assert ref.n_blocks == attempted
    assert ref.n_pass_top == attempted


def test_flagged_columns_never_move_even_inside_a_refined_chain():
    cfg = CONFIG()
    D, bnd, edg = _slab_case()
    flagged = np.zeros(edg.flags.size, dtype=bool)
    flagged[[70, 71]] = True          # mid-chain, so interpolation covers them
    edg.flags[flagged] = 4
    edg2, ref = refine_edges(D, edg, bnd, cfg)

    assert ref.n_pass_top == 10       # 14/16 anchors still clears the quorum
    assert not ref.refined_top[flagged].any()
    assert not ref.refined_bot[flagged].any()
    assert np.isnan(ref.o_top[flagged]).all()
    np.testing.assert_array_equal(edg2.y_top[flagged], edg.y_top[flagged])
    np.testing.assert_array_equal(edg2.y_bot[flagged], edg.y_bot[flagged])
    np.testing.assert_array_equal(edg2.diameter[flagged], edg.diameter[flagged])
    assert ref.refined_top[72]        # ... while their neighbours did move
    assert edg2.y_top[72] != edg.y_top[72]


def test_no_extrapolation_beyond_the_outermost_passing_centres():
    cfg = CONFIG()
    D, bnd, edg = _slab_case()
    edg2, ref = refine_edges(D, edg, bnd, cfg)
    # block centres are 8 .. 152; nothing outside that span may move
    assert not ref.refined_top[:8].any() and not ref.refined_top[153:].any()
    np.testing.assert_array_equal(edg2.y_top[:8], edg.y_top[:8])
    np.testing.assert_array_equal(edg2.y_bot[153:], edg.y_bot[153:])
    assert ref.refined_top[8] and ref.refined_top[152]


def test_all_fits_failing_leaves_the_edges_untouched():
    cfg = replace(CONFIG(), refine_sigma_min=100.0)  # no fit can clear the gate
    D, bnd, edg = _slab_case()
    edg2, ref = refine_edges(D, edg, bnd, cfg)

    assert ref.n_blocks == 10
    assert ref.n_pass_top == 0 and ref.n_pass_bot == 0
    assert not ref.refined_top.any() and not ref.refined_bot.any()
    np.testing.assert_array_equal(edg2.y_top, edg.y_top)
    np.testing.assert_array_equal(edg2.y_bot, edg.y_bot)
    np.testing.assert_array_equal(edg2.diameter, edg.diameter)


def test_border_clipped_columns_are_dropped_before_the_quorum():
    """A wall whose outward window leaves the image is not fitted at all.

    Clipping the sample rows instead would silently fit replicated border rows.
    The other wall, whose window is entirely inside, must still be refined.
    """
    cfg = CONFIG()
    D, bnd, edg = _slab_case(top=20.0, bot=100.0)  # top anchor at row 17, out = 35
    assert 17.0 - cfg.refine_out < 0.0
    edg2, ref = refine_edges(D, edg, bnd, cfg)

    assert ref.n_pass_top == 0
    assert not ref.refined_top.any()
    np.testing.assert_array_equal(edg2.y_top, edg.y_top)
    assert ref.n_pass_bot == 10
    assert ref.refined_bot.any()
    assert np.allclose(edg2.y_bot[ref.refined_bot], 100.0, atol=0.3)
    # the block was still attempted -- the bottom wall carried it
    assert ref.n_blocks == 10


def test_block_profile_returns_none_below_the_quorum():
    """The quorum is over the block's columns, counted after the border drop."""
    cfg = CONFIG()
    D, _bnd, edg = _slab_case()
    t = _t_grid(cfg.refine_out, 28.0)
    cols = np.arange(0, 16)
    ya = edg.y_top.astype(np.float32).copy()
    assert _block_profile(D, ya, cols, t, 1.0, 1.0) is not None
    ya[:5] = np.nan  # 11 of 16 usable -> below 70%
    assert _block_profile(D, ya, cols, t, 1.0, 1.0) is None


# --------------------------------------------------------------------------- #
# group D -- full pipeline integration                                        #
# --------------------------------------------------------------------------- #
TRUE_W = 60.0  # rendered perpendicular width of _blurred_fibre (px)
D_NOISE = 0.03  # RGB noise giving a real-image-like z-map contrast (~13 z)


def _median_diameter_px(rgb: np.ndarray, cfg: CONFIG) -> float:
    mr = compute_measurement(rgb, cfg, "synthetic")
    assert mr.res.valid.any()
    return float(np.nanmedian(mr.res.diameter_raw))


def test_refined_diameter_is_independent_of_the_rendered_blur():
    """The headline claim: same fibre, three focus settings, one diameter.

    The legacy fixed level walks down the wall as the blur grows, so its
    measured width grows with sigma; the erf midpoint does not.
    """
    on, off = CONFIG(), replace(CONFIG(), refine_on=False)
    sigmas = (2.0, 6.0, 12.0)
    imgs = {s: _blurred_fibre(sigma_px=s, noise=D_NOISE) for s in sigmas}

    refined = [_median_diameter_px(imgs[s], on) for s in sigmas]
    legacy = [_median_diameter_px(imgs[s], off) for s in sigmas]

    spread_refined = max(refined) - min(refined)
    spread_legacy = max(legacy) - min(legacy)
    assert spread_refined <= 1.0, f"refined widths drifted with blur: {refined}"
    assert spread_legacy > 2.0, f"legacy widths were unexpectedly stable: {legacy}"
    # and the refined width lands on the rendered one
    for s, v in zip(sigmas, refined):
        assert abs(v - TRUE_W) <= 2.0, f"sigma={s}: refined {v:.2f} vs true {TRUE_W}"


@pytest.mark.parametrize("angle", [0.0, 20.0, 35.0])
def test_refined_diameter_is_tilt_invariant(angle):
    """Mirrors test_edges_tilt.py:96-112 with refinement on: the perpendicular
    frame must survive the extra stage (2% pairwise tolerance)."""
    cfg = CONFIG()
    ref = _median_diameter_px(_blurred_fibre(sigma_px=4.0, noise=D_NOISE), cfg)
    got = _median_diameter_px(
        _blurred_fibre(sigma_px=4.0, angle_deg=angle, noise=D_NOISE), cfg
    )
    assert got == pytest.approx(ref, rel=0.02), f"angle {angle}: {got:.2f} vs {ref:.2f}"


@pytest.mark.parametrize("sigma", [2.0, 6.0, 12.0])
def test_sigma_map_recovers_the_rendered_blur(sigma):
    """sigma(x) is the focus map the GUI plots -- it must mean what it says."""
    mr = compute_measurement(_blurred_fibre(sigma_px=sigma, noise=D_NOISE), CONFIG(), "syn")
    refine_meta = mr.meta["refine"]
    for key in ("median_sigma_top", "median_sigma_bot"):
        assert refine_meta[key] == pytest.approx(sigma, rel=0.2), (
            f"{key} = {refine_meta[key]:.2f}, rendered {sigma}"
        )
    assert refine_meta["coverage_top"] >= 0.9
    assert refine_meta["coverage_bot"] >= 0.9


def test_unrefined_columns_keep_their_legacy_edges_and_flags():
    """Refinement is strictly additive: what it does not touch stays byte-equal.

    (QC *validity* may still differ -- the handful of unrefined columns at the
    ends of the chain now sit ~3 px off their refined neighbours and the
    rolling-MAD filter drops them -- so this checks the guarantee refine.py
    actually makes: edges and flags of unrefined columns are untouched.)
    """
    rgb = _blurred_fibre(sigma_px=4.0, noise=D_NOISE)
    mr_on = compute_measurement(rgb, CONFIG(), "syn")
    mr_off = compute_measurement(rgb, replace(CONFIG(), refine_on=False), "syn")

    unrefined = ~(mr_on.ref.refined_top | mr_on.ref.refined_bot)
    assert unrefined.any() and not unrefined.all()
    for name in ("y_top", "y_bot", "diameter"):
        np.testing.assert_array_equal(
            getattr(mr_on.edg, name)[unrefined], getattr(mr_off.edg, name)[unrefined],
            err_msg=f"{name} changed on an unrefined column",
        )
    np.testing.assert_array_equal(mr_on.edg.flags, mr_off.edg.flags)
    assert mr_on.res.coverage >= mr_off.res.coverage - 0.05
