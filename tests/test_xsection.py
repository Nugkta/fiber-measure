"""Unit tests for the pure cross-section math module (xsection.py).

All tests operate on hand-built arrays — no images, no I/O.
"""

from __future__ import annotations

import numpy as np
import pytest

from fibrecv.config import CONFIG
from fibrecv.xsection import (
    NOMINAL_ANGLES_DEG,
    PartStack,
    build_part_stack,
    fit_ellipse_projections,
    hexagon_area,
    hexagon_area_expected,
    pair_differences,
    predict_anisotropy,
    predict_phi_transfer,
    split_half_area,
)


# ---------------------------------------------------------------- helpers


def ellipse_widths(a: float, b: float, phi_deg: float,
                   theta_deg=NOMINAL_ANGLES_DEG) -> np.ndarray:
    """Exact projection widths of an (a, b, phi) ellipse at the given angles."""
    th = np.radians(np.asarray(theta_deg, dtype=float))
    ph = np.radians(phi_deg)
    return 2.0 * np.sqrt(a ** 2 * np.cos(th - ph) ** 2
                         + b ** 2 * np.sin(th - ph) ** 2)


def tile_W(w6: np.ndarray, n: int = 4) -> np.ndarray:
    return np.repeat(w6[:, None], n, axis=1).astype(float)


def _clip_poly(poly, nx, ny, h):
    """Sutherland–Hodgman clip of polygon by half-plane nx*x+ny*y <= h."""
    out = []
    m = len(poly)
    for i in range(m):
        p, q = poly[i], poly[(i + 1) % m]
        dp = nx * p[0] + ny * p[1] - h
        dq = nx * q[0] + ny * q[1] - h
        if dp <= 0:
            out.append(p)
        if (dp < 0 < dq) or (dq < 0 < dp):
            t = dp / (dp - dq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def clip_hexagon_area(h0: float, h1: float, h2: float) -> float:
    """Brute-force area of the 3-slab intersection (normals at 0/60/120 deg)."""
    poly = [(-1e5, -1e5), (1e5, -1e5), (1e5, 1e5), (-1e5, 1e5)]
    for d, h in enumerate((h0, h1, h2)):
        for s in (1.0, -1.0):
            nx = s * np.cos(np.radians(60.0 * d))
            ny = s * np.sin(np.radians(60.0 * d))
            poly = _clip_poly(poly, nx, ny, h)
    area = 0.0
    for i in range(len(poly)):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % len(poly)]
        area += x1 * y2 - x2 * y1
    return abs(area) / 2.0


def phi_dist(p: float, q: float) -> float:
    """Distance between axis orientations (period 180 deg)."""
    d = abs((p - q) % 180.0)
    return min(d, 180.0 - d)


# ---------------------------------------------------- ellipse fit: exact


@pytest.mark.parametrize("phi", [0.0, 30.0, 89.0, 90.0, 91.0, 179.0])
@pytest.mark.parametrize("ratio", [1.05, 1.3, 2.0])
def test_exact_recovery(phi, ratio):
    b = 80.0
    a = ratio * b
    W = tile_W(ellipse_widths(a, b, phi))
    fit = fit_ellipse_projections(W)
    assert fit.valid.all()
    assert np.allclose(fit.a, a, rtol=1e-9)
    assert np.allclose(fit.b, b, rtol=1e-9)
    assert np.allclose(fit.area, np.pi * a * b, rtol=1e-9)
    assert all(phi_dist(p, phi) < 1e-6 for p in fit.phi_deg)
    assert np.allclose(fit.rms_resid, 0.0, atol=1e-7)
    assert (fit.n_angles == 6).all()


def test_circle_phi_nan_area_exact():
    W = tile_W(ellipse_widths(90.0, 90.0, 0.0))
    fit = fit_ellipse_projections(W)
    assert fit.valid.all()
    assert np.isnan(fit.phi_deg).all()
    assert np.allclose(fit.a, 90.0) and np.allclose(fit.b, 90.0)
    assert np.allclose(fit.area, np.pi * 90.0 ** 2, rtol=1e-9)


def test_a_geq_b_always_under_noise():
    rng = np.random.default_rng(0)
    W = tile_W(ellipse_widths(100.0, 95.0, 20.0), n=500)
    W += rng.normal(0.0, 3.0, W.shape)
    fit = fit_ellipse_projections(W)
    ok = fit.valid
    assert ok.any()
    assert (fit.a[ok] >= fit.b[ok] - 1e-12).all()


# ------------------------------------------------- missing-data handling


def test_one_angle_nan_still_valid():
    W = tile_W(ellipse_widths(100.0, 80.0, 35.0), n=3)
    W[2, :] = np.nan
    fit = fit_ellipse_projections(W)
    assert fit.valid.all()
    assert (fit.n_angles == 5).all()
    assert np.allclose(fit.a, 100.0, rtol=1e-9)  # 5 exact points, exact fit
    assert np.isnan(fit.resid[2]).all()


def test_full_pair_nan_invalid():
    W = tile_W(ellipse_widths(100.0, 80.0, 35.0), n=3)
    W[1, :] = np.nan
    W[4, :] = np.nan  # direction 2 fully gone
    fit = fit_ellipse_projections(W)
    assert not fit.valid.any()
    assert np.isnan(fit.a).all() and np.isnan(fit.area).all()
    assert (fit.n_angles == 4).all()


def test_three_angles_one_per_pair_exact_zero_resid():
    W = tile_W(ellipse_widths(100.0, 80.0, 35.0), n=3)
    W[3:, :] = np.nan  # keep a1, a2, a3 only
    fit = fit_ellipse_projections(W)
    assert fit.valid.all()
    assert (fit.n_angles == 3).all()
    assert np.allclose(fit.a, 100.0, rtol=1e-9)
    finite = np.isfinite(fit.resid)
    assert np.allclose(fit.resid[finite], 0.0, atol=1e-7)


def test_negative_b2_invalid_not_clamped():
    # widths (100, 10, 10) at one angle per direction -> c0 < R -> b^2 < 0
    W = np.full((6, 2), np.nan)
    W[0] = 100.0
    W[1] = 10.0
    W[2] = 10.0
    fit = fit_ellipse_projections(W)
    assert not fit.valid.any()
    assert np.isnan(fit.b).all() and np.isnan(fit.area).all()
    assert (fit.n_angles == 3).all()


def test_mixed_column_patterns_in_one_call():
    w6 = ellipse_widths(100.0, 80.0, 35.0)
    W = tile_W(w6, n=4)
    W[2, 1] = np.nan                # col 1: 5 angles
    W[1, 2] = np.nan
    W[4, 2] = np.nan                # col 2: direction 2 gone -> invalid
    W[:, 3] = np.nan                # col 3: empty
    fit = fit_ellipse_projections(W)
    assert fit.valid.tolist() == [True, True, False, False]
    assert fit.n_angles.tolist() == [6, 5, 4, 0]
    assert np.allclose(fit.a[:2], 100.0, rtol=1e-9)


# ------------------------------------------------------- Monte-Carlo bias


def test_montecarlo_unbiased_well_separated():
    rng = np.random.default_rng(1)
    a, b, phi = 110.0, 80.0, 50.0
    W = tile_W(ellipse_widths(a, b, phi), n=3000)
    W += rng.normal(0.0, 1.0, W.shape)
    fit = fit_ellipse_projections(W)
    ok = fit.valid
    assert ok.mean() > 0.99
    assert abs(np.mean(fit.a[ok]) - a) < 0.2
    assert abs(np.mean(fit.b[ok]) - b) < 0.2
    assert abs(np.mean(fit.area[ok]) / (np.pi * a * b) - 1.0) < 0.005
    ph = np.radians(fit.phi_deg[ok])
    phi_mean = np.degrees(np.arctan2(np.mean(np.sin(2 * ph)),
                                     np.mean(np.cos(2 * ph))) / 2)
    assert phi_dist(phi_mean, phi) < 0.5


def test_montecarlo_circle_ratio_bias_positive_small():
    rng = np.random.default_rng(2)
    r = 95.0
    W = tile_W(ellipse_widths(r, r, 0.0), n=3000)
    W += rng.normal(0.0, 1.0, W.shape)
    fit = fit_ellipse_projections(W)
    ok = fit.valid
    ratio = fit.a[ok] / fit.b[ok]
    # Rician: the estimated anisotropy R-hat is biased upward at circularity
    assert 1.0 < np.mean(ratio) < 1.05
    assert abs(np.mean(fit.area[ok]) / (np.pi * r * r) - 1.0) < 0.005


# ----------------------------------------------------------- hexagon area


def test_hexagon_regular():
    W = np.full((6, 2), 100.0)  # all direction widths equal -> h = 50
    area, degen = hexagon_area(W)
    assert not degen.any()
    assert np.allclose(area, 2 * np.sqrt(3) * 50.0 ** 2, rtol=1e-12)


def test_hexagon_circle_anchor():
    r = 77.0
    W = tile_W(ellipse_widths(r, r, 0.0))
    area, degen = hexagon_area(W)
    assert np.allclose(np.pi * r * r / area, np.pi / (2 * np.sqrt(3)), rtol=1e-12)


def test_hexagon_matches_bruteforce_clip_random():
    rng = np.random.default_rng(3)
    for _ in range(50):
        h = rng.uniform(50.0, 100.0, 3)
        if not all(h[d] < h[(d + 1) % 3] + h[(d + 2) % 3] for d in range(3)):
            continue
        W = np.repeat((2 * h)[:, None], 2, axis=1)
        W = np.vstack([W, W])  # 6 rows: pairs repeat direction widths
        area, degen = hexagon_area(W)
        assert not degen.any()
        ref = clip_hexagon_area(*h)
        assert np.allclose(area, ref, rtol=1e-9)


def test_hexagon_degenerate_nan():
    # h = (10, 6, 4): h0 == h1 + h2 -> slab 0 non-binding at the boundary
    for h0 in (10.0, 10.1):
        W = np.repeat(np.array([2 * h0, 12.0, 8.0])[:, None], 2, axis=1)
        W = np.vstack([W, W])
        area, degen = hexagon_area(W)
        assert degen.all()
        assert np.isnan(area).all()
    # strictly inside -> valid and matches the clip
    W = np.repeat(np.array([19.0, 12.0, 8.0])[:, None], 2, axis=1)
    W = np.vstack([W, W])
    area, degen = hexagon_area(W)
    assert not degen.any()
    assert np.allclose(area, clip_hexagon_area(9.5, 6.0, 4.0), rtol=1e-9)


def test_hexagon_missing_direction_nan():
    W = np.full((6, 2), 100.0)
    W[2] = np.nan
    W[5] = np.nan
    area, degen = hexagon_area(W)
    assert np.isnan(area).all()


def test_hexagon_area_expected_matches_clip_on_ellipse():
    for a, b, phi in [(100.0, 80.0, 35.0), (90.0, 60.0, 0.0), (70.0, 65.0, 120.0)]:
        h = ellipse_widths(a, b, phi)[:3] / 2.0
        expected = hexagon_area_expected(np.array([a]), np.array([b]),
                                         np.array([phi]))
        assert np.allclose(expected, clip_hexagon_area(*h), rtol=1e-9)
        # and consistent with hexagon_area on the exact projection widths
        area, _ = hexagon_area(tile_W(ellipse_widths(a, b, phi), n=1))
        assert np.allclose(expected, area, rtol=1e-9)


def test_hexagon_area_expected_circle_phi_nan():
    r = 88.0
    expected = hexagon_area_expected(np.array([r]), np.array([r]),
                                     np.array([np.nan]))
    assert np.allclose(expected, 2 * np.sqrt(3) * r * r, rtol=1e-12)


# --------------------------------------------- split halves + pair diffs


def test_split_half_area_exact_halves_agree():
    W = tile_W(ellipse_widths(100.0, 80.0, 35.0), n=3)
    lo, hi = split_half_area(W)
    assert np.allclose(lo, np.pi * 100.0 * 80.0, rtol=1e-9)
    assert np.allclose(hi, np.pi * 100.0 * 80.0, rtol=1e-9)


def test_split_half_area_detects_half_disagreement():
    W = tile_W(ellipse_widths(100.0, 80.0, 35.0), n=3)
    W[3:, :] = tile_W(ellipse_widths(104.0, 84.0, 35.0), n=3)[3:, :]
    lo, hi = split_half_area(W)
    assert np.allclose(lo, np.pi * 100.0 * 80.0, rtol=1e-9)
    assert np.allclose(hi, np.pi * 104.0 * 84.0, rtol=1e-9)


def test_split_half_area_missing_angle_nan():
    W = tile_W(ellipse_widths(100.0, 80.0, 35.0), n=2)
    W[1, 0] = np.nan  # first half incomplete at col 0
    lo, hi = split_half_area(W)
    assert np.isnan(lo[0]) and np.isfinite(hi[0])
    assert np.isfinite(lo[1])


def test_pair_differences():
    W = tile_W(np.array([10.0, 20.0, 30.0, 11.0, 22.0, 33.0]), n=2)
    d = pair_differences(W)
    assert d.shape == (3, 2)
    assert np.allclose(d[:, 0], [-1.0, -2.0, -3.0])


# --------------------------------------------------- prediction criteria


def test_predict_anisotropy_elliptical_directional_wins():
    rng = np.random.default_rng(4)
    W = tile_W(ellipse_widths(110.0, 80.0, 35.0), n=2000)
    W += rng.normal(0.0, 1.0, W.shape)
    dir_err, circ_err = predict_anisotropy(W)
    rms_dir = np.sqrt(np.nanmean(dir_err ** 2))
    rms_circ = np.sqrt(np.nanmean(circ_err ** 2))
    assert rms_dir < 2.0          # ~ noise (sqrt(2) * sigma)
    assert rms_circ > 3 * rms_dir


def test_predict_anisotropy_circle_directional_does_not_win():
    """On a true circle the 5-sample isotropic mean must beat the 1-sample
    directional predictor: rms_dir ~ sqrt(2)*sigma, rms_circ ~ sqrt(1.2)*sigma.
    """
    rng = np.random.default_rng(5)
    W = tile_W(ellipse_widths(95.0, 95.0, 0.0), n=2000)
    W += rng.normal(0.0, 1.0, W.shape)
    dir_err, circ_err = predict_anisotropy(W)
    rms_dir = np.sqrt(np.nanmean(dir_err ** 2))
    rms_circ = np.sqrt(np.nanmean(circ_err ** 2))
    assert abs(rms_dir - np.sqrt(2.0)) < 0.1
    assert abs(rms_circ - np.sqrt(1.2)) < 0.1
    assert rms_circ < rms_dir


def test_predict_anisotropy_nan_partner_excluded():
    W = tile_W(ellipse_widths(110.0, 80.0, 35.0), n=4)
    W[3, 1] = np.nan  # partner of angle 0 missing at col 1
    dir_err, circ_err = predict_anisotropy(W)
    assert np.isnan(dir_err[0, 1]) and np.isnan(circ_err[0, 1])
    assert np.isnan(dir_err[3, 1]) and np.isnan(circ_err[3, 1])
    assert np.isfinite(dir_err[1, 1])  # other directions unaffected


def test_predict_phi_transfer_exact():
    W = tile_W(ellipse_widths(110.0, 80.0, 35.0), n=8)
    err_ell, err_circ, phi_hat = predict_phi_transfer(W)
    assert phi_dist(phi_hat, 35.0) < 1e-6
    even = np.zeros(8, bool)
    even[0::2] = True
    assert np.allclose(err_ell[:, even], 0.0, atol=1e-7)
    assert np.isnan(err_ell[:, ~even]).all()
    # the isotropic baseline is genuinely worse on an elliptical section
    assert np.sqrt(np.nanmean(err_circ ** 2)) > 5.0


def test_predict_phi_transfer_noisy_beats_circle():
    rng = np.random.default_rng(6)
    W = tile_W(ellipse_widths(110.0, 80.0, 35.0), n=2000)
    W += rng.normal(0.0, 1.0, W.shape)
    err_ell, err_circ, phi_hat = predict_phi_transfer(W)
    assert phi_dist(phi_hat, 35.0) < 2.0
    assert np.sqrt(np.nanmean(err_ell ** 2)) < 0.5 * np.sqrt(np.nanmean(err_circ ** 2))


def test_predict_phi_transfer_bisector_direction_skipped_not_garbage():
    # phi = 30 deg bisects the kept pair {0, 60} used when direction 3
    # (120 deg) is dropped: both kept z values coincide, the (c0, R) solve is
    # rank 1, and an unguarded lstsq returns a minimum-norm answer that is
    # ~29 px wrong on EXACT widths. Degenerate columns must be NaN, not emitted.
    W = tile_W(ellipse_widths(110.0, 80.0, 30.0), n=8)
    err_ell, err_circ, phi_hat = predict_phi_transfer(W)
    assert phi_dist(phi_hat, 30.0) < 1e-6
    even = np.zeros(8, bool)
    even[0::2] = True
    assert np.isnan(err_ell[[2, 5]][:, even]).all()
    # the two non-degenerate directions stay exact
    assert np.allclose(err_ell[[0, 1, 3, 4]][:, even], 0.0, atol=1e-7)


def test_predict_phi_transfer_near_bisector_guard_band():
    even = np.zeros(8, bool)
    even[0::2] = True
    # 2 deg off the bisector: z-spread ~0.12 < _MIN_Z_SPREAD -> still skipped
    W = tile_W(ellipse_widths(110.0, 80.0, 32.0), n=8)
    err_ell, _, _ = predict_phi_transfer(W)
    assert np.isnan(err_ell[[2, 5]][:, even]).all()
    # 10 deg off: z-spread ~0.59 -> solved, and exactly
    W2 = tile_W(ellipse_widths(110.0, 80.0, 40.0), n=8)
    err_ell2, _, _ = predict_phi_transfer(W2)
    assert np.allclose(err_ell2[[2, 5]][:, even], 0.0, atol=1e-7)


def test_predict_phi_transfer_unsolvable_returns_nan_phi():
    # every odd column misses direction 3 entirely -> no solvable phi fit;
    # phi_hat must be NaN (never a silent 0.0) and both error arrays all-NaN
    W = tile_W(ellipse_widths(110.0, 80.0, 35.0), n=8)
    W[[2, 5], 1::2] = np.nan
    err_ell, err_circ, phi_hat = predict_phi_transfer(W)
    assert np.isnan(phi_hat)
    assert np.isnan(err_ell).all()
    assert np.isnan(err_circ).all()


# ------------------------------------------------------- build_part_stack


# Aperiodic reference width curve. A periodic fixture aliases the correlation
# peak and a global ramp biases it toward zero lag (window-mean artifacts), so
# use pure smooth filtered noise, interpolated for float shifts.
_CURVE_X = np.arange(-2000, 5000)
_curve_rng = np.random.default_rng(42)
_raw = _curve_rng.normal(0.0, 1.0, _CURVE_X.size)
_kernel = np.exp(-0.5 * (np.arange(-60, 61) / 20.0) ** 2)
_smooth = np.convolve(_raw, _kernel / _kernel.sum(), mode="same")
_CURVE = 190.0 + 50.0 * _smooth


def _structured_width(x: np.ndarray) -> np.ndarray:
    return np.interp(x, _CURVE_X, _CURVE)


def _make_profiles(deltas: dict[int, float],
                   spans: dict[int, tuple[int, int]] | None = None,
                   flat: set[int] = frozenset(),
                   seed: int = 0):
    """Profiles whose content is w_ref(x - delta) on per-angle spans."""
    rng = np.random.default_rng(seed)
    profiles = {}
    for a, delta in deltas.items():
        x0, x1 = (spans or {}).get(a, (60, 2400))
        x = np.arange(x0, x1, dtype=float)
        if a in flat:
            w = np.full(x.size, 190.0)
        else:
            w = _structured_width(x - delta) + rng.normal(0.0, 0.15, x.size)
        profiles[a] = {"x": x, "w": w}
    return profiles


def test_build_part_stack_recovers_shifts_unequal_x0():
    deltas = {1: 0.0, 2: 40.0, 3: -70.0, 4: 15.0, 5: -220.0, 6: 310.0}
    spans = {1: (60, 2400), 2: (200, 2380), 3: (0, 2100),
             4: (150, 2500), 5: (300, 2450), 6: (100, 2300)}
    profiles = _make_profiles(deltas, spans)
    cfg = CONFIG()
    st = build_part_stack(3, 2, profiles, cfg)
    assert isinstance(st, PartStack)
    assert st.fiber == 3 and st.part == 2
    assert st.W.shape[0] == 6 and st.W.shape[1] == st.x.size
    by_angle = {s["angle"]: s for s in st.shifts}
    for a, delta in deltas.items():
        s = by_angle[a]
        assert not s["uncertain"], f"a{a} flagged uncertain"
        # content w_ref(x - delta) sits delta to the right -> shift back = -delta.
        # Tolerance 1.25 px: the NaN-fill constant extension biases the
        # correlation peak by up to ~1.1 px when spans differ (the synthetic
        # end-to-end pin is +/-1.5 px; ~1 px misalignment of a smooth width
        # profile is negligible for the fit).
        expected = -(delta - deltas[1])
        assert abs(s["shift_px"] - expected) <= 1.25, (a, s["shift_px"], expected)
    # aligned rows must all reproduce the reference curve on the grid
    ref = _structured_width(st.x)
    for k in range(6):
        row = st.W[k]
        finite = np.isfinite(row)
        assert finite.sum() > 1500
        assert np.nanmedian(np.abs(row[finite] - ref[finite])) < 0.5


def test_build_part_stack_missing_angle_row_nan():
    deltas = {1: 0.0, 2: 20.0, 4: -10.0, 5: 0.0, 6: 5.0}  # a3 absent
    profiles = _make_profiles(deltas)
    st = build_part_stack(1, 1, profiles, CONFIG())
    assert np.isnan(st.W[2]).all()
    by_angle = {s["angle"]: s for s in st.shifts}
    assert not by_angle[3]["present"]


def test_build_part_stack_flat_profile_uncertain():
    deltas = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0}
    profiles = _make_profiles(deltas, flat={5})
    st = build_part_stack(1, 1, profiles, CONFIG())
    by_angle = {s["angle"]: s for s in st.shifts}
    assert by_angle[5]["uncertain"]
    assert by_angle[5]["shift_px"] == 0.0


def test_build_part_stack_xsec_min_corr_wiring():
    """S1: the gate must read cfg.xsec_min_corr, not cfg.min_corr."""
    deltas = {1: 0.0, 2: 25.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0}
    profiles = _make_profiles(deltas)
    # absurdly high xsec gate -> everything uncertain, zero shifts
    st = build_part_stack(1, 1, profiles, CONFIG(xsec_min_corr=1.1))
    assert all(s["uncertain"] for s in st.shifts if s["angle"] != 1)
    assert all(s["shift_px"] == 0.0 for s in st.shifts)
    # high min_corr but permissive xsec_min_corr -> NOT gated
    st2 = build_part_stack(1, 1, profiles,
                           CONFIG(min_corr=1.1, xsec_min_corr=0.0))
    by_angle = {s["angle"]: s for s in st2.shifts}
    assert not by_angle[2]["uncertain"]
    assert abs(by_angle[2]["shift_px"] + 25.0) <= 0.5


def test_build_part_stack_max_shift_uses_xsec_bound():
    # a2's content is 600 px away: recoverable only if xsec_max_shift(800) is used
    deltas = {1: 0.0, 2: 600.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0}
    profiles = _make_profiles(deltas)
    cfg = CONFIG()
    assert cfg.max_shift == 400  # the aggregate-stage default stays untouched
    st = build_part_stack(1, 1, profiles, cfg)
    by_angle = {s["angle"]: s for s in st.shifts}
    assert abs(by_angle[2]["shift_px"] + 600.0) <= 0.5
