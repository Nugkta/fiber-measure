"""Pure multi-angle helper tests: profile extraction, preview, headless batch.

Fixtures (`_bright_fibre`, `multiangle_folder`) are module level so Task 3 can
reuse them for Streamlit AppTest cases without redefining anything.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import imageio.v3 as iio
import numpy as np
import pytest

from fibrecv.compute import compute_measurement
from fibrecv.config import CONFIG
from fibrecv.gui_app import (
    _MA_UM_PER_PX_DEFAULT,
    _profile_from_mr,
    _profiles_from_results,
    multiangle_preview,
    run_batch,
    run_multiangle_batch,
)
from fibrecv.xsection import NOMINAL_ANGLES_DEG

A_TRUE, B_TRUE, PHI_TRUE = 35.0, 25.0, 20.0


def ellipse_width(angle_index: int, a: float = A_TRUE, b: float = B_TRUE,
                  phi: float = PHI_TRUE) -> float:
    """Full width w(theta) at one nominal angle (w^2 = c0+c1cos2t+c2sin2t)."""
    th = np.radians(NOMINAL_ANGLES_DEG[angle_index - 1])
    ph = np.radians(phi)
    return float(2 * np.sqrt(a ** 2 * np.cos(th - ph) ** 2
                             + b ** 2 * np.sin(th - ph) ** 2))


def _bright_fibre(width_px: float, W: int = 500, H: int = 200,
                  seed: int = 0, noise: float = 0.0) -> np.ndarray:
    """Dark background, flat bright band of constant ``width_px`` (no
    x-modulation) -- the minimal fixture for pure a/b/phi recovery checks.
    ``noise=0`` (default) keeps every column numerically identical, which is
    what the recovery/n_uncertain tests below rely on being deterministic.
    """
    y = np.arange(H, dtype=float)[:, None]
    dist = np.abs(y - H / 2.0)
    t = np.clip((width_px / 2.0 + 1.0 - dist) / 2.0, 0.0, 1.0)
    t = np.repeat(t, W, axis=1)[..., None]
    bg = np.array([0.15, 0.14, 0.13])
    fg = np.array([0.85, 0.84, 0.83])
    img = bg + (fg - bg) * t
    if noise:
        rng = np.random.default_rng(seed)
        img = img + rng.normal(0.0, noise, img.shape)
    return np.clip(img, 0.0, 1.0).astype(np.float32)


def _measure_angle(angle: int, a: float = A_TRUE, b: float = B_TRUE,
                   phi: float = PHI_TRUE) -> object:
    """Run the real bright-mode pipeline on one angle's flat-width fixture."""
    width = ellipse_width(angle, a, b, phi)
    rgb = _bright_fibre(width, seed=angle)
    cfg = CONFIG(feature_mode="bright")
    return compute_measurement(rgb, cfg, name=f"C1_01_a{angle}_part1")


@pytest.fixture(scope="module")
def six_angle_mrs() -> dict:
    """MeasureResults for all six angles of the a=35/b=25/phi=20 fixture."""
    mrs = {}
    for a in range(1, 7):
        mr = _measure_angle(a)
        assert mr.meta["coverage"] >= 0.8, f"angle {a}: bad coverage"
        mrs[a] = mr
    return mrs


@pytest.fixture
def multiangle_folder(tmp_path: Path) -> Path:
    """C1 fiber 1 part 1, all six angles of the a=35/b=25/phi=20 fixture,
    plus a scale-bar twin ("s" suffix, skipped by discover_multiangle) and a
    stray non-multiangle file (skipped too)."""
    for a in range(1, 7):
        arr = (_bright_fibre(ellipse_width(a)) * 255).astype(np.uint8)
        iio.imwrite(tmp_path / f"C1_01_a{a}_part1.tiff", arr)
    arr = (_bright_fibre(ellipse_width(1), seed=99) * 255).astype(np.uint8)
    iio.imwrite(tmp_path / "C1_01_a1_part1s.tiff", arr)  # scale-bar twin
    arr = (_bright_fibre(40.0, seed=42) * 255).astype(np.uint8)
    iio.imwrite(tmp_path / "notes.png", arr)  # stray, not multiangle-parseable
    return tmp_path


# --------------------------------------------------------------------- #
# _profile_from_mr
# --------------------------------------------------------------------- #

def test_profile_from_mr_absolute_x_and_span_sliced():
    """Directly exercises the audit-critical math with a fake MeasureResult
    (rgb=None, like the cached folder path): x is ABSOLUTE (from bnd, never
    from mr.rgb), and the full-image-width valid/diameter_smooth arrays are
    sliced by the SAME span before masking."""
    W = 20
    valid = np.zeros(W, dtype=bool)
    valid[5:15] = True
    valid[10] = False  # one invalid column inside the span (absolute x=10)
    diameter_smooth = np.arange(W, dtype=float)
    mr = SimpleNamespace(
        rgb=None, D=None,
        bnd=SimpleNamespace(x0=5, x1=14),
        res=SimpleNamespace(valid=valid, diameter_smooth=diameter_smooth),
    )
    prof = _profile_from_mr(mr)
    assert prof["x"][0] == 5
    assert len(prof["x"]) == len(prof["w"]) == 10
    assert np.isnan(prof["w"][5])              # x=10 -> span-relative index 5
    assert prof["w"][0] == pytest.approx(5.0)
    assert prof["w"][-1] == pytest.approx(14.0)


def test_profile_from_mr_real_pipeline(six_angle_mrs):
    mr = six_angle_mrs[1]
    prof = _profile_from_mr(mr)
    assert prof["x"][0] == mr.bnd.x0
    assert len(prof["x"]) == len(prof["w"]) == mr.bnd.x1 - mr.bnd.x0 + 1
    assert np.isfinite(prof["w"]).any()


# --------------------------------------------------------------------- #
# _profiles_from_results
# --------------------------------------------------------------------- #

def test_profiles_from_results_all_pass_qc(six_angle_mrs):
    cfg = CONFIG(feature_mode="bright")
    profiles, excluded = _profiles_from_results(six_angle_mrs, cfg)
    assert set(profiles) == {1, 2, 3, 4, 5, 6}
    assert excluded == {}


def test_profiles_from_results_excludes_band_mismatch(six_angle_mrs):
    cfg = CONFIG(feature_mode="bright")
    mrs = dict(six_angle_mrs)
    mrs[5] = SimpleNamespace(res=SimpleNamespace(
        band_mismatch=True, coverage=mrs[5].res.coverage,
        diameter_smooth=mrs[5].res.diameter_smooth, valid=mrs[5].res.valid,
        anomaly=mrs[5].res.anomaly), bnd=mrs[5].bnd)
    profiles, excluded = _profiles_from_results(mrs, cfg)
    assert set(profiles) == {1, 2, 3, 4, 6}
    assert 5 in excluded
    assert "band mismatch" in excluded[5]  # friendly label, not raw flag name
    assert "band_mismatch" not in excluded[5]


# --------------------------------------------------------------------- #
# multiangle_preview
# --------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def six_angle_profiles(six_angle_mrs) -> dict:
    return {a: _profile_from_mr(mr) for a, mr in six_angle_mrs.items()}


def test_multiangle_preview_recovers_known_ellipse(six_angle_profiles):
    cfg = CONFIG(feature_mode="bright")
    mp = multiangle_preview(six_angle_profiles, cfg)
    assert mp.present == {1, 2, 3, 4, 5, 6}
    assert mp.n_directions == 3
    assert mp.fittable is True
    assert mp.stack is not None and mp.fit is not None
    assert mp.med["a"] == pytest.approx(A_TRUE, abs=1.5)
    assert mp.med["b"] == pytest.approx(B_TRUE, abs=1.5)
    assert mp.med["phi"] == pytest.approx(PHI_TRUE, abs=2.0)


def test_multiangle_preview_missing_angle_stays_fittable(six_angle_profiles):
    cfg = CONFIG(feature_mode="bright")
    profiles = {a: p for a, p in six_angle_profiles.items() if a != 5}
    mp = multiangle_preview(profiles, cfg)
    assert mp.present == {1, 2, 3, 4, 6}
    assert mp.fittable is True

    # n_uncertain must count only PRESENT-and-uncertain shifts: build_part_stack
    # marks the MISSING angle (5) uncertain=True too, and a naive sum over all
    # shifts would double-report it -- assert the two counts actually differ.
    raw_uncertain = sum(1 for s in mp.stack.shifts if s["uncertain"])
    angle5 = next(s for s in mp.stack.shifts if s["angle"] == 5)
    assert angle5["present"] is False and angle5["uncertain"] is True
    assert mp.med["n_uncertain"] == raw_uncertain - 1
    assert mp.med["n_uncertain"] == sum(
        1 for s in mp.stack.shifts if s["present"] and s["uncertain"])


def test_multiangle_preview_two_same_direction_not_fittable(six_angle_profiles):
    profiles = {a: p for a, p in six_angle_profiles.items() if a in (1, 4)}
    cfg = CONFIG(feature_mode="bright")
    mp = multiangle_preview(profiles, cfg)
    assert mp.present == {1, 4}
    assert mp.n_directions == 1     # angles 1 and 4 are a 180-degree pair
    assert mp.fittable is False


def test_multiangle_preview_empty_profiles():
    cfg = CONFIG(feature_mode="bright")
    mp = multiangle_preview({}, cfg)
    assert mp.stack is None
    assert mp.fit is None
    assert mp.present == set()
    assert mp.n_directions == 0
    assert mp.fittable is False
    assert mp.med == {}


# --------------------------------------------------------------------- #
# run_batch(aggregate=)
# --------------------------------------------------------------------- #

def test_run_batch_aggregate_false_skips_stage2(tmp_path):
    """aggregate=False skips run_aggregate.main entirely -- the multi-angle
    batch relies on this since stage-2 group registration is meaningless for
    angle data (run_xsection is its own aggregation stage)."""
    cfg = CONFIG(feature_mode="bright")
    out_root = tmp_path / "out"
    img_path = tmp_path / "C1_01_a1_part1.tiff"
    iio.imwrite(img_path, (_bright_fibre(ellipse_width(1)) * 255).astype(np.uint8))
    master, results = run_batch([img_path], out_root, cfg, jobs=1, aggregate=False)
    assert master.empty
    assert len(results) == 1
    assert not (out_root / "summary" / "master_summary.csv").exists()
    assert (out_root / "summary" / "run_config.json").exists()  # per-image side still runs


# --------------------------------------------------------------------- #
# run_multiangle_batch (headless)
# --------------------------------------------------------------------- #

def test_run_multiangle_batch_headless(multiangle_folder, tmp_path):
    cfg = CONFIG(feature_mode="bright")
    out_root = tmp_path / "out"
    summary, results, rc = run_multiangle_batch(
        multiangle_folder, "C1", out_root, cfg, _MA_UM_PER_PX_DEFAULT, jobs=1)
    assert rc == 0
    assert len(results) == 6  # scale-bar twin + notes.png excluded by discover_multiangle
    assert not any("error" in r for r in results)
    assert len(summary) == 1
    assert summary.iloc[0]["fiber"] == 1
    assert summary.iloc[0]["um_per_px"] == pytest.approx(_MA_UM_PER_PX_DEFAULT)


def test_run_multiangle_batch_requires_condition_string(multiangle_folder, tmp_path):
    cfg = CONFIG(feature_mode="bright")
    with pytest.raises(ValueError):
        run_multiangle_batch(multiangle_folder, None, tmp_path / "out", cfg,
                             _MA_UM_PER_PX_DEFAULT, jobs=1)
