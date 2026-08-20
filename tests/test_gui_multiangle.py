"""Multi-angle tests: pure helpers (profiles, preview, headless batch) plus
the Streamlit AppTest cases for the multi-angle analysis mode.

Fixtures (`_bright_fibre`, `multiangle_folder`) are module level and shared by
both halves.
"""

from __future__ import annotations

import re
from pathlib import Path
from types import SimpleNamespace

import imageio.v3 as iio
import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

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

APP = str(Path(__file__).resolve().parents[1] / "src" / "fibrecv" / "gui_app.py")
MULTIANGLE_MODE = "Multi-angle cross-section"

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


def _write_angles(folder: Path, angles) -> Path:
    """Write one C1 fiber-1 part-1 image per requested angle."""
    for a in angles:
        arr = (_bright_fibre(ellipse_width(a)) * 255).astype(np.uint8)
        iio.imwrite(folder / f"C1_01_a{a}_part1.tiff", arr)
    return folder


@pytest.fixture
def multiangle_folder(tmp_path: Path) -> Path:
    """C1 fiber 1 part 1, all six angles of the a=35/b=25/phi=20 fixture,
    plus a scale-bar twin ("s" suffix, skipped by discover_multiangle) and a
    stray non-multiangle file (skipped too)."""
    _write_angles(tmp_path, range(1, 7))
    arr = (_bright_fibre(ellipse_width(1), seed=99) * 255).astype(np.uint8)
    iio.imwrite(tmp_path / "C1_01_a1_part1s.tiff", arr)  # scale-bar twin
    arr = (_bright_fibre(40.0, seed=42) * 255).astype(np.uint8)
    iio.imwrite(tmp_path / "notes.png", arr)  # stray, not multiangle-parseable
    return tmp_path


@pytest.fixture
def multiangle_folder_minus_a5(tmp_path: Path) -> Path:
    """Five angles: a5 missing, so the three projection directions are still
    covered (a1/a4 = 0, a2 = 60, a3/a6 = 120) but the split half a4-a6 is not."""
    return _write_angles(tmp_path, [1, 2, 3, 4, 6])


@pytest.fixture
def multiangle_folder_one_direction(tmp_path: Path) -> Path:
    """a1 + a4 only — a 180-degree pair, i.e. ONE projection direction."""
    return _write_angles(tmp_path, [1, 4])


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


# --------------------------------------------------------------------- #
# GUI: multi-angle analysis mode (streamlit AppTest)
# --------------------------------------------------------------------- #

def _mode_radio(at: AppTest):
    """The analysis-mode radio, selected BY KEY (never by index: the smoke
    test's at.sidebar.text_input[0] contract shows how brittle index-based
    selectors are, and _load_reps' own source radio carries no key)."""
    return next(r for r in at.sidebar.radio if r.key == "analysis_mode")


def _in_multiangle(folder: Path | None = None, timeout: int = 180) -> AppTest:
    """Boot the app, switch to multi-angle mode, optionally point it at a
    folder (the shared "Image folder" field, still text_input[0])."""
    at = AppTest.from_file(APP, default_timeout=timeout).run()
    assert not at.exception
    _mode_radio(at).set_value(MULTIANGLE_MODE).run()
    assert not at.exception
    if folder is not None:
        at.sidebar.text_input[0].set_value(str(folder)).run()
        assert not at.exception
    return at


def test_multiangle_mode_applies_bright_defaults():
    """Entering multi-angle mode couples the detector to bright mode once."""
    at = AppTest.from_file(APP, default_timeout=60).run()
    assert at.session_state["cfg_dict"]["feature_mode"] == "desat"

    _mode_radio(at).set_value(MULTIANGLE_MODE).run()
    assert not at.exception
    cfg = at.session_state["cfg_dict"]
    assert cfg["feature_mode"] == "bright"
    assert cfg["edge_frac"] == pytest.approx(0.30)   # BRIGHT_DEFAULTS
    assert cfg["k_band"] == pytest.approx(6.0)

    # and only on the TRANSITION: tuning a knob and re-running must not be
    # stomped back to the bright defaults on every interaction
    edge_z = next(s for s in at.sidebar.slider if s.key and "edge_z" in s.key)
    edge_z.set_value(6.5)
    next(b for b in at.sidebar.button if b.label == "Apply").click()
    at.run()
    assert not at.exception
    assert at.session_state["cfg_dict"]["edge_z"] == pytest.approx(6.5)


def test_reset_in_multiangle_mode_keeps_bright():
    """Reset is mode-aware: in multi-angle mode its base is the bright
    baseline, not the desat dataclass defaults (which would silently
    mis-detect every angle image)."""
    at = _in_multiangle(timeout=60)
    next(b for b in at.sidebar.button
         if b.label == "Reset to defaults").click()
    at.run()
    assert not at.exception
    cfg = at.session_state["cfg_dict"]
    assert cfg["feature_mode"] == "bright"
    assert cfg["edge_frac"] == pytest.approx(0.30)
    assert cfg["k_band"] == pytest.approx(6.0)

    # the same Reset in replicate mode still lands on the desat defaults
    _mode_radio(at).set_value("Replicates").run()
    next(b for b in at.sidebar.button
         if b.label == "Reset to defaults").click()
    at.run()
    assert not at.exception
    assert at.session_state["cfg_dict"]["feature_mode"] == "desat"
    assert at.session_state["cfg_dict"]["edge_frac"] == pytest.approx(0.65)


def test_multiangle_scale_survives_mode_round_trip():
    """A hand-calibrated µm/px must survive leaving and re-entering the mode.

    Streamlit drops a keyed widget's session entry on any run where the widget
    is not rendered, so without a plain mirror key the scale would silently
    revert to the C1 default — and it drives every µm number in the mode AND
    the batch's --scale-source.
    """
    at = _in_multiangle(timeout=60)
    scale = next(n for n in at.sidebar.number_input if n.key == "ma_um_per_px")
    assert scale.value == pytest.approx(_MA_UM_PER_PX_DEFAULT)
    scale.set_value(0.250000).run()
    assert not at.exception

    _mode_radio(at).set_value("Replicates").run()
    _mode_radio(at).set_value(MULTIANGLE_MODE).run()
    assert not at.exception
    scale = next(n for n in at.sidebar.number_input if n.key == "ma_um_per_px")
    assert scale.value == pytest.approx(0.250000)


def test_multiangle_six_angles_render_and_fit(multiangle_folder: Path):
    at = _in_multiangle(multiangle_folder)

    # card 01: one tab per angle, always six
    tab_labels = [t.label for t in at.tabs]
    for a in range(1, 7):
        assert any(lbl.endswith(f"a{a}") for lbl in tab_labels), tab_labels

    labels = [m.label for m in at.metric]
    assert "median width" in labels          # per-angle row (µm/px, not ppu)
    assert "median area" in labels           # cross-section row
    assert "axis ratio a/b" in labels
    assert "orientation φ" in labels
    assert "fit rms residual" in labels

    ratio = float(next(m.value for m in at.metric
                       if m.label == "axis ratio a/b"))
    assert ratio == pytest.approx(A_TRUE / B_TRUE, rel=0.10)

    # six angles present -> the split-half error is defined and shown
    area = next(m.value for m in at.metric if m.label == "median area")
    assert "µm²" in area and "±—" not in area

    # per-angle QC table
    tables = [df.value for df in at.dataframe]
    qc = [t for t in tables if "corr peak" in t.columns]
    assert qc, "per-angle alignment table not rendered"
    assert len(qc[0]) == 6

    captions = " | ".join(str(c.value) for c in at.caption)
    assert "µm/px (manual)" in captions
    # alignment shifts render signed at the "px2" kind's 2 dp, not 1 dp: a
    # sub-pixel shift is real information about how well the angles line up
    assert re.search(r"a1: [+-]\d+\.\d{2} px", captions), captions


def test_multiangle_missing_angle_warns_and_drops_split_half(
        multiangle_folder_minus_a5: Path):
    at = _in_multiangle(multiangle_folder_minus_a5)
    warnings_txt = " | ".join(str(w.value) for w in at.warning)
    assert "5 of 6 angles" in warnings_txt, warnings_txt
    assert "a5" in warnings_txt, warnings_txt
    assert not at.error, [e.value for e in at.error]

    # the fit still runs, but the split-half uncertainty is undefined
    area = next(m.value for m in at.metric if m.label == "median area")
    assert "±—" in area, area
    captions = " | ".join(str(c.value) for c in at.caption)
    assert "split-half uncertainty needs all six angles" in captions


def test_multiangle_batch_button_runs_both_stages(multiangle_folder: Path):
    """Card 03's one button drives measure + run_xsection and renders the
    summary (the helper itself is covered headlessly above; this covers the
    wiring: progress callback, status column, download button)."""
    at = _in_multiangle(multiangle_folder, timeout=600)
    out = multiangle_folder / "batch_out"
    next(t for t in at.text_input
         if t.label == "Output folder").set_value(str(out)).run()
    next(n for n in at.number_input
         if n.label == "parallel jobs").set_value(1).run()
    next(b for b in at.button
         if b.label.startswith("Run multi-angle")).click().run()
    assert not at.exception
    assert (out / "summary" / "xsection_summary.csv").exists()
    assert any("fibre row" in str(s.value) for s in at.success)
    tables = [df.value for df in at.dataframe]
    assert any("status" in t.columns for t in tables), [list(t) for t in tables]


def test_multiangle_upload_source_has_no_folder_upload(multiangle_folder: Path):
    """The upload source offers hand-picked files only: a folder chooser
    pointed at a real condition directory would push tens of GB of TIFFs
    through the browser websocket. The batch button needs a real directory
    too, so it is disabled here."""
    at = _in_multiangle(multiangle_folder)
    assert not any("upload a whole folder" in cb.label
                   for cb in at.sidebar.checkbox)
    assert next(b for b in at.button
                if b.label.startswith("Run multi-angle")).disabled is False

    source = next(r for r in at.sidebar.radio if r.label == "multi-angle source")
    source.set_value("Upload").run()
    assert not at.exception
    assert not any("upload a whole folder" in cb.label
                   for cb in at.sidebar.checkbox)
    # no uploads yet -> nothing is loaded and the cards do not render
    assert not any(b.label.startswith("Run multi-angle") for b in at.button)
    assert any("upload the six angle images" in str(i.value) for i in at.info)


def test_multiangle_one_direction_cannot_fit(
        multiangle_folder_one_direction: Path):
    at = _in_multiangle(multiangle_folder_one_direction)
    errors = " | ".join(str(e.value) for e in at.error)
    assert "Cannot fit" in errors, errors
    assert "1 of 3 projection directions" in errors, errors

    # no ellipse metrics, but the aligned plot and the QC table still render
    labels = [m.label for m in at.metric]
    assert "median area" not in labels
    tables = [df.value for df in at.dataframe]
    assert any("corr peak" in t.columns for t in tables)
