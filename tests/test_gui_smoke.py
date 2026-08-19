"""Headless GUI smoke tests (streamlit AppTest): boot, full flow, manual edits."""

from __future__ import annotations

from pathlib import Path

import imageio.v3 as iio
import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).resolve().parents[1] / "src" / "fibrecv" / "gui_app.py")


def _synthetic_fibre(W: int = 800, H: int = 300, seed: int = 0) -> np.ndarray:
    """Saturated pink background with a pale desaturated full-width band."""
    rng = np.random.default_rng(seed)
    img = np.empty((H, W, 3), dtype=np.float32)
    img[:] = (0.95, 0.45, 0.65)
    img[130:170] = (0.92, 0.90, 0.91)
    img += rng.normal(0.0, 0.01, img.shape).astype(np.float32)
    return np.clip(img, 0.0, 1.0)


def _step_fibre(W: int = 800, H: int = 300, seed: int = 0) -> np.ndarray:
    """Like _synthetic_fibre but the band steps 40 -> 70 px tall at x=400
    (a ~55% level shift, safely above the calibrated step_frac=0.25)."""
    rng = np.random.default_rng(seed)
    img = np.empty((H, W, 3), dtype=np.float32)
    img[:] = (0.95, 0.45, 0.65)
    img[130:170, :400] = (0.92, 0.90, 0.91)
    img[115:185, 400:] = (0.92, 0.90, 0.91)
    img += rng.normal(0.0, 0.01, img.shape).astype(np.float32)
    return np.clip(img, 0.0, 1.0)


@pytest.fixture
def image_folder(tmp_path: Path) -> Path:
    for repn in (1, 2):
        arr = (_synthetic_fibre(seed=repn) * 255).astype(np.uint8)
        iio.imwrite(tmp_path / f"test 1_1_{repn}.png", arr)
    return tmp_path


@pytest.fixture
def anomaly_folder(tmp_path: Path) -> Path:
    """rep 1 clean, rep 2 with a mid-image diameter step."""
    arr = (_synthetic_fibre(seed=1) * 255).astype(np.uint8)
    iio.imwrite(tmp_path / "test 1_1_1.png", arr)
    arr = (_step_fibre(seed=2) * 255).astype(np.uint8)
    iio.imwrite(tmp_path / "test 1_1_2.png", arr)
    return tmp_path


def test_boot_without_data():
    at = AppTest.from_file(APP, default_timeout=60).run()
    assert not at.exception


def test_full_flow_with_manual_edits(image_folder: Path):
    at = AppTest.from_file(APP, default_timeout=180)
    at.run()
    at.sidebar.text_input[0].set_value(str(image_folder)).run()
    assert not at.exception
    labels = [m.label for m in at.metric]
    assert "group mean Ø" in labels
    assert "between-replicate std" in labels

    # per-replicate metric row: each image shows its own mean and spread
    assert "mean Ø" in labels
    assert "along-fibre std" in labels

    # group panel: per-image stats table with one row per replicate
    tables = [df.value for df in at.dataframe]
    stats = [t for t in tables if "mean Ø (µm)" in t.columns]
    assert stats, "per-image stats table not rendered"
    assert len(stats[0]) == 2
    assert np.isfinite(stats[0]["std (µm)"]).all()

    # inject a manual edit (what a click would store) and rerun: the apply
    # choke point in main() must run it through apply_manual_edits cleanly
    at.session_state["manual_edits"] = {
        "test 1_1_1": {
            "top": [[(100.0, 128.0), (300.0, 128.0)],
                    [(500.0, 127.0), (600.0, 127.0)]],
            "bot": [],
            "nudge_top": 0.0,
            "nudge_bot": 0.0,
        }
    }
    at.run()
    assert not at.exception
    captions = " | ".join(str(c.value) for c in at.caption)
    assert "edited columns" in captions
    img_captions = " | ".join(
        i.caption for el in at.get("imgs") for i in el.proto.imgs)
    assert "magenta = manual edit" in img_captions


def test_anomaly_flags_surface_in_gui(anomaly_folder: Path):
    at = AppTest.from_file(APP, default_timeout=180)
    at.run()
    at.sidebar.text_input[0].set_value(str(anomaly_folder)).run()
    assert not at.exception

    # the anomaly knobs live in their own collapsed sidebar expander
    exp_labels = [e.label for e in at.sidebar.expander]
    assert "Anomaly flags" in exp_labels, exp_labels

    # rep 2's diameter step earns a warning prefix on its tab label
    tab_labels = [t.label for t in at.tabs]
    assert any(lbl.startswith("⚠") for lbl in tab_labels), tab_labels
    assert any(not lbl.startswith("⚠") for lbl in tab_labels), tab_labels

    # flags badge on rep 2 names the anomaly
    # human-readable label in the badge, never the raw snake_case flag name
    flags_vals = [m.value for m in at.metric if m.label == "flags"]
    assert any("diameter step" in v for v in flags_vals), flags_vals
    assert not any("diameter_step" in v for v in flags_vals), flags_vals

    # per-image stats table gained the anomalies column (readable labels too)
    tables = [df.value for df in at.dataframe]
    stats = [t for t in tables if "anomalies" in t.columns]
    assert stats, "anomalies column missing from per-image stats"
    assert stats[0]["anomalies"].str.contains("diameter step").any()

    # advisory by default: nothing dropped from registration
    captions = " | ".join(str(c.value) for c in at.caption)
    assert "QC-dropped" not in captions

    # flip the exclude switch through the REAL widget path (checkbox inside
    # the expander inside the form, then Apply): rep 2 must be dropped with
    # an anomaly reason
    exclude_cb = next(cb for cb in at.sidebar.checkbox
                      if cb.key and "anomaly_exclude" in cb.key)
    exclude_cb.set_value(True)
    next(b for b in at.sidebar.button if b.label == "Apply").click()
    at.run()
    assert not at.exception
    assert at.session_state["cfg_dict"]["anomaly_exclude"] is True
    captions = " | ".join(str(c.value) for c in at.caption)
    assert "QC-dropped" in captions
    assert "test 1_1_2" in captions and "anomaly" in captions


def test_mode_switch_applies_calibrated_defaults():
    """Switching feature_mode through the real widget path re-applies that
    mode's calibrated edge_frac/k_band (bright: BRIGHT_DEFAULTS 0.30/6.0,
    desat: dataclass 0.65/4.0) while non-coupled user tuning survives."""
    at = AppTest.from_file(APP, default_timeout=60).run()
    assert not at.exception

    cfg = at.session_state["cfg_dict"]
    assert cfg["feature_mode"] == "desat"
    assert cfg["edge_frac"] == pytest.approx(0.65)
    assert cfg["k_band"] == pytest.approx(4.0)

    # tune a non-coupled knob, then switch to bright and Apply
    edge_z = next(s for s in at.sidebar.slider if s.key and "edge_z" in s.key)
    edge_z.set_value(5.0)
    mode_sb = next(sb for sb in at.sidebar.selectbox
                   if sb.key and "feature_mode" in sb.key)
    mode_sb.set_value("bright")
    next(b for b in at.sidebar.button if b.label == "Apply").click()
    at.run()
    assert not at.exception
    cfg = at.session_state["cfg_dict"]
    assert cfg["feature_mode"] == "bright"
    assert cfg["edge_frac"] == pytest.approx(0.30)  # BRIGHT_DEFAULTS
    assert cfg["k_band"] == pytest.approx(6.0)      # hidden knob rides along
    assert cfg["edge_z"] == pytest.approx(5.0)      # user tuning survives

    # switch back: desat calibrated defaults return (widgets re-keyed by the
    # form_version bump, so re-fetch them)
    mode_sb = next(sb for sb in at.sidebar.selectbox
                   if sb.key and "feature_mode" in sb.key)
    mode_sb.set_value("desat")
    next(b for b in at.sidebar.button if b.label == "Apply").click()
    at.run()
    assert not at.exception
    cfg = at.session_state["cfg_dict"]
    assert cfg["feature_mode"] == "desat"
    assert cfg["edge_frac"] == pytest.approx(0.65)
    assert cfg["k_band"] == pytest.approx(4.0)
    assert cfg["edge_z"] == pytest.approx(5.0)
