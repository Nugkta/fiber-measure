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


@pytest.fixture
def image_folder(tmp_path: Path) -> Path:
    for repn in (1, 2):
        arr = (_synthetic_fibre(seed=repn) * 255).astype(np.uint8)
        iio.imwrite(tmp_path / f"test 1_1_{repn}.png", arr)
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
