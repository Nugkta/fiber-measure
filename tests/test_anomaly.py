"""Unit tests for image-level and group-level anomaly detection."""

from __future__ import annotations

import json
from dataclasses import replace

import numpy as np

from fibrecv.anomaly import (
    AnomalyResult,
    detect_image_anomalies,
    detect_replicate_outliers,
    exclusion_reason,
)
from fibrecv.config import CONFIG


def _straight_fibre(W=300, x0=5, x1=294, diameter=100.0, center=150.0):
    """Flat horizontal fibre with every span column valid."""
    y_top = np.full(W, center - diameter / 2.0)
    y_bot = np.full(W, center + diameter / 2.0)
    smooth = np.full(W, np.nan)
    smooth[x0:x1 + 1] = diameter
    valid = np.zeros(W, dtype=bool)
    valid[x0:x1 + 1] = True
    return y_top, y_bot, smooth, valid, x0, x1


def _detect(y_top, y_bot, smooth, valid, x0, x1, slope=0.0, cfg=None):
    return detect_image_anomalies(
        y_top, y_bot, smooth, valid, x0, x1, slope, cfg or CONFIG()
    )


# --- config defaults -------------------------------------------------------

def test_config_defaults():
    cfg = CONFIG()
    assert cfg.jump_thresh_px == 10.0
    assert cfg.gap_frac == 0.10
    assert cfg.step_frac == 0.05
    assert cfg.step_window_px == 100
    assert cfg.rep_dev_frac == 0.25
    assert cfg.anomaly_exclude is False


# --- edge_jump -------------------------------------------------------------

def test_jump_on_top_edge_triggers():
    y_top, y_bot, smooth, valid, x0, x1 = _straight_fibre()
    y_top[150:] -= 15.0  # 15 px jump at column 150
    res = _detect(y_top, y_bot, smooth, valid, x0, x1)
    assert "edge_jump" in res.flags
    assert 150 in res.jump_cols
    assert res.max_jump_px >= 15.0
    assert res.n_jumps >= 1


def test_small_jump_is_clean():
    y_top, y_bot, smooth, valid, x0, x1 = _straight_fibre()
    y_top[150:] -= 5.0  # below jump_thresh_px=10
    res = _detect(y_top, y_bot, smooth, valid, x0, x1)
    assert "edge_jump" not in res.flags
    assert res.n_jumps == 0
    assert res.jump_cols == []


def test_bottom_edge_is_checked_too():
    y_top, y_bot, smooth, valid, x0, x1 = _straight_fibre()
    y_bot[200:] += 12.0
    res = _detect(y_top, y_bot, smooth, valid, x0, x1)
    assert "edge_jump" in res.flags
    assert 200 in res.jump_cols


def test_tilted_fibre_with_gap_not_false_flagged():
    """Diffing across a 41-col invalid gap on a slope-0.5 fibre gives a raw
    Δy of ~21 px; the slope-detrended residual must be ~0."""
    y_top, y_bot, smooth, valid, x0, x1 = _straight_fibre()
    slope = 0.5
    x = np.arange(y_top.size)
    y_top = y_top + slope * x
    y_bot = y_bot + slope * x
    valid[100:141] = False
    res = _detect(y_top, y_bot, smooth, valid, x0, x1, slope=slope)
    assert "edge_jump" not in res.flags
    assert res.max_jump_px < 1.0


def test_all_invalid_does_not_crash():
    y_top, y_bot, smooth, valid, x0, x1 = _straight_fibre()
    valid[:] = False
    smooth[:] = np.nan
    res = _detect(y_top, y_bot, smooth, valid, x0, x1)
    assert res.longest_gap_frac == 1.0
    assert "large_gap" in res.flags
    assert res.n_jumps == 0


# --- large_gap -------------------------------------------------------------

def test_gap_40_of_290_triggers():
    y_top, y_bot, smooth, valid, x0, x1 = _straight_fibre()
    valid[100:140] = False  # 40 / 290 span cols = 13.8% > 10%
    res = _detect(y_top, y_bot, smooth, valid, x0, x1)
    assert "large_gap" in res.flags
    assert res.longest_gap_frac > 0.10
    assert res.gap_start_col == 100


def test_gap_20_of_290_is_clean():
    y_top, y_bot, smooth, valid, x0, x1 = _straight_fibre()
    valid[100:120] = False  # 20 / 290 = 6.9%
    res = _detect(y_top, y_bot, smooth, valid, x0, x1)
    assert "large_gap" not in res.flags
    assert 0.0 < res.longest_gap_frac < 0.10


def test_no_gap_reports_zero():
    y_top, y_bot, smooth, valid, x0, x1 = _straight_fibre()
    res = _detect(y_top, y_bot, smooth, valid, x0, x1)
    assert res.longest_gap_frac == 0.0
    assert res.gap_start_col is None


# --- diameter_step ---------------------------------------------------------

def _stepped_fibre(W=800, step_at=400, d0=100.0, d1=112.0):
    y_top = np.full(W, 100.0)
    y_bot = y_top + d0
    y_bot[step_at:] = 100.0 + d1
    smooth = np.where(np.arange(W) < step_at, d0, d1)
    valid = np.ones(W, dtype=bool)
    return y_top, y_bot, smooth.astype(float), valid, 0, W - 1


def test_diameter_step_12pct_triggers():
    y_top, y_bot, smooth, valid, x0, x1 = _stepped_fibre()
    res = _detect(y_top, y_bot, smooth, valid, x0, x1)
    assert "diameter_step" in res.flags
    assert res.step_frac > 0.05
    assert abs(res.step_col - 400) <= 10


def test_linear_taper_is_clean():
    y_top, y_bot, smooth, valid, x0, x1 = _stepped_fibre()
    smooth = np.linspace(100.0, 110.0, smooth.size)  # 10% taper over 800 px
    res = _detect(y_top, y_top + smooth, smooth, valid, x0, x1)
    assert "diameter_step" not in res.flags


def test_short_span_skips_step_detector():
    y_top, y_bot, smooth, valid, x0, x1 = _straight_fibre(W=300, x0=5, x1=154)
    res = _detect(y_top, y_bot, smooth, valid, x0, x1)  # span 150 < 2*100
    assert res.step_frac is None or np.isnan(res.step_frac)
    assert "diameter_step" not in res.flags


# --- replicate_outlier -----------------------------------------------------

def test_replicate_outlier_detected():
    cfg = CONFIG()
    devs, outliers = detect_replicate_outliers({1: 40.0, 2: 41.0, 3: 60.0}, cfg)
    assert outliers == {3}
    assert abs(devs[3] - 19.0 / 41.0) < 1e-9
    assert devs[1] < cfg.rep_dev_frac


def test_replicate_outlier_needs_three():
    cfg = CONFIG()
    devs, outliers = detect_replicate_outliers({1: 40.0, 2: 60.0}, cfg)
    assert outliers == set()
    assert devs == {}


def test_replicate_outlier_drops_none_medians():
    cfg = CONFIG()
    devs, outliers = detect_replicate_outliers(
        {1: 40.0, 2: 41.0, 3: None, 4: 60.0}, cfg
    )
    assert outliers == {4}
    assert 3 not in devs


# --- exclusion_reason ------------------------------------------------------

def test_exclusion_priority_band_mismatch_first():
    cfg = replace(CONFIG(), anomaly_exclude=True)
    reason = exclusion_reason(True, 0.2, ["large_gap"], cfg)
    assert reason == "band_mismatch"


def test_exclusion_coverage_wording():
    cfg = CONFIG()
    assert exclusion_reason(False, 0.4, [], cfg) == "coverage 40% < 50%"
    assert exclusion_reason(False, 0.9, [], cfg) is None
    assert exclusion_reason(False, None, [], cfg) is None


def test_exclusion_anomaly_only_when_enabled():
    cfg = CONFIG()
    assert exclusion_reason(False, 0.9, ["large_gap"], cfg) is None
    cfg_ex = replace(cfg, anomaly_exclude=True)
    assert exclusion_reason(False, 0.9, ["large_gap"], cfg_ex) == "anomaly: large_gap"
    assert exclusion_reason(
        False, 0.9, ["edge_jump", "large_gap"], cfg_ex
    ) == "anomaly: edge_jump, large_gap"


def test_replicate_outlier_never_excludes():
    cfg = replace(CONFIG(), anomaly_exclude=True)
    assert exclusion_reason(False, 0.9, ["replicate_outlier"], cfg) is None


# --- meta JSON (compute + manual_edit) -------------------------------------

_META_KEYS = {
    "flags", "jump_cols", "max_jump_px", "n_jumps",
    "longest_gap_frac", "gap_start_col", "step_frac", "step_col",
}


def test_compute_meta_carries_anomaly_dict():
    from fibrecv.compute import compute_measurement

    rng = np.random.default_rng(0)
    rgb = np.empty((300, 800, 3), dtype=np.float32)
    rgb[:] = (0.95, 0.45, 0.65)          # saturated pink background
    rgb[130:170] = (0.92, 0.90, 0.91)    # pale desaturated fibre band
    rgb += rng.normal(0.0, 0.01, rgb.shape).astype(np.float32)
    rgb = np.clip(rgb, 0.0, 1.0)

    mr = compute_measurement(rgb, CONFIG(), name="test 1_1_1")
    assert set(mr.meta["anomaly"].keys()) == _META_KEYS
    json.dumps(mr.meta)  # whole meta must stay JSON-safe


def _mr_with_jump(cfg):
    """MeasureResult whose top edge has a 15-px 20-col dent (edge_jump)."""
    from fibrecv.band import BandResult
    from fibrecv.compute import MeasureResult
    from fibrecv.edges import EdgeResult
    from fibrecv.qc import run_qc

    W = 300
    y_top = np.full(W, 100.0)
    y_bot = np.full(W, 200.0)
    y_top[140:160] = 85.0  # 15-px jump down at 140 and back up at 160
    edg = EdgeResult(
        y_top=y_top, y_bot=y_bot, diameter=y_bot - y_top,
        amp=np.full(W, 10.0), y_core=np.full(W, 150.0),
        flags=np.zeros(W, dtype=np.int64), half_window=80,
    )
    bnd = BandResult(
        mask=np.zeros((300, W), dtype=bool), c_fit=np.full(W, 150.0),
        slope=0.0, intercept=150.0, band_half=50.0, x0=5, x1=294,
        centroid=np.full(W, 150.0), low_confidence=False, n_components=1,
    )
    res = run_qc(edg, bnd, cfg)
    diameter_um = np.where(res.valid, res.diameter_raw / cfg.ppu, np.nan)
    meta = {"name": "test 1_1_1", "coverage": res.coverage,
            "anomaly": res.anomaly.as_dict()}
    return MeasureResult(
        rgb=None, D=None, bnd=bnd, edg=edg, res=res, diameter_um=diameter_um,
        name="test 1_1_1", group="1_1", replicate=1, meta=meta,
    )


def test_manual_edit_clears_anomaly_in_res_and_meta():
    from fibrecv.manual_edit import apply_manual_edits, empty_edits

    cfg = CONFIG()
    mr = _mr_with_jump(cfg)
    assert "edge_jump" in mr.res.anomaly.flags
    assert "edge_jump" in mr.meta["anomaly"]["flags"]

    edits = empty_edits()
    # redraw the dented range back to 100; anchors extend past the dent so
    # the ramp's diluted ends land on already-correct columns
    edits["top"] = [[(120.0, 100.0), (180.0, 100.0)]]
    new_mr, _, _ = apply_manual_edits(mr, edits, cfg)
    assert "edge_jump" not in new_mr.res.anomaly.flags
    assert "edge_jump" not in new_mr.meta["anomaly"]["flags"]
    json.dumps(new_mr.meta)


# --- as_dict ---------------------------------------------------------------

def test_as_dict_json_safe_and_capped():
    y_top, y_bot, smooth, valid, x0, x1 = _straight_fibre(W=600, x0=0, x1=599)
    # 30 separate jumps -> jump_cols must be capped at 20 in the dict
    for i, c in enumerate(range(20, 580, 18)):
        y_top[c] += 15.0 if i % 2 == 0 else -15.0
    res = _detect(y_top, y_bot, smooth, valid, x0, x1)
    d = res.as_dict()
    assert set(d.keys()) == {
        "flags", "jump_cols", "max_jump_px", "n_jumps",
        "longest_gap_frac", "gap_start_col", "step_frac", "step_col",
    }
    assert len(d["jump_cols"]) <= 20
    assert d["n_jumps"] > 20  # true count survives the cap
    json.dumps(d)  # must not raise
    assert all(isinstance(c, int) for c in d["jump_cols"])
    assert isinstance(d["max_jump_px"], float)


def _qc_fixture(diameter_px=100.0, band_half=50.0, slope=0.0, W=300):
    """EdgeResult/BandResult pair mirroring tests/test_edges_wall.py."""
    from fibrecv.band import BandResult
    from fibrecv.edges import EdgeResult

    x = np.arange(W, dtype=float)
    center = 150.0 + slope * x
    y_top = center - diameter_px / 2.0
    y_bot = center + diameter_px / 2.0
    edg = EdgeResult(
        y_top=y_top, y_bot=y_bot, diameter=y_bot - y_top,
        amp=np.full(W, 10.0), y_core=center.copy(),
        flags=np.zeros(W, dtype=np.int64), half_window=120,
    )
    bnd = BandResult(
        mask=np.zeros((300, W), dtype=bool), c_fit=center.copy(),
        slope=slope, intercept=150.0, band_half=band_half, x0=5, x1=294,
        centroid=center.copy(), low_confidence=False, n_components=1,
    )
    return edg, bnd


def test_qc_clean_fibre_has_no_anomalies():
    from fibrecv.qc import run_qc

    edg, bnd = _qc_fixture()
    res = run_qc(edg, bnd, CONFIG())
    assert res.anomaly.flags == []


def test_qc_injected_jump_flagged_but_stays_valid():
    from fibrecv.qc import run_qc

    edg, bnd = _qc_fixture()
    edg.y_top[150:160] -= 15.0
    edg.diameter[150:160] += 15.0
    res = run_qc(edg, bnd, CONFIG())
    assert "edge_jump" in res.anomaly.flags
    assert 150 in res.anomaly.jump_cols
    # advisory: the jump must not invalidate columns around it
    assert res.valid[148:162].all()


def test_qc_tilted_gap_uses_refit_slope():
    """A 41-col gap on a slope-0.5 fibre: the detrend must use the refit
    centerline slope, so the across-gap step is not a false edge_jump."""
    from fibrecv.qc import run_qc

    edg, bnd = _qc_fixture(slope=0.5)
    edg.flags[100:141] = 1  # invalidate a 41-col run
    res = run_qc(edg, bnd, CONFIG())
    assert "edge_jump" not in res.anomaly.flags
    assert "large_gap" in res.anomaly.flags


def test_as_dict_nan_becomes_none():
    y_top, y_bot, smooth, valid, x0, x1 = _straight_fibre(W=300, x0=5, x1=154)
    res = _detect(y_top, y_bot, smooth, valid, x0, x1)  # span too short for step
    d = res.as_dict()
    assert d["step_frac"] is None
    assert d["step_col"] is None
    text = json.dumps(d)
    assert "NaN" not in text
