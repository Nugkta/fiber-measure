"""Unit tests for fibrecv.manual_edit (pure manual boundary correction)."""

from __future__ import annotations

import json
from collections import Counter

import numpy as np
import pytest

from fibrecv.band import BandResult
from fibrecv.compute import MeasureResult
from fibrecv.config import CONFIG
from fibrecv.edges import FLAG_BAD_GRAD, FLAG_OK, EdgeResult
from fibrecv.manual_edit import (
    apply_manual_edits,
    corrected_boundary,
    corrected_segments,
    display_to_native,
    empty_edits,
    has_edits,
)
from fibrecv.qc import run_qc

W, H = 300, 300
BAD = slice(120, 161)  # corrupted columns: NaN boundaries + FLAG_BAD_GRAD


def _band() -> BandResult:
    return BandResult(
        mask=np.zeros((H, W), dtype=bool),
        c_fit=np.full(W, 150.0),
        slope=0.0,
        intercept=150.0,
        band_half=50.0,
        x0=5,
        x1=294,
        centroid=np.full(W, 150.0),
        low_confidence=False,
        n_components=1,
    )


def _edges() -> EdgeResult:
    y_top = np.full(W, 100.0)
    y_bot = np.full(W, 200.0)
    flags = np.zeros(W, dtype=np.int64)
    y_top[BAD] = np.nan
    y_bot[BAD] = np.nan
    flags[BAD] = FLAG_BAD_GRAD
    return EdgeResult(
        y_top=y_top,
        y_bot=y_bot,
        diameter=y_bot - y_top,
        amp=np.full(W, 10.0),
        y_core=np.full(W, 150.0),
        flags=flags,
        half_window=80,
    )


def _mr(cfg: CONFIG) -> MeasureResult:
    bnd, edg = _band(), _edges()
    res = run_qc(edg, bnd, cfg)
    diameter_um = np.where(res.valid, res.diameter_raw / cfg.ppu, np.nan)
    span = slice(bnd.x0, bnd.x1 + 1)
    flag_counts = Counter(int(f) for f in res.reason[span])
    meta = {
        "name": "test 1_1_1",
        "coverage": res.coverage,
        "n_valid": int(res.valid.sum()),
        "low_confidence": bool(res.low_confidence),
        "band_mismatch": bool(res.band_mismatch),
        "flag_counts": {str(k): int(v) for k, v in flag_counts.items()},
        "median_diameter_um": float(np.nanmedian(diameter_um)),
        "params": cfg.as_dict(),
    }
    return MeasureResult(
        rgb=None, D=None, bnd=bnd, edg=edg, res=res, diameter_um=diameter_um,
        name="test 1_1_1", group="1_1", replicate=1, meta=meta,
    )


# --------------------------------------------------------------------------- #
# corrected_boundary                                                           #
# --------------------------------------------------------------------------- #
def test_two_anchor_interpolation():
    y = np.full(W, 100.0)
    anchors = [(50.0, 110.0), (80.0, 113.0)]
    y_new, edited = corrected_boundary(y, anchors, ramp=1)
    cols = np.arange(50, 81)
    expect = np.interp(cols, [50.0, 80.0], [110.0, 113.0])
    assert np.allclose(y_new[cols], expect)
    want_mask = np.zeros(W, dtype=bool)
    want_mask[cols] = True
    assert np.array_equal(edited, want_mask)
    # everything outside untouched
    assert np.allclose(y_new[:50], 100.0) and np.allclose(y_new[81:], 100.0)


def test_end_blend_continuity():
    y = np.full(W, 100.0)
    anchors = [(50.0, 110.0), (80.0, 110.0)]
    y_new, _ = corrected_boundary(y, anchors)  # default ramp
    half = min(25, (31 + 1) // 2)
    # at the very ends the correction is diluted to ~1/half of the jump
    assert abs(y_new[50] - 100.0) <= 10.0 / half + 1e-9
    assert abs(y_new[80] - 100.0) <= 10.0 / half + 1e-9
    # middle of the range reaches the full target
    assert y_new[65] == pytest.approx(110.0)
    # outside untouched
    assert y_new[49] == 100.0 and y_new[81] == 100.0


def test_single_anchor_exact():
    y = np.full(W, 100.0)
    y_new, edited = corrected_boundary(y, [(60.0, 95.0)])
    assert y_new[60] == pytest.approx(95.0)
    assert edited.sum() == 1 and edited[60]
    assert y_new[59] == 100.0 and y_new[61] == 100.0


def test_nan_region_replaced():
    y = np.full(W, 100.0)
    y[BAD] = np.nan
    y_new, edited = corrected_boundary(y, [(120.0, 102.0), (160.0, 102.0)])
    assert np.all(np.isfinite(y_new[BAD]))
    assert np.allclose(y_new[BAD], 102.0)
    assert edited[BAD].all()


def test_anchor_x_clamped():
    y = np.full(W, 100.0)
    y_new, edited = corrected_boundary(y, [(250.0, 105.0), (1e6, 105.0)], ramp=1)
    assert edited[250:W].all()
    assert not edited[:250].any()
    assert np.allclose(y_new[250:W], 105.0)


def test_empty_anchors_noop():
    y = np.full(W, 100.0)
    y_new, edited = corrected_boundary(y, [])
    assert np.array_equal(y_new, y) and not edited.any()


# --------------------------------------------------------------------------- #
# corrected_segments (independent anchor sets)                                 #
# --------------------------------------------------------------------------- #
def test_disjoint_segments_leave_gap_untouched():
    y = np.full(W, 100.0)
    segments = [
        [(20.0, 98.0), (60.0, 98.0)],
        [(220.0, 103.0), (260.0, 103.0)],
    ]
    y_new, edited = corrected_segments(y, segments, ramp=1)
    assert np.allclose(y_new[20:61], 98.0)
    assert np.allclose(y_new[220:261], 103.0)
    # the long gap between the two sets is NOT redrawn
    assert np.allclose(y_new[61:220], 100.0)
    assert not edited[61:220].any()
    assert edited[20:61].all() and edited[220:261].all()


def test_overlapping_segments_later_wins():
    y = np.full(W, 100.0)
    segments = [
        [(50.0, 110.0), (100.0, 110.0)],
        [(80.0, 120.0), (130.0, 120.0)],
    ]
    y_new, edited = corrected_segments(y, segments, ramp=1)
    assert np.allclose(y_new[50:80], 110.0)
    assert np.allclose(y_new[80:131], 120.0)
    assert edited[50:131].all()


def test_empty_segments_noop():
    y = np.full(W, 100.0)
    y_new, edited = corrected_segments(y, [], ramp=1)
    assert np.array_equal(y_new, y) and not edited.any()
    y_new, edited = corrected_segments(y, [[]], ramp=1)  # empty set skipped
    assert np.array_equal(y_new, y) and not edited.any()


# --------------------------------------------------------------------------- #
# has_edits / empty_edits                                                      #
# --------------------------------------------------------------------------- #
def test_has_edits():
    assert not has_edits(None)
    assert not has_edits(empty_edits())
    e = empty_edits()
    e["top"].append([])         # an empty set alone is not an edit
    assert not has_edits(e)
    e["top"][0].append((1.0, 2.0))
    assert has_edits(e)
    e2 = empty_edits()
    e2["nudge_bot"] = -1.5
    assert has_edits(e2)


# --------------------------------------------------------------------------- #
# apply_manual_edits                                                           #
# --------------------------------------------------------------------------- #
def test_apply_purity():
    cfg = CONFIG()
    mr = _mr(cfg)
    snap_top = mr.edg.y_top.copy()
    snap_bot = mr.edg.y_bot.copy()
    snap_flags = mr.edg.flags.copy()
    snap_valid = mr.res.valid.copy()
    snap_meta = dict(mr.meta)
    edits = empty_edits()
    edits["top"].append([(120.0, 100.0), (160.0, 100.0)])
    edits["bot"].append([(120.0, 200.0), (160.0, 200.0)])
    edits["nudge_top"] = -1.0
    apply_manual_edits(mr, edits, cfg)
    assert np.array_equal(mr.edg.y_top, snap_top, equal_nan=True)
    assert np.array_equal(mr.edg.y_bot, snap_bot, equal_nan=True)
    assert np.array_equal(mr.edg.flags, snap_flags)
    assert np.array_equal(mr.res.valid, snap_valid)
    assert mr.meta == snap_meta


def test_validity_and_coverage_recovery():
    cfg = CONFIG()
    mr = _mr(cfg)
    cov_before = mr.res.coverage
    assert not mr.res.valid[BAD].any()
    edits = empty_edits()
    edits["top"].append([(120.0, 100.0), (160.0, 100.0)])
    edits["bot"].append([(120.0, 200.0), (160.0, 200.0)])
    new_mr, edited_top, edited_bot = apply_manual_edits(mr, edits, cfg)
    assert new_mr.res.valid[BAD].all()
    assert new_mr.res.coverage > cov_before
    assert np.allclose(new_mr.diameter_um[BAD], 100.0 / cfg.ppu)
    assert edited_top[BAD].all() and edited_bot[BAD].all()
    # meta refreshed + JSON-safe provenance
    assert new_mr.meta["coverage"] == pytest.approx(new_mr.res.coverage)
    assert new_mr.meta["n_valid"] == int(new_mr.res.valid.sum())
    assert "manual_edit" in new_mr.meta
    json.dumps(new_mr.meta["manual_edit"])
    assert new_mr.meta["manual_edit"]["top_sets"] == [[[120.0, 100.0], [160.0, 100.0]]]


def test_one_sided_edit_stays_invalid():
    cfg = CONFIG()
    mr = _mr(cfg)
    edits = empty_edits()
    edits["top"].append([(120.0, 100.0), (160.0, 100.0)])  # bottom stays NaN
    new_mr, _, _ = apply_manual_edits(mr, edits, cfg)
    assert not new_mr.res.valid[BAD].any()
    assert (new_mr.edg.flags[BAD] == FLAG_BAD_GRAD).all()


def test_nudge_shifts_diameter():
    cfg = CONFIG()
    mr = _mr(cfg)
    edits = empty_edits()
    edits["nudge_top"] = -2.0  # top moves up 2 px -> diameter grows by 2
    new_mr, edited_top, _ = apply_manual_edits(mr, edits, cfg)
    good = new_mr.res.valid & mr.res.valid
    assert good.any()
    assert np.allclose(new_mr.res.diameter_raw[good], 102.0)
    assert (new_mr.edg.flags[good] == FLAG_OK).all()
    assert edited_top[good].all()


# --------------------------------------------------------------------------- #
# display_to_native                                                            #
# --------------------------------------------------------------------------- #
def test_display_to_native():
    # asymmetric stretch: sx = 1200/600 = 2, sy = 300/100 = 3, crop offset 400
    x, y = display_to_native(300.0, 50.0, disp_w=600, disp_h=100,
                             crop_w=1200, crop_h=300,
                             y_offset=400, img_w=1200, img_h=1000)
    assert x == pytest.approx(600.0)
    assert y == pytest.approx(550.0)
    # clamping to the native image bounds
    x2, y2 = display_to_native(10_000.0, 10_000.0, disp_w=600, disp_h=100,
                               crop_w=1200, crop_h=300,
                               y_offset=400, img_w=1200, img_h=1000)
    assert x2 == 1200 - 1 and y2 == 1000 - 1
