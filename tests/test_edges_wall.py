"""Unit tests for gap-tolerant wall finding and the wide-band QC guard."""

from __future__ import annotations

import numpy as np

from fibrecv.band import BandResult
from fibrecv.config import CONFIG
from fibrecv.edges import EdgeResult, _find_wall
from fibrecv.qc import run_qc


def _plateaued_ramp():
    """Profile with a defocus-soft outer wall broken by two flat plateaus.

    Background 0 for indices 0..59; three rising fragments of 1.8 z each
    (below rise_min=2.0) separated by 11-px plateaus, total rise 5.4 z, core
    plateau from ~118 on. The strict contiguous-run rule rejects every
    fragment for rise; bridging the plateaus heals the run.
    """
    d = np.zeros(160)
    level = 0.0
    for i in range(60, 118):
        if 72 <= i <= 82 or 95 <= i <= 105:
            pass                  # plateau: no rise
        else:
            level += 0.15         # rising fragments: 12 px x 0.15 = 1.8 z each
        d[i] = level
    d[118:] = level
    g = np.gradient(d)
    return d, g


def test_legacy_gap0_fragments_the_ramp():
    d, g = _plateaued_ramp()
    cfg = CONFIG()
    wall = _find_wall(d, g, 0, 130, +1, cfg, gap=0)
    # every contiguous fragment rises < rise_min -> no wall at all
    assert wall is None


def test_gap_bridges_plateaus():
    d, g = _plateaued_ramp()
    cfg = CONFIG()
    wall = _find_wall(d, g, 0, 130, +1, cfg, gap=12)
    assert wall is not None
    i_outer, i_inner = wall
    assert i_outer <= 65            # run starts at the true outer wall
    assert d[i_inner] - d[i_outer] >= cfg.rise_min


def test_gap_does_not_bridge_falling_profile():
    """A falling stretch (not a plateau) must still terminate the run."""
    d = np.zeros(160)
    d[60:80] = np.linspace(0, 1.8, 20)      # small outer blip, rise < rise_min
    d[80:100] = np.linspace(1.8, 0.2, 20)   # falls back down
    d[100:130] = np.linspace(0.2, 6.0, 30)  # true wall further in
    d[130:] = 6.0
    g = np.gradient(d)
    cfg = CONFIG()
    wall = _find_wall(d, g, 0, 140, +1, cfg, gap=12)
    assert wall is not None
    i_outer, _ = wall
    assert i_outer >= 95            # blip rejected: run did not bridge the fall


def _qc_fixture(diameter_px: float, band_half: float):
    W = 300
    y_top = np.full(W, 150.0 - diameter_px / 2.0)
    y_bot = np.full(W, 150.0 + diameter_px / 2.0)
    edg = EdgeResult(
        y_top=y_top, y_bot=y_bot, diameter=y_bot - y_top,
        amp=np.full(W, 10.0), y_core=np.full(W, 150.0),
        flags=np.zeros(W, dtype=np.int64), half_window=120,
    )
    bnd = BandResult(
        mask=np.zeros((300, W), dtype=bool), c_fit=np.full(W, 150.0),
        slope=0.0, intercept=150.0, band_half=band_half, x0=5, x1=294,
        centroid=np.full(W, 150.0), low_confidence=False, n_components=1,
    )
    return edg, bnd


def test_qc_flags_too_wide_vs_band():
    cfg = CONFIG()
    edg, bnd = _qc_fixture(diameter_px=200.0, band_half=50.0)  # ratio 2.0
    res = run_qc(edg, bnd, cfg)
    assert res.band_mismatch
    assert res.low_confidence


def test_qc_accepts_consistent_diameter():
    cfg = CONFIG()
    edg, bnd = _qc_fixture(diameter_px=100.0, band_half=50.0)  # ratio 1.0
    res = run_qc(edg, bnd, cfg)
    assert not res.band_mismatch
