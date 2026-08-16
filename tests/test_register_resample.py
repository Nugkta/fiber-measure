"""Tests for register.resample_to_grid (shared by aggregate + xsection)."""

from __future__ import annotations

import numpy as np

from fibrecv.config import CONFIG
from fibrecv.register import register_sample, resample_to_grid


def test_resample_to_grid_basic():
    # replicate 0: x 0..4 values 10..14 ; replicate 1 shifted by +0.5
    a0 = (np.arange(5, dtype=float), np.array([10., 11., 12., 13., 14.]))
    a1 = (np.arange(5, dtype=float) + 0.5, np.array([20., 21., 22., 23., 24.]))
    grid, stack = resample_to_grid([a0, a1])
    assert grid.tolist() == [0, 1, 2, 3, 4, 5]
    assert stack.shape == (2, 6)
    # replicate 0 exact on its own grid, NaN at 5
    assert np.allclose(stack[0, :5], a0[1]) and np.isnan(stack[0, 5])
    # replicate 1 midpoints on integer grid, NaN outside its 0.5..4.5 range
    assert np.isnan(stack[1, 0]) and np.isnan(stack[1, 5])
    assert np.allclose(stack[1, 1:5], [20.5, 21.5, 22.5, 23.5])


def test_resample_to_grid_interior_gap_bridged_by_default():
    # legacy behaviour (stage 2): interior NaN runs are linearly bridged
    x = np.arange(11, dtype=float)
    d = x.copy() + 10.0
    d[4:8] = np.nan
    grid, stack = resample_to_grid([(x, d)])
    assert np.allclose(stack[0], x + 10.0)  # gap filled by interpolation


def test_resample_to_grid_max_gap_masks_interior_gap():
    x = np.arange(11, dtype=float)
    d = x.copy() + 10.0
    d[4:8] = np.nan  # finite neighbours at x=3 and x=8: gap of 5 px
    grid, stack = resample_to_grid([(x, d)], max_gap=1.5)
    assert np.isnan(stack[0, 4:8]).all()  # strictly inside the gap: NaN
    assert np.allclose(stack[0, :4], d[:4])  # edges + finite samples intact
    assert np.allclose(stack[0, 8:], d[8:])


def test_resample_to_grid_max_gap_allows_subpixel_spacing():
    # consecutive finite samples 1 px apart (integer columns + sub-pixel
    # shift) must still interpolate normally under max_gap=1.5
    x = np.arange(6, dtype=float) + 0.4
    d = np.array([10., 11., 12., 13., 14., 15.])
    grid, stack = resample_to_grid([(x, d)], max_gap=1.5)
    inside = (grid >= x[0]) & (grid <= x[-1])
    assert np.isfinite(stack[0, inside]).all()


def test_resample_to_grid_all_nan_degrades():
    aligned = [(np.arange(4, dtype=float), np.full(4, np.nan)),
               (np.arange(4, dtype=float), np.full(4, np.nan))]
    grid, stack = resample_to_grid(aligned)
    assert grid.size == 0
    assert stack.shape == (2, 0)


def test_register_sample_all_nan_profiles_no_crash():
    profiles = [
        {"x": np.arange(10), "diameter_px_raw": np.full(10, np.nan),
         "diameter_px_smooth": np.full(10, np.nan), "valid": np.zeros(10, bool),
         "replicate": r}
        for r in (1, 2)
    ]
    table, shifts, summary = register_sample(profiles, CONFIG())
    assert table["x_aligned_px"].size == 0
    assert np.isnan(summary["mean_um"])
    assert len(shifts) == 2
