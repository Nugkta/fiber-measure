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
