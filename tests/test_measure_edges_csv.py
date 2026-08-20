"""Per-image CSV must carry per-column edge positions (additive columns)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from fibrecv import run_aggregate
from fibrecv.compute import compute_measurement
from fibrecv.config import CONFIG
from fibrecv.measure import write_measurement

from test_features_modes import _pink_fibre

EXPECTED_LEADING = ["x_px", "diameter_px_raw", "diameter_px_smooth",
                    "diameter_um", "valid", "interpolated"]


def _measure_to(tmp_path, name: str, seed: int) -> pd.DataFrame:
    rgb = _pink_fibre(seed=seed)
    cfg = CONFIG()
    mr = compute_measurement(rgb, cfg, name)
    write_measurement(rgb, mr, cfg, tmp_path)
    return pd.read_csv(tmp_path / "per_image" / "csv" / f"{name}_profile.csv")


def test_csv_has_edge_columns_and_existing_schema(tmp_path):
    df = _measure_to(tmp_path, "masp2 3_1_1", seed=0)
    # existing columns unchanged, in order; new ones appended
    assert list(df.columns)[:6] == EXPECTED_LEADING
    assert "y_top_px" in df.columns and "y_bot_px" in df.columns

    valid = df["valid"].astype(bool).to_numpy()
    y_top = df["y_top_px"].to_numpy()
    y_bot = df["y_bot_px"].to_numpy()
    # valid columns: finite edges, bottom below top, consistent with diameter
    assert np.isfinite(y_top[valid]).all() and np.isfinite(y_bot[valid]).all()
    assert (y_bot[valid] > y_top[valid]).all()
    # invalid columns: NaN edges
    if (~valid).any():
        assert np.isnan(y_top[~valid]).all() and np.isnan(y_bot[~valid]).all()


def test_run_aggregate_accepts_csv_with_edge_columns(tmp_path):
    for rep, seed in [(1, 0), (2, 1)]:
        _measure_to(tmp_path, f"masp2 3_1_{rep}", seed=seed)
    rc = run_aggregate.main(["--out", str(tmp_path), "--groups", "3_1"])
    assert rc == 0
    assert (tmp_path / "summary" / "master_summary.csv").exists()
