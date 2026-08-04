"""Aggregation-level anomaly outputs: per_image_summary.csv + exclusion flags."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fibrecv.config import CONFIG
from fibrecv.run_aggregate import main

PPU = CONFIG().ppu

EXPECTED_COLS = [
    "name", "group", "replicate", "median_diameter_um", "coverage",
    "anomaly_flags", "max_jump_px", "longest_gap_frac", "step_frac",
    "rep_dev_frac", "excluded", "excluded_reason",
]


def _clean_anomaly() -> dict:
    return {"flags": [], "jump_cols": [], "max_jump_px": 0.0, "n_jumps": 0,
            "longest_gap_frac": 0.0, "gap_start_col": None,
            "step_frac": 0.01, "step_col": 200}


def _write_replicate(out: Path, rep: int, diameter_px: float,
                     anomaly: dict | None, drop_anomaly_key: bool = False) -> None:
    """Write the run_measure output contract for one replicate directly."""
    csv_dir = out / "per_image" / "csv"
    meta_dir = out / "per_image" / "diagnostics"
    csv_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)

    x = np.arange(400, dtype=float)
    wiggle = 2.0 * np.sin(x / 30.0)  # gives registration a correlation signal
    d = np.full(400, diameter_px) + wiggle
    base = f"test 1_1_{rep}"
    pd.DataFrame({
        "x_px": x,
        "diameter_px_raw": d,
        "diameter_px_smooth": d,
        "valid": np.ones(400, dtype=bool),
    }).to_csv(csv_dir / f"{base}_profile.csv", index=False)

    meta = {
        "name": base, "group": "1_1", "replicate": rep,
        "coverage": 0.95, "low_confidence": False, "band_mismatch": False,
        "median_diameter_um": diameter_px / PPU,
    }
    if not drop_anomaly_key:
        meta["anomaly"] = anomaly if anomaly is not None else _clean_anomaly()
    with open(meta_dir / f"{base}_meta.json", "w") as fh:
        json.dump(meta, fh)


@pytest.fixture
def out_root(tmp_path: Path) -> Path:
    """3 replicates: rep1 clean, rep2 large_gap, rep3 60%-deviant median."""
    gap = _clean_anomaly()
    gap.update({"flags": ["large_gap"], "longest_gap_frac": 0.2,
                "gap_start_col": 120})
    _write_replicate(tmp_path, 1, 100.0, None)
    _write_replicate(tmp_path, 2, 100.0, gap)
    _write_replicate(tmp_path, 3, 160.0, None)
    return tmp_path


def _read_summary(out: Path) -> pd.DataFrame:
    df = pd.read_csv(out / "summary" / "per_image_summary.csv")
    df["anomaly_flags"] = df["anomaly_flags"].fillna("")
    df["excluded_reason"] = df["excluded_reason"].fillna("")
    return df.set_index("name")


def test_advisory_mode_keeps_everything(out_root: Path):
    assert main(["--out", str(out_root), "--all"]) == 0
    df = _read_summary(out_root)
    assert list(df.reset_index().columns) == EXPECTED_COLS
    assert len(df) == 3
    assert not df["excluded"].any()

    assert df.loc["test 1_1_2", "anomaly_flags"] == "large_gap"
    assert df.loc["test 1_1_3", "anomaly_flags"] == "replicate_outlier"
    assert abs(df.loc["test 1_1_3", "rep_dev_frac"] - 0.6) < 0.01
    assert df.loc["test 1_1_1", "rep_dev_frac"] < 0.05

    master = pd.read_csv(out_root / "summary" / "master_summary.csv")
    assert int(master.loc[0, "n_replicates_used"]) == 3


def test_anomaly_exclude_drops_flagged_image(out_root: Path):
    assert main(["--out", str(out_root), "--all", "--anomaly-exclude"]) == 0
    df = _read_summary(out_root)

    assert bool(df.loc["test 1_1_2", "excluded"])
    assert df.loc["test 1_1_2", "excluded_reason"] == "anomaly: large_gap"
    # replicate_outlier never excludes; deviant rep 3 stays in
    assert not bool(df.loc["test 1_1_3", "excluded"])
    # excluded rep 2 still participates in the group-median computation
    assert abs(df.loc["test 1_1_3", "rep_dev_frac"] - 0.6) < 0.01

    master = pd.read_csv(out_root / "summary" / "master_summary.csv")
    assert int(master.loc[0, "n_replicates_used"]) == 2


def test_summary_written_even_when_everything_excluded(tmp_path: Path):
    """All replicates excluded -> no groups to register (exit 1), but the
    per-image summary must still exist: it is the audit of WHY."""
    gap = _clean_anomaly()
    gap.update({"flags": ["large_gap"], "longest_gap_frac": 0.3})
    for rep in (1, 2, 3):
        _write_replicate(tmp_path, rep, 100.0, dict(gap))
    assert main(["--out", str(tmp_path), "--all", "--anomaly-exclude"]) == 1
    df = _read_summary(tmp_path)
    assert len(df) == 3
    assert df["excluded"].all()
    assert (df["excluded_reason"] == "anomaly: large_gap").all()
    assert not (tmp_path / "summary" / "master_summary.csv").exists()


def test_old_meta_without_anomaly_key(tmp_path: Path):
    for rep in (1, 2, 3):
        _write_replicate(tmp_path, rep, 100.0, None, drop_anomaly_key=True)
    assert main(["--out", str(tmp_path), "--all"]) == 0
    df = _read_summary(tmp_path)
    assert len(df) == 3
    assert (df["anomaly_flags"] == "").all()
    assert not df["excluded"].any()
