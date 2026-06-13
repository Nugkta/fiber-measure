"""Tests for tensile (stress-strain) ingestion and metric computation.

The synthetic curves use a perfectly linear-elastic-then-fracture ("triangular")
stress-strain response so the expected Young's modulus, tensile strength and
toughness have closed forms we can assert against, independent of the
implementation.
"""

import numpy as np
import pandas as pd
import pytest

from fibrecv.config import CONFIG
from fibrecv.tensile import (
    build_matrix, compute_tensile, discover_tensile, discover_tensile_files,
    parse_tensile_name, read_trace)


class _FakeUpload:
    """Minimal stand-in for a Streamlit UploadedFile (``.name`` + ``getvalue``)."""

    def __init__(self, name, data=b""):
        self.name = name
        self._data = data

    def getvalue(self):
        return self._data


# --------------------------------------------------------------------------- #
# Helpers                                                                      #
# --------------------------------------------------------------------------- #
def _write_ta_csv(path, disp_mm, load_n, force_col="Load 3"):
    """Write a minimal TA-Instruments-WinTest-style CSV (metadata + data block)."""
    lines = [
        r"C:\TA Instruments\WinTest\TestData\fake.CSV",
        "",
        "Created ,10/06/2026,12:29:56",
        "",
        "Channels",
        "Disp,Load,Load 2,Load 3,",
        "Units",
        "mm,N,N,N,",
        "",
        "UserID:, username",
        "",
        f'"Points","Elapsed Time","Disp","{force_col}",',
        '"","Sec","mm","N",',
    ]
    # data rows; Points counter resets every 100 (mimics the real per-scan reset)
    for i, (d, ld) in enumerate(zip(disp_mm, load_n)):
        lines.append(f"{i % 100 + 1},{i * 0.01:.4f},{d:.5f},{ld:.6f},")
    path.write_text("\r\n".join(lines), encoding="latin-1")
    return path


def _triangular(d_um=10.0, L0_mm=10.0, E_pa=1.0e9, strain_b=0.05, n=51):
    """Linear-elastic rise to (strain_b, E*strain_b) then a sharp fracture drop.

    Returns (df, area_m2, fmax_n) with closed-form metrics:
        modulus = E_pa, strength = E*strain_b, toughness = 0.5*E*strain_b**2.
    """
    area = np.pi * (d_um * 1e-6 / 2.0) ** 2
    disp_rise = np.linspace(0.0, strain_b * L0_mm, n)
    load_rise = E_pa * (disp_rise / L0_mm) * area
    fmax = load_rise[-1]
    disp_drop = strain_b * L0_mm + np.array([0.005, 0.010, 0.015])
    load_drop = np.array([0.10, 0.02, 0.0]) * fmax   # all < break_drop_frac*Fmax
    df = pd.DataFrame({
        "disp_mm": np.concatenate([disp_rise, disp_drop]),
        "load_n": np.concatenate([load_rise, load_drop]),
    })
    return df, area, fmax


# --------------------------------------------------------------------------- #
# Filename parsing / discovery                                                 #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name,group,kind", [
    ("masp2 10_1 06102026 122956_std.CSV", "10_1", "std"),
    ("masp2 10_1 06102026 122956_rdr.CSV", "10_1", "rdr"),
    ("masp2 3_5 06112026 031644_std.CSV", "3_5", "std"),
    ("masp2 7_10 06112026 021020_std.CSV", "7_10", "std"),
])
def test_parse_tensile_name(name, group, kind):
    assert parse_tensile_name(name) == (group, kind)


@pytest.mark.parametrize("name", ["notes.txt", "calibration only.CSV"])
def test_parse_tensile_name_rejects(name):
    with pytest.raises(ValueError):
        parse_tensile_name(name)


def test_discover_tensile_prefers_std(tmp_path):
    df, _, _ = _triangular()
    _write_ta_csv(tmp_path / "masp2 3_1 06112026 030624_std.CSV",
                  df["disp_mm"], df["load_n"])
    _write_ta_csv(tmp_path / "masp2 3_1 06112026 030624_rdr.CSV",
                  df["disp_mm"], df["load_n"])
    _write_ta_csv(tmp_path / "masp2 4_2 06112026 022059_rdr.CSV",
                  df["disp_mm"], df["load_n"])
    found = discover_tensile(tmp_path)
    assert set(found) == {"3_1", "4_2"}
    assert found["3_1"].name.endswith("_std.CSV")   # std wins over rdr
    assert found["4_2"].name.endswith("_rdr.CSV")   # rdr-only still discovered


# --------------------------------------------------------------------------- #
# Trace reading                                                                #
# --------------------------------------------------------------------------- #
def test_read_trace_parses_data_block(tmp_path):
    disp = np.array([0.015, 0.016, 0.017, 0.020])
    load = np.array([0.024, 0.023, 0.026, 0.050])
    p = _write_ta_csv(tmp_path / "masp2 5_1 06112026 024206_std.CSV", disp, load)
    df = read_trace(p)
    assert list(df.columns) == ["disp_mm", "load_n"]
    assert len(df) == 4
    np.testing.assert_allclose(df["disp_mm"].to_numpy(), disp, atol=1e-9)
    np.testing.assert_allclose(df["load_n"].to_numpy(), load, atol=1e-9)


def test_read_trace_from_bytes_and_upload(tmp_path):
    disp = np.array([0.01, 0.02, 0.03])
    load = np.array([0.10, 0.20, 0.30])
    p = _write_ta_csv(tmp_path / "masp2 6_1 06112026 010000_std.CSV", disp, load)
    raw = p.read_bytes()
    # raw bytes (filename required to pick the reader)
    df_b = read_trace(raw, name=p.name)
    np.testing.assert_allclose(df_b["load_n"].to_numpy(), load, atol=1e-9)
    # file-like upload (getvalue + .name), as drag-and-drop would provide
    df_u = read_trace(_FakeUpload(p.name, raw))
    np.testing.assert_allclose(df_u["disp_mm"].to_numpy(), disp, atol=1e-9)


def test_discover_tensile_files_prefers_std():
    ups = [
        _FakeUpload("masp2 3_1 06112026 030624_std.CSV"),
        _FakeUpload("masp2 3_1 06112026 030624_rdr.CSV"),
        _FakeUpload("masp2 4_2 06112026 022059_rdr.CSV"),
        _FakeUpload("notes.txt"),                 # unparseable -> skipped
    ]
    found = discover_tensile_files(ups)
    assert set(found) == {"3_1", "4_2"}
    assert found["3_1"].name.endswith("_std.CSV")
    assert found["4_2"].name.endswith("_rdr.CSV")


def test_build_matrix_accepts_uploads(tmp_path):
    df, _, _ = _triangular()
    raw = _write_ta_csv(tmp_path / "masp2 3_1 06112026 030624_std.CSV",
                        df["disp_mm"], df["load_n"]).read_bytes()
    tensile = {"3_1": _FakeUpload("masp2 3_1 06112026 030624_std.CSV", raw)}
    mat = build_matrix({"3_1": 10.0}, tensile, CONFIG())
    row = mat.iloc[0]
    assert row["group"] == "3_1" and row["flag"] == ""
    assert np.isfinite(row["youngs_modulus_GPa"])


def test_read_trace_force_column_fallback(tmp_path):
    disp = np.array([0.01, 0.02, 0.03])
    load = np.array([0.1, 0.2, 0.3])
    p = _write_ta_csv(tmp_path / "masp2 6_1 06112026 010000_std.CSV", disp, load,
                      force_col="Load")          # no "Load 3" -> fall back to "Load"
    df = read_trace(p)
    np.testing.assert_allclose(df["load_n"].to_numpy(), load, atol=1e-9)


# --------------------------------------------------------------------------- #
# Metric computation                                                           #
# --------------------------------------------------------------------------- #
def test_compute_tensile_recovers_known_metrics():
    df, area, fmax = _triangular(d_um=10.0, L0_mm=10.0, E_pa=1.0e9, strain_b=0.05)
    res = compute_tensile(df, diameter_um=10.0, gauge_length_mm=10.0, cfg=CONFIG())

    assert res.fmax_n == pytest.approx(fmax, rel=1e-6)
    assert res.area_m2 == pytest.approx(area, rel=1e-9)
    assert res.tensile_strength_pa == pytest.approx(0.05 * 1.0e9, rel=1e-3)  # E*strain_b
    assert res.youngs_modulus_pa == pytest.approx(1.0e9, rel=0.02)
    assert res.toughness_j_m3 == pytest.approx(0.5 * 1.0e9 * 0.05 ** 2, rel=0.10)
    assert res.strain_at_break == pytest.approx(0.0505, abs=2e-3)
    assert res.extension_at_break_mm == pytest.approx(0.505, abs=1e-2)
    # fracture marked just past the peak, not at the noisy tail
    assert df["disp_mm"].iloc[res.break_index] == pytest.approx(0.505, abs=1e-2)


def test_compute_tensile_without_diameter_gives_force_only():
    df, _, fmax = _triangular()
    res = compute_tensile(df, diameter_um=None, gauge_length_mm=10.0, cfg=CONFIG())
    assert res.fmax_n == pytest.approx(fmax, rel=1e-6)        # force metrics still work
    assert res.extension_at_break_mm == pytest.approx(0.505, abs=1e-2)
    assert np.isnan(res.tensile_strength_pa)                  # need area for stress
    assert np.isnan(res.youngs_modulus_pa)
    assert np.isnan(res.toughness_j_m3)


# --------------------------------------------------------------------------- #
# Fracture detection on awkward real-world trace shapes                        #
# --------------------------------------------------------------------------- #
def _df(disp, load):
    return pd.DataFrame({"disp_mm": np.asarray(disp, float),
                         "load_n": np.asarray(load, float)})


def test_fracture_caught_despite_high_residual_load():
    """A brittle snap that leaves a high residual (grip friction) must still be
    found -- the old 'load < 20% of Fmax' rule missed these (see fibre 10_6)."""
    fmax = 0.10
    disp_rise = np.linspace(0.0, 0.5, 500)
    load_rise = fmax * disp_rise / 0.5
    disp_tail = np.linspace(0.505, 2.0, 400)
    load_tail = np.full_like(disp_tail, 0.35 * fmax)  # holds at 35% -> never < 20%
    res = compute_tensile(
        _df(np.r_[disp_rise, disp_tail], np.r_[load_rise, load_tail]),
        diameter_um=10.0, gauge_length_mm=10.0, cfg=CONFIG())
    assert "no_fracture_drop" not in res.flags          # the snap IS detected
    assert res.fmax_n == pytest.approx(fmax, rel=1e-3)
    assert res.strain_at_break == pytest.approx(0.05, abs=5e-3)  # at the snap, not the tail


def test_post_fracture_recoil_spike_is_not_the_breaking_force():
    """A post-test recoil spike larger than the true peak must not become Fmax,
    and the real rupture (not the spike) sets strain-at-break (see fibre 4_3)."""
    fmax = 0.10
    disp_rise = np.linspace(0.0, 0.5, 500)
    load_rise = fmax * disp_rise / 0.5
    disp_flat = np.linspace(0.505, 0.9, 200)
    load_flat = np.zeros_like(disp_flat)               # snapped, load ~0
    disp_spk = np.array([0.95, 0.951, 0.952, 1.0])
    load_spk = np.array([2 * fmax, 2 * fmax, 2 * fmax, 0.0])  # recoil > true peak
    res = compute_tensile(
        _df(np.r_[disp_rise, disp_flat, disp_spk],
            np.r_[load_rise, load_flat, load_spk]),
        diameter_um=10.0, gauge_length_mm=10.0, cfg=CONFIG())
    assert res.fmax_n == pytest.approx(fmax, rel=1e-3)  # the bump, not 2*fmax
    assert "artifact_after_break" in res.flags
    assert res.strain_at_break == pytest.approx(0.05, abs=5e-3)


def test_manual_break_index_overrides_detector():
    """A user-chosen break sample pins the fracture, and every downstream metric
    is taken up to it (the GUI's manual break point)."""
    df, _, _ = _triangular(d_um=10.0, L0_mm=10.0, E_pa=1.0e9, strain_b=0.05)
    res = compute_tensile(df, diameter_um=10.0, gauge_length_mm=10.0, cfg=CONFIG(),
                          break_index=25)  # mid-rise, well before the real snap
    assert "manual_break" in res.flags
    assert res.break_index == 25
    assert res.strain_at_break == pytest.approx(df["disp_mm"].iloc[25] / 10.0, abs=1e-9)
    assert res.fmax_n == pytest.approx(float(df["load_n"].iloc[:26].max()), rel=1e-9)


def test_build_matrix_honours_manual_break(tmp_path):
    df, _, _ = _triangular()
    p = _write_ta_csv(tmp_path / "masp2 3_1 06112026 030624_std.CSV",
                      df["disp_mm"], df["load_n"])
    auto = build_matrix({"3_1": 10.0}, {"3_1": p}, CONFIG())
    manual = build_matrix({"3_1": 10.0}, {"3_1": p}, CONFIG(), breaks={"3_1": 20})
    assert manual.iloc[0]["strain_at_break_pct"] < auto.iloc[0]["strain_at_break_pct"]
    assert "manual_break" in manual.iloc[0]["notes"]


def test_ductile_decline_stays_no_fracture_drop():
    """A gradual post-peak decline that never collapses is left unbroken and
    flagged (see fibre 10_2) rather than snapped at a spurious point."""
    fmax = 0.10
    disp_rise = np.linspace(0.0, 0.5, 500)
    load_rise = fmax * disp_rise / 0.5
    disp_dec = np.linspace(0.501, 3.0, 2000)
    load_dec = np.linspace(fmax, 0.3 * fmax, 2000)     # slow, stays above 20%
    res = compute_tensile(
        _df(np.r_[disp_rise, disp_dec], np.r_[load_rise, load_dec]),
        diameter_um=10.0, gauge_length_mm=10.0, cfg=CONFIG())
    assert "no_fracture_drop" in res.flags
    assert res.fmax_n == pytest.approx(fmax, rel=1e-3)
    assert res.strain_at_break == pytest.approx(0.30, abs=1e-2)  # runs to the end


# --------------------------------------------------------------------------- #
# Matrix assembly                                                              #
# --------------------------------------------------------------------------- #
def test_build_matrix_matches_and_flags(tmp_path):
    df, _, _ = _triangular()
    p31 = _write_ta_csv(tmp_path / "masp2 3_1 06112026 030624_std.CSV",
                        df["disp_mm"], df["load_n"])
    p99 = _write_ta_csv(tmp_path / "masp2 9_9 06102026 014137_std.CSV",
                        df["disp_mm"], df["load_n"])
    diameters = {"3_1": 10.0, "3_2": 12.0}      # 3_2 has an image but no curve
    tensile = {"3_1": p31, "9_9": p99}          # 9_9 has a curve but no image
    mat = build_matrix(diameters, tensile, CONFIG())

    by_group = {r["group"]: r for _, r in mat.iterrows()}
    assert set(by_group) == {"3_1", "3_2", "9_9"}
    assert by_group["3_1"]["flag"] == ""                          # fully matched
    assert np.isfinite(by_group["3_1"]["youngs_modulus_GPa"])
    assert by_group["3_2"]["flag"] == "unmatched_tensile"        # diameter, no curve
    assert not np.isfinite(by_group["3_2"]["fmax_N"])
    assert by_group["9_9"]["flag"] == "unmatched_image"          # curve, no diameter
    assert np.isfinite(by_group["9_9"]["fmax_N"])                # force metric present
    assert not np.isfinite(by_group["9_9"]["youngs_modulus_GPa"])
