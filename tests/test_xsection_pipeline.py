"""End-to-end tests for the run_xsection CLI on a fabricated mini-tree."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from fibrecv import run_xsection
from fibrecv.xsection import NOMINAL_ANGLES_DEG

A_TRUE, B_TRUE, PHI_TRUE = 100.0, 80.0, 35.0
UM_ITEMS = 0.38892367379682846
UM_CAMERA = 0.22

EXPECTED_PART_COLS = [
    "x_px", "w_a1_px", "w_a2_px", "w_a3_px", "w_a4_px", "w_a5_px", "w_a6_px",
    "n_angles", "a_px", "b_px", "phi_deg", "rms_resid_px", "a_um", "b_um",
    "area_um2", "area_err_um2", "area_hex_um2", "hex_ratio",
    "hex_ratio_expected", "valid",
]

EXPECTED_SUMMARY_COLS = [
    "fiber", "n_parts", "n_positions", "um_per_px", "A_mean_um2", "A_harm_um2",
    "A_min_um2", "A_min_part", "A_min_x_px", "axis_ratio_median",
    "axis_ratio_iqr", "phi_med_deg", "hex_ratio_median", "pair_dw_frac_median",
    "A_circle_um2", "area_ratio", "uniformity", "n_uncertain_shifts",
    "n_saturated_shifts", "part_rms_med_max_px", "low_confidence",
]


def ellipse_width(angle_index: int, a=A_TRUE, b=B_TRUE, phi=PHI_TRUE) -> float:
    th = np.radians(NOMINAL_ANGLES_DEG[angle_index - 1])
    ph = np.radians(phi)
    return float(2 * np.sqrt(a ** 2 * np.cos(th - ph) ** 2
                             + b ** 2 * np.sin(th - ph) ** 2))


def _write_sidecar(data_root, image_name, items_m=UM_ITEMS * 1e-6):
    xml = f"""<?xml version="1.0" encoding="utf-8"?>
<ImageMetadata>
  <NominalMagnification>10</NominalMagnification>
  <CameraAdapterMagnification>1</CameraAdapterMagnification>
  <OptovarMagnification>1</OptovarMagnification>
  <Scaling>
    <AutoScaling><CameraPixelDistance>2.2,2.2</CameraPixelDistance></AutoScaling>
    <Items>
      <Distance Id="X"><Value>{items_m!r}</Value></Distance>
      <Distance Id="Y"><Value>{items_m!r}</Value></Distance>
    </Items>
  </Scaling>
</ImageMetadata>
"""
    (data_root / f"{image_name}_metadata.xml").write_text(xml, encoding="utf-8")


def _write_image_profile(out_root, base, width, x0=100, x1=400,
                         band_mismatch=False, coverage=1.0):
    csv_dir = out_root / "per_image" / "csv"
    meta_dir = out_root / "per_image" / "diagnostics"
    csv_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)
    x = np.arange(x0, x1, dtype=float)
    df = pd.DataFrame({
        "x_px": x,
        "diameter_px_raw": np.full(x.size, width),
        "diameter_px_smooth": np.full(x.size, width),
        "diameter_um": np.full(x.size, width / 2.5712),
        "valid": np.ones(x.size, dtype=int),
        "interpolated": np.zeros(x.size, dtype=int),
        "y_top_px": np.full(x.size, 500.0 - width / 2),
        "y_bot_px": np.full(x.size, 500.0 + width / 2),
    })
    df.to_csv(csv_dir / f"{base}_profile.csv", index=False)
    meta = {"coverage": coverage, "band_mismatch": band_mismatch,
            "low_confidence": False, "anomaly": {"flags": []}}
    with open(meta_dir / f"{base}_meta.json", "w") as fh:
        json.dump(meta, fh)


@pytest.fixture()
def mini_tree(tmp_path):
    out_root = tmp_path / "out"
    data_root = tmp_path / "data"
    data_root.mkdir()
    # fiber 1 part 1: all six angles of the known ellipse
    for a in range(1, 7):
        base = f"C1_01_a{a}_part1"
        _write_image_profile(out_root, base, ellipse_width(a))
        _write_sidecar(data_root, f"{base}.tiff")
    return out_root, data_root


def test_pipeline_known_ellipse(mini_tree):
    out_root, data_root = mini_tree
    rc = run_xsection.main([
        "--out", str(out_root), "--data-root", str(data_root),
        "--condition", "C1", "--scale-source", "items",
    ])
    assert rc == 0

    part_csv = out_root / "per_part" / "csv" / "xsec_C1_01_p1.csv"
    df = pd.read_csv(part_csv)
    assert list(df.columns) == EXPECTED_PART_COLS
    v = df[df["valid"] == 1]
    assert len(v) > 250
    assert np.allclose(v["a_px"], A_TRUE, rtol=1e-6)
    assert np.allclose(v["b_px"], B_TRUE, rtol=1e-6)
    assert np.allclose(v["phi_deg"], PHI_TRUE, atol=1e-4)
    assert np.allclose(v["n_angles"], 6)
    assert np.allclose(v["area_um2"],
                       np.pi * A_TRUE * B_TRUE * UM_ITEMS ** 2, rtol=1e-6)
    assert np.allclose(v["area_err_um2"], 0.0, atol=1e-6)
    # hexagon QC pair: measured ratio equals expected-for-fit ratio exactly
    assert np.allclose(v["hex_ratio"], v["hex_ratio_expected"], rtol=1e-9)
    assert (v["hex_ratio"] < np.pi / (2 * np.sqrt(3))).all()  # ellipse < circle anchor

    assert (out_root / "per_part" / "plots" / "xsec_C1_01_p1.png").exists()
    shifts = json.loads(
        (out_root / "per_part" / "shifts" / "xsec_C1_01_p1.json").read_text())
    assert shifts["um_per_px_resolved"] == pytest.approx(UM_ITEMS, rel=1e-9)
    assert shifts["um_per_px_camera"] == pytest.approx(UM_CAMERA, rel=1e-6)
    assert len(shifts["shifts"]) == 6

    summary = pd.read_csv(out_root / "summary" / "xsection_summary.csv")
    assert list(summary.columns) == EXPECTED_SUMMARY_COLS
    assert len(summary) == 1
    row = summary.iloc[0]
    assert row["fiber"] == 1
    assert row["n_parts"] == 1
    A_um = np.pi * A_TRUE * B_TRUE * UM_ITEMS ** 2
    assert row["A_mean_um2"] == pytest.approx(A_um, rel=1e-6)
    assert row["A_harm_um2"] == pytest.approx(A_um, rel=1e-6)
    assert row["A_min_um2"] == pytest.approx(A_um, rel=1e-6)
    assert row["axis_ratio_median"] == pytest.approx(A_TRUE / B_TRUE, rel=1e-6)
    assert row["phi_med_deg"] == pytest.approx(PHI_TRUE, abs=1e-3)
    # circle comparison: d_bar = grand mean width over the six angles
    d_bar = np.mean([ellipse_width(a) for a in range(1, 7)]) * UM_ITEMS
    assert row["A_circle_um2"] == pytest.approx(np.pi * (d_bar / 2) ** 2, rel=1e-6)
    assert row["area_ratio"] == pytest.approx(A_um / (np.pi * (d_bar / 2) ** 2),
                                              rel=1e-6)
    assert row["uniformity"] == pytest.approx(1.0, rel=1e-6)

    assert (out_root / "summary" / "xsection_angle_residuals.csv").exists()
    cfg_json = json.loads(
        (out_root / "summary" / "xsection_run_config.json").read_text())
    assert cfg_json["scale_source"] == "items"


def test_pipeline_numeric_scale_source(mini_tree):
    out_root, data_root = mini_tree
    rc = run_xsection.main([
        "--out", str(out_root), "--data-root", str(data_root),
        "--scale-source", "0.5",
    ])
    assert rc == 0
    df = pd.read_csv(out_root / "per_part" / "csv" / "xsec_C1_01_p1.csv")
    v = df[df["valid"] == 1]
    assert np.allclose(v["area_um2"], np.pi * A_TRUE * B_TRUE * 0.25, rtol=1e-6)


def test_pipeline_missing_angle_and_excluded_image(mini_tree):
    out_root, data_root = mini_tree
    # part 2: a3 never measured; a5 measured but band_mismatch -> excluded
    for a in [1, 2, 4, 5, 6]:
        base = f"C1_01_a{a}_part2"
        _write_image_profile(out_root, base, ellipse_width(a),
                             band_mismatch=(a == 5))
        _write_sidecar(data_root, f"{base}.tiff")
    rc = run_xsection.main([
        "--out", str(out_root), "--data-root", str(data_root),
        "--scale-source", "items",
    ])
    assert rc == 0
    df = pd.read_csv(out_root / "per_part" / "csv" / "xsec_C1_01_p2.csv")
    assert df["w_a3_px"].isna().all()
    assert df["w_a5_px"].isna().all()
    v = df[df["valid"] == 1]
    assert len(v) > 250  # d3 still covered by a6 -> fit remains valid
    assert np.allclose(v["a_px"], A_TRUE, rtol=1e-6)
    summary = pd.read_csv(out_root / "summary" / "xsection_summary.csv")
    assert summary.iloc[0]["n_parts"] == 2


def test_pipeline_disagreeing_scales_abort(tmp_path):
    out_root = tmp_path / "out"
    data_root = tmp_path / "data"
    data_root.mkdir()
    for a in range(1, 7):
        base = f"C1_01_a{a}_part1"
        _write_image_profile(out_root, base, ellipse_width(a))
        items = UM_ITEMS * 1e-6 * (1.01 if a == 4 else 1.0)  # 1% off on a4
        _write_sidecar(data_root, f"{base}.tiff", items_m=items)
    rc = run_xsection.main([
        "--out", str(out_root), "--data-root", str(data_root),
        "--scale-source", "items",
    ])
    assert rc != 0


def test_pipeline_fibers_filter(mini_tree):
    out_root, data_root = mini_tree
    for a in range(1, 7):
        base = f"C1_02_a{a}_part1"
        _write_image_profile(out_root, base, ellipse_width(a))
        _write_sidecar(data_root, f"{base}.tiff")
    rc = run_xsection.main([
        "--out", str(out_root), "--data-root", str(data_root),
        "--scale-source", "items", "--fibers", "2",
    ])
    assert rc == 0
    assert (out_root / "per_part" / "csv" / "xsec_C1_02_p1.csv").exists()
    assert not (out_root / "per_part" / "csv" / "xsec_C1_01_p1.csv").exists()
    summary = pd.read_csv(out_root / "summary" / "xsection_summary.csv")
    assert summary["fiber"].tolist() == [2]


def test_pipeline_validation_writer(mini_tree):
    out_root, data_root = mini_tree
    rc = run_xsection.main([
        "--out", str(out_root), "--data-root", str(data_root),
        "--scale-source", "items", "--validation",
    ])
    assert rc == 0
    val = pd.read_csv(out_root / "summary" / "xsection_validation.csv")
    assert set(val["kind"]) == {"anisotropy", "phi_transfer"}
    ani = val[val["kind"] == "anisotropy"]
    # exact ellipse: directional prediction is exact, circle baseline is not
    assert (ani["rmse_directional"] < 1e-6).all()
    assert (ani["rmse_circle"] > 1.0).all()


def test_pipeline_truncated_sidecar_rc2_not_traceback(mini_tree):
    # ET.ParseError subclasses SyntaxError — it must land on the rc=2 scale
    # path like every other bad-sidecar failure, not escape as a traceback
    out_root, data_root = mini_tree
    (data_root / "C1_01_a3_part1.tiff_metadata.xml").write_text(
        '<?xml version="1.0"?><ImageMetada', encoding="utf-8")
    rc = run_xsection.main([
        "--out", str(out_root), "--data-root", str(data_root),
        "--scale-source", "items",
    ])
    assert rc == 2


def test_pipeline_garbled_profile_csv_skipped(mini_tree, capsys):
    # a truncated artifact from a crashed run_measure is skipped with a
    # warning (image treated as absent), never a KeyError crash
    out_root, data_root = mini_tree
    csv = out_root / "per_image" / "csv" / "C1_01_a3_part1_profile.csv"
    csv.write_text("this,is\nnot,a,profile\n")
    rc = run_xsection.main([
        "--out", str(out_root), "--data-root", str(data_root),
        "--scale-source", "items",
    ])
    assert rc == 0
    assert "skipping unreadable profile" in capsys.readouterr().out
    part = pd.read_csv(out_root / "per_part" / "csv" / "xsec_C1_01_p1.csv")
    assert part["w_a3_px"].isna().all()
