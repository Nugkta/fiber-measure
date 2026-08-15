"""Tests for the Zeiss XML scale reader (scale.py)."""

from pathlib import Path

import pytest

from fibrecv.scale import ScaleInfo, read_scale, resolve_um_per_px, sidecar_for

_REAL_SIDECAR = Path(
    "/Users/stan/Documents/UOM/spins/multiangle/C1/C1_01_a1_part1.tiff_metadata.xml")


def _write_xml(path: Path, *, cam="2.2,2.2", dist_x="3.8892367379682846E-07",
               dist_y="3.8892367379682846E-07", nominal="10", adapter="1",
               optovar="1", bom=True, drop_camera=False, drop_items=False,
               cam_tag="CameraPixelDistance") -> Path:
    items = "" if drop_items else f"""
    <Items>
      <Pixel Id="Z"><Value>1</Value></Pixel>
      <Distance Id="X"><Value>{dist_x}</Value><DefaultUnitFormat>µm</DefaultUnitFormat></Distance>
      <Distance Id="Y"><Value>{dist_y}</Value><DefaultUnitFormat>µm</DefaultUnitFormat></Distance>
    </Items>"""
    camera = "" if drop_camera else f"""
    <AutoScaling>
      <CameraName>AxioCam ERc5s</CameraName>
      <CameraAdapterMagnification>{adapter}</CameraAdapterMagnification>
      <{cam_tag}>{cam}</{cam_tag}>
    </AutoScaling>"""
    text = f"""<?xml version="1.0" encoding="utf-8"?>
<ImageMetadata>
  <Information>
    <Instrument>
      <Objectives>
        <Objective><NominalMagnification>{nominal}</NominalMagnification></Objective>
      </Objectives>
    </Instrument>
  </Information>
  <OptovarMagnification>{optovar}</OptovarMagnification>
  <Scaling>{camera}{items}
  </Scaling>
</ImageMetadata>
"""
    encoding = "utf-8-sig" if bom else "utf-8"
    path.write_text(text, encoding=encoding)
    return path


def test_sidecar_for_appends_to_full_filename(tmp_path):
    img = tmp_path / "C1_01_a1_part1.tiff"
    assert sidecar_for(img) == tmp_path / "C1_01_a1_part1.tiff_metadata.xml"


def test_read_scale_both_fields(tmp_path):
    xml = _write_xml(tmp_path / "a.xml")
    info = read_scale(xml)
    assert isinstance(info, ScaleInfo)
    assert info.camera_pixel_distance_um == pytest.approx(2.2)
    assert info.total_magnification == pytest.approx(10.0)
    assert info.um_per_px_camera == pytest.approx(0.22)
    assert info.um_per_px_items == pytest.approx(0.38892367, rel=1e-6)
    assert info.source_path == xml


def test_read_scale_accepts_plural_camera_tag(tmp_path):
    xml = _write_xml(tmp_path / "a.xml", cam_tag="CameraPixelDistances")
    assert read_scale(xml).um_per_px_camera == pytest.approx(0.22)


def test_read_scale_no_bom(tmp_path):
    xml = _write_xml(tmp_path / "a.xml", bom=False)
    assert read_scale(xml).um_per_px_camera == pytest.approx(0.22)


def test_read_scale_anisotropic_camera_rejected(tmp_path):
    xml = _write_xml(tmp_path / "a.xml", cam="2.2,3.3")
    with pytest.raises(ValueError):
        read_scale(xml)


def test_read_scale_anisotropic_items_rejected(tmp_path):
    xml = _write_xml(tmp_path / "a.xml", dist_y="9.9E-07")
    with pytest.raises(ValueError):
        read_scale(xml)


def test_read_scale_missing_camera_rejected(tmp_path):
    xml = _write_xml(tmp_path / "a.xml", drop_camera=True)
    with pytest.raises(ValueError):
        read_scale(xml)


def test_read_scale_missing_items_rejected(tmp_path):
    xml = _write_xml(tmp_path / "a.xml", drop_items=True)
    with pytest.raises(ValueError):
        read_scale(xml)


def test_resolve_um_per_px(tmp_path):
    info = read_scale(_write_xml(tmp_path / "a.xml"))
    assert resolve_um_per_px(info, "camera") == pytest.approx(0.22)
    assert resolve_um_per_px(info, "items") == pytest.approx(0.38892367, rel=1e-6)
    assert resolve_um_per_px(info, "0.25") == pytest.approx(0.25)
    assert resolve_um_per_px(info, 0.25) == pytest.approx(0.25)
    with pytest.raises(ValueError):
        resolve_um_per_px(info, "bogus")
    with pytest.raises(ValueError):
        resolve_um_per_px(info, "-1.0")


@pytest.mark.skipif(not _REAL_SIDECAR.exists(), reason="C1 data not on this machine")
def test_read_scale_real_sidecar():
    info = read_scale(_REAL_SIDECAR)
    assert info.um_per_px_camera == pytest.approx(0.22, rel=1e-3)
    assert info.um_per_px_items == pytest.approx(0.38892367, rel=1e-4)
    assert info.total_magnification == pytest.approx(10.0)
