"""Zeiss XML sidecar scale reader.

Dependencies
------------
Stdlib only: ``xml.etree.ElementTree``, ``dataclasses``, ``pathlib``.

Inputs
------
- ``*_metadata.xml`` sidecars written by Zeiss ZEN next to each TIFF, holding
  two independent µm/px candidates: the camera sensor pitch
  (``CameraPixelDistance[s]``, µm, divided by total magnification) and the
  precomputed ``Scaling/Items/Distance`` values (metres per pixel).

Output
------
- ``sidecar_for(image_path)`` -> path of the sidecar for an image.
- ``read_scale(xml_path)`` -> ``ScaleInfo`` with both µm/px candidates.
- ``resolve_um_per_px(info, source)`` -> float µm/px for source
  ``"camera"`` | ``"items"`` | numeric literal.

Pos
---
Third-stage helper: ``run_xsection.py`` converts px to µm exclusively through
this module (never ``config.ppu``). The two XML fields disagree by 1.77x on
C1; the scale-bar cross-check (labbook 02) adjudicates which one is real.

Once I am updated, update my header comments and folder's md.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

# The camera pixel pitch tag appears as CameraPixelDistance (singular, what
# the C1 sidecars actually contain) in some ZEN exports and
# CameraPixelDistances (plural) in others — accept both.
_CAMERA_TAGS = ("CameraPixelDistance", "CameraPixelDistances")
_MAG_TAGS = ("NominalMagnification", "CameraAdapterMagnification",
             "OptovarMagnification")


@dataclass(frozen=True)
class ScaleInfo:
    """Both µm/px candidates read from one Zeiss sidecar."""

    um_per_px_camera: float
    um_per_px_items: float
    camera_pixel_distance_um: float
    total_magnification: float
    source_path: Path


def sidecar_for(image_path: str | Path) -> Path:
    """Sidecar path for an image: ``_metadata.xml`` appended to the FULL name
    (``C1_01_a1_part1.tiff`` -> ``C1_01_a1_part1.tiff_metadata.xml``)."""
    image_path = Path(image_path)
    return image_path.with_name(image_path.name + "_metadata.xml")


def _find_text(root: ET.Element, tag: str, xml_path: Path) -> str:
    node = root.find(f".//{tag}")
    if node is None or node.text is None:
        raise ValueError(f"missing <{tag}> in {xml_path}")
    return node.text.strip()


def read_scale(xml_path: str | Path) -> ScaleInfo:
    """Read both µm/px candidates from a Zeiss sidecar.

    Raises ``ValueError`` on a missing field or an anisotropic (X != Y)
    pixel in either candidate.
    """
    xml_path = Path(xml_path)
    root = ET.parse(xml_path).getroot()

    # Candidate 1: camera sensor pitch / total optical magnification.
    for tag in _CAMERA_TAGS:
        node = root.find(f".//{tag}")
        if node is not None and node.text:
            pair = node.text.strip().split(",")
            break
    else:
        raise ValueError(
            f"missing <{'|'.join(_CAMERA_TAGS)}> in {xml_path}")
    values = [float(v) for v in pair]
    if len(values) != 2 or values[0] != values[1]:
        raise ValueError(
            f"anisotropic camera pixel distance {pair} in {xml_path}")
    cam_um = values[0]

    magnification = 1.0
    for tag in _MAG_TAGS:
        magnification *= float(_find_text(root, tag, xml_path))
    if magnification <= 0:
        raise ValueError(f"non-positive magnification in {xml_path}")

    # Candidate 2: precomputed Scaling/Items distances (metres per pixel).
    dist = {}
    for axis in ("X", "Y"):
        node = root.find(f".//Scaling/Items/Distance[@Id='{axis}']/Value")
        if node is None or node.text is None:
            raise ValueError(
                f"missing Scaling/Items/Distance[{axis}] in {xml_path}")
        dist[axis] = float(node.text)
    if dist["X"] != dist["Y"]:
        raise ValueError(
            f"anisotropic Items distance {dist} in {xml_path}")

    return ScaleInfo(
        um_per_px_camera=cam_um / magnification,
        um_per_px_items=dist["X"] * 1e6,
        camera_pixel_distance_um=cam_um,
        total_magnification=magnification,
        source_path=xml_path,
    )


def resolve_um_per_px(info: ScaleInfo, source: str | float) -> float:
    """Resolve a ``--scale-source`` value to µm/px.

    ``"camera"`` and ``"items"`` pick the corresponding XML candidate; any
    numeric literal (str or float) is used directly. Raises ``ValueError``
    for anything else or a non-positive value.
    """
    if source == "camera":
        return info.um_per_px_camera
    if source == "items":
        return info.um_per_px_items
    try:
        value = float(source)
    except (TypeError, ValueError):
        raise ValueError(f"unrecognised scale source: {source!r}") from None
    if value <= 0:
        raise ValueError(f"non-positive um/px: {source!r}")
    return value
