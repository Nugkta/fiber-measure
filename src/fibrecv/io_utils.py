"""Image discovery, filename parsing and RGB loading.

Dependencies
------------
``numpy``, ``imageio.v3`` (JPEG decode), ``pathlib``, ``re``.

Inputs
------
- A directory of ``masp2 A_B_C.jpg`` microscopy images (2560x1920 RGB) plus
  sidecar ``*.jpg_metadata.xml`` files that must be ignored.
- Glob / group selectors from the CLI.

Output
------
- ``discover_images(root, glob)`` -> sorted list of ``Path`` to .jpg files.
- ``parse_name(path)`` -> ``(group, replicate)`` where group is ``"A_B"`` and
  replicate is the integer ``C``.
- ``load_rgb(path)`` -> float32 array in [0, 1], shape (H, W, 3).

Pos
---
Entry layer of the per-image pipeline. Used by ``measure.py`` (load + parse) and
by both CLIs (discovery + grouping). Knows the dataset's naming convention.
"""

from __future__ import annotations

import re
from pathlib import Path

import imageio.v3 as iio
import numpy as np

# "masp2 10_5_2.jpg" -> A=10, B=5, C=2 ; tolerant of extra spaces.
_NAME_RE = re.compile(r"masp2\s+(\d+)_(\d+)_(\d+)$", re.IGNORECASE)


def discover_images(root: str | Path, glob: str = "masp2 *_*.jpg") -> list[Path]:
    """Return sorted .jpg images under ``root`` matching ``glob``.

    The sidecar ``*.jpg_metadata.xml`` files share the ``masp2 *`` prefix, so we
    explicitly keep only paths whose suffix is exactly ``.jpg``.
    """
    root = Path(root)
    paths = [p for p in root.glob(glob) if p.suffix.lower() == ".jpg"]
    return sorted(paths, key=lambda p: _sort_key(p))


def parse_name(path: str | Path) -> tuple[str, int]:
    """Parse ``masp2 A_B_C.jpg`` -> (``"A_B"``, ``C``).

    Raises ``ValueError`` if the stem does not match the expected pattern.
    """
    stem = Path(path).stem
    m = _NAME_RE.match(stem.strip())
    if not m:
        raise ValueError(f"unrecognised image name: {path!r}")
    a, b, c = m.group(1), m.group(2), m.group(3)
    return f"{a}_{b}", int(c)


def _sort_key(path: Path):
    """Natural sort by (A, B, C) integers so 3_1_2 < 3_1_10 < 10_1_1."""
    try:
        group, rep = parse_name(path)
        a, b = group.split("_")
        return (int(a), int(b), rep)
    except ValueError:
        return (1 << 30, 1 << 30, 1 << 30, path.name)


def load_rgb(path: str | Path) -> np.ndarray:
    """Load a JPEG as float32 RGB in [0, 1], shape (H, W, 3).

    Drops an alpha channel if present and promotes greyscale to 3 channels.
    """
    arr = iio.imread(path)
    arr = np.asarray(arr)
    if arr.ndim == 2:
        arr = np.repeat(arr[:, :, None], 3, axis=2)
    if arr.shape[2] == 4:
        arr = arr[:, :, :3]
    if arr.dtype != np.float32 and arr.dtype != np.float64:
        arr = arr.astype(np.float32) / 255.0
    else:
        arr = arr.astype(np.float32)
        if arr.max() > 1.0:
            arr = arr / 255.0
    return arr
