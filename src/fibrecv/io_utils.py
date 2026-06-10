"""Image discovery, filename parsing and RGB loading.

Dependencies
------------
``numpy``, ``imageio.v3`` (JPEG decode), ``pathlib``, ``re``.

Inputs
------
- A directory of microscopy images whose stems end in numbers (e.g.
  ``masp2 A_B_C.jpg``, ``3-1-2.png``, ``IMG_0123.jpg``) plus sidecar
  ``*.jpg_metadata.xml`` files that must be ignored.
- Glob / group selectors from the CLI.

Output
------
- ``discover_images(root, glob)`` -> sorted list of ``Path`` to .jpg files.
- ``parse_name(path)`` -> ``(group, replicate)`` from the trailing run of
  numbers in the stem (last number = replicate, the rest = group).
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

# Separators that may appear between name tokens: space _ - . ( ) [ ]
_SEP_RE = re.compile(r"[\s_\-.()\[\]]+")
_DIGITS_RE = re.compile(r"(\d+)")


def natural_key(text: str):
    """Sort key interleaving numeric and text runs: '3_1' < '3_3' < '10_5'."""
    return tuple((0, int(t)) if t.isdigit() else (1, t.lower())
                 for t in _DIGITS_RE.split(text) if t)


def discover_images(root: str | Path, glob: str = "masp2 *_*.jpg") -> list[Path]:
    """Return sorted .jpg images under ``root`` matching ``glob``.

    The sidecar ``*.jpg_metadata.xml`` files share the ``masp2 *`` prefix, so we
    explicitly keep only paths whose suffix is exactly ``.jpg``.
    """
    root = Path(root)
    paths = [p for p in root.glob(glob) if p.suffix.lower() == ".jpg"]
    return sorted(paths, key=lambda p: _sort_key(p))


def parse_name(path: str | Path) -> tuple[str, int]:
    """Parse an image name into ``(group, replicate)`` via trailing numbers.

    The stem is split on spaces/underscores/dashes/dots/brackets and the
    trailing run of integer tokens drives the result:

    - two or more trailing integers: the last is the replicate, the rest
      joined with ``_`` form the group ('masp2 10_5_2' -> ('10_5', 2),
      '3-1-2' -> ('3_1', 2), 'sampleA 10_5_2' -> ('10_5', 2));
    - exactly one trailing integer: it is the replicate and the text prefix
      is the group ('IMG_0123' -> ('IMG', 123)).

    Raises ``ValueError`` if the stem does not end in an integer token, or a
    single trailing integer has no text prefix (e.g. '3.jpg').
    """
    stem = Path(path).stem.strip()
    tokens = [t for t in _SEP_RE.split(stem) if t]
    i = len(tokens)
    while i > 0 and tokens[i - 1].isdigit():
        i -= 1
    nums = tokens[i:]
    if not nums:
        raise ValueError(f"unrecognised image name: {path!r}")
    replicate = int(nums[-1])
    if len(nums) >= 2:
        return "_".join(nums[:-1]), replicate
    prefix = " ".join(tokens[:i])
    if not prefix:
        raise ValueError(f"unrecognised image name: {path!r}")
    return prefix, replicate


def _sort_key(path: Path):
    """Sort parseable names by (group, replicate); the rest by filename."""
    try:
        group, rep = parse_name(path)
        return (0, natural_key(group), rep)
    except ValueError:
        return (1, natural_key(Path(path).name), 0)


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
