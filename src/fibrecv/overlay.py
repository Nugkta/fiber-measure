"""Draw the boundary overlay PNG (user-requested review artifact).

Dependencies
------------
``numpy``, ``imageio.v3`` (lossless PNG write).

Inputs
------
- ``rgb``: original float RGB image in [0, 1], shape (H, W, 3).
- ``y_top`` / ``y_bot`` / ``c_fit``: per-column boundary + centerline arrays.
- ``valid``: bool mask of accepted columns (boundaries only drawn where valid).
- ``out_path``: destination PNG.

Output
------
- ``render_overlay(...)`` -> the full-resolution (H, W, 3) uint8 RGB array with
  the boundaries drawn on (pure; the GUI shows it directly via ``st.image``).
- ``draw_overlay(...)`` renders that array and writes it as a lossless PNG.
Boundaries: top in cyan, bottom in yellow, fitted centerline as a dashed grey
line, drawn directly onto the pixel grid (no matplotlib rescaling) so edges can
be checked pixel-for-pixel.

Pos
---
Side-output of ``measure.py`` and the live preview of ``gui_app.py``. This is the
primary human review artifact for confirming boundaries hug the true fibre edge
and exclude shadow/halo.
"""

from __future__ import annotations

from pathlib import Path

import imageio.v3 as iio
import numpy as np

CYAN = np.array([0, 255, 255], dtype=np.uint8)
YELLOW = np.array([255, 255, 0], dtype=np.uint8)
GREY = np.array([180, 180, 180], dtype=np.uint8)


def _stamp(img: np.ndarray, x: int, y: float, color: np.ndarray, thick: int = 1) -> None:
    """Colour a small vertical mark of half-width ``thick`` at column x, row y."""
    if not np.isfinite(y):
        return
    H = img.shape[0]
    yc = int(round(y))
    lo, hi = max(0, yc - thick), min(H, yc + thick + 1)
    img[lo:hi, x] = color


def render_overlay(
    rgb: np.ndarray,
    y_top: np.ndarray,
    y_bot: np.ndarray,
    c_fit: np.ndarray,
    valid: np.ndarray,
    x0: int,
    x1: int,
    thick: int = 1,
) -> np.ndarray:
    """Draw the boundaries onto a copy of ``rgb`` and return the uint8 array.

    Pure: no disk I/O. Returns an (H, W, 3) uint8 RGB image with the top boundary
    in cyan, bottom in yellow and a dashed grey centerline, ready for either
    ``st.image`` (GUI) or ``imwrite`` (CLI via ``draw_overlay``).
    """
    img = (np.clip(rgb, 0, 1) * 255).astype(np.uint8).copy()
    W = img.shape[1]
    for x in range(max(0, x0), min(W, x1 + 1)):
        # dashed centerline (every other 8-px run) for context
        if (x // 8) % 2 == 0:
            _stamp(img, x, c_fit[x], GREY, thick=0)
        if x < valid.size and valid[x]:
            _stamp(img, x, y_top[x], CYAN, thick=thick)
            _stamp(img, x, y_bot[x], YELLOW, thick=thick)
    return img


def draw_overlay(
    rgb: np.ndarray,
    y_top: np.ndarray,
    y_bot: np.ndarray,
    c_fit: np.ndarray,
    valid: np.ndarray,
    out_path: str | Path,
    x0: int,
    x1: int,
    thick: int = 1,
) -> None:
    """Render the boundary overlay and save it as a lossless PNG."""
    img = render_overlay(rgb, y_top, y_bot, c_fit, valid, x0, x1, thick=thick)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    iio.imwrite(out_path, img)
