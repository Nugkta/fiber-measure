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
- ``draw_perp_chords(...)`` -> green measurement chords (used by the above).
Boundaries: top in cyan, bottom in yellow, fitted centerline as a dashed grey
line, drawn directly onto the pixel grid (no matplotlib rescaling) so edges can
be checked pixel-for-pixel. When a tilt ``slope`` is given, a few green chords
show the exact perpendicular diameter being reported: each starts on the top
boundary and runs perpendicular to the fibre axis for
``(y_bot - y_top) * cos(theta)`` px, so its far end must land on the bottom
boundary — a gap, overshoot or wrong angle exposes a tilt bug at a glance.

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

from .band import tilt_geometry

CYAN = np.array([0, 255, 255], dtype=np.uint8)
YELLOW = np.array([255, 255, 0], dtype=np.uint8)
GREY = np.array([180, 180, 180], dtype=np.uint8)
MAGENTA = np.array([255, 0, 255], dtype=np.uint8)  # manually corrected columns
WHITE = np.array([255, 255, 255], dtype=np.uint8)
GREEN = np.array([0, 230, 60], dtype=np.uint8)  # measured perpendicular chords


def _stamp(img: np.ndarray, x: int, y: float, color: np.ndarray, thick: int = 1) -> None:
    """Colour a small vertical mark of half-width ``thick`` at column x, row y."""
    if not np.isfinite(y):
        return
    H = img.shape[0]
    yc = int(round(y))
    lo, hi = max(0, yc - thick), min(H, yc + thick + 1)
    img[lo:hi, x] = color


def _draw_segment(
    img: np.ndarray,
    xa: float,
    ya: float,
    xb: float,
    yb: float,
    color: np.ndarray,
    brush: int = 2,
) -> None:
    """Stamp the straight segment (xa, ya)->(xb, yb), clipped to the image.

    Dense point sampling with a square brush of side ``2*brush - 1`` px — no
    matplotlib, so the overlay stays pixel-exact.
    """
    H, W = img.shape[:2]
    n = int(2 * max(abs(xb - xa), abs(yb - ya))) + 2
    xs = np.rint(np.linspace(xa, xb, n)).astype(int)
    ys = np.rint(np.linspace(ya, yb, n)).astype(int)
    for dy in range(-brush + 1, brush):
        for dx in range(-brush + 1, brush):
            xk, yk = xs + dx, ys + dy
            ok = (xk >= 0) & (xk < W) & (yk >= 0) & (yk < H)
            img[yk[ok], xk[ok]] = color


def draw_perp_chords(
    img: np.ndarray,
    y_top: np.ndarray,
    y_bot: np.ndarray,
    valid: np.ndarray,
    slope: float,
    x0: int,
    x1: int,
    n_chords: int = 7,
    color: np.ndarray = GREEN,
    tick: int = 6,
) -> None:
    """Draw the measured perpendicular diameter at a few columns, in place.

    At ``n_chords`` evenly spaced valid columns, a chord starts on the top
    boundary and runs perpendicular to the (clamped) fibre axis for exactly the
    reported length ``(y_bot - y_top) * cos(theta)``. Geometry check by eye:
    the far end must land on the bottom boundary and the chord must cross the
    fibre at a right angle — a gap, overshoot or skewed angle means the tilt
    correction is wrong. ``tick``-px end bars run along the wall direction so
    the touch-down points are easy to judge. Chords use the automatic band
    tilt; manually edited columns re-fit their own slope for the number they
    report, so a small end-gap on heavily edited stretches is expected.
    """
    if n_chords <= 0:
        return
    m, cth = tilt_geometry(slope)
    sth = m * cth  # sin(theta), signed
    hi = min(img.shape[1], x1 + 1, valid.size, y_top.size, y_bot.size)
    xs = np.arange(max(0, x0), hi)
    if xs.size == 0:
        return
    cand = xs[valid[xs] & np.isfinite(y_top[xs]) & np.isfinite(y_bot[xs])]
    if cand.size == 0:
        return
    # interior sample positions: n_chords picks strictly between the span ends
    idx = np.unique(np.round(np.linspace(0, cand.size - 1, n_chords + 2)
                             ).astype(int)[1:-1])
    for x in cand[idx]:
        length = (y_bot[x] - y_top[x]) * cth
        xa, ya = float(x), float(y_top[x])
        xb, yb = xa - length * sth, ya + length * cth
        _draw_segment(img, xa, ya, xb, yb, color)
        if tick > 0:
            for px, py in ((xa, ya), (xb, yb)):
                _draw_segment(img, px - tick * cth, py - tick * sth,
                              px + tick * cth, py + tick * sth, color)


def render_overlay(
    rgb: np.ndarray,
    y_top: np.ndarray,
    y_bot: np.ndarray,
    c_fit: np.ndarray,
    valid: np.ndarray,
    x0: int,
    x1: int,
    thick: int = 1,
    *,
    edited_top: np.ndarray | None = None,
    edited_bot: np.ndarray | None = None,
    slope: float | None = None,
    n_chords: int = 7,
) -> np.ndarray:
    """Draw the boundaries onto a copy of ``rgb`` and return the uint8 array.

    Pure: no disk I/O. Returns an (H, W, 3) uint8 RGB image with the top boundary
    in cyan, bottom in yellow and a dashed grey centerline, ready for either
    ``st.image`` (GUI) or ``imwrite`` (CLI via ``draw_overlay``). Columns marked
    in the optional ``edited_top``/``edited_bot`` bool masks (manual GUI
    corrections) are drawn in magenta instead. When ``slope`` is given,
    ``n_chords`` green perpendicular chords show the measured diameter
    (``draw_perp_chords``).
    """
    img = (np.clip(rgb, 0, 1) * 255).astype(np.uint8).copy()
    W = img.shape[1]
    for x in range(max(0, x0), min(W, x1 + 1)):
        # dashed centerline (every other 8-px run) for context
        if (x // 8) % 2 == 0:
            _stamp(img, x, c_fit[x], GREY, thick=0)
        if x < valid.size and valid[x]:
            top_c = MAGENTA if (edited_top is not None and edited_top[x]) else CYAN
            bot_c = MAGENTA if (edited_bot is not None and edited_bot[x]) else YELLOW
            _stamp(img, x, y_top[x], top_c, thick=thick)
            _stamp(img, x, y_bot[x], bot_c, thick=thick)
    if slope is not None:
        draw_perp_chords(img, y_top, y_bot, valid, slope, x0, x1,
                         n_chords=n_chords)
    return img


def mark_anchors(
    img: np.ndarray,
    anchors,
    color: np.ndarray = WHITE,
    size: int = 4,
) -> None:
    """Stamp a small filled square at each (x, y) anchor, in place, clipped."""
    H, W = img.shape[:2]
    for ax, ay in anchors:
        if not (np.isfinite(ax) and np.isfinite(ay)):
            continue
        xc, yc = int(round(ax)), int(round(ay))
        x_lo, x_hi = max(0, xc - size), min(W, xc + size + 1)
        y_lo, y_hi = max(0, yc - size), min(H, yc + size + 1)
        if x_lo < x_hi and y_lo < y_hi:
            img[y_lo:y_hi, x_lo:x_hi] = color


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
    slope: float | None = None,
) -> None:
    """Render the boundary overlay and save it as a lossless PNG."""
    img = render_overlay(rgb, y_top, y_bot, c_fit, valid, x0, x1, thick=thick,
                         slope=slope)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    iio.imwrite(out_path, img)
