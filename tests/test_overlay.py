"""Unit tests for the overlay's perpendicular measurement chords."""

from __future__ import annotations

import numpy as np

from fibrecv.overlay import GREEN, draw_perp_chords, render_overlay


def _tilted_walls(H=300, W=400, slope=0.4, y0=60.0, dv=120.0):
    """Parallel wall lines y_top = y0 + slope*x, y_bot = y_top + dv."""
    x = np.arange(W, dtype=float)
    y_top = y0 + slope * x
    y_bot = y_top + dv
    valid = np.ones(W, dtype=bool)
    return y_top, y_bot, valid


def _green_pixels(img):
    ys, xs = np.where(np.all(img == GREEN, axis=-1))
    return xs, ys


def test_chord_end_lands_on_bottom_wall():
    """The chord starts on the top wall, has the reported perpendicular
    length, and its far end touches the bottom wall line."""
    slope, y0, dv = 0.4, 60.0, 120.0
    y_top, y_bot, valid = _tilted_walls(slope=slope, y0=y0, dv=dv)
    img = np.zeros((300, 400, 3), dtype=np.uint8)
    draw_perp_chords(img, y_top, y_bot, valid, slope, 0, 399,
                     n_chords=1, tick=0)

    xs, ys = _green_pixels(img)
    assert xs.size > 0
    top = (xs[np.argmin(ys)], ys.min())
    bot = (xs[np.argmax(ys)], ys.max())

    # start on the top wall at the sampled column
    assert abs(top[1] - (y0 + slope * top[0])) <= 2.0
    # far end on the bottom wall line (the bug-spotting property)
    assert abs(bot[1] - (y0 + dv + slope * bot[0])) <= 2.0
    # pixel length equals the reported perpendicular diameter
    cth = 1.0 / np.sqrt(1.0 + slope * slope)
    length = np.hypot(bot[0] - top[0], bot[1] - top[1])
    assert abs(length - dv * cth) <= 2.5
    # chord runs perpendicular to the axis: direction (dx,dy) ~ (-sin,cos)
    dx, dy = bot[0] - top[0], bot[1] - top[1]
    assert abs(dx / dy + slope) <= 0.05


def test_render_overlay_chords_opt_in():
    """No slope -> no green chords; slope given -> chords drawn."""
    y_top, y_bot, valid = _tilted_walls()
    rgb = np.zeros((300, 400, 3), dtype=float)
    c_fit = (y_top + y_bot) / 2.0

    plain = render_overlay(rgb, y_top, y_bot, c_fit, valid, 0, 399)
    withc = render_overlay(rgb, y_top, y_bot, c_fit, valid, 0, 399, slope=0.4)
    assert _green_pixels(plain)[0].size == 0
    assert _green_pixels(withc)[0].size > 0


def test_chords_skip_invalid_and_nan_columns():
    slope, y0 = 0.4, 60.0
    y_top, y_bot, valid = _tilted_walls(slope=slope, y0=y0)
    valid[100:300] = False
    y_top[:50] = np.nan
    img = np.zeros((300, 400, 3), dtype=np.uint8)
    draw_perp_chords(img, y_top, y_bot, valid, slope, 0, 399, n_chords=5)
    xs, ys = _green_pixels(img)
    assert xs.size > 0
    # a slanted chord may sweep across invalid columns, but it must START on
    # the top wall at a valid column: no wall-touching green inside the
    # invalid span (8-px margin covers the end ticks)
    on_wall = np.abs(ys - (y0 + slope * xs)) <= 1.5
    assert not np.any((xs[on_wall] >= 108) & (xs[on_wall] <= 292))

    # fully invalid span draws nothing and does not raise
    img2 = np.zeros((300, 400, 3), dtype=np.uint8)
    draw_perp_chords(img2, y_top, y_bot, np.zeros(400, bool), 0.4, 0, 399)
    assert _green_pixels(img2)[0].size == 0
