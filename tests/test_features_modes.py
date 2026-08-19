"""Tests for the feature_mode dispatch: "desat" (MasP2) vs "bright" (C1)."""

from __future__ import annotations

import argparse
import math

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from skimage.color import rgb2hsv

from fibrecv.compute import compute_measurement
from fibrecv.config import CONFIG
from fibrecv.features import estimate_bg, rgb_to_desaturation

TRUE_W = 60.0


def _pink_fibre(W: int = 400, H: int = 300, seed: int = 0) -> np.ndarray:
    """Pale fibre on saturated pink (the MasP2 look, as in test_edges_tilt)."""
    y, x = np.mgrid[0:H, 0:W].astype(np.float64)
    dist = np.abs(y - H / 2)
    t = np.clip((TRUE_W / 2 + 1.0 - dist) / 2.0, 0.0, 1.0)[..., None]
    bg = np.array([0.95, 0.45, 0.65])
    fg = np.array([0.92, 0.90, 0.91])
    img = bg + (fg - bg) * t
    rng = np.random.default_rng(seed)
    img = img + rng.normal(0.0, 0.006, img.shape)
    return np.clip(img, 0.0, 1.0).astype(np.float32)


def _bright_fibre(fg_level: float = 0.85, noise: float = 0.006,
                  angle_deg: float = 3.0, W: int = 800, H: int = 700,
                  seed: int = 0) -> np.ndarray:
    """Bright fibre on dark grey background (the C1 look).

    Geometry copied from test_edges_tilt._inclined_fibre; only the colours
    differ: near-neutral dark bg (~0.20) and bright near-neutral fibre.
    """
    m = math.tan(math.radians(angle_deg))
    y, x = np.mgrid[0:H, 0:W].astype(np.float64)
    dist = np.abs(y - (m * (x - W / 2) + H / 2)) / math.sqrt(1 + m * m)
    t = np.clip((TRUE_W / 2 + 1.0 - dist) / 2.0, 0.0, 1.0)[..., None]
    bg = np.array([0.21, 0.20, 0.19])   # dark grey, slight warm cast like C1
    fg = np.array([fg_level, fg_level - 0.01, fg_level - 0.02])
    img = bg + (fg - bg) * t
    rng = np.random.default_rng(seed)
    img = img + rng.normal(0.0, noise, img.shape)
    return np.clip(img, 0.0, 1.0).astype(np.float32)


def test_desat_mode_bit_identical_to_original_formula():
    """feature_mode="desat" (the default) must reproduce the original
    desaturation computation exactly — MasP2 behaviour is frozen."""
    rgb = _pink_fibre()
    cfg = CONFIG()
    assert cfg.feature_mode == "desat"
    D, S, s_bg, mad = rgb_to_desaturation(rgb, cfg)
    # original formula, inline
    S_ref = rgb2hsv(rgb)[:, :, 1].astype(np.float32)
    s_bg_ref, mad_ref = estimate_bg(S_ref, cfg)
    D_ref = ((s_bg_ref - S_ref) / (cfg.mad_scale * mad_ref + cfg.eps)).astype(np.float32)
    assert s_bg == s_bg_ref and mad == mad_ref
    assert_array_equal(S, S_ref)
    assert_array_equal(D, D_ref)


def test_bright_mode_formula():
    rgb = _bright_fibre()
    cfg = CONFIG(feature_mode="bright")
    D, V, v_bg, mad = rgb_to_desaturation(rgb, cfg)
    V_ref = np.median(rgb, axis=2).astype(np.float32)
    v_bg_ref, mad_ref = estimate_bg(V_ref, cfg)
    D_ref = ((V_ref - v_bg_ref) / (cfg.mad_scale * mad_ref + cfg.eps)).astype(np.float32)
    assert v_bg == v_bg_ref and mad == mad_ref
    assert_array_equal(V, V_ref)
    assert_array_equal(D, D_ref)
    # fibre core strongly positive, background ~0
    assert np.median(D[345:355, :]) > 8.0
    assert abs(np.median(D[:40, :])) < 1.0


def test_invalid_mode_raises():
    with pytest.raises(ValueError):
        rgb_to_desaturation(_bright_fibre(), CONFIG(feature_mode="hsv"))


def test_bright_fibre_end_to_end_bright_mode():
    """C1 bright end: realistic noise puts the fibre at z ~ 20.

    The fixed edge_z=4 level sits at ~20% of the sigma_y-smoothed shoulder,
    which biases each wall outward by ~2 px/side (amplitude-dependent; the
    synthetic end-to-end study measures and corrects this bias as delta).
    Here we only pin coverage and a coarse width window.
    """
    rgb = _bright_fibre(noise=0.03)
    mr = compute_measurement(rgb, CONFIG(feature_mode="bright"), "synthetic")
    assert mr.meta["coverage"] >= 0.9
    med = float(np.nanmedian(mr.res.diameter_raw))
    assert abs(med - TRUE_W) < 12.0


def test_dim_fibre_end_to_end_bright_mode():
    """C1 dark end: fibre ~0.50 on bg ~0.20 with realistic noise (z ~ 8)."""
    rgb = _bright_fibre(fg_level=0.50, noise=0.03)
    mr = compute_measurement(rgb, CONFIG(feature_mode="bright"), "synthetic")
    assert mr.meta["coverage"] >= 0.9
    med = float(np.nanmedian(mr.res.diameter_raw))
    assert abs(med - TRUE_W) < 8.0


def test_bright_fibre_has_no_desat_signal():
    """Under the default desat mode the C1-style image must find ~nothing."""
    rgb = _bright_fibre()
    mr = compute_measurement(rgb, CONFIG(), "synthetic")
    assert mr.meta["coverage"] < 0.2


# --- CLI mode presets (study 03) -------------------------------------------

def _args(**kw):
    """Namespace with every build_config flag None except those given."""
    from fibrecv.run_measure import build_config
    fields = ["feature_mode", "ppu", "edge_z", "edge_frac", "edge_cap", "k_band",
              "min_width", "sigma_y", "wcol", "guard", "amin", "reject_dev",
              "margin", "min_coverage", "max_shift", "slope_min", "slope_rel",
              "rise_min"]
    ns = argparse.Namespace(**{f: None for f in fields})
    for k, v in kw.items():
        setattr(ns, k, v)
    return build_config(ns)


def test_desat_cli_keeps_dataclass_defaults():
    """No --feature-mode: the desat-calibrated CONFIG defaults are untouched."""
    cfg = _args()
    base = CONFIG()
    assert cfg.feature_mode == "desat"
    assert cfg.edge_frac == base.edge_frac == 0.65
    assert cfg.k_band == base.k_band == 4.0


def test_bright_cli_applies_calibrated_presets():
    """--feature-mode bright swaps in the study-03 calibrated knobs.

    Both are load-bearing: edge_frac is the PRIMARY relative threshold under
    the bright formula (not a faint-wall cap), and k_band must rise because the
    median z-map inflates z-scores -- at 4.0 the defocus halo balloons the
    coarse band and raises false band_mismatch.
    """
    cfg = _args(feature_mode="bright")
    assert cfg.feature_mode == "bright"
    assert cfg.edge_frac == 0.30
    assert cfg.k_band == 6.0


def test_explicit_flags_override_bright_presets():
    """An explicit flag always wins over the mode preset."""
    cfg = _args(feature_mode="bright", edge_frac=0.42, k_band=9.0)
    assert cfg.edge_frac == 0.42
    assert cfg.k_band == 9.0
