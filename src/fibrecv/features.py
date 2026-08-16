"""Desaturation feature map -- the load-bearing signal for fibre detection.

Dependencies
------------
``numpy``, ``skimage.color.rgb2hsv``.

Inputs
------
- ``rgb``: float RGB image in [0, 1], shape (H, W, 3).
- ``CONFIG`` for the margin fraction and numerical constants.

Output
------
- ``estimate_bg(S, cfg)`` -> ``(s_bg, mad)`` robust background saturation stats
  from the top+bottom margin rows.
- ``rgb_to_desaturation(rgb, cfg)`` -> ``(D, S, s_bg, mad)``, dispatching on
  ``cfg.feature_mode`` (unknown mode -> ValueError). ``"desat"``: ``D`` is the
  robust z-like desaturation map ``(s_bg - S) / (mad_scale*MAD + eps)``.
  ``"bright"``: same machinery on brightness ``V = rgb.max(axis=2)`` with the
  sign flipped (``(V - v_bg) / ...``) for bright-on-dark fibres. Either way
  ``D`` is large-positive inside the fibre, ~0 in background, and
  self-normalising per image so faint and dark-background cases are comparable.

Pos
---
Second stage of the per-image pipeline (after io_utils.load_rgb). Feeds
``band.py`` (coarse mask + centerline) and ``edges.py`` (per-column profiles).
Saturation is the discriminating feature for MasP2 (``"desat"``, calibrated);
brightness is for the multi-angle C1 set (``"bright"``).
"""

from __future__ import annotations

import numpy as np
from skimage.color import rgb2hsv

from .config import CONFIG


def estimate_bg(S: np.ndarray, cfg: CONFIG) -> tuple[float, float]:
    """Robust background saturation from the top and bottom margin rows.

    Uses the outer ``cfg.margin`` fraction of rows at top and bottom (assumed to
    be background, since the fibre runs roughly horizontally through the middle).
    Returns the median saturation ``s_bg`` and its MAD (median absolute
    deviation), both scalars.
    """
    h = S.shape[0]
    m = max(1, int(round(cfg.margin * h)))
    margin = np.concatenate([S[:m, :].ravel(), S[-m:, :].ravel()])
    s_bg = float(np.median(margin))
    mad = float(np.median(np.abs(margin - s_bg)))
    return s_bg, mad


def rgb_to_desaturation(rgb: np.ndarray, cfg: CONFIG) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Compute the feature z-map ``D`` from an RGB image (mode-dispatched).

    ``cfg.feature_mode == "desat"`` (default, MasP2): RGB->HSV, take
    saturation ``S``; estimate background ``(s_bg, mad)`` from the margins;
    build ``D = (s_bg - S) / (mad_scale*MAD + eps)``.

    ``cfg.feature_mode == "bright"`` (C1, bright fibre on dark bg): the same
    margin-row median/MAD machinery applied to brightness ``V = max(R,G,B)``,
    with the sign flipped: ``D = (V - v_bg) / (mad_scale*MAD + eps)``.

    Returns ``(D, F, f_bg, mad)`` where ``F`` is the feature channel (S or V)
    and ``f_bg`` its background level. NB downstream meta stores ``f_bg`` under
    the key ``"bg_S"`` even in bright mode (semantic overload; the active mode
    is recoverable from ``meta["params"]["feature_mode"]``). ``D`` is float32
    with the same H x W shape. Raises ``ValueError`` on an unknown mode.
    """
    if cfg.feature_mode == "desat":
        hsv = rgb2hsv(rgb)
        S = hsv[:, :, 1].astype(np.float32)
        s_bg, mad = estimate_bg(S, cfg)
        denom = cfg.mad_scale * mad + cfg.eps
        D = ((s_bg - S) / denom).astype(np.float32)
        return D, S, s_bg, mad
    if cfg.feature_mode == "bright":
        V = rgb.max(axis=2).astype(np.float32)
        v_bg, mad = estimate_bg(V, cfg)
        denom = cfg.mad_scale * mad + cfg.eps
        D = ((V - v_bg) / denom).astype(np.float32)
        return D, V, v_bg, mad
    raise ValueError(f"unknown feature_mode: {cfg.feature_mode!r}")
