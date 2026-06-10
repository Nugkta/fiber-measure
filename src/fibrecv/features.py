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
- ``rgb_to_desaturation(rgb, cfg)`` -> ``(D, S, s_bg, mad)`` where ``D`` is the
  robust z-like desaturation map ``(s_bg - S) / (mad_scale*MAD + eps)``:
  large-positive inside the (desaturated) fibre, ~0 in background, and
  self-normalising per image so faint and dark-background cases are comparable.

Pos
---
Second stage of the per-image pipeline (after io_utils.load_rgb). Feeds
``band.py`` (coarse mask + centerline) and ``edges.py`` (per-column profiles).
Saturation -- not brightness -- is the discriminating feature (calibrated).
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
    """Compute the desaturation z-map ``D`` from an RGB image.

    Steps: RGB->HSV, take saturation ``S``; estimate background ``(s_bg, mad)``
    from the margins; build ``D = (s_bg - S) / (mad_scale*MAD + eps)``.

    Returns ``(D, S, s_bg, mad)``. ``D`` is float32 with the same H x W shape.
    """
    hsv = rgb2hsv(rgb)
    S = hsv[:, :, 1].astype(np.float32)
    s_bg, mad = estimate_bg(S, cfg)
    denom = cfg.mad_scale * mad + cfg.eps
    D = ((s_bg - S) / denom).astype(np.float32)
    return D, S, s_bg, mad
