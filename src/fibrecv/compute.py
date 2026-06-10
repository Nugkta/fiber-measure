"""Pure per-image compute core: features -> band -> edges -> qc, no file I/O.

Dependencies
------------
``numpy`` plus the pipeline modules ``io_utils``, ``features``, ``band``,
``edges``, ``qc``. Deliberately imports nothing that writes to disk.

Inputs
------
- ``rgb``: float RGB image in [0, 1], shape (H, W, 3) (already decoded).
- A ``CONFIG`` carrying every tunable parameter.
- Optional ``name`` (the image stem, e.g. ``"masp2 3_1_2"``) used only to parse
  the ``group``/``replicate`` labels and stamp the meta dict.

Output
------
``compute_measurement(rgb, cfg, name=None)`` -> ``MeasureResult`` bundling the
desaturation map ``D``, the ``BandResult``/``EdgeResult``/``QCResult`` objects,
the per-column ``diameter_um`` array, the parsed ``name``/``group``/``replicate``
and the diagnostics ``meta`` dict -- all in memory, nothing written.

Pos
---
The shared compute heart. ``measure.measure_image`` calls it and then writes the
CSV/plot/overlay/meta artifacts; the Streamlit GUI calls it to redraw boundaries
on every parameter change without touching disk. Splitting compute from
artifact-writing is what makes both paths possible from one code path.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from . import band as band_mod
from . import edges as edges_mod
from . import features as feat_mod
from . import io_utils
from . import qc as qc_mod
from .band import BandResult
from .config import CONFIG
from .edges import EdgeResult
from .qc import QCResult


@dataclass
class MeasureResult:
    """In-memory result of measuring one image (no artifacts written)."""

    rgb: np.ndarray            # original float RGB in [0, 1], (H, W, 3)
    D: np.ndarray             # desaturation z-map, (H, W)
    bnd: BandResult           # band localisation + centerline
    edg: EdgeResult           # per-column boundaries
    res: QCResult             # cleaned profile + coverage/flags
    diameter_um: np.ndarray   # per-column diameter in microns, NaN where invalid
    name: str | None          # image stem, or None for ad-hoc (uploaded) images
    group: str | None         # parsed "A_B" group, or None if name unparseable
    replicate: int | None     # parsed replicate C, or None if name unparseable
    meta: dict = field(default_factory=dict)  # diagnostics (identical to the meta JSON)


def compute_measurement(rgb: np.ndarray, cfg: CONFIG, name: str | None = None) -> MeasureResult:
    """Run features -> band -> edges -> qc in memory and assemble the meta dict.

    This is the exact computation that used to live inline in
    ``measure.measure_image`` (lines ~68-123), lifted out verbatim so the CLI and
    the GUI share one code path. No files are read or written here.
    """
    H, W = rgb.shape[:2]
    D, S, s_bg, mad = feat_mod.rgb_to_desaturation(rgb, cfg)
    bnd = band_mod.locate_band(D, cfg)
    edg = edges_mod.detect_edges(D, bnd, cfg)
    res = qc_mod.run_qc(edg, bnd, cfg)

    diameter_um = np.where(res.valid, res.diameter_raw / cfg.ppu, np.nan)

    group: str | None = None
    replicate: int | None = None
    if name is not None:
        try:
            group, replicate = io_utils.parse_name(name)
        except ValueError:
            group, replicate = None, None

    # diagnostics meta -- identical content/order to the JSON written by the CLI
    span = slice(bnd.x0, bnd.x1 + 1)
    flag_counts = Counter(int(f) for f in res.reason[span])
    meta = {
        "name": name,
        "group": group,
        "replicate": replicate,
        "image_shape": [int(H), int(W)],
        "bg_S": s_bg,
        "MAD": mad,
        "tilt_slope": bnd.slope,
        "band_half_px": bnd.band_half,
        "band_span": [int(bnd.x0), int(bnd.x1)],
        "half_window_px": int(edg.half_window),
        "n_components": int(bnd.n_components),
        "coverage": res.coverage,
        "n_valid": int(res.valid.sum()),
        "n_span": int(bnd.x1 - bnd.x0 + 1),
        "low_confidence": bool(res.low_confidence),
        "band_mismatch": bool(res.band_mismatch),
        "flag_counts": {str(k): int(v) for k, v in flag_counts.items()},
        "median_diameter_um": float(np.nanmedian(diameter_um)) if res.valid.any() else None,
        "params": cfg.as_dict(),
    }

    return MeasureResult(
        rgb=rgb,
        D=D,
        bnd=bnd,
        edg=edg,
        res=res,
        diameter_um=diameter_um,
        name=name,
        group=group,
        replicate=replicate,
        meta=meta,
    )
