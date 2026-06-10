"""Per-image orchestration: features -> band -> edges -> qc -> artifacts.

Dependencies
------------
``numpy``, ``pandas``, ``matplotlib`` (Agg backend), ``io_utils`` (load/parse),
``compute`` (the pure detection core) and ``overlay`` (boundary PNG). The
features/band/edges/qc stages are reached via ``compute.compute_measurement``.

Inputs
------
- One image path (``masp2 A_B_C.jpg``).
- A ``CONFIG`` and the output root directory.

Output
------
Writes four artifacts per image and returns a small summary dict:
  - ``overlays/<name>_overlay.png``        boundary check image
  - ``per_image/csv/<name>_profile.csv``   x_px, diameter_px_raw/smooth, diameter_um, valid, interpolated
  - ``per_image/plots/<name>_profile.png`` diameter-vs-position plot
  - ``per_image/diagnostics/<name>_meta.json`` bg_S, MAD, params, coverage, tilt, flags

Pos
---
The heart of the per-image pipeline; called by ``run_measure.py`` (one call per
image, parallelised). Its CSV output is the input contract for ``run_aggregate``.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from . import io_utils  # noqa: E402
from .compute import compute_measurement  # noqa: E402
from .config import CONFIG  # noqa: E402
from .overlay import draw_overlay  # noqa: E402


def _out_paths(out_root: Path, name: str) -> dict[str, Path]:
    return {
        "overlay": out_root / "overlays" / f"{name}_overlay.png",
        "csv": out_root / "per_image" / "csv" / f"{name}_profile.csv",
        "plot": out_root / "per_image" / "plots" / f"{name}_profile.png",
        "meta": out_root / "per_image" / "diagnostics" / f"{name}_meta.json",
    }


def measure_image(path: str | Path, cfg: CONFIG, out_root: str | Path) -> dict:
    """Measure one image end-to-end and write its CSV/plot/overlay/meta.

    Computation is delegated to ``compute.compute_measurement`` (pure, no I/O)
    and artifact-writing to ``write_measurement``; this function only resolves the
    path and loads the image, so its outputs are byte-identical to the
    pre-refactor implementation.
    """
    path = Path(path)
    rgb = io_utils.load_rgb(path)
    mr = compute_measurement(rgb, cfg, name=path.stem)
    return write_measurement(rgb, mr, cfg, out_root)


def write_measurement(rgb, mr, cfg: CONFIG, out_root: str | Path) -> dict:
    """Write the four per-image artifacts for an already-computed measurement.

    Shared by ``measure_image`` (CLI) and the GUI's group export so both write
    byte-identical CSV/plot/overlay/meta. ``mr.name`` supplies the output stem.
    """
    out_root = Path(out_root)
    name = mr.name
    paths = _out_paths(out_root, name)
    for p in paths.values():
        p.parent.mkdir(parents=True, exist_ok=True)

    bnd, edg, res = mr.bnd, mr.edg, mr.res
    x = np.arange(rgb.shape[1])

    # --- per-image CSV (restricted to the band span for compactness) ---
    span = slice(bnd.x0, bnd.x1 + 1)
    df = pd.DataFrame(
        {
            "x_px": x[span],
            "diameter_px_raw": res.diameter_raw[span],
            "diameter_px_smooth": res.diameter_smooth[span],
            "diameter_um": mr.diameter_um[span],
            "valid": res.valid[span].astype(int),
            "interpolated": res.interpolated[span].astype(int),
        }
    )
    df.to_csv(paths["csv"], index=False)

    # --- overlay PNG ---
    draw_overlay(
        rgb, edg.y_top, edg.y_bot, bnd.c_fit, res.valid,
        paths["overlay"], bnd.x0, bnd.x1, thick=1,
    )

    # --- profile plot ---
    _plot_profile(df, name, res.coverage, cfg, paths["plot"])

    # --- diagnostics meta (assembled in compute_measurement) ---
    with open(paths["meta"], "w") as fh:
        json.dump(mr.meta, fh, indent=2)

    return {
        "name": name,
        "group": mr.group,
        "replicate": mr.replicate,
        "coverage": res.coverage,
        "low_confidence": res.low_confidence,
        "median_diameter_um": mr.meta["median_diameter_um"],
        "tilt_slope": bnd.slope,
    }


def _plot_profile(df: pd.DataFrame, name: str, coverage: float, cfg: CONFIG, out: Path) -> None:
    """Diameter-vs-position scatter (raw) + smoothed line, in microns."""
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(df["x_px"], df["diameter_um"], ".", ms=2, alpha=0.4, color="tab:gray", label="raw")
    sm_um = df["diameter_px_smooth"] / cfg.ppu
    ax.plot(df["x_px"], sm_um, "-", lw=1.3, color="tab:blue", label="smooth")
    ax.set_xlabel("x position (px)")
    ax.set_ylabel("diameter (µm)")
    ax.set_title(f"{name}   coverage={coverage:.0%}")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)
