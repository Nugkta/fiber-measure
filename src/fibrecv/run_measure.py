"""CLI: measure a glob of images in parallel (slow per-image stage).

Dependencies
------------
``numpy``, ``pandas`` (indirect), ``concurrent.futures.ProcessPoolExecutor``,
the ``measure`` + ``io_utils`` + ``config`` modules.

Inputs
------
Command-line flags: image selector (``--glob`` / ``--groups`` / ``--all``),
``--root`` (image dir), ``--out`` (output root), ``--jobs`` and every CONFIG
parameter (``--ppu``, ``--edge-frac`` strictness knob, ``--k-band`` ...).

Output
------
Per-image artifacts under ``<out>/overlays`` and ``<out>/per_image/...`` plus a
provenance snapshot ``<out>/summary/run_config.json`` and ``<out>/summary/run_log.txt``.

Pos
---
First entrypoint of the two-stage workflow. Heavy CPU stage -> run via ``srun``
on the HPC cluster, never on the login node. Its CSVs feed ``run_aggregate.py``.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
from pathlib import Path

from .config import CONFIG
from .io_utils import discover_images, parse_name
from .measure import measure_image

DEFAULT_ROOT = "/net/scratch/j56806hx/spins-cv/Images MasP2"
DEFAULT_OUT = "/net/scratch/j56806hx/spins-cv/output"


def build_config(args: argparse.Namespace) -> CONFIG:
    """Override CONFIG defaults from CLI flags (only those explicitly set)."""
    cfg = CONFIG()
    overrides = {
        "ppu": args.ppu,
        "edge_z": args.edge_z,
        "edge_frac": args.edge_frac,
        "k_band": args.k_band,
        "min_width": args.min_width,
        "sigma_y": args.sigma_y,
        "wcol": args.wcol,
        "guard": args.guard,
        "amin": args.amin,
        "reject_dev": args.reject_dev,
        "margin": args.margin,
        "min_coverage": args.min_coverage,
        "max_shift": args.max_shift,
        "slope_min": args.slope_min,
        "slope_rel": args.slope_rel,
        "rise_min": args.rise_min,
    }
    overrides = {k: v for k, v in overrides.items() if v is not None}
    return replace(cfg, **overrides)


def select_images(args: argparse.Namespace) -> list[Path]:
    """Resolve the image selector flags to a sorted, de-duplicated path list."""
    root = Path(args.root)
    if args.glob:
        return discover_images(root, args.glob)
    paths = discover_images(root)
    if args.groups:
        wanted = set(args.groups)
        kept: list[Path] = []
        for p in paths:
            try:
                group, _ = parse_name(p)
            except ValueError:
                continue
            if group in wanted:
                kept.append(p)
        return kept
    # default / --all
    return paths


def _worker(path_str: str, cfg: CONFIG, out_root: str) -> dict:
    """Pool worker: measure one image, capturing any error as a result dict."""
    try:
        return measure_image(path_str, cfg, out_root)
    except Exception as exc:  # noqa: BLE001 - report, never crash the pool
        return {"name": Path(path_str).stem, "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def _lib_versions() -> dict:
    import imageio
    import matplotlib
    import numpy
    import pandas
    import PIL
    import scipy
    import skimage
    return {
        "python": platform.python_version(),
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "scikit-image": skimage.__version__,
        "matplotlib": matplotlib.__version__,
        "pandas": pandas.__version__,
        "pillow": PIL.__version__,
        "imageio": imageio.__version__,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Measure fibre diameter profiles from images.")
    ap.add_argument("--root", default=DEFAULT_ROOT, help="image directory")
    ap.add_argument("--out", default=DEFAULT_OUT, help="output root directory")
    sel = ap.add_mutually_exclusive_group()
    sel.add_argument("--glob", help='image glob, e.g. "masp2 3_1_*.jpg"')
    sel.add_argument("--groups", nargs="+", help='A_B groups, e.g. 3_1 10_5')
    sel.add_argument("--all", action="store_true", help="all masp2 *_*.jpg images")
    ap.add_argument("--jobs", type=int, default=4, help="parallel processes")
    # CONFIG overrides (None -> keep calibrated default)
    ap.add_argument("--ppu", type=float, default=None)
    ap.add_argument("--edge-z", dest="edge_z", type=float, default=None,
                    help="strictness knob (absolute z above bg); higher -> tighter")
    ap.add_argument("--edge-frac", dest="edge_frac", type=float, default=None,
                    help="relative cap on the edge level for faint fibres")
    ap.add_argument("--k-band", dest="k_band", type=float, default=None)
    ap.add_argument("--min-width", dest="min_width", type=float, default=None)
    ap.add_argument("--sigma-y", dest="sigma_y", type=float, default=None)
    ap.add_argument("--wcol", type=int, default=None)
    ap.add_argument("--guard", type=int, default=None)
    ap.add_argument("--amin", type=float, default=None)
    ap.add_argument("--reject-dev", dest="reject_dev", type=float, default=None)
    ap.add_argument("--margin", type=float, default=None)
    ap.add_argument("--min-coverage", dest="min_coverage", type=float, default=None)
    ap.add_argument("--max-shift", dest="max_shift", type=int, default=None)
    ap.add_argument("--slope-min", dest="slope_min", type=float, default=None,
                    help="absolute wall-slope floor (z/px)")
    ap.add_argument("--slope-rel", dest="slope_rel", type=float, default=None,
                    help="wall slope as fraction of side max slope")
    ap.add_argument("--rise-min", dest="rise_min", type=float, default=None,
                    help="minimum z-rise of a wall run")
    args = ap.parse_args(argv)

    cfg = build_config(args)
    out_root = Path(args.out)
    (out_root / "summary").mkdir(parents=True, exist_ok=True)

    images = select_images(args)
    if not images:
        print("No images matched the selector.", file=sys.stderr)
        return 1

    # provenance snapshot
    with open(out_root / "summary" / "run_config.json", "w") as fh:
        json.dump({"params": cfg.as_dict(), "versions": _lib_versions(),
                   "n_images": len(images), "root": str(args.root)}, fh, indent=2)

    print(f"Measuring {len(images)} images with {args.jobs} jobs "
          f"(edge_z={cfg.edge_z}, wcol={cfg.wcol}) -> {out_root}")
    results = []
    if args.jobs <= 1:
        for p in images:
            results.append(_worker(str(p), cfg, str(out_root)))
            _echo(results[-1])
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(_worker, str(p), cfg, str(out_root)): p for p in images}
            for fut in as_completed(futs):
                r = fut.result()
                results.append(r)
                _echo(r)

    # run log
    errors = [r for r in results if "error" in r]
    low = [r for r in results if r.get("low_confidence")]
    log_path = out_root / "summary" / "run_log.txt"
    with open(log_path, "w") as fh:
        fh.write(f"images={len(images)} ok={len(results)-len(errors)} "
                 f"errors={len(errors)} low_confidence={len(low)}\n\n")
        for r in sorted(results, key=lambda d: d.get("name", "")):
            if "error" in r:
                fh.write(f"[ERROR] {r['name']}: {r['error']}\n")
            else:
                fh.write(f"{r['name']}: cov={r['coverage']:.0%} "
                         f"med={r['median_diameter_um']} "
                         f"tilt={r['tilt_slope']:.4f} "
                         f"{'LOWCONF' if r['low_confidence'] else ''}\n")
    print(f"Done. ok={len(results)-len(errors)} errors={len(errors)} "
          f"low_confidence={len(low)}. Log -> {log_path}")
    for r in errors:
        print(f"  [ERROR] {r['name']}: {r['error']}", file=sys.stderr)
    return 0


def _echo(r: dict) -> None:
    if "error" in r:
        print(f"  [ERROR] {r['name']}: {r['error']}")
    else:
        tag = " LOWCONF" if r.get("low_confidence") else ""
        print(f"  {r['name']}: cov={r['coverage']:.0%} med={r['median_diameter_um']}{tag}")


if __name__ == "__main__":
    raise SystemExit(main())
