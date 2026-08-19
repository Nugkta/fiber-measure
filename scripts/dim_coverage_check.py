"""Coverage vs k_band on the hardest C1 images, with a DATA-selected hard set.

Why this exists
---------------
The study's "dimmest / hardest images" coverage check had no saved script and
no stated selection rule, so it could not be told apart from an outcome-selected
set (picking the images that happened to behave). Here the hard set is chosen by
a statistic computed from the pixels alone, before any edge detection runs, and
the selection rule is recorded in the output.

Selection criterion (DATA-selected, not outcome-selected)
---------------------------------------------------------
For every C1 image, using only ``features.estimate_bg`` on the brightness
channel ``V = median(R, G, B)``:

    cnr = (percentile(V, 99) - v_bg) / (mad_scale * MAD_bg + eps)

i.e. how many robust background sigmas the brightest tissue sits above the
background level. It uses no band, no edges, no coverage and no diameter, so it
cannot be contaminated by the outcome being measured. The hard set is the
LOWEST QUARTILE of ``cnr`` (dimmest fibre against its own background noise).
``--criterion mad`` instead takes the HIGHEST quartile of the background MAD
(noisiest background); both are recorded so the choice is auditable.

Input
-----
- C1 images: ``/Users/stan/Documents/UOM/spins/multiangle/C1/C1_*_a*_part*.tiff``
  (the ``*s.tiff`` scale images and the ``_metadata.xml`` sidecars are excluded,
  matching the run's own file filter).

Output
------
``scripts/dim_coverage.json`` (``--out``): the criterion and its quartile cut,
the selected image list with scores, and for each ``k_band`` in [4,5,6,7,8]
(feature_mode=bright, edge_frac=0.30, edge_cap=0.50) the min / mean / median /
p05 coverage and the median width over that set -- plus a per-image
width-invariance summary across k, which is what "widths identical at every
k_band" actually needs to be checked against.

Pos
---
Study-03 reproducibility artifact. Makes the k_band robustness claim
reproducible and states the selection provenance explicitly.

Usage
-----
    uv run python scripts/dim_coverage_check.py
    uv run python scripts/dim_coverage_check.py --criterion mad --jobs 6

Reminder: once I am updated, update my header comments and the folder's md.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from fibrecv.compute import compute_measurement
from fibrecv.config import CONFIG
from fibrecv.features import estimate_bg
from fibrecv.io_utils import load_rgb

C1_ROOT = Path("/Users/stan/Documents/UOM/spins/multiangle/C1")
K_VALUES = [4.0, 5.0, 6.0, 7.0, 8.0]


def list_images(root: Path) -> list[Path]:
    """The 450 measurement images (same filter the run uses)."""
    return sorted(p for p in root.glob("C1_*_a*_part*.tiff")
                  if not p.name.endswith("s.tiff") and "_metadata" not in p.name)


def score_image(path: Path) -> dict:
    """Background-only difficulty score. No band, no edges, no QC."""
    cfg = CONFIG(feature_mode="bright")
    rgb = load_rgb(str(path))
    V = np.median(rgb, axis=2).astype(np.float32)
    v_bg, mad = estimate_bg(V, cfg)
    p99 = float(np.percentile(V, 99))
    cnr = (p99 - v_bg) / (cfg.mad_scale * mad + cfg.eps)
    return {"stem": path.stem, "path": str(path), "v_bg": float(v_bg),
            "mad_bg": float(mad), "p99_V": p99, "cnr": float(cnr)}


def measure_k(job: tuple[str, float]) -> dict:
    """Measure one image at one k_band with the shipped bright config."""
    path, k = job
    cfg = CONFIG(feature_mode="bright", edge_frac=0.30, edge_cap=0.50, k_band=k)
    rgb = load_rgb(path)
    mr = compute_measurement(rgb, cfg, Path(path).stem)
    valid = bool(mr.res.valid.any())
    return {
        "stem": Path(path).stem, "k_band": k,
        "coverage": float(mr.meta["coverage"]),
        "median_width_px": float(np.nanmedian(mr.res.diameter_raw)) if valid else float("nan"),
        "low_confidence": bool(mr.meta["low_confidence"]),
        "band_mismatch": bool(mr.meta["band_mismatch"]),
        "band_half_px": float(mr.bnd.band_half),
    }


def _agg(vals: list[float]) -> dict:
    a = np.asarray([v for v in vals if np.isfinite(v)], dtype=float)
    if a.size == 0:
        return {"n": 0}
    return {"n": int(a.size), "min": float(a.min()), "mean": float(a.mean()),
            "median": float(np.median(a)), "p05": float(np.percentile(a, 5)),
            "max": float(a.max())}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", type=Path, default=C1_ROOT)
    ap.add_argument("--out", type=Path, default=Path("scripts/dim_coverage.json"))
    ap.add_argument("--criterion", choices=["cnr", "mad"], default="cnr")
    ap.add_argument("--quantile", type=float, default=0.25)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--max-images", type=int, default=None,
                    help="cap the hard set (kept deterministic: the hardest N)")
    args = ap.parse_args(argv)

    rev = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                         text=True).stdout.strip()

    imgs = list_images(args.root)
    print(f"scoring {len(imgs)} images ({args.jobs} jobs)...")
    sys.stdout.flush()
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        scores = list(ex.map(score_image, imgs, chunksize=4))

    if args.criterion == "cnr":
        key, ascending = "cnr", True          # lowest quartile of contrast
        cut = float(np.quantile([s["cnr"] for s in scores], args.quantile))
        hard = [s for s in scores if s["cnr"] <= cut]
        rule = (f"lowest {args.quantile:.0%} quantile of cnr "
                f"= (p99(V) - v_bg) / (1.4826*MAD_bg); cut = {cut:.4f}")
    else:
        key, ascending = "mad_bg", False      # highest quartile of background MAD
        cut = float(np.quantile([s["mad_bg"] for s in scores], 1 - args.quantile))
        hard = [s for s in scores if s["mad_bg"] >= cut]
        rule = (f"highest {args.quantile:.0%} quantile of background MAD; "
                f"cut = {cut:.6f}")

    hard.sort(key=lambda s: s[key], reverse=not ascending)
    if args.max_images:
        hard = hard[:args.max_images]
    print(f"hard set: {len(hard)} / {len(imgs)} images  ({rule})")
    sys.stdout.flush()

    jobs = [(s["path"], k) for k in K_VALUES for s in hard]
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        rows = list(ex.map(measure_k, jobs, chunksize=2))

    by_k = {}
    print(f"\n{'k_band':>7} {'cov_min':>8} {'cov_mean':>9} {'cov_p05':>8} "
          f"{'w_med':>8} {'lowconf':>8} {'mismatch':>9}")
    for k in K_VALUES:
        sel = [r for r in rows if r["k_band"] == k]
        cov = _agg([r["coverage"] for r in sel])
        wid = _agg([r["median_width_px"] for r in sel])
        by_k[str(k)] = {
            "n_images": len(sel),
            "coverage": cov,
            "median_width_px": wid,
            "n_low_confidence": int(sum(r["low_confidence"] for r in sel)),
            "n_band_mismatch": int(sum(r["band_mismatch"] for r in sel)),
        }
        print(f"{k:>7.1f} {cov['min']:>8.4f} {cov['mean']:>9.4f} {cov['p05']:>8.4f} "
              f"{wid['median']:>8.3f} {by_k[str(k)]['n_low_confidence']:>8d} "
              f"{by_k[str(k)]['n_band_mismatch']:>9d}")

    # ---- per-image width invariance across k (tests "identical at every k_band") ----
    per_img = {}
    for r in rows:
        per_img.setdefault(r["stem"], {})[r["k_band"]] = r["median_width_px"]
    spreads, cov_spreads = [], []
    for stem, d in per_img.items():
        w = np.array([d[k] for k in K_VALUES], dtype=float)
        if np.isfinite(w).all():
            spreads.append((stem, float(w.max() - w.min())))
    cov_by_img = {}
    for r in rows:
        cov_by_img.setdefault(r["stem"], {})[r["k_band"]] = r["coverage"]
    for stem, d in cov_by_img.items():
        c = np.array([d[k] for k in K_VALUES], dtype=float)
        cov_spreads.append((stem, float(c.max() - c.min())))
    sv = np.array([s for _, s in spreads]) if spreads else np.array([np.nan])
    cv = np.array([s for _, s in cov_spreads])
    worst_w = max(spreads, key=lambda t: t[1]) if spreads else (None, float("nan"))
    worst_c = max(cov_spreads, key=lambda t: t[1])

    invariance = {
        "definition": "per image, (max over k of median width) - (min over k of median width)",
        "n_images": int(sv.size),
        "width_spread_px": {"median": float(np.median(sv)), "mean": float(sv.mean()),
                            "p95": float(np.percentile(sv, 95)), "max": float(sv.max())},
        "worst_image_width": {"stem": worst_w[0], "spread_px": worst_w[1]},
        "n_images_width_exactly_identical": int((sv == 0.0).sum()),
        "coverage_spread": {"median": float(np.median(cv)), "mean": float(cv.mean()),
                            "p95": float(np.percentile(cv, 95)), "max": float(cv.max())},
        "worst_image_coverage": {"stem": worst_c[0], "spread": worst_c[1]},
    }
    print(f"\nwidth spread across k: median={invariance['width_spread_px']['median']:.4f} "
          f"p95={invariance['width_spread_px']['p95']:.4f} "
          f"max={invariance['width_spread_px']['max']:.4f} px "
          f"({invariance['n_images_width_exactly_identical']}/{invariance['n_images']} "
          f"exactly identical)")

    out = {
        "kind": "dim_coverage_check",
        "git_rev": rev,
        "root": str(args.root),
        "n_images_total": len(imgs),
        "selection": {
            "outcome_selected": False,
            "data_selected": True,
            "criterion": args.criterion,
            "rule": rule,
            "cut_value": cut,
            "computed_from": "pixels only: V = median(R,G,B), background level and "
                             "MAD from features.estimate_bg on the top/bottom margin "
                             "rows, plus percentile(V, 99). No band, no edge "
                             "detection, no coverage, no diameter enters the score.",
            "n_selected": len(hard),
        },
        "measure_config": "CONFIG(feature_mode='bright', edge_frac=0.30, "
                          "edge_cap=0.50, k_band=<swept>)",
        "k_values": K_VALUES,
        "by_k": by_k,
        "width_invariance_across_k": invariance,
        "hard_set": hard,
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
