"""Synthetic control WITH a chromatic fringe -- the median-vs-max discriminator.

Why this exists
---------------
The pre-existing synthetic control (``scripts/calibrate_edge_frac.py::_render_circle``)
paints the *same* geometry into R, G and B. It therefore contains no chromatic
fringe at all, so ``max(R,G,B)`` and ``median(R,G,B)`` see the same feature map
and the control cannot test the z-map change -- it only tests the threshold
change. (It is also 500x1 px: ``W_IMG=1200`` is declared but never used, so the
delta is measured on a *single column*.) This control fixes both: full-width
images, and the R channel's edge displaced outward by ``disp`` px to emulate
lateral chromatic aberration.

Input
-----
Nothing on disk. ``fibrecv.compute`` / ``fibrecv.config`` / ``fibrecv.features``
from the worktree whose ``.venv`` is active.

Output
------
``scripts/synthetic_fringe.json`` (``--out``): every (displacement, fringe
variant, config) cell with median measured width, per-side delta, coverage,
plus the mechanism-attribution deltas.

How the OLD pipeline is emulated WITHOUT editing ``src``
-------------------------------------------------------
Two independent knobs are needed: the z-map, and the edge-level formula.

1. **z-map** -- ``fibrecv.features.rgb_to_desaturation`` is monkeypatched (on
   the module object that ``compute`` looks it up on) with a function that
   builds ``V`` as either ``max(R,G,B)`` (old) or ``median(R,G,B)`` (new), using
   the identical background/MAD machinery. ``cfg.feature_mode`` no longer
   selects the feature; the ``--zmap`` argument does.

2. **edge-level formula** -- because the patch decouples the feature from
   ``cfg.feature_mode``, ``feature_mode`` is free to select the *formula* in
   ``edges._side_edge``. The base branch used ``level = base + min(edge_z,
   edge_frac*A)`` unconditionally, which is byte-for-byte the surviving
   ``desat`` branch on HEAD. So ``feature_mode="desat"`` + ``edge_z=4.0`` +
   ``edge_frac=0.65`` + patched max z-map reproduces the old pipeline exactly,
   and ``feature_mode="bright"`` selects the new clamped formula.

``--native`` skips the patch entirely and runs ``CONFIG(feature_mode="bright",
edge_z=4.0, edge_frac=0.65, k_band=4.0)``. Run with ``--native`` from the
*base* worktree and the numbers must match config (i) measured here; that is
the emulation's audit trail.

Pos
---
Study-03 reproducibility artifact. Answers: does a chromatic fringe actually
separate median from max, and how much of the reported synthetic improvement is
the z-map versus the threshold?

Usage
-----
    uv run python scripts/synthetic_fringe_control.py
    (cd <base worktree> && uv run python <this file> --native --out /tmp/native.json)

Reminder: once I am updated, update my header comments and the folder's md.
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np

import fibrecv.features as feat_mod
from fibrecv.compute import compute_measurement
from fibrecv.config import CONFIG
from fibrecv.features import estimate_bg

H_IMG = 500
W_IMG = 1200
CIRCLE_D = 190.0          # TRUE full width, set by the un-displaced G/B edges
NOISE = 0.03
SEED = 200

_ORIG_FEATURE = feat_mod.rgb_to_desaturation
_ZMAP = {"mode": "median"}


def _patched_feature(rgb: np.ndarray, cfg: CONFIG):
    """Brightness z-map with the R/G/B reducer chosen by ``_ZMAP`` (not cfg)."""
    if _ZMAP["mode"] == "max":
        V = rgb.max(axis=2).astype(np.float32)
    elif _ZMAP["mode"] == "median":
        V = np.median(rgb, axis=2).astype(np.float32)
    else:                                     # pragma: no cover
        raise ValueError(_ZMAP["mode"])
    v_bg, mad = estimate_bg(V, cfg)
    D = ((V - v_bg) / (cfg.mad_scale * mad + cfg.eps)).astype(np.float32)
    return D, V, v_bg, mad


def render_fringe(disp: float = 0.0, fringe_channels: tuple[int, ...] = (0,),
                  seed: int = SEED, noise: float = NOISE,
                  width: int = W_IMG) -> np.ndarray:
    """Horizontal circular fibre whose ``fringe_channels`` edges sit ``disp`` px out.

    Same profile as ``calibrate_edge_frac._render_circle`` (a 2 px linear ramp
    whose 50% point is the true boundary), but (a) full width and (b) the
    listed channels get a half-width of ``h + disp``, so the outer ``disp`` px
    ring is bright in those channels only -- exactly what lateral chromatic
    aberration produces. ``disp=0`` reproduces the achromatic control.
    """
    h = CIRCLE_D / 2.0
    y = np.arange(H_IMG, dtype=float)[:, None]
    dist = np.abs(y - H_IMG / 2.0)
    bg = np.array([0.21, 0.20, 0.19])
    fg = np.array([0.85, 0.84, 0.83])
    img = np.empty((H_IMG, width, 3), dtype=float)
    for c in range(3):
        hc = h + (disp if c in fringe_channels else 0.0)
        t = np.clip((hc + 1.0 - dist) / 2.0, 0.0, 1.0)        # (H, 1)
        img[:, :, c] = bg[c] + (fg[c] - bg[c]) * t
    rng = np.random.default_rng(seed)
    img = img + rng.normal(0.0, noise, img.shape)
    return np.clip(img, 0.0, 1.0).astype(np.float32)


def measure(rgb: np.ndarray, cfg: CONFIG, zmap: str | None) -> dict:
    """Measure one image; ``zmap=None`` -> no patch (native package behaviour)."""
    if zmap is None:
        mr = compute_measurement(rgb, cfg, "fringe")
    else:
        _ZMAP["mode"] = zmap
        feat_mod.rgb_to_desaturation = _patched_feature
        try:
            mr = compute_measurement(rgb, cfg, "fringe")
        finally:
            feat_mod.rgb_to_desaturation = _ORIG_FEATURE
    med = float(np.nanmedian(mr.res.diameter_raw))
    return {
        "median_width_px": med,
        "delta_per_side_px": (med - CIRCLE_D) / 2.0,
        "coverage": float(mr.meta["coverage"]),
        "band_half_px": float(mr.bnd.band_half),
    }


# label -> (zmap, CONFIG, human description)
def build_configs() -> dict[str, tuple[str, CONFIG, str]]:
    return {
        "i_old_max_ez4": (
            "max",
            CONFIG(feature_mode="desat", edge_z=4.0, edge_frac=0.65, k_band=4.0),
            "OLD pipeline: max(R,G,B) z-map, level = bg + min(edge_z=4, 0.65*A), k_band=4",
        ),
        "ii_median_frac065": (
            "median",
            CONFIG(feature_mode="bright", edge_z=4.0, edge_frac=0.65,
                   edge_cap=0.50, k_band=4.0),
            "median z-map, dataclass-default edge_frac=0.65 (effective 0.50 after edge_cap), k_band=4",
        ),
        "iii_median_frac030_shipped": (
            "median",
            CONFIG(feature_mode="bright", edge_z=4.0, edge_frac=0.30,
                   edge_cap=0.50, k_band=6.0),
            "SHIPPED CLI config: median z-map, edge_frac=0.30, edge_cap=0.50, k_band=6",
        ),
        "iv_max_frac030": (
            "max",
            CONFIG(feature_mode="bright", edge_z=4.0, edge_frac=0.30,
                   edge_cap=0.50, k_band=6.0),
            "max z-map with the NEW threshold (isolates z-map from threshold vs iii)",
        ),
        "ii_median_frac065_k6": (
            "median",
            CONFIG(feature_mode="bright", edge_z=4.0, edge_frac=0.65,
                   edge_cap=0.50, k_band=6.0),
            "control: same as (ii) but k_band=6 -- shows k_band is inert on synthetic",
        ),
    }


VARIANTS = {"R_only": (0,), "R_and_B": (0, 2)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=Path("scripts/synthetic_fringe.json"))
    ap.add_argument("--disps", type=float, nargs="+", default=[0.0, 1.0, 2.0, 3.0])
    ap.add_argument("--widths", type=int, nargs="+", default=[1200, 1],
                    help="1200 = proper image; 1 = the shipped calibrate_edge_frac geometry")
    ap.add_argument("--native", action="store_true",
                    help="no monkeypatch: run CONFIG(feature_mode='bright', edge_z=4, "
                         "edge_frac=0.65, k_band=4) as the package defines it "
                         "(on the BASE worktree this IS the old pipeline)")
    ap.add_argument("--validate-native", type=Path, default=None,
                    metavar="NATIVE_JSON",
                    help="a --native record produced on the BASE worktree; its "
                         "rows are compared cell-by-cell against config (i) to "
                         "prove the monkeypatch emulation is exact")
    args = ap.parse_args(argv)

    rev = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                         text=True).stdout.strip()
    branch = subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"],
                            capture_output=True, text=True).stdout.strip()

    if args.native:
        cfgs = {"native_bright_ez4_frac065_k4": (
            None, CONFIG(feature_mode="bright", edge_z=4.0, edge_frac=0.65, k_band=4.0),
            "native package behaviour, no monkeypatch")}
    else:
        cfgs = build_configs()

    rows = []
    print(f"{'width':>6} {'variant':>9} {'disp':>5} {'config':>28} "
          f"{'width_px':>9} {'delta':>8} {'cov':>6}")
    print("-" * 82)
    for width in args.widths:
        for vname, chans in VARIANTS.items():
            for disp in args.disps:
                if disp == 0.0 and vname != "R_only":
                    continue          # disp=0 is identical for every variant
                rgb = render_fringe(disp=disp, fringe_channels=chans, width=width)
                for label, (zmap, cfg, desc) in cfgs.items():
                    r = measure(rgb, cfg, zmap)
                    row = {"width_px_image": width, "fringe_variant": vname,
                           "fringe_channels": list(chans), "disp_px": disp,
                           "config": label, "config_desc": desc,
                           "zmap": zmap or "native", **r}
                    rows.append(row)
                    print(f"{width:>6d} {vname:>9} {disp:>5.1f} {label:>28} "
                          f"{r['median_width_px']:>9.3f} "
                          f"{r['delta_per_side_px']:>+8.3f} {r['coverage']:>6.3f}")

    # ---- mechanism attribution (full-width images, R_only fringe) ----
    def get(width, variant, disp, label):
        for r in rows:
            if (r["width_px_image"] == width and r["fringe_variant"] == variant
                    and r["disp_px"] == disp and r["config"] == label):
                return r["delta_per_side_px"]
        return None

    attribution = []
    if not args.native:
        for width in args.widths:
            for vname in VARIANTS:
                for disp in args.disps:
                    v = vname if disp != 0.0 else "R_only"
                    if disp == 0.0 and vname != "R_only":
                        continue
                    d_old = get(width, v, disp, "i_old_max_ez4")
                    d_ship = get(width, v, disp, "iii_median_frac030_shipped")
                    d_max30 = get(width, v, disp, "iv_max_frac030")
                    d_065 = get(width, v, disp, "ii_median_frac065")
                    if None in (d_old, d_ship, d_max30, d_065):
                        continue
                    attribution.append({
                        "width_px_image": width, "fringe_variant": v, "disp_px": disp,
                        "delta_old_i": d_old,
                        "delta_shipped_iii": d_ship,
                        "delta_max_frac030_iv": d_max30,
                        "delta_median_frac065_ii": d_065,
                        "total_change_i_to_iii": d_ship - d_old,
                        "zmap_only_iv_minus_iii": d_max30 - d_ship,
                        "threshold_only_i_minus_iv": d_max30 - d_old,
                        "fringe_separates_median_from_max_px": d_max30 - d_ship,
                    })

    # ---- emulation audit: config (i) here vs the base branch running natively ----
    emulation_check = None
    if args.validate_native and not args.native:
        nat = json.loads(args.validate_native.read_text())
        cells, worst, all_exact = [], 0.0, True
        for nr in nat["rows"]:
            mine = get(nr["width_px_image"], nr["fringe_variant"], nr["disp_px"],
                       "i_old_max_ez4")
            if mine is None:
                continue
            d = abs(mine - nr["delta_per_side_px"])
            worst = max(worst, d)
            all_exact &= (mine == nr["delta_per_side_px"])
            cells.append({"width_px_image": nr["width_px_image"],
                          "fringe_variant": nr["fringe_variant"],
                          "disp_px": nr["disp_px"],
                          "emulated_delta_i": mine,
                          "native_base_delta": nr["delta_per_side_px"],
                          "abs_diff_px": d})
        emulation_check = {
            "native_source": str(args.validate_native),
            "native_git_branch": nat.get("git_branch"),
            "native_git_rev": nat.get("git_rev"),
            "n_cells": len(cells),
            "max_abs_diff_px": worst,
            "bitwise_equal_all_cells": bool(all_exact),
            "verdict": ("config (i) reproduces the base branch EXACTLY"
                        if all_exact else
                        f"config (i) differs from the base branch by up to {worst:g} px"),
            "cells": cells,
        }
        print(f"\nemulation check vs base branch: {emulation_check['verdict']}")

    out = {
        "kind": "synthetic_fringe_control",
        "emulation_check": emulation_check,
        "git_branch": branch, "git_rev": rev,
        "native_mode": bool(args.native),
        "true_width_px": CIRCLE_D,
        "noise_sigma": NOISE, "seed": SEED,
        "notes": [
            "delta_per_side_px = (median measured width - 190.0) / 2.",
            "The R channel's edge is displaced OUTWARD by disp px; the true "
            "boundary stays at 190 px (set by the un-displaced G/B edges).",
            "config (i) emulates the base branch: patched max(R,G,B) z-map + "
            "feature_mode='desat', whose level formula min(edge_z, edge_frac*A) "
            "is byte-identical to the base branch's unconditional formula.",
            "width_px_image=1 reproduces the geometry actually used by "
            "scripts/calibrate_edge_frac.py (W_IMG=1200 is declared but unused "
            "there, so its images are 500x1).",
        ],
        "configs": {k: v[2] for k, v in cfgs.items()},
        "rows": rows,
        "attribution": attribution,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
