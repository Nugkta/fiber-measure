"""How stable is the synthetic per-side delta across noise seeds?

Why this exists
---------------
The study quotes the synthetic edge bias to three significant figures
("3.90 -> -0.51 px") from a SINGLE noise realisation. A single draw carries no
error bar, so this re-measures delta over many seeds and reports mean / sd /
min / max, on both synthetic controls that exist in the repo.

Input
-----
- ``scripts/calibrate_edge_frac.py::_render_circle`` -- the single-image control
  the study actually used (note: it renders a 500x1 image; ``W_IMG=1200`` is
  declared but never used there).
- ``scripts/synthetic_fringe_control.py::render_fringe`` -- the same circle at
  full width, for a like-for-like full-image comparison.
- ``tests/test_xsection_synthetic.py::_render`` + ``fibrecv.xsection`` -- the
  6-angle stacked control, whose ``delta`` fixture is where "3.90 / -0.51" comes
  from.

Output
------
``scripts/delta_seeds.json`` (``--out``): per (control, config) the full list of
per-seed deltas plus mean, sd (ddof=1), min, max, median and range.

Configs measured
----------------
- ``shipped``   feature_mode=bright, edge_frac=0.30, edge_cap=0.50, k_band=6.0
- ``default``   feature_mode=bright, dataclass defaults (edge_frac=0.65,
                edge_cap=0.50, k_band=4.0) -- ``CONFIG(feature_mode="bright")``
- ``old``       the base-branch pipeline, emulated exactly via the monkeypatch
                documented in ``synthetic_fringe_control`` (max(R,G,B) z-map +
                the ``min(edge_z, edge_frac*A)`` level formula). Included so the
                "before" number carries an error bar too.

Pos
---
Study-03 reproducibility artifact. Puts a dispersion on every quoted synthetic
delta so the report can stop reporting 3 s.f. from one draw.

Usage
-----
    uv run python scripts/delta_seed_stability.py
    uv run python scripts/delta_seed_stability.py --n-seeds 20 --skip-stack

Reminder: once I am updated, update my header comments and the folder's md.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from fibrecv.compute import compute_measurement
from fibrecv.config import CONFIG
from fibrecv.xsection import NOMINAL_ANGLES_DEG, build_part_stack

SCRIPTS = Path(__file__).resolve().parent
REPO = SCRIPTS.parent

sys.path.insert(0, str(SCRIPTS))
from calibrate_edge_frac import CIRCLE_D, _render_circle          # noqa: E402
from synthetic_fringe_control import (                             # noqa: E402
    _ORIG_FEATURE, _ZMAP, _patched_feature, render_fringe,
)
import fibrecv.features as feat_mod                                # noqa: E402


def _load_test_module():
    """Import ``tests/test_xsection_synthetic.py`` by path (it is not a package)."""
    path = REPO / "tests" / "test_xsection_synthetic.py"
    spec = importlib.util.spec_from_file_location("_test_xsec_synth", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# label -> (zmap or None, CONFIG, description)
CONFIGS = {
    "shipped": ("median",
                CONFIG(feature_mode="bright", edge_z=4.0, edge_frac=0.30,
                       edge_cap=0.50, k_band=6.0),
                "SHIPPED CLI config: median z-map, edge_frac=0.30, edge_cap=0.50, k_band=6.0"),
    "default": ("median",
                CONFIG(feature_mode="bright"),
                "dataclass default: median z-map, edge_frac=0.65 (effective 0.50 "
                "after edge_cap), k_band=4.0"),
    "old": ("max",
            CONFIG(feature_mode="desat", edge_z=4.0, edge_frac=0.65, k_band=4.0),
            "OLD base-branch pipeline, emulated: max(R,G,B) z-map + "
            "level = bg + min(edge_z=4, 0.65*A), k_band=4.0"),
}


def _with_zmap(zmap: str, fn):
    """Run ``fn()`` with the feature reducer forced to ``zmap``."""
    _ZMAP["mode"] = zmap
    feat_mod.rgb_to_desaturation = _patched_feature
    try:
        return fn()
    finally:
        feat_mod.rgb_to_desaturation = _ORIG_FEATURE


def delta_single(cfg: CONFIG, zmap: str, seed: int, full_width: bool) -> float:
    """Per-side delta on the one-image circular control."""
    rgb = (render_fringe(disp=0.0, seed=seed, width=1200) if full_width
           else _render_circle(seed=seed))

    def run():
        mr = compute_measurement(rgb, cfg, "circle")
        return (float(np.nanmedian(mr.res.diameter_raw)) - CIRCLE_D) / 2.0

    return _with_zmap(zmap, run)


def delta_stack(tmod, cfg: CONFIG, zmap: str, seed0: int) -> float:
    """Per-side delta on the 6-angle stacked circular control (the test's fixture)."""
    def run():
        profiles = {}
        for a in range(1, 7):
            rgb = tmod._render(NOMINAL_ANGLES_DEG[a - 1], tmod.CIRCLE_D,
                               tmod.CIRCLE_D, 0.0, tmod.SHIFTS[a], seed=seed0 + a)
            mr = compute_measurement(rgb, cfg, f"circle_a{a}")
            sp = slice(mr.bnd.x0, mr.bnd.x1 + 1)
            profiles[a] = {
                "x": np.arange(rgb.shape[1], dtype=float)[sp],
                "w": np.where(mr.res.valid[sp], mr.res.diameter_smooth[sp], np.nan),
            }
        st = build_part_stack(1, 1, profiles, cfg)
        return (float(np.nanmedian(st.W)) - tmod.CIRCLE_D) / 2.0

    return _with_zmap(zmap, run)


def summarise(vals: list[float]) -> dict:
    a = np.asarray(vals, dtype=float)
    return {
        "n": int(a.size),
        "mean_px": float(a.mean()),
        "sd_px": float(a.std(ddof=1)) if a.size > 1 else 0.0,
        "min_px": float(a.min()),
        "max_px": float(a.max()),
        "median_px": float(np.median(a)),
        "range_px": float(a.max() - a.min()),
        "values_px": [float(v) for v in a],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=Path("scripts/delta_seeds.json"))
    ap.add_argument("--n-seeds", type=int, default=12,
                    help="number of noise seeds per control (>=10)")
    ap.add_argument("--skip-stack", action="store_true",
                    help="skip the 6-angle control (the slow one)")
    args = ap.parse_args(argv)
    if args.n_seeds < 10:
        raise SystemExit("--n-seeds must be >= 10")

    rev = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                         text=True).stdout.strip()
    branch = subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"],
                            capture_output=True, text=True).stdout.strip()

    single_seeds = [200 + i for i in range(args.n_seeds)]
    stack_seeds = [200 + 100 * i for i in range(args.n_seeds)]

    controls: dict[str, dict] = {}

    for ctrl, full in (("single_image_500x1_as_used_by_calibrate_edge_frac", False),
                       ("single_image_500x1200_full_width", True)):
        controls[ctrl] = {"seeds": single_seeds, "configs": {}}
        for label, (zmap, cfg, desc) in CONFIGS.items():
            vals = [delta_single(cfg, zmap, s, full) for s in single_seeds]
            controls[ctrl]["configs"][label] = {"desc": desc, **summarise(vals)}
            s = controls[ctrl]["configs"][label]
            print(f"{ctrl:52s} {label:8s} mean={s['mean_px']:+7.3f} "
                  f"sd={s['sd_px']:6.3f} min={s['min_px']:+7.3f} max={s['max_px']:+7.3f}")
            sys.stdout.flush()

    if not args.skip_stack:
        tmod = _load_test_module()
        ctrl = "six_angle_stack_as_used_by_test_xsection_synthetic"
        controls[ctrl] = {"seeds": stack_seeds, "configs": {}}
        for label, (zmap, cfg, desc) in CONFIGS.items():
            vals = [delta_stack(tmod, cfg, zmap, s) for s in stack_seeds]
            controls[ctrl]["configs"][label] = {"desc": desc, **summarise(vals)}
            s = controls[ctrl]["configs"][label]
            print(f"{ctrl:52s} {label:8s} mean={s['mean_px']:+7.3f} "
                  f"sd={s['sd_px']:6.3f} min={s['min_px']:+7.3f} max={s['max_px']:+7.3f}")
            sys.stdout.flush()

    out = {
        "kind": "delta_seed_stability",
        "git_branch": branch, "git_rev": rev,
        "true_width_px": CIRCLE_D,
        "n_seeds": args.n_seeds,
        "notes": [
            "delta = (median measured width - 190.0) / 2, per side, in px.",
            "The 'single_image_500x1' control is the geometry the study actually "
            "used: calibrate_edge_frac._render_circle declares W_IMG=1200 but "
            "never uses it, so its images are 500 rows x 1 column.",
            "The 6-angle control is the fixture behind the quoted 3.90 / -0.51.",
            "'old' is the base branch emulated by monkeypatch; the emulation is "
            "verified bit-exact against the base worktree in "
            "scripts/synthetic_fringe.json -> emulation_check.",
        ],
        "controls": controls,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
