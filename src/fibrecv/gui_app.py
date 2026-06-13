"""Local Streamlit GUI for tuning, previewing, batch-processing and exporting.

Dependencies
------------
``streamlit``, ``numpy``, ``pandas``, ``matplotlib`` (Agg), ``imageio`` plus the
reused ``fibrecv`` package: ``compute`` (pure detection core), ``config``,
``io_utils``, ``overlay`` (``render_overlay``), ``register``, ``measure``
(``write_measurement``) and ``run_measure``/``run_aggregate`` (batch + aggregate).
Imports are absolute (``fibrecv.*``) because Streamlit runs this file as a script.

Inputs
------
- Images from a local folder OR any number of uploaded files; both are
  auto-grouped via ``parse_name``'s trailing-numbers rule, with unparseable
  names collected in an "ungrouped" bucket.
- The three boundary knobs (``edge_z``/``edge_frac``/``wcol``) plus the ``ppu``
  calibration, edited in a sidebar form and applied on demand; all other
  ``CONFIG`` fields stay at the validated defaults.
- An output-folder path for export/batch.

Output
------
- Live, in-memory preview: full-res boundary overlays, per-replicate diameter
  profiles, and a registered mean+/-std group curve -- all redrawn when the user
  changes parameters and clicks Apply (no disk writes for preview).
- Manual boundary correction: per replicate, the user can click anchor points
  on a zoomed strip (or nudge a whole line) to redraw the detected top/bottom
  boundary where detection fails. Points are grouped into independent sets,
  one per corrected stretch, so two far-apart fixes are never joined across
  the gap. Corrections (``manual_edit``) re-run QC and flow into the profile
  plot, group registration and export (drawn in magenta). Batch runs recompute
  from disk and never see these edits.
- On request: the standard fibrecv output tree (overlays/, per_image/*,
  per_sample/*, summary/master_summary.csv, run_config.json) written locally for
  the current group or for a whole folder (in-process batch with a progress bar).

Pos
---
Thin front-end over the validated pipeline; adds no detection logic. Designed to
run on the user's local Mac/Windows machine (``fibrecv-gui`` or
``streamlit run src/fibrecv/gui_app.py``), reading a copied images folder. The
heavy compute reuses ``compute.compute_measurement`` so preview == CLI output.
"""

from __future__ import annotations

import io
import json
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import imageio.v3 as iio  # noqa: E402
import streamlit as st  # noqa: E402
import streamlit.components.v1 as components  # noqa: E402

from fibrecv import run_aggregate  # noqa: E402
from fibrecv.compute import compute_measurement  # noqa: E402
from fibrecv.config import CONFIG  # noqa: E402
from fibrecv.io_utils import (  # noqa: E402
    IMAGE_SUFFIXES, discover_images, natural_key, parse_name)
from fibrecv.io_utils import load_rgb as _io_load_rgb  # noqa: E402
from fibrecv.manual_edit import (  # noqa: E402
    apply_manual_edits, display_to_native, empty_edits, has_edits)
from fibrecv.measure import write_measurement  # noqa: E402
from fibrecv.overlay import GREY, WHITE, mark_anchors, render_overlay  # noqa: E402
from fibrecv.register import register_sample  # noqa: E402
from fibrecv.run_measure import _lib_versions, _worker  # noqa: E402
from fibrecv.tensile import (  # noqa: E402
    build_matrix, compute_tensile, discover_tensile, discover_tensile_files,
    read_trace)
from streamlit_image_coordinates import streamlit_image_coordinates  # noqa: E402

DEFAULTS = CONFIG()  # never mutated; the source of widget defaults + reset target

# --- visible parameters: the three knobs that move the detected boundary,
# plus the ppu calibration. Everything else in CONFIG stays at the validated
# defaults (CLI keeps full control). spec = (name, kind, help, step, lo, hi,
# fmt). ---
PARAM_SPECS: list[tuple] = [
    ("edge_z", "slider",
     "Where on the fibre wall the boundary line is drawn (the main knob). The "
     "detector finds the steep brightness wall at each fibre edge and draws "
     "the line where the signal has risen edge_z units above the background "
     "just outside that wall. Higher = the line sits higher up the wall, "
     "further inside the fibre = thinner reading; lower = further out = "
     "thicker reading. If the line cuts into the fibre (too thin), lower it; "
     "if it sits in the shadow outside the fibre (too thick), raise it. "
     "Recommended 3-5, default 4.0.",
     0.5, 1.0, 12.0, "%.1f"),
    ("edge_frac", "float",
     "Protection for faint fibres. A pale, low-contrast fibre has a weak wall "
     "that may never rise edge_z units above background; this caps the "
     "crossing level at edge_frac of the wall's own height so the line stays "
     "on the wall instead of drifting outward. It only kicks in when the wall "
     "is weaker than edge_z. Leave at 0.65 unless a very faint fibre loses "
     "its boundary — then lower it slightly.",
     0.05, 0.0, 1.0, "%.2f"),
    ("wcol", "int",
     "Horizontal smoothing width (pixels). If the boundary line is jittery or "
     "ragged, raise it (61-81) for a smoother, more stable line; lower it "
     "(15-25) to preserve fine thickness variation. Too high flattens real "
     "variation. Default 41.",
     1, 1, 201, None),
    ("ppu", "float",
     "Calibration: camera pixels per micron. Diameters in µm = pixels / ppu, "
     "so this scales every µm number in the app and the exports (pixel values "
     "are unaffected). Measure it with a stage micrometer for your microscope "
     "+ camera combination. Default 1.3680 (the original calibrated setup).",
     0.001, 0.1, 10.0, "%.4f"),
]

# names of int-typed visible fields (so widgets return int, not float)
_INT_FIELDS = {name for (name, kind, *_rest) in PARAM_SPECS if kind == "int"}


# --------------------------------------------------------------------------- #
# Config <-> cache-key helpers                                                 #
# --------------------------------------------------------------------------- #
def _cfg_items(cfg_dict: dict) -> tuple:
    """Hashable, order-stable view of a config dict for cache keys."""
    return tuple(sorted(cfg_dict.items()))


def _cfg_from_items(cfg_items: tuple) -> CONFIG:
    """Rebuild a CONFIG from a cache-key tuple (coercing int fields)."""
    d = dict(cfg_items)
    for k in _INT_FIELDS:
        if k in d:
            d[k] = int(round(d[k]))
    return replace(CONFIG(), **d)


def _rgb_from_bytes(data: bytes) -> np.ndarray:
    """Decode uploaded image bytes to float32 RGB in [0, 1] (mirrors load_rgb)."""
    arr = np.asarray(iio.imread(io.BytesIO(data)))
    if arr.ndim == 2:
        arr = np.repeat(arr[:, :, None], 3, axis=2)
    if arr.shape[2] == 4:
        arr = arr[:, :, :3]
    if arr.dtype != np.float32 and arr.dtype != np.float64:
        arr = arr.astype(np.float32) / 255.0
    else:
        arr = arr.astype(np.float32)
        if arr.max() > 1.0:
            arr = arr / 255.0
    return arr


# --------------------------------------------------------------------------- #
# Cached loaders + compute (keyed so unchanged inputs return instantly)        #
# --------------------------------------------------------------------------- #
@st.cache_data(show_spinner=False, max_entries=12)
def _cached_rgb_from_path(path: str, mtime: float) -> np.ndarray:
    """Decode a JPEG once; ``mtime`` invalidates the cache if the file changes."""
    return _io_load_rgb(path)


@st.cache_data(show_spinner=False, max_entries=24)
def _cached_compute_path(path: str, mtime: float, cfg_items: tuple):
    """Compute one folder image; cached on (path, mtime, params). rgb/D dropped."""
    cfg = _cfg_from_items(cfg_items)
    rgb = _cached_rgb_from_path(path, mtime)
    mr = compute_measurement(rgb, cfg, name=Path(path).stem)
    mr.rgb = None   # the big arrays are not needed downstream; keep the cache light
    mr.D = None
    return mr


@st.cache_data(show_spinner=False, max_entries=8)
def _cached_discover_tensile(folder: str) -> dict[str, str]:
    """``discover_tensile`` keyed by folder string (Path -> str so it caches)."""
    return {g: str(p) for g, p in discover_tensile(folder).items()}


@st.cache_data(show_spinner=False, max_entries=12)
def _cached_compute_upload(file_key: str, data: bytes, cfg_items: tuple):
    """Compute one uploaded image; returns (mr, rgb) since uploads keep no path."""
    cfg = _cfg_from_items(cfg_items)
    rgb = _rgb_from_bytes(data)
    mr = compute_measurement(rgb, cfg, name=file_key)
    mr.D = None
    return mr, rgb


UNGROUPED = "ungrouped"


def _group_sort_key(group: str):
    """Natural sort for group labels ('3_1' < '3_3' < '10_5'); ungrouped last."""
    return (group == UNGROUPED, natural_key(group))


def _group_by_name(items: list, key) -> dict[str, list]:
    """Bucket items by their parse_name group; unparseable names -> UNGROUPED."""
    groups: dict[str, list] = {}
    for it in items:
        try:
            g, _ = parse_name(key(it))
        except ValueError:
            g = UNGROUPED
        groups.setdefault(g, []).append(it)
    return groups


def _sorted_reps(items: list, key) -> list:
    """Sort one group's items by replicate number; unparseable last, by name."""
    def _k(it):
        try:
            return (0, parse_name(key(it))[1], key(it))
        except ValueError:
            return (1, 0, key(it))
    return sorted(items, key=_k)


# --------------------------------------------------------------------------- #
# Pure export / batch logic (no Streamlit calls -> headlessly testable)        #
# --------------------------------------------------------------------------- #
def export_group(reps: list[dict], out_root: str | Path, cfg: CONFIG) -> str:
    """Write the standard output tree for one already-computed group.

    ``reps`` is a list of ``{"name", "rgb", "mr"}`` dicts. Writes per-image
    artifacts (via the shared ``write_measurement``) then reuses
    ``run_aggregate.main`` to build per_sample/* and the master_summary row for
    this group, so the export is identical to the two-stage CLI. Returns the
    group label. Raises ``ValueError`` if the group is not name-parseable.
    """
    out_root = Path(out_root)
    group = next((r["mr"].group for r in reps if r["mr"].group is not None), None)
    if group is None:
        raise ValueError(
            "Cannot export: image names must end in numbers "
            "(e.g. 'name 3_1_2.jpg') to derive a group label."
        )
    for rep in reps:
        write_measurement(rep["rgb"], rep["mr"], cfg, out_root)
    run_aggregate.main([
        "--out", str(out_root), "--groups", group,
        "--ppu", str(cfg.ppu), "--max-shift", str(cfg.max_shift),
        "--min-corr", str(cfg.min_corr), "--min-coverage", str(cfg.min_coverage),
    ])
    return group


def export_all_groups(grouped_reps: dict[str, list[dict]], out_root: str | Path,
                      cfg: CONFIG) -> list[str]:
    """Write the output tree for *every* group in one session, not just the loaded one.

    ``grouped_reps`` maps group label -> its already-computed reps. Writes each
    image's artifacts via ``write_measurement`` (so manual edits carried on the
    loaded group's reps are honoured), then runs a single ``run_aggregate`` pass
    over all of them so ``summary/master_summary.csv`` has one row per group.
    Returns the list of exported group labels.
    """
    out_root = Path(out_root)
    groups: list[str] = []
    for reps in grouped_reps.values():
        group = next((r["mr"].group for r in reps if r["mr"].group is not None), None)
        if group is None:
            continue
        for rep in reps:
            write_measurement(rep["rgb"], rep["mr"], cfg, out_root)
        groups.append(group)
    if groups:
        run_aggregate.main([
            "--out", str(out_root), "--groups", *groups,
            "--ppu", str(cfg.ppu), "--max-shift", str(cfg.max_shift),
            "--min-corr", str(cfg.min_corr), "--min-coverage", str(cfg.min_coverage),
        ])
    return groups


def run_batch(
    image_paths: list[Path],
    out_root: str | Path,
    cfg: CONFIG,
    jobs: int,
    progress_cb: Callable[[float, dict], None] | None = None,
) -> tuple[pd.DataFrame, list[dict]]:
    """Measure every image in-process, then aggregate the whole folder.

    Tries a ``ProcessPoolExecutor`` (reusing the picklable ``run_measure._worker``)
    and falls back to sequential if the pool cannot start (Windows ``spawn``
    safety); both paths report progress through ``progress_cb(frac, result)``.
    Writes the full output tree (per_image/*, overlays/, per_sample/*,
    summary/master_summary.csv, run_config.json) and returns
    ``(master_summary_df, per_image_results)``.
    """
    out_root = Path(out_root)
    (out_root / "summary").mkdir(parents=True, exist_ok=True)
    n = len(image_paths)

    # provenance snapshot (mirrors run_measure.main)
    with open(out_root / "summary" / "run_config.json", "w") as fh:
        json.dump({"params": cfg.as_dict(), "versions": _lib_versions(),
                   "n_images": n, "root": str(image_paths[0].parent) if n else ""},
                  fh, indent=2)

    results: list[dict] | None = None
    if jobs and jobs > 1:
        try:
            collected: list[dict] = []
            with ProcessPoolExecutor(max_workers=int(jobs)) as ex:
                futs = {ex.submit(_worker, str(p), cfg, str(out_root)): p
                        for p in image_paths}
                for i, fut in enumerate(as_completed(futs)):
                    r = fut.result()
                    collected.append(r)
                    if progress_cb:
                        progress_cb((i + 1) / n, r)
            results = collected
        except Exception:  # noqa: BLE001 - any pool start failure -> sequential
            results = None

    if results is None:
        results = []
        for i, p in enumerate(image_paths):
            r = _worker(str(p), cfg, str(out_root))
            results.append(r)
            if progress_cb:
                progress_cb((i + 1) / n, r)

    # run log (mirrors run_measure.main)
    errors = [r for r in results if "error" in r]
    low = [r for r in results if r.get("low_confidence")]
    with open(out_root / "summary" / "run_log.txt", "w") as fh:
        fh.write(f"images={n} ok={len(results) - len(errors)} "
                 f"errors={len(errors)} low_confidence={len(low)}\n\n")
        for r in sorted(results, key=lambda d: d.get("name", "")):
            if "error" in r:
                fh.write(f"[ERROR] {r['name']}: {r['error']}\n")
            else:
                fh.write(f"{r['name']}: cov={r['coverage']:.0%} "
                         f"med={r['median_diameter_um']} tilt={r['tilt_slope']:.4f} "
                         f"{'LOWCONF' if r['low_confidence'] else ''}\n")

    # aggregate the whole folder (reuses the validated CLI aggregator)
    run_aggregate.main([
        "--out", str(out_root), "--all",
        "--ppu", str(cfg.ppu), "--max-shift", str(cfg.max_shift),
        "--min-corr", str(cfg.min_corr), "--min-coverage", str(cfg.min_coverage),
    ])
    master_path = out_root / "summary" / "master_summary.csv"
    master = pd.read_csv(master_path) if master_path.exists() else pd.DataFrame()
    return master, results


def _group_mean_um(reps: list[dict], cfg: CONFIG) -> float | None:
    """Registered mean diameter (µm) for one group's reps, or ``None`` if it
    cannot be measured. Same profile-building + registration as ``_render_group``,
    so the matrix matches the on-screen group panel."""
    profiles = []
    for rep in reps:
        mr = rep["mr"]
        res, bnd = mr.res, mr.bnd
        if res.band_mismatch or res.coverage < cfg.min_coverage:
            continue
        W = rep["rgb"].shape[1]
        span = slice(bnd.x0, bnd.x1 + 1)
        profiles.append({
            "replicate": mr.replicate if mr.replicate is not None else rep["idx"] + 1,
            "coverage": res.coverage,
            "x": np.arange(W)[span].astype(float),
            "diameter_px_raw": res.diameter_raw[span].astype(float),
            "diameter_px_smooth": res.diameter_smooth[span].astype(float),
            "valid": res.valid[span].astype(bool),
        })
    if not profiles:
        return None
    profiles.sort(key=lambda p: p["replicate"])
    try:
        _, _, summary = register_sample(profiles, cfg)
    except Exception:  # noqa: BLE001
        return None
    mean_um = summary.get("mean_um")
    return float(mean_um) if mean_um is not None and np.isfinite(mean_um) else None


def _diameters_from_uploads(image_uploads, cfg_items: tuple) -> dict[str, float]:
    """Measure every uploaded image in-memory -> ``{group: mean diameter µm}``.

    The folder path uses ``run_batch``; this is its upload twin so the tensile
    matrix works with drag-and-drop images too. Reuses ``_cached_compute_upload``
    (so the loaded group is not recomputed) and the same registration as the
    group panel. ``cfg_items`` is the image config (not the tensile one) to keep
    the cache keys identical to the live view.
    """
    cfg = _cfg_from_items(cfg_items)
    diameters: dict[str, float] = {}
    groups = _group_by_name(list(image_uploads), key=lambda u: u.name)
    for g, items in groups.items():
        if g == UNGROUPED:
            continue
        reps = []
        for i, up in enumerate(_sorted_reps(items, key=lambda u: u.name)):
            mr, rgb = _cached_compute_upload(Path(up.name).stem, up.getvalue(),
                                             cfg_items)
            reps.append({"name": Path(up.name).stem, "rgb": rgb, "mr": mr, "idx": i})
        mean_um = _group_mean_um(reps, cfg)
        if mean_um is not None:
            diameters[g] = mean_um
    return diameters


def _grouped_reps_from_uploads(image_uploads, cfg_items: tuple,
                               loaded_reps: list[dict] | None,
                               group_label: str | None) -> dict[str, list[dict]]:
    """Build ``{group: reps}`` for *every* uploaded group (for an all-groups export).

    Reuses ``_cached_compute_upload`` (so nothing is recomputed needlessly) and
    keeps the on-screen group's already-loaded reps so its manual edits are
    preserved; other groups are measured fresh.
    """
    grouped = _group_by_name(list(image_uploads), key=lambda u: u.name)
    out: dict[str, list[dict]] = {}
    for g, items in grouped.items():
        if g == UNGROUPED:
            continue
        if g == group_label and loaded_reps:
            out[g] = loaded_reps
            continue
        reps = []
        for i, up in enumerate(_sorted_reps(items, key=lambda u: u.name)):
            mr, rgb = _cached_compute_upload(Path(up.name).stem, up.getvalue(),
                                             cfg_items)
            reps.append({"name": Path(up.name).stem, "rgb": rgb, "mr": mr, "idx": i})
        out[g] = reps
    return out


# --------------------------------------------------------------------------- #
# Plot builders                                                               #
# --------------------------------------------------------------------------- #
def _profile_fig(mr, rgb, cfg: CONFIG):
    """Per-replicate diameter-vs-position figure (raw points + smoothed line, µm)."""
    bnd, res = mr.bnd, mr.res
    span = slice(bnd.x0, bnd.x1 + 1)
    x = np.arange(rgb.shape[1])[span]
    raw_um = mr.diameter_um[span]
    sm_um = res.diameter_smooth[span] / cfg.ppu
    fig, ax = plt.subplots(figsize=(9, 3))
    ax.plot(x, raw_um, ".", ms=2, alpha=0.4, color="tab:gray", label="raw")
    ax.plot(x, sm_um, "-", lw=1.3, color="tab:blue", label="smooth")
    ax.set_xlabel("x position (px)")
    ax.set_ylabel("diameter (µm)")
    ax.set_title(f"{mr.name}   coverage={res.coverage:.0%}")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def _group_fig(table: dict, group_label: str, rep_curves: list[tuple]):
    """Aligned replicate curves behind the registered mean +/- std band.

    ``rep_curves`` is a list of ``(replicate, x_aligned, diameter_um)`` tuples,
    one per replicate, drawn thin/faded so the user can see what the mean and
    std are built from.
    """
    fig, ax = plt.subplots(figsize=(9, 3))
    cmap = plt.get_cmap("tab10")
    for i, (rep, rx, ry) in enumerate(rep_curves):
        ax.plot(rx, ry, "-", lw=0.8, alpha=0.45, color=cmap(i % 10),
                label=f"rep {rep}")
    x = table["x_aligned_px"]
    mean = table["mean_um"]
    std = table["std_um"]
    ax.plot(x, mean, "-", lw=2.0, color="tab:blue", zorder=5,
            label="mean of replicates")
    band = np.where(np.isfinite(std), std, 0.0)
    ax.fill_between(x, mean - band, mean + band, alpha=0.25, color="tab:blue",
                    label="±std across replicates")
    ax.set_xlabel("aligned x position (px)")
    ax.set_ylabel("diameter (µm)")
    ax.set_title(f"sample {group_label} — replicates aligned and averaged "
                 f"(n_reps={int(np.nanmax(table['n']))})")
    ax.legend(loc="best", fontsize=7, ncols=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def _tensile_fig(res):
    """Stress-strain curve with the modulus fit, toughness area and break point.

    Mirrors ``_group_fig``'s style (Agg, 9x3, grid, tight_layout). Plots stress
    (MPa) vs strain (%) up to fracture as the main line, a faded post-break
    tail, the shaded toughness area, the steepest-slope modulus fit line over
    its own segment and the marked break point. When no diameter matched the
    fibre (all-NaN stress) it falls back to the raw force/displacement trace so
    the curve is still informative.
    """
    fig, ax = plt.subplots(figsize=(9, 3))
    group = res.group or "fibre"
    brk = int(res.break_index)

    def _focus(x_used, y_used, has_tail):
        """Bound the view to the pre-break curve so a post-break recoil spike in
        the faded tail cannot dominate the axes (a short stub still shows)."""
        xu = np.asarray(x_used, float); yu = np.asarray(y_used, float)
        xu = xu[np.isfinite(xu)]; yu = yu[np.isfinite(yu)]
        if not xu.size or not yu.size:
            return
        ymin, ymax = min(0.0, float(yu.min())), float(yu.max())
        pad = 0.10 * ((ymax - ymin) or 1.0)
        ax.set_ylim(ymin - 0.3 * pad, ymax + pad)
        xlo, xhi = min(0.0, float(xu.min())), float(xu.max())
        if has_tail:
            xhi += 0.30 * ((xhi - xlo) or 1.0)
        ax.set_xlim(xlo, xhi)

    # no matched diameter -> stress is all NaN; show force vs displacement instead
    if not np.isfinite(res.stress_pa).any():
        x = res.disp_mm
        y = res.load_n
        ax.plot(x[:brk + 1], y[:brk + 1], "-", lw=1.6, color="tab:blue",
                label="load")
        if brk + 1 < x.size:
            ax.plot(x[brk:], y[brk:], "-", lw=0.8, alpha=0.35, color="tab:blue")
        ax.scatter([x[brk]], [y[brk]], s=30, color="tab:red", zorder=6,
                   label=f"break (Fmax = {res.fmax_n:.3f} N)")
        _focus(x[:brk + 1], y[:brk + 1], brk + 1 < x.size)
        ax.set_xlabel("displacement (mm)")
        ax.set_ylabel("load (N)")
        ax.set_title(f"{group} — stress–strain "
                     f"(no matched diameter — force/displacement only)")
        ax.legend(loc="best", fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        return fig

    x = res.strain * 100.0          # %
    y = res.stress_pa / 1e6         # MPa
    ax.plot(x[:brk + 1], y[:brk + 1], "-", lw=1.6, color="tab:blue",
            label="stress–strain")
    if brk + 1 < x.size:
        ax.plot(x[brk:], y[brk:], "-", lw=0.8, alpha=0.35, color="tab:blue",
                label="post-break (unused)")

    # toughness = area under the curve up to fracture
    ax.fill_between(x[:brk + 1], 0.0, y[:brk + 1], alpha=0.15,
                    color="tab:blue", label="toughness (area)")

    # steepest initial slope = Young's modulus, drawn over its fit segment
    fit = res.modulus_fit or {}
    if fit and np.isfinite(fit.get("slope", np.nan)):
        s_lo, s_hi = fit["strain_lo"], fit["strain_hi"]
        s_line = np.array([s_lo, s_hi], dtype=float)
        stress_fit = fit["slope"] * s_line + fit["intercept"]
        ax.plot(s_line * 100.0, stress_fit / 1e6, "--", lw=1.4,
                color="tab:orange",
                label=f"E = {res.youngs_modulus_pa / 1e9:.2f} GPa")

    # break point
    ax.scatter([x[brk]], [y[brk]], s=30, color="tab:red", zorder=6,
               label=f"break (Fmax = {res.fmax_n:.3f} N)")

    _focus(x[:brk + 1], y[:brk + 1], brk + 1 < x.size)
    ax.set_xlabel("strain (%)")
    ax.set_ylabel("stress (MPa)")
    ax.set_title(f"{group} — stress–strain")
    ax.legend(loc="best", fontsize=7)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------- #
# Sidebar: parameter form                                                     #
# --------------------------------------------------------------------------- #
def _param_form() -> None:
    """Render the 3-knob parameter form; updates session_state on Apply/Reset."""
    applied = st.session_state.cfg_dict
    ver = st.session_state.form_version

    with st.sidebar.form("params", clear_on_submit=False):
        st.markdown("**Parameters** — edit, then click **Apply** to re-render.")
        new_vals: dict = {}
        for (name, kind, help_txt, step, lo, hi, fmt) in PARAM_SPECS:
            if name == "ppu":
                st.markdown("**Calibration**")
            key = f"p_{name}_v{ver}"
            cur = applied[name]
            if kind == "slider":
                new_vals[name] = st.slider(
                    name, min_value=float(lo), max_value=float(hi),
                    value=float(cur), step=float(step), help=help_txt, key=key)
            elif kind == "int":
                new_vals[name] = int(st.number_input(
                    name, min_value=int(lo), max_value=int(hi),
                    value=int(cur), step=int(step), help=help_txt, key=key))
            else:  # float
                kwargs = dict(min_value=float(lo), max_value=float(hi),
                              value=float(cur), step=float(step), help=help_txt, key=key)
                if fmt:
                    kwargs["format"] = fmt
                new_vals[name] = float(st.number_input(name, **kwargs))
        c1, c2 = st.columns(2)
        apply = c1.form_submit_button("Apply", type="primary", width="stretch")
        reset = c2.form_submit_button("Reset to defaults", width="stretch")

    if reset:
        st.session_state.cfg_dict = DEFAULTS.as_dict()
        st.session_state.form_version += 1
        st.rerun()
    if apply:
        # merge: the three knobs override a full defaults dict, so hidden
        # fields always carry the validated values
        st.session_state.cfg_dict = {**DEFAULTS.as_dict(), **new_vals}
        st.rerun()  # re-run top-to-bottom so reps recompute with the new params


# --------------------------------------------------------------------------- #
# Sidebar: tensile controls                                                   #
# --------------------------------------------------------------------------- #
def _tensile_controls() -> dict:
    """Render the sidebar tensile section; return the fibre→source map + params.

    Like the image data source, tensile data can come from a local folder or
    drag-and-drop uploads. Either way the result is ``tmap`` = ``{group: source}``
    (a ``Path`` for the folder, an uploaded file for uploads), both of which feed
    the polymorphic ``read_trace``. The tester logs only crosshead displacement
    and force, so the gauge length sets the strain scale and the modulus window
    the auto-fit width; these feed a tensile-specific config in ``main()`` (the
    diameter knobs are untouched). Returns
    ``{"tmap", "folder", "gauge_length_mm", "modulus_window"}``.
    """
    st.sidebar.markdown("### Tensile data")
    source = st.sidebar.radio("tensile source", ["Local folder", "Upload"],
                              horizontal=True, label_visibility="collapsed",
                              key="tensile_source")
    tmap: dict[str, object] = {}
    folder: str | None = None
    if source == "Local folder":
        folder = st.sidebar.text_input(
            "Tensile data folder", value=st.session_state.get("tensile_folder", ""))
        st.session_state.tensile_folder = folder
        if folder and Path(folder).is_dir():
            try:
                tmap = {g: Path(p)
                        for g, p in _cached_discover_tensile(folder).items()}
            except Exception as exc:  # noqa: BLE001
                st.sidebar.warning(f"Could not scan tensile folder: {exc}")
        elif folder:
            st.sidebar.warning("Enter a valid tensile folder path.")
    else:
        folder_mode = st.sidebar.checkbox(
            "📁 upload a whole folder", key="tensile_folder_mode",
            help="Make Browse open a folder chooser; non-tensile files are ignored.")
        ups = st.sidebar.file_uploader(
            "Upload tensile files",
            type=None if folder_mode else ["csv", "xls", "xlsx"],
            accept_multiple_files=True, key="tensile_uploads")
        if folder_mode:
            _enable_folder_upload("Upload tensile files")
        if ups:
            tmap = discover_tensile_files(ups)
    if tmap:
        st.sidebar.caption(f"{len(tmap)} tensile fibre(s) matched.")

    gauge_length_mm = float(st.sidebar.number_input(
        "Gauge length L₀ (mm)", min_value=0.1, value=float(DEFAULTS.gauge_length_mm),
        step=1.0,
        help="Grip separation; strain = displacement / L₀. The tester records "
             "only displacement, so this sets the strain scale."))
    modulus_window = float(st.sidebar.number_input(
        "Modulus window (fraction)", min_value=0.02, max_value=0.5,
        value=float(DEFAULTS.modulus_window), step=0.01,
        help="Width of the sliding linear fit used to auto-detect the steepest "
             "initial slope (Young's modulus), as a fraction of the rising "
             "region."))
    return {
        "tmap": tmap,
        "folder": folder or None,
        "gauge_length_mm": gauge_length_mm,
        "modulus_window": modulus_window,
    }


# --------------------------------------------------------------------------- #
# Sidebar: data source -> list of replicate dicts                             #
# --------------------------------------------------------------------------- #
def _enable_folder_upload(label: str) -> None:
    """Let the ``st.file_uploader`` whose label contains ``label`` take a folder.

    Sets ``webkitdirectory`` on that uploader's hidden ``<input type=file>`` from a
    zero-height html component (same-origin, reaching the parent document). The
    Browse button then opens a folder chooser and the browser hands back every
    file inside; the caller keeps only the ones with the right extension.
    """
    js = """
    <script>
    const WANT = %s;
    const doc = window.parent.document;
    function apply() {
      doc.querySelectorAll('[data-testid="stFileUploader"]').forEach(u => {
        if ((u.innerText || "").includes(WANT)) {
          const inp = u.querySelector('input[type="file"]');
          if (inp) {
            inp.setAttribute("webkitdirectory", "");
            inp.setAttribute("directory", "");
          }
        }
      });
    }
    apply(); setTimeout(apply, 150); setTimeout(apply, 500);
    </script>
    """ % json.dumps(label)
    with st.sidebar:
        components.html(js, height=0, width=0)


def _load_reps(cfg_items: tuple) -> tuple[list[dict], str | None, str | None]:
    """Resolve the data-source controls to a list of replicate dicts.

    Returns ``(reps, group_label, folder)``. Each rep is
    ``{"name", "rgb", "mr", "idx"}``. ``folder`` is the scanned folder path
    (for batch) or None for the upload source. Both sources group images with
    the shared ``parse_name`` rule; unparseable names land in an "ungrouped"
    bucket instead of being hidden.
    """
    st.sidebar.markdown("### Data source")
    source = st.sidebar.radio("source", ["Local folder", "Upload"],
                              horizontal=True, label_visibility="collapsed")
    reps: list[dict] = []
    group_label: str | None = None
    folder: str | None = None
    st.session_state["image_uploads"] = []  # all uploaded images (for matrix export)

    if source == "Local folder":
        folder = st.sidebar.text_input("Image folder", value=st.session_state.get(
            "folder", ""))
        st.session_state.folder = folder
        if not folder or not Path(folder).is_dir():
            st.sidebar.warning("Enter a valid local folder path.")
            return reps, None, None
        paths = discover_images(folder)
        if not paths:
            st.sidebar.warning("No image files found in this folder.")
            return reps, None, folder
        groups = _group_by_name(paths, key=lambda p: p.name)
        keys = sorted(groups, key=_group_sort_key)
        group_label = st.sidebar.selectbox("Group", keys)
        st.sidebar.caption(f"{len(paths)} images, {len(keys)} groups")
        for i, p in enumerate(_sorted_reps(groups[group_label],
                                           key=lambda p: p.name)):
            mtime = Path(p).stat().st_mtime
            mr = _cached_compute_path(str(p), mtime, cfg_items)
            rgb = _cached_rgb_from_path(str(p), mtime)
            reps.append({"name": Path(p).stem, "rgb": rgb, "mr": mr, "idx": i})
    else:
        folder_mode = st.sidebar.checkbox(
            "📁 upload a whole folder", key="img_folder_mode",
            help="Make Browse open a folder chooser and ingest every image inside. "
                 "You can also drag a folder onto the box either way.")
        uploads = st.sidebar.file_uploader(
            "Upload images",
            type=None if folder_mode else ["jpg", "jpeg", "png", "tif", "tiff", "bmp"],
            accept_multiple_files=True)
        if folder_mode:
            _enable_folder_upload("Upload images")
            uploads = [u for u in (uploads or [])
                       if Path(u.name).suffix.lower() in IMAGE_SUFFIXES]
        if uploads:
            st.session_state["image_uploads"] = list(uploads)
            groups = _group_by_name(list(uploads), key=lambda u: u.name)
            keys = sorted(groups, key=_group_sort_key)
            group_label = (st.sidebar.selectbox("Group", keys)
                           if len(keys) > 1 else keys[0])
            st.sidebar.caption(f"{len(uploads)} files, {len(keys)} groups")
            for i, up in enumerate(_sorted_reps(groups[group_label],
                                                key=lambda u: u.name)):
                data = up.getvalue()
                stem = Path(up.name).stem
                mr, rgb = _cached_compute_upload(stem, data, cfg_items)
                reps.append({"name": stem, "rgb": rgb, "mr": mr, "idx": i})
            if group_label == UNGROUPED:
                group_label = None  # header falls back to "Replicates (uploaded)"

    return reps, group_label, folder


# --------------------------------------------------------------------------- #
# Main-area renderers                                                         #
# --------------------------------------------------------------------------- #
def _render_replicate(rep: dict, cfg: CONFIG) -> None:
    mr, rgb = rep["mr"], rep["rgb"]
    bnd, edg, res = mr.bnd, mr.edg, mr.res
    edited_top, edited_bot = rep.get("edited_top"), rep.get("edited_bot")
    overlay = render_overlay(rgb, edg.y_top, edg.y_bot, bnd.c_fit, res.valid,
                             bnd.x0, bnd.x1, thick=1,
                             edited_top=edited_top, edited_bot=edited_bot)
    caption = f"{mr.name} — cyan top / yellow bottom / dashed centerline"
    if edited_top is not None or edited_bot is not None:
        caption += " / magenta = manual edit"
    st.image(overlay, width="stretch", caption=caption)

    med = mr.meta["median_diameter_um"]
    d = mr.diameter_um
    have_d = bool(np.isfinite(d).any())
    mean_um = float(np.nanmean(d)) if have_d else None
    std_um = float(np.nanstd(d)) if have_d else None
    flags = []
    if res.low_confidence:
        flags.append("low_confidence")
    if res.band_mismatch:
        flags.append("band_mismatch")
    c = st.columns(6)
    c[0].metric("mean Ø", f"{mean_um:.2f} µm" if mean_um is not None else "—",
                help="Mean diameter of this image, averaged along the fibre "
                     "(valid columns only).")
    c[1].metric("median Ø", f"{med:.2f} µm" if med is not None else "—")
    c[2].metric("along-fibre std", f"{std_um:.2f} µm" if std_um is not None else "—",
                help="Std of the diameter along this image's fibre — thickness "
                     "variation within the picture, not disagreement between "
                     "replicates (that is the group panel's std).")
    c[3].metric("coverage", f"{res.coverage:.0%}")
    c[4].metric("tilt slope", f"{bnd.slope:.4f}")
    c[5].metric("flags", ", ".join(flags) if flags else "none")

    _render_edit_expander(rep)

    fig = _profile_fig(mr, rgb, cfg)
    st.pyplot(fig)
    plt.close(fig)


def _render_edit_expander(rep: dict) -> None:
    """Manual boundary correction: click anchors on a zoomed strip, or nudge.

    Anchors and nudges live in ``st.session_state.manual_edits[name]`` in
    native full-image pixel coordinates, so they survive parameter changes;
    ``main()`` applies them via ``apply_manual_edits`` right after loading, so
    by the time this renders, ``rep["mr"]`` already reflects the edits.
    """
    mr, rgb, name = rep["mr"], rep["rgb"], rep["name"]
    bnd, edg, res = mr.bnd, mr.edg, mr.res
    H, W = rgb.shape[:2]
    edits = st.session_state.manual_edits.get(name) or empty_edits()

    with st.expander("Edit boundaries (manual correction)"):
        side = st.radio("Line to edit", ["top (cyan)", "bottom (yellow)"],
                        horizontal=True, key=f"side_{name}")
        side_key = "top" if side.startswith("top") else "bot"
        st.caption(
            "Click 2+ points along the true edge in the zoomed strip below; "
            "the active set's line is redrawn through its points and blended "
            "into the detected line at the ends (a single click edits one "
            "column). Each set corrects one stretch independently — start a "
            "new set to fix another stretch without the two being joined. "
            "Corrections are drawn in magenta and feed the profile plot, the "
            "group statistics and the export. Points outside the detected "
            "band span stay invalid.")

        # anchor sets of the side being edited; the radio's key carries a
        # version counter (same pattern as the param form) so 'Start new set'
        # and deletions can re-default it to the newest set
        sets_cur = edits[side_key]
        ver_key = f"setver_{name}"
        ver = st.session_state.setdefault(ver_key, 0)
        ctl = st.columns([3, 1])
        if sets_cur:
            active = ctl[0].radio(
                "Set to extend (new clicks are added to it)",
                options=list(range(len(sets_cur))),
                format_func=lambda i: f"set {i + 1} ({len(sets_cur[i])} pts)",
                horizontal=True, index=len(sets_cur) - 1,
                key=f"actset_{side_key}_{name}_v{ver}")
        else:
            active = None
            ctl[0].caption("No sets yet — the first click starts set 1.")
        if ctl[1].button("Start new set", key=f"newset_{name}",
                         disabled=bool(sets_cur) and not sets_cur[-1]):
            ed = st.session_state.manual_edits.setdefault(name, empty_edits())
            ed[side_key].append([])
            st.session_state[ver_key] = ver + 1
            st.rerun()

        # zoomed strip: full width, vertically cropped around the band so
        # clicks have much better y-resolution than on the full image
        span = slice(bnd.x0, bnd.x1 + 1)
        c_span = bnd.c_fit[span]
        pad = 12
        y_lo = max(0, int(np.floor(np.nanmin(c_span) - edg.half_window - pad)))
        y_hi = min(H, int(np.ceil(np.nanmax(c_span) + edg.half_window + pad)) + 1)
        overlay_full = render_overlay(
            rgb, edg.y_top, edg.y_bot, bnd.c_fit, res.valid, bnd.x0, bnd.x1,
            thick=1, edited_top=rep.get("edited_top"),
            edited_bot=rep.get("edited_bot"))
        for i, anchor_set in enumerate(sets_cur):
            mark_anchors(overlay_full, anchor_set,
                         color=WHITE if i == active else GREY)
        crop = overlay_full[y_lo:y_hi]

        value = streamlit_image_coordinates(crop, width=1200, height=360,
                                            key=f"clk_{name}")
        # the component re-emits the last click on every rerun -> dedupe on
        # its unix_time stamp, otherwise anchors would silently duplicate
        if value and value.get("unix_time") != st.session_state.get(f"last_click_{name}"):
            st.session_state[f"last_click_{name}"] = value["unix_time"]
            x_nat, y_nat = display_to_native(
                value["x"], value["y"], value["width"], value["height"],
                crop.shape[1], crop.shape[0], y_lo, W, H)
            ed = st.session_state.manual_edits.setdefault(name, empty_edits())
            sets_ed = ed[side_key]
            if not sets_ed:
                sets_ed.append([])
            idx = (active if active is not None and active < len(sets_ed)
                   else len(sets_ed) - 1)
            sets_ed[idx].append((x_nat, y_nat))
            st.rerun()

        if sets_cur:
            parts = []
            for i, s in enumerate(sets_cur):
                if s:
                    xs = [p[0] for p in s]
                    parts.append(f"set {i + 1}: {len(s)} pts "
                                 f"(x {min(xs):.0f}–{max(xs):.0f})")
                else:
                    parts.append(f"set {i + 1}: empty")
            st.caption(" · ".join(parts))
            if active is not None and sets_cur[active]:
                st.caption(f"set {active + 1} points: " + ", ".join(
                    f"({ax:.0f}, {ay:.0f})" for ax, ay in sets_cur[active]))
        b = st.columns(3)
        has_active_pts = active is not None and bool(sets_cur[active])
        if b[0].button("Undo last point", key=f"undo_{name}",
                       disabled=not has_active_pts):
            ed = st.session_state.manual_edits.get(name)
            if ed and active is not None and active < len(ed[side_key]):
                if ed[side_key][active]:
                    ed[side_key][active].pop()
                if not ed[side_key][active]:
                    ed[side_key].pop(active)
                    st.session_state[ver_key] = ver + 1
            st.rerun()
        if b[1].button("Delete set", key=f"delset_{name}",
                       disabled=active is None):
            ed = st.session_state.manual_edits.get(name)
            if ed and active is not None and active < len(ed[side_key]):
                ed[side_key].pop(active)
                st.session_state[ver_key] = ver + 1
            st.rerun()
        if b[2].button("Clear all edits", key=f"clrall_{name}",
                       disabled=not has_edits(st.session_state.manual_edits.get(name))):
            st.session_state.manual_edits.pop(name, None)
            st.session_state[ver_key] = ver + 1
            st.rerun()

        n1, n2 = st.columns(2)
        new_nt = n1.number_input("Nudge top line (px, + = down)",
                                 value=float(edits["nudge_top"]), step=0.5,
                                 key=f"nudgetop_{name}")
        new_nb = n2.number_input("Nudge bottom line (px, + = down)",
                                 value=float(edits["nudge_bot"]), step=0.5,
                                 key=f"nudgebot_{name}")
        if new_nt != edits["nudge_top"] or new_nb != edits["nudge_bot"]:
            ed = st.session_state.manual_edits.setdefault(name, empty_edits())
            ed["nudge_top"], ed["nudge_bot"] = float(new_nt), float(new_nb)
            st.rerun()

        et, eb = rep.get("edited_top"), rep.get("edited_bot")
        if et is not None or eb is not None:
            edited_cols = np.zeros(W, dtype=bool)
            if et is not None:
                edited_cols |= et
            if eb is not None:
                edited_cols |= eb
            st.caption(f"edited columns: {int(edited_cols.sum())}; "
                       f"valid after re-QC: {int((edited_cols & res.valid).sum())}")


def _render_per_image_stats(reps: list[dict]) -> None:
    """Compact per-image table so individual stats sit next to the group stats.

    Stats are taken along each picture's own fibre (valid columns, before
    registration), so they answer "how thick is this picture's fibre and how
    much does it vary along its length" — distinct from the group panel's
    between-replicate numbers.
    """
    rows = []
    for rep in reps:
        mr = rep["mr"]
        d = mr.diameter_um
        ok = bool(np.isfinite(d).any())
        rows.append({
            "image": mr.name,
            "mean Ø (µm)": float(np.nanmean(d)) if ok else np.nan,
            "std (µm)": float(np.nanstd(d)) if ok else np.nan,
            "median Ø (µm)": float(np.nanmedian(d)) if ok else np.nan,
            "coverage": mr.res.coverage,
        })
    st.markdown("**Per-image statistics** — each picture's own fibre, before "
                "registration")
    st.dataframe(
        pd.DataFrame(rows), width="stretch", hide_index=True,
        column_config={
            "mean Ø (µm)": st.column_config.NumberColumn(format="%.2f"),
            "std (µm)": st.column_config.NumberColumn(
                format="%.2f",
                help="Thickness variation along this picture's fibre."),
            "median Ø (µm)": st.column_config.NumberColumn(format="%.2f"),
            "coverage": st.column_config.NumberColumn(format="percent"),
        })


def _render_group(reps: list[dict], cfg: CONFIG, group_label: str | None,
                  tmap: dict) -> None:
    st.subheader("Group panel — registered mean ± std")
    profiles, dropped = [], []
    for rep in reps:
        mr = rep["mr"]
        res, bnd = mr.res, mr.bnd
        if res.band_mismatch:
            dropped.append(f"{mr.name} (band_mismatch)")
            continue
        if res.coverage < cfg.min_coverage:
            dropped.append(f"{mr.name} (coverage {res.coverage:.0%} < {cfg.min_coverage:.0%})")
            continue
        W = rep["rgb"].shape[1]
        span = slice(bnd.x0, bnd.x1 + 1)
        profiles.append({
            "replicate": mr.replicate if mr.replicate is not None else rep["idx"] + 1,
            "coverage": res.coverage,
            "x": np.arange(W)[span].astype(float),
            "diameter_px_raw": res.diameter_raw[span].astype(float),
            "diameter_px_smooth": res.diameter_smooth[span].astype(float),
            "valid": res.valid[span].astype(bool),
        })
    if dropped:
        st.caption("QC-dropped from registration: " + "; ".join(dropped))
    if not profiles:
        st.warning("No replicate passed QC (coverage / band_mismatch); nothing to register.")
        _render_per_image_stats(reps)
        return
    # sort with the same key register_sample uses, so zip(profiles, shifts)
    # below is order-aligned (replicate-keyed dicts would silently collide for
    # ungrouped uploads that share the idx+1 fallback number)
    profiles.sort(key=lambda p: p["replicate"])
    try:
        table, shifts, summary = register_sample(profiles, cfg)
    except Exception as exc:  # noqa: BLE001
        st.error(f"Registration failed: {exc}")
        _render_per_image_stats(reps)
        return
    assert all(s["replicate"] == p["replicate"] for p, s in zip(profiles, shifts))

    rep_curves = [(p["replicate"], p["x"] + s["shift_px"],
                   p["diameter_px_smooth"] / cfg.ppu)
                  for p, s in zip(profiles, shifts)]
    fig = _group_fig(table, group_label or "uploaded", rep_curves)
    st.pyplot(fig)
    plt.close(fig)

    shift_txt = ", ".join(f"rep {s['replicate']}: {s['shift_px']:+.0f} px"
                          for s in shifts)
    st.caption(
        "Each replicate's profile is shifted horizontally (cross-correlation "
        "against the first replicate) so the same physical stretch of fibre "
        "lines up before averaging — the photos never frame the fibre "
        f"identically. Applied shifts: {shift_txt}. Thin lines = individual "
        "replicates after alignment; bold line = pointwise mean across "
        "replicates; band = ±1 std.")

    cv = summary["cv"]
    c = st.columns(6)
    c[0].metric("group mean Ø", f"{summary['mean_um']:.2f} µm",
                help="Mean diameter across the replicates of this group, "
                     "averaged along the aligned overlap region (where all "
                     "replicates are present after registration).")
    c[1].metric("between-replicate std", f"{summary['std_um']:.2f} µm",
                help="At each aligned position, the std of the diameter across "
                     "replicates; this is its average along the fibre. It "
                     "measures replicate-to-replicate disagreement, not "
                     "thickness variation along the fibre.")
    c[2].metric("CV", f"{cv:.3f}" if np.isfinite(cv) else "—",
                help="between-replicate std / group mean (dimensionless).")
    c[3].metric("n reps used", f"{summary['n_replicates_used']}")
    c[4].metric("overlap", f"{summary['overlap_px']} px",
                help="Number of aligned columns where every replicate has data.")
    c[5].metric("registration",
                "uncertain" if summary["registration_uncertain"] else "ok",
                help="'uncertain' = the cross-correlation peak was below "
                     "min_corr for at least one replicate, so its shift was "
                     "reset to 0.")

    _render_per_image_stats(reps)
    st.caption(
        "Per-image mean/median/std are computed along each picture's own "
        "fibre (valid columns, no registration): std there = thickness "
        "variation within one picture. The group numbers above instead "
        "average the registered replicates, so their std = disagreement "
        "between replicates.")

    _render_tensile(group_label, summary["mean_um"], cfg, tmap)


def _manual_break_control(df, mean_um, cfg: CONFIG, group_label: str,
                          auto, breaks: dict) -> object:
    """Checkbox + strain slider to override the auto-detected fracture.

    Returns the ``TensileResult`` to display: ``auto`` when manual mode is off,
    or a recompute pinned to the user-chosen sample when on. Records the chosen
    sample index in ``breaks[group_label]`` (or clears it) so the export matrix
    honours the same break.
    """
    on = st.checkbox(
        "Set break point manually", key=f"tbreak_on_{group_label}",
        help="Override the auto-detected fracture; drag to the strain where the "
             "fibre actually breaks. The export uses this break too.")
    strain = np.asarray(auto.strain, dtype=float)
    finite = strain[np.isfinite(strain)] * 100.0
    lo = round(float(finite.min()), 3) if finite.size else 0.0
    hi = round(float(finite.max()), 3) if finite.size else 0.0
    if not on or finite.size < 2 or hi <= lo:
        breaks.pop(group_label, None)
        if on:
            st.caption("Not enough distinct points to set a manual break.")
        return auto

    auto_pct = float(auto.strain_at_break * 100.0)
    step = max(round((hi - lo) / 200.0, 3), 0.01)
    skey = f"tbreak_strain_{group_label}"
    if skey not in st.session_state or not (lo <= st.session_state[skey] <= hi):
        st.session_state[skey] = float(min(max(auto_pct, lo), hi))
    sel = st.slider("break at strain (%)", min_value=lo, max_value=hi,
                    step=step, key=skey)

    diffs = np.abs(strain - sel / 100.0)
    diffs[~np.isfinite(diffs)] = np.inf
    brk = int(np.argmin(diffs))
    breaks[group_label] = brk
    st.caption(f"Manual break at {sel:.2f}% strain (auto was {auto_pct:.2f}%). "
               "Untick to restore automatic detection.")
    return compute_tensile(df, diameter_um=mean_um,
                           gauge_length_mm=cfg.gauge_length_mm, cfg=cfg,
                           break_index=brk)


def _render_tensile(group_label: str | None, mean_um: float, cfg: CONFIG,
                    tmap: dict) -> None:
    """Stress-strain subsection for the loaded fibre, joined by its group key.

    ``tmap`` is the resolved ``{group: source}`` map from the sidebar (folder or
    uploads). Uses the group's registered mean diameter for the cross-section,
    so stress = force / area is tied to this fibre's own measurement. Degrades
    gracefully: no data -> a hint; no matched fibre -> a caption; any read/parse
    error -> ``st.error`` rather than a crashed page.
    """
    if not tmap or group_label is None:
        st.caption("Set a tensile data folder or upload tensile files in the "
                   "sidebar to see stress–strain.")
        return
    if group_label not in tmap:
        st.caption(f"No tensile file matched fibre {group_label}.")
        return

    try:
        df = read_trace(tmap[group_label])
        auto = compute_tensile(df, diameter_um=mean_um,
                               gauge_length_mm=cfg.gauge_length_mm, cfg=cfg)
        st.subheader("Tensile (stress–strain)")

        # manual break point: override the auto-detected fracture by dragging it
        # along the strain axis; the choice is stored per fibre so the export
        # matrix uses the same break the user picked here.
        breaks = st.session_state.setdefault("tensile_breaks", {})
        res = _manual_break_control(df, mean_um, cfg, group_label, auto, breaks)

        fig = _tensile_fig(res)
        st.pyplot(fig)
        plt.close(fig)

        # Six metrics in six narrow columns: shrink the value font so figures
        # like "161.00 mN" are not clipped with an ellipsis. Scoped to this
        # container's key so the diameter metrics elsewhere keep their size.
        st.markdown(
            "<style>.st-key-tensile_metrics [data-testid='stMetricValue']"
            "{font-size:1.1rem;}.st-key-tensile_metrics "
            "[data-testid='stMetricLabel']{font-size:0.8rem;}</style>",
            unsafe_allow_html=True)
        with st.container(key="tensile_metrics"):
            c = st.columns(6)
            c[0].metric("breaking force", f"{res.fmax_n * 1000:.2f} mN",
                        help="Peak of the load trace (Fmax).")
            c[1].metric("tensile strength",
                        f"{res.tensile_strength_pa / 1e6:.1f} MPa"
                        if np.isfinite(res.tensile_strength_pa) else "—",
                        help="Fmax / cross-sectional area.")
            c[2].metric("extension at break",
                        f"{res.extension_at_break_mm:.3f} mm")
            c[3].metric("strain at break", f"{res.strain_at_break * 100:.2f} %")
            c[4].metric("Young's modulus",
                        f"{res.youngs_modulus_pa / 1e9:.2f} GPa"
                        if np.isfinite(res.youngs_modulus_pa) else "—")
            c[5].metric("toughness",
                        f"{res.toughness_j_m3 / 1e6:.2f} MJ/m³"
                        if np.isfinite(res.toughness_j_m3) else "—")

        area_um2 = res.area_m2 * 1e12 if np.isfinite(res.area_m2) else np.nan
        d_txt = (f"{res.diameter_um:.2f} µm" if res.diameter_um is not None
                 else "—")
        a_txt = f"{area_um2:.1f} µm²" if np.isfinite(area_um2) else "—"
        flag_txt = ("; flags: " + ", ".join(res.flags)) if res.flags else ""
        st.caption(f"Diameter used: {d_txt}; cross-section area: {a_txt}"
                   f"{flag_txt}.")
        st.caption(
            "stress = force / area (area from this fibre's mean measured "
            "diameter), strain = displacement / L₀, modulus = steepest initial "
            "slope of the curve, toughness = shaded area under the curve.")
    except Exception as exc:  # noqa: BLE001
        st.error(f"Tensile analysis failed for {group_label}: {exc}")


def _render_export_batch(reps: list[dict], cfg: CONFIG, group_label: str | None,
                         folder: str | None, tmap: dict,
                         cfg_items: tuple) -> None:
    st.divider()
    st.subheader("Export & batch")
    out_folder = st.text_input("Output folder",
                               value=st.session_state.get("out_folder", "./fibrecv_output"))
    st.session_state.out_folder = out_folder

    col_exp, col_batch = st.columns(2)

    edited_names = sorted(n for n, e in st.session_state.get("manual_edits", {}).items()
                          if has_edits(e))

    with col_exp:
        st.markdown("**Export current group**")
        st.caption("Writes overlays/, per_image/* and per_sample/* for the loaded group.")
        if edited_names:
            st.info("Manual edits active for: " + ", ".join(edited_names)
                    + " — they ARE included in this export.")
        if st.button("Export current group", disabled=not reps, width="stretch"):
            try:
                with st.spinner("Writing output tree…"):
                    g = export_group(reps, out_folder, cfg)
                st.success(f"Exported group {g} → {out_folder}")
            except Exception as exc:  # noqa: BLE001
                st.error(f"Export failed: {exc}")

    image_uploads = st.session_state.get("image_uploads") or []
    folder_ok = folder is not None and Path(folder or "").is_dir()

    with col_batch:
        st.markdown("**Export all groups**")
        st.caption("Measures every group in this session — the whole folder, or "
                   "all uploaded images — and writes the full output tree + "
                   "master_summary.csv (one row per group).")
        jobs = (st.number_input("parallel jobs", min_value=1, max_value=16,
                                 value=4, step=1) if folder_ok else 4)
        if edited_names and folder_ok:
            st.warning("Folder export recomputes every image from disk and "
                       "IGNORES the manual edits made here.")
        elif edited_names:
            st.info("Manual edits on the loaded group ARE kept; other groups are "
                    "measured fresh.")
        can_all = folder_ok or bool(image_uploads)
        if not can_all:
            st.caption("Needs a *Local folder* source or uploaded images.")
        if st.button("Export all groups", disabled=not can_all, width="stretch"):
            try:
                if folder_ok:
                    images = discover_images(folder)
                    if not images:
                        st.warning("No image files found in the folder.")
                        return
                    n = len(images)
                    prog = st.progress(0.0, text=f"Measuring 0/{n}…")

                    def _cb(frac: float, _r: dict) -> None:
                        prog.progress(frac, text=f"Measuring {int(round(frac * n))}/{n}…")

                    with st.spinner("Exporting all groups…"):
                        master, results = run_batch(images, out_folder, cfg, int(jobs), _cb)
                    prog.empty()
                    n_err = sum(1 for r in results if "error" in r)
                    st.success(f"Exported {len(results)} images, {n_err} errors "
                               f"→ {out_folder}")
                else:
                    grouped = _grouped_reps_from_uploads(image_uploads, cfg_items,
                                                         reps, group_label)
                    with st.spinner(f"Exporting {len(grouped)} groups "
                                    f"({len(image_uploads)} images)…"):
                        exported = export_all_groups(grouped, out_folder, cfg)
                    mpath = Path(out_folder) / "summary" / "master_summary.csv"
                    master = pd.read_csv(mpath) if mpath.exists() else pd.DataFrame()
                    st.success(f"Exported {len(exported)} groups → {out_folder}")
                if not master.empty:
                    st.dataframe(master, width="stretch")
                    st.download_button(
                        "Download master_summary.csv",
                        master.to_csv(index=False).encode(),
                        file_name="master_summary.csv", mime="text/csv")
                else:
                    st.info("No groups aggregated (check filenames / coverage).")
            except Exception as exc:  # noqa: BLE001
                st.error(f"Export all groups failed: {exc}")
                st.code(traceback.format_exc())

    st.divider()
    st.markdown("**Export tensile matrix (all fibres)**")
    st.caption("Measures every image's diameter and joins each fibre's tensile "
               "metrics into one row per fibre.")
    img_ok = folder_ok or bool(image_uploads)
    tns_ok = bool(tmap)
    if not (img_ok and tns_ok):
        st.caption("Needs an image source (Local folder or uploaded images) and "
                   "tensile data (folder or uploads) set in the sidebar.")
    if st.button("Build & export tensile matrix",
                 disabled=not (img_ok and tns_ok), width="stretch"):
        try:
            if folder_ok:
                # disk batch: writes the full output tree + master_summary
                images = discover_images(folder)
                if not images:
                    st.warning("No image files found in the folder.")
                    return
                n = len(images)
                prog = st.progress(0.0, text=f"Measuring 0/{n}…")

                def _mcb(frac: float, _r: dict) -> None:
                    prog.progress(frac, text=f"Measuring {int(round(frac * n))}/{n}…")

                with st.spinner("Measuring all images…"):
                    master, _ = run_batch(images, out_folder, cfg, int(jobs), _mcb)
                prog.empty()
                if master.empty:
                    st.info("No groups aggregated (check filenames / coverage).")
                    return
                diameters = dict(zip(master["group"].astype(str),
                                     master["mean_um"].astype(float)))
            else:
                # uploaded images: measure in-memory (no disk batch)
                with st.spinner(f"Measuring {len(image_uploads)} uploaded images…"):
                    diameters = _diameters_from_uploads(image_uploads, cfg_items)
                if not diameters:
                    st.info("No image groups could be measured from the uploads "
                            "(check filenames / coverage).")
                    return

            breaks = st.session_state.get("tensile_breaks", {})
            matrix = build_matrix(diameters, tmap, cfg, breaks=breaks)
            out_path = Path(out_folder) / "summary" / "tensile_matrix.csv"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            matrix.to_csv(out_path, index=False)
            st.success(f"Wrote {out_path}")
            st.dataframe(matrix, width="stretch")
            st.download_button(
                "Download tensile_matrix.csv",
                matrix.to_csv(index=False).encode(),
                file_name="tensile_matrix.csv", mime="text/csv")
        except Exception as exc:  # noqa: BLE001
            st.error(f"Tensile matrix export failed: {exc}")
            st.code(traceback.format_exc())


# --------------------------------------------------------------------------- #
# App entry point                                                             #
# --------------------------------------------------------------------------- #
def main() -> None:
    st.set_page_config(page_title="fibrecv — fibre diameter GUI", layout="wide")
    st.session_state.setdefault("cfg_dict", DEFAULTS.as_dict())
    st.session_state.setdefault("form_version", 0)
    st.session_state.setdefault("manual_edits", {})  # image name -> edits dict

    st.title("fibrecv — fibre diameter detection")
    st.caption("Local preview / tuning / batch / export over the validated pipeline. "
               f"Strictness knob edge_z = {st.session_state.cfg_dict['edge_z']}.")

    cfg_items = _cfg_items(st.session_state.cfg_dict)
    cfg = _cfg_from_items(cfg_items)

    # sidebar
    reps, group_label, folder = _load_reps(cfg_items)
    _param_form()
    tensile = _tensile_controls()

    # tensile-specific config: the diameter knobs stay as tuned; only the
    # strain scale and modulus-fit width come from the tensile controls
    tcfg = replace(cfg, gauge_length_mm=tensile["gauge_length_mm"],
                   modulus_window=tensile["modulus_window"])

    # apply manual boundary edits once, right here: every downstream consumer
    # (overlay, profile plot, group registration, export) reads rep["mr"], so
    # corrections flow everywhere; the cached MeasureResult is never mutated
    for rep in reps:
        edits = st.session_state.manual_edits.get(rep["name"])
        if has_edits(edits):
            rep["mr"], rep["edited_top"], rep["edited_bot"] = \
                apply_manual_edits(rep["mr"], edits, cfg)

    # main area
    if not reps:
        st.info("Pick a folder + group, or upload images, to begin.")
        return

    st.subheader(f"Replicates — group {group_label}" if group_label else "Replicates (uploaded)")
    tabs = st.tabs([f"_{r['mr'].replicate}" if r["mr"].replicate is not None else r["name"]
                    for r in reps])
    for tab, rep in zip(tabs, reps):
        with tab:
            _render_replicate(rep, cfg)

    _render_group(reps, tcfg, group_label, tensile["tmap"])
    _render_export_batch(reps, tcfg, group_label, folder, tensile["tmap"], cfg_items)


if __name__ == "__main__":
    main()
