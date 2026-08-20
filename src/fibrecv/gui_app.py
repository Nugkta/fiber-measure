"""Local Streamlit GUI for tuning, previewing, batch-processing and exporting.

Dependencies
------------
``streamlit``, ``numpy``, ``pandas``, ``matplotlib`` (Agg), ``imageio`` plus the
reused ``fibrecv`` package: ``compute`` (pure detection core), ``config``,
``io_utils`` (incl. ``discover_multiangle``/``parse_multiangle_name``),
``overlay`` (``render_overlay``), ``register``, ``measure``
(``write_measurement``), ``xsection`` (``build_part_stack``,
``fit_ellipse_projections``, split-half/pair helpers) and
``run_measure``/``run_aggregate``/``run_xsection`` (batch + aggregate stages).
Imports are absolute (``fibrecv.*``) because Streamlit runs this file as a script.
The sibling ``.streamlit/config.toml`` supplies the Streamlit-engine half of the
theme (base colours, radii, ``toolbarMode``); this file's own ``_CSS`` constant
covers what that engine cannot reach (header/chip markup, section cards, the
jump menu) -- the two are one visual design and must be read together.

Inputs
------
- The analysis-mode radio (``analysis_mode``: Replicates | Multi-angle
  cross-section) picking which of the two main-area layouts renders; entering
  multi-angle mode once re-applies the calibrated bright defaults.
- Images from a local folder OR any number of uploaded files; both are
  auto-grouped via ``parse_name``'s trailing-numbers rule, with unparseable
  names collected in an "ungrouped" bucket. In multi-angle mode the same two
  sources are grouped by ``parse_multiangle_name`` instead (condition / fibre /
  part selectors), with a manual µm/px scale replacing the ppu calibration.
- The image-mode selector (``feature_mode``: desat | bright — switching modes
  re-applies that mode's calibrated ``edge_frac``/``k_band`` defaults, bright
  via ``run_measure.BRIGHT_DEFAULTS``), the three boundary knobs
  (``edge_z``/``edge_frac``/``wcol``), the ``ppu`` calibration, and the
  anomaly-flag thresholds (``jump_thresh_px``/``gap_frac``/``step_frac``/
  ``rep_dev_frac``) plus the ``anomaly_exclude`` checkbox, edited in a sidebar
  form and applied on demand; all other ``CONFIG`` fields stay at the
  validated defaults.
- An output-folder path for export/batch.

Output
------
- Live, in-memory preview: full-res boundary overlays, per-replicate diameter
  profiles, and a registered mean+/-std group curve -- all redrawn when the user
  changes parameters and clicks Apply (no disk writes for preview). Anomalous
  measurements are surfaced as a "⚠" replicate-tab prefix, the anomaly names in
  the amber flags badge and the per-image stats' ``anomalies`` column; with
  ``anomaly_exclude`` on they drop the replicate from registration via the same
  ``anomaly.exclusion_reason`` policy the CLI aggregator uses.
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
- In multi-angle mode instead: the six angle views of one (fibre, part) as
  tabs, their aligned width stack, the per-position ellipse cross-section
  (area +/- split-half error, axis ratio, orientation, rms residual, drawn
  cross-section) and a one-button batch, scoped to the selected condition
  (every fibre/part under it, not just the one previewed), that runs the
  per-image measure pass and ``run_xsection`` back to back with a numeric
  ``--scale-source`` (the sidebar's manual µm/px field, never ``cfg.ppu``).
- A "clean lab" light UI over that data: a slim violet-accent header + state
  chip (``_render_header``), numbered card sections -- 01 Replicates,
  02 Group panel, 03 Tensile, 04 Export & batch (or, in multi-angle mode,
  01 Angles, 02 Cross-section, 03 Batch & export) -- each an
  ``st.container(key="card_...")`` framed by ``_CSS``, with a fixed jump menu
  (``_render_jump_menu``, hidden below 1200px) linking to their anchors.

Pos
---
Thin front-end over the validated pipeline; adds no detection logic. Designed to
run on the user's local Mac/Windows machine (``fibrecv-gui`` or
``streamlit run src/fibrecv/gui_app.py``), reading a copied images folder. The
heavy compute reuses ``compute.compute_measurement`` so preview == CLI output.
The styling layer is self-contained and additive over that data path: ``_CSS``
(injected once via ``_inject_css()`` at the top of ``main()``) styles the
header/chip, section cards and jump menu by targeting ``data-testid``/
``.st-key-*`` hooks on ordinary widgets rather than replacing them (metric
cards are real ``st.metric`` underneath, kept AppTest-inspectable); ``_fmt``
centralises the one value-formatting policy (precision + unit suffix per
``kind``, ``None``/NaN -> "—") used by every metric, caption and the header
chip; ``_styled_fig``/``_ACCENT``/``_REP_CYCLE`` give every matplotlib figure
the same muted, transparent-background look. None of this touches
``CONFIG``, session-state keys, ``compute_measurement`` or the export/batch
code paths.
"""

from __future__ import annotations

import html
import io
import json
import re
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, replace
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

from fibrecv import run_aggregate, run_xsection  # noqa: E402
from fibrecv.anomaly import detect_replicate_outliers, exclusion_reason  # noqa: E402
from fibrecv.compute import compute_measurement  # noqa: E402
from fibrecv.config import CONFIG  # noqa: E402
from fibrecv.io_utils import (  # noqa: E402
    IMAGE_SUFFIXES, discover_images, discover_multiangle, natural_key,
    parse_multiangle_name, parse_name)
from fibrecv.io_utils import load_rgb as _io_load_rgb  # noqa: E402
from fibrecv.manual_edit import (  # noqa: E402
    apply_manual_edits, display_to_native, empty_edits, has_edits)
from fibrecv.measure import write_measurement  # noqa: E402
from fibrecv.overlay import GREY, WHITE, mark_anchors, render_overlay  # noqa: E402
from fibrecv.register import register_sample  # noqa: E402
from fibrecv.run_measure import BRIGHT_DEFAULTS, _lib_versions, _worker  # noqa: E402
from fibrecv.tensile import (  # noqa: E402
    build_matrix, compute_tensile, discover_tensile, discover_tensile_files,
    read_trace)
from fibrecv.xsection import (  # noqa: E402
    PartStack, XsecFit, build_part_stack, fit_ellipse_projections,
    pair_differences, split_half_area)
from streamlit_image_coordinates import streamlit_image_coordinates  # noqa: E402

DEFAULTS = CONFIG()  # never mutated; the source of widget defaults + reset target

# --- visible parameters: the three knobs that move the detected boundary,
# plus the ppu calibration. Everything else in CONFIG stays at the validated
# defaults (CLI keeps full control). spec = (name, label, kind, help, step,
# lo, hi, fmt); ``name`` stays the CONFIG field key, ``label`` is the display
# text only. ---
PARAM_SPECS: list[tuple] = [
    ("edge_z", "Boundary tightness (edge_z)", "slider",
     "Where on the fibre wall the boundary line is drawn. In desat mode this "
     "is the main knob: higher sits further inside the fibre (thinner "
     "reading), lower sits further outside (thicker reading). If the line "
     "cuts into the fibre, lower it; if it sits in the shadow outside the "
     "fibre, raise it. Recommended 3-5, default 4.0. In bright mode it is "
     "only the noise floor under edge_frac — leave it at 4.0 there.",
     0.5, 1.0, 12.0, "%.1f"),
    ("edge_frac", "Relative threshold (edge_frac)", "float",
     "Fraction of the wall's own height where the boundary is placed. "
     "Bright mode: this is the main knob — the line sits where the wall "
     "reaches edge_frac of its full height (calibrated 0.30; higher = "
     "thinner reading). Desat mode: only a faint-fibre guard — it caps the "
     "crossing level for walls too weak to reach edge_z above background; "
     "leave at 0.65 unless a very faint fibre loses its boundary. Switching "
     "the image mode resets this to that mode's calibrated default.",
     0.05, 0.0, 1.0, "%.2f"),
    ("wcol", "Smoothing width (wcol)", "int",
     "Horizontal smoothing width in pixels. If the boundary line is jittery "
     "or ragged, raise it (61-81) for a smoother line; lower it (15-25) to "
     "preserve fine thickness variation — too high flattens real variation. "
     "Default 41.",
     1, 1, 201, None),
    ("ppu", "Pixels per micron (ppu)", "float",
     "Calibration: camera pixels per micron; diameters in µm = pixels / ppu, "
     "so this scales every µm number in the app and exports (pixel values "
     "are unaffected). Measure it with a stage micrometer for your "
     "microscope + camera combination. Default 1.3680 (the original "
     "calibrated setup). Not used in Multi-angle cross-section mode, which "
     "takes its scale from the sidebar's µm/px field.",
     0.001, 0.1, 10.0, "%.4f"),
    ("jump_thresh_px", "Edge jump threshold (px)", "float",
     "edge_jump anomaly: a detected edge that moves more than this many "
     "pixels between neighbouring measured columns (tilt-corrected) is "
     "flagged — usually a boundary latching onto a reflection or shadow. "
     "Default 20 (calibrated on the MasP2 set).",
     1.0, 1.0, 100.0, "%.0f"),
    ("gap_frac", "Missing-stretch fraction", "float",
     "large_gap anomaly: flag the image when the longest stretch of "
     "unmeasurable columns exceeds this fraction of the fibre span — long "
     "gaps usually mean focus or lighting problems. Default 0.10.",
     0.01, 0.01, 1.0, "%.2f"),
    ("step_frac", "Diameter step fraction", "float",
     "diameter_step anomaly: flag when the diameter level shifts by more "
     "than this fraction of the image's median diameter between adjacent "
     "windows — a sudden persistent shift, not gradual taper. Default 0.25 "
     "(calibrated on the MasP2 set; smooth real variation reaches ~0.20).",
     0.01, 0.01, 1.0, "%.2f"),
    ("rep_dev_frac", "Replicate deviation fraction", "float",
     "replicate_outlier: flag a replicate whose median diameter deviates "
     "from the group median by more than this fraction. Advisory only — it "
     "never excludes, since diameter genuinely varies along a fibre. "
     "Default 0.25.",
     0.05, 0.05, 2.0, "%.2f"),
    ("anomaly_exclude", "Exclude flagged images from group stats", "bool",
     "When on, images with edge_jump / large_gap / diameter_step anomalies "
     "are dropped from registration and the group mean, like band_mismatch. "
     "Off = anomalies are advisory badges only. replicate_outlier never "
     "excludes either way.",
     None, None, None, None),
]

# names of int-typed visible fields (so widgets return int, not float)
_INT_FIELDS = {name for (name, label, kind, *_rest) in PARAM_SPECS if kind == "int"}

# feature_mode selector: internal value -> display label. Bright mode starts
# from the study-03 calibrated BRIGHT_DEFAULTS (edge_frac 0.30 / k_band 6.0),
# mirroring run_measure.build_config; desat keeps the dataclass defaults.
_MODE_LABELS: dict[str, str] = {
    "desat": "Pale fibre on saturated background (desat)",
    "bright": "Bright fibre on dark background (bright)",
}


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


# internal flag name -> human-readable display label. The raw snake_case names
# stay in the meta JSON / per_image_summary.csv (stable machine schema); every
# user-facing surface (badge, stats table, dropped caption) shows these labels.
_FLAG_LABELS: dict[str, str] = {
    "low_confidence": "low confidence",
    "band_mismatch": "band mismatch",
    "edge_jump": "edge jump",
    "large_gap": "large gap",
    "diameter_step": "diameter step",
    "replicate_outlier": "deviant replicate",
}

# the flags badge's help tooltip: every flag this app can raise, in one place
_FLAGS_HELP = (
    "Everything this app can warn about for one image:\n"
    "- **low confidence** — too few columns could be measured (low coverage) "
    "or the band detection itself is doubtful; treat this image's numbers "
    "with care.\n"
    "- **band mismatch** — the measured diameter disagrees badly with the "
    "coarse band width, so the detector probably locked onto a reflection or "
    "shadow instead of the fibre walls. Always excluded from group stats.\n"
    "- **edge jump** — a boundary line jumps suddenly between neighbouring "
    "columns, usually because it briefly grabbed a reflection or shadow.\n"
    "- **large gap** — a long stretch of the fibre could not be measured at "
    "all (focus or lighting problems).\n"
    "- **diameter step** — the diameter shifts to a different level "
    "mid-image instead of varying smoothly; often a focus change or an "
    "overlapping object.\n"
    "- **deviant replicate** — this image's median diameter sits far from "
    "its group's median. Advisory only: diameter genuinely varies along a "
    "fibre, so this never excludes the image.\n\n"
    "The last four are advisory photo-quality warnings: consider re-taking "
    "the photo, or fix the boundary under *Edit boundaries*. They exclude an "
    "image from the group stats only when *Exclude flagged images* is on "
    "(deviant replicate never excludes)."
)


def _flag_labels(flags: list[str]) -> list[str]:
    """Map internal flag names to their display labels (unknown -> as-is)."""
    return [_FLAG_LABELS.get(f, f) for f in flags]


def _friendly_reason(reason: str) -> str:
    """Rewrite raw flag names inside an exclusion reason for display."""
    for raw, label in _FLAG_LABELS.items():
        reason = reason.replace(raw, label)
    return reason


# kind -> (format spec, unit suffix); the single source of display precision
# for every metric/caption number in the app, so the same quantity is never
# shown at two different precisions.
_FMT_KINDS: dict[str, tuple[str, str]] = {
    "um":     ("{:.2f}", " µm"),
    "um2":    ("{:.1f}", " µm²"),
    "cv":     ("{:.3f}", ""),
    "coverage": ("{:.0%}", ""),
    "px":     ("{:.0f}", " px"),
    # px2 = the same unit at 2 dp, for small pixel quantities where the
    # 0-dp "px" kind would hide the value's relation to a threshold (an rms
    # residual of 6.27 px vs the 6.0 px flag would both render as "6 px")
    "px2":    ("{:.2f}", " px"),
    "deg":    ("{:.1f}", "°"),
    "ppu":    ("{:.4f}", ""),
    "tilt":   ("{:.4f}", ""),
    "mm":     ("{:.3f}", " mm"),
    "mN":     ("{:.2f}", " mN"),
    "GPa":    ("{:.2f}", " GPa"),
    "MJ/m3":  ("{:.2f}", " MJ/m³"),
    "MPa":    ("{:.1f}", " MPa"),
    "one_dp": ("{:.1f}", ""),
    "pct100": ("{:.2f}", " %"),
    "int":    ("{:.0f}", ""),
}


def _fmt(value, kind: str, dash: str = "—") -> str:
    """Format ``value`` per ``kind``'s entry in ``_FMT_KINDS``; ``None``/NaN -> ``dash``.

    Centralises the precision + unit suffix used across metric cards, table
    captions and the header chip, so a missing value always renders as the
    same em-dash instead of "nan" or an inconsistent decimal count.
    """
    if value is None:
        return dash
    try:
        fval = float(value)
    except (TypeError, ValueError):
        return dash
    if not np.isfinite(fval):
        return dash
    fmt, suffix = _FMT_KINDS[kind]
    return fmt.format(fval) + suffix


def _safe_key(name: str) -> str:
    """Sanitise a free-text name (e.g. an image filename) into a valid
    Streamlit container key: alphanumerics, ``_`` and ``-`` only. Only used to
    build container keys — never changes ``rep["name"]`` or any session-state
    key derived from it."""
    return re.sub(r"[^A-Za-z0-9_-]", "_", name)


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
        "--rep-dev-frac", str(cfg.rep_dev_frac),
    ] + (["--anomaly-exclude"] if cfg.anomaly_exclude else []))
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
            "--rep-dev-frac", str(cfg.rep_dev_frac),
        ] + (["--anomaly-exclude"] if cfg.anomaly_exclude else []))
    return groups


def run_batch(
    image_paths: list[Path],
    out_root: str | Path,
    cfg: CONFIG,
    jobs: int,
    progress_cb: Callable[[float, dict], None] | None = None,
    aggregate: bool = True,
) -> tuple[pd.DataFrame, list[dict]]:
    """Measure every image in-process, then (optionally) aggregate the folder.

    Tries a ``ProcessPoolExecutor`` (reusing the picklable ``run_measure._worker``)
    and falls back to sequential if the pool cannot start (Windows ``spawn``
    safety); both paths report progress through ``progress_cb(frac, result)``.
    Writes the per-image output tree (per_image/*, overlays/, run_config.json)
    always; with ``aggregate=True`` (default) also runs ``run_aggregate.main``
    to build per_sample/* and summary/master_summary.csv, returning
    ``(master_summary_df, per_image_results)``. ``aggregate=False`` skips that
    stage-2 pass entirely -- the multi-angle batch uses this, since stage 2's
    group registration is meaningless for angle data (``run_xsection`` is its
    own aggregation stage). The DataFrame read-back below is unconditional,
    so with ``aggregate=False`` the returned DataFrame is whatever
    ``summary/master_summary.csv`` already holds in ``out_root`` -- empty for
    a fresh out_root, but possibly stale rows left over from a previous
    ``aggregate=True`` run against the same out_root. Callers using
    ``aggregate=False`` should ignore the returned DataFrame.
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
    if aggregate:
        run_aggregate.main([
            "--out", str(out_root), "--all",
            "--ppu", str(cfg.ppu), "--max-shift", str(cfg.max_shift),
            "--min-corr", str(cfg.min_corr), "--min-coverage", str(cfg.min_coverage),
            "--rep-dev-frac", str(cfg.rep_dev_frac),
        ] + (["--anomaly-exclude"] if cfg.anomaly_exclude else []))
    master_path = out_root / "summary" / "master_summary.csv"
    master = pd.read_csv(master_path) if master_path.exists() else pd.DataFrame()
    return master, results


# --------------------------------------------------------------------------- #
# Pure multi-angle preview / batch logic (no Streamlit calls -> headlessly    #
# testable). Task 3 wires these into a new UI section; nothing here renders.  #
# --------------------------------------------------------------------------- #

#: Default µm/px offered for the manual multi-angle scale field -- the C1
#: "Scaling/Items" value resolved and used throughout study 03 (labbook 03).
_MA_UM_PER_PX_DEFAULT = 0.388924


def _profile_from_mr(mr) -> dict[str, np.ndarray]:
    """Build one angle's ``{"x", "w"}`` profile for ``xsection.build_part_stack``.

    Mirrors ``write_measurement``'s profile CSV contract (``x_px`` /
    ``diameter_px_smooth`` / ``valid`` -- the only three columns
    ``run_xsection`` reads back) directly from an in-memory ``MeasureResult``,
    without ever touching ``mr.rgb`` (it is ``None`` on the cached folder
    path -- see ``_cached_compute_path``). ``x`` is ABSOLUTE image pixels
    (``bnd.x0 .. bnd.x1`` inclusive, never relative to the span); ``valid``
    and ``diameter_smooth`` are FULL-image-width arrays on ``mr.res``, so both
    must be sliced by the same ``span`` before combining, exactly like
    ``_render_group``/``_group_mean_um`` already do (their own ``span =
    slice(bnd.x0, bnd.x1 + 1)`` + masked-profile blocks, below in this file).
    """
    bnd, res = mr.bnd, mr.res
    x = np.arange(bnd.x0, bnd.x1 + 1, dtype=float)
    span = slice(bnd.x0, bnd.x1 + 1)
    w = np.where(res.valid[span], res.diameter_smooth[span], np.nan).astype(float)
    return {"x": x, "w": w}


def _profiles_from_results(mrs: dict[int, object], cfg: CONFIG
                           ) -> tuple[dict[int, dict], dict[int, str]]:
    """Split per-angle ``MeasureResult``s into included profiles + drop reasons.

    ``mrs`` maps angle (1..6) -> ``MeasureResult``. Applies the exact same QC
    policy ``_render_group`` uses (``anomaly.exclusion_reason`` on
    ``band_mismatch``/``coverage``/``anomaly.flags``), so a GUI multi-angle
    preview agrees with the CLI's ``run_xsection`` (which gained the same
    ``--anomaly-exclude`` flag in Task 1). Returns ``(profiles, excluded)``:
    ``profiles`` maps angle -> ``_profile_from_mr(mr)`` for angles that passed
    QC; ``excluded`` maps every dropped angle to a user-facing
    (``_friendly_reason``) string.
    """
    profiles: dict[int, dict] = {}
    excluded: dict[int, str] = {}
    for angle, mr in mrs.items():
        res = mr.res
        reason = exclusion_reason(res.band_mismatch, res.coverage,
                                  res.anomaly.flags, cfg)
        if reason is not None:
            excluded[angle] = _friendly_reason(reason)
            continue
        profiles[angle] = _profile_from_mr(mr)
    return profiles, excluded


def _ma_circular_phi_median(phi_deg: np.ndarray) -> float:
    """Orientation average with period 180 deg (vector mean of 2*phi).

    Local copy of ``run_xsection._circular_phi_median`` (private there, not
    part of this task's import list) so this module stays independent of
    ``run_xsection``'s internals -- only its public ``main`` is imported.
    """
    ph = np.radians(np.asarray(phi_deg, dtype=float))
    ph = ph[np.isfinite(ph)]
    if ph.size == 0:
        return float("nan")
    return float(np.degrees(np.arctan2(np.mean(np.sin(2 * ph)),
                                       np.mean(np.cos(2 * ph))) / 2) % 180.0)


@dataclass
class MultiAnglePreview:
    """Headless multi-angle cross-section preview for one (fiber, part).

    All geometry stays in pixels (no micron conversion -- the caller applies
    ``um_per_px``, mirroring ``xsection.py``'s split from ``run_xsection.py``).
    """

    stack: PartStack | None       # None when profiles was empty
    fit: XsecFit | None           # None when profiles was empty
    present: set                  # angles (1..6) present in the input profiles
    n_directions: int             # distinct projection directions among present
    fittable: bool                # n_directions == 3 (necessary, not sufficient
                                  # -- a fittable part can still have 0 valid
                                  # columns after the per-column 3-direction rule)
    med: dict                     # summary stats (empty dict when profiles empty)
                                  # keys: a/b/ratio/phi/area/area_err/rms/
                                  # dw_frac/valid_frac/n_uncertain/n_saturated/
                                  # n_links


def multiangle_preview(profiles: dict[int, dict], cfg: CONFIG) -> MultiAnglePreview:
    """Pure multi-angle cross-section preview from already-QC'd profiles.

    Mirrors ``run_xsection.py``'s per-part pipeline: ``build_part_stack`` ->
    ``fit_ellipse_projections`` -> ``split_half_area``/``pair_differences``,
    including its ``w_dir`` nanmean denominator and RuntimeWarning
    suppression for ``dw_frac``
    (run_xsection.py:360-366). ``profiles`` maps angle (1..6) -> ``{"x", "w"}``
    absolute-px profiles, e.g. the ``included`` half of
    ``_profiles_from_results``'s return. Empty ``profiles`` short-circuits to
    ``stack=None``/``fit=None``/``med={}`` (``build_part_stack`` raises on an
    empty dict, so this case must never reach it).

    ``present``/``n_directions``/``fittable`` are computed even when empty
    (all empty/zero/False) so the caller can render "0 of 6 angles" copy
    without a None-check. ``fiber``/``part`` are fixed placeholders (0, 1) --
    this is a single-part preview, not tied to any real (fiber, part) label.
    """
    present = set(profiles)
    n_directions = len({(a - 1) % 3 for a in present})
    fittable = n_directions == 3
    if not profiles:
        return MultiAnglePreview(stack=None, fit=None, present=present,
                                 n_directions=n_directions, fittable=fittable,
                                 med={})

    stack = build_part_stack(0, 1, profiles, cfg)
    fit = fit_ellipse_projections(stack.W)
    A_lo, A_hi = split_half_area(stack.W)
    pair_dw = pair_differences(stack.W)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        with np.errstate(invalid="ignore"):
            area_err_px2 = np.abs(A_lo - A_hi) / 2.0
            w_dir = np.nanmean(np.stack([stack.W[:3], stack.W[3:]]), axis=0)
            dw_frac = np.abs(pair_dw) / w_dir

    valid = fit.valid
    # n_uncertain counts only PRESENT angles marked uncertain: build_part_stack
    # also marks every MISSING angle uncertain=True (shifts.append at
    # xsection.py:390-392), so a naive `sum(uncertain)` would double-report
    # angles that are already absent from `present` -- mirrors the CLI's
    # `n_unc = sum(1 for s in st.shifts if s["present"] and s["uncertain"])`
    # at run_xsection.py:367-368.
    n_uncertain = sum(1 for s in stack.shifts if s["present"] and s["uncertain"])
    n_saturated = sum(1 for s in stack.shifts if s.get("saturated"))
    # n_links mirrors run_xsection.py:370/392 (`sum(present) - 1`, floored at
    # 0): the number of cross-angle correlation LINKS actually estimated, i.e.
    # one fewer than the count of present angles (the reference angle has no
    # link to itself). len(present) is exactly that present-count.
    n_links = max(len(present) - 1, 0)

    Wf = np.isfinite(stack.W)
    n_measurable = int(((Wf[0] | Wf[3]) & (Wf[1] | Wf[4]) & (Wf[2] | Wf[5])).sum())
    n_pos = int(valid.sum())
    valid_frac = n_pos / n_measurable if n_measurable else 0.0

    dwf = dw_frac[np.isfinite(dw_frac)]

    # nanmedian/nanmean below can legitimately see an all-NaN slice (e.g. a
    # missing angle removes a whole direction from one split-half fit, so
    # area_err_px2 is NaN at every column even though the full 6-angle fit
    # stays valid there) -- suppress numpy's advisory RuntimeWarning for that
    # case, same as xsection.py:165-167 / run_xsection.py:361-366.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        med = {
            "a": float(np.nanmedian(fit.a[valid])) if valid.any() else float("nan"),
            "b": float(np.nanmedian(fit.b[valid])) if valid.any() else float("nan"),
            "ratio": float(np.nanmedian((fit.a / fit.b)[valid])) if valid.any() else float("nan"),
            "phi": _ma_circular_phi_median(fit.phi_deg),
            "area": float(np.nanmedian(fit.area[valid])) if valid.any() else float("nan"),
            "area_err": float(np.nanmedian(area_err_px2[valid])) if valid.any() else float("nan"),
            "rms": float(np.nanmedian(fit.rms_resid[valid])) if valid.any() else float("nan"),
            "dw_frac": float(np.median(dwf)) if dwf.size else float("nan"),
        }
    med.update({
        "valid_frac": valid_frac,
        "n_uncertain": n_uncertain,
        "n_saturated": n_saturated,
        "n_links": n_links,
    })

    return MultiAnglePreview(stack=stack, fit=fit, present=present,
                             n_directions=n_directions, fittable=fittable,
                             med=med)


def run_multiangle_batch(
    folder: str | Path,
    condition: str | None,
    out_root: str | Path,
    cfg: CONFIG,
    um_per_px: float,
    jobs: int,
    progress_cb: Callable[[float, dict], None] | None = None,
) -> tuple[pd.DataFrame, list[dict], int]:
    """Measure a whole multi-angle folder, then run the ``run_xsection`` stage.

    ``discover_multiangle(folder, condition)`` finds every ``(fiber, part)``
    -> ``{angle: path}`` group (scale-bar twins and unparseable names are
    already skipped there); the paths are flattened and measured in-process
    via ``run_batch(..., aggregate=False)`` (stage-2 group registration is
    meaningless for angle data -- ``run_xsection`` is its own aggregation
    stage). Then runs ``run_xsection.main`` with a NUMERIC ``--scale-source``
    (``um_per_px``, e.g. ``0.388924`` -- bypasses the XML sidecars entirely,
    per Task 1) and ``--anomaly-exclude`` mirrored from ``cfg.anomaly_exclude``.
    ``--condition`` is ALWAYS passed explicitly: ``run_xsection``'s
    ``by_fiber`` keys only on fiber number (run_xsection.py:261/379), so a
    folder mixing conditions would silently merge them without this filter.

    ``condition`` must be a concrete string -- raises ``ValueError`` on
    ``None`` rather than silently degrading (e.g. to "measure everything");
    the caller (a later task) is responsible for resolving a single detected
    condition, or asking the user, before calling this.

    Returns ``(summary_df, results, rc)``: ``summary_df`` is
    ``summary/xsection_summary.csv`` read back (empty DataFrame if
    ``run_xsection`` failed or wrote nothing), ``results`` is ``run_batch``'s
    per-image result dicts, and ``rc`` is ``run_xsection.main``'s return code
    (0 on success).
    """
    if condition is None:
        raise ValueError(
            "run_multiangle_batch requires a concrete condition string "
            "(discover_multiangle can find more than one; resolve to a "
            "single condition before calling this, never pass None)."
        )
    folder = Path(folder)
    out_root = Path(out_root)
    groups = discover_multiangle(folder, condition)
    image_paths = [p for angles in groups.values() for p in angles.values()]

    _master, results = run_batch(image_paths, out_root, cfg, jobs, progress_cb,
                                 aggregate=False)

    rc = run_xsection.main([
        "--out", str(out_root), "--data-root", str(folder),
        "--condition", condition, "--scale-source", str(um_per_px),
    ] + (["--anomaly-exclude"] if cfg.anomaly_exclude else []))

    summary_path = out_root / "summary" / "xsection_summary.csv"
    summary_df = pd.read_csv(summary_path) if summary_path.exists() else pd.DataFrame()
    return summary_df, results, rc


def _group_mean_um(reps: list[dict], cfg: CONFIG) -> float | None:
    """Registered mean diameter (µm) for one group's reps, or ``None`` if it
    cannot be measured. Same profile-building + registration as ``_render_group``,
    so the matrix matches the on-screen group panel."""
    profiles = []
    for rep in reps:
        mr = rep["mr"]
        res, bnd = mr.res, mr.bnd
        if exclusion_reason(res.band_mismatch, res.coverage,
                            res.anomaly.flags, cfg) is not None:
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
_ACCENT = "#660099"
_REP_CYCLE = ("#A78BFA", "#F9A8D4", "#93C5FD", "#FCD34D")
_MUTED = "#64748B"
_GRID = "#E2E8F0"


def _styled_fig(figsize: tuple[float, float] = (9, 3)):
    """Build a Figure/Axes pair sharing the app's muted, transparent plot style.

    Top/right spines are dropped; the remaining spines, ticks and axis labels
    use the muted slate (``_MUTED``); the grid is a faint slate (``_GRID``);
    the figure and axes backgrounds are transparent so the plot sits directly
    on the surrounding section card instead of a white rectangle.
    """
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_alpha(0.0)
    ax.set_facecolor("none")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(_MUTED)
    ax.spines["bottom"].set_color(_MUTED)
    ax.tick_params(colors=_MUTED)
    ax.xaxis.label.set_color(_MUTED)
    ax.yaxis.label.set_color(_MUTED)
    ax.grid(color=_GRID, lw=0.8, alpha=0.3)
    return fig, ax


def _profile_fig(mr, rgb, cfg: CONFIG, um_per_px: float | None = None):
    """Per-replicate diameter-vs-position figure (raw points + smoothed line, µm).

    By default both series are scaled by ``cfg.ppu`` (pixels per micron, so µm
    = px / ppu). ``um_per_px`` overrides that with the reciprocal calibration
    used by the multi-angle mode (µm = px * um_per_px); it must replace BOTH
    series, since ``mr.diameter_um`` is itself already ppu-derived and mixing
    the two would draw the raw points and the smooth line on different scales.
    """
    bnd, res = mr.bnd, mr.res
    span = slice(bnd.x0, bnd.x1 + 1)
    x = np.arange(rgb.shape[1])[span]
    if um_per_px is None:
        raw_um = mr.diameter_um[span]
        sm_um = res.diameter_smooth[span] / cfg.ppu
    else:
        raw_um = np.where(res.valid[span],
                          res.diameter_raw[span] * um_per_px, np.nan)
        sm_um = res.diameter_smooth[span] * um_per_px
    fig, ax = _styled_fig()
    ax.plot(x, raw_um, ".", ms=2, alpha=0.4, color="#94A3B8", label="raw")
    ax.plot(x, sm_um, "-", lw=1.3, color=_ACCENT, label="smooth")
    ax.set_xlabel("x position (px)")
    ax.set_ylabel("diameter (µm)")
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
    fig, ax = _styled_fig()
    for i, (rep, rx, ry) in enumerate(rep_curves):
        ax.plot(rx, ry, "-", lw=0.8, alpha=0.45, color=_REP_CYCLE[i % 4],
                label=f"rep {rep}")
    x = table["x_aligned_px"]
    mean = table["mean_um"]
    std = table["std_um"]
    ax.plot(x, mean, "-", lw=2.0, color=_ACCENT, zorder=5,
            label="mean of replicates")
    band = np.where(np.isfinite(std), std, 0.0)
    ax.fill_between(x, mean - band, mean + band, alpha=0.25, color=_ACCENT,
                    label="±std across replicates")
    ax.set_xlabel("aligned x position (px)")
    ax.set_ylabel("diameter (µm)")
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
    fig, ax = _styled_fig()
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
        ax.plot(x[:brk + 1], y[:brk + 1], "-", lw=1.6, color=_ACCENT,
                label="load")
        if brk + 1 < x.size:
            ax.plot(x[brk:], y[brk:], "-", lw=0.8, alpha=0.35, color=_ACCENT)
        ax.scatter([x[brk]], [y[brk]], s=30, color="#DC2626", zorder=6,
                   label="break point")
        _focus(x[:brk + 1], y[:brk + 1], brk + 1 < x.size)
        ax.set_xlabel("displacement (mm)")
        ax.set_ylabel("load (N)")
        ax.legend(loc="best", fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        return fig

    x = res.strain * 100.0          # %
    y = res.stress_pa / 1e6         # MPa
    ax.plot(x[:brk + 1], y[:brk + 1], "-", lw=1.6, color=_ACCENT,
            label="stress–strain")
    if brk + 1 < x.size:
        ax.plot(x[brk:], y[brk:], "-", lw=0.8, alpha=0.35, color=_ACCENT,
                label="post-break (unused)")

    # toughness = area under the curve up to fracture
    ax.fill_between(x[:brk + 1], 0.0, y[:brk + 1], alpha=0.15,
                    color=_ACCENT, label="toughness (area)")

    # steepest initial slope = Young's modulus, drawn over its fit segment
    fit = res.modulus_fit or {}
    if fit and np.isfinite(fit.get("slope", np.nan)):
        s_lo, s_hi = fit["strain_lo"], fit["strain_hi"]
        s_line = np.array([s_lo, s_hi], dtype=float)
        stress_fit = fit["slope"] * s_line + fit["intercept"]
        ax.plot(s_line * 100.0, stress_fit / 1e6, "--", lw=1.4,
                color="#D97706",
                label=f"E = {_fmt(res.youngs_modulus_pa / 1e9, 'GPa')}")

    # break point
    ax.scatter([x[brk]], [y[brk]], s=30, color="#DC2626", zorder=6,
               label="break point")

    _focus(x[:brk + 1], y[:brk + 1], brk + 1 < x.size)
    ax.set_xlabel("strain (%)")
    ax.set_ylabel("stress (MPa)")
    ax.legend(loc="best", fontsize=7)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


# six-angle colour cycle for the multi-angle plots: one hue per angle, with
# each 180-degree pair (a1/a4, a2/a5, a3/a6) sharing a hue family so a pair
# that disagrees is visible as two shades of the same colour
_ANGLE_CYCLE = ("#660099", "#0EA5E9", "#059669",
                "#C084FC", "#7DD3FC", "#6EE7B7")


def _xsec_stack_fig(stack, fit, um_per_px: float):
    """The six aligned width profiles (µm) of one part, on the shared x grid.

    Columns where no ellipse could be fitted (fewer than three projection
    directions measurable after alignment) are shaded, so a metric computed
    over "valid" columns can be read against how much of the part that was.
    """
    fig, ax = _styled_fig()
    x = stack.x
    for k in range(6):
        w = stack.W[k]
        if not np.isfinite(w).any():
            continue
        ax.plot(x, w * um_per_px, "-", lw=1.0, color=_ANGLE_CYCLE[k],
                label=f"a{k + 1}")
    invalid = ~fit.valid
    if invalid.any():
        ax.fill_between(x, 0, 1, where=invalid,
                        transform=ax.get_xaxis_transform(), step="mid",
                        color="#94A3B8", alpha=0.15, lw=0,
                        label="no 3-direction fit")
    ax.set_xlabel("aligned x position (px, reference-angle frame)")
    ax.set_ylabel("width (µm)")
    ax.legend(loc="best", fontsize=7, ncols=4)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def _xsec_ellipse_fig(med: dict, um_per_px: float):
    """The median cross-section drawn to scale: ellipse, equal-area circle, φ.

    ``med`` is ``MultiAnglePreview.med`` (px); everything is converted here.
    The dashed circle has the same area as the ellipse, so the gap between
    the two curves *is* the anisotropy the axis ratio reports. Returns
    ``None`` when the median axes are not finite (nothing to draw).
    """
    a = float(med.get("a", np.nan)) * um_per_px
    b = float(med.get("b", np.nan)) * um_per_px
    if not (np.isfinite(a) and np.isfinite(b)):
        return None
    phi_deg = float(med.get("phi", np.nan))
    ph = np.radians(phi_deg if np.isfinite(phi_deg) else 0.0)

    fig, ax = _styled_fig((3.8, 3.8))
    t = np.linspace(0.0, 2 * np.pi, 361)
    ex = a * np.cos(t) * np.cos(ph) - b * np.sin(t) * np.sin(ph)
    ey = a * np.cos(t) * np.sin(ph) + b * np.sin(t) * np.cos(ph)
    ax.fill(ex, ey, color=_ACCENT, alpha=0.18)
    ax.plot(ex, ey, "-", lw=1.6, color=_ACCENT, label="median ellipse")
    r = float(np.sqrt(a * b))          # equal-area circle
    ax.plot(r * np.cos(t), r * np.sin(t), "--", lw=1.0, color=_MUTED,
            label="equal-area circle")
    if np.isfinite(phi_deg):
        ax.annotate("", xy=(a * np.cos(ph), a * np.sin(ph)), xytext=(0.0, 0.0),
                    arrowprops=dict(arrowstyle="->", color="#D97706", lw=1.4))
        ax.text(0.55 * a * np.cos(ph), 0.55 * a * np.sin(ph),
                f" φ {phi_deg:.0f}°", color="#D97706", fontsize=8)
    lim = 1.15 * max(a, r)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("µm")
    ax.set_ylabel("µm")
    ax.legend(loc="upper right", fontsize=7)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def _xsec_area_fig(stack, fit, um_per_px: float, show_err: bool):
    """Per-column ellipse area along the part, with its bounds (mirrors the
    CLI's ``run_xsection._plot_part`` lower panel).

    The split-half ± band is recomputed here from the stack
    (``MultiAnglePreview`` carries only the medians, not the per-column
    arrays); it is a pure function over ``stack.W``, so this cannot disagree
    with the preview's numbers. ``show_err=False`` drops the band, which is
    all-NaN whenever an angle is missing (each split half needs all three
    directions on its own).
    """
    A_lo, A_hi = split_half_area(stack.W)
    um2 = um_per_px ** 2
    area = np.where(fit.valid, fit.area, np.nan) * um2

    fig, ax = _styled_fig()
    ax.plot(stack.x, area, "-", lw=1.4, color=_ACCENT, label="ellipse area")
    if show_err:
        err = np.abs(A_lo - A_hi) / 2.0 * um2
        band = np.where(np.isfinite(err), err, 0.0)
        ax.fill_between(stack.x, area - band, area + band, alpha=0.22,
                        color=_ACCENT, label="±split-half err")
    ax.set_xlabel("aligned x position (px, reference-angle frame)")
    ax.set_ylabel("area (µm²)")
    ax.legend(loc="best", fontsize=7)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------- #
# Sidebar: parameter form                                                     #
# --------------------------------------------------------------------------- #
#: the two analysis modes (radio option strings, also the session values of
#: the "analysis_mode" key)
_MODE_REPS = "Replicates"
_MODE_MULTIANGLE = "Multi-angle cross-section"


def _analysis_mode() -> str:
    """Sidebar analysis-mode radio; couples multi-angle to bright mode ONCE.

    Called at the top of ``main()``, before the ``cfg_items`` snapshot, so a
    mode switch is reflected in the same run (no extra rerun). The C1
    multi-angle set is bright-on-dark, so ENTERING multi-angle mode re-applies
    the calibrated bright defaults; that must happen only on the transition,
    not on every rerun, or the user's own tuning would be stomped on each
    interaction. Since the radio carries ``key="analysis_mode"``, its session
    value is already the NEW value by the time this returns, so the transition
    is detected against a private shadow key instead. Leaving the mode does not
    auto-revert (the user may well want to keep bright).
    """
    st.sidebar.markdown('<div class="fcv-side-label">Analysis</div>',
                        unsafe_allow_html=True)
    mode = st.sidebar.radio(
        "analysis mode", [_MODE_REPS, _MODE_MULTIANGLE],
        label_visibility="collapsed", key="analysis_mode",
        help="Replicates: several photos of the same fibre, averaged into one "
             "diameter profile. Multi-angle cross-section: the six rotated "
             "views (a1..a6) of one fibre part, combined into a per-position "
             "ellipse cross-section.")
    prev = st.session_state.get("_prev_analysis_mode")
    if prev != mode and mode == _MODE_MULTIANGLE:
        st.session_state.cfg_dict = {**st.session_state.cfg_dict,
                                     "feature_mode": "bright",
                                     **BRIGHT_DEFAULTS}
        st.session_state.form_version += 1
    st.session_state["_prev_analysis_mode"] = mode
    return mode


def _param_form(multiangle: bool = False) -> None:
    """Render the 3-knob parameter form; updates session_state on Apply/Reset.

    ``multiangle`` makes Reset mode-aware: the multi-angle images are
    bright-on-dark, so resetting there must land on the bright baseline
    (``feature_mode`` + ``BRIGHT_DEFAULTS``) rather than silently dropping the
    session back to the desat defaults, which would leave every angle
    mis-detected until the user noticed the mode selector had flipped.
    """
    applied = st.session_state.cfg_dict
    ver = st.session_state.form_version

    st.sidebar.markdown('<div class="fcv-side-label">Detection</div>',
                        unsafe_allow_html=True)
    with st.sidebar.form("params", clear_on_submit=False):
        st.markdown("**Parameters** — edit, then click **Apply** to re-render.")
        new_vals: dict = {}

        mode_options = list(_MODE_LABELS)
        cur_mode = applied.get("feature_mode", "desat")
        new_vals["feature_mode"] = st.selectbox(
            "Image mode",
            mode_options,
            index=(mode_options.index(cur_mode)
                   if cur_mode in mode_options else 0),
            format_func=lambda m: _MODE_LABELS[m],
            help="Which feature map the detector runs on. desat: HSV "
                 "desaturation z-map (MasP2 look — pale fibre on a saturated "
                 "background). bright: median-RGB brightness z-map (C1 look — "
                 "bright fibre on a dark background). Applying a mode switch "
                 "re-applies that mode's calibrated defaults for edge_frac "
                 "and k_band (bright: 0.30 / 6.0, desat: 0.65 / 4.0).",
            key=f"p_feature_mode_v{ver}")

        def _render_spec(spec: tuple) -> None:
            name, label, kind, help_txt, step, lo, hi, fmt = spec
            key = f"p_{name}_v{ver}"
            cur = applied[name]
            if kind == "slider":
                new_vals[name] = st.slider(
                    label, min_value=float(lo), max_value=float(hi),
                    value=float(cur), step=float(step), help=help_txt, key=key)
            elif kind == "int":
                new_vals[name] = int(st.number_input(
                    label, min_value=int(lo), max_value=int(hi),
                    value=int(cur), step=int(step), help=help_txt, key=key))
            elif kind == "bool":
                new_vals[name] = bool(st.checkbox(
                    label, value=bool(cur), help=help_txt, key=key))
            else:  # float
                kwargs = dict(min_value=float(lo), max_value=float(hi),
                              value=float(cur), step=float(step), help=help_txt, key=key)
                if fmt:
                    kwargs["format"] = fmt
                new_vals[name] = float(st.number_input(label, **kwargs))

        # anomaly knobs (the tail of PARAM_SPECS) fold into their own expander
        split = next((i for i, s in enumerate(PARAM_SPECS)
                      if s[0] == "jump_thresh_px"), len(PARAM_SPECS))
        for spec in PARAM_SPECS[:split]:
            if spec[0] == "ppu":
                st.markdown("**Calibration**")
            _render_spec(spec)
        if split < len(PARAM_SPECS):
            with st.expander("Anomaly flags", expanded=False):
                for spec in PARAM_SPECS[split:]:
                    _render_spec(spec)
        c1, c2 = st.columns(2)
        apply = c1.form_submit_button("Apply", type="primary", width="stretch")
        reset = c2.form_submit_button("Reset to defaults", width="stretch")

    if reset:
        base = DEFAULTS.as_dict()
        if multiangle:
            base.update({"feature_mode": "bright", **BRIGHT_DEFAULTS})
        st.session_state.cfg_dict = base
        st.session_state.form_version += 1
        st.rerun()
    if apply:
        # merge: the visible knobs override a full defaults dict, so hidden
        # fields always carry the validated values. Bright mode starts from
        # BRIGHT_DEFAULTS (study 03) so the hidden k_band rides along,
        # mirroring run_measure.build_config.
        mode = new_vals.get("feature_mode", "desat")
        base = DEFAULTS.as_dict()
        if mode == "bright":
            base.update(BRIGHT_DEFAULTS)
        if mode != applied.get("feature_mode", "desat"):
            # mode switch: the mode-coupled visible knob jumps to the new
            # mode's calibrated default rather than keeping the stale widget
            # value; the version bump refreshes the widgets to show it
            new_vals.pop("edge_frac", None)
            st.session_state.form_version += 1
        st.session_state.cfg_dict = {**base, **new_vals}
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
    with st.sidebar.expander("Tensile", expanded=False):
        source = st.radio("tensile source", ["Local folder", "Upload"],
                          horizontal=True, label_visibility="collapsed",
                          key="tensile_source")
        tmap: dict[str, object] = {}
        folder: str | None = None
        if source == "Local folder":
            folder = st.text_input(
                "Tensile data folder", value=st.session_state.get("tensile_folder", ""))
            st.session_state.tensile_folder = folder
            if folder and Path(folder).is_dir():
                try:
                    tmap = {g: Path(p)
                            for g, p in _cached_discover_tensile(folder).items()}
                except Exception as exc:  # noqa: BLE001
                    st.warning(f"Could not scan tensile folder: {exc}")
            elif folder:
                st.warning("Enter a valid tensile folder path.")
        else:
            folder_mode = st.checkbox(
                "📁 upload a whole folder", key="tensile_folder_mode",
                help="Make Browse open a folder chooser; non-tensile files are ignored.")
            ups = st.file_uploader(
                "Upload tensile files",
                type=None if folder_mode else ["csv", "xls", "xlsx"],
                accept_multiple_files=True, key="tensile_uploads")
            if folder_mode:
                _enable_folder_upload("Upload tensile files")
            if ups:
                tmap = discover_tensile_files(ups)
        if tmap:
            st.caption(f"{len(tmap)} tensile fibre(s) matched.")

        gauge_length_mm = float(st.number_input(
            "Gauge length L₀ (mm)", min_value=0.1, max_value=1000.0,
            value=float(DEFAULTS.gauge_length_mm),
            step=1.0, format="%.1f",
            help="Grip separation; strain = displacement / L₀. The tester records "
                 "only displacement, so this sets the strain scale."))
        modulus_window = float(st.number_input(
            "Modulus window (fraction)", min_value=0.02, max_value=0.5,
            value=float(DEFAULTS.modulus_window), step=0.01, format="%.2f",
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

    A collapsed sidebar expander re-mounts the uploader's DOM on re-expand, which
    drops the one-shot attribute. A ``MutationObserver`` on the parent document
    re-applies it whenever a matching uploader (re)appears, so it survives both
    client-side expander toggles (no Python rerun) and Streamlit reruns. The
    observer is created once per label (guarded by a flag on ``window.parent``)
    and cleaned up on iframe teardown (``pagehide`` — ``unload`` is deprecated
    and Chrome may never fire it, leaving the flag stuck) so a later
    re-creation of this component (e.g. the folder-mode checkbox toggled off
    then on) gets a fresh observer instead of being silently blocked by a
    stale flag.
    """
    js = """
    <script>
    (function () {
      const WANT = %s;
      const doc = window.parent.document;
      const FLAG = "__fcv_folder_observer_" + WANT;

      function apply() {
        doc.querySelectorAll('[data-testid="stFileUploader"]').forEach(u => {
          if ((u.innerText || "").includes(WANT)) {
            const inp = u.querySelector('input[type="file"]');
            if (inp && !inp.hasAttribute("webkitdirectory")) {
              inp.setAttribute("webkitdirectory", "");
              inp.setAttribute("directory", "");
            }
          }
        });
      }

      apply(); setTimeout(apply, 150); setTimeout(apply, 500);

      if (!window.parent[FLAG]) {
        window.parent[FLAG] = true;
        const observer = new MutationObserver(apply);
        observer.observe(doc.body, {childList: true, subtree: true});
        window.addEventListener("pagehide", function () {
          observer.disconnect();
          window.parent[FLAG] = false;
        });
      }
    })();
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
    st.sidebar.markdown('<div class="fcv-side-label">Data</div>',
                        unsafe_allow_html=True)
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


def _ma_classify(items: list, key) -> tuple[list[tuple], int, int]:
    """Split a file listing into multi-angle images, scale-bar twins and rest.

    Returns ``(parsed, n_scalebar, n_other)`` where ``parsed`` is a list of
    ``(MultiAngleKey, item)`` for plain (non-scale-bar) images whose stem
    matches the strict ``<cond>_<fiber>_a<angle>_part<part>`` convention.
    Angles outside 1..6 count as "other": the whole cross-section stage is
    built on the six nominal rotations, so an ``a7`` would be silently ignored
    downstream rather than measured. Unlike ``discover_multiangle`` (which
    accepts any parseable angle and lets a missing/extra angle simply be an
    absent or unused dict key), this function actively rejects angles outside
    1..6 so the UI can report them as skipped rather than silently drop them.
    """
    parsed: list[tuple] = []
    n_scalebar = n_other = 0
    for it in items:
        try:
            k = parse_multiangle_name(key(it))
        except ValueError:
            n_other += 1
            continue
        if k.scalebar:
            n_scalebar += 1
            continue
        if not 1 <= k.angle <= 6:
            n_other += 1
            continue
        parsed.append((k, it))
    return parsed, n_scalebar, n_other


def _load_multiangle(cfg_items: tuple) -> dict:
    """Resolve the multi-angle data-source controls to one (fiber, part) set.

    The multi-angle twin of ``_load_reps``: same "Image folder" text input
    (so it stays ``text_input[0]`` and shares ``st.session_state.folder``),
    but the images are grouped by ``parse_multiangle_name`` instead of
    ``parse_name``, and the selection narrows to a single condition / fibre /
    part whose ≤6 angle images are then computed with the shared caches.

    The upload source deliberately offers NO "upload a whole folder" checkbox
    (unlike replicate mode): a real multi-angle condition directory is tens of
    GB of TIFFs and pointing a folder chooser at it would push all of it
    through the browser websocket. Hand-picked files only.

    Returns ``{"angles": {angle: {"name", "rgb", "mr"}}, "condition", "fiber",
    "part", "folder", "um_per_px", "n_images", "n_scalebar", "n_other"}``;
    ``folder`` is None for the upload source (the batch card needs a real
    directory), and ``condition``/``fiber``/``part`` are None when nothing
    could be loaded.
    """
    st.sidebar.markdown('<div class="fcv-side-label">Data</div>',
                        unsafe_allow_html=True)
    source = st.sidebar.radio("multi-angle source", ["Local folder", "Upload"],
                              horizontal=True, label_visibility="collapsed")
    out: dict = {"angles": {}, "condition": None, "fiber": None, "part": None,
                 "folder": None, "n_images": 0, "n_scalebar": 0, "n_other": 0}

    parsed: list[tuple] = []
    is_upload = source != "Local folder"
    if not is_upload:
        folder = st.sidebar.text_input(
            "Image folder", value=st.session_state.get("folder", ""))
        st.session_state.folder = folder
        if folder and Path(folder).is_dir():
            out["folder"] = folder
            parsed, n_sb, n_other = _ma_classify(discover_images(folder),
                                                 key=lambda p: p.name)
            out.update(n_scalebar=n_sb, n_other=n_other)
        elif folder:
            st.sidebar.warning("Enter a valid local folder path.")
    else:
        uploads = st.sidebar.file_uploader(
            "Upload multi-angle images",
            type=["jpg", "jpeg", "png", "tif", "tiff", "bmp"],
            accept_multiple_files=True, key="ma_uploads")
        parsed, n_sb, n_other = _ma_classify(list(uploads or []),
                                             key=lambda u: u.name)
        out.update(n_scalebar=n_sb, n_other=n_other)

    out["n_images"] = len(parsed)
    if parsed:
        st.sidebar.caption(
            f"{out['n_images']} multi-angle image(s) · {out['n_scalebar']} "
            f"scale-bar twin(s) skipped · {out['n_other']} other file(s) ignored")
        conditions = sorted({k.condition for k, _ in parsed})
        condition = (st.sidebar.selectbox("Condition", conditions)
                     if len(conditions) > 1 else conditions[0])
        in_cond = [(k, it) for k, it in parsed if k.condition == condition]
        fibers = sorted({k.fiber for k, _ in in_cond})
        fiber = st.sidebar.selectbox("Fibre", fibers,
                                     format_func=lambda f: f"{f:02d}")
        parts = sorted({k.part for k, _ in in_cond if k.fiber == fiber})
        part = st.sidebar.selectbox("Part", parts)
        chosen = {k.angle: it for k, it in in_cond
                  if k.fiber == fiber and k.part == part}
        out.update(condition=condition, fiber=fiber, part=part)

        angles: dict[int, dict] = {}
        for angle, item in sorted(chosen.items()):
            if is_upload:
                stem = Path(item.name).stem
                mr, rgb = _cached_compute_upload(stem, item.getvalue(), cfg_items)
            else:
                mtime = Path(item).stat().st_mtime
                stem = Path(item).stem
                mr = _cached_compute_path(str(item), mtime, cfg_items)
                rgb = _cached_rgb_from_path(str(item), mtime)
            angles[angle] = {"name": stem, "rgb": rgb, "mr": mr}
        out["angles"] = angles
    elif out["folder"] is not None or out["n_other"] or out["n_scalebar"]:
        # something WAS offered and none of it is usable -- an empty uploader
        # (nothing offered yet) must stay quiet
        st.sidebar.warning(
            "No multi-angle filenames found (expected "
            "`<cond>_<fibre>_a<angle>_part<part>`, e.g. C1_01_a1_part1.tiff). "
            "Switch to Replicates mode for ordinary replicate images.")

    # scale: manual and deliberately OUTSIDE the parameter form -- it is not a
    # CONFIG field, no compute cache key includes it, and every µm number in
    # this mode is derived from it at render time, so changing it must not
    # trigger a recompute (an Apply would).
    #
    # The widget's own value is seeded from (and mirrored back into) the plain
    # session key "ma_um_per_px_value": Streamlit DROPS a keyed widget's
    # session entry on any run where the widget is not rendered, so a trip
    # through replicate mode would otherwise silently revert a hand-calibrated
    # scale to the C1 default -- and this number drives every µm figure here
    # plus the batch's --scale-source. Same mirroring trick as
    # ``st.session_state.folder`` above; the mirror key is deliberately a
    # DIFFERENT name from the widget key, since assigning to the widget's own
    # key while also passing ``value=`` is what Streamlit warns about.
    out["um_per_px"] = float(st.sidebar.number_input(
        "Scale (µm per pixel)", min_value=0.000001, max_value=100.0,
        value=float(st.session_state.get("ma_um_per_px_value",
                                         _MA_UM_PER_PX_DEFAULT)),
        step=0.001, format="%.6f", key="ma_um_per_px",
        help="Microns per pixel for the multi-angle images. Default "
             f"{_MA_UM_PER_PX_DEFAULT} is the C1 objective's resolved "
             "Scaling/Items value (study 03). The sidebar's ppu calibration is "
             "NOT used in this mode — every µm number here comes from this "
             "field, and the batch passes it to run_xsection as "
             "--scale-source."))
    st.session_state.ma_um_per_px_value = out["um_per_px"]
    return out


# --------------------------------------------------------------------------- #
# Theme: CSS, slim header, jump menu                                          #
# --------------------------------------------------------------------------- #
# Most colours/radii live in .streamlit/config.toml (Streamlit's own theme
# engine); this constant covers what that engine cannot reach: the custom
# header/chip markup, section-card framing around the existing
# st.container(key=...) pattern, the fixed jump nav, and a chrome backstop.
_CSS = """
<style>
/* ---- chrome backstop: config.toml's toolbarMode="minimal" hides most of
   this already; kept as a belt-and-braces override ---- */
#MainMenu, footer, [data-testid="stDecoration"],
[data-testid="stAppDeployButton"] { display: none !important; }

/* ---- slim header + state chip (see _render_header) ---- */
.fcv-header {
    display: flex;
    align-items: baseline;
    flex-wrap: wrap;
    gap: 0.6rem 1rem;
    padding-bottom: 0.6rem;
    margin-bottom: 1rem;
    border-bottom: 1px solid #E2E8F0;
}
.fcv-header h1 {
    font-size: 1.4rem;
    font-weight: 700;
    color: #0F172A;
    margin: 0;
}
.fcv-header .fcv-sub {
    font-size: 0.85rem;
    color: #64748B;
}
.fcv-chip {
    margin-left: auto;
    padding: 0.2rem 0.7rem;
    border-radius: 999px;
    background: #F3E8FF;
    color: #660099;
    font-size: 0.8rem;
    font-weight: 600;
    white-space: nowrap;
}

/* ---- sidebar group captions (DATA / DETECTION), rendered via
   st.sidebar.markdown('<div class="fcv-side-label">...</div>'); the Tensile
   group uses a native st.sidebar.expander instead, so it needs no caption ---- */
.fcv-side-label {
    font-size: 0.72rem;
    font-weight: 700;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    color: #64748B;
    margin: 0.2rem 0 0.35rem;
}

/* ---- multi-angle coverage chips (_render_xsection): one per angle a1..a6.
   ok = in the fit, miss = no image for that angle, warn = QC-excluded or an
   uncertain/saturated alignment shift ---- */
.fcv-achips {
    display: flex;
    flex-wrap: wrap;
    gap: 0.4rem;
    margin: 0.1rem 0 0.9rem;
}
.fcv-achip {
    padding: 0.15rem 0.65rem;
    border-radius: 999px;
    border: 1px solid transparent;
    font-size: 0.78rem;
    font-weight: 600;
    white-space: nowrap;
}
.fcv-achip.ok   { background: #DCFCE7; color: #15803D; border-color: #BBF7D0; }
.fcv-achip.miss { background: #F1F5F9; color: #64748B; border-color: #E2E8F0; }
.fcv-achip.warn { background: #FEF3C7; color: #B45309; border-color: #FDE68A; }

/* ---- section cards: st.container(key="card_replicates"/"card_group"/
   "card_tensile"/"card_export", and the multi-angle mode's "card_angles"/
   "card_xsec"/"card_xsec_export") -> class st-key-card_<name>; matched by
   prefix so this one rule covers every card in both modes ---- */
[class*="st-key-card_"] {
    background: #FFFFFF;
    border: 1px solid #E2E8F0;
    border-radius: 0.75rem;
    padding: 1.1rem 1.4rem 1.4rem;
    margin-bottom: 1.2rem;
}

/* ---- tensile metrics: shrink the value font so figures like "161.00 mN"
   don't ellipsis in six narrow columns (moved from the inline <style> that
   used to sit next to st.container(key="tensile_metrics")) ---- */
.st-key-tensile_metrics [data-testid="stMetricValue"] { font-size: 1.1rem; }
.st-key-tensile_metrics [data-testid="stMetricLabel"] { font-size: 0.8rem; }

/* ---- metric cards: st.container(key="metrics_rep_<name>"/"metrics_group")
   -> class st-key-metrics_<name>; matched by prefix so this one rule covers
   both. Each real st.metric widget inside becomes a bordered white card;
   label is small/uppercase/grey, value is tabular-nums and contained
   (min-width:0 so it never grows its column). Overflow/white-space/text-
   overflow are deliberately left unset here (rather than nowrap+ellipsis)
   so the shared stMetricValue wrap rule below is the one in effect: a value
   like "62.43 µm" must always be fully readable, so it wraps to a second
   line instead of being clipped. */
[class*="st-key-metrics_"] [data-testid="stMetric"] {
    background: #FFFFFF;
    border: 1px solid #E2E8F0;
    border-radius: 0.6rem;
    padding: 0.6rem 0.8rem;
    min-width: 0;
    overflow: hidden;
}
[class*="st-key-metrics_"] [data-testid="stMetricLabel"] {
    font-size: 0.72rem;
    font-weight: 700;
    letter-spacing: 0.04em;
    text-transform: uppercase;
    color: #64748B;
}
[class*="st-key-metrics_"] [data-testid="stMetricValue"] {
    font-variant-numeric: tabular-nums;
    min-width: 0;
}
/* primary metric (first column of the row) reads as the headline number */
[class*="st-key-metrics_"] [data-testid="stColumn"]:first-of-type
    [data-testid="stMetricValue"] { color: #660099; }

/* ---- flags / registration status: the single metric wrapped in
   st.container(key="status_ok_*"/"status_warn_*") gets its value coloured
   green (good: no flags / registration ok) or amber (needs attention) ---- */
[class*="st-key-status_ok"] [data-testid="stMetricValue"] { color: #15803D; }
[class*="st-key-status_warn"] [data-testid="stMetricValue"] { color: #B45309; }

/* ---- wrap, don't clip: lets any metric value wrap to a second line instead
   of being silently ellipsised/clipped when its column is narrow (e.g.
   "62.43 µm" at 1280px). Streamlit renders the value text in a <p> inside
   stMetricValue and ships its own emotion-generated ".st-emotion-cache-* p"
   rule (specificity 0,1,1: one class + the p type) hard-coding
   nowrap/hidden/ellipsis on that <p> -- overflow/text-overflow/white-space
   are not inherited, so overriding stMetricValue (the div) alone never
   reaches the actual text node. The " p" selector here matches that
   specificity so ours (declared later) wins in the cascade. ---- */
[data-testid="stMetricValue"],
[data-testid="stMetricValue"] p {
    overflow: visible;
    white-space: normal;
    text-overflow: clip;
}

/* ---- jump menu: fixed right-hand nav, hidden below 1200px ---- */
.fcv-jump {
    position: fixed;
    top: 5.5rem;
    right: 1rem;
    z-index: 999;
    display: flex;
    flex-direction: column;
    gap: 0.3rem;
    padding: 0.6rem 0.9rem;
    background: #FFFFFF;
    border: 1px solid #E2E8F0;
    border-radius: 0.75rem;
    font-size: 0.8rem;
}
.fcv-jump a { color: #660099; text-decoration: none; }
.fcv-jump a:hover { text-decoration: underline; }
@media (max-width: 1200px) {
    .fcv-jump { display: none; }
}
/* ---- reserve the nav's column wherever it is shown (>=1200px) so the fixed
   .fcv-jump never sits over page content: the header chip's right edge at
   wide viewports, or a plot's corner at narrower ones. Padding goes on the
   block container (covers every child, incl. .fcv-header and st.pyplot
   figures) sized to clear the nav's own width (~1rem right offset + ~0.9rem
   padding each side + its longest link, "04 Export & batch") plus a gap ---- */
@media (min-width: 1200px) {
    [data-testid="stMainBlockContainer"] { padding-right: 11rem; }
}

/* anchored headings (card subheaders) land below the fixed header on jump */
h1, h2, h3 { scroll-margin-top: 4.5rem; }
</style>
"""


def _inject_css() -> None:
    """Emit ``_CSS`` once per run; called first thing in ``main()``."""
    st.markdown(_CSS, unsafe_allow_html=True)


def _render_header(group_label: str | None, n_reps: int, edge_z: float,
                   feature_mode: str = "desat", noun: str = "replicate",
                   label_prefix: str = "group") -> None:
    """Slim header + state chip, replacing ``st.title``/``st.caption``.

    Must run after the sidebar builders (``_load_reps`` etc.) since the chip
    needs ``n_reps``/``group_label``. The no-data branch is keyed on
    ``n_reps == 0``, not on ``group_label`` truthiness: unparseable-name
    uploads land in the "ungrouped" bucket with ``group_label=None`` while
    ``reps`` is still populated, so that case must still show a truthful
    "loaded" chip rather than "no data loaded". Contains the literal
    "fibrecv" the external screenshot harness waits on (``text=fibrecv``).

    ``noun`` names what is being counted ("replicate" / "angle") and
    ``label_prefix`` what ``group_label`` is ("group", or "" in multi-angle
    mode, where the label is already a fibre + part identity and "group
    C1_01 part1" would be plain wrong).
    """
    if n_reps == 0:
        chip = "no data loaded"
    else:
        reps_txt = f"{n_reps} {noun}" + ("" if n_reps == 1 else "s")
        edge_z_txt = _fmt(edge_z, "one_dp")
        tail = f"{feature_mode} · edge_z {edge_z_txt}"
        prefix = f"{label_prefix} " if label_prefix else ""
        chip = (f"{prefix}{group_label} · {reps_txt} · {tail}"
                if group_label else
                f"{reps_txt} (ungrouped) · {tail}")
    st.markdown(
        '<div class="fcv-header">'
        '<h1>fibrecv — fibre diameter detection</h1>'
        '<span class="fcv-sub">Local preview / tuning / batch / export over '
        'the validated pipeline.</span>'
        f'<span class="fcv-chip">{html.escape(chip)}</span>'
        '</div>',
        unsafe_allow_html=True)


#: replicate mode's card anchors, the default _render_jump_menu contents
_JUMP_REPLICATES = (("replicates", "01 Replicates"),
                    ("group-panel", "02 Group panel"),
                    ("tensile", "03 Tensile"),
                    ("export", "04 Export & batch"))


def _render_jump_menu(items: list[tuple[str, str]] | None = None) -> None:
    """Fixed right-hand nav to the card anchors; only called once data
    is loaded. ``items`` is a list of ``(anchor, label)`` pairs, defaulting to
    replicate mode's four cards. If a target anchor did not render (e.g. card
    03 is skipped because no group mean is available), its link is simply
    inert — the nav itself still degrades gracefully rather than erroring."""
    links = "".join(f'<a href="#{html.escape(anchor)}">{html.escape(label)}</a>'
                    for anchor, label in (items or _JUMP_REPLICATES))
    st.markdown(f'<nav class="fcv-jump">{links}</nav>', unsafe_allow_html=True)


# --------------------------------------------------------------------------- #
# Main-area renderers                                                         #
# --------------------------------------------------------------------------- #
def _render_replicate(rep: dict, cfg: CONFIG) -> None:
    mr, rgb = rep["mr"], rep["rgb"]
    bnd, edg, res = mr.bnd, mr.edg, mr.res
    edited_top, edited_bot = rep.get("edited_top"), rep.get("edited_bot")
    overlay = render_overlay(rgb, edg.y_top, edg.y_bot, bnd.c_fit, res.valid,
                             bnd.x0, bnd.x1, thick=1,
                             edited_top=edited_top, edited_bot=edited_bot,
                             slope=bnd.slope)
    caption = (f"{mr.name} — cyan top / yellow bottom / dashed centerline / "
               f"green = measured perpendicular diameter")
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
    flags.extend(res.anomaly.flags)
    if rep.get("rep_outlier"):
        flags.append("replicate_outlier")
    flags_txt = ", ".join(_flag_labels(flags)) if flags else "none"
    # rep["idx"] (position within the group) disambiguates names that collide
    # after sanitisation, e.g. "masp2 1_1_1.png" and "masp2_1_1_1.png" both
    # reduce to "masp2_1_1_1_png" via _safe_key alone -> StreamlitDuplicateElementKey.
    safe_name = f"{rep['idx']}_{_safe_key(rep['name'])}"
    with st.container(key=f"metrics_rep_{safe_name}"):
        c = st.columns(6)
        c[0].metric("mean Ø", _fmt(mean_um, "um"),
                    help="Mean diameter of this image, averaged along the fibre "
                         "(valid columns only).")
        c[1].metric("median Ø", _fmt(med, "um"))
        c[2].metric("along-fibre std", _fmt(std_um, "um"),
                    help="Std of the diameter along this image's fibre — thickness "
                         "variation within the picture, not disagreement between "
                         "replicates (that is the group panel's std).")
        c[3].metric("coverage", _fmt(res.coverage, "coverage"))
        c[4].metric("tilt slope", _fmt(bnd.slope, "tilt"))
        status_key = (f"status_ok_flags_{safe_name}" if not flags
                      else f"status_warn_flags_{safe_name}")
        with c[5]:
            with st.container(key=status_key):
                st.metric("flags", flags_txt, help=_FLAGS_HELP)

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
            edited_bot=rep.get("edited_bot"), slope=bnd.slope)
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
                                 min_value=-200.0, max_value=200.0,
                                 value=float(edits["nudge_top"]), step=0.5,
                                 format="%.1f", key=f"nudgetop_{name}")
        new_nb = n2.number_input("Nudge bottom line (px, + = down)",
                                 min_value=-200.0, max_value=200.0,
                                 value=float(edits["nudge_bot"]), step=0.5,
                                 format="%.1f", key=f"nudgebot_{name}")
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
        anomalies = _flag_labels(mr.res.anomaly.flags)
        if rep.get("rep_outlier"):
            anomalies.append(_FLAG_LABELS["replicate_outlier"])
        rows.append({
            "image": mr.name,
            "mean Ø (µm)": float(np.nanmean(d)) if ok else np.nan,
            "std (µm)": float(np.nanstd(d)) if ok else np.nan,
            "median Ø (µm)": float(np.nanmedian(d)) if ok else np.nan,
            "coverage (%)": mr.res.coverage * 100.0,
            "anomalies": "; ".join(anomalies) if anomalies else "—",
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
            "coverage (%)": st.column_config.NumberColumn(format="%.1f%%"),
        })


def _render_group(reps: list[dict], cfg: CONFIG,
                  group_label: str | None) -> float | None:
    """Registered mean ± std group panel; renders inside card 02.

    Returns the group's registered mean diameter (µm) so ``main()`` can gate
    card 03 (Tensile) on it — ``None`` on either early-return path (no
    replicate passed QC, or registration raised), matching exactly what used
    to gate the inline ``_render_tensile`` call this function made itself.
    The "Group panel" heading now lives in main()'s card-02 subheader; this
    function no longer renders its own.
    """
    profiles, dropped = [], []
    for rep in reps:
        mr = rep["mr"]
        res, bnd = mr.res, mr.bnd
        reason = exclusion_reason(res.band_mismatch, res.coverage,
                                  res.anomaly.flags, cfg)
        if reason is not None:
            dropped.append(f"{mr.name} ({_friendly_reason(reason)})")
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
        return None
    # sort with the same key register_sample uses, so zip(profiles, shifts)
    # below is order-aligned (replicate-keyed dicts would silently collide for
    # ungrouped uploads that share the idx+1 fallback number)
    profiles.sort(key=lambda p: p["replicate"])
    try:
        table, shifts, summary = register_sample(profiles, cfg)
    except Exception as exc:  # noqa: BLE001
        st.error(f"Registration failed: {exc}")
        _render_per_image_stats(reps)
        return None
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
    reg_uncertain = summary["registration_uncertain"]
    with st.container(key="metrics_group"):
        c = st.columns(6)
        c[0].metric("group mean Ø", _fmt(summary["mean_um"], "um"),
                    help="Mean diameter across the replicates of this group, "
                         "averaged along the aligned overlap region (where all "
                         "replicates are present after registration).")
        c[1].metric("between-replicate std", _fmt(summary["std_um"], "um"),
                    help="At each aligned position, the std of the diameter across "
                         "replicates; this is its average along the fibre. It "
                         "measures replicate-to-replicate disagreement, not "
                         "thickness variation along the fibre.")
        c[2].metric("CV", _fmt(cv, "cv"),
                    help="between-replicate std / group mean (dimensionless).")
        c[3].metric("reps used", _fmt(summary["n_replicates_used"], "int"))
        c[4].metric("overlap", _fmt(summary["overlap_px"], "px"),
                    help="Number of aligned columns where every replicate has data.")
        status_key = ("status_warn_registration" if reg_uncertain
                      else "status_ok_registration")
        with c[5]:
            with st.container(key=status_key):
                st.metric("registration", "uncertain" if reg_uncertain else "ok",
                          help="'uncertain' = the cross-correlation peak was below "
                               "min_corr for at least one replicate, so its shift "
                               "was reset to 0.")

    _render_per_image_stats(reps)
    st.caption(
        "Per-image mean/median/std are computed along each picture's own "
        "fibre (valid columns, no registration): std there = thickness "
        "variation within one picture. The group numbers above instead "
        "average the registered replicates, so their std = disagreement "
        "between replicates.")

    return summary["mean_um"]


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
        fallback = not np.isfinite(res.stress_pa).any()
        if fallback:
            st.caption("No matched diameter — showing force vs displacement.")

        # Six metrics in six narrow columns: the value font is shrunk (scoped
        # to this container's key, see _CSS's ".st-key-tensile_metrics" rules)
        # so figures like "161.00 mN" are not clipped with an ellipsis.
        with st.container(key="tensile_metrics"):
            c = st.columns(6)
            c[0].metric("breaking force", _fmt(res.fmax_n * 1000, "mN"),
                        help="Maximum load reached before the break (Fmax).")
            c[1].metric("tensile strength",
                        _fmt(res.tensile_strength_pa / 1e6, "MPa"),
                        help="Fmax / cross-sectional area.")
            c[2].metric("extension at break",
                        _fmt(res.extension_at_break_mm, "mm"))
            c[3].metric("strain at break",
                        _fmt(res.strain_at_break * 100, "pct100"))
            c[4].metric("Young's modulus",
                        _fmt(res.youngs_modulus_pa / 1e9, "GPa"))
            c[5].metric("toughness",
                        _fmt(res.toughness_j_m3 / 1e6, "MJ/m3"))

        area_um2 = res.area_m2 * 1e12 if np.isfinite(res.area_m2) else np.nan
        d_txt = _fmt(res.diameter_um, "um")
        a_txt = _fmt(area_um2, "um2")
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
    """Export/batch subsection; renders inside card 04 (main() owns the
    "04 Export & batch" heading, so the leading divider + subheader this
    function used to render itself are gone)."""
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
# Main-area renderers: multi-angle mode                                       #
# --------------------------------------------------------------------------- #
#: which nominal projection direction each angle measures (a1/a4 = 0 deg,
#: a2/a5 = 60 deg, a3/a6 = 120 deg) -- the ellipse fit needs all three
_MA_DIRECTIONS = "a1/a4 = 0°, a2/a5 = 60°, a3/a6 = 120°"


def _render_angle_tab(entry: dict, angle: int, cfg: CONFIG, um_per_px: float,
                      reason: str | None, shift: dict | None) -> None:
    """One angle's tab in card 01: overlay, metrics, optional editor, profile.

    Every µm number here comes from ``um_per_px`` (the sidebar's manual
    multi-angle scale), never from ``cfg.ppu`` — which is why this cannot
    reuse ``_render_replicate``. ``reason`` is this angle's QC exclusion
    string (None if it passed) and ``shift`` its entry from
    ``PartStack.shifts`` (None before a stack exists).
    """
    mr, rgb = entry["mr"], entry["rgb"]
    bnd, edg, res = mr.bnd, mr.edg, mr.res
    edited_top, edited_bot = entry.get("edited_top"), entry.get("edited_bot")
    overlay = render_overlay(rgb, edg.y_top, edg.y_bot, bnd.c_fit, res.valid,
                             bnd.x0, bnd.x1, thick=1,
                             edited_top=edited_top, edited_bot=edited_bot,
                             slope=bnd.slope)
    caption = (f"{entry['name']} — cyan top / yellow bottom / dashed "
               "centerline / green = measured perpendicular width")
    if edited_top is not None or edited_bot is not None:
        caption += " / magenta = manual edit"
    st.image(overlay, width="stretch", caption=caption)

    if reason:
        st.warning(f"Excluded from the cross-section fit — {reason}. Fix the "
                   "boundary below, or the fit runs without this angle.")
    if shift and shift["present"] and (shift["uncertain"] or shift["saturated"]):
        st.caption("⚠ this angle's alignment against the reference angle is "
                   "uncertain (weak or boundary-clamped correlation peak), so "
                   "it was stacked with a zero shift — see the per-angle table "
                   "in card 02.")

    span = slice(bnd.x0, bnd.x1 + 1)
    d_px = np.where(res.valid[span], res.diameter_raw[span], np.nan)
    have_d = bool(np.isfinite(d_px).any())
    med_px = float(np.nanmedian(d_px)) if have_d else None
    med_um = med_px * um_per_px if med_px is not None else None

    flags = []
    if res.low_confidence:
        flags.append("low_confidence")
    if res.band_mismatch:
        flags.append("band_mismatch")
    flags.extend(res.anomaly.flags)
    flags_txt = ", ".join(_flag_labels(flags)) if flags else "none"

    with st.container(key=f"metrics_angle_a{angle}"):
        c = st.columns(5)
        c[0].metric("median width", _fmt(med_um, "um"),
                    help="Median width of this view, over its measured "
                         "columns, converted with the sidebar's manual µm/px "
                         "scale (the ppu calibration is unused in this mode).")
        c[1].metric("median width (px)", _fmt(med_px, "px2"))
        c[2].metric("coverage", _fmt(res.coverage, "coverage"))
        c[3].metric("tilt slope", _fmt(bnd.slope, "tilt"))
        status_key = (f"status_ok_flags_a{angle}" if not flags
                      else f"status_warn_flags_a{angle}")
        with c[4]:
            with st.container(key=status_key):
                st.metric("flags", flags_txt, help=_FLAGS_HELP)

    # the editor is gated behind a checkbox: st.tabs runs all six tab bodies
    # on every rerun, so six always-open editors would mean six extra
    # full-res overlay renders plus six click-component iframes per
    # interaction on 2560x1920 images
    if st.checkbox("Edit boundaries", key=f"ma_edit_on_{entry['name']}",
                   help="Open the manual boundary editor for this angle. Kept "
                        "off by default because all six tabs render on every "
                        "interaction."):
        _render_edit_expander(entry)

    fig = _profile_fig(mr, rgb, cfg, um_per_px=um_per_px)
    st.pyplot(fig)
    plt.close(fig)


def _render_xsection(pv: MultiAnglePreview, excluded: dict[int, str],
                     um_per_px: float, cfg: CONFIG) -> None:
    """Card 02: angle coverage, aligned stack, ellipse metrics and QC table.

    Degrades in three steps rather than all-or-nothing: no angle passed QC ->
    an error listing the reasons; fewer than three projection directions, or
    no column where all three are measurable -> the aligned plot and the QC
    table still render, but no ellipse metrics; missing angles with the fit
    still possible -> everything renders, minus the split-half uncertainty
    (which needs each half to be a complete 3-direction fit of its own).
    """
    loaded = set(pv.present) | set(excluded)
    shifts = ({s["angle"]: s for s in pv.stack.shifts}
              if pv.stack is not None else {})

    # --- 1. coverage chips -------------------------------------------------
    missing: list[int] = []
    chips = []
    for a in range(1, 7):
        s = shifts.get(a, {})
        if a not in loaded:
            cls, note = "miss", "no image"
            missing.append(a)
        elif a in excluded:
            cls, note = "warn", "QC-excluded"
        elif s.get("present") and (s.get("uncertain") or s.get("saturated")):
            cls, note = "warn", "uncertain shift"
        else:
            cls, note = "ok", "in fit"
        chips.append(f'<span class="fcv-achip {cls}">a{a} · {note}</span>')
    st.markdown('<div class="fcv-achips">' + "".join(chips) + "</div>",
                unsafe_allow_html=True)

    if pv.stack is None or pv.fit is None:
        reasons = "; ".join(f"a{a}: {excluded[a]}" for a in sorted(excluded))
        st.error("Cannot fit a cross-section: no angle passed QC"
                 + (f" ({reasons})." if reasons else
                    " — no angle image is loaded for this part."))
        return

    med = pv.med
    n_present = len(pv.present)
    has_columns = bool(pv.fit.valid.any())
    if not pv.fittable:
        st.error(
            f"Cannot fit: only {pv.n_directions} of 3 projection directions "
            f"covered ({_MA_DIRECTIONS}). The aligned widths and the per-angle "
            "table below still describe what was measured, but a per-position "
            "ellipse needs at least one angle from each direction.")
    else:
        if n_present < 6:
            miss_txt = (", ".join(f"a{a}" for a in missing) or "none")
            exc_txt = (", ".join(f"a{a}" for a in sorted(excluded)) or "none")
            st.warning(
                f"Fitting with {n_present} of 6 angles — no image: {miss_txt}; "
                f"QC-excluded: {exc_txt}. The ellipse is still exactly "
                "determined, but the 180° pairs no longer cross-check it.")
        if not has_columns:
            st.error(
                "No column has all three projection directions after "
                "alignment, so no ellipse could be fitted anywhere along this "
                "part — the angles overlap too little, or their measured "
                "stretches do not.")

    # --- 2. aligned width stack -------------------------------------------
    fig = _xsec_stack_fig(pv.stack, pv.fit, um_per_px)
    st.pyplot(fig)
    plt.close(fig)
    # 2 dp to match the "px2" _FMT_KINDS precision (a sub-pixel shift is real
    # information here); the sign is kept, which _fmt cannot express
    shift_txt = ", ".join(f"a{s['angle']}: {s['shift_px']:+.2f} px"
                          for s in pv.stack.shifts if s["present"])
    st.caption(
        "Each angle's width profile is shifted along x (cross-correlation "
        "against the lowest present angle) so the same physical position is "
        f"compared across views. Applied shifts: {shift_txt}.")

    if pv.fittable and has_columns:
        # --- 3. ellipse metrics -------------------------------------------
        area_um2 = med["area"] * um_per_px ** 2
        err_um2 = med["area_err"] * um_per_px ** 2
        # split_half_area fits a1-a3 and a4-a6 independently, and each half
        # needs all three directions of its own -- one missing angle makes the
        # error NaN at every column, so the band is not merely absent, it is
        # undefined; say so instead of printing "± —" without explanation
        full_six = n_present == 6
        err_txt = _fmt(err_um2, "um2") if full_six else "—"
        with st.container(key="metrics_xsec"):
            c = st.columns(6)
            c[0].metric("median area", f"{_fmt(area_um2, 'um2')} ±{err_txt}",
                        help="Median over columns of the per-column ellipse "
                             "area, ± half the split-half (a1-a3 vs a4-a6) "
                             "difference. The batch CSV's A_mean_um2 is a MEAN "
                             "over columns of the whole fibre, so the two "
                             "numbers are close but never identical.")
            c[1].metric("axis ratio a/b", _fmt(med["ratio"], "cv"),
                        help="Median per-column ratio of the semi-major to the "
                             "semi-minor axis; 1.0 = circular.")
            c[2].metric("orientation φ", _fmt(med["phi"], "deg"),
                        help="Major-axis angle, circular median with period "
                             "180°. NaN columns (near-circular) are skipped.")
            c[3].metric("semi-axis a", _fmt(med["a"] * um_per_px, "um"))
            c[4].metric("semi-axis b", _fmt(med["b"] * um_per_px, "um"))
            rms = med["rms"]
            rms_warn = bool(np.isfinite(rms) and rms > cfg.xsec_rms_flag_px)
            with c[5]:
                with st.container(key=("status_warn_xsec_rms" if rms_warn
                                       else "status_ok_xsec_rms")):
                    st.metric("fit rms residual", _fmt(rms, "px2"),
                              help="Median per-column rms of the measured "
                                   "widths against the fitted ellipse. Above "
                                   f"{cfg.xsec_rms_flag_px:.0f} px the batch "
                                   "flags the fibre low-confidence.")
        if not full_six:
            st.caption("split-half uncertainty needs all six angles — with "
                       f"{n_present} present it is undefined at every column, "
                       "so the ± figure and the area band are omitted.")

        # --- 4. cross-section drawing + area along the part ----------------
        col_e, col_a = st.columns([1, 2])
        with col_e:
            efig = _xsec_ellipse_fig(med, um_per_px)
            if efig is not None:
                st.pyplot(efig)
                plt.close(efig)
                st.caption("Median cross-section drawn to scale.")
        with col_a:
            afig = _xsec_area_fig(pv.stack, pv.fit, um_per_px, show_err=full_six)
            st.pyplot(afig)
            plt.close(afig)
            st.caption("Per-column ellipse area with the ±split-half "
                       "uncertainty band.")

    # --- 5. per-angle QC table --------------------------------------------
    rows = []
    for a in range(1, 7):
        s = shifts.get(a, {})
        present = bool(s.get("present"))
        r = pv.fit.resid[a - 1]
        finite = np.isfinite(r)
        rows.append({
            "angle": f"a{a}",
            "in fit": present,
            "shift (px)": float(s["shift_px"]) if present else np.nan,
            "corr peak": float(s["corr_peak"]) if present else np.nan,
            "uncertain": bool(present and s.get("uncertain")),
            "saturated": bool(s.get("saturated")),
            "median signed residual (px)": (float(np.median(r[finite]))
                                            if finite.any() else np.nan),
            "columns": int(finite.sum()),
        })
    st.markdown("**Per-angle alignment & residuals**")
    st.dataframe(
        pd.DataFrame(rows), width="stretch", hide_index=True,
        column_config={
            # same 2 dp as the "px2" kind / the shift caption above
            "shift (px)": st.column_config.NumberColumn(format="%.2f"),
            "corr peak": st.column_config.NumberColumn(format="%.3f"),
            "median signed residual (px)": st.column_config.NumberColumn(
                format="%.2f",
                help="Systematically positive or negative here means this "
                     "view reads consistently wider or narrower than the "
                     "fitted ellipse predicts."),
        })

    # --- 6. secondary numbers ---------------------------------------------
    n_links = int(med.get("n_links", 0))
    st.caption(
        f"pair Δw/w {_fmt(med.get('dw_frac'), 'cv')} · valid columns "
        f"{_fmt(med.get('valid_frac'), 'coverage')} of the measurable ones · "
        f"uncertain shifts {int(med.get('n_uncertain', 0))}/{n_links} "
        f"({int(med.get('n_saturated', 0))} saturated) · "
        f"Scale: {um_per_px:.4f} µm/px (manual).")


def _render_xsec_batch(ma: dict, cfg: CONFIG) -> None:
    """Card 03: run the whole folder through measure + cross-sections.

    One button for both stages (``run_multiangle_batch``): the per-image
    measure pass reports through the progress bar, ``run_xsection`` runs
    behind the spinner. Folder source only — the batch reads the images back
    off disk by path, which uploads have none of.
    """
    out_folder = st.text_input(
        "Output folder",
        value=st.session_state.get("out_folder", "./fibrecv_output"))
    st.session_state.out_folder = out_folder

    folder = ma["folder"]
    condition = ma["condition"]
    folder_ok = folder is not None and Path(folder).is_dir()
    can_run = bool(folder_ok and condition)
    jobs = int(st.number_input("parallel jobs", min_value=1, max_value=16,
                               value=4, step=1))

    st.caption(
        "Tune the parameters on a few images in the preview above first — this "
        "runs the whole folder once with the current parameters (roughly 20 "
        "minutes for 450 images), measuring every angle image and then fitting "
        "one cross-section per (fibre, part).")
    st.caption(
        f"Scoped to condition {condition or '—'}; images of other conditions "
        "in the same folder are skipped. The batch recomputes every image from "
        "disk and IGNORES the manual edits made in card 01. Point it at a "
        "fresh output folder — it reuses per_image/* from whatever is already "
        "there.")
    if not can_run:
        st.caption("Needs a *Local folder* source with multi-angle filenames.")

    if st.button("Run multi-angle batch (measure + cross-sections)",
                 disabled=not can_run, type="primary", width="stretch"):
        try:
            groups = discover_multiangle(folder, condition)
            n = sum(len(angles) for angles in groups.values())
            if not n:
                st.warning("No multi-angle images found for this condition.")
                return
            prog = st.progress(0.0, text=f"Measuring 0/{n}…")

            def _cb(frac: float, _r: dict) -> None:
                prog.progress(frac,
                              text=f"Measuring {int(round(frac * n))}/{n}…")

            with st.spinner(f"Measuring {n} images, then fitting "
                            f"{len(groups)} cross-section(s)…"):
                summary, results, rc = run_multiangle_batch(
                    folder, condition, out_folder, cfg, ma["um_per_px"],
                    jobs, _cb)
            prog.empty()
            if rc != 0:
                st.error(f"Cross-section stage failed (exit code {rc}) — the "
                         "per-image measurements were still written to "
                         f"{out_folder}.")
                return
            n_err = sum(1 for r in results if "error" in r)
            st.success(f"Measured {len(results)} images ({n_err} errors), "
                       f"{len(summary)} fibre row(s) → {out_folder}")
            if summary.empty:
                st.info("No fibre could be cross-sectioned (check the angle "
                        "coverage and the per-image QC).")
                return
            table = summary.copy()
            if "low_confidence" in table.columns:
                low = table["low_confidence"].astype(bool)
                table.insert(0, "status",
                             np.where(low, "⚠ low confidence", "ok"))
            else:
                low = None
            # dataframe + download must render before the low-confidence
            # fibre-list lookup below: a schema surprise in that lookup (e.g.
            # a renamed/missing "fiber" column) must not raise and mask a
            # batch that actually succeeded and already wrote its files.
            st.dataframe(table, width="stretch", hide_index=True)
            st.download_button(
                "Download xsection_summary.csv",
                summary.to_csv(index=False).encode(),
                file_name="xsection_summary.csv", mime="text/csv")
            try:
                if low is not None:
                    flagged = [f"{int(f):02d}" for f in table.loc[low, "fiber"]]
                    if flagged:
                        st.warning(
                            "Low-confidence fibre(s): " + ", ".join(flagged)
                            + " — too few valid columns, too many uncertain "
                              "shifts, or an rms residual above "
                              f"{cfg.xsec_rms_flag_px:.0f} px. Check their "
                              "angle coverage before using these areas.")
            except Exception:  # noqa: BLE001 - degrade to skipping the
                # warning; the dataframe and download above already rendered.
                pass
        except Exception as exc:  # noqa: BLE001
            st.error(f"Multi-angle batch failed: {exc}")
            st.code(traceback.format_exc())


def _run_multiangle_mode(cfg_items: tuple, cfg: CONFIG) -> None:
    """The multi-angle branch of ``main()``: sidebar, header, three cards."""
    ma = _load_multiangle(cfg_items)
    _param_form(multiangle=True)
    angles: dict[int, dict] = ma["angles"]
    um_per_px = ma["um_per_px"]

    # manual boundary edits, applied once per angle right here -- same choke
    # point as the replicate flow, so corrections reach the overlay, the
    # profile plot and (through _profiles_from_results) the ellipse fit. The
    # replicate-outlier pass has no meaning across angles and is skipped.
    for entry in angles.values():
        edits = st.session_state.manual_edits.get(entry["name"])
        if has_edits(edits):
            entry["mr"], entry["edited_top"], entry["edited_bot"] = \
                apply_manual_edits(entry["mr"], edits, cfg)

    part_label = (f"{ma['condition']}_{ma['fiber']:02d} part{ma['part']}"
                  if ma["condition"] is not None else None)
    _render_header(part_label, len(angles), cfg.edge_z, cfg.feature_mode,
                   noun="angle", label_prefix="")

    if not angles:
        st.info("Pick a folder and a fibre, or upload the six angle images of "
                "one part, to begin.")
        return

    profiles, excluded = _profiles_from_results(
        {a: e["mr"] for a, e in angles.items()}, cfg)
    pv = multiangle_preview(profiles, cfg)

    _render_jump_menu([("angles", "01 Angles"),
                       ("xsection", "02 Cross-section"),
                       ("xsec-export", "03 Batch & export")])

    shifts = ({s["angle"]: s for s in pv.stack.shifts}
              if pv.stack is not None else {})

    def _tab_label(a: int) -> str:
        if a not in angles:
            return f"a{a}"
        s = shifts.get(a, {})
        warn = (a in excluded
                or bool(angles[a]["mr"].res.anomaly.flags)
                or bool(s.get("present")
                        and (s.get("uncertain") or s.get("saturated"))))
        return f"⚠ a{a}" if warn else f"a{a}"

    with st.container(key="card_angles"):
        st.subheader(f"01 Angles — {part_label}", anchor="angles")
        st.caption("The six rotated views of this fibre part; widths use the "
                   "sidebar's manual µm/px scale, not the ppu calibration.")
        tabs = st.tabs([_tab_label(a) for a in range(1, 7)])
        for a, tab in zip(range(1, 7), tabs):
            with tab:
                if a not in angles:
                    st.info(f"No image for angle a{a} of this part.")
                    continue
                _render_angle_tab(angles[a], a, cfg, um_per_px,
                                  excluded.get(a), shifts.get(a))

    with st.container(key="card_xsec"):
        st.subheader("02 Cross-section", anchor="xsection")
        st.caption("Per-position ellipse from the aligned angle widths "
                   f"({_MA_DIRECTIONS}).")
        _render_xsection(pv, excluded, um_per_px, cfg)

    with st.container(key="card_xsec_export"):
        st.subheader("03 Batch & export", anchor="xsec-export")
        _render_xsec_batch(ma, cfg)


# --------------------------------------------------------------------------- #
# App entry point                                                             #
# --------------------------------------------------------------------------- #
def main() -> None:
    st.set_page_config(page_title="fibrecv — fibre diameter GUI", layout="wide")
    st.session_state.setdefault("cfg_dict", DEFAULTS.as_dict())
    st.session_state.setdefault("form_version", 0)
    st.session_state.setdefault("manual_edits", {})  # image name -> edits dict

    _inject_css()

    # the analysis-mode radio runs BEFORE the cfg snapshot: entering
    # multi-angle mode re-applies the bright defaults, and that must be in
    # cfg_items for this same run rather than one rerun late
    mode = _analysis_mode()

    cfg_items = _cfg_items(st.session_state.cfg_dict)
    cfg = _cfg_from_items(cfg_items)

    if mode == _MODE_MULTIANGLE:
        _run_multiangle_mode(cfg_items, cfg)
        return

    # sidebar — _load_reps must run first: the smoke test's
    # at.sidebar.text_input[0] is the "Image folder" field it renders
    reps, group_label, folder = _load_reps(cfg_items)
    _param_form()
    tensile = _tensile_controls()

    # header renders after the sidebar builders since its state chip needs
    # the group/replicate state they just produced
    _render_header(group_label, len(reps), cfg.edge_z, cfg.feature_mode)

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

    # group-level replicate_outlier check (advisory: badges/table only, never
    # exclusion) — after manual edits so corrected medians feed the comparison;
    # keyed by idx, not name: uploaded files can share a stem
    rep_devs, rep_outliers = detect_replicate_outliers(
        {rep["idx"]: rep["mr"].meta.get("median_diameter_um") for rep in reps},
        cfg)
    for rep in reps:
        rep["rep_dev"] = rep_devs.get(rep["idx"])
        rep["rep_outlier"] = rep["idx"] in rep_outliers

    # main area
    if not reps:
        st.info("Pick a folder + group, or upload images, to begin.")
        return

    _render_jump_menu()

    with st.container(key="card_replicates"):
        st.subheader(
            f"01 Replicates — group {group_label}" if group_label
            else "01 Replicates (uploaded)", anchor="replicates")
        def _tab_label(r: dict) -> str:
            base = (f"Rep {r['mr'].replicate}"
                    if r["mr"].replicate is not None else r["name"])
            warn = r["mr"].res.anomaly.flags or r.get("rep_outlier")
            return f"⚠ {base}" if warn else base

        tabs = st.tabs([_tab_label(r) for r in reps])
        for tab, rep in zip(tabs, reps):
            with tab:
                _render_replicate(rep, cfg)

    with st.container(key="card_group"):
        st.subheader("02 Group panel", anchor="group-panel")
        st.caption("Registered mean ± std across the group's aligned replicates.")
        mean_um = _render_group(reps, tcfg, group_label)

    # card 03 (Tensile) only when the group panel produced a mean — identical
    # gating to the pre-redesign code, where _render_tensile was called from
    # inside _render_group only on its non-early-return path
    if mean_um is not None:
        with st.container(key="card_tensile"):
            st.subheader("03 Tensile", anchor="tensile")
            _render_tensile(group_label, mean_um, tcfg, tensile["tmap"])

    with st.container(key="card_export"):
        st.subheader("04 Export & batch", anchor="export")
        _render_export_batch(reps, tcfg, group_label, folder, tensile["tmap"], cfg_items)


if __name__ == "__main__":
    main()
