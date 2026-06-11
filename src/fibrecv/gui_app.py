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

from fibrecv import run_aggregate  # noqa: E402
from fibrecv.compute import compute_measurement  # noqa: E402
from fibrecv.config import CONFIG  # noqa: E402
from fibrecv.io_utils import discover_images, natural_key, parse_name  # noqa: E402
from fibrecv.io_utils import load_rgb as _io_load_rgb  # noqa: E402
from fibrecv.manual_edit import (  # noqa: E402
    apply_manual_edits, display_to_native, empty_edits, has_edits)
from fibrecv.measure import write_measurement  # noqa: E402
from fibrecv.overlay import GREY, WHITE, mark_anchors, render_overlay  # noqa: E402
from fibrecv.register import register_sample  # noqa: E402
from fibrecv.run_measure import _lib_versions, _worker  # noqa: E402
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
# Sidebar: data source -> list of replicate dicts                             #
# --------------------------------------------------------------------------- #
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
        uploads = st.sidebar.file_uploader(
            "Upload images", type=["jpg", "jpeg", "png", "tif", "tiff", "bmp"],
            accept_multiple_files=True)
        if uploads:
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
    flags = []
    if res.low_confidence:
        flags.append("low_confidence")
    if res.band_mismatch:
        flags.append("band_mismatch")
    c = st.columns(4)
    c[0].metric("median Ø", f"{med:.2f} µm" if med is not None else "—")
    c[1].metric("coverage", f"{res.coverage:.0%}")
    c[2].metric("tilt slope", f"{bnd.slope:.4f}")
    c[3].metric("flags", ", ".join(flags) if flags else "none")

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


def _render_group(reps: list[dict], cfg: CONFIG, group_label: str | None) -> None:
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
        return
    # sort with the same key register_sample uses, so zip(profiles, shifts)
    # below is order-aligned (replicate-keyed dicts would silently collide for
    # ungrouped uploads that share the idx+1 fallback number)
    profiles.sort(key=lambda p: p["replicate"])
    try:
        table, shifts, summary = register_sample(profiles, cfg)
    except Exception as exc:  # noqa: BLE001
        st.error(f"Registration failed: {exc}")
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


def _render_export_batch(reps: list[dict], cfg: CONFIG, group_label: str | None,
                         folder: str | None) -> None:
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

    with col_batch:
        st.markdown("**Run batch (whole folder)**")
        if edited_names:
            st.warning("Run batch recomputes every image from disk and "
                       "IGNORES the manual edits made here.")
        jobs = st.number_input("parallel jobs", min_value=1, max_value=16, value=4, step=1)
        disabled = folder is None or not Path(folder or "").is_dir()
        if disabled:
            st.caption("Batch needs the *Local folder* source.")
        if st.button("Run batch", disabled=disabled, width="stretch"):
            images = discover_images(folder)
            if not images:
                st.warning("No image files found in the folder.")
                return
            prog = st.progress(0.0, text=f"Measuring 0/{len(images)}…")
            n = len(images)

            def _cb(frac: float, _r: dict) -> None:
                prog.progress(frac, text=f"Measuring {int(round(frac * n))}/{n}…")

            try:
                with st.spinner("Batch running…"):
                    master, results = run_batch(images, out_folder, cfg, int(jobs), _cb)
                prog.empty()
                n_err = sum(1 for r in results if "error" in r)
                st.success(f"Batch done: {len(results)} images, {n_err} errors → {out_folder}")
                if not master.empty:
                    st.dataframe(master, width="stretch")
                    st.download_button(
                        "Download master_summary.csv",
                        master.to_csv(index=False).encode(),
                        file_name="master_summary.csv", mime="text/csv")
                else:
                    st.info("No groups aggregated (check filenames / coverage).")
            except Exception as exc:  # noqa: BLE001
                prog.empty()
                st.error(f"Batch failed: {exc}")
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

    _render_group(reps, cfg, group_label)
    _render_export_batch(reps, cfg, group_label, folder)


if __name__ == "__main__":
    main()
