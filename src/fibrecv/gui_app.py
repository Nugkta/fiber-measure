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
- The three boundary knobs (``edge_z``/``edge_frac``/``wcol``), edited in a
  sidebar form and applied on demand; all other ``CONFIG`` fields stay at the
  validated defaults.
- An output-folder path for export/batch.

Output
------
- Live, in-memory preview: full-res boundary overlays, per-replicate diameter
  profiles, and a registered mean+/-std group curve -- all redrawn when the user
  changes parameters and clicks Apply (no disk writes for preview).
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
from fibrecv.measure import write_measurement  # noqa: E402
from fibrecv.overlay import render_overlay  # noqa: E402
from fibrecv.register import register_sample  # noqa: E402
from fibrecv.run_measure import DEFAULT_ROOT, _lib_versions, _worker  # noqa: E402

DEFAULTS = CONFIG()  # never mutated; the source of widget defaults + reset target

# --- visible parameters: the three knobs that move the detected boundary.
# Everything else in CONFIG stays at the validated defaults (CLI keeps full
# control). spec = (name, kind, help, step, lo, hi, fmt). ---
PARAM_SPECS: list[tuple] = [
    ("edge_z", "slider",
     "边界松紧（最重要）。调大 → 边界向纤维内侧收紧；调小 → 边界向外放松。"
     "如果画出的线被拽进纤维内部（比实际细），就调小；"
     "如果线跑到纤维外面的阴影上（比实际粗），就调大。推荐 3–5，默认 4.0。",
     0.5, 1.0, 12.0, "%.1f"),
    ("edge_frac", "float",
     "淡纤维保护。纤维颜色很淡、对比度低时防止边界跑飞。"
     "一般保持默认 0.65 不用动；只有很淡的纤维检测不到边时才微调。",
     0.05, 0.0, 1.0, "%.2f"),
    ("wcol", "int",
     "横向平滑宽度（像素）。边界线抖动、锯齿、毛刺严重时调大（如 31、51），"
     "线会更平滑稳定；调太大会抹平真实的粗细变化。默认 15。",
     1, 1, 201, None),
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


def _group_fig(table: dict, group_label: str):
    """Registered mean diameter curve with a +/-std shaded band."""
    x = table["x_aligned_px"]
    mean = table["mean_um"]
    std = table["std_um"]
    fig, ax = plt.subplots(figsize=(9, 3))
    ax.plot(x, mean, "-", lw=1.5, color="tab:blue", label="mean")
    band = np.where(np.isfinite(std), std, 0.0)
    ax.fill_between(x, mean - band, mean + band, alpha=0.25, color="tab:blue", label="±std")
    ax.set_xlabel("aligned x position (px)")
    ax.set_ylabel("diameter (µm)")
    ax.set_title(f"sample {group_label}  (n_reps={int(np.nanmax(table['n']))})")
    ax.legend(loc="best", fontsize=8)
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
            "folder", DEFAULT_ROOT))
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
    overlay = render_overlay(rgb, edg.y_top, edg.y_bot, bnd.c_fit, res.valid,
                             bnd.x0, bnd.x1, thick=1)
    st.image(overlay, width="stretch",
             caption=f"{mr.name} — cyan top / yellow bottom / dashed centerline")

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

    fig = _profile_fig(mr, rgb, cfg)
    st.pyplot(fig)
    plt.close(fig)


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
    try:
        table, shifts, summary = register_sample(profiles, cfg)
    except Exception as exc:  # noqa: BLE001
        st.error(f"Registration failed: {exc}")
        return

    fig = _group_fig(table, group_label or "uploaded")
    st.pyplot(fig)
    plt.close(fig)

    cv = summary["cv"]
    c = st.columns(6)
    c[0].metric("mean Ø", f"{summary['mean_um']:.2f} µm")
    c[1].metric("std", f"{summary['std_um']:.2f} µm")
    c[2].metric("CV", f"{cv:.3f}" if np.isfinite(cv) else "—")
    c[3].metric("n reps used", f"{summary['n_replicates_used']}")
    c[4].metric("overlap", f"{summary['overlap_px']} px")
    c[5].metric("registration", "uncertain" if summary["registration_uncertain"] else "ok")


def _render_export_batch(reps: list[dict], cfg: CONFIG, group_label: str | None,
                         folder: str | None) -> None:
    st.divider()
    st.subheader("Export & batch")
    out_folder = st.text_input("Output folder",
                               value=st.session_state.get("out_folder", "./fibrecv_output"))
    st.session_state.out_folder = out_folder

    col_exp, col_batch = st.columns(2)

    with col_exp:
        st.markdown("**Export current group**")
        st.caption("Writes overlays/, per_image/* and per_sample/* for the loaded group.")
        if st.button("Export current group", disabled=not reps, width="stretch"):
            try:
                with st.spinner("Writing output tree…"):
                    g = export_group(reps, out_folder, cfg)
                st.success(f"Exported group {g} → {out_folder}")
            except Exception as exc:  # noqa: BLE001
                st.error(f"Export failed: {exc}")

    with col_batch:
        st.markdown("**Run batch (whole folder)**")
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

    st.title("fibrecv — fibre diameter detection")
    st.caption("Local preview / tuning / batch / export over the validated pipeline. "
               f"Strictness knob edge_z = {st.session_state.cfg_dict['edge_z']}.")

    cfg_items = _cfg_items(st.session_state.cfg_dict)
    cfg = _cfg_from_items(cfg_items)

    # sidebar
    reps, group_label, folder = _load_reps(cfg_items)
    _param_form()

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
