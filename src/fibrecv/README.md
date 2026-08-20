# fibrecv/

Once my folder has changes, please update me.

Per-image detection pipeline (features → band → edges → qc) driven by
`compute`/`measure`/`run_measure`, replicate registration (`register`/
`run_aggregate`), multi-angle cross-sections (`xsection`/`run_xsection`/
`scale`), tensile ingestion (`tensile`), and the front ends (three CLIs +
`gui_app.py`/`gui_launch.py` Streamlit GUI) sharing one compute core.

## Files

- `__init__.py` — current — package docstring naming the three-stage
  pipeline (`run_measure` → `run_aggregate` | `run_xsection`) and
  `__version__`.
- `config.py` — current — the `CONFIG` dataclass carrying every tunable
  detection/registration/xsection/tensile parameter (incl. `feature_mode`,
  `xsec_min_corr`/`xsec_max_shift`) and the `px_to_um` helper; imported by
  every other module.
- `io_utils.py` — current — image discovery (`discover_images`), filename
  parsing into (group, replicate) (`parse_name`), the strict multi-angle
  convention (`MultiAngleKey`, `parse_multiangle_name`, `multiangle_group`,
  `discover_multiangle`), and RGB loading (`load_rgb`, with a pillow
  fallback for LZW TIFFs).
- `features.py` — current — RGB → feature z-map (`rgb_to_desaturation`),
  mode-dispatched: `"desat"` (MasP2 pale-on-pink saturation) or `"bright"`
  (C1 bright-on-dark brightness), same margin-median/MAD machinery.
- `scale.py` — current — Zeiss XML sidecar reader (`read_scale`,
  `sidecar_for`, `resolve_um_per_px`): both µm/px candidates (camera pitch
  vs Scaling/Items) for the xsection stage's µm conversion.
- `xsection.py` — current — pure multi-angle cross-section math:
  `build_part_stack` (cross-angle alignment), `fit_ellipse_projections`
  (w²-space linear ellipse fit), `split_half_area` uncertainty,
  `predict_anisotropy`/`predict_phi_transfer` validation criteria — see
  `docs/features/03_multiangle_xsection.md`.
- `band.py` — current — locates the single full-width fibre band and fits
  its centerline/tilt (`BandResult`, `tilt_geometry`).
- `edges.py` — current — per-column tight-inner-edge detection producing
  sub-pixel top/bottom boundaries and the tilt-corrected perpendicular
  diameter (`EdgeResult`, `detect_edges`).
- `anomaly.py` — current — pure anomaly detectors: per-image edge_jump /
  large_gap / diameter_step (`detect_image_anomalies` → `AnomalyResult`),
  group-level `detect_replicate_outliers`, and the shared `exclusion_reason`
  drop policy used by both `run_aggregate` and the GUI.
- `qc.py` — current — outlier rejection, smoothing and coverage gating on
  the raw edge diameters (`QCResult`), plus the advisory anomaly pass
  (`QCResult.anomaly`, detectors in `anomaly.py`).
- `compute.py` — current — pure per-image compute core chaining
  features → band → edges → qc with no file I/O (`compute_measurement`);
  shared by the CLI and the GUI so preview == CLI output.
- `manual_edit.py` — current — headless manual boundary correction:
  anchors/nudges → corrected edges + re-run QC (`apply_manual_edits`); no
  Streamlit dependency.
- `overlay.py` — current — draws the boundary-overlay PNG (top/bottom
  edges, centerline, green perpendicular measurement chords)
  (`render_overlay`, `draw_overlay`, `draw_perp_chords`).
- `measure.py` — current — per-image orchestration: runs `compute`, writes
  the overlay/CSV/plot/meta artifacts for one image (`write_measurement`).
- `register.py` — current — aligns a group's replicate profiles into a
  registered mean ± std curve (`register_sample`).
- `run_measure.py` — current — CLI, stage 1: measures a glob of images in
  parallel (`ProcessPoolExecutor`), writing per-image artifacts.
- `run_aggregate.py` — current — CLI, stage 2: groups per-image profiles by
  sample and builds the registered averages + `master_summary.csv` +
  `per_image_summary.csv` (every image with its anomaly flags and
  exclusion reason; `--anomaly-exclude` turns flags into exclusions).
- `run_xsection.py` — current — CLI, stage 3 (multi-angle sets): aligns the
  six per-angle width profiles of each (fiber, part), fits per-column
  ellipses, converts to µm via the XML sidecar scale
  (`--scale-source`; a numeric literal bypasses the sidecars entirely, for
  manual-scale GUI runs), optional `--anomaly-exclude` (mirrors
  `run_aggregate`'s flag), and writes per-part CSV/plot/shifts +
  `xsection_summary.csv`/`xsection_angle_residuals.csv`
  (+`xsection_validation.csv` with `--validation`).
- `tensile.py` — current — tensile (stress-strain) CSV/Excel ingestion and
  single-fibre metric computation (modulus, toughness, break point) from a
  fibre's measured mean diameter.
- `gui_app.py` — current — local Streamlit GUI: tuning/preview/manual-edit/
  batch/export over the same `compute`/`register`/`tensile` core, plus the
  "clean lab" visual layer (`_CSS`/`_inject_css`, header+chip, numbered
  card sections with a jump menu, `_fmt`, `_styled_fig`) — see
  `docs/features/02_gui-redesign.md`. Sidebar image-mode selector
  (desat | bright); a mode switch re-applies that mode's calibrated
  `edge_frac`/`k_band` defaults (bright via `run_measure.BRIGHT_DEFAULTS`).
  Sidebar analysis-mode radio (`_analysis_mode`) switching the whole main area
  between replicate mode and a multi-angle cross-section mode: `_load_multiangle`
  (condition/fibre/part selectors + manual µm/px scale), 01 Angles (six angle
  tabs, `_render_angle_tab`), 02 Cross-section (`_render_xsection`,
  `_xsec_stack_fig`/`_xsec_ellipse_fig`/`_xsec_area_fig`) and 03 Batch & export
  (`_render_xsec_batch`). Underneath, pure headless multi-angle helpers over the
  `xsection`/`run_xsection` stage-3 pipeline (`_profile_from_mr`,
  `_profiles_from_results`, `multiangle_preview` -> `MultiAnglePreview`,
  `run_multiangle_batch`) plus `run_batch(aggregate=)` to skip stage-2
  registration for angle data.
- `gui_launch.py` — current — `fibrecv-gui` console-script entry point; a
  thin subprocess wrapper around `streamlit run gui_app.py`.
