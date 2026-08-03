# fibrecv/

Once my folder has changes, please update me.

Per-image detection pipeline (features → band → edges → refine → qc) driven by
`compute`/`measure`/`run_measure`, replicate registration (`register`/
`run_aggregate`), tensile ingestion (`tensile`), and the two front ends
(`run_measure.py`/`run_aggregate.py` CLIs, `gui_app.py`/`gui_launch.py`
Streamlit GUI) that call the same core so CLI and GUI outputs match.

## Files

- `__init__.py` — current — package docstring naming the two-stage
  pipeline (`run_measure` → `run_aggregate`) and `__version__`.
- `config.py` — current — the `CONFIG` dataclass carrying every tunable
  detection/registration/tensile parameter and the `px_to_um` helper;
  imported by every other module.
- `io_utils.py` — current — image discovery (`discover_images`), filename
  parsing into (group, replicate) (`parse_name`), and RGB loading
  (`load_rgb`).
- `features.py` — current — RGB → desaturation z-map (`rgb_to_desaturation`),
  the self-normalising per-image signal detection is built on.
- `band.py` — current — locates the single full-width fibre band and fits
  its centerline/tilt (`BandResult`, `tilt_geometry`).
- `edges.py` — current — per-column tight-inner-edge detection producing
  sub-pixel top/bottom boundaries and the tilt-corrected perpendicular
  diameter (`EdgeResult`, `detect_edges`).
- `refine.py` — current — erf edge-refinement stage: refits each detected
  wall as a Gaussian-blurred step over blocks of columns and shifts the
  boundary to the fitted midpoint, sub-pixel and blur-invariant; falls back
  to the legacy edge per block/column when a fit doesn't pass its gates, and
  is bit-identical to the pre-refinement edges when `refine_on=False`
  (`RefineResult`, `refine_edges`) — see
  `docs/features/03_esf_edge_refinement.md`.
- `qc.py` — current — outlier rejection, smoothing and coverage gating on
  the raw edge diameters (`QCResult`).
- `compute.py` — current — pure per-image compute core chaining
  features → band → edges → refine → qc with no file I/O (`compute_measurement`);
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
  parallel (`ProcessPoolExecutor`), writing per-image artifacts; exposes every
  `CONFIG` knob as a flag, including `--refine`/`--no-refine` for the erf
  edge-refinement stage.
- `run_aggregate.py` — current — CLI, stage 2: groups per-image profiles by
  sample and builds the registered averages + `master_summary.csv`.
- `tensile.py` — current — tensile (stress-strain) CSV/Excel ingestion and
  single-fibre metric computation (modulus, toughness, break point) from a
  fibre's measured mean diameter.
- `gui_app.py` — current — local Streamlit GUI: tuning/preview/manual-edit/
  batch/export over the same `compute`/`register`/`tensile` core (sidebar
  params include the `refine_on` erf edge-refinement toggle, and the
  per-replicate profile plot overlays a sigma(x) fit-quality trace when
  refinement ran), plus the "clean lab" visual layer (`_CSS`/`_inject_css`,
  header+chip, four numbered card sections with a jump menu, `_fmt`,
  `_styled_fig`) — see `docs/features/02_gui-redesign.md`.
- `gui_launch.py` — current — `fibrecv-gui` console-script entry point; a
  thin subprocess wrapper around `streamlit run gui_app.py`.
