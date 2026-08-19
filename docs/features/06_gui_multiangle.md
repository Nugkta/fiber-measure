# GUI multi-angle cross-section mode
Added: 2026-08-19
Code: `src/fibrecv/gui_app.py`, `src/fibrecv/run_xsection.py`, `src/fibrecv/xsection.py`, `tests/test_gui_multiangle.py`

## What it does
Adds a second analysis mode to the Streamlit GUI, next to the existing replicate
workflow, that wires the CLI's multi-angle cross-section stage (feature
`03_multiangle_xsection`) into interactive preview. A sidebar radio
("Replicates" / "Multi-angle cross-section") switches the whole main area.
In multi-angle mode the user points at a folder (or uploads files) named
`<condition>_<fibre>_a<angle>_part<part>.tiff`, picks a condition/fibre/part,
and sees three card sections: **01 Angles** (six tabs, one per rotation image,
each with its overlay, per-angle metrics in µm, and an optional manual
boundary editor), **02 Cross-section** (the six aligned width profiles, the
fitted ellipse's area/axis-ratio/orientation with a split-half uncertainty,
a drawn cross-section, an area-vs-position plot, and a per-angle QC table),
and **03 Batch & export** (one button that measures every image for the whole
condition and runs the cross-section fit, with a progress bar and a
downloadable `xsection_summary.csv`). Entering the mode auto-switches to the
`bright` feature-mode defaults (C1 images are bright-on-dark), since the
cross-section fit is only meaningful with a bright-mode boundary.

## Design choices
- **Manual µm/px scale instead of `ppu`** — chosen because the multi-angle
  pipeline's true scale comes from the C1 microscope's XML sidecar
  (`Scaling/Items`), not the GUI's `ppu` calibration field, and the GUI has no
  sidecar reader wired to its preview path. A `number_input` (default
  `0.388924`, the C1 study-03 resolved `Scaling/Items` value, matching the
  `s03d` run's `0.38892367`) sits outside the parameter form so changing it
  never triggers a recompute — every µm figure in this mode is derived from
  it at render time. Rejected: reading the real sidecar in the GUI, which
  would need per-image sidecar discovery logic the preview path does not
  otherwise need. The `ppu` field from the parameter form stays visible in
  this mode (the form is shared with replicate mode) but is unused —
  documented in the field's own help text and in tooltips throughout the
  mode.
- **Numeric `--scale-source` bypasses XML sidecars entirely** — `run_xsection`'s
  `_part_scale` now tries `float(scale_source)` first; on success it returns
  the manual value directly with no sidecar reads at all, which also means it
  skips the 0.1% inter-angle sidecar-agreement check the `"camera"`/`"items"`
  string paths still enforce. This is a genuine behaviour change on the CLI's
  numeric path (a user who already passed a hardcoded number as
  `--scale-source` before this change was still silently required to have
  matching sidecars present; now they are not read at all), accepted because
  the GUI's batch card always calls it with a numeric literal and reading
  sidecars for a value the user already supplied by hand would be pure
  overhead. Non-positive values raise `ValueError` (CLI exit code 2).
- **≥3 of 3 projection directions required to fit** — the six angles measure
  three physical directions twice each (a1/a4 = 0°, a2/a5 = 60°, a3/a6 =
  120°); the per-column ellipse fit is a 3-parameter LSQ and needs all three,
  so `n_directions = len({(angle-1) % 3 for angle in present})` must equal 3
  before a fit is attempted (`fittable`, in `MultiAnglePreview`). This is
  necessary but not sufficient: even a fittable part can end up with 0 valid
  columns if, after cross-angle alignment, no single column has all three
  directions present at once — the GUI checks `fit.valid.any()` separately
  and shows a distinct error for that case.
- **Split-half uncertainty needs all six angles, not just three directions**
  — `area_err` is `|A(a1,a2,a3) - A(a4,a5,a6)| / 2`, two independent
  3-direction fits compared against each other. With fewer than six angles
  present the GUI still runs the main fit and shows the area, axis ratio, and
  orientation, but the uncertainty band is unavailable (rendered as `±—`) and
  a caption explains why, rather than silently omitting the number or fitting
  it from an incomplete pair.
- **Batch is condition-scoped and always passes `--condition` explicitly** —
  `run_xsection`'s internal grouping keys results only on fiber number, so a
  folder holding more than one condition (e.g. C1 and C2 side by side) would
  silently merge their fibers under the same key without an explicit filter.
  The batch card always resolves to one concrete condition string before
  running (the button is disabled otherwise) and processes every fibre/part
  under that condition, not just the currently previewed one — a full run is
  scoped to "the whole condition", explained in the card's own caption
  ("roughly 20 minutes for 450 images").
- **Batch skips stage-2 replicate registration** — `run_batch` gained an
  `aggregate: bool = True` parameter; the multi-angle batch calls it with
  `aggregate=False` because stage 2 (grouping by replicate number into a mean
  ± std curve) is meaningless for angle images, where `run_xsection` is
  itself the aggregation stage. Existing replicate-mode call sites are
  unaffected (they don't pass the new argument).
- **Per-tab manual editor gated behind a checkbox** — `st.tabs` renders all
  six tab bodies on every Streamlit rerun regardless of which tab is visible,
  so an always-open editor in each of the six angle tabs would mean six extra
  full-resolution overlay renders plus six click-capture iframes on every
  interaction. Each tab's editor is instead behind its own "Edit boundaries"
  checkbox, off by default.

## Algorithm details
1. **Loading** (`_load_multiangle`): folder or upload source is parsed with
   `parse_multiangle_name` into condition/fibre/part groups; scale-bar name
   twins (the `...part1s.tiff` convention) and non-matching filenames are
   counted and reported but excluded. The ≤6 chosen angle images are run
   through the same per-image compute cache the replicate mode uses
   (`_cached_compute_path` / `_cached_compute_upload`), so preview numbers
   match the CLI's stage-1 output exactly.
2. **QC split** (`_profiles_from_results`): each angle's `MeasureResult` is
   passed through the same `anomaly.exclusion_reason` policy the CLI
   aggregator uses; angles that fail QC are excluded from the fit and shown
   with a friendly reason string instead.
3. **Preview fit** (`multiangle_preview` → `MultiAnglePreview`): included
   angles' profiles go through `xsection.build_part_stack` (cross-angle
   correlation alignment) → `fit_ellipse_projections` (per-column w²-space
   ellipse LSQ) → `hexagon_area`/`hexagon_area_expected` → `split_half_area`
   → `pair_differences` — the same call chain `run_xsection` uses per part.
   Summary scalars (`med`) are medians (area, axis lengths, ratio, rms) or
   circular medians (orientation φ) over the fit's valid columns; `hex_ratio`
   and `dw_frac` are medians over their own finite values, independent of the
   ellipse-fit validity mask.
4. **Rendering**: card 01 shows each angle's overlay, µm metrics (using the
   manual scale, never `ppu`), and its own diameter-vs-position plot. Card 02
   shows the aligned width stack, the ellipse metrics with the split-half
   band when all six angles are present, a drawn cross-section (equal-aspect
   ellipse + equal-area circle reference), an area-vs-position plot, and a
   per-angle QC table (shift px, correlation peak, uncertain/saturated
   flags, median signed residual). Card 03's batch button calls
   `run_multiangle_batch`: an in-process `run_batch(..., aggregate=False)`
   measure pass (progress bar), followed by
   `run_xsection.main(["--scale-source", str(um_per_px), "--condition",
   condition, ...])`, then reads `summary/xsection_summary.csv` back into a
   dataframe for display and download.
5. **Angle geometry** — the three projection directions (`a1/a4=0°`,
   `a2/a5=60°`, `a3/a6=120°`) come from `xsection.NOMINAL_ANGLES_DEG`, the
   same constant the CLI's ellipse fit uses; the GUI assumes it without
   re-deriving it.

## Caveats
- **Batch recomputes from disk and ignores manual boundary edits** — any
  correction made in card 01's per-tab editor only affects the live preview;
  clicking the batch button in card 03 re-measures every image from disk
  through `run_measure`'s worker path, which knows nothing about the edits.
- **Batch is condition-scoped** — a folder with multiple conditions requires
  one batch run per condition; the button always passes `--condition`
  explicitly rather than ever measuring "everything".
- **The non-numeric scale path still hardcodes `.tiff`** — `_part_scale`'s
  sidecar-reading branch (used for `--scale-source camera`/`items`, not the
  GUI's manual-scale path) builds each image's expected filename as
  `f"{cond}_{fiber:02d}_a{a}_part{part}.tiff"`; a C1 dataset with a different
  image extension would fail sidecar lookup on that path. The GUI's own
  numeric-scale path is unaffected since it never reads sidecars.
- **60°-spacing is an unverified working assumption** — inherited from the
  underlying `xsection` module (`NOMINAL_ANGLES_DEG`); if the true angle step
  is not exactly 60°, both the CLI and this GUI mode need `theta_deg`
  corrected and re-run, and the GUI has no way to detect a wrong spacing on
  its own.
- **Memory/latency for six full-resolution tabs** — C1 images are 2560×1920;
  `st.tabs` renders all six tab bodies (overlay + profile plot each) on every
  interaction even though only one is visible at a time, so switching a
  parameter can be noticeably slower here than in replicate mode with fewer
  images loaded at once. The per-tab editor checkbox (see Design choices)
  mitigates the worst case but does not remove the base six-image cost.
- **`ppu` is visible but unused in this mode** — the parameter form is
  shared with replicate mode, so the `ppu` calibration field still renders
  in the sidebar; it has no effect on anything shown in multi-angle mode,
  which is documented in the field's own help text but easy to miss.
- **"median area" (preview) vs "A_mean_um2" (batch CSV) are different
  statistics** — card 02's headline area metric is `nanmedian` over the
  fitted columns; the batch's `xsection_summary.csv` reports `A_mean_um2`,
  the mean over the same columns. The two are usually close but not
  identical, and the preview metric's own help text says so explicitly.
- **A fittable part can still yield 0 valid columns** — having all three
  projection directions present (`fittable=True`) only guarantees the fit is
  *attempted*; after cross-angle alignment, no single x-position may retain
  all three directions at once, in which case the GUI shows a distinct "no
  column has all three projection directions" error rather than a normal
  result.
