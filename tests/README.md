# tests/

Once my folder has changes, please update me.

Pytest suite for `fibrecv` — pure-array unit tests plus rendered-image
end-to-end tests; no external data needed (the one real-sidecar check is
skipif-guarded). Run with `uv run pytest -q` from the repo root.

## Files

- `test_io_utils.py` — current — MasP2 name parsing/discovery, natural sort,
  `load_rgb` (LZW pillow fallback + double-failure diagnostics).
- `test_multiangle_names.py` — current — C1 multi-angle name parsing,
  (group, replicate) adapter, directory discovery.
- `test_scale.py` — current — Zeiss XML scale reader: both camera-tag
  spellings, magnification product, X≠Y and missing-node errors, resolver.
- `test_features_modes.py` — current — feature-mode dispatch: `"desat"`
  bit-identity regression, `"bright"` end-to-end coverage (bright + dim
  fixtures), invalid-mode error.
- `test_edges_tilt.py` / `test_edges_wall.py` — current — per-column edge
  detection and tilt-corrected diameter on synthetic fibres.
- `test_measure_edges_csv.py` — current — per-image CSV contract: additive
  `y_top_px`/`y_bot_px` columns, existing schema untouched.
- `test_register_resample.py` — current — `resample_to_grid`: grid/NaN
  policy, all-NaN degrade, interior-gap bridging (legacy) vs `max_gap`
  masking (xsection).
- `test_aggregate_summary.py` — current — stage-2 replicate registration and
  master summary.
- `test_xsection.py` — current — pure cross-section math: w²-space ellipse
  fit, hexagon closed form vs brute-force clip, split-half, anisotropy +
  φ-transfer predictors (incl. bisector rank guard), `build_part_stack`
  (shift recovery, saturation flag, interior gaps, config wiring).
- `test_xsection_pipeline.py` — current — `run_xsection` CLI on a fabricated
  known-ellipse tree: schemas, summaries, missing-angle/excluded-image/
  bad-scale/garbled-input paths, `--fibers` filter, validation writer.
- `test_xsection_synthetic.py` — current — rendered multi-angle images
  end-to-end: edge-bias measurement, ellipse/φ/area recovery pins, shift
  recovery, hexagon QC, area_err asymmetry response.
- `test_anomaly.py` — current — advisory anomaly flags + exclusion switch.
- `test_overlay.py` — current — overlay rendering smoke.
- `test_manual_edit.py` — current — manual boundary-edit round-trip.
- `test_gui_smoke.py` — current — Streamlit GUI smoke (expander widget path).
- `test_gui_multiangle.py` — current — pure `gui_app.py` multi-angle helpers:
  `_profile_from_mr`/`_profiles_from_results` (absolute-x span slicing, QC
  split), `multiangle_preview` (known-ellipse recovery, missing-angle
  `n_uncertain` rule, non-fittable direction case), `run_batch(aggregate=)`
  and headless `run_multiangle_batch`. Module-level `_bright_fibre`/
  `multiangle_folder` fixtures are reused by a later task's AppTest cases.
- `test_tensile.py` — current — stress–strain analysis stage.
