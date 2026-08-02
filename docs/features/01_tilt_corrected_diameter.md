# Tilt-corrected (perpendicular) diameter
Added: 2026-07-31
Code: `src/fibrecv/edges.py` (`_axis_average`, `_vshift`, `detect_edges`),
`src/fibrecv/band.py` (`tilt_geometry`, `TAN_TILT_MAX`), `src/fibrecv/qc.py`
(band-mismatch ratio), `src/fibrecv/manual_edit.py` (edited-column diameter);
tests in `tests/test_edges_tilt.py`; visual explainer in
`others/tilt_corrected_diameter/` (open `algorithm_visualisation.html`).

## What it does

The diameter used to be reported as `y_bot − y_top`, the vertical distance
between the detected walls. On a fibre inclined at angle θ that is wrong twice:
the vertical chord is a slanted cut inflated by 1/cos θ (+15.5 % at 30°), and —
about twice as large — the 41-column neighbourhood average cut straight across
the tilted walls, smearing them into ramps so the boundary crossings landed too
far out on both sides (+37 % total at 30°, session-validated on synthetic
fibres of known width). The feature makes detection tilt-invariant and reports
the width perpendicular to the fibre axis. The measured value is a clean
replacement — same CSV columns, meta keys and summary tables — flat to within
2 % from 0° to 45° of tilt, and bit-for-bit identical to the old pipeline for a
horizontal fibre.

## Design choices

- **Compensate the vertical scan rather than rotate the image.** Scanning
  stays per-column and `y_top`/`y_bot` stay sub-pixel rows in original image
  coordinates, so overlays, manual-edit clicks and QC keep their meaning;
  rotation would resample the image and shift the calibrated signal statistics.
- **Tilt from the existing global Theil–Sen centerline fit, not local gradient
  orientation.** Local gradient directions are exactly what this imagery
  corrupts (specular rim, iridescence banding, shadow ramps — the documented
  failure modes), while the mask-centroid + Theil–Sen route integrates over the
  whole band and tolerates ~29 % outlier columns. The cost is a straight-fibre
  assumption whose error is second-order (≈ δ·tan θ for local deviation δ); at
  the real-data tilt range (median 0.9°, max 6.7°) that is < 0.5 %.
- **Axis-following average as shear → box filter → unshear** (three vectorised
  passes, cost independent of `wcol`, exact identity at zero slope). Rejected: a
  per-column Python loop (slow) and a 2-D oriented kernel (changes the
  calibrated smoothing behaviour).
- **One shared, clamped tilt helper.** `band.tilt_geometry` is the single
  source of (slope, cos θ) for `edges`, `qc` and `manual_edit`, so the three
  consumers cannot drift apart; |slope| is clamped to tan 60°, bounding every
  1/cos θ compensation at 2× even for a pathological band fit.
- **Clean replacement, no config switch or legacy column** (user decision).
  Every downstream consumer reads `EdgeResult.diameter` / `QCResult.diameter_raw`,
  so only the producers changed.
- **Manual edits re-estimate tilt from the edited boundaries** (Theil–Sen on
  the edited midline, `manual_edit.py`), because manual editing exists
  precisely for images where the automatic band fit failed — its stale slope
  cannot be trusted.

## Algorithm details

Tilt comes from the band stage: per-column centroids of the coarse
desaturation mask, robust Theil–Sen line fit → `band.slope`
(`band.py:centerline_fit`), then `tilt_geometry` gives the clamped slope m and
cos θ = 1/√(1+m²). Four compensations in `edges.detect_edges`, all exact
identities at θ = 0:

1. **Axis-following neighbourhood average** (`_axis_average`): shear each
   column by m·(x − x_c) (`_vshift`, bilinear in y), horizontal
   `uniform_filter1d` of width `wcol` (default 41), unshear. The window follows
   the fibre axis instead of cutting across the walls; "nearest" padding then
   extends the profile *along the axis* at image borders.
2. **Vertical smoothing rescaled**: Gaussian σ = `sigma_y`/cos θ, so smoothing
   keeps its calibrated perpendicular width.
3. **Gradient and guard rescaled**: the vertical gradient is divided by cos θ
   so the calibrated wall gates (`slope_min`, `slope_rel`, `slope_cap`,
   `rise_min`) see perpendicular-frame slopes; the `guard` row count is
   multiplied by 1/cos θ so the wall-local background zone keeps its
   perpendicular size.
4. **Perpendicular projection**: `diameter = (y_bot − y_top) · cos θ`.

Downstream: QC's band-mismatch check converts band thickness to the same
perpendicular frame (`2 · band_half · cos θ`, `qc.py`) before comparing;
manual-edit columns use `(y_bot − y_top) · cos θ` with θ re-fitted from the
edited boundaries. Validated by `tests/test_edges_tilt.py`: invariance ≤ 2 %
relative across 0–45°, identity at zero slope, exact-width even-`wcol`
behaviour, and clamping.

## Caveats

- **Straight-fibre assumption.** One global θ per image; a bent fibre gives a
  local error ≈ δ·tan θ. Negligible at the observed tilt range, but a visibly
  curved fibre at high tilt is outside the model. Check
  `band.centroid − band.c_fit` residuals if in doubt; the upgrade path is a
  smoothed per-column centerline fed to the same shear machinery.
- **Validated to 45°, clamped at 60°.** Between 45° and 60° values are still
  produced but sit outside the validated envelope; beyond 60° the clamp
  deliberately under-corrects rather than letting compensations explode.
- **Axial resolution.** The 41-column average low-passes the diameter profile
  along the fibre: necks/defects shorter than ~`wcol` px are attenuated, and
  neighbouring columns are correlated over that span (fewer independent samples
  than columns).
- **Tilt quality depends on the band mask.** A `low_confidence` band (no
  component spanning ≥ 85 % of the width) can bias the slope; non-finite slopes
  are treated as horizontal.
- **Absolute width vs invariance.** The boundary level intentionally sits
  partway up the smoothed wall (`edge_z` above the wall-local base), which on
  the synthetic fixture reads systematically wide — the validated property is
  tilt *invariance*, not absolute accuracy.
