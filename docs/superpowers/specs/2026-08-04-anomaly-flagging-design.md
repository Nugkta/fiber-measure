# Anomaly Flagging — Design Spec (feature 04_anomaly_flags)

Date: 2026-08-04
Status: approved (brainstormed with user)

## Problem

Measurements sometimes come out visibly wrong because of problems in the photos themselves
(bad focus, lighting, fibre placement): detected edges jump suddenly, long stretches of the
profile are missing, or the diameter shifts to a new level mid-image. Study 01 showed the
dominant error source is the photos/dataset, not the edge algorithm. Existing QC only checks
pointwise centre deviation, local rolling-MAD outliers, coverage, and band mismatch — none of
these image-level failure modes are surfaced. This feature flags anomalous measurements in
both the GUI and the output files so bad photos can be spotted and re-taken.

Naming convention: in image `753`, `7` = condition, `5` = sample #5 under that condition,
`3` = position #3 on that sample. "Replicates" are different positions on the same physical
sample — diameter genuinely varies along a fibre, hence the generous, never-excluding
replicate-outlier threshold.

## Decisions

1. **Four detectors**:
   - `edge_jump` — detrended jump between consecutive *valid* columns on either edge
     (residual = Δy − slope·Δx, refit centerline slope), threshold `jump_thresh_px`.
   - `large_gap` — longest invalid run within the fibre span as a fraction of span length,
     threshold `gap_frac`.
   - `diameter_step` — two-window (`step_window_px`) median scan over the smoothed diameter;
     max |Δmedian| / global nanmedian, threshold `step_frac`. Span < 2·window → skipped (NaN).
   - `replicate_outlier` — group-level: |median − group-median| / group-median >
     `rep_dev_frac`; needs ≥3 finite medians; NEVER excludes.
2. **Advisory by default**: new `anomaly_exclude` config switch (default False). When on,
   image-level anomalies exclude a replicate from registration (like band_mismatch), in both
   CLI and GUI via one shared `exclusion_reason` helper. Priority:
   band_mismatch → coverage → anomaly (only when `anomaly_exclude`) → None.
3. **Outputs**: `"anomaly"` sub-dict in per-image meta JSON + new
   `summary/per_image_summary.csv` (every image, excluded or not). profile CSV and
   master_summary.csv unchanged.
4. **Architecture**: pure module `src/fibrecv/anomaly.py` (numpy + config only); per-image
   detection called at the end of `qc.run_qc` (`QCResult` gains an `anomaly` field);
   replicate_outlier computed at aggregation time, shared by CLI and GUI. No bits are added
   to `EdgeResult.flags` (those gate validity).
5. **GUI**: reuse ok/warn flags badge; `⚠ ` tab-label prefix; `anomalies` column in the
   group per-image stats table. No profile-plot/overlay markers (deselected).
6. **Params**: `jump_thresh_px=10.0`, `gap_frac=0.10`, `step_frac=0.05`, `rep_dev_frac=0.25`
   exposed in the GUI param panel plus an `anomaly_exclude` checkbox; `step_window_px=100`
   config-only.

## Compatibility

- Old output trees have metas without `"anomaly"` — always read via
  `meta.get("anomaly") or {}`.
- `manual_edit.apply_manual_edits` rebuilds meta explicitly — must refresh `"anomaly"` there.
- `run_aggregate` gains `--anomaly-exclude` (BooleanOptionalAction, default None) and
  `--rep-dev-frac`; the GUI threads them at all three export call sites.
