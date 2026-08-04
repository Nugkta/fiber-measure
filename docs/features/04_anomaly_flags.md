# Anomaly flagging
Added: 2026-08-04
Code: `src/fibrecv/anomaly.py` (new), `src/fibrecv/qc.py`, `src/fibrecv/compute.py`,
`src/fibrecv/manual_edit.py`, `src/fibrecv/run_aggregate.py`, `src/fibrecv/gui_app.py`,
`src/fibrecv/config.py`

## What it does
Measurements sometimes come out visibly wrong because of problems in the photos
themselves (bad focus, lighting, fibre placement) — study 01 showed the dataset,
not the edge algorithm, is the dominant error source. This feature detects three
image-level symptoms (a detected edge that jumps suddenly, a long unmeasurable
stretch, a mid-image shift of the diameter level) plus one group-level symptom
(a replicate whose median diameter sits far from its group), and surfaces them
everywhere a user looks: a `⚠` prefix on the replicate's GUI tab, the anomaly
names in the amber flags badge and in a new `anomalies` column of the per-image
stats table, an `"anomaly"` sub-dict in every per-image meta JSON, and a new
`summary/per_image_summary.csv` with one row per image (excluded ones included).
Flags are advisory by default; an `anomaly_exclude` switch (GUI checkbox /
`--anomaly-exclude` CLI flag) makes image-level anomalies drop the replicate
from registration, exactly like `band_mismatch` does.

## Design choices
- **Advisory by default** — a flag marks a photo for review/retake; it only
  affects the numbers when the user opts in via `anomaly_exclude`. Rejected:
  auto-excluding, because study-01 style datasets would silently lose most
  replicates.
- **`replicate_outlier` never excludes** — replicates are different positions
  along the same physical sample (naming `7_5_3` = condition 7, sample 5,
  position 3) and diameter genuinely varies along a fibre, so a generous 25%
  deviation threshold flags it for attention only.
- **Pure detector module** (`anomaly.py`, numpy + config only) called from the
  end of `run_qc` — the GUI's manual-edit path re-runs QC, so corrected
  boundaries clear their flags with no extra wiring. Rejected: a separate
  post-processing pass, which manual edits would have bypassed.
- **No new bits in `EdgeResult.flags`** — those gate per-column validity
  (`flags == FLAG_OK`); anomaly evidence lives in its own `AnomalyResult` so
  detection numbers are untouched.
- **One shared `exclusion_reason` policy** used by both `run_aggregate` and the
  GUI group panel, replacing two previously-duplicated band_mismatch/coverage
  checks — the CLI and GUI can no longer disagree about who is excluded.
- **Existing outputs unchanged** — profile CSVs and `master_summary.csv` keep
  their schemas; new data goes to the meta JSON and the new
  `per_image_summary.csv`. Old output trees (metas without `"anomaly"`) load
  fine. Rejected: profile-plot overlay markers (deselected during design).

## Algorithm details
Per image, `detect_image_anomalies(y_top, y_bot, diameter_smooth, valid, x0,
x1, slope, cfg)` runs at the end of `qc.run_qc` and attaches an
`AnomalyResult` to `QCResult` (serialised by `as_dict()` into the meta JSON;
`jump_cols` capped at 20 entries, NaN → null):

- **edge_jump** — for each edge separately, diff the edge row over consecutive
  *valid* columns and detrend with the refit centerline slope
  (`residual = Δy − slope·Δx`, the same Theil-Sen slope QC fits for the
  centre-deviation check; `band.slope` fallback when < 10 valid centres).
  `|residual| > jump_thresh_px` (default 10) flags the later column of the
  pair. Diffing valid neighbours means a gap cannot hide a jump; detrending
  keeps a tilted fibre from false-flagging across that gap.
- **large_gap** — longest run of invalid columns inside the band span, as a
  fraction of span length; `> gap_frac` (default 0.10) flags. All-invalid
  gives `longest_gap_frac = 1.0`.
- **diameter_step** — two adjacent `step_window_px` (default 100) windows
  slide along `diameter_smooth` (stride `window/25`); the shift is
  `|median(right) − median(left)| / global median`. The maximum forms a
  plateau around a genuine level shift, so the reported `step_col` is the
  plateau centre; `> step_frac` (default 0.05) flags. Spans shorter than two
  windows skip the check (`step_frac = NaN`). A linear taper moves both
  window medians together and stays below threshold.
- **replicate_outlier** — at aggregation (`run_aggregate.main`) and in the GUI
  (`main()`, after manual edits), `detect_replicate_outliers` compares each
  replicate's `median_diameter_um` to the group median; deviation
  `> rep_dev_frac` (default 0.25) flags. Needs ≥ 3 finite medians; excluded
  replicates still participate in the group median.
- **Exclusion policy** — `exclusion_reason(band_mismatch, coverage, flags,
  cfg)` returns the first matching reason with priority band_mismatch →
  coverage → anomaly (the last only when `cfg.anomaly_exclude`), or `None` to
  keep the replicate.

`summary/per_image_summary.csv` columns: `name, group, replicate,
median_diameter_um, coverage, anomaly_flags` (`;`-joined), `max_jump_px,
longest_gap_frac, step_frac, rep_dev_frac, excluded, excluded_reason`.
It is deliberately a tree-wide audit: it covers every parseable image found,
even under a `--groups` selection (where `master_summary.csv` covers only the
selected groups), and it is written even when exclusions leave nothing to
register — precisely then it is the record of why.

## Caveats
- Thresholds were chosen by design reasoning, not calibrated on the full
  dataset; expect to tune `jump_thresh_px`/`step_frac` once real bad photos
  are run through them.
- `diameter_step` needs a span of at least `2·step_window_px` (200 px default)
  — short fibres skip the check silently (`step_frac` = null in the meta).
- `replicate_outlier` is silent for groups with fewer than 3 usable medians,
  and its deviation is computed from `median_diameter_um` in the meta — old
  trees without that field skip it.
- A genuine sharp physical feature (e.g. a bead on the fibre) can legitimately
  trigger `edge_jump`/`diameter_step`; the flags mark "look at this photo",
  not "this photo is wrong".
- `edge_jump` compares consecutive valid columns; a jump that coincides
  exactly with the start of a long invalid gap on a *curved* fibre can be
  under- or over-detrended since a single global slope is used.
