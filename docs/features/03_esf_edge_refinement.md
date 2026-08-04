# Erf edge refinement
Added: 2026-08-03
Code: `src/fibrecv/refine.py` (new stage), `src/fibrecv/config.py`
(`refine_*` CONFIG fields), `src/fibrecv/compute.py` (pipeline wiring +
`meta["refine"]`), `src/fibrecv/run_measure.py` (`--refine`/`--no-refine`
CLI flag), `src/fibrecv/gui_app.py` (`refine_on` toggle, σ(x) profile-plot
trace); tests in `tests/test_refine.py`.

## What it does

Optical blur (focus-dependent PSF, σ roughly 5–15 px) smears the true
air/silk boundary into a ramp, and the existing wall finder (`edges.py`)
places the edge at a fixed low threshold near the base of that ramp — a
threshold crossing that drifts with focus, giving an image-dependent
diameter bias of several px per side (measured +2.3 / −3.4 / +3.4 px on
three MasP2 images, see `docs/labbooks/01_esf_edge_consistency.md`). Erf
edge refinement adds a new stage between `edges` and `qc`: it refits each
detected wall, in blocks of columns, as a Gaussian-blurred step and moves
the boundary to the fitted 50% midpoint, which by the knife-edge principle
is invariant to the blur width. Where a block's fit doesn't pass its
quality gates, the column keeps its original (legacy) edge unchanged. The
stage is **off by default** (`refine_on=False`) — the full-set A/B in
`docs/report/01_esf_edge_consistency.md` found the real-image benefit
unproven, so as of the owner decision on 2026-08-04 it is opt-in via
`--refine` / the GUI checkbox pending a focus-sweep study. It is
bit-identical to the old pipeline when off, so it is both the current
shipped default and its own A/B control.

## Design choices

- **Work entirely in the perpendicular frame.** The profile coordinate `t`
  is a perpendicular distance from the legacy edge, sampled at
  `y = anchor ± t/cth` using `m, cth = band.tilt_geometry(bnd.slope)` — the
  same single tilt source `edges.py` and `qc.py` already use (`band.py:60`).
  The inward window (`inside = clip(refine_in_frac · band_half · cth, 8,
  refine_in_max)`), the fitted `t0`/`sigma`, and the gates that judge them
  are all perpendicular px, so gates apply unconverted; only the final
  applied shift is converted back to vertical (`t0/cth`). Every conversion
  is an exact identity at zero tilt. Rejected: re-deriving tilt locally per
  wall, which risked drifting from the `edges`/`qc` frame.
- **No pre-smoothing of the desaturation map `D`.** Smoothing `D` before
  fitting would inflate the fitted blur sigma and defeat the measurement.
  Noise is handled instead by fitting the block-mean profile (a block of
  `refine_block` columns averaged, ~86–127 samples along `t`) with least
  squares — the study's own erf fits on raw channels already ran at
  3.5–6.7% residual, so this was validated headroom, not a new risk.
- **`RefineResult` carries `o_top`/`o_bot`** (the applied per-column
  perpendicular offsets), beyond what the design spec's dataclass strictly
  required — needed to compute meta's `median_abs_t0_top/bot` and useful
  directly in tests.
- **Offsets apply only to anchor columns** — finite `y_top` & `y_bot` AND
  `flags == FLAG_OK` — and only within a chain's interpolation coverage;
  flagged columns never move, so refinement cannot paper over a column
  `edges`/`qc` already distrusts.
- **Bit-identical when `refine_on=False`.** `refine_edges` returns the SAME
  `EdgeResult` object (not a copy) and an empty `RefineResult`; every
  numeric array is untouched. Meta necessarily still gains a
  `refine.enabled=false` sub-dict and the new `refine_*` CONFIG keys appear
  in `params`, but nothing numeric about the measurement changes.
- **Sigma and residual are interpolated the same way as the offsets** —
  between passing block centres — so the σ(x) diagnostic trace and the
  offset field share one interpolation rule instead of two.
- **Blocks are consecutive, full-width, from `bnd.x0`.** A trailing partial
  block is skipped rather than padded; block centres are integers; nothing
  is extrapolated beyond a chain's outermost passing centre — an
  unsupported stretch of fibre is left at its legacy edge rather than
  guessed.
- **Border-clipped columns are dropped before the 70% quorum check**, not
  after — a column whose sample window would leave the image never
  contributes a replicated-border sample to a block's mean profile, and
  never counts toward the quorum either.
- **`refine_on=True` at landing** (the design spec's default); `--no-refine`
  / `refine_on=False` was the deliberate A/B control for the Task 5
  validation, not a fallback for uncertainty. **Flipped to `refine_on=False`
  on 2026-08-04** (owner decision, after Task 5 returned a null result on
  real data) — `--refine` / `refine_on=True` is now the opt-in arm instead.
- **float32 throughout** (NEP-50 promotion rules), so the refined
  `diameter` matches the dtype/rounding semantics of `edges.py`'s own
  diameter computation (`edges.py:359`) — no numeric drift from a dtype
  mismatch between the legacy and refined values.
- **Residual gate at 0.15, not the 0.08 the single-image prototype
  suggested.** On the full 141-image MasP2 set, 0.08 refined only 76.4% of
  anchor columns (median image 81.9%), under the design spec's ~80% bar. The
  cause is not block noise — widening the block barely moves coverage (block
  16 → 66.9%, 32 → 70.3%, 64 → 68.0% on a 12-image probe) — but systematic
  model mismatch from specular stripes and shadow ramps. 0.15 lifts coverage
  to 89.6% (median image 95.7%). Where the tight gate already had a
  representative sample, the added blocks agree with it — `masp2 10_1_1` top
  wall, 65.3% → 71.6% coverage: median σ 11.654 → 11.680 px, median |t₀|
  5.658 → 5.245 px. Where it did not (the same image's bottom wall, 23.1% →
  99.1%), the fit population is replaced rather than extended and does move
  (σ 9.075 → 7.746 px, |t₀| 4.945 → 3.527 px) — but a gate rejecting three
  quarters of a wall was not delivering a defensible estimate there in the
  first place. The cost is bounded and measured: a fit admitted at residual
  0.08 is within ~0.7 px of the true step, one admitted at 0.15 within
  ~1.5 px — still smaller than the 2.3–3.4 px per-side bias of the legacy edge
  the block would otherwise fall back to. Rejected alternative: shrinking
  `refine_in_max`, which lifts coverage further but truncates the erf's inner
  plateau, biasing σ and the midpoint of blocks that were *already passing* —
  a change to the estimate itself, where `relmax` leaves every passing block's
  fit untouched and only adds more of them (see
  `docs/report/01_esf_edge_consistency.md` §3.1).
- **Perf budget: < 1 s added per 2560-px image** (~320 fits/image, block
  width 16) — met, see Algorithm details.
- **Manual edits run after refine and are not re-refined** — see Caveats.

## Algorithm details

Runs between `edges.detect_edges` and `qc.run_qc` inside
`compute.compute_measurement` (`compute.py:80`). For each wall (top:
`sign=+1`, bottom: `sign=-1`, `+t` always points into the fibre):

1. **Erf model** (`refine.py:_erf_model`):
   `D(t) = a + (b − a) · 0.5 · [1 + erf((t − t0) / (σ√2))]`, a
   Gaussian-blurred step of height `b − a`, midpoint `t0`, blur width `σ`.
   By the knife-edge principle `t0` is invariant to `σ`; the legacy
   fixed-level crossing (`edges.py`) is not — it drifts with focus.
2. **Blocking** (`refine.py:refine_edges`): consecutive `refine_block`-column
   (default 16) runs from `bnd.x0`; trailing partial dropped. Each block's
   profile is the mean of its columns' aligned samples (`_block_profile`),
   sampled bilinearly along `t ∈ [−refine_out, +inside]` in 0.5 px steps
   (`_T_STEP`) — `refine_out` default 35 px outward (background side),
   `inside = clip(refine_in_frac · band_half · cth, 8, refine_in_max)` px
   inward (default `refine_in_frac=0.8`, `refine_in_max=28`, floor 8) —
   ~86–127 samples per fit. The inward cap keeps the fit window short of the
   specular core inside a well-lit fibre.
3. **Quorum** (`_QUORUM=0.7`): columns whose sample window would leave the
   image are dropped first; a block is only fit if ≥ 70% of its columns
   survive.
4. **Fit** (`_fit_block`, `scipy.optimize.curve_fit`): initial guess from
   the outer/inner 6-sample means and the sample nearest the level
   midpoint, σ0=3.0, loose bounds (σ floor 0.05). A block's fit is accepted
   only when all four gates pass: `b − a > 0` (genuine rising step);
   `rms_residual / (b − a) < refine_relmax` (default 0.15, raised from 0.08
   by the study-01 M5 tuning — see Design choices);
   `refine_sigma_min ≤ σ ≤ refine_sigma_max` (default 0.8–20 px);
   `|t0| ≤ refine_maxshift` (default 12 px).
5. **Interpolation** (`_interp_side`): passing blocks' `t0`/σ/residual are
   `np.interp`-spread over the columns between their centres; a chain
   breaks where the gap between consecutive passing centres exceeds
   `refine_gap_blocks` (default 2) blocks, and nothing is extrapolated past
   a chain's outermost centre — unsupported columns are NaN and keep their
   legacy edge.
6. **Apply** (`refine_edges`): on anchor columns within coverage,
   `y_top += o_top/cth`, `y_bot -= o_bot/cth` (perpendicular offset
   converted to vertical), then `diameter = (y_bot − y_top) · cth` is
   recomputed — `qc.run_qc` trusts `EdgeResult.diameter` and never
   recomputes it itself, so skipping this step would silently discard the
   refinement.
7. **Diagnostics**: `RefineResult` carries `refined_top/bot`,
   `sigma_top/bot`, `resid_top/bot`, `o_top/bot` (all NaN off-coverage) plus
   `n_blocks`/`n_pass_top`/`n_pass_bot`; `compute.compute_measurement` folds
   these into `meta["refine"]` — `enabled`, `n_blocks`, `n_pass_top/bot`,
   `coverage_top/bot` (refined / anchor columns), `median_sigma_top/bot`,
   `median_abs_t0_top/bot` (None-safe when nothing was refined).

The GUI (`gui_app.py`) exposes `refine_on` as a sidebar checkbox and, when
`ref` carries finite sigma values, overlays a σ(x) blur-width trace on a
twin axis of the per-replicate diameter plot; the CLI (`run_measure.py`)
exposes `--refine`/`--no-refine` (`argparse.BooleanOptionalAction`).

Performance: ~320 fits per 2560-px image at `refine_block=16`; measured
0.32–0.59 s per image on `masp2 10_1_1.jpg` (< 1 s budget).

## Caveats

- **Partial coverage blends estimators — a measurement-validity caveat.**
  The reported `diameter`/`diameter_um`, and therefore the per-image median
  that reaches `qc`, is the erf-refined value *only* on columns a passing
  block chain covers; every other column keeps its unchanged legacy
  fixed-threshold edge. On a wall with low refine coverage the median is
  therefore still mostly the legacy estimator even though `refine_on=True`
  — this is the plan-mandated fallback (no half-converged fit is ever
  substituted for a full legacy edge), but it means "refinement is on"
  does not imply "the reported number is the erf estimate" for a given
  column or image. Single-image observation (`masp2 10_1_1.jpg`,
  `docs/labbooks/01_esf_edge_consistency.md`): top-wall coverage 65.3%,
  bottom-wall coverage 23.1% — the bottom wall's specular core stripe sits
  inside the inward fit window and trips the residual gate on most blocks,
  so the bottom-wall median there is still mostly the legacy edge. Check
  `meta["refine"]["coverage_top"/"coverage_bot"]` before treating a
  per-image median as "refined." Across the whole MasP2 set at the tuned
  `refine_relmax=0.15`, coverage averages 89.6% (median image 95.7%), but the
  p10 image is at 59.6% on the top wall and 79.8% on the bottom — so
  low-coverage walls do still exist and are still mostly legacy.
- **Manual edits run after refine and override it, with no re-refine.**
  `manual_edit.py` corrects `EdgeResult` boundaries downstream of
  `refine_edges`; a manually edited stretch is not re-fit as a blurred
  step, so its diameter is whatever the user's clicked points imply, not an
  erf midpoint — expected, since manual editing exists for cases where
  automatic detection (refined or not) failed.
- **Thin fibres / specular cores can push the inward window onto the wrong
  feature.** `inside` is capped at `refine_in_max` (default 28 px, floor
  8 px) specifically to keep the fit window short of a well-lit fibre's
  specular core, but on a thin fibre, or one with a bright core close to
  the wall, the window can still see the core or opposite wall and trip the
  residual or shift gates — the block then fails and falls back to the
  legacy edge rather than fitting an incorrect blurred step (the coverage
  numbers above are exactly this failure mode).
- **Interpolation seams and along-fibre ripple.** A chain break (gap >
  `refine_gap_blocks`) leaves a hard boundary between refined and legacy
  columns, and within a chain the offset field is piecewise linear between
  block centres. `qc`'s rolling-MAD outlier check does not object — the
  outlier flag count actually *fell* 14% across MasP2 — but the detrended
  along-fibre noise `std(raw − smooth)` rose from 0.681 µm to 0.846 µm
  (median image 0.342 → 0.491), improving in only 8 of 141 images. If you
  need per-column fidelity along the fibre, this stage currently costs you
  some.
- **It does not (yet) deliver replicate consistency, its stated goal.** The
  full-set A/B (`docs/report/01_esf_edge_consistency.md`) found the
  within-group between-replicate std of the per-image mean diameter falls in
  23 of 46 groups and rises in 23 (sign test p = 1.00; mean 11.14 → 10.85 µm).
  The focus-dependent bias the stage targets is real and measurable
  (+1.87 µm per px of fitted σ, within-group), but refinement removes only
  ~11% of it, and on MasP2 the replicate spread is dominated by the sample
  (the three replicates of a group differ in median diameter by 22.9% at the
  median group), not by edge placement. Treat the stage as a better-motivated
  edge estimator with a validated fit model, not as a consistency fix.
