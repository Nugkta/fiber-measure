# Multi-angle cross-section measurement — design spec

Date: 2026-08-15 · Study: `02_multiangle_xsection` · Status: approved in chat, pre-implementation
Owner decisions recorded from the deep-planning session of 2026-08-15.

## Motivation

`fibrecv` currently feeds `tensile.py` a single area `A = π(d/2)²` from the group-mean diameter
(`tensile.py:296`) — a circular-cross-section assumption. Young's modulus and toughness are both
∝ 1/A, so any area error propagates 1:1 into the headline mechanical properties. The new C1
dataset images every fiber from six rotation angles precisely so that the cross-section can be
measured instead of assumed. This spec defines the method that turns the six projections into
per-position cross-section areas and per-fiber summaries.

## Data

`/Users/stan/Documents/UOM/spins/multiangle/C1/`: 15 fibers × 6 angles × 5 parts; per shot a
plain TIFF (2560×1920 RGB), an "s" twin with burned-in scale bar, and a Zeiss ZEN XML sidecar.
Established (labbook 02, 2026-08-15): angles cover a full revolution in 60° steps; a4/a5/a6 are
the 180° flips of a1/a2/a3 → each position has 3 independent directions × 2 repeats. Parts are
different positions along the fiber (non-overlapping windows; each ≈ 0.56 mm of a 10 mm gauge,
~30% coverage). Scale from metadata: CameraPixelDistance 2.2 µm / 10× → 0.22 µm/px
(≈ 4.545 px/µm); the hardcoded `ppu = 1.3680` does NOT apply to C1. No tensile data for C1 yet.

Known filename defects to fix by rename: `C1_13_a6_part34.tiff_metadata.xml` → `part4` metadata;
`C1_09_a6_par5[s].tiff` → `part5[s]`.

## Method

1. **Width profiles** — reuse `compute_measurement` (desaturation z-map → per-column sub-pixel
   edges → QC) per image to get w(x) per (fiber, angle, part). New C1 name parsing supplies the
   (condition, fiber, angle, part) key.
2. **Scale** — read µm/px per image from the Zeiss sidecar (CameraPixelDistance /
   TotalMagnification); one-off cross-validation against the "s" scale-bar images.
3. **Angle registration** — within a (fiber, part), align the six w(x) profiles by
   cross-correlation (stage repositioning shifts; prototype showed ±≲400 px works).
4. **Per-position cross-section** — at each aligned x, fit the 3-parameter ellipse projection
   w(θ) = √(a²cos²(θ+φ) + b²sin²(θ+φ)) to the six widths at nominal θ = 0…300°;
   A(x) = πab/4, plus axis ratio and orientation. The 2 repeats per direction give a per-position
   uncertainty. QC companion: area of the circumscribed hexagon from the 3 directions' support
   lines — a model-free upper bound; for a true ellipse A_ellipse/A_hex ≈ 0.91, systematic
   deviation flags non-elliptical sections. Silhouettes cannot see concavities; all results are
   convex-hull areas (stated limitation).
5. **Per-fiber summaries** — from pooled A(x) over all parts: A_mean (arithmetic; toughness),
   A_harm (harmonic; Young's modulus — series compliance weights thin segments), A_min (+
   position; tensile strength / weakest link). Headline stays A_mean for continuity. The ratio
   A_min/A_mean doubles as a thickness-uniformity index.
6. **Angle-error handling (staged)** — v1 fits with nominal angles and reports per-angle residual
   diagnostics (free by-product). Only if residuals show systematic per-angle offsets (motivated
   by ~10% width mismatch inside some 180° pairs) does v2 add per-fiber angle offsets (5 free
   params, constrained by thousands of positions) to the joint fit.

## Validation (fixed before running)

- **Synthetic**: images of known elliptical cross-sections through the real pipeline; fitted
  a, b, φ, A must be unbiased within stated tolerance. Proves the implementation.
- **Real data, primary criterion**: leave-one-direction-out — fit ellipse on 2 of 3 directions,
  predict the held-out width; compare prediction error against the circular model (mean of used
  widths), paired over 15 fibers. Ellipse must beat circle significantly. Noise floor: 180°-pair
  repeatability. Proves the model.
- Report A_ellipse vs π(d/2)² differences descriptively (effect size, no pass line).
- Negative result (ellipse ⊁ circle) still gets written up; if batch ellipticity is low
  (axis ratio < 1.05), stratify by ellipticity before concluding.

## Build

All development in git worktree `worktree-multiangle-xsection` (study-01 pattern); merge decision
at study close. Touch points identified: `io_utils.parse_name` (fails on C1 names today;
aggregate silently skips), new scale reader, new xsection fit module, two-level grouping
(fiber → angle → part) distinct from the current replicate model, `tensile.py` join deferred.
Tests: parser unit tests, fitter unit tests, synthetic end-to-end assertion of unbiasedness.
Done = validation results written into `docs/report/02_multiangle_xsection.md` + feature doc +
labbook/metalabbook updates.

## Decision record

| Decision | What was decided | Why | Rejected, and why | What would overturn it |
|---|---|---|---|---|
| Deliverable endpoint | Method + per-fiber A only; tensile join deferred | C1 tensile not measured yet | Full chain now — no data to join | C1 tensile data arriving |
| Cross-section model | Ellipse fit primary + circumscribed-hexagon upper bound as QC | 3 directions determine 3 params; model-free guard needed absent ground truth; bounded numbers | Ellipse-only (no guard); hexagon-only (coarse, correction itself shape-dependent) | A_ellipse/A_hex systematically ≠ ~0.91 → sections non-elliptical → promote corrected hull |
| Area summaries | Export A_mean, A_harm, A_min(+pos); headline A_mean | Each mechanical property needs a different physical aggregate (E→harm, strength→min, toughness→mean) | Single-metric options — strength biased on necked fibers | Community/collaborator convention mandating one metric for comparability |
| Rotation-angle error | v1 nominal 0–300° + per-angle residual diagnostics; v2 per-fiber offset fitting only if flagged | ~10% width mismatch in some 180° pairs suggests manual rotation error; diagnostics are free | Free angles from the start — risks absorbing focus bias into angles | Diagnostics showing pair mismatch is noise-level (drop v2) or clearly focus-driven (model bias, not angle) |
| Validation | Synthetic (implementation) + real-data leave-one-direction-out vs circle (model, primary) + 180°-pair noise floor | Study-01 lesson: synthetic success ≠ real-data success | Synthetic-only; real-only | Low-ellipticity batch making the test insensitive → stratify by axis ratio |
| Code location | Everything in worktree `worktree-multiangle-xsection` | Owner decision; study-01 precedent, main stays frozen | Split fixes-to-main / method-to-study (my proposal) | Merge split at study close |
| Study registration | `02_multiangle_xsection` opened 2026-08-15 with angle result recorded as pending confirmation | Angle identification is real progress; labbook records it as awaiting collaborator | Waiting for confirmation first | Collaborator contradicting the 60° step → revise labbook + refit with actual angles |
| Scale source (self-decided) | Zeiss XML metadata authoritative; "s" images one-off cross-check only | Metadata is exact; scale-bar OCR is a workaround | Scale-bar reading as primary | Metadata vs scale-bar disagreement in the cross-check |
| Filename typos (self-decided) | Rename `part34` XML and `par5` pair in place | Two isolated typos, unambiguous targets | Special-casing in parser — noise forever | — |
