# Multi-angle cross-section measurement (third pipeline stage)
Added: 2026-08-15 · Revised: 2026-08-16 (post-review fixes: saturation flag,
φ-transfer rank guard, interior-gap masking, rms confidence gate,
`xsec_max_shift` 800→2000)
Code: `src/fibrecv/xsection.py`, `src/fibrecv/run_xsection.py`, `src/fibrecv/scale.py`,
`src/fibrecv/io_utils.py` (multi-angle parsing + LZW fallback), `src/fibrecv/features.py`
(bright feature mode), `src/fibrecv/measure.py` (edge columns), `src/fibrecv/register.py`
(`resample_to_grid`)

## What it does
Turns the six-rotation-angle image sets of the C1 protocol (15 fibers × 6 angles × 5 parts)
into per-position cross-sections instead of assuming circularity. For every (fiber, part) it
aligns the six per-angle width profiles onto a common axis, fits a per-column ellipse to the
projection widths, converts to µm via the Zeiss XML sidecar scale, and exports per-fiber area
statistics (`A_mean`/`A_harm`/`A_min`) that a tensile
join can consume in place of `π(d/2)²`. Run as
`uv run python -m fibrecv.run_xsection --out <run> --data-root <C1 dir> --scale-source items
[--validation]` after a `run_measure --feature-mode bright` pass.

## Design choices
- **w²-space linear fit** — the projection width of an ellipse obeys
  `w²(θ) = c0 + c1·cos2θ + c2·sin2θ`, exactly linear in 3 parameters, so a per-column LSQ has
  no iteration and no near-circular degeneracy: `a²=(c0+R)/4`, `b²=(c0−R)/4`, `R=√(c1²+c2²)`
  guarantees a ≥ b, and φ degrades to NaN (not a wild angle) when `R < 1e-6·c0`. A nonlinear
  (a, b, φ) fit was rejected for exactly those degeneracies.
- **Column validity, never interpolation** — a column is fitted only if all 3 projection
  directions have ≥1 finite width AND the implied b² is strictly positive; negative b² marks
  the column invalid rather than clamping, so impossible width triples are counted, not hidden.
- **Split-half uncertainty** — `area_err = |A₁₂₃ − A₄₅₆|/2` from two complete 3-direction
  exact fits; conservative by construction because it absorbs the 180° focus asymmetry (the
  dominant real error, ~8–17 px RMSE vs 4.6 px pair noise floor on C1).
- **A_min on an 11-px rolling median** — the weakest-link statistic cannot be set by a single
  bad column.
- **Anisotropy test instead of literal leave-one-direction-out** — dropping one direction
  leaves a rank-2 design for 3 parameters (θ and θ+180° share a design row), and in
  leave-one-repeat-out the LSQ passes exactly through any single-observation direction, so the
  "ellipse prediction" of a dropped width is identically its 180° partner. The honest primary
  criterion is therefore the leak-free anisotropy test (partner width vs mean-of-other-five),
  with the ellipse form judged separately by the RMS residual and a φ-transfer test. The
  φ-transfer two-direction solve is rank 2 only when the kept directions differ in
  `z = cos 2(θ−φ̂)` — when φ̂ (near-)bisects the kept pair both z values coincide, so such
  columns are skipped under `_MIN_Z_SPREAD = 0.2` (~±3.3° band per bisector) instead of
  emitting a degenerate minimum-norm answer, and `phi_hat` is NaN when no odd column solves.
- **Cross-angle alignment gates `xsec_min_corr=0.2` / `xsec_max_shift=2000`** — recalibrated
  on the full C1 set (2026-08-16): real inter-angle stage repositioning spans nearly the whole
  2560 px frame (|lag| p95 = 1741 px, visually verified), and ~2% of links are genuinely
  unalignable and fall back to zero shift + `uncertain`. A correlation peak ON the ±bound is
  flagged `saturated` and treated as uncertain — the clamped lag carries no evidence (the
  pre-review code accepted such lags as confident, corrupting 23/75 parts at the old 800 px
  bound). The aggregate-stage `min_corr`/`max_shift` are untouched; the xsection values are
  wired in via `dataclasses.replace` inside `build_part_stack`.
- **Interior dropouts stay NaN** — `resample_to_grid(..., max_gap=1.5)` masks grid points
  inside any gap between consecutive finite samples wider than 1.5 px, so a dropout inside one
  angle's profile is honest missing data instead of interpolated fabrication counted in
  `n_angles` (stage 2 keeps the legacy bridging via the `None` default).
- **rms confidence gate `xsec_rms_flag_px=6`** — a part whose median per-column fit residual
  exceeds 6 px marks its fiber `low_confidence`. Calibrated on C1 @ bound 2000 (p50 = 3.7 px,
  p90 = 6.6): every part above ~6.4 px shows one angle image whose widths sit 5–15 µm off the
  other five (verified visually) — the failure mode `rms_resid` was computed for but never
  gated on. `valid_frac` counts only *measurable* columns (all 3 directions present), so
  well-aligned parts are not penalised for the union grid widening with large shifts.
- **Scale from the XML sidecar, never `cfg.ppu`** — two conflicting fields exist
  (`CameraPixelDistance`/magnification → 0.220 µm/px vs `Scaling/Items` → 0.38892 µm/px,
  1.77×). Burned-in 100 µm scale bars adjudicated **items** (measured 0.3876 µm/px, −0.34%);
  both candidates are recorded in every shifts JSON and the six sidecars of a part must agree
  to 0.1% or the run aborts.

## Algorithm details
1. **Per-image measurement** (`run_measure --feature-mode bright`): C1 fibers are bright on a
   dark background with no saturation contrast, so `features.rgb_to_desaturation` dispatches on
   `CONFIG.feature_mode`: `"bright"` builds the z-map from `V = max(R,G,B)` with the same
   margin-row median/MAD machinery (`D = (V−v_bg)/(1.4826·MAD)`); everything downstream
   (band/edges/qc/tilt) is feature-agnostic. Per-column `y_top_px`/`y_bot_px` are persisted as
   additive CSV columns for focus/asymmetry diagnostics.
2. **Stacking** (`xsection.build_part_stack`): each angle's span-restricted profile is
   NaN-padded onto the shared absolute x grid *before* correlation (differing per-image x0
   would otherwise inject a silent misalignment), then `register.estimate_shift` aligns it to
   the lowest present angle under the xsec gates — boundary-saturated peaks → zero shift +
   `saturated`+`uncertain`; `register.resample_to_grid(..., max_gap=1.5)` produces the
   `W (6, N)` stack with interior dropouts kept NaN. Missing/excluded angles become NaN rows.
3. **Fit** (`fit_ellipse_projections`): per-column 3-parameter LSQ in w²-space, vectorised by
   grouping columns on their 2⁶ finite-mask patterns; residuals and per-column validity as
   above. Angles are `NOMINAL_ANGLES_DEG = 0..300°` in 60° steps — a parameter everywhere, so
   a corrected spacing only requires re-running.
4. **Outputs** (`run_xsection`): per-part CSV (aligned widths + a/b/φ/area columns), plot,
   shifts JSON (incl. per-link `saturated` flags + `n_saturated`); per-fiber
   `xsection_summary.csv` (A_mean/A_harm/A_min(+part,+x), axis-ratio stats, φ circular mean,
   pair Δw fraction, A_circle=π(d̄/2)² comparison, `n_uncertain_shifts`,
   `n_saturated_shifts`, `part_rms_med_max_px`, `low_confidence`);
   `xsection_angle_residuals.csv` (v2-gate evidence); `xsection_validation.csv` with
   `--validation` (anisotropy + φ-transfer RMSEs). Garbled per-image CSVs/meta JSONs are
   skipped with a `[WARN]` (image treated as absent); a truncated sidecar aborts rc=2.

## Caveats
- **Silhouettes see the convex hull only** — a concave (e.g. kidney-shaped) section is
  invisible to any projection method; the ellipse fit assumes convexity.
- **Nominal 60° spacing is a working assumption** (180° pairing is proven; collaborator
  confirmation pending). If the true step differs, re-run with corrected `theta_deg`.
- **Ellipse form is unverifiable from 3 directions** — the fit is saturated, so held-out
  widths cannot test ellipticity directly; the RMS residual and φ-transfer are the
  model-form evidence. On C1 the φ-transfer fails because φ genuinely varies along fibers
  (twist), so per-part fixed-φ predictions are poor — a protocol limit, not a fit bug.
- **Focus bias entangles with ellipticity** — per-image focus states shift widths by ~2–10 px;
  the split-half `area_err` absorbs it conservatively, and within-part φ coherence can be
  focus-induced (cross-part coherence is the artifact-immune evidence).
- **Edge bias δ** — the fixed edge_z level sits on the smoothed shoulder, biasing each wall
  outward by an amplitude-dependent ~2–4 px/side (δ=3.9 px at the C1 bright end in synthetic
  replay). Ratios and φ are unaffected to first order; absolute areas inherit ~2δ per axis.
- **A_min is still an extreme statistic** — the rolling median guards single columns and,
  after the alignment fix, end-clustering is gone (0/15 fibers place A_min near a span end);
  but A_min remains the output most sensitive to alignment and single-image quality — consume
  it together with `low_confidence`/`n_saturated_shifts`.
- **Six C1 fibers are low-confidence** (01, 07, 08, 09, 13, 15): each has at least one angle
  image whose widths sit 5–15 µm off the other five (bad single shot — focus/edge state).
  The rms gate flags them; their A_min/uniformity should not be quoted without the flag.
