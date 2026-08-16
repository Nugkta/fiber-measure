Hypothesis: Per-position ellipse fitting of the six 60°-spaced projections predicts held-out
direction widths better than the circular assumption, and yields cross-section areas
systematically different from π(d/2)².
Status: done
Started: 2026-08-15

## 2026-08-16
**What I did** — Close-out (Task 13). Feature doc `docs/features/03_multiangle_xsection.md`;
study report `docs/report/02_multiangle_xsection.md` (intro/method/results/discussion, synthetic
recovery table, Wilcoxon stats, six figures `docs/report/02_xsec_*.png`); `src/fibrecv/README.md`
+ root `README.md` updated to the three-stage pipeline (on branch `worktree-multiangle-xsection`);
metalabbook row → done. Full suite green on the branch (202 tests).

**What I observed** — Final deliverables: `fibrecv_output/multiangle_c1/` in the worktree holds
75 per-part CSVs/plots/shifts, `xsection_summary.csv` (15 rows), `xsection_angle_residuals.csv`,
`xsection_validation.csv`. Headline numbers in the report; hypothesis verdict split — shape half
supported (all fibers elliptical, ratio med ≈1.13, φ persists along fibers, Fisher p=1.4e-4),
mechanical half largely negated (A_mean within ±3% of π(d̄/2)²; A_min/A_mean 0.71–0.92 is the
quantity the circular pipeline actually misses).

**What I think it means** — Study closed. The practical product is the third pipeline stage +
per-fiber A_mean/A_harm/A_min; the circular area assumption is acceptable for C1-like fibers,
but uniformity (A_min) is new information. Confidence: high on the code (synthetic proof +
208 tests), medium-high on the science (focus systematics are the accuracy ceiling). Next
leverage: one-focus-protocol imaging, collaborator angle confirmation, C1 tensile join.

## 2026-08-15
**CORRECTION — scale adjudication (Task 4 hard gate): the `Scaling/Items` field is the true
scale, NOT the camera-pitch value.** The earlier reading (this labbook, planning entry below:
"CameraPixelDistance 2.2 µm, TotalMagnification 10× → 0.22 µm/px") is wrong for these TIFFs.
Burned-in 100 µm scale bars on three "s" twins spread across the dataset (C1_01_a1_part1s,
C1_08_a3_part2s, C1_15_a6_part5s; crops kept at
`fibrecv_output/multiangle_c1/scale_check/crop_*.png` in the worktree) all measure 258.0 px
cap-center-to-cap-center → 0.3876 µm/px, which is −0.34% from `Scaling/Items/Distance` =
0.38892 µm/px and +76% from camera 2.2/10 = 0.220 µm/px. Decision per spec overturn condition:
`--scale-source items` for all later runs; fibers are therefore ~72–80 µm in diameter
(185–205 px), not 41–45 µm. Sidecar note: the camera tag in these files is
`CameraPixelDistance` (singular), value "2.2,2.2"; `scale.py` accepts both spellings.

**What I did (implementation, Tasks 2–9)** — LZW pillow fallback in `load_rgb` (commit 2dc2085);
`scale.py` Zeiss reader (bf7a182; NB real tag is singular `CameraPixelDistance`, both spellings
accepted); scale-bar hard gate (see CORRECTION above); `feature_mode="bright"` V-based z-map
(7e57995, zero knob tuning needed); additive `y_top_px`/`y_bot_px` CSV columns (fc68872);
`resample_to_grid` extraction + all-NaN guard (d207e90); measurement pilot (fibers 01+08) then
full 450-image batch `run_measure --feature-mode bright --ppu 2.5712 --jobs 6 --out
fibrecv_output/multiangle_c1` (all other knobs default); `xsection.py` pure math module + 48
unit tests (a4b4740).

**What I observed (Tasks 5–9)** —
- Probe gate: C1_01_a1_part1 cov 0.91 / 212 px, C1_08_a3_part2 cov 1.00 / 117.6 px,
  C1_15_a6_part5 cov 0.96 / 166 px; overlays on true walls (fibers 08/15 are genuinely thinner
  than the 185–205 px planning range, which came from fibers 01–06 only).
- Full batch audit: 450/450 ok, 0 errors, 0 low-confidence, coverage ≥0.8 on 449/450 (only
  C1_04_a1_part1 at 0.76 — defocused left end, measured columns clean), noisiest-MAD-quartile
  median coverage 1.00 (no dim-fiber degradation). Median diameters 45.7–84.0 µm (p50 63.9) at
  0.38892 µm/px.
- Cross-angle alignment study (pilot, 180°-pairs): at `max_shift=400` 16/30 pair links pinned at
  the bound; at 800 the p25 corr rises 0.228→0.559 and pinning drops to 3/30 — real inter-angle
  stage repositioning reaches ~±800 px (~±310 µm). a1-referenced peaks at ms=800: p10 0.194,
  median 0.562. Pair-mean chaining tested and rejected (cross-link median 0.432→0.484, not worth
  deviating from the plan's ref=a1 design). Frozen: `xsec_max_shift=800` (new knob, evidence
  above), `xsec_min_corr=0.2` (zeroes the genuinely unalignable ~10% — near-flat profiles).
- Plan deviation (recorded): build_part_stack unit-test shift tolerance relaxed 0.5→1.25 px —
  the NaN-fill constant extension biases the correlation peak by up to ~1.1 px when per-angle
  spans differ; the synthetic end-to-end pin stays ±1.5 px as planned.

**What I did (Task 12, real-data run + validation)** — Full third-stage run:
`run_xsection --out fibrecv_output/multiangle_c1 --data-root .../C1 --condition C1
--scale-source items --validation` (75 parts, 15 fibers). Scratchpad analysis: per-fiber paired
Wilcoxon on the anisotropy test, σ_pair noise floor, cross-part φ coherence permutation test,
v2 angle-residual gate, shifts-on-vs-zero sensitivity (fibers 01+08), A_min end-clustering.
Six figures → `docs/report/02_xsec_*.png`.

**What I observed (Task 12)** —
- Summary: axis_ratio_median 1.094–1.245 across the 15 fibers (grand median ≈1.13); area_ratio
  (A_mean/A_circle) 0.972–1.011; uncertain shift links 27/375 (7.2%), pinned at ±800 px 29/375.
- Primary anisotropy test: directional (180°-partner) beats the isotropic 5-width mean in 9/15
  fibers; paired Wilcoxon p=0.076; median relative RMSE reduction 6.5% (−14%..+36%). Pooled
  RMSEs 7.6–17.3 px vs pair noise floor σ_pair = 4.56 px (1.77 µm) — per-image systematics
  (focus/edge bias), not sensor noise, dominate both predictors. NB the directional predictor
  carries a 1-sample-vs-5-sample variance handicap (√2σ vs √1.2σ), so beating the mean at all
  requires genuine direction dependence; implied direction modulation ~4% of width.
- φ-transfer (model-form) test: NEGATIVE — rmse_ell 35.0 px vs circle 13.6 px. Cause: φ is not
  constant — within-part φ concentration median 0.70 (p25 0.51), and part-mean φ moves along
  fibers. A per-part global φ therefore does not transfer to held-out directions.
- Cross-part φ coherence (immune to per-image focus artifacts, which re-randomise between
  parts): 5/15 fibers individually p<0.05, Fisher combined p=1.4e-4 — section orientation
  genuinely persists along fibers at the population level (impossible for circular sections).
- Hexagon QC: median measured hex_ratio 0.9060 vs expected-for-fit 0.9052 — the convex-bound
  geometry is fully consistent with the fitted ellipses.
- v1→v2 gate: CLOSED — "nominal angles suffice". Per-angle median signed residuals across
  fibers all ≤0.6 px vs the 9.1 px (2× floor) threshold; the one strong per-fiber pattern
  (fiber 08: a1–a3 negative, a4–a6 positive, ±9 px) is a half-turn (focus-session) split, not
  an angle-offset signature, and is absorbed by the split-half area_err.
- Sensitivity: forced-zero shifts change per-part median area by 2.9% median (14.2% worst,
  08_p4) and INFLATE the median axis ratio 1.141→1.197 — alignment removes spurious anisotropy,
  as designed. A_min positions: 5/15 within 10% of a span end (mild clustering; edge-trim
  consideration noted for the report).

**What I think it means (Task 12)** — The circular assumption is measurably wrong in shape
(all 15 fibers fit ratio >1.09; orientation persists along fibers) but nearly right in area
(A_mean within ±3% of π(d̄/2)² — mathematically expected, since mean-width circles nearly
area-match moderate ellipses). The headline anisotropy evidence is marginal at the per-fiber
Wilcoxon level (p=0.076) because per-image focus systematics (~8–17 px) swamp the ~4% direction
signal; the cross-part φ coherence test (Fisher p=1.4e-4) and hexagon consistency carry the
non-circularity claim instead. φ varies along fibers (twist), so a fixed-φ transfer fails — an
honest model-form limitation to report, not a fit bug (synthetic φ-transfer passes at constant
φ). Confidence: medium-high. Evidence that would change this: collaborator confirming a
different angle step (re-run Tasks 9–12), or C1 tensile data showing area errors matter more
than the ±3% found here.

**What I did (Tasks 10–11)** — `run_xsection.py` third-stage CLI (795ce71; per-part CSV/plot/
shifts JSON, per-fiber summary, angle residuals, `--validation` writer; 6 pipeline tests on a
fabricated known-ellipse tree). Synthetic end-to-end validation `tests/test_xsection_synthetic.py`
(7137a51): rendered bright-on-dark elliptical fibres (1200×500, 6 angles, injected shifts ±110 px,
noise 0.03) through the REAL measure pipeline → stack → fit.

**What I observed (Task 11, recovery table)** — per-side edge bias δ = 3.90 px (circle control,
z≈20 — matches the predicted fixed-z-on-smoothed-shoulder bias); circle ratio_med 1.0009
(pin ≤1.03); ellipse (a=210, b=165, φ=35°): a_med 217.8 vs δ-corrected truth 217.8, b_med 172.9
vs 172.8, φ_med 34.97°, area +0.10% (pin ±3%), valid 89.4%. Two pin deviations recorded:
(1) shift pin 1.5→2.5 px — the correlation's constant-extension edge fill biases the peak
~|shift|/span (2.0 px at 110/1200); width impact ~0.2 px, negligible, bounded later by the
shifts-on-vs-zero sensitivity check. (2) The ×2-noise area_err scaling pin is unattainable —
pipeline smoothing crushes per-column noise so the synthetic split-half error is a ~0.04%-of-area
systematic floor; replaced by a direct mechanism test (inject +2 px/side on the second half-turn
→ reported area_err within 6% of the analytic value, tolerance 25%).

**What I did (implementation, Task 1)** — Started executing the approved plan
(`docs/superpowers/plans/2026-08-15-multiangle-xsection-impl.md`, committed on branch
`worktree-multiangle-xsection` @ 9151b82). Data hygiene in `/Users/stan/Documents/UOM/spins/multiangle/C1/`:
renamed `C1_13_a6_part34.tiff_metadata.xml → C1_13_a6_part4.tiff_metadata.xml` (verified the target
was absent first); renamed the four `C1_09_a6_par5*` files → `C1_09_a6_part5*` (.tiff, .tiff_metadata.xml,
s.tiff, s.tiff_metadata.xml); moved the NUL-corrupted sync-conflict duplicate
`C1_07_a5_part3.tiff_metadata-C-LOSX3JVJ9.xml` into `C1/_quarantine/` (the real part3 sidecar is intact).
Added the strict C1 parsing API to `src/fibrecv/io_utils.py` (`MultiAngleKey`, `parse_multiangle_name`,
`multiangle_group`, `discover_multiangle`) + `tests/test_multiangle_names.py` (commit f19f222).

**What I observed (Task 1)** — Post-hygiene counts: 900 tiffs (450 plain + 450 s), 899 `*_metadata.xml`,
every plain image has a sidecar (`C1_14_a5_part1s.tiff` lacking one is known-benign — s-sidecars are
never read). Real-dir discovery smoke: 75 (fiber, part) keys × 6 angles = 450 images, zero missing,
fibers 1–15, parts 1–5. Full test suite 120 passed.

**What I did** — Surveyed the new multi-angle dataset `/Users/stan/Documents/UOM/spins/multiangle/C1/`.
Ran per-column edge extraction (threshold on grayscale, longest-run top/bottom edges) on fibers
01–06, all 6 angles × 5 parts, via throwaway scratchpad scripts (angle_structure.py /
angle_structure2.py, not kept). Three tests for the angle structure: aligned cross-correlation of
width-deviation profiles, mirror test (top-edge roughness of a_i vs negated bottom-edge roughness
of a_j), and ellipse-projection fit residuals under 30° vs 60° spacing. Extracted Zeiss ZEN
acquisition timestamps and scaling metadata. Ran a deep-planning design session; all seven design
decisions closed (see spec `docs/superpowers/specs/2026-08-15-multiangle-xsection-design.md`).

**What I observed** —
- Dataset: 15 fibers × 6 angles × 5 parts, each shot as plain + "s" (scale-bar) TIFF (2560×1920,
  8-bit RGB) + Zeiss XML sidecar. 450 measurement images total. Shot order angle-major (all 5
  parts at one angle, then rotate); ~1.5–2 h per fiber. Two filename anomalies:
  `C1_13_a6_part34.tiff_metadata.xml` (orphan XML; `part4`'s own metadata missing) and
  `C1_09_a6_par5[s].tiff` (typo).
- Mirror test (fibers 1–6, 30 parts): lag-3 pairs (a1,a4),(a2,a5),(a3,a6) show top-vs-negated-
  bottom roughness correlation +0.443/+0.472/+0.399 (reverse direction +0.355/+0.464/+0.410),
  all other pairings ~+0.29 (same-side baseline +0.31/+0.30/+0.29 at lags 1/2/3).
- Width-deviation profile correlation (fibers 1–3): lag-3 pairs +0.628/+0.619/+0.708; all other
  pairs ≤ +0.563.
- Ellipse-projection fit on per-part median widths (fibers 1–6, 30 parts): RMS residual median
  4.24 px under 60° spacing vs 5.32 px under 30° spacing.
- Some 180° pairs differ in median width by ~10% (fiber 01: a1 185.8 px vs a4 203.2 px).
- Zeiss metadata: StageRotation=NaN (rotation not recorded). CameraPixelDistance 2.2 µm,
  TotalMagnification 10× → 0.22 µm/px ≈ 4.545 px/µm; hardcoded `config.ppu=1.3680` px/µm does
  NOT apply to C1 (3.3× off). Fiber widths 185–205 px ≈ 41–45 µm at metadata scale.
- No tensile data for C1 exists yet (`MaSp2 10_06_26/` covers the old MasP2 fibers only).

**What I think it means** — The six angles cover a full revolution in 60° steps
(0/60/120/180/240/300°); a4/a5/a6 are the 180° flips of a1/a2/a3, so each cross-section has
3 independent directions × 2 repeats. Confidence: high for the 180° pairing (the mirror
signature is unique to it); medium-high for exactly-60° spacing (even steps is the natural
protocol and fits best, but exact within-half-turn spacing is not fully identifiable from
silhouettes — collaborator to confirm). The ~10% width mismatch within 180° pairs means manual
rotation error and/or focus-dependent edge bias; the fitter must expose per-angle residual
diagnostics. Evidence that would change this: collaborator stating a different step size, or
per-angle residuals showing the pairing is coincidental.

**Next** — Owner confirms 60° step with collaborator. Implementation starts in worktree
`worktree-multiangle-xsection` per the spec: scale from metadata, C1 name parsing, per-position
ellipse fit + circumscribed-hexagon QC bound, A_mean/A_harm/A_min export, synthetic +
leave-one-direction-out validation.
