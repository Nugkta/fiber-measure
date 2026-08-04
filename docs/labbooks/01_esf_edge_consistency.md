Hypothesis: Refitting each detected wall as a Gaussian-blurred step (erf) and moving the boundary to the fitted 50% midpoint removes the focus-dependent edge-placement bias, measurably reducing replicate-to-replicate diameter spread on MasP2 versus the current fixed-threshold edge.
Status: active
Started: 2026-08-03

## 2026-08-04
**What I did** — (later in the day) Merged `main` into the branch
(`git merge main`, commit `d6addb3`): main had gained the anomaly-flag
feature line (12 commits, `anomaly.py` + aggregate/GUI wiring). Resolved 5
conflicts (`compute.py` meta, `gui_app.py` docstring + checkbox branch,
three READMEs) by keeping both feature lines side by side. Then mined the
existing A/B outputs (`fibrecv_output/ab_refine_off` vs `ab_refine_on`) with
a scratchpad script to pick showcase cases for the owner's visual acceptance:
biggest per-image |Δdiameter|, biggest group-std wins/losses, largest fitted
σ, lowest refine coverage.
(Earlier) Reviewed the Task 5 M5 A/B outcome (null result on real
MasP2 data) and decided `refine_on` should not default to on until a proper
focus-sweep study validates it. Flipped `CONFIG.refine_on` `True` → `False`
in `src/fibrecv/config.py` and swept every doc/help-text spot that stated or
implied the default was on: `run_measure.py` `--refine` help text,
`GUI_README.md` checkbox description, `README.md` one-liner,
`docs/features/03_esf_edge_refinement.md` (behaviour paragraph + design-choice
bullet), `docs/superpowers/specs/2026-08-03-esf-edge-refinement-design.md`
(new revision footnote, historical body left untouched), and
`docs/report/01_esf_edge_consistency.md` §5 (appended the decision sentence,
no numbers/verdicts changed). Made 14 tests in `tests/test_refine.py`
explicit about `refine_on` (they exercised the on-path via the bare `CONFIG()`
default; now set `refine_on=True` explicitly since it is no longer the
default) — no assertions weakened. Ran the full suite. Later the same day,
the project owner clarified the MasP2 image-naming semantics: `X_Y_Z.jpg` is
condition (`X`) / sample-fibre number (`Y`) / position along that fibre
(`Z`), so a group `X_Y` is **one fibre** and its 2–3 images are different
positions along it, not repeat shots of the same spot. Corrected
`docs/report/01_esf_edge_consistency.md` (the dataset/method section, the
`registration_uncertain` interpretation, the §3.4 heterogeneity discussion,
and the Discussion/Limitations framing) accordingly — no numbers, table
values, or PASS/FAIL verdicts touched.
**What I observed** — Post-merge suite: 155 passed (124 ours + 31 from the
anomaly line). Acceptance-case shortlist from the A/B outputs: biggest movers
`7_3_3` (−8.0 µm, σ≈10), `10_8_3` (−7.6), `4_6_3` (−7.0); group-std wins
`4_6` (6.13→2.42 µm), `7_1` (7.53→3.63); regressions `7_4` (4.40→6.96 µm),
`10_4` (4.36→6.72); heavy-defocus/low-coverage `9_7_*` (σ≈16, coverage
0.29–0.80). Earlier: `uv run pytest -q` → 124 passed (same count as before
the change). No test needed a loosened assertion, only an explicit
`refine_on=True`/`replace(CONFIG(), refine_on=True, ...)` on tests that were
implicitly relying on the old default. Naming semantics (owner, 2026-08-04,
facts as stated): `X_Y_Z.jpg` → X = experimental condition, Y = sample
(fibre) number under that condition, Z = position along that fibre. A group
(`X_Y`) is therefore one fibre imaged at 2–3 different positions, not the
same spot re-photographed.
**What I think it means** — Confidence: high. The flip is purely a
configuration/documentation change — the algorithm itself is untouched and
still passes every synthetic ground-truth test — so refinement remains fully
available and correct, just opt-in (`--refine` / GUI checkbox) rather than
on by default, pending the focus-sweep study recommended in the report. The
naming correction withdraws the report's `registration_uncertain` reading
(33–36/46 groups) as an acquisition defect — that check assumed a group's
images share close to the same field of view, which by the true semantics
they never do, so it was never a meaningful diagnostic here. It does not
change the study's verdict: within-group spread now mixes real along-fibre
variation with measurement/acquisition noise, and per the owner's stated
expectation that along-fibre variation is small ("a bit of difference but not
much"), the ~23% median spread is still believed to be dominated by noise
upstream of the edge estimator, not by taper — so the primary metric's noise
floor still swamps the ~1–3 µm refine effect, and all acceptance verdicts and
numbers in the report are unchanged. Confidence: medium-high that noise
dominates over real taper (matches the owner's stated expectation, but
unquantified without new data); high that the registration-based reading is
withdrawn (a direct consequence of the corrected semantics). Evidence that
would change this reading: the focus-sweep / marked-segment validation
already recommended in the report, which would separate real along-fibre
variation from measurement noise directly.
**Next** — Design and run the focus-sweep validation (same marked segment,
deliberate focus settings) recommended in `docs/report/01_esf_edge_consistency.md`
§5/§6 before reconsidering the default.

## 2026-08-03
**What I did** — Assessed the "sharp air/silk step blurred by the optical PSF" idea on real data before any implementation. Wrote two throwaway analysis scripts (session scratchpad, since removed per repo policy): (1) ran the existing pipeline (`fibrecv.compute.compute_measurement`, default `CONFIG`) on `masp2 10_1_1.jpg`, `10_5_1.jpg`, `10_10_1.jpg` from `/Users/stan/Documents/UOM/spins/Images MasP2`, built edge-aligned 64-column mean profiles around each detected wall (aligned per column on the detected edge to avoid tilt smear; window −35 px outside → min(0.8·band_half, 28) px inside), and fitted the 4-parameter erf model `I(y) = a + (b−a)·½[1+erf((y−y0)/(√2σ))]`; (2) repeated the fits on three channels — saturation S, linearised R−G, linearised G — and compared fitted midpoints y0 pairwise on blocks where both fits had relative residual < 8%.

Implemented the design (Tasks 1–4 of `docs/superpowers/plans/2026-08-03-esf-edge-refinement-impl.md`): new `src/fibrecv/refine.py` stage inserted between `edges` and `qc` (block-wise erf fit of each wall over the perpendicular profile `t`, gated on step sign/residual/sigma-range/max-shift, interpolated between passing block centres, applied only to anchor columns); `refine_*` CONFIG fields (`config.py`); pipeline wiring + `meta["refine"]` diagnostics (`compute.py`); `--refine`/`--no-refine` CLI flag (`run_measure.py`); `refine_on` GUI checkbox + σ(x) profile-plot trace (`gui_app.py`); `tests/test_refine.py` groups A–D (identity/wiring, 1-D fit+gates, offset-field interpolation, full-pipeline integration). Eyeballed the refined output on `masp2 10_1_1.jpg` (read-only input; output image/plot to the session scratchpad only, then deleted) as the Task 2 gate. Commit range `ec70e0c..9b3fc7d`. Docs sync (this entry, `docs/features/03_esf_edge_refinement.md`, READMEs, header/docstring updates) is Task 4.

Ran the M5 A/B validation on the full set. Dataset `/Users/stan/Documents/UOM/spins/Images MasP2`: 144 JPEGs, of which 141 parse into 47 groups × 3 replicates (`3_1`–`3_5`, `4_2`–`4_8`, `7_1`–`7_10`, `8_1`–`8_10`, `9_5`–`9_9`, `10_1`–`10_10`); `teste`/`teste2`/`teste3` carry no group in their names and are excluded from group statistics. Three arms, each `uv run python -m fibrecv.run_measure --root "/Users/stan/Documents/UOM/spins/Images MasP2" --all --jobs 6 {--no-refine|--refine} --out fibrecv_output/<arm>` followed by `uv run python -m fibrecv.run_aggregate --out fibrecv_output/<arm> --all`: `ab_refine_off` (control), `ab_refine_on` (final tuned defaults), `ab_refine_on_rel008` (the as-landed gate, kept as the tuning before-picture). Three throwaway analysis scripts in the session scratchpad — a cached-pipeline knob sweep, the A/B analysis, the figure — all removed per repo policy. Tuned `CONFIG.refine_relmax` 0.08 → 0.15 in `src/fibrecv/config.py` as the sanctioned coverage fix, re-calibrated the two contamination-dose gate tests in `tests/test_refine.py` and added a midpoint-error-bound test, then re-ran BOTH arms with the final value so their `run_config.json` snapshots match. Wrote `docs/report/01_esf_edge_consistency.md` and `docs/report/01_esf_ab_summary.png`.

**What I observed** —
- erf fits on S: residual RMS / amplitude median 3.5–6.7% per image; blur σ median 7.3–9.3 px, p10–p90 ≈ 5.5–14.7 px, varying both between images and along a single fibre. Median diameters 60–120 px, so the ±2σ transition zone is a large fraction of the radius.
- Offset of the current threshold edge from the fitted 50% midpoint (per side, + = midpoint inside): **+2.3 px (10_1_1), −3.4 px (10_5_1), +3.4 px (10_10_1)** on S — i.e. an image-dependent diameter bias of up to ~7 px between images (G-channel fits gave per-side offsets up to +5.6 px).
- Channel comparison (blocks passing <8% residual): S passed 33/62, 50/67, 58/63 blocks; R−G passed 62/62, 55/67 but collapsed to 22/63 on the thin-fibre image (iridescent, non-spectrally-neutral specular rim suspected); G passed 51/62, 51/67, 16/63 (specular/shadow contamination). Fitted midpoints disagree between channels: median y0 differences 2–8 px with p10–p90 spreads up to ±20 px, and the S−G difference correlates with σ (r ≈ −0.5 to −0.7).
- Evidence figure: `docs/labbooks/01_esf_fits.png` (example aligned profiles + erf fits + normalised residual clouds for the three images, G channel).

Implementation-stage numbers (single-image, `masp2 10_1_1.jpg`, `refine_on=True` defaults; NOT the systematic A/B validation — that is Task 5):
- Refine coverage (refined / anchor columns): top wall 65.3%, bottom wall 23.1%.
- Median fitted blur σ: 11.65 px (top), 9.08 px (bottom).
- Median |applied offset| (|t0|): 5.66 px (top), 4.94 px (bottom).
- Median diameter: 110.16 px / 80.53 µm with `refine_on=False`, vs 115.23 px / 84.23 µm with `refine_on=True`.
- Test suite: 69 → 119 tests (50 new, `tests/test_refine.py` groups A–D).
- Perf: 0.32–0.59 s added per 2560-px image (< 1 s budget).

M5 A/B numbers (141 images, 46 of 47 groups usable in both arms — `7_7` has two `band_mismatch` replicates):
- **Coverage vs knobs.** As landed (`refine_relmax=0.08`) refinement covered 76.4% of anchor columns on average (median image 81.9%, 15.6% of images below 50%); the design spec's bar is ~80%. Block-width probe at that gate (12 images, both walls pooled): block 16 → 66.9%, block 32 → 70.3%, block 64 → 68.0%. Residual-gate sweep over all 141 images: 0.10 → 82.9%, 0.12 → 86.6%, **0.15 (adopted) → 89.6%** (median image 95.7%, 81.6% of images ≥ 80%, 2.8% below 50%), 0.20 → 91.0%, 0.25 → 91.6%. `refine_in_max` probe 28 → 20 → 16 → 12 (12 images, gate held at 0.08): coverage 66.9 → 79.0 → 85.5 → 89.3%, `masp2 10_1_1` median fitted σ 10.8 → 8.1 → 6.5 px, its bottom-wall median offset +4.9 → +2.0 → +0.5 px. Synthetic contaminated profiles (specular bump / shadow ramp of controlled amplitude): admitted midpoint error ≤ ~0.7 px at residual 0.08, ≤ ~1.5 px at 0.15; the pilot measured the legacy per-side bias at 2.3–3.4 px. Suite after the change: **124 passed**.
- **Effect of the gate change on the fits themselves** (`masp2 10_1_1`, `meta["refine"]` of the 0.08 arm vs the 0.15 arm). Top wall: coverage 65.3% → 71.6%, median σ 11.654 → 11.680 px, median |t0| 5.658 → 5.245 px. Bottom wall: coverage 23.1% → 99.1%, median σ 9.075 → 7.746 px, median |t0| 4.945 → 3.527 px.
- **Fits (tuned arm).** Per-image median fitted σ = 8.53 px (p10 5.62, p90 12.26), matching the pilot's 7.3–9.3 px medians and 5.5–14.7 px p10–p90. Median applied offset 2.46 px per side (p10 1.03, p90 4.78).
- **Primary metric** (within-group between-replicate std of the per-image mean diameter, µm): OFF mean 11.143 / median 9.004; ON mean 10.851 / median 9.435. **23 groups fall, 23 rise; exact two-sided sign test p = 1.000.** Median Δ +0.018 µm, median |Δ| 1.458 µm, worst fall −4.080 (`8_10`), worst rise +2.615 (`4_8`); two groups rise by more than their own off-arm noise (`10_5` 0.005 → 0.915, `10_7` 0.804 → 3.358). 6 of the 46 groups have only 2 usable replicates (`8_5`, `8_7`, `9_9`, `10_2`, `10_5`, `10_7`), including both of those two. Robustness on per-image medians: 25/46 improved (p = 0.659), mean 11.445 → 10.961. Stratified: off-std ≥ 5 µm → 18/33 improved (14.37 → 13.81 µm); off-std < 5 µm → 5/13 (2.94 → 3.34 µm). The gate value is irrelevant to this: improved-out-of-46 = 21 / 23 / 23 / 23 / 22 / 22 at relmax 0.08 / 0.10 / 0.12 / 0.15 / 0.20 / 0.25, median group std 9.004 off vs 9.31 / 9.51 / 9.49 / 9.44 / 9.39 / 9.39.
- **Replicate heterogeneity.** The median spread of the three replicate median diameters within a group is 22.9% of their mean, and 16 of 46 groups exceed 30% (`3_5` = 33.2 / 58.9 / 79.2 µm; `8_8` = 57.1 / 45.2 / 100.5 µm; `10_10` = 43.7 / 49.8 / 62.4 µm). `run_aggregate` reports `registration_uncertain=True` for 36 of 46 groups. Same primary metric on well-matched subsets: replicate medians within 10% (n = 10) → 3/10 improved (4.291 → 4.864 µm); within 20% (n = 20) → 8/20 (5.559 → 5.631); registration certain in both arms (n = 10, pooled pointwise std) → 5/10 (11.494 → 10.710).
- **Focus dependence.** Within-group centred, the measured diameter rises with the image's fitted blur σ: OFF slope **+1.87 µm per px of σ** (r = +0.162, n = 141), ON **+1.67 µm/px** (r = +0.146). Per-image mean-diameter change ON−OFF: mean −0.238 µm, sd 3.301, p10 −4.797, p90 +3.845, max |Δ| 8.043 µm; its correlation with σ is r = −0.044.
- **Secondaries.** Detrended along-fibre noise `std(raw−smooth)` got **worse**: 0.681 → 0.846 µm mean (median image 0.342 → 0.491), improving in only 8 of 141 images. `run_aggregate` pooled pointwise between-replicate std 13.733 → 13.377 µm mean, 24/46 improved (p = 0.883). QC unharmed and marginally better: `band_mismatch` 7 = 7, `low_confidence` 9 = 9, `coverage<50%` 3 = 3, mean QC coverage 95.38% → 95.50%, valid columns 340 746 → 341 152 (+406), `FLAG_ROLL_OUTLIER` 2 739 → 2 345 (−14%), `FLAG_CENTER_DEV` 12 478 → 12 466. Runtime over 144 images at `--jobs 6`: 1 min 56 s off vs 2 min 18 s on.
- **Not measured.** The planned secondary "σ(x) maps match visually defocused stretches" was not run; the diameter-vs-σ regression above was used in its place.
- Figure `docs/report/01_esf_ab_summary.png`; full write-up `docs/report/01_esf_edge_consistency.md`.

**What I think it means** — The blurred-step model is empirically valid (few-% residuals), and the current fixed-low-threshold edge placement carries a focus-dependent, image-dependent bias of several px per side that the σ-invariant erf midpoint would remove — this is the mechanism behind replicate disagreement. The edge is chromatic: different channels place it 2–8 px apart with focus-dependent disagreement, so a "true geometric edge" is not recoverable to better than a few px in these images; a consistent operational edge is the achievable goal, and S (= the pipeline's D map, an affine rescale) is the right channel to define it on — it is the only channel that stayed usable across all three images. Confidence: medium-high on the mechanism; the hypothesis's replicate-spread claim is untested until the A/B validation runs. Evidence that would change this reading: an A/B run showing no reduction in between-replicate spread.

The implementation-stage numbers are a sanity check, not evidence for the hypothesis: on this one image, refinement shifts the reported diameter by +4.7% (110.16 → 115.23 px) and fitted σ (9–12 px) sits well inside the range the earlier assessment measured (5.5–14.7 px), so the fit is behaving as designed. The bottom wall's low coverage (23.1%) is a specific, understood failure mode — its specular core stripe sits inside the inward fit window (capped at 28 px) and trips the residual gate — not a general problem with the method; but it is a live instance of the caveat now documented in `docs/features/03_esf_edge_refinement.md`: at low coverage, the median diameter reaching `qc` is still mostly the legacy estimator, not the erf one, so a wall like this bottom one gets only partial benefit from refinement as configured. Confidence: high on the numbers themselves (single deterministic run), low on what they imply for the hypothesis — one image, no replicate comparison. Evidence that would change this reading: the Task 5 A/B run.

M5, on the tuning first: the plan's suspicion that the 16-column blocks were noisier than the pilot's 64-column ones is **wrong** — widening the block barely moves coverage, so the residual is systematic model mismatch (specular stripe inside the window, shadow ramp outside), not noise that averaging can remove. That is why raising the residual gate is the knob that works. I read the `masp2 10_1_1` before/after as follows: on the top wall, where 0.08 already had a two-thirds sample, the extra blocks agree with the ones the tight gate accepted (σ moves 0.03 px, median offset 0.4 px), so the looser gate is safe there; on the bottom wall it is not a test at all — 23.1% coverage is not a representative sample, and going to 99.1% replaces the fit population rather than extending it. Nothing says the 23% subset was the more trustworthy of the two, only that a gate rejecting three quarters of a wall leaves no estimate worth defending. `refine_in_max` was rejected for a different reason: a shorter window refits blocks that were *already passing* through less data, changing the estimate itself, whereas `relmax` leaves every passing block's fit untouched and only adds more. Coverage bought by biasing good fits is not coverage worth having. The `FLAG_ROLL_OUTLIER` count falling 14% says the interpolation-seam risk the plan flagged did not materialise.

M5, on the result: **the hypothesis is not supported on MasP2, and this experiment could not have supported it.** The three replicates of a group plainly do not agree about the fibre itself — 22.9% median spread is not an edge-placement effect — so the primary metric's noise floor (~9–11 µm, set by replicates that are demonstrably not the same piece of fibre) is 3–10× the effect the erf midpoint can remove (~1.8 µm per side), so 23/46 at p = 1.000 is what a null result looks like when the measurement is dominated by the sample rather than the estimator. Acceptance against the spec: coverage **PASS** (89.6% mean, 95.7% median image), "primary falls in a clear majority" **FAIL** (23/46, p = 1.000), "rises nowhere beyond old noise" **FAIL** (2 groups). That second failure is the weakest of the three findings and I would not defend it hard: both offending groups are 2-replicate with a near-zero baseline (0.005 and 0.804 µm), so almost any change would breach it; restricted to 3-replicate groups, criterion 2 passes. Criterion 1 is the one that decides the verdict. What the run *did* establish is positive and worth keeping: the blurred-step model holds at scale (90% of anchor columns yield an accepted fit, σ matching the pilot's hand fits), and the focus-dependent bias is real and now quantified on 141 images at +1.87 µm per px of blur σ — a ~12 µm systematic swing over the observed σ range. The erf midpoint as implemented takes that only to +1.67 µm/px, i.e. removes ~11% of the focus dependence when the σ-invariant midpoint was supposed to remove essentially all of it, so a second mechanism dominates the residual: the ~10% of columns still on the legacy edge, and/or the wall-to-wall specular/chromatic asymmetry the pilot already saw (channel midpoints 2–8 px apart). The stage also costs along-fibre smoothness (+0.17 µm detrended noise) from the block-and-interpolate offset field. Confidence: **high** that refinement does not improve replicate consistency on this dataset; **low** that it would not help where replicates really are the same segment — untested. Evidence that would change this reading: a deliberate focus sweep of one fixed, marked segment, with σ as the independent variable and "diameter-vs-σ slope goes to zero" as the acceptance test — a design that needs no replicate assumption and isolates the estimator from the sample. One caveat I am carrying knowingly: I substituted that diameter-vs-σ regression for the planned "σ(x) matches visible defocus" eyeball check. The regression is the stronger test of whether σ still biases the measurement, but it cannot catch a σ(x) map that is smoothly mis-scaled — self-consistent and well-fit but systematically wrong — which is exactly what the eyeball check was for. So the +1.87 µm/px slope is trustworthy as a relative before/after, and only weakly grounded as an absolute focus scale.

**Next** — Design approved in chat (approach: erf-refinement stage after the existing wall finder, fallback to current edge on fit failure). Write the design spec, implement `refine` stage, then run the A/B validation over the full MasP2 folder: metric = within-group between-replicate std of mean diameter, old vs new.

Design spec, implementation (Tasks 1–3) and docs sync (Task 4) are done. Next: Task 5 — run `--refine` vs `--no-refine` over the full MasP2 folder, aggregate both, and compute the primary metric (within-group between-replicate std of per-image mean diameter, off vs on) plus secondary diagnostics (coverage %, median σ, detrended along-fibre noise). Status stays `active` until that validation lands.

Task 5 has now landed with a null result, so the study is **not** being closed on my own authority — status stays `active` pending three decisions: (1) acquire a proper replicate set (same marked segment, deliberate focus sweep) before any further edge-algorithm validation, since MasP2 cannot resolve a 1–3 µm estimator effect; (2) decide whether `refine_on` should stay `True` by default — it is roughly neutral on batch numbers (mean between-replicate std 11.14 → 10.85 µm, mean per-image diameter change −0.24 µm) but costs along-fibre smoothness and ~0.15 s/image; (3) if the mechanism is pursued, chase the residual +1.67 µm/px — joint two-wall fitting, and whether the two walls' offsets are correlated (adds to the diameter error) or anti-correlated (cancels).
