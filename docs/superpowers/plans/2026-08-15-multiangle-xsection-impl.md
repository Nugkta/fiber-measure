# Implementation plan — Study 02: multi-angle cross-section (`02_multiangle_xsection`)

## Context

`fibrecv` measures fiber diameters from single-view microscopy and feeds `tensile.py` a circular
area `A = π(d/2)²`; Young's modulus and toughness are ∝ 1/A, so the circular assumption's error
propagates 1:1 into the headline properties. The new C1 dataset
(`/Users/stan/Documents/UOM/spins/multiangle/C1/`: 15 fibers × 6 angles × 5 parts, TIFF + "s"
scale twin + Zeiss XML per shot) images each fiber from six rotation angles — 180° pairing proven
(a4/a5/a6 = flips of a1/a2/a3 → 3 independent directions × 2 repeats per position); exactly-60°
spacing is the working assumption pending collaborator confirmation (see Risks) — so the
cross-section can be measured instead of assumed.

Design is fixed by the approved spec `docs/superpowers/specs/2026-08-15-multiangle-xsection-design.md`
(decision record inside): ellipse fit primary + circumscribed-hexagon upper bound QC; per-fiber
A_mean/A_harm/A_min; v1 nominal angles + residual diagnostics (v2 offsets only if flagged); scale
from Zeiss XML with scale-bar cross-check; validation = synthetic unbiasedness + real-data
held-out prediction vs circle with 180°-pair noise floor; tensile join deferred (no C1 tensile
data yet). All work on branch `worktree-multiangle-xsection` in the existing worktree
`/Users/stan/Documents/UOM/spins/fiber-measure/.claude/worktrees/multiangle-xsection` (clean at
main 4be80c2); docs go to main at study close (study-01 convention). On execution start, copy
this plan to `docs/superpowers/plans/2026-08-15-multiangle-xsection-impl.md` and commit it.

## Plan-time discoveries (verified, load-bearing)

1. **LZW decode failure** — `io_utils.load_rgb` → imageio → tifffile raises
   `ValueError: COMPRESSION.LZW requires the 'imagecodecs' package` on every C1 TIFF;
   `iio.imread(path, plugin="pillow")` decodes fine. Fixed in Task 2, no new dependency.
2. **Desaturation feature has zero signal on C1** — `compute_measurement` returns coverage 0.00
   on C1 probes. C1 is a bright fiber (V≈0.80–0.88) on dark gray bg (V≈0.18–0.24) with no
   saturation contrast (S≈0.06–0.12 both sides). Needs a config-selected brightness z-map mode
   (Task 5); everything downstream of the z-map (band/edges/qc/tilt) is feature-agnostic.
3. **Two conflicting scale fields in the Zeiss XML** (re-verified directly):
   `CameraPixelDistances 2.2,2.2` (NB: plural tag) with 10× magnification → 0.220 µm/px, but
   `Scaling/Items/Distance[X] = 3.8892e-7 m` → 0.3889 µm/px — 1.768× apart (3.1× in area).
   The "s" scale-bar cross-check (Task 4) is a HARD GATE that adjudicates before any µm output.
4. **Review-verified data quirks** (adversarial review 2026-08-15): 900 tiffs + 899 sidecars —
   `C1_14_a5_part1s.tiff` has no sidecar (benign; s-sidecars never read); one NUL-corrupted
   sync-conflict duplicate `C1_07_a5_part3.tiff_metadata-C-LOSX3JVJ9.xml` (real part3 sidecar
   intact). Fiber brightness spans ~2×: fiber 01 V≈0.80 (z≈20) but e.g. C1_08_a3_part2 V≈0.51
   (z≈8) — only 2× above `k_band=4`/`edge_z=4`.

## Global constraints

- No new dependencies; `uv run` for everything; run all commands from the worktree root.
- `io_utils.parse_name` and its 2-tuple untouched (6 call sites destructure positionally); C1
  parsing is a NEW `parse_multiangle_name` + adapter. MasP2 behavior stays bit-identical
  (`feature_mode="desat"` default); per-image CSV changes are additive columns only; the
  existing test suite stays green after every task.
- All fitting math in px; µm conversion once, at the xsection stage, from the per-image XML
  scale (never `cfg.ppu`); error loudly if the six sidecars of a (fiber, part) disagree (>0.1%).
- `xsection.py` is pure (no I/O), mirroring the compute/measure split; all I/O in `run_xsection.py`.
- TDD per task: failing tests → implement → gate green → commit (`feat(scope): …`).
- Throwaway analysis scripts in the session scratchpad only; data renames recorded in the labbook.

## Tasks

### Task 1 — C1 data hygiene + name parsing + discovery
Files: `src/fibrecv/io_utils.py`; new `tests/test_multiangle_names.py`.
- Rename in the data dir (outside repo; verify `C1_13_a6_part4.tiff_metadata.xml` absent first):
  `C1_13_a6_part34.tiff_metadata.xml → C1_13_a6_part4.tiff_metadata.xml`; the four
  `C1_09_a6_par5*` files → `part5*`. Also move the NUL-corrupted sync-conflict duplicate
  `C1_07_a5_part3.tiff_metadata-C-LOSX3JVJ9.xml` into a `_quarantine/` subfolder (real part3
  sidecar is intact). Record all moves in the labbook.
  Post-checks: 900 tiffs (450 plain + 450 s), 899 `*_metadata.xml` in place, and every one of
  the 450 PLAIN images has an intact sidecar (`C1_14_a5_part1s.tiff` lacking one is known-benign).
- New API (regex `^(?P<cond>[A-Za-z]+\d*)_(?P<fiber>\d{1,3})_a(?P<angle>\d)_part(?P<part>\d)(?P<s>s?)$`):
  - `MultiAngleKey(condition, fiber, angle, part, scalebar)` frozen dataclass
  - `parse_multiangle_name(path) -> MultiAngleKey` (ValueError on no match; strict — no `par5` special case)
  - `multiangle_group(key) -> (f"{cond}_{fiber:02d}_a{angle}", part)` adapter to (group, replicate) shape
  - `discover_multiangle(root, condition=None) -> dict[(fiber, part), dict[angle, Path]]`
    (plain images only; s twins/XML excluded; missing angle = absent key, not error)
- Tests: accept/reject parametrize table; adapter; discovery on a fabricated tmp_path tree;
  existing `test_io_utils.py` untouched and green.
- Gate: real-dir smoke prints 75 (fiber, part) keys × 6 angles = 450, zero missing.

### Task 2 — LZW TIFF loading
Files: `io_utils.load_rgb`; extend `tests/test_io_utils.py`.
- Test: Pillow-written LZW TIFF in tmp_path → `load_rgb` returns float32 (H,W,3) in [0,1].
- Implement: on `iio.imread(path)` exception, retry `iio.imread(path, plugin="pillow")`
  (generic path first → existing JPEG behavior byte-identical).
- Gate: real C1 tiff loads → (1920, 2560, 3).

### Task 3 — Zeiss XML scale reader
Files: new `src/fibrecv/scale.py`; new `tests/test_scale.py`. Stdlib ElementTree; `.//` searches.
- API: `ScaleInfo(um_per_px_camera, um_per_px_items, camera_pixel_distance_um, total_magnification, source_path)`;
  `sidecar_for(image_path)` (append `_metadata.xml` to the FULL filename);
  `read_scale(xml_path)`; `resolve_um_per_px(info, source)` with source ∈ {"camera","items",numeric literal}.
- Field names (review-verified against real sidecars): the tag is **`CameraPixelDistances`**
  (plural; value `"2.2,2.2"`) — a singular search finds nothing. Magnification = product of
  `NominalMagnification` × `CameraAdapterMagnification` × `OptovarMagnification` (= 10×1×1 here),
  NOT the eyepiece `TotalMagnification` (numerically identical here but wrong provenance).
- Tests: minimal XML fixtures (UTF-8 BOM, both fields; X≠Y → ValueError; missing node → ValueError).
  Real-sidecar spot-check (expect 0.220 / 0.38892) runs as a scratchpad command or a test guarded
  by `pytest.mark.skipif(not Path(...).exists())` — the committed suite must pass off this machine.

### Task 4 — Scale-bar cross-check (HARD GATE, scratchpad, no committed code)
- Scratchpad script: crop the burned-in scale bar from 3 "s" twins spread across fibers
  (e.g. C1_01_a1_part1s, C1_08_a3_part2s, C1_15_a6_part5s), read the µm label visually from the
  saved crops, measure bar length in px (bright/dark run threshold), compute µm/px each.
- Decision: whichever XML field matches within ~2% becomes `--scale-source` for all later runs
  (if "items" wins, add a dated correction entry to labbook 02 — the spec's recorded overturn
  condition). Neither matches → STOP, escalate to owner. Record numbers + crop paths in the
  labbook BEFORE Task 8's full batch. Expected: 0.22 wins (2.2 µm is the AxioCam ERc5s native
  pixel; TIFFs are native-resolution 2560×1920) — but do not assume.

### Task 5 — Brightness feature mode
Files: `config.py`, `features.py`, `run_measure.py`; new `tests/test_features_modes.py`.
- `config.feature_mode: str = "desat"` (`"desat"` MasP2 | `"bright"` C1); `--feature-mode` CLI flag.
- `features.py`: `"bright"` branch uses V = `rgb.max(axis=2)` with the SAME margin-row background
  median + MAD machinery: `D = (V − v_bg)/(mad_scale·MAD_v)`; return contract unchanged.
- Tests: bit-identical regression guard for `"desat"` on the existing pink fixture (assert_array_equal);
  `_bright_fibre()` fixture (adapted from `tests/test_edges_tilt.py:19-36 _inclined_fibre`,
  bright ≈0.85 on dark ≈0.20) reaches coverage ≥0.9 end-to-end under `"bright"` and ~0 under
  `"desat"`; a DIM variant (fiber ≈0.50 on ≈0.20 → z≈8, the C1 dark end per review) must also
  reach coverage ≥0.9; invalid mode → ValueError. Note in the docstring that `meta["bg_S"]`
  holds v_bg under bright mode (semantic overload, mode recoverable from `params`).
- Gate: the three C1 probe images reach coverage ≥ ~0.9, width ~185–205 px. Knob tuning (if any)
  via CLI-exposed flags only; record values in labbook.

### Task 6 — Persist per-column edges in the per-image CSV
Files: `measure.py` (`write_measurement`); small new test.
- Add `y_top_px`, `y_bot_px` (span-restricted, NaN-invalid) as ADDITIVE columns; all existing
  columns unchanged; `run_aggregate.main` still exits 0 on a MasP2-style-named fixture tree
  (on a C1-only tree it returns 1 "No groups matched" by design — don't use C1 names here).
- Rationale: widths already persist per angle (all the fit needs); edges enable Task 12's
  focus/asymmetry diagnostics and any v2 work.

### Task 7 — Shared resample helper
Files: `register.py`; unit test.
- Extract `register.py:123-138` verbatim into `resample_to_grid(aligned) -> (grid, stack)`;
  `register_sample` calls it; existing aggregate tests are the no-change guard. Add an all-NaN
  input guard + test (current inline code crashes on `min()` of an empty generator — a fully
  invalid (fiber, part) must degrade gracefully, not crash).

### Task 8 — Measurement batch (pilot → full 450)
No new code. From the worktree root, `<PPU>`/`<SCALE-SOURCE>` per Task 4 (ppu cosmetic only —
xsection re-derives µm from XML):
```bash
uv run python -m fibrecv.run_measure \
  --root /Users/stan/Documents/UOM/spins/multiangle/C1 \
  --glob "C1_01_a*_part[0-9].tiff" \
  --feature-mode bright --ppu <PPU> --jobs 6 \
  --out fibrecv_output/multiangle_c1
```
- Pilot gate (fiber 01, 30 images — the BRIGHT end; also pilot one dim fiber, e.g. fiber 08):
  0 errors; coverage ≥0.8 on ≥90%; EYEBALL ≥6 overlays (one per angle) — boundaries on true
  walls. During the pilot also record the real cross-angle `corr_peak` distribution (feeds the
  final `xsec_min_corr` value). Tune → delete out dir → rerun; freeze flags in labbook.
- Full batch: glob `"C1_*_a*_part[0-9].tiff"` (exactly the 450 plain images). ~10–20 min.
  (If Task 11's synthetic proof later changes a knob, rerunning the batch costs ~20 min — accepted.)
- Audit (scratchpad): run_log errors = 0; coverage distribution STRATIFIED by image brightness
  (the ~2× V spread — dim images are the failure risk); LOWCONF list short and understood.

### Task 9 — xsection math module (the numerical heart; strict TDD)
Files: new `src/fibrecv/xsection.py` (pure: numpy + `register.estimate_shift`/`resample_to_grid`);
new `tests/test_xsection.py` (~25 tests on hand-built arrays, no images).
- **Parameterization** (kills both degeneracies): fit in squared-width space, exactly linear:
  `w²(θ) = c0 + c1·cos2θ + c2·sin2θ`; `R = √(c1²+c2²)`; `a² = c0+R`, `b² = c0−R` (a ≥ b by
  construction), `φ = wrap(−atan2(c2,c1)/2, [0°,180°))`, NaN when near-circular (R < 1e-6·c0).
  3-param linear LSQ per column over ≤6 points; vectorize by grouping columns on the 2⁶
  finite-mask patterns. Column valid iff all 3 direction pairs have ≥1 finite width AND b² > 0
  (negative → invalid + counted, never clamped).
- API: `NOMINAL_ANGLES_DEG`; `PartStack(fiber, part, x, W(6,N), shifts)`;
  `build_part_stack(fiber, part, profiles, cfg)` (ref = lowest angle). **B2 fix (review):**
  per-image CSVs are span-restricted with per-image x0 — resample every profile onto the
  absolute full-width x grid (NaN-padded from its own `x_px` column) BEFORE `estimate_shift`,
  else differing x0 injects a silent (x0_ref − x0_other) misalignment. **S1 fix:** wire the new
  `cfg.xsec_min_corr = 0.2` explicitly via `replace(cfg, min_corr=cfg.xsec_min_corr)` around the
  `estimate_shift` call (the threshold is read inside as `cfg.min_corr`), with a test asserting
  the 0.2 gate actually takes effect. (0.2 is provisional — re-measure the real `corr_peak`
  distribution during the Task 8 pilot before freezing.)
  `XsecFit(a, b, phi_deg, area, resid(6,N), rms_resid, n_angles, valid)`;
  `fit_ellipse_projections(W, theta_deg)`; `hexagon_area(W)`; `hexagon_area_expected(a, b, phi_deg)`;
  `split_half_area(W, theta_deg)`; `pair_differences(W)`;
  `predict_anisotropy(W)` → (directional-predictor errors, circle-baseline errors) and
  `predict_phi_transfer(W, theta_deg)` (Task 12 criteria).
- **Hexagon**: direction widths = nanmean of each 180° pair, h_d = w_d/2;
  `A_hex = (2/√3)·(2(h₀h₁+h₁h₂+h₂h₀) − (h₀²+h₁²+h₂²))`; valid iff all 3 directions present and
  h_d < h_e + h_f for all rotations (else NaN + `hex_degenerate` — NaN-ing is REQUIRED: beyond
  degeneracy the closed form under-states the true area and would silently break the upper-bound
  property; state this in the feature doc). Circle anchor: ratio π/(2√3) = 0.9069 — **exact for
  circles only** (review-verified: a true ellipse at axis ratio 1.5 gives 0.883–0.893, at 2.0
  gives 0.827–0.876, φ-dependent). QC must therefore compare measured hex_ratio against the
  value EXPECTED for the fitted (a, b, φ) — add `hexagon_area_expected(a, b, phi_deg)` (support
  widths of the fitted ellipse → same hexagon formula) so the flag is
  `hex_ratio / hex_ratio_expected ≠ 1`, not `hex_ratio ≠ 0.9069`.
- **Per-position uncertainty**: `area_err = |A₁₂₃ − A₄₅₆|/2` from `split_half_area` (each half is
  a complete 3-direction exact fit); doubles as the pair noise floor; documented as conservative
  (absorbs 180° focus asymmetry).
- Tests: exact recovery to ~1e-9 over φ ∈ {0,30,89,90,91,179}° × ratio ∈ {1.0,1.05,1.3,2.0};
  a ≥ b always; circle → φ NaN; missing-data table (one angle NaN → valid; full pair NaN →
  invalid; 3 angles one-per-pair → valid, resid ≡ 0); Monte-Carlo under noise — SCOPED (review
  S5): assert area ≈ unbiased and a/b/φ unbiased only for well-separated axes; at circularity
  R̂ is Rician-biased upward so assert the axis ratio bias is POSITIVE and small, don't assert
  zero; hexagon vs brute-force half-plane clip (shoelace) on random triples, and
  `hexagon_area_expected` matches the clip on exact ellipse support widths; degeneracy guard at
  h₀ = h₁+h₂; build_part_stack recovers injected shifts ±0.5 px INCLUDING profiles with unequal
  per-angle x0 (the B2 regression test), flat profile → uncertain flag.

### Task 10 — `run_xsection.py` CLI + output schema
Files: new `src/fibrecv/run_xsection.py`; `__init__.py` stub print updated; new
`tests/test_xsection_pipeline.py`. No pyproject change (matches run_measure/run_aggregate convention).
```bash
uv run python -m fibrecv.run_xsection \
  --out fibrecv_output/multiangle_c1 \
  --data-root /Users/stan/Documents/UOM/spins/multiangle/C1 \
  [--condition C1] [--fibers ...] [--scale-source camera|items|<um_per_px>] \
  [--min-corr 0.2] [--validation]
```
Reads `per_image/csv/*_profile.csv` (+ meta JSONs; excluded images → NaN angle rows), keys via
`parse_multiangle_name`, scale via `scale.sidecar_for` (error if the 6 sidecars disagree >0.1%).
Outputs:
```
per_part/csv/xsec_C1_<ff>_p<p>.csv    x_px, w_a1..w_a6_px, n_angles, a_px, b_px, phi_deg,
                                      rms_resid_px, a_um, b_um, area_um2, area_err_um2,
                                      area_hex_um2, hex_ratio, hex_ratio_expected, valid
per_part/plots/xsec_C1_<ff>_p<p>.png  6 aligned width profiles / A(x) ± err + hex bound
per_part/shifts/xsec_C1_<ff>_p<p>.json  shift/corr/uncertain per angle, both XML scales + resolved
summary/xsection_summary.csv          per fiber: n_parts, n_positions, um_per_px, A_mean_um2,
                                      A_harm_um2, A_min_um2(+part,+x), axis_ratio_median/iqr,
                                      phi_med_deg, hex_ratio_median, pair_dw_frac_median,
                                      A_circle_um2, area_ratio, uniformity, n_uncertain_shifts,
                                      low_confidence
summary/xsection_angle_residuals.csv  per (fiber, angle) residual diagnostics (v2 gate evidence)
summary/xsection_run_config.json      provenance
```
A_min taken on an 11-px rolling-median-filtered A(x) (single bad column can't set the weakest
link); `A_circle_um2 = π(d̄/2)²` with d̄ = grand mean width (tensile.py convention, descriptive
comparison); A_harm = n/Σ(1/A). Tests: fabricated mini-tree from a known ellipse → exit 0, exact
schemas, hand-computed summaries match; missing angle / excluded image / disagreeing scales /
`--fibers` filter covered.

### Task 11 — Synthetic end-to-end validation (implementation proof; committed tests)
Files: new `tests/test_xsection_synthetic.py`.
`_elliptical_fibre_images(...)`: bright-on-dark render (Task 5 geometry) with per-column
half-width from a(x), b(x), φ; mild taper; per-angle x-shifts (±≤120 px) and independent noise
seeds; all 6 nominal angles (θ/θ+180 same geometry, different seeds). W=1200, H=500 for speed.
Assertions (initial pins; record-and-tighten per study-01 pattern — ratio and φ stay centered):
circular control a=b=190 → ratio ≤1.03, measure edge bias δ; elliptical a=210, b=165, φ=35° →
φ ±4°, ratio ±0.05 of the δ-corrected truth (a+2δ)/(b+2δ), area ±3% of π(a+2δ)(b+2δ)/4; shifts
±1.5 px; hex_ratio within ±0.02 of `hexagon_area_expected` for the rendered ellipse (NOT the
circle anchor 0.9069 — at ratio 1.27 the true-ellipse value is ~0.90, review S4); area_err
scales ~linearly with injected noise (×2 → ×1.5–3).
Gate: labbook entry with δ and recovered-vs-true table (feeds the report).

### Task 12 — Real-data run + held-out validation (model proof)
Files: `predict_anisotropy` + `predict_phi_transfer` in `xsection.py` + unit tests (synthetic W:
strongly elliptical → directional error ≈ noise while circle error ≫; circular → the two are
statistically indistinguishable); `--validation` writer in `run_xsection.py`; analysis in
scratchpad; figures into `docs/report/`.
- **Criterion resolution (document in report; sharpened by adversarial review B1):** the spec's
  literal "fit on 2 of 3 directions" is rank-deficient (θ and θ+180° share a design row → rank 2
  for 3 params; the held-out direction is unconstrained). Moreover, in leave-one-repeat-out the
  LSQ fit passes EXACTLY through any single-observation direction, so the "ellipse prediction"
  of a dropped width is identically its 180° partner — the ellipse machinery adds nothing to
  that prediction, and its error is by construction the pair difference that defines σ_pair.
  Honest framing, implemented as:
  - **Primary — anisotropy test (leak-free):** for each (fiber, part, angle k) with a finite
    partner w_{k+3} (positions with NaN partner are EXCLUDED — rank-2 otherwise), the
    directional predictor is w_{k+3}; the isotropic baseline is the mean of the other 5 widths.
    Per-fiber RMSE, paired Wilcoxon across 15 fibers, p<0.05 + median relative error reduction.
    This tests "width is direction-dependent (non-circular section)" — which is the load-bearing
    claim for replacing π(d/2)². The noise floor σ_pair = median(|w_k − w_{k+3}|)/√2 is plotted
    against the CIRCLE RMSE only (plotting it against the directional RMSE would compare a
    quantity with itself).
  - **Model-form evidence (ellipse specifically):** with 3 directions the ellipse form is
    saturated and cannot be validated by held-out widths — the report must say so. Ellipse
    adequacy is instead judged by (a) hex_ratio vs `hexagon_area_expected` consistency and
    (b) the **secondary φ-transfer test**: φ fitted on odd grid columns (all 6 angles), then at
    even columns drop BOTH repeats of one direction, solve (c0, R) from the remaining 2
    directions with φ fixed, predict the dropped direction (now well-posed, rank 2 for 2
    unknowns); reported with the φ-transfer declared.
  - Spec fallbacks: negative result written up; median axis ratio <1.05 → stratify by ellipticity.
- Full run: `run_xsection … --scale-source <per Task 4> --validation` → `summary/xsection_validation.csv`.
- Scratchpad figures → `docs/report/`: `02_xsec_anisotropy_paired.png` (headline: per-fiber
  paired RMSE, directional vs circle, noise floor drawn against the circle RMSE only),
  `02_xsec_angle_residuals.png`, `02_xsec_hex_ratio.png` (measured vs expected-for-fit),
  `02_xsec_example_part.png`, `02_xsec_area_vs_circle.png`, `02_xsec_axis_ratio.png`.
- **v1→v2 gate:** open a v2 angle-offset follow-up only if some angle's median signed residual
  > ~2× pair noise floor consistently across fibers AND the pattern doesn't track edge asymmetry
  (`y_bot−y_top` from Task 6) instead. Otherwise record "nominal angles suffice".
- Sensitivity checks: fits with shifts on vs forced-zero on 2 fibers; watch shifts pinned at ±400;
  check whether A_min positions cluster at part-span ends (defocus/vignette → edge-trim paragraph).

### Task 13 — Docs + close-out
- `docs/features/03_multiangle_xsection.md` (dirs hold 01, 02, 04 → next UNUSED per CLAUDE.md
  is 03; re-verify at execution; CLAUDE.md template; Design choices: w²-space linear fit,
  split-half uncertainty, A_min median guard, anisotropy-test-vs-literal-LODO rank argument,
  xsec_min_corr, hexagon NaN-on-degeneracy as an upper-bound requirement; Caveats:
  convex-hull-only silhouettes, nominal angles, focus-bias entanglement, scale-field
  adjudication, ellipse form unverifiable from 3 directions).
- `src/fibrecv/README.md` + root `README.md`: three-stage pipeline + new module bullets; header docstrings.
- `docs/report/02_multiangle_xsection.md`: intro/method/results/discussion, synthetic table,
  Wilcoxon stats, six figures.
- Labbook prepend-entries throughout; metalabbook row → done + re-sort at close.
- Full suite green; then superpowers:finishing-a-development-branch — merge decision is the
  owner's (code on branch; docs on main).

## Risks (watch actively)

- **The 60° step is an assumption, not a confirmed fact** (labbook: 180° pairing proven, exact
  spacing medium-high confidence, collaborator confirmation pending — the spec's recorded
  overturn condition). If the within-half-turn spacing differs, the fit design matrix AND the
  hexagon's 60°-apart normals both change: keep `theta_deg` a parameter everywhere (it already
  is), and re-run Tasks 9–12 with the confirmed angles if the answer differs.
- **Scale ambiguity** (1.768× → 3.1× in area): Task 4 hard gate; both XML values carried in every
  shifts JSON for audit.
- **C1 optics vs MasP2-calibrated knobs, and the ~2× brightness spread**: fiber 01 has huge
  headroom (z≈20) but the dim end (e.g. C1_08_a3_part2, z≈8) sits only 2× above the z=4 gates —
  Task 8 pilots BOTH a bright and a dim fiber, the audit stratifies coverage by brightness, and
  Task 5/11 fixtures include the dim variant.
- **Cross-angle alignment weaker than replicate alignment** (best corr 0.63–0.71): misalignment
  biases axis ratio toward 1; xsec_min_corr gate + uncertain flags + shifts-on-vs-zero sensitivity.
- **NaN patterns**: columns valid only with all 3 direction pairs; partial columns invalid and
  counted, never interpolated across angles.
- **A_min is an extreme statistic**: rolling-median guard + end-clustering check.
- **Focus bias vs ellipticity vs angle error**: separated by per-angle residuals + edge-asymmetry
  cross-check before any v2 decision; split-half area_err absorbs it conservatively.
- **180° pairs are NOT x-mirrored** (rotation about the fiber's own axis preserves x) — no
  profile reversal in alignment.
- **Don't run `run_aggregate` on C1** — parse_name skips every C1 name by design; run_xsection
  is the aggregation stage here.

## Verification

- After every task: `uv run pytest -q` green in the worktree (existing suite + new tests).
- Task 4 gate: scale-bar µm/px matches exactly one XML field within ~2%.
- Task 8 gate: 450/450 measured, ≥ ~95% coverage ≥0.8, overlays eyeballed per angle.
- Task 11 gate: synthetic recovery table within pins (φ ±4°, area ±3% δ-corrected).
- End-to-end: `run_xsection --validation` produces per-part CSVs (75), xsection_summary.csv
  (15 rows), validation CSV; headline figure shows ellipse vs circle RMSE per fiber against the
  noise floor; labbook + report updated.

## Critical files

- `src/fibrecv/io_utils.py` — C1 parsing, discovery, LZW loading (Tasks 1–2)
- `src/fibrecv/scale.py` (new) — Zeiss XML scale reader (Task 3)
- `src/fibrecv/features.py` — bright-fibre feature mode (Task 5)
- `src/fibrecv/measure.py` — y_top/y_bot persistence (Task 6)
- `src/fibrecv/register.py` — estimate_shift reuse + resample_to_grid extraction (Task 7)
- `src/fibrecv/xsection.py` (new) — alignment stack, ellipse + hexagon + LORO (Tasks 9, 12)
- `src/fibrecv/run_xsection.py` (new) — third-stage CLI (Task 10)
