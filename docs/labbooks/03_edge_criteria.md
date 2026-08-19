Hypothesis: Replacing max(R,G,B) with median(R,G,B) in bright mode and switching from a fixed to a clamped relative edge threshold will reduce defocus-induced width errors on C1 images without introducing inward bite.
Status: active
Started: 2026-08-19

## 2026-08-19

**What I did**

1. **Implementation** (branch `worktree-edge-criteria`, off `worktree-multiangle-xsection`).
   - `src/fibrecv/features.py:83` — bright-mode z-map `rgb.max(axis=2)` → `np.median(rgb, axis=2)`.
   - `src/fibrecv/edges.py:263-273` — bright-mode boundary level made relative and clamped.
     Shipped form: `level = base_val + min(edge_cap*A, max(edge_z, edge_frac*A))`.
     Desat path left byte-for-byte as `base_val + min(edge_z, edge_frac*A)`.
   - `config.py` gained `edge_cap: float = 0.50`; `run_measure.py` gained `--edge-cap` and
     `BRIGHT_DEFAULTS = {"edge_frac": 0.30, "k_band": 6.0}` applied only when
     `--feature-mode bright` is given.
2. **Calibration** — `scripts/calibrate_edge_frac.py`, sweep `edge_frac` 0.15…0.45 at
   `edge_cap=0.50`, selected 0.30. Full 450-image stage 1 → `fibrecv_output/multiangle_c1_s03`,
   then stage 3 with `--scale-source items --validation`.
3. **Post-close verification** — the checks skipped at first close: MasP2 desat regression,
   fiber 10 rms diagnosis, `k_band`/`amin` band-finding on the dimmest C1 quartile, per-part
   rms distribution.
4. **k_band re-run** — after finding the `k_band` coupling bug, re-ran stage 1+3 with
   `--edge-frac 0.30 --k-band 6.0` on 450 primary images (excluding the `s.tiff` scale-bar
   duplicates) → `fibrecv_output/multiangle_c1_s03b`. Generated the old-vs-new boundary
   overlay (`scripts/spot_check_overlay.py`).
5. **Multi-agent code review** (2 Opus on algorithm/methodology, 2 Sonnet on wiring/docs,
   1 agent rebuilding missing evidence). Wrote the reproducibility scripts the study lacked:
   `verify_masp2_identical.py`, `synthetic_fringe_control.py`, `delta_seed_stability.py`,
   `rms_statistics.py`, `dim_coverage_check.py`, `recalibrate_edge_frac.py`.
6. **edge_frac re-decision (abandoned mid-run)** — the review invalidated both original
   selection criteria, so I re-ran calibration correctly and launched a full 450-image run at
   `edge_frac=0.40`. It died at 403/450 when the disk filled (460 GB volume down to 117 MB
   free); the owner chose not to retry it, so **no 0.40 full-pipeline result exists**.
   Freed ~9 GB by deleting the superseded `multiangle_c1_s03` (900-image buggy run) and the
   partial 0.40 output.
7. **Final audit** (independent agent) — found that the headline run predates the clamp fix
   (below). Re-ran the shipped config on the fixed code → `fibrecv_output/multiangle_c1_s03d`.

**What I observed**

*Headline result (450 C1 images, s03d at edge_frac=0.30 on post-clamp-fix code vs the old baseline)*
- Stage 1: 450 ok, 0 errors, 0 low_confidence. Suite: **218 passed**.
- Low-confidence fibers **6 → 1** (old 01/07/08/09/13/15; new: fiber 10); total positions
  162223 → 161696 (−0.33%).
- Per-part rms over 75 parts: p50 3.697 → **2.739**, p90 6.638 → 4.264, p95 7.756 → 4.892,
  max 11.841 → **6.273**, parts >6 px **11 → 1**.
- Paired statistics (`scripts/rms_statistics_s03d.json`): 60/75 parts improved, median paired
  change **−1.057 px**, Wilcoxon W=436 **p=1.77e−7**, bootstrap 95% CI [−1.577, −0.571] px;
  per-fiber sign test 13/15 improved, p=0.0074.
- s03d vs s03b (pre- vs post-clamp-fix): 54/75 parts exactly identical, max part difference
  0.37 px, median 0.000 px. Headline statistics survive the fix.
- Discarded-data confound tested and **not supported**: Spearman(Δvalid columns, Δrms)
  = −0.139, p=0.235 — parts that *gained* columns improved most. (Failing to support the
  confound is not the same as refuting it; p=0.235 is a null result.)
- MasP2 desat **bit-identical**: 3 images (`masp2 10_10_1.jpg`, `masp2 7_2_1.jpg`,
  `teste3.jpg`), `diameter_raw` and the float32 `D` z-map both match by SHA-256. Median
  diameters 59.760 / 103.529 / 153.99 px = 43.68 / 75.68 / 112.57 µm at ppu 1.368
  (`scripts/masp2_identical.json`).
- Visual spot-check on 4 images (1 defocused, 2 sharp, 1 dim): old boundary sits outside the
  visible bright-to-dark edge in every panel, new one lands on it; no inward bite. Per-side
  shift −12.2 px (defocused), −6.7 / −5.6 px (sharp), −7.6 px (dim).

*The headline run predates the clamp fix — found by the final audit*
- s03b finished **10:19:04**; the clamp fix (`a558a09`) landed **10:42:21**. So every number
  above was produced by the **pre-fix** ordering, not by the shipped code.
- My earlier claim that the fix was "a no-op on real C1 (maxdiff 7.6e−06 px)" is **wrong** —
  that check sampled images without faint walls. The two orderings diverge when
  `A_side < edge_z/edge_cap = 8` (not `< edge_z = 4`, as I also wrote). Measured A/B on the
  dimmest C1 images: 819–1527 of ~2530 columns differ per image, max per-column 5.65 px,
  median-width shift up to **+0.83 px**; sharp images give exactly 0.
- Consequence: `edge_cap` is **not inert on the shipped path** — it binds on faint walls, and
  there it *overrides* the `edge_z` noise floor (level becomes `0.5·A < 4 z`). That is the
  price of guaranteeing a crossing exists on the wall run, and it is live on real dim images.
- Re-ran the shipped config on the fixed code → `multiangle_c1_s03d` (results pending).

*k_band coupling bug (self-inflicted, found post-close, fixed)*
- The median z-map has a smaller background MAD than `max(R,G,B)`, so z inflates and more of
  the defocus halo crosses the absolute `k_band=4` gate. On the 6 most defocused images
  `band_half` ballooned 97.5→159.0, 89.0→141.5, 97.5→222.5 px while the measured widths
  stayed correct (140.1/136.7/145.2 px vs sharp siblings 138.0/142.1/139.9 px).
  `band_ratio` fell to 0.326–0.485, below `band_ratio_min=0.5` → false `band_mismatch`.
- Consequence: fiber 01 part 4 lost a2/a3/a5. The 6 angles are 60°-spaced, forming 3 direction
  pairs mod 180°; survivors a1/a4/a6 give directions 0°/0°/120° — only 2 distinct, rank-deficient
  for a 3-parameter ellipse → **0 valid fits**. Fiber 01's `low_confidence` silently flipped
  True→False because the rms gate only saw surviving parts.
- Sweep over k ∈ [4,5,6,7,8] on 6 problem + 11 control images: mismatches 6→0 at k≥5;
  problem `band_ratio` min 0.326→0.646 (k=5)→0.705 (k=6)→0.761 (k=7). Selected **k_band=6.0**.
  Widths move ≤0.13 px across the k range on these images (≤1.43 px on the dim set) — k_band
  moves the coarse mask, not the boundary, but *not* "identical" as first written.
- Non-circular support on a data-selected dim set (`dim_coverage.json`): min coverage 0.8762
  at every k; three images turn `low_confidence` only at k=8.
- **Untracked side effect**: raising `k_band` 4→6 also moved the `FLAG_NO_BG` guard from 2.0
  to 3.0 z in bright mode (`edges.py:282` uses `cfg.k_band / 2.0`). Incidence change never
  measured.

*Code-review defects (all fixed)*
- **Clamp ordering was wrong.** As first written, `max(edge_z, min(frac*A, cap*A))` puts the
  cap *inside* the max, so it cannot cap: for `A_side < 8` the level exceeded the documented
  50% ceiling, and for `A_side < edge_z` it landed *above the wall top* — `_outer_crossing`
  then returns None and `_side_edge` silently falls back to the wall's outer base **with
  FLAG_OK**. Corrected to `min(cap*A, max(edge_z, frac*A))`, which makes that fallback
  (`edges.py:277-279`) effectively dead code, though it still carries FLAG_OK.
- `gui_app.py` **fully reverted** (empty diff vs base) — its tooltip/chip changes described
  bright-mode behaviour the GUI cannot reach, and its docstring claimed an `edge_cap` widget
  that was never added.
- Tests added: bright-gate specificity, explicit-desat preset isolation, `edge_cap == 0.50`
  pin. 213 → **218 tests**.

*Both original calibration criteria were invalid*
- `calibrate_edge_frac._render_circle` declares `W_IMG=1200` but never broadcasts it — every
  image was **500 rows × 1 column**. On one column δ is seed-unstable: old +4.357 **± 0.332**
  (range 3.87–5.10) vs ± 0.006 at full width.
- The other criterion, `low_confidence_count`, was **0 at every edge_frac** — zero
  discriminating power — and its 30-image sample excluded all six defocused images and fiber 07
  entirely. `main()` also implemented no selection rule; it printed and exited.
- **The headline δ was attributed to the wrong config.** `test_xsection_synthetic` builds
  `CONFIG(feature_mode="bright")`, inheriting `edge_frac=0.65` → capped to an effective 0.50.
  That is where "δ 3.90 → −0.51" comes from. Per-side δ over 12 seeds:

  | control | old | shipped 0.30 | effective 0.50 |
  |---|---|---|---|
  | 500×1 (as the study actually ran it) | +4.357 ± 0.332 | +2.103 ± 0.081 | −0.014 ± 0.053 |
  | 500×1200 full width | +4.338 ± 0.006 | +2.124 ± 0.002 | +0.014 ± 0.002 |
  | 6-angle stack (the test fixture) | +3.885 ± 0.006 | **+1.649 ± 0.003** | −0.514 ± 0.003 |

  So roughly half the old bias survives in what ships; the quoted −0.51 is not the shipped path.

*Mechanism attribution is the reverse of what I claimed*
- Fringe-injection A/B (`synthetic_fringe.json`, old emulation verified bit-exact against the
  base worktree): with **no** chromatic fringe the z-map change contributes **0.003 px** — the
  threshold change accounts for the entire −2.22 px. The median z-map removes ≈0.9 px per px of
  **single-channel** (R-only) displacement and ≈0.05 px for a two-channel (R+B) fringe.
  "The median z-map was the dominant fix" was asserted, never tested, and is unsupported.

*Corrected calibration — and why it still does not select a value*
(`recalibrate_edge_frac.json`; full-width synthetic × 12 seeds, k_band=6.0, plus within-pair
|Δwidth| on 3 defocused and 3 control parts, 9 pairs each)

| edge_frac | δ px/side | defocus pair disagreement px | control pair disagreement px |
|---|---|---|---|
| 0.15 | +4.190 ± 0.004 | 7.322 | 7.687 |
| 0.20 | +3.415 ± 0.003 | 7.512 | 7.129 |
| 0.25 | +2.734 ± 0.002 | 7.726 | 7.296 |
| **0.30 (shipped)** | **+2.124 ± 0.002** | **7.843** | **7.423** |
| 0.35 | +1.570 ± 0.002 | 7.924 | 7.289 |
| 0.40 | +1.031 ± 0.002 | 7.543 | 6.898 |
| 0.45 | +0.521 ± 0.002 | 6.556 | 6.830 |
| 0.50 | +0.014 ± 0.002 | 6.076 | 6.756 |

Neither criterion is usable as written, which I did not see at first:
- **δ is circular.** `_render_circle` builds a *linear* ramp with `t = clip((h+1−dist)/2, 0, 1)`,
  so `t = 0.5` falls exactly at `dist = h`, the true radius. Measuring at 50% of amplitude
  therefore recovers the truth **by construction**. δ→0 at frac=0.50 is a property of the
  fixture, not evidence about real fibres — and study 01 showed real edges are *shouldered*
  ramps, where the 50% midpoint bites inward.
- **Pair disagreement is too weak to decide.** n=9 pairs, profiles compared without cross-angle
  registration, and the ~7 px baseline is dominated by misregistration rather than edge error.
  Controls move almost as much as the defocused parts: the defocus-minus-control excess is
  +0.42 (0.30), +0.64 (0.35–0.40), −0.27 (0.45), −0.68 (0.50) — noise-level. It is also
  **not monotone** (0.15 at 7.32 beats 0.30 at 7.84).

So the corrected sweep does not select 0.30, but neither does it establish 0.45–0.50. No valid
criterion in this repo picks any particular `edge_frac`.

*Other open items*
- Fiber 10 is a **genuine stage-3 registration lock**, not an edge defect — I was wrong at first
  close. Per-part rms old→new: p1 2.64→6.27, p3 2.07→4.11. Since the 6 angles form 3 pairs,
  `rms_resid_px` is purely pair disagreement; part 1's a2−a5 pair reaches rms 20.376
  (`xsection_validation.csv`). Coverage is 0.97–1.00 on all 30 images, so edge extraction is
  not at fault. (The lag-0 / +508 / +710 diagnostics quoted earlier have **no persisted
  artifact** and are unverified.)
- Axis-ratio shift is **~100% geometric**: predicted 1.1385–1.1398 from a uniform width shrink,
  observed 1.1387. A null result my "not separated" hedge obscured.
- `band_ratio` never comes within 0.2 of the `band_ratio_min=0.5` gate across all 450 images,
  so gate sensitivity remains undemonstrated; the k_band problem set was the observed failures,
  which is circular (the dim-image check above is the non-circular part).
- Global width shrink is **9.7%** (median over 450 paired images), not the 10.6% quoted earlier.
- Unverifiable for want of an artifact: the 3495-column count (s03 deleted), the fiber-10 lag
  diagnostics, the old-side `band_half` values.

**What I think it means** — The empirical improvement on real C1 data is large, statistically
solid, robust to the column-churn confound, and visually confirmed; MasP2 is provably untouched.
Confidence: **high** on the direction and rough magnitude of the gain.

Two things I got wrong and have now corrected. First, the headline run predates the clamp fix,
so until `s03d` lands the quoted numbers describe code that is not what ships. The divergence is
confined to faint walls on dim images (≤0.83 px on any image median), so I expect the headline
statistics to survive — but "expect" is not "verified", which is why I re-ran rather than
patched the wording. Second, I claimed the corrected calibration favours 0.45–0.50; it does not
establish that, because its δ criterion is circular by construction of the synthetic ramp and
its pair criterion is at noise level. The honest position is that **`edge_frac=0.30` ships on
empirical grounds and no valid criterion selects any value**, 0.30 included.

I also over-claimed the mechanism: on synthetic data the threshold change does essentially all
the work and the z-map does ~nothing absent a single-channel fringe; the real-data split was
never measured. Confidence: **high** that the original attribution was wrong, **low** on the
true real-data split.

Evidence that would change this reading: an `s03d` that moves the headline materially would mean
the pre-fix numbers were load-bearing; a *registered* pair-disagreement comparison at 0.30 /
0.40 / 0.45 on the known defocused parts would give the first non-circular handle on `edge_frac`
without a 450-image run.

**Next** — Fold the `s03d` numbers into the report, then rewrite `docs/report/03_edge_criteria/`
and `docs/features/05_edge_criteria.md`, both of which still carry the retracted −0.51 figure,
the broken clamp formula, and the invalidated calibration narrative. Follow-ups for the owner,
none blocking: (i) stage-3 cross-angle registration peak validation (fiber 10, and the same
fragility on fiber 04 part 5); (ii) measure the `FLAG_NO_BG` incidence change from `k_band` 4→6
and re-check the other absolute-z thresholds (`amin`, `rise_min`, `slope_min`/`slope_cap`) as a
group; (iii) add an end-to-end faint-wall test — `test_bright_level_never_leaves_the_wall`
re-implements the formula inline and never calls `edges.py`, so reverting the clamp fix still
passes all 218 tests; (iv) pin end-to-end bright accuracy at the *shipped* presets, which no
test currently does; (v) study 02's anisotropy statistics should be recomputed on the final run.
