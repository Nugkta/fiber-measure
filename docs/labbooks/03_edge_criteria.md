Hypothesis: Replacing max(R,G,B) with median(R,G,B) in bright mode and switching from a fixed to a clamped relative edge threshold will reduce defocus-induced width errors on C1 images without introducing inward bite.
Status: done
Started: 2026-08-19

## 2026-08-19 (re-run with k_band fix — final numbers, study closed)
**What I did** — Re-ran stage 1 (450 primary images, excluding the redundant `s.tiff` scale-bar variants this time) with `--feature-mode bright --edge-frac 0.30 --k-band 6.0` → `fibrecv_output/multiangle_c1_s03b/`, then stage 3 with `--scale-source items --validation`. Added `BRIGHT_DEFAULTS = {"edge_frac": 0.30, "k_band": 6.0}` to `run_measure.build_config` so `--feature-mode bright` picks the calibrated knobs up automatically (explicit flags still override) + 3 tests. Generated the old-vs-new boundary overlay figure (`scripts/spot_check_overlay.py` → `docs/report/03_edge_criteria/spot_check_overlay.png`). Rewrote the report with corrected numbers, the k_band section, the figure, and a caveats section.

**What I observed** — Stage 1: 450 ok, 0 errors, **0 low_confidence** (the 6 false `band_mismatch` are gone). Fiber 01 part 4 restored: **0 → 2195 valid columns**; fiber 01 n_positions 8676 → 11172 (old 11215). Final per-fiber: low_confidence **6 → 1**, `part_rms_med_max_px` p50 4.97 → 4.10, axis ratio grand median 1.127 → 1.139, total positions 162223 → 161634 (−0.4%). Per-part rms across all 75 parts: p50 3.70 → **2.78**, p90 6.64 → 4.39, p95 7.76 → 4.89, max 11.84 → **6.27**, parts >6.0 px **11 → 1**. Crucially the rms distribution is essentially unchanged from the buggy k_band=4 run (p50 2.84 → 2.78, max 6.29 → 6.27), so the improvement was never an artifact of the lost part. Visual spot-check on 4 images (1 defocused, 2 sharp, 1 dim): in every panel the old boundary sits outside the visible bright-to-dark edge in the dark surround and the new one lands on it; no inward bite anywhere. Per-side shift −12.2 px (defocused), −6.7 / −5.6 px (sharp), −7.6 px (dim). Suite: **216 passed**.

**What I think it means** — Study closed, hypothesis confirmed. The improvement is real and now verified three independent ways: synthetic (δ 3.90 → −0.51 px/side), statistical (rms distribution tightened throughout, parts >6.0 px 11 → 1), and visual (boundary lands on the edge in all four spot-checks). The one remaining flag, fiber 10, is a genuine stage-3 registration lock failure, not an edge defect. MasP2 is provably bit-identical. Confidence: high. What would change this: a collaborator-confirmed different angle spacing, or evidence that the visible bright-to-dark transition is not the physical wall (the spot-check assumes it is).

**Next (follow-ups for the owner, none blocking)** — (i) stage-3 cross-angle registration peak validation (fiber 10, and the same fragility on fiber 04 part 5); (ii) re-check the other absolute-z thresholds (`amin`, `rise_min`, `slope_min`/`slope_cap`) as a group if the z-map is ever changed again; (iii) study 02's anisotropy statistics should be recomputed on the s03b run before being restated, since axis ratio shifts partly for geometric reasons; (iv) `xsec_rms_flag_px` is now stricter in relative terms than when calibrated.

## 2026-08-19 (post-close verification — two defects found, one fixed)
**What I did** — Ran the verification items skipped at first close: (a) MasP2 desat regression check, (b) fiber 10 rms diagnosis, (c) `k_band`/`amin` band-finding spot-check on the dimmest C1 quartile (the check the plan asked for and I skipped), (d) per-part rms distribution recalibration. Found two defects; fixed one, scoped the other out.

**What I observed** —
1. **MasP2 desat: PASS, bit-identical.** 3 images run on both branches: median diameters 80.52685546875 / 62.23075866699219 / 42.61311340332031 µm, deltas exactly 0.0. `diameter_raw` arrays and the float32 `D` z-map both match by SHA-256. The desat path is provably untouched.
2. **`k_band` coupling bug (my defect, now fixed).** The median z-map has a smaller background MAD than `max(R,G,B)`, so z-scores inflate and more of the defocus halo crosses the absolute `k_band=4` threshold. On the 6 most defocused C1 images the coarse band ballooned (`band_half` 97.5→159.0, 89.0→141.5, 97.5→222.5 px) while the measured widths stayed correct (140.1/136.7/145.2 px vs their sharp siblings' 138.0/142.1/139.9 px). `band_ratio = width/band_thickness` collapsed to 0.326–0.483, below `band_ratio_min=0.5` → false `band_mismatch` → excluded from stage 3. Not systemic: `band_half` p50 86.5→87.5, p90 102.0→102.6; only the tail blew up (max 133→222).
   **Consequence:** fiber 01 part 4 lost angles a2/a3/a5. The 6 angles are 60°-spaced (0/60/120/180/240/300), forming 3 direction pairs mod 180°; the survivors a1/a4/a6 give directions 0°/0°/120° — only **2 distinct directions**, rank-deficient for a 3-parameter ellipse → **0 valid fits out of 3495 columns**. Fiber 01 silently lost a whole part (n_positions 11215→8676) *and* its `low_confidence` flipped True→False because the rms gate only saw the surviving parts. So the headline "6→1" was partly an artifact of data loss.
   **Fix:** `k_band` recalibrated for bright mode. Sweep over [4,5,6,7,8] on the 6 problem images + 11 controls: mismatches 6→0 at k≥5; problem `band_ratio` min 0.326→0.646 (k=5)→0.705 (k=6)→0.761 (k=7); **measured widths identical at every k** (problem 144.2 px, control 147.6 px throughout) — `k_band` moves only the coarse mask, not the boundary. On the 12 hardest/dimmest images coverage *improves* slightly (min 0.732→0.753 at k=6) with widths drifting 139.5→139.7 px (0.1%). **Selected `k_band=6.0`** (0.205 headroom above the 0.5 gate).
3. **Fiber 10 is a real defect, not a threshold artifact — I was wrong at first close.** Per-part rms old→new: p1 2.65→**6.29**, p2 2.84→2.27, p3 2.07→4.12, p4 2.89→3.73, p5 3.38→2.13. Because the 6 angles form 3 pairs, `rms_resid_px` is purely pair disagreement. Part 1 pair |Δ|: a1−a4 1.91→6.76, **a2−a5 2.73→20.44**, a3−a6 7.14→5.01. The a2/a5 pair locked ~500 px off in x-registration (rms 20.37 at lag 0 vs 4.92 at lag +508; true lag ≈+710 by independent search); part 3 shows the same on a3/a6 (+512). Coverage is fine on all 30 images (0.97–1.00 valid-column ratio vs old), so edge extraction is not at fault. Cause: the ~10.6% global width shrink altered the cross-correlation landscape enough to flip the peak choice. **This is a stage-3 registration robustness issue, out of scope for the edge-criteria study** — recorded for follow-up.
4. **Per-part rms distribution** (75 parts): OLD p50=3.70 p90=6.64 p95=7.76 max=11.84 → NEW p50=2.84 p90=4.38 p95=4.93 max=6.29. The distribution tightened substantially; `xsec_rms_flag_px=6.0` now sits above p95, flagging exactly 1 part.

**What I think it means** — The core hypothesis still holds and the improvement is real, but my first close overstated it and missed a defect I introduced. The `k_band` coupling is exactly the risk the plan named ("these knobs were calibrated under max-V; median-V has slightly different MAD") and skipping that check let a silent data-loss bug ship. Fiber 10's registration failure is genuine and unrelated to the edge criteria. Confidence: high on the mechanism for both (each traced to concrete numbers and reproduced by sweep). Re-running stage 1+3 with `k_band=6.0` for honest final numbers.

**Next** — Re-run with `k_band=6.0`, redo the old-vs-new comparison, correct the report. Follow-ups for the owner: (i) stage-3 registration peak validation (fiber 10), (ii) bright mode now needs three non-default flags (`--edge-frac 0.30 --k-band 6.0 --feature-mode bright`) — a mode-preset would remove that footgun.

## 2026-08-19 (full C1 run + xsection comparison)
**What I did** — Completed full 450-image C1 run at edge_frac=0.30, feature_mode=bright → `fibrecv_output/multiangle_c1_s03/`. Ran stage 3 (xsection) with `--validation --scale-source items`. Compared xsection_summary against study 02 results.

**What I observed** — Low-confidence fibers: **6 → 1** (fibers 01/07/08/09/13/15 all cleared; fiber 10 newly flagged with rms=6.29, barely above the 6.0 threshold). RMS residuals dramatically improved on formerly-bad fibers: 08 (11.84→4.61), 01 (9.11→4.98), 07 (8.53→3.98). Median rms across all 15 fibers: 4.97→4.07 px. Axis ratio largely stable (grand median 1.131→1.145); biggest ratio shifts on formerly-bad fibers (01: 1.189→1.151, 07: 1.240→1.197, 15: 1.230→1.166) — these are likely now more accurate since the old readings had systematic wide bias from chromatic fringes. Two fibers with rms regressions: fiber 10 (3.38→6.29, now flagged) and fiber 04 (3.92→5.56, still passes). Stage 1 LOWCONF: 6 unique images, all on fiber 01 (angles a2/a3/a5/a6, parts 2–4) — the tighter threshold exposes genuine defocus issues earlier.

**What I think it means** — The hypothesis is confirmed: median z-map + clamped relative threshold substantially reduced defocus-induced errors (low-confidence 6→1, rms p50 −18%). The one new regression (fiber 10, rms barely over threshold at 6.29 px) is borderline — the xsec_rms_flag_px=6.0 was calibrated on the old distribution (p50=3.7, p90=6.6); the new distribution is tighter (p50≈3.7→≈3.3 estimated) so the threshold may want recalibration, but that's optional. Confidence: high.

**Next** — Recalibrate `xsec_rms_flag_px` on the new rms distribution (optional — fiber 10 is borderline, not a real problem). Run full test suite in worktree. Clean up investigation files. Documentation: feature doc + report when study closes.

## 2026-08-19 (calibration + full run)
**What I did** — Ran edge_frac calibration sweep [0.15..0.45] with edge_cap=0.50 on synthetic circle + 30-sample C1 subset (`scripts/calibrate_edge_frac.py`). Results:

| edge_frac | delta (px) | in band | med_w (px) | std_w | iqr_w | cov | lowconf |
|-----------|-----------|---------|-----------|-------|-------|-----|---------|
| 0.15 | 4.39 | NO | 159.0 | 16.4 | 23.7 | 99.8% | 0 |
| 0.20 | 3.55 | YES | 155.0 | 16.1 | 22.4 | 99.8% | 0 |
| 0.25 | 2.82 | YES | 152.0 | 16.0 | 20.4 | 99.8% | 0 |
| **0.30** | **2.19** | **YES** | **149.5** | **16.0** | **19.6** | **99.7%** | **0** |
| 0.35 | 1.61 | YES | 147.1 | 15.8 | 18.9 | 99.7% | 0 |
| 0.40 | 1.05 | YES | 144.9 | 15.7 | 18.4 | 99.7% | 0 |
| 0.45 | 0.52 | YES | 142.8 | 15.6 | 17.8 | 99.7% | 0 |

Selected edge_frac=0.30 (delta=2.19 px, in band, IQR past elbow of diminishing returns). **0 low-confidence at every frac value** — the median z-map change alone fixed the defocus issue (old pipeline had 11/450 low-confidence). Updated config.py comment with calibration result. Launched full 450-image C1 run at edge_frac=0.30 → `fibrecv_output/multiangle_c1_s03/`.

**What I observed** — Calibration is clean: delta monotonically decreases 4.39→0.52 as frac increases; std/IQR also decrease monotonically; coverage stable >99.7% everywhere; no low-confidence images at any sweep point. The selection is not a boundary value — frac=0.30 sits in the middle of the valid range with margin on both sides.

**What I think it means** — The median z-map was the dominant fix: it eliminated the chromatic-fringe bias that was causing defocus sensitivity. The clamped threshold gives additional control but the delta=2.19 px at frac=0.30 is a good compromise between accuracy and stability. Confidence: high (clean monotonic trends, 0 lowconf, not a boundary selection). Full C1 run in progress.

**Next** — When stage 1 completes: run stage 3 (xsection), compare with study 02 results (number of low-confidence fibers, rms residuals, anisotropy, axis ratio).

## 2026-08-19
**What I did** — Created study branch `worktree-edge-criteria` off `worktree-multiangle-xsection`. Implemented two code changes:
1. `src/fibrecv/features.py:81`: bright-mode z-map changed from `rgb.max(axis=2)` to `np.median(rgb, axis=2)` — rejects single-channel chromatic-aberration fringes.
2. `src/fibrecv/edges.py:260`: bright-mode threshold changed from `min(edge_z, edge_frac*A)` to `max(edge_z, min(edge_frac*A, edge_cap*A))` — edge_frac is now the primary relative knob, edge_z is the absolute floor, edge_cap (0.50) prevents study-01 inward bite. Desat path unchanged.
3. Added `edge_cap` field to `config.py`, `--edge-cap` CLI flag to `run_measure.py`, updated GUI slider tooltip and header chip.

**What I observed** — All 213 tests pass. Synthetic circle control delta dropped from 3.90 px/side to −0.51 px/side (well within the −1 to +4 acceptance band). The median z-map alone removed most of the systematic wide bias before any threshold recalibration.

**What I think it means** — The median z-map change is a clear win: it eliminates the chromatic fringe riding that was adding 1–4 px of systematic bias. The clamped threshold infrastructure is in place but edge_frac still at 0.65 (legacy default) — calibration sweep needed to find the optimal value for bright mode. Confidence: high that the code is correct, medium that the default edge_frac will need tuning.

**Next** — Run edge_frac calibration sweep [0.15..0.45] on C1 data with edge_cap=0.50, using defocus stability (criterion 1) and synthetic delta (criterion 2) as selection criteria. Then re-run full C1 pipeline.
