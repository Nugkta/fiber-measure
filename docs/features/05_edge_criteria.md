# Edge criteria upgrade (median RGB + clamped relative threshold)
Added: 2026-08-19
Code: `src/fibrecv/features.py`, `src/fibrecv/edges.py`, `src/fibrecv/config.py`, `src/fibrecv/run_measure.py`

## What it does
Replaces the bright-mode z-map from `max(R,G,B)` to `median(R,G,B)` and changes the edge threshold formula from a fixed absolute level to a clamped relative threshold. The median z-map rejects single-channel chromatic-aberration fringes (synthetic-verified contribution ~0.9 px per px of single-channel displacement; real-data contribution not isolated). The clamped relative threshold makes the edge crossing scale with wall amplitude while preventing convergence on 50% (which causes inward bite on cylindrical fibres). On synthetic data without a chromatic fringe, the threshold change accounts for essentially all of the improvement (~2.22 px) and the z-map contributes ~0.003 px. Because the median z-map shifts the z-scale (smaller background MAD → inflated z-scores), `k_band` is recalibrated for bright mode as part of the same change. The desat mode (MasP2) path is bit-identical (verified by SHA-256 on the `D` z-map and `diameter_raw`).

## Design choices
- **Median over max**: median of 3 channels selects the middle-ranked channel, rejecting a single-channel outlier where chromatic aberration concentrates. Against a two-channel fringe the median rides it too (~0.05 px per px of displacement). Rejected alternatives: min (halves the faint-wall guard, degrading coverage on defocused images in both modes — per early-stage audit, though no artifact was saved), mean (still pulled by the outlier), max (the original, rides the fringe).
- **Two-sided clamp, not simple relative**: `min(edge_cap*A, max(edge_z, edge_frac*A))` gives three independent controls. `edge_z` is the absolute floor protecting faint walls from going below noise. `edge_frac` is the primary relative knob. `edge_cap` prevents convergence on 50% where study 01 showed inward bite, and guarantees a crossing always exists on the wall run (`level ≤ base + cap*A < wall top`). On faint walls (`A_side < edge_z/edge_cap = 8 z`), the cap overrides the `edge_z` floor — the price of the crossing guarantee. Rejected: unconstrained calibration (converges on 50%), single relative fraction (no floor for faint walls).
- **Mode-gated threshold**: the new formula only applies when `feature_mode == "bright"`. Desat mode keeps the original `min(edge_z, edge_frac*A)` formula unchanged, so MasP2 behaviour is frozen.
- **edge_frac=0.30 selected empirically**: the original calibration sweep (0.15–0.45) used two selection criteria that were both invalid — the synthetic circle rendered at 500×1 pixels (a broadcast bug), and `low_confidence_count` was 0 at every value (zero discriminating power). A corrected calibration also cannot rescue the choice: the synthetic δ criterion is circular (the ramp's 50% point equals the true radius by construction, so δ→0 at frac=0.50 is a property of the fixture, not a fact about real fibres), and the pair-disagreement criterion sits at noise level (n=9, unregistered profiles, controls move almost as much as defocused parts). No valid criterion selects any specific `edge_frac`; 0.30 ships because the full 450-image C1 run at that value shows large, statistically significant improvement (Wilcoxon p=2e−7, 60/75 parts improved, bootstrap 95% CI [−1.58, −0.57] px on per-part rms change).
- **k_band=6.0 for bright mode**: the median z-map has a smaller background MAD, so z-scores inflate and defocus halo crosses the absolute `k_band=4.0`, ballooning the coarse band mask and raising false `band_mismatch`. Swept [4..8]: mismatches 6→0 at k≥5, widths move ≤0.13 px on the calibration set (≤1.43 px on the dimmest images), dim-image coverage slightly improved. 6.0 gives 0.205 headroom above `band_ratio_min`. Rejected: lowering `band_ratio_min` (the false positives sat at 0.326–0.485, below the gate, while the rest of the population stays well above 0.5, so no alternative threshold separates them).
- **Presets in the CLI, not the dataclass**: `CONFIG` defaults stay desat-calibrated so MasP2 cannot regress; `--feature-mode bright` applies `BRIGHT_DEFAULTS` in `run_measure.build_config`, and explicit flags still override. Rejected: changing the dataclass defaults (would silently alter MasP2), and per-mode dataclass factories (larger API change than the study warranted).

## Algorithm details
In `features.py`, bright-mode z-map computation (line 83):
```
V = np.median(rgb, axis=2).astype(np.float32)
```
The median of 3 values equals the middle value, rejecting a single-channel outlier.

In `edges.py`, the threshold level computation (`_side_edge`, lines 263–273):
```
if cfg.feature_mode == "bright":
    level = base_val + min(cfg.edge_cap * A_side,
                           max(cfg.edge_z, cfg.edge_frac * A_side))
else:
    level = base_val + min(cfg.edge_z, cfg.edge_frac * A_side)
```

The cap must be **outermost** (`min(cap*A, max(…))`). With it inside the max (`max(edge_z, min(frac*A, cap*A))`), a faint wall where `A_side < edge_z` gets a level above the wall top — `_outer_crossing` returns None and `_side_edge` silently falls back to the wall's outer base with FLAG_OK.

For C1 bright images (typical A_side 40–70 z):
- At frac=0.30, cap=0.50: level = base + min(0.50×50, max(4, 0.30×50)) = base + min(25, 15) = base + **15 z**

For a faint wall (A_side=5, below the `edge_z/edge_cap = 8 z` crossover):
- level = base + min(0.50×5, max(4, 0.30×5)) = base + min(2.5, max(4, 1.5)) = base + min(2.5, 4) = base + **2.5 z** (cap wins, overriding the noise floor). This is by design: it guarantees a crossing exists on the wall run, at the cost of measuring at a sub-noise level on the faintest walls.

Config fields in `config.py`: `edge_frac` (dataclass default 0.65 = desat legacy), `edge_cap` (0.50), `edge_z` (4.0), `k_band` (dataclass default 4.0 = desat). CLI flags: `--edge-frac`, `--edge-cap`, `--edge-z`, `--k-band`.

Bright-mode presets live in `run_measure.py`:
```
BRIGHT_DEFAULTS = {"edge_frac": 0.30, "k_band": 6.0}
```
applied in `build_config` when `--feature-mode bright` is given, before explicit flag overrides. The GUI (`gui_app.py`) exposes the same presets through a sidebar `Image mode` selectbox: applying a mode switch starts from that mode's calibrated defaults (bright: `BRIGHT_DEFAULTS`; desat: dataclass) and resets the visible `edge_frac` widget accordingly, while non-coupled knobs survive.

## Caveats
- **The presets apply only through the two front ends (CLI `build_config`, GUI mode selector).** Code that constructs `CONFIG(feature_mode="bright")` directly (tests, notebooks, scripts) gets the desat-calibrated `edge_frac=0.65` and `k_band=4.0`, which are both wrong for bright images — `k_band=4.0` in particular reintroduces the false `band_mismatch` on defocused images. Set both explicitly in that case.
- **`edge_cap` is not inert.** At the shipped `edge_frac=0.30`, the cap binds whenever `A_side < edge_z/edge_cap = 8 z`. On dim C1 images this affects ~800–1500 of ~2500 columns per image, shifting medians up to ~0.83 px. The cap overrides the `edge_z` noise floor on these faint walls — the crossing is guaranteed to exist, but the measurement sits at a sub-noise level. A follow-up could flag such columns rather than silently measuring them.
- **Other absolute-z thresholds were not individually re-swept.** `amin=3.0`, `rise_min=2.0`, `slope_min`, `slope_cap` are all absolute in z units and carry the same coupling to the z-scale that `k_band` did. They showed no symptoms across 450 C1 images, but a future z-map change should re-check them as a group.
- **`k_band` raise moves the NO_BG guard.** The `FLAG_NO_BG` threshold at `edges.py:282` is `cfg.k_band / 2.0`, so raising `k_band` from 4.0 to 6.0 moved it from 2.0 to 3.0 z in bright mode. The incidence change was never measured.
- `xsec_rms_flag_px=6.0` was calibrated against the old, wider rms distribution (p50=3.70, p90=6.64). Under the new distribution (p50=2.74, p90=4.26) it sits above p95 — it now flags only more extreme relative outliers, though it still isolates exactly the one genuinely defective part (fiber 10, a stage-3 registration failure).
- `np.median` on axis=2 with 3 channels returns the middle value as float64; the `.astype(np.float32)` cast is required.
- The width shrink is not uniform across images (−5.6 to −12.2 px/side in the spot-check), so any downstream quantity that is a ratio of widths — notably `axis_ratio` — shifts partly for geometric rather than physical reasons. The observed axis-ratio shift is ~100% geometric.
- **No valid criterion selects `edge_frac`.** The shipped 0.30 rests on empirical real-data performance. The corrected synthetic-δ criterion is circular (favours 0.50 by construction of the fixture) and the pair-disagreement criterion is noise-dominated. A properly registered pair-disagreement comparison or a full-pipeline run at alternative values would give the first non-circular evidence.
