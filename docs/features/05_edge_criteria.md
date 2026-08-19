# Edge criteria upgrade (median RGB + clamped relative threshold)
Added: 2026-08-19
Code: `src/fibrecv/features.py`, `src/fibrecv/edges.py`, `src/fibrecv/config.py`, `src/fibrecv/run_measure.py`

## What it does
Replaces the bright-mode z-map from `max(R,G,B)` to `median(R,G,B)` and changes the edge threshold formula from a fixed absolute level to a clamped relative threshold. This eliminates the systematic 1-4 px wide bias caused by single-channel chromatic-aberration fringes riding the max channel, and makes the edge crossing scale with wall amplitude while preventing convergence on 50% (which causes inward bite on cylindrical fibres). Because the median z-map shifts the z-scale, `k_band` is recalibrated for bright mode as part of the same change. The desat mode (MasP2) path is bit-identical (verified by SHA-256 on the `D` z-map and `diameter_raw`).

## Design choices
- **Median over max**: median of 3 channels rejects the single outlier channel where chromatic aberration concentrates. Rejected alternatives: min (deletes faint-wall cap, crashes 46% of MasP2 defocused walls per audit), mean (still pulled by the outlier), max (the original, rides the fringe).
- **Two-sided clamp, not simple relative**: `max(edge_z, min(edge_frac*A, edge_cap*A))` gives three independent controls. `edge_z` is the absolute floor protecting faint walls from going below noise. `edge_frac` is the primary relative knob. `edge_cap` prevents convergence on 50% where study 01 showed inward bite. Rejected: unconstrained calibration (converges on 50%), single relative fraction (no floor for faint walls).
- **Mode-gated threshold**: the new formula only applies when `feature_mode == "bright"`. Desat mode keeps the original `min(edge_z, edge_frac*A)` formula unchanged, so MasP2 behaviour is frozen.
- **edge_frac=0.30 calibrated**: swept [0.15..0.45] on synthetic circle + 30-sample C1 subset. At 0.30: delta=2.19 px (in [-1,+4] acceptance band), 0 low-confidence, IQR past diminishing-returns elbow.
- **k_band=6.0 for bright mode**: the median z-map has a smaller background MAD, so z-scores inflate and defocus halo crosses the absolute `k_band=4.0`, ballooning the coarse band mask and raising false `band_mismatch`. Swept [4..8]: mismatches 6→0 at k≥5, measured widths identical at every k (k_band moves only the coarse mask, never the boundary), dim-image coverage slightly improved. 6.0 gives 0.205 headroom above `band_ratio_min`. Rejected: lowering `band_ratio_min` (the healthy population sits at p05=0.767 while the false positives were 0.326–0.483, so no threshold separates them without disabling the check).
- **Presets in the CLI, not the dataclass**: `CONFIG` defaults stay desat-calibrated so MasP2 cannot regress; `--feature-mode bright` applies `BRIGHT_DEFAULTS` in `run_measure.build_config`, and explicit flags still override. Rejected: changing the dataclass defaults (would silently alter MasP2), and per-mode dataclass factories (larger API change than the study warranted).

## Algorithm details
In `features.py`, bright-mode z-map computation (line 81):
```
V = np.median(rgb, axis=2).astype(np.float32)
```
The median of 3 values equals the middle value, rejecting the single channel where chromatic aberration concentrates.

In `edges.py`, the threshold level computation (`_side_edge`, line 260):
```
if cfg.feature_mode == "bright":
    level = base_val + max(cfg.edge_z, min(cfg.edge_frac * A_side, cfg.edge_cap * A_side))
else:
    level = base_val + min(cfg.edge_z, cfg.edge_frac * A_side)
```

For C1 bright images (typical A_side 40-70 z):
- At frac=0.30, cap=0.50: level = base + max(4, min(0.30*50, 0.50*50)) = base + max(4, 15) = base + 15 z
- For a hypothetical faint wall (A_side=5): level = base + max(4, min(1.5, 2.5)) = base + 4 z (floor wins)

Config fields in `config.py`: `edge_frac` (dataclass default 0.65 = desat legacy), `edge_cap` (0.50), `edge_z` (4.0), `k_band` (dataclass default 4.0 = desat). CLI flags: `--edge-frac`, `--edge-cap`, `--edge-z`, `--k-band`.

Bright-mode presets live in `run_measure.py`:
```
BRIGHT_DEFAULTS = {"edge_frac": 0.30, "k_band": 6.0}
```
applied in `build_config` when `--feature-mode bright` is given, before explicit flag overrides.

## Caveats
- **The presets are CLI-only.** Code that constructs `CONFIG(feature_mode="bright")` directly (tests, notebooks, any future GUI bright mode) gets the desat-calibrated `edge_frac=0.65` and `k_band=4.0`, which are both wrong for bright images — `k_band=4.0` in particular reintroduces the false `band_mismatch` on defocused images. Set both explicitly in that case.
- **Other absolute-z thresholds were not individually re-swept.** `amin=3.0`, `rise_min=2.0`, `slope_min`, `slope_cap` are all absolute in z units and carry the same coupling to the z-scale that `k_band` did. They showed no symptoms across 450 C1 images, but a future z-map change should re-check them as a group.
- `xsec_rms_flag_px=6.0` was calibrated against the old, wider rms distribution (p50=3.70, p90=6.64). Under the new distribution (p50=2.78, p90=4.39) it sits above p95 rather than near p90 — stricter in relative terms, though it still isolates exactly the one genuinely defective part.
- `np.median` on axis=2 with 3 channels returns the middle value as float64; the `.astype(np.float32)` cast is required.
- The width shrink is not uniform across images (−5.6 to −12.2 px/side in the spot-check), so any downstream quantity that is a ratio of widths — notably `axis_ratio` — shifts partly for geometric rather than physical reasons.
