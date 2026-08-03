# ESF Edge Refinement — Design Spec

Date: 2026-08-03
Status: approved in chat, pending spec review
Study: `docs/labbooks/01_esf_edge_consistency.md` (hypothesis + empirical basis)

## Problem

The optical setup blurs the sharp air/silk boundary: measured edge profiles are
S-shaped ramps with Gaussian-blur width σ ≈ 5–15 px, varying between images and
along a single fibre (defocus). The current recipe places each boundary at a
fixed low threshold near the base of the wall ramp
(`base + min(edge_z, edge_frac·A_side)` in `edges._side_edge`). The distance
from ramp base to true edge grows with σ, so the measured diameter carries a
focus-dependent bias — measured at +2.3 / −3.4 / +3.4 px **per side** (S
channel) on three MasP2 images, i.e. up to ~7 px image-to-image diameter drift
for 60–120 px fibres. This is a driver of replicate disagreement.

For any symmetric PSF, the 50% midpoint of the blurred profile coincides with
the true step position (knife-edge principle) and is invariant to σ. Fitting
the erf model to each wall ramp and moving the boundary to the fitted midpoint
removes the focus-dependent bias.

## Goal and non-goals

**Goal:** replicate-to-replicate consistency — the same fibre measured at
different focus should give the same diameter.

**Non-goals (explicitly out of scope):**
- Absolute µm accuracy. The edge is chromatic (channel-dependent midpoints
  2–8 px apart, focus-correlated), so a "true geometric edge" is not
  recoverable to better than a few px; absolute calibration needs a
  known-diameter standard (future study).
- Joint two-edge fitting for thin fibres where the two blur ramps overlap.
  Only if the A/B validation shows elevated failure on thin-fibre groups.
- Acquisition changes (PNG instead of JPEG, backlit illumination) — report
  recommendations only.

## Empirical basis (study 01, 2026-08-03)

- erf fits on the saturation channel: 3.5–6.7% residual/amplitude on
  `masp2 10_1_1 / 10_5_1 / 10_10_1`.
- Channel comparison: S was the only channel passing the <8%-residual gate
  across all three images (33/62, 50/67, 58/63 blocks); linearised R−G and G
  collapse near the iridescent specular rim on the thin-fibre image.
- Fitting is done on the pipeline's `D` z-map, which is an affine rescale of
  S per image — affine transforms are absorbed by the erf's `a, b` parameters,
  so `D` and `S` give identical midpoints.

## Architecture

New module `src/fibrecv/refine.py`, called from `compute.compute_measurement`
between `edges.detect_edges` and `qc.run_qc`:

```
D, S, … = features.rgb_to_desaturation(rgb, cfg)
bnd     = band.locate_band(D, cfg)
edg     = edges.detect_edges(D, bnd, cfg)
edg, ref = refine.refine_edges(D, edg, bnd, cfg)   # NEW (identity if cfg.refine_on is False)
res     = qc.run_qc(edg, bnd, cfg)
```

API:

```python
@dataclass
class RefineResult:
    refined_top: np.ndarray   # bool (W,) column's top edge was moved
    refined_bot: np.ndarray   # bool (W,)
    sigma_top: np.ndarray     # float (W,) fitted blur sigma, NaN where unrefined
    sigma_bot: np.ndarray     # float (W,)
    resid_top: np.ndarray     # float (W,) relative fit residual, NaN where unrefined
    resid_bot: np.ndarray     # float (W,)
    n_blocks: int             # blocks attempted (per side)
    n_pass_top: int           # blocks passing all gates
    n_pass_bot: int

def refine_edges(D, edg, bnd, cfg) -> tuple[EdgeResult, RefineResult]
```

`edges.py` detection logic is untouched. `EdgeResult.flags` is **not** modified:
`qc.run_qc` requires `flags == FLAG_OK` for validity and already reuses bits
32/64 for its own reasons, so unrefined columns must not gain a flag bit —
they keep today's edge and today's validity. Refinement coverage lives only in
`RefineResult` and the meta dict.

## Algorithm

Per side (top wall, bottom wall) independently:

1. **Anchor columns**: columns with finite `y_top`/`y_bot` and
   `flags == FLAG_OK`.
2. **Blocks**: non-overlapping runs of `refine_block` (default 16) columns
   within the band span; a block is attempted if ≥ 70% of its columns are
   anchor columns.
3. **Aligned mean profile**: sample `D[y, x]` at
   `y = anchor[x] + t` for offsets `t ∈ [−refine_out, +inside]` in steps of
   0.5 px (bilinear in y), average over the block's anchor columns.
   `inside = clip(refine_in_frac · band_half, 8, refine_in_max)`.
   Per-column alignment (not row-average) so fibre tilt does not smear the
   ramp. Sign convention: `+t` points into the fibre for both sides.
4. **Fit** the 4-parameter model by least squares
   (`scipy.optimize.curve_fit`):
   `I(t) = a + (b−a)·½·[1 + erf((t − t₀)/(√2·σ))]`,
   init: `a` = mean of outermost 6 samples, `b` = mean of innermost 6,
   `t₀` = level-midpoint sample, `σ` = 3.
5. **Gates** (all must pass, else the block is rejected):
   - relative residual `rms(fit − data)/|b−a| < refine_relmax`
   - `refine_sigma_min ≤ σ ≤ refine_sigma_max`
   - `|t₀| ≤ refine_maxshift`
   - `b − a > 0` (desaturation rises into the fibre)
6. **Offset field**: each passing block yields offset `t₀` at its centre
   column. Per-column offsets are linearly interpolated across x between
   passing block centres; runs of more than `refine_gap_blocks` consecutive
   failed blocks are left unrefined (interpolation does not bridge them), as
   are columns outside the outermost passing centres (no extrapolation).
7. **Apply**: `y_top_new = y_top + o_top(x)` (top edge moves down into the
   fibre for positive t₀), `y_bot_new = y_bot − o_bot(x)`; diameter is
   recomputed as `(y_bot_new − y_top_new)·cos(tilt)` exactly as in
   `detect_edges`. Non-anchor columns keep their NaN/flagged state.

With `cfg.refine_on == False`, `refine_edges` returns the inputs unchanged
plus an all-False `RefineResult` — bit-identical to today's pipeline.

## Config additions (`config.py`)

| Field | Default | Meaning |
|---|---|---|
| `refine_on` | `True` | master switch (also the A/B switch) |
| `refine_block` | 16 | columns per fitting block |
| `refine_out` | 35 | px sampled outside the anchor (background side) |
| `refine_in_frac` | 0.8 | inside extent as fraction of `band_half` |
| `refine_in_max` | 28 | hard cap on inside extent (px), keeps the specular rim out |
| `refine_relmax` | 0.08 | max relative fit residual |
| `refine_sigma_min` | 0.8 | min accepted σ (px) |
| `refine_sigma_max` | 20.0 | max accepted σ (px) |
| `refine_maxshift` | 12.0 | max |t₀| offset from anchor (px) |
| `refine_gap_blocks` | 2 | max consecutive failed blocks bridged by interpolation |

Defaults derive from the study-01 measurements (σ range 5–15 px, observed
offsets ≤ ~6 px/side, residuals 3–7%).

> **Revised after M5 (2026-08-03):** `refine_relmax` was raised 0.08 → **0.15**
> when the full-set validation showed 0.08 refined only 76% of anchor columns,
> under this spec's own ≥80% coverage bar. Rationale, alternatives tried and
> the bound on the admitted midpoint error are in
> `docs/report/01_esf_edge_consistency.md` §3.1. The other nine defaults are
> unchanged.

## Diagnostics

- meta JSON gains a `refine` sub-dict: enabled flag, per-side block pass
  counts, refined-column coverage %, median σ and median |t₀| per side.
- GUI: one checkbox for `refine_on` in the boundary-knobs card; the
  diagnostics figure gains a σ(x) trace (the "focus map") when refinement ran.
  No other GUI changes.

## Performance

~`(W/16)` blocks × 2 sides ≈ 320 fits per 2560-px image, each ≤ ~140 samples
and 4 parameters: well under 1 s, negligible next to the existing pipeline.

## Validation protocol (main experiment of study 01)

Run the full MasP2 folder (10 groups × 3 replicates) twice —
`refine_on=False` vs `True`, all other config identical:

- **Primary metric:** within-group between-replicate std of the mean
  diameter. Expectation: consistent reduction across groups.
- **Secondary:** std of the detrended per-column diameter profile
  (along-fibre noise), refined coverage %, and a sanity check that σ(x) maps
  match visually defocused stretches.
- Record both runs' numbers in the labbook; on completion write
  `docs/report/01_esf_edge_consistency.md` (intro/method/results/discussion).

Acceptance for the feature: primary metric falls in a clear majority of
groups and rises in none by more than its old between-replicate noise;
coverage ≥ ~80% of valid columns on typical images.

## Testing

- Unit (synthetic): render Gaussian-blurred step profiles across σ ∈ [1, 15]
  with noise; refined edge must recover the true step within 0.3 px and be
  independent of σ (the current threshold edge shows a σ-proportional drift —
  regression-test the contrast).
- Unit (gates): specular-bump-contaminated and shadow-ramp profiles must be
  rejected by the residual gate and fall back to the anchor edge.
- Integration: `refine_on=False` reproduces today's outputs bit-identically
  on a golden image; `refine_on=True` leaves flags/validity of unrefined
  columns unchanged.

## Documentation obligations (on landing)

`docs/features/03_esf_edge_refinement.md` (template per CLAUDE.md), one-line
README mention, `refine.py` header comments (Input/Output/Pos), folder README
updates (`src/fibrecv/`, `docs/features/`), labbook + metalabbook updates.
