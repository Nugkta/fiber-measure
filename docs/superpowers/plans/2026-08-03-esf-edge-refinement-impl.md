# Implementation Plan — ESF Edge Refinement (study 01)

Status: approved 2026-08-03 — executing Tasks 1→5 (milestones M1→M5).

## Context

Optical blur (PSF, σ ≈ 5–15 px, focus-dependent) smears the air/silk boundary;
the current recipe places edges at a fixed low threshold near the wall-ramp
base, giving a focus-dependent diameter bias of several px per side (measured
+2.3 / −3.4 / +3.4 px on three MasP2 images). The erf midpoint of the blurred
profile is invariant to σ (knife-edge principle). New stage `refine.py` refits
each detected wall as a Gaussian-blurred step on the D z-map and moves the
boundary to the fitted 50% point; falls back to the current edge when the fit
fails. Goal: replicate-to-replicate consistency.

- Spec (authoritative): `docs/superpowers/specs/2026-08-03-esf-edge-refinement-design.md`
- Study/labbook: `docs/labbooks/01_esf_edge_consistency.md`

## Global Constraints

1. **Perpendicular frame throughout**: profile coordinate `t` is perpendicular
   distance; sample at `y = anchor ± t/cth` (`m, cth = band.tilt_geometry(bnd.slope)`
   — the single tilt source, band.py:60); `inside = clip(refine_in_frac ·
   band_half · cth, 8, refine_in_max)`; fitted t0/σ are perpendicular px so
   gates apply unconverted; applied vertical offset = `t0/cth`; identities at 0 tilt.
2. **No pre-smoothing of D** (would inflate fitted σ; block mean + LSQ over
   ~86–127 samples is the noise handling; study fitted raw channels at 3.5–6.7% resid).
3. **RefineResult adds `o_top`/`o_bot`** (applied per-column offsets) beyond the
   spec dataclass — needed for meta's median |t0|; aids tests.
4. Offsets applied **only to anchor columns** (finite y_top & y_bot AND
   flags==FLAG_OK) within interpolation coverage; flagged columns never move.
5. **Bit-identical when off** = all numeric arrays (same `EdgeResult` object
   returned); meta necessarily gains `refine.enabled=false` + new params keys.
6. σ/resid interpolated between passing block centres like the offsets.
7. Blocks: consecutive full `refine_block` runs from `bnd.x0`; trailing partial
   skipped; integer block centres; no extrapolation beyond outermost passing centres.
8. Columns whose window leaves the image are dropped **before** the 70% quorum
   check (no replicated-border contamination).
9. `refine_on=True` at landing (spec default); `--no-refine` is the A/B control.
10. Manual edits run after refine and override it (no re-refine) — Caveats item.
11. float32 throughout (NEP-50) so diameter matches edges.py:359 semantics.
12. Perf budget: < 1 s added per 2560-px image (~320 fits/image).
13. `uv run pytest` green after every task; the 69 existing tests must not
    change (except one added GUI smoke assertion in Task 3).

## Task 1 — M1: config + identity wiring

Files: `src/fibrecv/config.py`, `src/fibrecv/refine.py` (new skeleton),
`src/fibrecv/compute.py`, `tests/test_refine.py` (new, group A).

### config.py
New `# --- erf edge refinement ---` section after the window-sizing section
(~line 95), house comment style, with defaults:
```
refine_on=True, refine_block=16, refine_out=35.0, refine_in_frac=0.8,
refine_in_max=28.0, refine_relmax=0.08, refine_sigma_min=0.8,
refine_sigma_max=20.0, refine_maxshift=12.0, refine_gap_blocks=2
```

### refine.py skeleton
House-style header docstring (Dependencies/Inputs/Output/"Erf model & gates"/
Pos — "Fifth stage" — study existing modules e.g. edges.py/qc.py for the
exact house format). Dataclass:
```python
@dataclass
class RefineResult:
    refined_top: np.ndarray  # bool (W,)
    refined_bot: np.ndarray
    sigma_top: np.ndarray    # float32 (W,) perpendicular px, NaN unrefined
    sigma_bot: np.ndarray
    resid_top: np.ndarray    # float32 (W,) relative residual
    resid_bot: np.ndarray
    o_top: np.ndarray        # float32 (W,) applied offset (+into fibre), NaN unrefined
    o_bot: np.ndarray
    n_blocks: int            # attempted (quorum met; same both sides)
    n_pass_top: int
    n_pass_bot: int
    @classmethod
    def empty(cls, W): ...
```
`refine_edges(D, edg, bnd, cfg)`: when `cfg.refine_on` is False, return
`(edg, RefineResult.empty(W))` with `edg` the SAME object (identity, not a
copy). In this task the on-path is a placeholder: also return
`(edg, RefineResult.empty(W))` (no-op refinement, n_blocks=0) — the real
algorithm lands in Task 2.

### compute.py
- Insert the refine call between current lines 74/75 (i.e. after edges, before
  qc).
- `MeasureResult.ref: RefineResult | None = None` — the field goes **after
  `meta`** (a default is required — tests/test_manual_edit.py:80 constructs
  MeasureResult directly without it).
- meta gains a `"refine"` sub-dict placed between `"flag_counts"` and
  `"median_diameter_um"` (precedent for nested meta dict: manual_edit.py:239-255)
  with keys: `enabled`, `n_blocks`, `n_pass_top`, `n_pass_bot`,
  `coverage_top`, `coverage_bot` (refined columns / anchor columns, 0.0 if no
  anchors), `median_sigma_top`, `median_sigma_bot` (None-safe when nothing
  refined), `median_abs_t0_top`, `median_abs_t0_bot` (None-safe). All values
  JSON-safe (plain Python types).

### tests/test_refine.py — group A (identity/wiring)
Shared generator `_blurred_fibre(sigma_px, angle_deg=0, width=60)` adapted
from `_inclined_fibre` (tests/test_edges_tilt.py:19-37) but with a true erf
edge blend and fg=(0.95, 0.90, 0.91) so saturation is exactly affine in the
blend. Tests:
- off-path (`refine_on=False`) bit-identical: `np.testing.assert_array_equal`
  on all numeric arrays of the edge result AND `is` identity on the EdgeResult
  object returned by refine_edges.
- MeasureResult can still be constructed without `ref` (default None).
- meta `"refine"` sub-dict exists, is placed between `"flag_counts"` and
  `"median_diameter_um"` in key order, and `json.dumps(meta)` succeeds.

### Gate
`uv run pytest` fully green (69 existing + new group A).

## Task 2 — M2: the algorithm

Files: `src/fibrecv/refine.py` (full implementation), `tests/test_refine.py`
(groups B–D). No other production files change.

### refine.py helpers
- `_erf_model(t, a, b, t0, sigma)` — Gaussian-blurred step:
  `a + (b - a) * 0.5 * (1 + erf((t - t0) / (sigma * sqrt(2))))`.
- `_block_profile(D, y_anchor, cols, t_grid, sign, cth)` — vectorised bilinear
  gather (pattern: `edges._vshift`, edges.py:121-131) sampling
  `y = anchor + sign * t/cth` per column, then mean over the block's columns;
  returns None below the 70% column quorum. Columns whose sample window leaves
  the image are dropped BEFORE the quorum check (no replicated-border
  contamination).
- `_fit_block(t, prof, cfg)` — `scipy.optimize.curve_fit`, p0 = [mean of
  outer-6 samples, mean of inner-6 samples, t at the sample closest to the
  level midpoint, 3.0], loose trf bounds (σ > 0.05, t0 inside the window),
  maxfev=2000, try/except → fail. Gates (all must pass): b − a > 0,
  rms_residual/(b − a) < refine_relmax, σ ∈ [refine_sigma_min,
  refine_sigma_max], |t0| ≤ refine_maxshift.
- `_interp_side(...)` — offsets (and σ, resid) interpolated between passing
  block centres via `np.interp` per chain; chains split where the gap between
  consecutive passing centres exceeds refine_gap_blocks blocks; NaN outside
  chains; no extrapolation beyond the outermost passing centres.

### refine_edges on-path
Anchor columns = finite y_top & y_bot AND flags==FLAG_OK. Blocks: consecutive
full `refine_block`-column runs starting at `bnd.x0`; trailing partial block
skipped; integer block centres. Profile window: outward `refine_out` px,
inward `inside = clip(refine_in_frac * band_half * cth, 8, refine_in_max)`
(perpendicular px). Tilt frame: `m, cth = band.tilt_geometry(bnd.slope)`;
sample at `y = anchor ± t/cth`; fitted t0/σ are perpendicular px (gates apply
unconverted); applied vertical offset = `t0/cth`. Apply on refined columns
only: `y_top += o_top/cth`, `y_bot -= o_bot/cth`, then **recompute
`diameter = (y_bot − y_top) * cth`** (qc.py:93 trusts `edg.diameter` and never
recomputes — forgetting this is the one hard bug). Return
`dataclasses.replace(edg, y_top=…, y_bot=…, diameter=…)`;
`flags`/`amp`/`y_core`/`half_window` untouched. float32 throughout. No
pre-smoothing of D. ~320 fits/image, < 1 s added.

### tests groups B–D
- **B fit/gates (1-D, on synthetic profiles)**: midpoint recovery ≤ 0.3 px
  across σ ∈ {1,2,4,8,12,15} × t0 ∈ {−4,0,+5} with noise amplitude
  (b−a)/40; **contrast test** — a legacy fixed-threshold crossing drifts by
  px-scale amounts as σ goes 1→15 while the erf midpoint stays ≤ 0.3 px;
  rejection tests: specular-bump profile, shadow-ramp profile, inverted step
  (b−a ≤ 0), σ out of range, |t0| > maxshift.
- **C offset field**: interpolation across ≤ 2-block holes; > 2-block holes
  and chain ends stay NaN; block below 70% quorum skipped; flagged/unrefined
  columns untouched (array-equal on those columns); all-fits-fail returns
  input y arrays unchanged; border-clipped columns excluded before quorum.
- **D integration (full pipeline on `_blurred_fibre` images)**: refined median
  diameter independent of σ ∈ {2,6,12} (≤ 1 px spread; legacy spread > 2 px);
  tilt invariance at {0°, 20°, 35°} within 2% (mirrors
  test_edges_tilt.py:96-112); σ-map recovers the rendered blur ± 20%.

### Gate
`uv run pytest` fully green + eyeball one real MasP2 image (read-only input,
output image/plot to the session scratchpad only, then delete).

## Task 3 — M3: CLI + GUI

Files: `src/fibrecv/run_measure.py`, `src/fibrecv/gui_app.py`,
existing GUI smoke test file (one added assertion).

### run_measure.py
`--refine/--no-refine` via `argparse.BooleanOptionalAction`, `default=None`
(must survive the overrides None-filter at run_measure.py:43-65); add the
entry to the overrides dict mapping onto `refine_on`.

### gui_app.py
- `PARAM_SPECS` gains a bool entry for `refine_on` placed **before** the
  `ppu` entry (the Calibration divider triggers on name=='ppu', line 679).
- New `kind == "bool"` checkbox branch in the `_param_form` loop (lines
  678-696), with a `.get` fallback for pre-upgrade session state.
- `_profile_fig` (lines 539-554) gains a twinx σ(x) trace (style consistent
  with `_styled_fig` + `_MUTED`/`_ACCENT`), gated on
  `getattr(mr, "ref", None)` being present with finite σ values.
- One added smoke-test assertion: the refine checkbox renders before the
  Calibration divider.

### Gate
Suite green; `uv run python -m fibrecv.run_measure --help` shows the flag;
Streamlit AppTest boots; meta `refine.enabled` matches flag polarity on a
real image.

## Task 4 — M4: docs

No production code. Changes:
- qc.py header Pos: "Fifth stage" → "Sixth stage".
- The pipeline arrow chain `features -> band -> edges -> qc` gains `refine`
  everywhere it appears (src/fibrecv/README.md:5, compute.py and measure.py
  docstrings — grep for `edges -> qc` to find all).
- src/fibrecv/README.md: file bullet for refine.py.
- `docs/features/03_esf_edge_refinement.md` — CLAUDE.md feature template
  (What it does / Design choices / Algorithm details / Caveats), with the
  Global Constraints decision log restated under Design choices; Caveat:
  manual edits run after refine and override it (no re-refine).
- docs/features/README.md bullet; root README.md one-liner.
- GUI_README.md: "Four knobs" → "Five" + σ-trace mention.
- Labbook `docs/labbooks/01_esf_edge_consistency.md`: extend today's entry
  (same-date extension, newest-first rule); `docs/metalabbook.md` Last
  Updated refresh + re-sort; docs/labbooks/README.md if file list changed.

### Gate
`grep -rn "edges -> qc\|Fifth stage\|Four knobs"` returns only intentional
hits; suite still green.

## Task 5 — M5: validation (study 01 main experiment)

No production code changes. Steps:
- Run A/B on the full MasP2 set:
  `uv run python -m fibrecv.run_measure --root "/Users/stan/Documents/UOM/spins/Images MasP2" --all --jobs 6 --no-refine --out fibrecv_output/ab_refine_off`
  and the `--refine` twin to `fibrecv_output/ab_refine_on`; run_aggregate on
  each.
- Ephemeral analysis script in the session scratchpad (deleted after):
  primary metric = within-group between-replicate std of per-image mean
  diameter, off vs on; secondary = detrended along-fibre noise, refine
  coverage %, median σ, σ(x) vs visible defocus. Acceptance: primary falls in
  a clear majority of groups, rises nowhere beyond old noise; coverage ≥ ~80%.
  Tuning knobs if under-covered: `refine_relmax` first (16-col blocks are
  noisier than the study's 64-col), then `refine_block`.
- Labbook entry with the numbers; report
  `docs/report/01_esf_edge_consistency.md` (intro, method, results,
  discussion); study status → done in both labbook header and metalabbook row.

## Risks to watch

- Thin fibres: inside-window floor (8 px) can see the opposite wall → gates
  reject → fallback (spec non-goal; watch failure rate in Task 5).
- Interpolation seams can trip qc rolling-MAD on a few columns — check
  `flag_counts` deltas in Task 5.
- σ≈15 fits: p0 tails not saturated at window ends — maxfev=2000 headroom.
- GUI stale sessions/caches: `.get` fallback + `getattr(mr, "ref", None)`.
- Keep float32 (NEP-50) so diameter matches edges.py:359 semantics.
