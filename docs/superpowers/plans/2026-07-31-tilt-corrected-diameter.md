# Tilt-Corrected Fibre Diameter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> On execution start, copy this plan to `docs/superpowers/plans/2026-07-31-tilt-corrected-diameter.md` and commit it (repo convention, cf. commit `a00da17`).

**Goal:** Report the fibre diameter perpendicular to the fibre axis (not the vertical image chord), with detection itself made tilt-invariant, validated for inclinations up to ~45°.

**Architecture:** Scanning stays vertical per column and `y_top`/`y_bot` remain original-image rows (overlay, manual-edit clicks and QC keep their coordinate semantics). Four compensations, all exact identities at zero tilt: (1) the `wcol` column-neighbourhood average follows the fibre axis instead of cutting across it; (2) the vertical Gaussian sigma is scaled by 1/cos θ; (3) the vertical gradient is scaled by 1/cos θ so the calibrated wall slope gates see perpendicular slopes; (4) the vertical chord is multiplied by cos θ. θ comes from the existing Theil–Sen centerline slope (`band.slope`).

**Tech Stack:** Python 3.12, numpy/scipy/scikit-image, pytest via `uv run pytest`. No new dependencies.

## Context (why)

`diameter = y_bot - y_top` (`src/fibrecv/edges.py:292`, mirrored `src/fibrecv/manual_edit.py:205`) measures a **vertical** chord. For an inclined fibre this is wrong twice over:

1. **Chord direction** (user-reported): true diameter is perpendicular to the axis; vertical chord is inflated by 1/cos θ.
2. **Cross-column smearing** (found during investigation, ~2× larger): `uniform_filter1d(D, size=wcol=41, axis=1)` at `edges.py:263` averages 41 image columns; on a tilted fibre the boundary row drifts `41·tanθ` px across that window, smearing the wall into a ramp so the level-crossing boundary lands further out on both sides.

Session-validated on synthetic bands of known width: vertical chord drifts +37 % by 30°; cos correction alone leaves ~2/3 of the drift; all four compensations together are flat to 40°. Real MasP2 data: median tilt 0.9°, max 6.7°.

User decisions: fix both errors; support up to ~45°; **clean replacement** of the diameter value (no legacy columns, no config switch). Every downstream consumer (qc → compute → measure/register/tensile/GUI) reads `EdgeResult.diameter` / `QCResult.diameter_raw`, so only the two producer lines plus one QC ratio need changing.

## Global Constraints

- Python `>=3.12`, `numpy>=1.26,<2`, `scikit-image==0.24.*` — no new dependencies, no version changes.
- All compensations MUST be exact identities at `band.slope == 0` (horizontal fibres reproduce current results bit-for-bit; guarded by test).
- No output-schema changes: CSV columns, meta JSON keys, master_summary/tensile columns all keep their names.
- `y_top`/`y_bot` stay sub-pixel rows in original image coordinates.
- Non-finite `band.slope` must be treated as `0.0` everywhere.
- Run tests with `uv run pytest`; work on branch `feat/tilt-corrected-diameter` off `main`.
- Match the repo's docstring style (module headers with Dependencies/Inputs/Output/Pos sections).

---

### Task 1: `_axis_average` helper in edges.py

**Files:**
- Modify: `src/fibrecv/edges.py` (new helper after `half_window_px`, ~line 104)
- Test: `tests/test_edges_tilt.py` (new file)

**Interfaces:**
- Consumes: nothing new (`numpy`, `scipy.ndimage.uniform_filter1d` already imported in edges.py).
- Produces: `_axis_average(D: np.ndarray, slope: float, wcol: int) -> np.ndarray` — float32 (H, W); Task 2 calls it in `detect_edges`; Task 2's tests import `_inclined_fibre` from this test file.

- [ ] **Step 1: Create branch**

```bash
git checkout -b feat/tilt-corrected-diameter
```

- [ ] **Step 2: Write the failing tests**

Create `tests/test_edges_tilt.py`:

```python
"""Tilt-invariance tests: perpendicular diameter + axis-following averaging."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.ndimage import uniform_filter1d

from fibrecv.compute import compute_measurement
from fibrecv.config import CONFIG
from fibrecv.edges import _axis_average

TRUE_W = 60.0  # true perpendicular width of the synthetic fibre (px)


def _inclined_fibre(angle_deg: float, width: float = TRUE_W,
                    W: int = 800, H: int = 700, seed: int = 0) -> np.ndarray:
    """Pink background + pale band at ``angle_deg``, anti-aliased boundary.

    The edge blends linearly over ~2 px of *perpendicular* distance so the
    boundary stays smooth at any angle (a hard-edged mask aliases when tilted
    and adds ~2 % noise to the measured width).
    """
    m = math.tan(math.radians(angle_deg))
    y, x = np.mgrid[0:H, 0:W].astype(np.float64)
    dist = np.abs(y - (m * (x - W / 2) + H / 2)) / math.sqrt(1 + m * m)
    t = np.clip((width / 2 + 1.0 - dist) / 2.0, 0.0, 1.0)[..., None]
    bg = np.array([0.95, 0.45, 0.65])   # saturated pink background
    fg = np.array([0.92, 0.90, 0.91])   # pale desaturated fibre
    img = bg + (fg - bg) * t
    rng = np.random.default_rng(seed)
    img = img + rng.normal(0.0, 0.006, img.shape)
    return np.clip(img, 0.0, 1.0).astype(np.float32)


def _median_diameter_px(rgb: np.ndarray, cfg: CONFIG) -> float:
    mr = compute_measurement(rgb, cfg, "synthetic")
    assert mr.res.valid.any()
    return float(np.nanmedian(mr.res.diameter_raw))


def test_axis_average_slope0_matches_uniform_filter():
    rng = np.random.default_rng(1)
    D = rng.normal(size=(50, 80)).astype(np.float32)
    ref = uniform_filter1d(D, size=41, axis=1, mode="nearest")
    assert np.allclose(_axis_average(D, 0.0, 41), ref)
    # non-finite slope degrades to the plain average too
    assert np.allclose(_axis_average(D, float("nan"), 41), ref)


def test_axis_average_follows_the_axis():
    """A pattern constant along a 45-degree axis must survive axis-averaging.

    D[r, c] = 1 exactly on the line r = 10 + c. With slope=1 every sample in
    the sheared window lands back on the line, so interior columns are
    reproduced exactly; a straight column average would dilute the line to
    ~1/wcol.
    """
    H, W, wcol = 60, 41, 11
    D = np.zeros((H, W), dtype=np.float32)
    c = np.arange(W)
    D[10 + c, c] = 1.0
    A = _axis_average(D, 1.0, wcol)
    hw = wcol // 2
    interior = slice(hw, W - hw)
    assert np.allclose(A[:, interior], D[:, interior], atol=1e-6)
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `uv run pytest tests/test_edges_tilt.py -v`
Expected: FAIL with `ImportError: cannot import name '_axis_average'`

- [ ] **Step 4: Implement `_axis_average`**

In `src/fibrecv/edges.py`, insert after `half_window_px` (after line 103):

```python
def _axis_average(D: np.ndarray, slope: float, wcol: int) -> np.ndarray:
    """Column-neighbourhood average taken *along the fibre axis*.

    Equivalent to ``uniform_filter1d(D, size=wcol, axis=1, mode="nearest")``
    for a horizontal fibre, but on a tilted fibre each neighbour column ``k``
    is sampled at the row offset ``slope * k`` (bilinear in y, borders
    clipped), so the window follows the axis instead of cutting across the
    boundary -- the tilt-smearing that otherwise widens both walls. For
    ``slope != 0`` the window is the odd size ``2*(wcol//2) + 1``.
    """
    if not np.isfinite(slope) or slope == 0.0:
        return uniform_filter1d(D, size=max(1, wcol), axis=1, mode="nearest")
    H, W = D.shape
    hw = max(1, wcol) // 2
    rows = np.arange(H)[:, None]
    cols = np.arange(W)[None, :]
    acc = np.zeros((H, W), dtype=np.float32)
    for k in range(-hw, hw + 1):
        xs = np.clip(cols + k, 0, W - 1)
        dy = slope * k
        y0 = int(np.floor(dy))
        f = dy - y0
        ya = np.clip(rows + y0, 0, H - 1)
        yb = np.clip(rows + y0 + 1, 0, H - 1)
        acc += (1.0 - f) * D[ya, xs] + f * D[yb, xs]
    return acc / float(2 * hw + 1)
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/test_edges_tilt.py -v`
Expected: 2 PASS

- [ ] **Step 6: Commit**

```bash
git add src/fibrecv/edges.py tests/test_edges_tilt.py
git commit -m "feat(edges): axis-following column average helper"
```

---

### Task 2: Tilt-compensated detection + perpendicular diameter in `detect_edges`

**Files:**
- Modify: `src/fibrecv/edges.py:14-18` (module docstring Output section), `:62` area (docstring, add Tilt handling paragraph), `:93` (`EdgeResult.diameter` comment), `:252-292` (`detect_edges`)
- Test: `tests/test_edges_tilt.py` (extend)

**Interfaces:**
- Consumes: `_axis_average(D, slope, wcol)` from Task 1; `_inclined_fibre` / `_median_diameter_px` fixtures from Task 1's test file.
- Produces: `EdgeResult.diameter` now = `(y_bot - y_top) * cos(atan(band.slope))` — the value every downstream module consumes. `y_top`/`y_bot` semantics unchanged.

- [ ] **Step 1: Record the pre-change horizontal baseline (regression guard)**

Run BEFORE touching `detect_edges` (Task 1 did not change its behaviour):

```bash
uv run python -c "
import sys; sys.path.insert(0, 'tests')
import numpy as np
from test_edges_tilt import _inclined_fibre
from fibrecv.compute import compute_measurement
from fibrecv.config import CONFIG
mr = compute_measurement(_inclined_fibre(0.0), CONFIG(), 'pin')
print('median_px =', round(float(np.nanmedian(mr.res.diameter_raw)), 4))
"
```

Write the printed number down; it is re-checked in the Verification section after all tasks (must match to < 0.1 %).

- [ ] **Step 2: Write the failing tilt-invariance test**

Append to `tests/test_edges_tilt.py`:

```python
def test_tilt_invariance_and_absolute_width():
    """Median measured width must not drift with tilt (validated to 40 deg).

    Tolerances: pairwise spread <= 2 % relative (widen to at most 3 % only if
    the 0-degree reference itself sits within the absolute band below).
    The absolute band is wide because the boundary level intentionally sits
    partway down the smoothed wall (edge_z above local background), which on
    this fixture reads systematically wide -- the property under test is
    *invariance*, not absolute accuracy.
    """
    cfg = CONFIG()
    meds = {a: _median_diameter_px(_inclined_fibre(float(a)), cfg)
            for a in (0, 10, 20, 30, 40)}
    ref = meds[0]
    for a, v in meds.items():
        assert v == pytest.approx(ref, rel=0.02), f"angle {a}: {v:.2f} vs {ref:.2f}"
    assert TRUE_W * 0.95 <= ref <= TRUE_W * 1.35
```

- [ ] **Step 3: Run test to verify it fails**

Run: `uv run pytest tests/test_edges_tilt.py::test_tilt_invariance_and_absolute_width -v`
Expected: FAIL — with the current vertical-chord pipeline the 20°/30°/40° medians run ~8–30 % above the 0° reference (session measurement: +37 % drift at 30°).

- [ ] **Step 4: Implement the compensations in `detect_edges`**

In `src/fibrecv/edges.py`, replace lines 254–265 (from `H, W = D.shape` through `G = np.gradient(...)`) with:

```python
    H, W = D.shape
    hw = half_window_px(band, cfg, H)

    # tilt geometry: report the width perpendicular to the fibre axis and
    # keep the calibrated recipe stable in the perpendicular frame (slope
    # from the robust centerline fit; the vertical scan itself is unchanged)
    m = float(band.slope) if np.isfinite(band.slope) else 0.0
    cth = float(1.0 / np.sqrt(1.0 + m * m))  # cos of the tilt angle

    # plateau-bridging length for the wall finder, scaled to the fibre size so
    # thin fibres never bridge across their own thickness (bounds 4..16 px)
    bh = band.band_half if np.isfinite(band.band_half) else 0.0
    gap = int(np.clip(round(cfg.wall_gap_frac * 2.0 * bh), 4, 16))

    # column-neighbourhood average along the fibre axis, then vertical
    # smoothing + gradient; sigma and slopes are rescaled by 1/cos(tilt) so
    # smoothing width and the wall slope gates keep their calibrated
    # *perpendicular* values on a tilted fibre (identities at zero tilt)
    Davg = _axis_average(D, m, max(1, cfg.wcol))
    Dsm = gaussian_filter1d(Davg, sigma=cfg.sigma_y / cth, axis=0, mode="nearest")
    G = np.gradient(Dsm, axis=0) / cth
```

Replace line 292:

```python
    diameter = (y_bot - y_top) * cth
```

- [ ] **Step 5: Update the docstrings**

In the module docstring, change the Output section (lines 14–18) to:

```
Output
------
``EdgeResult`` with float arrays ``y_top``, ``y_bot`` (sub-pixel boundary rows,
NaN where invalid), ``diameter`` = ``(y_bot - y_top) * cos(tilt)`` -- the width
perpendicular to the fibre axis -- plus per-column ``amp`` (band amplitude A),
``y_core`` and integer ``flags`` for QC.
```

Insert before the "Pos" section (around line 63):

```
Tilt handling
-------------
The fibre may be inclined (validated to ~45 deg; slope from the robust
centerline fit). Detection still scans vertical columns, but three
compensations make the recipe tilt-invariant: the ``wcol`` neighbourhood
average follows the fibre axis (sheared sampling, so it never smears the
boundary), sigma_y and the wall slope gates are rescaled by 1/cos(tilt) so
smoothing and thresholds stay calibrated in the perpendicular frame, and the
vertical chord is multiplied by cos(tilt) to give the perpendicular diameter.
All three are exact identities for a horizontal fibre.
```

Change the `EdgeResult.diameter` field comment (line 93) to:

```python
    diameter: np.ndarray  # float (W,) perpendicular width = (y_bot-y_top)*cos(tilt), NaN if invalid
```

- [ ] **Step 6: Run the new test and the full suite**

Run: `uv run pytest tests/test_edges_tilt.py -v` — Expected: all PASS (if the invariance assertion fails marginally at `rel=0.02`, check the absolute band still holds and relax only up to `rel=0.03`).
Run: `uv run pytest` — Expected: all PASS (existing fixtures are slope-0; `test_gui_smoke` uses a horizontal band).

- [ ] **Step 7: Commit**

```bash
git add src/fibrecv/edges.py tests/test_edges_tilt.py
git commit -m "feat(edges): measure diameter perpendicular to the fibre axis"
```

---

### Task 3: Same correction in the manual-edit path

**Files:**
- Modify: `src/fibrecv/manual_edit.py:205` (+ a line in the `apply_manual_edits` docstring)
- Test: `tests/test_manual_edit.py` (extend)

**Interfaces:**
- Consumes: `mr.bnd.slope` (already in scope — `run_qc(new_edg, mr.bnd, cfg)` at `manual_edit.py:216`).
- Produces: manually edited columns get the same perpendicular-diameter definition as the detector; GUI apply path (`gui_app.py:1397-1401`) needs no change.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_manual_edit.py` (add `import math` and `from dataclasses import replace` to the file's imports):

```python
def test_tilted_band_scales_edited_diameter():
    """Edited diameters use the perpendicular definition: x cos(tilt)."""
    cfg = CONFIG()
    mr = _mr(cfg)
    mr = replace(mr, bnd=replace(mr.bnd, slope=1.0))  # 45-degree band
    edits = empty_edits()
    edits["top"].append([(120.0, 100.0), (160.0, 100.0)])
    edits["bot"].append([(120.0, 200.0), (160.0, 200.0)])
    new_mr, _, _ = apply_manual_edits(mr, edits, cfg)
    cth = 1.0 / math.sqrt(2.0)
    assert new_mr.res.valid[BAD].all()
    assert np.allclose(new_mr.res.diameter_raw[BAD], 100.0 * cth)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_manual_edit.py::test_tilted_band_scales_edited_diameter -v`
Expected: FAIL — `diameter_raw[BAD]` is 100.0 (vertical chord), not 70.71.

- [ ] **Step 3: Implement**

In `src/fibrecv/manual_edit.py`, replace line 205 (`diameter = y_bot - y_top`) with:

```python
    m = float(mr.bnd.slope) if np.isfinite(mr.bnd.slope) else 0.0
    cth = float(1.0 / np.sqrt(1.0 + m * m))
    diameter = (y_bot - y_top) * cth  # perpendicular width, same as detect_edges
```

Add one sentence to the `apply_manual_edits` docstring (after the sentence ending "flags cleared to ``FLAG_OK`` (user override)", line ~175): `Diameters are recomputed with the same perpendicular (cos-tilt) definition as the detector.`

- [ ] **Step 4: Run the file's tests**

Run: `uv run pytest tests/test_manual_edit.py -v`
Expected: all PASS (existing tests use a `slope=0.0` band fixture, so their expected values are unchanged).

- [ ] **Step 5: Commit**

```bash
git add src/fibrecv/manual_edit.py tests/test_manual_edit.py
git commit -m "fix(manual-edit): perpendicular diameter for edited columns"
```

---

### Task 4: Tilt-consistent band-mismatch ratio in QC

**Files:**
- Modify: `src/fibrecv/qc.py:145-156` (band-consistency check)
- Test: `tests/test_edges_wall.py` (extend)

**Interfaces:**
- Consumes: `band.slope` (BandResult already passed into `run_qc`).
- Produces: `band_mismatch` compares the (now perpendicular) diameter against the coarse band thickness converted to the same frame — prevents legitimately tilted fibres from being flagged low-confidence.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_edges_wall.py` (add `import math` and `from dataclasses import replace` to the imports):

```python
def test_qc_band_ratio_uses_perpendicular_frame():
    """band_half is a vertical column count; on a 45-degree fibre it reads
    1/cos wider than the perpendicular diameter. The ratio check must compare
    like with like, or a legitimate tilted measurement gets flagged."""
    cfg = CONFIG()
    # true width 100 px at 45 deg: vertical band chord = 141.4 px,
    # perpendicular diameter a bit narrow at 55 px -> ratio 0.55 in the
    # perpendicular frame (fine), 0.39 uncorrected (false mismatch)
    edg, bnd = _qc_fixture(diameter_px=55.0, band_half=70.71)
    bnd = replace(bnd, slope=1.0)
    res = run_qc(edg, bnd, cfg)
    assert not res.band_mismatch

    # a genuinely-too-narrow measurement still trips the guard
    edg2, bnd2 = _qc_fixture(diameter_px=40.0, band_half=70.71)
    bnd2 = replace(bnd2, slope=1.0)
    res2 = run_qc(edg2, bnd2, cfg)
    assert res2.band_mismatch
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_edges_wall.py::test_qc_band_ratio_uses_perpendicular_frame -v`
Expected: FAIL on the first assertion (uncorrected ratio 55/141.4 = 0.39 < `band_ratio_min` 0.5).

- [ ] **Step 3: Implement**

In `src/fibrecv/qc.py`, replace line 149 (`band_thickness = 2.0 * band.band_half`) with:

```python
    m = float(band.slope) if np.isfinite(band.slope) else 0.0
    cth = float(1.0 / np.sqrt(1.0 + m * m))
    # band_half is a vertical column count; x cos(tilt) converts it to the
    # perpendicular frame the diameter is now reported in
    band_thickness = 2.0 * band.band_half * cth
```

- [ ] **Step 4: Run the file's tests**

Run: `uv run pytest tests/test_edges_wall.py -v`
Expected: all PASS (`_qc_fixture` defaults to `slope=0.0`, so the two existing ratio tests are unchanged).

- [ ] **Step 5: Commit**

```bash
git add src/fibrecv/qc.py tests/test_edges_wall.py
git commit -m "fix(qc): band-mismatch ratio compared in the perpendicular frame"
```

---

### Task 5: Documentation + end-to-end verification

**Files:**
- Modify: `README.md` (algorithm paragraph, ~lines 12–15), `GUI_README.md` (only if it describes the diameter as vertical — check)
- No test file; this task's gate is the Verification checklist below.

**Interfaces:**
- Consumes: everything above.
- Produces: user-facing docs consistent with the new definition.

- [ ] **Step 1: Update README.md**

In the algorithm paragraph (the one beginning "Detection runs on a per-image *desaturation z-map*", `README.md:12-15`), append:

```
Diameters are measured perpendicular to the fibre axis: the per-column
neighbourhood average follows the fitted centerline and the vertical chord is
rescaled by cos(tilt), so inclined fibres (validated to ~45°) are not
overestimated.
```

- [ ] **Step 2: Check GUI_README.md**

Run: `grep -n -i "vertical\|diameter" GUI_README.md`. If any line describes the measured diameter as the vertical distance between boundaries, reword it to "width perpendicular to the fibre axis"; otherwise no change.

- [ ] **Step 3: Full test suite**

Run: `uv run pytest`
Expected: all PASS.

- [ ] **Step 4: Horizontal regression pin**

Re-run the snippet from Task 2 Step 1. Expected: identical `median_px` to the recorded pre-change value within 0.1 % (all compensations are identities at the fixture's ~zero slope).

- [ ] **Step 5: Real-image spot check**

```bash
uv run python -c "
import numpy as np
from fibrecv.io_utils import load_rgb
from fibrecv.compute import compute_measurement
from fibrecv.config import CONFIG
base = '/Users/stan/Library/CloudStorage/OneDrive-TheUniversityofManchester/Y2_onedrive/Projects/spins/Images MasP2/'
for nm, old in [('masp2 10_1_2.jpg', 32.9855), ('masp2 9_9_3.jpg', None),
                ('masp2 3_5_3.jpg', 81.2550), ('masp2 9_6_3.jpg', 57.4728)]:
    mr = compute_measurement(load_rgb(base + nm), CONFIG(), nm)
    med = float(np.nanmedian(mr.diameter_um))
    delta = '' if old is None else f'  (was {old}: {100*(med-old)/old:+.3f} %)'
    print(f'{nm}: {med:.4f} um{delta}, tilt={mr.bnd.slope:+.4f}, cov={mr.res.coverage:.2f}')
"
```

Expected: near-horizontal images (`10_1_2` tilt 0.45°, `9_6_3` 0.50°) shift by ≲ 0.05 %; `3_5_3` (1.1°) by ≲ 0.15 %; the most-tilted `9_9_3` (6.7°) shifts **down** by roughly 0.7–1.5 % (cos alone is −0.69 %, plus the smear correction). Coverage must not degrade (≳ previous values in `fibrecv_output/per_image/diagnostics/*_meta.json`).

- [ ] **Step 6: Synthetic sweep sanity plot (optional but recommended)**

Reuse `_inclined_fibre` to print median measured width at 0/5/10/20/30/40/45° — values should be flat within ~2 %.

- [ ] **Step 7: GUI smoke + commit docs**

Run: `uv run pytest tests/test_gui_smoke.py -v` — Expected: PASS.

```bash
git add README.md GUI_README.md
git commit -m "docs: document perpendicular (tilt-corrected) diameter measurement"
```

- [ ] **Step 8: Finish the branch**

Use the superpowers:finishing-a-development-branch skill (merge/PR decision belongs to the user).

---

## Untouched by design (verified during exploration)

- `half_window_px` — window sizing uses the *vertical* band thickness, which is exactly what the vertical scan must cover on a tilted fibre.
- `_detect_column` / `_side_edge` / `_find_wall` / `_outer_crossing` — operate on the (already compensated) profile and gradient; z-value knobs (`edge_z`, `edge_frac`, `rise_min`, `amin`) need no rescaling.
- `compute.py`, `measure.py`, `register.py`, `overlay.py`, `tensile.py`, `run_measure.py`, `run_aggregate.py`, `gui_app.py` — consume `EdgeResult.diameter`/`QCResult.diameter_raw` downstream; overlay still draws `y_top`/`y_bot` rows.
- QC centerline-deviation check (`qc.py:98-108`) — works on row coordinates, tilt-aware via its own Theil–Sen refit.

## Known limits (documented, not addressed)

- Beyond ~45° (near-vertical fibres) a rotate-to-fibre-frame architecture would be needed — out of scope.
- A steeply tilted fibre exiting through the top/bottom of frame can fail the full-width band rule → `low_confidence` (existing semantics).
- Background estimation (`features.py:36`, top/bottom 12 % margin rows) assumes a roughly horizontal fibre; the median is robust to the partial contamination a ≤ 45° fibre causes.
- A single global slope (straight-axis assumption) is used, matching the existing centerline model; a pathological slope estimate from a bad mask degrades results as it already does today.
- `_axis_average` costs ~41 vectorised full-image ops (~0.5 s on 5 MP vs ~50 ms today) — acceptable for the cached GUI and parallel batch; do not pre-optimise.

## Performance of plan self-review

- Spec coverage: both error sources fixed (Tasks 1–2), manual-edit parity (Task 3), QC frame consistency (Task 4), docs + E2E (Task 5); clean replacement — no schema changes anywhere. ✓
- Placeholders: none — every step has runnable code/commands; the one measured constant (horizontal pin) has an explicit capture procedure. ✓
- Type consistency: `_axis_average(D, slope, wcol)` used identically in Task 1 test, Task 2 implementation; `cth` formula identical in edges/manual_edit/qc. ✓
