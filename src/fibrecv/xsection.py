"""Pure multi-angle cross-section math: alignment stack, ellipse + hexagon.

Dependencies
------------
``numpy`` plus ``register.estimate_shift`` / ``register.resample_to_grid``
and ``config.CONFIG``. Deliberately no file I/O (mirrors the compute/measure
split); ``run_xsection.py`` owns all reading and writing.

Inputs
------
- Per-angle width profiles of one (fiber, part): ``{angle: {"x", "w"}}`` with
  ``x`` in ABSOLUTE image pixels (span-restricted is fine) and ``w`` the
  NaN-masked per-column diameter in px.
- ``CONFIG`` for the cross-angle alignment knobs ``xsec_min_corr`` /
  ``xsec_max_shift`` (wired into ``estimate_shift`` via ``dataclasses.replace``
  — the aggregate-stage ``min_corr``/``max_shift`` stay untouched).

Output
------
- ``build_part_stack`` -> ``PartStack`` with the aligned ``W (6, N)`` width
  stack (row k = angle a<k+1>, NaN where missing).
- ``fit_ellipse_projections(W, theta_deg)`` -> ``XsecFit`` per-column ellipse
  (fit in squared-width space, exactly linear; a >= b by construction; phi is
  the rotation angle at which the measured width is maximal, in [0, 180)).
- ``hexagon_area`` / ``hexagon_area_expected`` -> circumscribed-hexagon upper
  bound and its value expected for the fitted ellipse (QC ratio pair).
- ``split_half_area``, ``pair_differences`` -> 180-degree-pair uncertainty.
- ``predict_anisotropy``, ``predict_phi_transfer`` -> held-out validation
  errors (directional vs isotropic baseline; phi-transfer two-direction fit).

Pos
---
Third-stage compute heart. Consumed only by ``run_xsection.py`` and the test
suite. All math in px; micron conversion happens once, in the runner, from the
per-image Zeiss XML scale.

Once I am updated, update my header comments and folder's md.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, replace

import numpy as np

from .config import CONFIG
from .register import estimate_shift, resample_to_grid

#: Nominal rotation angles of a1..a6 (deg). 60-degree spacing is the working
#: assumption (labbook 02); every consumer takes theta_deg as a parameter so a
#: corrected spacing only requires re-running, not re-coding.
NOMINAL_ANGLES_DEG = np.array([0.0, 60.0, 120.0, 180.0, 240.0, 300.0])

_SQRT3 = np.sqrt(3.0)


@dataclass
class PartStack:
    """Aligned width profiles of one (fiber, part) across the six angles."""

    fiber: int
    part: int
    x: np.ndarray          # (N,) common grid, absolute px of the ref angle
    W: np.ndarray          # (6, N) width px, NaN = missing
    shifts: list           # per angle: angle/present/shift_px/corr_peak/uncertain


@dataclass
class XsecFit:
    """Per-column ellipse fit of the width stack."""

    a: np.ndarray          # (N,) semi-major axis, px (NaN invalid)
    b: np.ndarray          # (N,) semi-minor axis, px
    phi_deg: np.ndarray    # (N,) major-axis angle [0,180); NaN near-circular
    area: np.ndarray       # (N,) pi*a*b, px^2
    resid: np.ndarray      # (n_angles_in, N) w_obs - w_model, NaN unused
    rms_resid: np.ndarray  # (N,)
    n_angles: np.ndarray   # (N,) finite widths per column
    valid: np.ndarray      # (N,) bool


def _direction_ids(theta_deg: np.ndarray) -> np.ndarray:
    """Group angles by projection direction (theta mod 180 deg)."""
    mod = np.round(np.asarray(theta_deg, dtype=float) % 180.0, 6)
    _, ids = np.unique(mod, return_inverse=True)
    return ids


def _fit_coeffs(W: np.ndarray, theta_deg: np.ndarray
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Solve w^2 = c0 + c1*cos2t + c2*sin2t per column (vectorised by mask).

    Returns ``(C, solvable, n_angles)`` where ``C`` is (3, N) (NaN when the
    column lacks all three directions) and ``solvable`` marks columns with
    >= 3 distinct directions among their finite widths.
    """
    W = np.asarray(W, dtype=float)
    n_ang, n_col = W.shape
    theta = np.radians(np.asarray(theta_deg, dtype=float))
    design = np.column_stack(
        [np.ones(n_ang), np.cos(2 * theta), np.sin(2 * theta)])
    dir_ids = _direction_ids(theta_deg)
    n_dirs = int(dir_ids.max()) + 1

    finite = np.isfinite(W)
    n_angles = finite.sum(axis=0)
    C = np.full((3, n_col), np.nan)
    solvable = np.zeros(n_col, dtype=bool)

    code = (finite * (1 << np.arange(n_ang))[:, None]).sum(axis=0)
    for pattern in np.unique(code):
        rows = [k for k in range(n_ang) if pattern >> k & 1]
        cols = np.where(code == pattern)[0]
        if len({int(dir_ids[k]) for k in rows}) < n_dirs or n_dirs < 3:
            continue
        X = design[rows]
        Y = W[np.ix_(rows, cols)] ** 2
        coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
        C[:, cols] = coef
        solvable[cols] = True
    return C, solvable, n_angles


def fit_ellipse_projections(W: np.ndarray,
                            theta_deg: np.ndarray = NOMINAL_ANGLES_DEG
                            ) -> XsecFit:
    """Per-column ellipse from projection widths (squared-width linear LSQ).

    Model: ``w(theta)^2 = c0 + c1 cos 2theta + c2 sin 2theta`` with
    ``a^2 = (c0+R)/4``, ``b^2 = (c0-R)/4``, ``R = hypot(c1, c2)`` — so a >= b
    by construction and the near-circular case degrades to phi = NaN instead
    of a degenerate solve. A column is valid iff all three directions have at
    least one finite width AND the implied ``b^2`` is strictly positive
    (negative values mark the column invalid; they are never clamped).
    """
    W = np.asarray(W, dtype=float)
    n_ang, n_col = W.shape
    theta = np.radians(np.asarray(theta_deg, dtype=float))
    C, solvable, n_angles = _fit_coeffs(W, theta_deg)

    c0, c1, c2 = C
    with np.errstate(invalid="ignore"):
        R = np.hypot(c1, c2)
        a2 = (c0 + R) / 4.0
        b2 = (c0 - R) / 4.0
        valid = solvable & (b2 > 0)
        a = np.where(valid, np.sqrt(np.where(valid, a2, np.nan)), np.nan)
        b = np.where(valid, np.sqrt(np.where(valid, b2, np.nan)), np.nan)
        area = np.pi * a * b
        phi = np.degrees(np.arctan2(c2, c1) / 2.0) % 180.0
        phi = np.where(valid & (R >= 1e-6 * c0), phi, np.nan)

    resid = np.full((n_ang, n_col), np.nan)
    if valid.any():
        design = np.column_stack(
            [np.ones(n_ang), np.cos(2 * theta), np.sin(2 * theta)])
        with np.errstate(invalid="ignore"):
            model2 = design @ C[:, valid]
            model = np.sqrt(np.clip(model2, 0.0, None))
        obs = W[:, valid]
        r = obs - model
        r[~np.isfinite(obs)] = np.nan
        resid[:, valid] = r
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        rms_resid = np.sqrt(np.nanmean(resid ** 2, axis=0))

    return XsecFit(a=a, b=b, phi_deg=phi, area=area, resid=resid,
                   rms_resid=rms_resid, n_angles=n_angles, valid=valid)


def _hexagon_from_h(h: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Closed-form area of the 3-slab hexagon from support distances h (3, N).

    Valid only while every slab binds (h_d < h_e + h_f for all rotations);
    beyond that the closed form UNDER-states the true intersection area and
    would silently break the upper-bound property, so degenerate columns are
    NaN-ed and flagged instead.
    """
    h0, h1, h2 = h
    present = np.isfinite(h).all(axis=0)
    with np.errstate(invalid="ignore"):
        degen = present & ~((h0 < h1 + h2) & (h1 < h0 + h2) & (h2 < h0 + h1))
        area = (2.0 / _SQRT3) * (2.0 * (h0 * h1 + h1 * h2 + h2 * h0)
                                 - (h0 ** 2 + h1 ** 2 + h2 ** 2))
    area = np.where(present & ~degen, area, np.nan)
    return area, degen


def hexagon_area(W: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Circumscribed-hexagon upper bound from the three direction widths.

    Direction width = nanmean of each 180-degree pair (rows k, k+3); support
    distance ``h_d = w_d / 2``. Returns ``(area, degenerate)``; area is NaN
    where a direction is missing or the hexagon is degenerate.
    """
    W = np.asarray(W, dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        w_dir = np.nanmean(np.stack([W[:3], W[3:]]), axis=0)
    return _hexagon_from_h(w_dir / 2.0)


def hexagon_area_expected(a: np.ndarray, b: np.ndarray, phi_deg: np.ndarray
                          ) -> np.ndarray:
    """Hexagon area EXPECTED for a fitted (a, b, phi) ellipse.

    The 0.9069 circle anchor is exact for circles only; a true ellipse gives a
    phi-dependent lower ratio, so QC must compare the measured hexagon against
    this expected value, not against the anchor. NaN phi (near-circular) is
    evaluated at phi = 0, where the support widths are phi-independent anyway.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    phi = np.radians(np.where(np.isfinite(phi_deg), phi_deg, 0.0))
    theta = np.radians(NOMINAL_ANGLES_DEG[:3])[:, None]
    with np.errstate(invalid="ignore"):
        h = np.sqrt(a ** 2 * np.cos(theta - phi) ** 2
                    + b ** 2 * np.sin(theta - phi) ** 2)
    area, _ = _hexagon_from_h(h)
    return area


def split_half_area(W: np.ndarray,
                    theta_deg: np.ndarray = NOMINAL_ANGLES_DEG
                    ) -> tuple[np.ndarray, np.ndarray]:
    """Ellipse areas from the two 3-angle halves (a1-a3 vs a4-a6).

    Each half is a complete 3-direction exact fit, so
    ``area_err = |A_123 - A_456| / 2`` is a conservative per-position
    uncertainty that absorbs 180-degree focus asymmetry.
    """
    W = np.asarray(W, dtype=float)
    lo = W.copy()
    lo[3:] = np.nan
    hi = W.copy()
    hi[:3] = np.nan
    return (fit_ellipse_projections(lo, theta_deg).area,
            fit_ellipse_projections(hi, theta_deg).area)


def pair_differences(W: np.ndarray) -> np.ndarray:
    """(3, N) width differences of the 180-degree pairs: w_k - w_{k+3}."""
    W = np.asarray(W, dtype=float)
    return W[:3] - W[3:]


def predict_anisotropy(W: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Leak-free anisotropy test errors, both (6, N).

    For each (angle k, column) whose 180-degree partner is finite, the
    directional predictor of ``w_k`` is the partner width; the isotropic
    baseline is the mean of the other five widths. Positions with a NaN
    partner are excluded from BOTH errors (same mask, fair comparison).
    """
    W = np.asarray(W, dtype=float)
    n_ang, n_col = W.shape
    dir_err = np.full((n_ang, n_col), np.nan)
    circ_err = np.full((n_ang, n_col), np.nan)
    for k in range(n_ang):
        partner = (k + 3) % n_ang
        others = [j for j in range(n_ang) if j != k]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            base = np.nanmean(W[others], axis=0)
        ok = np.isfinite(W[k]) & np.isfinite(W[partner]) & np.isfinite(base)
        dir_err[k, ok] = W[k, ok] - W[partner, ok]
        circ_err[k, ok] = W[k, ok] - base[ok]
    return dir_err, circ_err


# Minimum spread of the frozen-phi regressor z = cos 2(theta - phi_hat) for
# the two-direction (c0, R) solve to be usable. Two kept directions give
# identical z whenever phi_hat bisects them (theta1 + theta2 = 2 phi_hat mod
# 180 deg) — the design matrix drops to rank 1 and lstsq silently returns a
# minimum-norm answer whose held-out prediction is garbage. The (c0, R) error
# amplification scales as 1/spread; 0.2 skips a ~±3.3 deg band around each
# bisector (for 60-deg-spaced kept pairs) instead of emitting it.
_MIN_Z_SPREAD = 0.2


def predict_phi_transfer(W: np.ndarray,
                         theta_deg: np.ndarray = NOMINAL_ANGLES_DEG
                         ) -> tuple[np.ndarray, np.ndarray, float]:
    """Phi-transfer test: held-out direction prediction with phi fixed.

    Phi is fitted on the ODD grid columns (all six angles, ellipticity-weighted
    coefficient average). At EVEN columns each direction d is dropped (both
    repeats) in turn and ``(c0, R)`` solved from the remaining two directions
    with phi frozen to predict the dropped widths. The solve is only rank 2
    when the kept directions actually differ in ``z = cos 2(theta - phi)`` —
    when ``phi_hat`` (near-)bisects the kept pair the two regressor values
    coincide, so such columns are skipped (NaN) under ``_MIN_Z_SPREAD``
    rather than solved degenerately. Returns ``(err_ellipse, err_circle,
    phi_deg)`` where the error arrays are (6, N), finite only at even columns
    with a well-posed solve; the circle baseline predicts the dropped width
    as the mean of the widths the ellipse solve used. ``phi_deg`` is NaN (and
    both error arrays all-NaN) when no odd column is solvable.
    """
    W = np.asarray(W, dtype=float)
    n_ang, n_col = W.shape
    theta = np.radians(np.asarray(theta_deg, dtype=float))
    dir_ids = _direction_ids(theta_deg)

    C, solvable, _ = _fit_coeffs(W[:, 1::2], theta_deg)
    if not solvable.any():
        return (np.full((n_ang, n_col), np.nan),
                np.full((n_ang, n_col), np.nan), float("nan"))
    s1 = np.nansum(C[1, solvable])
    s2 = np.nansum(C[2, solvable])
    phi_hat = float(np.degrees(np.arctan2(s2, s1) / 2.0) % 180.0)

    z = np.cos(2 * (theta - np.radians(phi_hat)))  # (n_ang,)
    even_idx = np.arange(0, n_col, 2)
    We = W[:, even_idx]
    finite = np.isfinite(We)
    err_ell = np.full((n_ang, n_col), np.nan)
    err_circ = np.full((n_ang, n_col), np.nan)

    code = (finite * (1 << np.arange(n_ang))[:, None]).sum(axis=0)
    n_dirs = int(dir_ids.max()) + 1
    for d in range(n_dirs):
        dropped = [k for k in range(n_ang) if dir_ids[k] == d]
        keep = [k for k in range(n_ang) if dir_ids[k] != d]
        for pattern in np.unique(code):
            rows = [k for k in keep if pattern >> k & 1]
            # rank-2 needs the kept rows to SPREAD in z, not merely to come
            # from two directions (identical z when phi_hat bisects the pair)
            if len(rows) < 2 or float(np.ptp(z[rows])) < _MIN_Z_SPREAD:
                continue
            cols = np.where(code == pattern)[0]
            X = np.column_stack([np.ones(len(rows)), z[rows]])
            coef, *_ = np.linalg.lstsq(X, We[np.ix_(rows, cols)] ** 2,
                                       rcond=None)
            base = We[np.ix_(rows, cols)].mean(axis=0)
            for t in dropped:
                pred2 = coef[0] + coef[1] * z[t]
                pred = np.sqrt(np.clip(pred2, 0.0, None))
                ok = np.isfinite(We[t, cols])
                tgt = even_idx[cols[ok]]
                err_ell[t, tgt] = We[t, cols[ok]] - pred[ok]
                err_circ[t, tgt] = We[t, cols[ok]] - base[ok]
    return err_ell, err_circ, phi_hat


def build_part_stack(fiber: int, part: int, profiles: dict, cfg: CONFIG
                     ) -> PartStack:
    """Align the per-angle profiles of one (fiber, part) into a width stack.

    ``profiles`` maps angle (1..6) to ``{"x": absolute px, "w": widths}``.
    Every profile is first NaN-padded onto the shared ABSOLUTE x grid (the
    per-image CSVs are span-restricted with differing x0 — correlating the
    span-relative arrays would inject a silent x0 misalignment), then shifted
    onto the lowest present angle via ``estimate_shift`` under the
    xsection-specific gates ``cfg.xsec_min_corr`` / ``cfg.xsec_max_shift``.
    A weak or absent correlation peak falls back to zero shift + uncertain;
    so does a peak sitting on the ``±xsec_max_shift`` search boundary
    (``saturated`` flag — the true lag may lie beyond the bound, so the
    clamped estimate is never applied).
    Missing angles become all-NaN rows, never errors. Interior measurement
    dropouts stay NaN in the stack (``max_gap=1.5``: source columns are
    integers + a sub-pixel shift, so consecutive finite samples sit 1 px
    apart and any wider gap means genuinely missing columns — bridging them
    would feed fabricated widths to the ellipse fit and overcount
    ``n_angles``).
    """
    if not profiles:
        raise ValueError(f"no profiles for fiber {fiber} part {part}")
    cfg2 = replace(cfg, min_corr=cfg.xsec_min_corr,
                   max_shift=cfg.xsec_max_shift)

    lo = min(int(np.floor(np.min(p["x"]))) for p in profiles.values())
    hi = max(int(np.ceil(np.max(p["x"]))) for p in profiles.values())
    n_grid = hi - lo + 1
    padded: dict[int, np.ndarray] = {}
    for a, p in profiles.items():
        full = np.full(n_grid, np.nan)
        xi = np.round(np.asarray(p["x"], dtype=float)).astype(int) - lo
        inside = (xi >= 0) & (xi < n_grid)
        full[xi[inside]] = np.asarray(p["w"], dtype=float)[inside]
        padded[a] = full
    grid0 = np.arange(lo, hi + 1, dtype=float)

    ref_angle = min(profiles)
    shifts: list[dict] = []
    aligned: list[tuple[np.ndarray, np.ndarray]] = []
    for a in range(1, 7):
        if a not in profiles:
            shifts.append({"angle": a, "present": False, "shift_px": 0.0,
                           "corr_peak": 0.0, "uncertain": True,
                           "saturated": False})
            aligned.append((np.array([float(lo)]), np.array([np.nan])))
            continue
        if a == ref_angle:
            sh, pk, unc, sat = 0.0, 1.0, False, False
        else:
            sh, pk, unc = estimate_shift(padded[ref_angle], padded[a], cfg2)
            # an argmax on the search boundary means the true peak may lie
            # BEYOND the bound — the clamped lag carries no evidence, so
            # treat it like a weak link (zero shift + uncertain), flagged
            sat = abs(sh) > cfg2.max_shift - 0.5
            if sat:
                sh, unc = 0.0, True
        shifts.append({"angle": a, "present": True, "shift_px": float(sh),
                       "corr_peak": float(pk), "uncertain": bool(unc),
                       "saturated": bool(sat)})
        aligned.append((grid0 + sh, padded[a]))

    grid, stack = resample_to_grid(aligned, max_gap=1.5)
    return PartStack(fiber=fiber, part=part, x=grid.astype(float),
                     W=stack, shifts=shifts)
