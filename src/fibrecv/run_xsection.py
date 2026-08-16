"""CLI: third stage — per-position cross-sections from multi-angle profiles.

Dependencies
------------
``numpy``, ``pandas``, ``matplotlib`` (Agg), the ``xsection`` + ``scale`` +
``io_utils`` + ``config`` modules. Reads the per-image CSV/meta artifacts of
``run_measure`` (multi-angle names) plus the Zeiss XML sidecars in the data dir.

Inputs
------
``--out`` (output root holding ``per_image/``), ``--data-root`` (image dir with
``*_metadata.xml`` sidecars), optional ``--condition`` / ``--fibers`` filters,
``--scale-source camera|items|<um_per_px>`` (µm conversion; the scale-bar cross
-check adjudicated "items" for C1 — labbook 02), ``--min-corr`` /
``--max-shift`` (cross-angle alignment gates -> ``xsec_min_corr`` /
``xsec_max_shift``), ``--validation``.

Output
------
Per (fiber, part): ``per_part/csv/xsec_<cond>_<ff>_p<p>.csv`` (aligned widths +
per-column ellipse/hexagon), ``per_part/plots/...png``, ``per_part/shifts/...json``
(alignment + both XML scales). Per run: ``summary/xsection_summary.csv`` (one
row per fiber), ``summary/xsection_angle_residuals.csv`` (v2-gate evidence),
``summary/xsection_run_config.json`` and, with ``--validation``,
``summary/xsection_validation.csv`` (anisotropy + phi-transfer errors).

Pos
---
Third entrypoint: run_measure -> (run_aggregate | run_xsection). This is the
aggregation stage for multi-angle (C1) data — ``run_aggregate`` skips those
names by design. All math is delegated to the pure ``xsection`` module; µm
conversion happens here, once, from the per-image XML scale (never cfg.ppu).

Once I am updated, update my header comments and folder's md.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .anomaly import exclusion_reason  # noqa: E402
from .config import CONFIG  # noqa: E402
from .io_utils import parse_multiangle_name  # noqa: E402
from .scale import read_scale, resolve_um_per_px, sidecar_for  # noqa: E402
from .xsection import (  # noqa: E402
    build_part_stack,
    fit_ellipse_projections,
    hexagon_area,
    hexagon_area_expected,
    pair_differences,
    predict_anisotropy,
    predict_phi_transfer,
    split_half_area,
)

DEFAULT_OUT = "./fibrecv_output"
_SCALE_TOL = 1e-3  # sidecars of one (fiber, part) must agree to 0.1%


def _load_multiangle_profiles(out_root: Path, cfg: CONFIG,
                              condition: str | None):
    """Bucket per-image profiles by (condition, fiber, part).

    Returns ``(parts, seen_angles)`` where ``parts[key][angle]`` holds the
    NaN-masked smooth width profile of an INCLUDED image and
    ``seen_angles[key]`` lists every measured angle (excluded ones too — their
    sidecars still participate in the scale-consistency check while their
    widths become NaN rows).
    """
    csv_dir = out_root / "per_image" / "csv"
    meta_dir = out_root / "per_image" / "diagnostics"
    parts: dict = defaultdict(dict)
    seen_angles: dict = defaultdict(set)
    for csv in sorted(csv_dir.glob("*_profile.csv")):
        base = csv.stem[: -len("_profile")]
        try:
            key = parse_multiangle_name(base)
        except ValueError:
            continue
        if key.scalebar:
            continue
        if condition is not None and key.condition != condition:
            continue
        kk = (key.condition, key.fiber, key.part)
        seen_angles[kk].add(key.angle)

        coverage = None
        band_mismatch = False
        flags: list = []
        meta_path = meta_dir / f"{base}_meta.json"
        if meta_path.exists():
            with open(meta_path) as fh:
                meta = json.load(fh)
            coverage = meta.get("coverage")
            band_mismatch = bool(meta.get("band_mismatch"))
            flags = list((meta.get("anomaly") or {}).get("flags") or [])
        if exclusion_reason(band_mismatch, coverage, flags, cfg) is not None:
            continue

        df = pd.read_csv(csv)
        w = np.where(df["valid"].to_numpy(bool),
                     df["diameter_px_smooth"].to_numpy(float), np.nan)
        parts[kk][key.angle] = {"x": df["x_px"].to_numpy(float), "w": w}
    return parts, seen_angles


def _part_scale(data_root: Path, cond: str, fiber: int, part: int,
                angles, scale_source):
    """Read every sidecar of the part; enforce 0.1% agreement; resolve µm/px."""
    infos = []
    for a in sorted(angles):
        img = data_root / f"{cond}_{fiber:02d}_a{a}_part{part}.tiff"
        infos.append(read_scale(sidecar_for(img)))
    for attr in ("um_per_px_camera", "um_per_px_items"):
        vals = np.array([getattr(i, attr) for i in infos])
        if vals.max() / vals.min() - 1.0 > _SCALE_TOL:
            raise ValueError(
                f"{cond}_{fiber:02d} part{part}: sidecars disagree on {attr} "
                f"({vals.min():.6g}..{vals.max():.6g})")
    return infos[0], float(resolve_um_per_px(infos[0], scale_source))


def _rolling_median(a: np.ndarray, window: int = 11) -> np.ndarray:
    return (pd.Series(a).rolling(window, center=True, min_periods=window // 2 + 1)
            .median().to_numpy())


def _circular_phi_median(phi_deg: np.ndarray) -> float:
    """Orientation average with period 180 deg (vector mean of 2*phi)."""
    ph = np.radians(phi_deg[np.isfinite(phi_deg)])
    if ph.size == 0:
        return float("nan")
    return float(np.degrees(np.arctan2(np.mean(np.sin(2 * ph)),
                                       np.mean(np.cos(2 * ph))) / 2) % 180.0)


def _rmse(err: np.ndarray) -> tuple[float, int]:
    finite = np.isfinite(err)
    n = int(finite.sum())
    if n == 0:
        return float("nan"), 0
    return float(np.sqrt(np.mean(err[finite] ** 2))), n


def _plot_part(df: pd.DataFrame, name: str, um: float, out: Path) -> None:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    for a in range(1, 7):
        ax1.plot(df["x_px"], df[f"w_a{a}_px"] * um, lw=0.8, label=f"a{a}")
    ax1.set_ylabel("aligned width (µm)")
    ax1.legend(loc="best", fontsize=7, ncol=6)
    ax1.grid(alpha=0.3)
    ax1.set_title(name)
    area = df["area_um2"]
    err = df["area_err_um2"].fillna(0.0)
    ax2.plot(df["x_px"], area, lw=1.2, color="tab:blue", label="ellipse area")
    ax2.fill_between(df["x_px"], area - err, area + err, alpha=0.25,
                     color="tab:blue", label="±split-half err")
    ax2.plot(df["x_px"], df["area_hex_um2"], lw=0.8, color="tab:orange",
             label="hexagon bound")
    ax2.set_xlabel("x (px, ref-angle frame)")
    ax2.set_ylabel("area (µm²)")
    ax2.legend(loc="best", fontsize=7)
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Per-position cross-sections from multi-angle profiles.")
    ap.add_argument("--out", default=DEFAULT_OUT,
                    help="output root (holds per_image/)")
    ap.add_argument("--data-root", required=True,
                    help="image directory with the *_metadata.xml sidecars")
    ap.add_argument("--condition", default=None,
                    help="condition label filter, e.g. C1")
    ap.add_argument("--fibers", nargs="+", type=int, default=None,
                    help="fiber numbers to process (default: all)")
    ap.add_argument("--scale-source", dest="scale_source", default="items",
                    help='µm/px source: "camera", "items" (default; scale-bar '
                         'adjudicated) or a numeric µm/px literal')
    ap.add_argument("--min-corr", dest="min_corr", type=float, default=None,
                    help="cross-angle corr-peak gate (xsec_min_corr)")
    ap.add_argument("--max-shift", dest="max_shift", type=int, default=None,
                    help="cross-angle lag bound px (xsec_max_shift)")
    ap.add_argument("--validation", action="store_true",
                    help="also write summary/xsection_validation.csv")
    args = ap.parse_args(argv)

    cfg = CONFIG()
    ov = {k: v for k, v in {
        "xsec_min_corr": args.min_corr, "xsec_max_shift": args.max_shift,
    }.items() if v is not None}
    cfg = replace(cfg, **ov)

    out_root = Path(args.out)
    data_root = Path(args.data_root)
    parts, seen_angles = _load_multiangle_profiles(out_root, cfg,
                                                   args.condition)
    if args.fibers:
        wanted = set(args.fibers)
        parts = {k: v for k, v in parts.items() if k[1] in wanted}
    if not parts:
        print("No multi-angle profiles matched.")
        return 1

    for d in ("per_part/csv", "per_part/plots", "per_part/shifts", "summary"):
        (out_root / d).mkdir(parents=True, exist_ok=True)

    by_fiber: dict = defaultdict(list)
    resid_rows: list[dict] = []
    val_rows: list[dict] = []

    for (cond, fiber, part) in sorted(parts):
        profiles = parts[(cond, fiber, part)]
        if not profiles:
            continue
        try:
            info, um = _part_scale(data_root, cond, fiber, part,
                                   seen_angles[(cond, fiber, part)],
                                   args.scale_source)
        except (ValueError, FileNotFoundError, OSError) as exc:
            print(f"[ERROR] scale: {exc}")
            return 2

        st = build_part_stack(fiber, part, profiles, cfg)
        fit = fit_ellipse_projections(st.W)
        hex_px2, hex_degen = hexagon_area(st.W)
        hex_exp_px2 = hexagon_area_expected(fit.a, fit.b, fit.phi_deg)
        A_lo, A_hi = split_half_area(st.W)
        with np.errstate(invalid="ignore"):
            area_err_px2 = np.abs(A_lo - A_hi) / 2.0
            hex_ratio = fit.area / hex_px2
            hex_ratio_expected = fit.area / hex_exp_px2

        name = f"xsec_{cond}_{fiber:02d}_p{part}"
        df = pd.DataFrame({
            "x_px": st.x,
            **{f"w_a{a}_px": st.W[a - 1] for a in range(1, 7)},
            "n_angles": fit.n_angles,
            "a_px": fit.a,
            "b_px": fit.b,
            "phi_deg": fit.phi_deg,
            "rms_resid_px": fit.rms_resid,
            "a_um": fit.a * um,
            "b_um": fit.b * um,
            "area_um2": fit.area * um ** 2,
            "area_err_um2": area_err_px2 * um ** 2,
            "area_hex_um2": hex_px2 * um ** 2,
            "hex_ratio": hex_ratio,
            "hex_ratio_expected": hex_ratio_expected,
            "valid": fit.valid.astype(int),
        })
        df.to_csv(out_root / "per_part" / "csv" / f"{name}.csv", index=False)
        _plot_part(df, name, um, out_root / "per_part" / "plots" / f"{name}.png")
        with open(out_root / "per_part" / "shifts" / f"{name}.json", "w") as fh:
            json.dump({
                "fiber": fiber, "part": part, "condition": cond,
                "shifts": st.shifts,
                "um_per_px_camera": info.um_per_px_camera,
                "um_per_px_items": info.um_per_px_items,
                "um_per_px_resolved": um,
                "n_columns": int(st.x.size),
                "n_valid": int(fit.valid.sum()),
                "n_hex_degenerate": int(hex_degen.sum()),
                "n_saturated": int(sum(1 for s in st.shifts
                                       if s.get("saturated"))),
            }, fh, indent=2)

        for k in range(6):
            r = fit.resid[k]
            finite = np.isfinite(r)
            resid_rows.append({
                "fiber": fiber, "angle": k + 1,
                "part": part,
                "n_cols": int(finite.sum()),
                "median_signed_resid_px": float(np.median(r[finite]))
                if finite.any() else float("nan"),
                "mad_resid_px": float(np.median(np.abs(
                    r[finite] - np.median(r[finite]))))
                if finite.any() else float("nan"),
                "corr_peak": st.shifts[k]["corr_peak"],
                "shift_px": st.shifts[k]["shift_px"],
                "uncertain": st.shifts[k]["uncertain"],
            })

        if args.validation:
            dir_err, circ_err = predict_anisotropy(st.W)
            for k in range(6):
                rd, n = _rmse(dir_err[k])
                rc_, _ = _rmse(circ_err[k])
                val_rows.append({
                    "fiber": fiber, "part": part, "kind": "anisotropy",
                    "index": k + 1, "n": n, "rmse_directional": rd,
                    "rmse_circle": rc_, "phi_used": float("nan")})
            t_ell, t_circ, phi_hat = predict_phi_transfer(st.W)
            for d in range(3):
                mask_rows = [d, d + 3]
                rd, n = _rmse(t_ell[mask_rows])
                rc_, _ = _rmse(t_circ[mask_rows])
                val_rows.append({
                    "fiber": fiber, "part": part, "kind": "phi_transfer",
                    "index": d + 1, "n": n, "rmse_directional": rd,
                    "rmse_circle": rc_, "phi_used": phi_hat})

        pair_dw = pair_differences(st.W)
        with np.errstate(invalid="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                w_dir = np.nanmean(np.stack([st.W[:3], st.W[3:]]), axis=0)
                dw_frac = np.abs(pair_dw) / w_dir
        n_unc = sum(1 for s in st.shifts
                    if s["present"] and s["uncertain"])
        n_sat = sum(1 for s in st.shifts if s.get("saturated"))
        n_links = sum(1 for s in st.shifts if s["present"]) - 1
        rms_med = (float(np.nanmedian(fit.rms_resid[fit.valid]))
                   if fit.valid.any() else float("nan"))
        by_fiber[fiber].append({
            "part": part, "um": um,
            "x": st.x, "W_um": st.W * um,
            "area_um2": df["area_um2"].to_numpy(),
            "valid": fit.valid,
            "ratio": fit.a / fit.b,
            "phi": fit.phi_deg,
            "hex_ratio": hex_ratio,
            "dw_frac": dw_frac,
            "n_uncertain": n_unc,
            "n_saturated": n_sat,
            "rms_med": rms_med,
            "n_links": max(n_links, 0),
        })
        print(f"  {name}: cols={st.x.size} valid={int(fit.valid.sum())} "
              f"uncertain_shifts={n_unc} saturated={n_sat}")

    # ---- per-fiber summary --------------------------------------------
    rows = []
    for fiber in sorted(by_fiber):
        ps = by_fiber[fiber]
        area = np.concatenate([p["area_um2"][p["valid"]] for p in ps])
        area = area[np.isfinite(area)]
        n_pos = int(area.size)
        A_mean = float(np.mean(area)) if n_pos else float("nan")
        A_harm = (float(n_pos / np.sum(1.0 / area))
                  if n_pos and (area > 0).all() else float("nan"))
        # weakest link on the 11-px rolling-median-filtered area profile
        A_min, min_part, min_x = float("nan"), -1, float("nan")
        for p in ps:
            filt = _rolling_median(
                np.where(p["valid"], p["area_um2"], np.nan))
            if np.isfinite(filt).any():
                i = int(np.nanargmin(filt))
                if not np.isfinite(A_min) or filt[i] < A_min:
                    A_min, min_part, min_x = float(filt[i]), p["part"], float(p["x"][i])
        ratio = np.concatenate([p["ratio"][p["valid"]] for p in ps])
        ratio = ratio[np.isfinite(ratio)]
        phi_all = np.concatenate([p["phi"] for p in ps])
        hexr = np.concatenate([p["hex_ratio"] for p in ps])
        hexr = hexr[np.isfinite(hexr)]
        dwf = np.concatenate([p["dw_frac"].ravel() for p in ps])
        dwf = dwf[np.isfinite(dwf)]
        W_um_all = np.concatenate([p["W_um"].ravel() for p in ps])
        d_bar = float(np.nanmean(W_um_all))
        A_circle = np.pi * (d_bar / 2.0) ** 2
        n_cols_total = int(sum(p["valid"].size for p in ps))
        n_unc = int(sum(p["n_uncertain"] for p in ps))
        n_sat = int(sum(p["n_saturated"] for p in ps))
        n_links = int(sum(p["n_links"] for p in ps))
        valid_frac = n_pos / n_cols_total if n_cols_total else 0.0
        part_rms = np.array([p["rms_med"] for p in ps])
        rms_max = (float(np.nanmax(part_rms))
                   if np.isfinite(part_rms).any() else float("nan"))
        rows.append({
            "fiber": fiber,
            "n_parts": len(ps),
            "n_positions": n_pos,
            "um_per_px": float(np.mean([p["um"] for p in ps])),
            "A_mean_um2": A_mean,
            "A_harm_um2": A_harm,
            "A_min_um2": A_min,
            "A_min_part": min_part,
            "A_min_x_px": min_x,
            "axis_ratio_median": float(np.median(ratio)) if ratio.size else float("nan"),
            "axis_ratio_iqr": float(np.percentile(ratio, 75)
                                    - np.percentile(ratio, 25)) if ratio.size else float("nan"),
            "phi_med_deg": _circular_phi_median(phi_all),
            "hex_ratio_median": float(np.median(hexr)) if hexr.size else float("nan"),
            "pair_dw_frac_median": float(np.median(dwf)) if dwf.size else float("nan"),
            "A_circle_um2": A_circle,
            "area_ratio": A_mean / A_circle if np.isfinite(A_mean) else float("nan"),
            "uniformity": A_min / A_mean if np.isfinite(A_min) and A_mean else float("nan"),
            "n_uncertain_shifts": n_unc,
            "n_saturated_shifts": n_sat,
            "part_rms_med_max_px": rms_max,
            "low_confidence": bool(
                valid_frac < 0.5
                or (n_links and n_unc / n_links > 1 / 3)
                or (np.isfinite(rms_max)
                    and rms_max > cfg.xsec_rms_flag_px)),
        })
        print(f"fiber {fiber:02d}: A_mean={A_mean:.1f}um2 "
              f"ratio_med={rows[-1]['axis_ratio_median']:.3f} "
              f"area_ratio={rows[-1]['area_ratio']:.3f} "
              f"uncertain={n_unc}/{n_links}")

    pd.DataFrame(rows).to_csv(
        out_root / "summary" / "xsection_summary.csv", index=False)

    # per (fiber, angle) residual diagnostics, aggregated over parts
    rr = pd.DataFrame(resid_rows)
    agg = (rr.groupby(["fiber", "angle"])
           .agg(n_cols=("n_cols", "sum"),
                median_signed_resid_px=("median_signed_resid_px", "median"),
                mad_resid_px=("mad_resid_px", "median"),
                mean_corr_peak=("corr_peak", "mean"),
                n_uncertain=("uncertain", "sum"))
           .reset_index())
    agg.to_csv(out_root / "summary" / "xsection_angle_residuals.csv",
               index=False)

    if args.validation:
        pd.DataFrame(val_rows).to_csv(
            out_root / "summary" / "xsection_validation.csv", index=False)

    with open(out_root / "summary" / "xsection_run_config.json", "w") as fh:
        json.dump({
            "params": cfg.as_dict(),
            "scale_source": str(args.scale_source),
            "data_root": str(data_root),
            "condition": args.condition,
            "fibers": args.fibers,
            "n_parts": len(parts),
        }, fh, indent=2)

    print(f"Wrote {len(rows)} fiber rows -> "
          f"{out_root / 'summary' / 'xsection_summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
