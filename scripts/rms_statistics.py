"""Paired old-vs-new statistics for the per-part ellipse-fit residual (rms).

Why this exists
---------------
The study quotes a per-part ellipse rms improvement with no test, no interval
and no accounting for the fact that the two runs do not keep the same columns.
A residual can fall simply because fewer / easier columns survived, so the
column-count change is measured alongside and correlated against the rms change.

Input
-----
- NEW run: ``<new-root>/per_part/csv/xsec_C1_*_p*.csv`` and
  ``<new-root>/summary/xsection_summary.csv``
  (default ``fibrecv_output/multiangle_c1_s03b`` in this worktree).
- OLD run: the same two paths under
  ``.../worktrees/multiangle-xsection/fibrecv_output/multiangle_c1``.
Parts are paired by (fiber, part) parsed from the filename.

Output
------
``scripts/rms_statistics.json`` (``--out``):
- per-part median ``rms_resid_px`` over ``valid``-true rows, old and new;
- n parts, n improved / worsened / tied;
- Wilcoxon signed-rank statistic and p-value on the paired change;
- median paired change + bootstrap 95% CI on that median (n_resamples, seed
  recorded);
- sign test (exact binomial) on the per-FIBER ``part_rms_med_max_px``
  improvement;
- per-part valid-column-count change: mean |delta| as a percent of part size,
  counts losing >10% / >25%, the worst part, and the Spearman correlation
  between column-count change and rms change.

Pos
---
Study-03 reproducibility artifact. Turns "rms improved" into a paired test with
an interval, and checks the improvement is not an artefact of dropped columns.

Usage
-----
    uv run python scripts/rms_statistics.py

Reminder: once I am updated, update my header comments and the folder's md.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr, wilcoxon

WORKTREES = Path("/Users/stan/Documents/UOM/spins/fiber-measure/.claude/worktrees")
NEW_ROOT = WORKTREES / "edge-criteria" / "fibrecv_output" / "multiangle_c1_s03b"
OLD_ROOT = WORKTREES / "multiangle-xsection" / "fibrecv_output" / "multiangle_c1"

PART_RE = re.compile(r"xsec_C1_(\d+)_p(\d+)\.csv$")


def load_parts(root: Path) -> dict[tuple[int, int], dict]:
    """Per-part median rms over valid rows + valid/total column counts."""
    out: dict[tuple[int, int], dict] = {}
    for p in sorted((root / "per_part" / "csv").glob("xsec_C1_*_p*.csv")):
        m = PART_RE.search(p.name)
        if not m:
            continue
        key = (int(m.group(1)), int(m.group(2)))
        df = pd.read_csv(p)
        valid = df["valid"].astype(bool)
        rms = pd.to_numeric(df.loc[valid, "rms_resid_px"], errors="coerce").dropna()
        out[key] = {
            "file": str(p),
            "n_rows_grid": int(len(df)),
            "n_valid_cols": int(valid.sum()),
            "n_rms_finite": int(len(rms)),
            "rms_med_px": float(rms.median()) if len(rms) else float("nan"),
        }
    return out


def bootstrap_median_ci(deltas: np.ndarray, n_resamples: int, seed: int,
                        alpha: float = 0.05) -> tuple[float, float]:
    """Percentile bootstrap CI on the median of the paired change (resample parts)."""
    rng = np.random.default_rng(seed)
    n = deltas.size
    idx = rng.integers(0, n, size=(n_resamples, n))
    meds = np.median(deltas[idx], axis=1)
    lo, hi = np.percentile(meds, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi)


def _subgroups(frac_old: np.ndarray, d_rms: np.ndarray) -> dict:
    """Is the rms gain concentrated in the parts that dropped columns?

    Splits parts by how many valid columns they lost relative to the old run and
    reports the rms change within each stratum, plus the Spearman correlation
    restricted to the losing parts (where the "fewer, easier columns" confound
    would live).
    """
    out = {}
    for name, mask in (("lost_gt_25pct", frac_old < -0.25),
                       ("lost_10_to_25pct", (frac_old < -0.10) & (frac_old >= -0.25)),
                       ("within_10pct", np.abs(frac_old) <= 0.10),
                       ("gained_gt_10pct", frac_old > 0.10)):
        sel = d_rms[mask]
        out[name] = {
            "n_parts": int(mask.sum()),
            "median_rms_change_px": float(np.median(sel)) if sel.size else None,
            "n_improved": int((sel < 0).sum()),
        }
    losers = frac_old < 0
    if losers.sum() >= 3:
        r, p = spearmanr(frac_old[losers], d_rms[losers])
        out["spearman_among_column_losers"] = {
            "n_parts": int(losers.sum()), "rho": float(r), "p_value": float(p),
            "note": "positive rho = the MORE columns a part lost, the MORE its rms fell "
                    "(i.e. the confound); negative/near-zero rho argues against it",
        }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--new-root", type=Path, default=NEW_ROOT)
    ap.add_argument("--old-root", type=Path, default=OLD_ROOT)
    ap.add_argument("--out", type=Path, default=Path("scripts/rms_statistics.json"))
    ap.add_argument("--n-resamples", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20260819)
    args = ap.parse_args(argv)

    rev = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                         text=True).stdout.strip()

    new, old = load_parts(args.new_root), load_parts(args.old_root)
    keys = sorted(set(new) & set(old))
    only_new, only_old = sorted(set(new) - set(old)), sorted(set(old) - set(new))

    per_part = []
    for k in keys:
        n, o = new[k], old[k]
        d_rms = n["rms_med_px"] - o["rms_med_px"]
        d_cols = n["n_valid_cols"] - o["n_valid_cols"]
        # two denominators, both reported: the part's union-grid length (new run)
        # and the OLD valid count (i.e. "what fraction of the columns we had")
        per_part.append({
            "fiber": k[0], "part": k[1],
            "rms_med_old_px": o["rms_med_px"], "rms_med_new_px": n["rms_med_px"],
            "rms_change_px": d_rms,
            "improved": bool(d_rms < 0),
            "n_valid_cols_old": o["n_valid_cols"], "n_valid_cols_new": n["n_valid_cols"],
            "n_rows_grid_old": o["n_rows_grid"], "n_rows_grid_new": n["n_rows_grid"],
            "d_valid_cols": d_cols,
            "d_cols_frac_of_old_valid": d_cols / o["n_valid_cols"] if o["n_valid_cols"] else float("nan"),
            "d_cols_frac_of_new_grid": d_cols / n["n_rows_grid"] if n["n_rows_grid"] else float("nan"),
        })

    df = pd.DataFrame(per_part)
    ok = df["rms_med_old_px"].notna() & df["rms_med_new_px"].notna()
    d = df.loc[ok, "rms_change_px"].to_numpy(dtype=float)

    w_stat, w_p = wilcoxon(df.loc[ok, "rms_med_new_px"], df.loc[ok, "rms_med_old_px"])
    med_change = float(np.median(d))
    ci_lo, ci_hi = bootstrap_median_ci(d, args.n_resamples, args.seed)

    # ---- per-fiber sign test on part_rms_med_max_px ----
    s_new = pd.read_csv(args.new_root / "summary" / "xsection_summary.csv").set_index("fiber")
    s_old = pd.read_csv(args.old_root / "summary" / "xsection_summary.csv").set_index("fiber")
    fibers = sorted(set(s_new.index) & set(s_old.index))
    fib_rows = []
    for f in fibers:
        o = float(s_old.loc[f, "part_rms_med_max_px"])
        n = float(s_new.loc[f, "part_rms_med_max_px"])
        fib_rows.append({"fiber": int(f), "old_px": o, "new_px": n,
                         "change_px": n - o, "improved": bool(n < o)})
    n_imp = sum(r["improved"] for r in fib_rows)
    n_tie = sum(r["change_px"] == 0.0 for r in fib_rows)
    n_eff = len(fib_rows) - n_tie
    bt = binomtest(n_imp, n_eff, 0.5, alternative="two-sided")

    # ---- column-count vs rms change ----
    dc = df.loc[ok, "d_valid_cols"].to_numpy(dtype=float)
    frac_old = df.loc[ok, "d_cols_frac_of_old_valid"].to_numpy(dtype=float)
    frac_grid = df.loc[ok, "d_cols_frac_of_new_grid"].to_numpy(dtype=float)
    rho, rho_p = spearmanr(dc, d)
    worst_i = int(np.argmin(frac_old))
    worst = df.loc[ok].iloc[worst_i]

    # ---- provenance: what actually differs between the two runs ----
    def _params(root: Path, fname: str) -> dict:
        p = root / "summary" / fname
        if not p.exists():
            return {}
        d = json.loads(p.read_text())
        return d.get("params", d)

    cfg_diff = {}
    for fname in ("run_config.json", "xsection_run_config.json"):
        o, n = _params(args.old_root, fname), _params(args.new_root, fname)
        diff = {k: [o.get(k, "<absent>"), n.get(k, "<absent>")]
                for k in sorted(set(o) | set(n))
                if o.get(k, "<absent>") != n.get(k, "<absent>")}
        cfg_diff[fname] = diff

    result = {
        "kind": "rms_statistics",
        "git_rev": rev,
        "new_root": str(args.new_root), "old_root": str(args.old_root),
        "run_config_diff_old_to_new": cfg_diff,
        "confound_note": "The A/B is not a single-knob change: edge_frac 0.65->0.30, "
                         "k_band 4.0->6.0, edge_cap added, and the bright z-map "
                         "max(R,G,B)->median(R,G,B) all move together. ppu also "
                         "differs (2.5712->1.368) but rms_resid_px is in PIXELS, so "
                         "it does not affect this comparison. The xsection stage "
                         "parameters (xsec_max_shift, xsec_min_corr, xsec_rms_flag_px) "
                         "are unchanged.",
        "pairing": "by (fiber, part) parsed from xsec_C1_<fiber>_p<part>.csv",
        "statistic": "per part: median of rms_resid_px over rows with valid == True",
        "n_parts_paired": int(ok.sum()),
        "parts_only_in_new": [list(k) for k in only_new],
        "parts_only_in_old": [list(k) for k in only_old],
        "rms": {
            "n_improved": int((d < 0).sum()),
            "n_worsened": int((d > 0).sum()),
            "n_tied": int((d == 0).sum()),
            "median_old_px": float(df.loc[ok, "rms_med_old_px"].median()),
            "median_new_px": float(df.loc[ok, "rms_med_new_px"].median()),
            "mean_old_px": float(df.loc[ok, "rms_med_old_px"].mean()),
            "mean_new_px": float(df.loc[ok, "rms_med_new_px"].mean()),
            "median_paired_change_px": med_change,
            "mean_paired_change_px": float(d.mean()),
            "wilcoxon_statistic": float(w_stat),
            "wilcoxon_p_value": float(w_p),
            "wilcoxon_note": "two-sided signed-rank on (new - old), scipy default "
                             "zero_method='wilcox', exact/approx chosen by scipy",
            "bootstrap_median_change_ci95_px": [ci_lo, ci_hi],
            "bootstrap_n_resamples": args.n_resamples,
            "bootstrap_seed": args.seed,
            "bootstrap_method": "percentile bootstrap, parts resampled with "
                                "replacement, statistic = median of the paired change",
        },
        "per_fiber_sign_test": {
            "field": "part_rms_med_max_px from summary/xsection_summary.csv",
            "n_fibers": len(fib_rows),
            "n_improved": int(n_imp),
            "n_worsened": int(n_eff - n_imp),
            "n_tied_excluded": int(n_tie),
            "binomial_p_value_two_sided": float(bt.pvalue),
            "rows": fib_rows,
        },
        "column_counts": {
            "definition": "d_valid_cols = (# rows with valid==True, new) - (same, old)",
            "mean_abs_delta_cols": float(np.abs(dc).mean()),
            "mean_abs_delta_pct_of_old_valid": float(np.abs(frac_old).mean() * 100),
            "mean_abs_delta_pct_of_new_grid": float(np.abs(frac_grid).mean() * 100),
            "mean_signed_delta_cols": float(dc.mean()),
            "median_signed_delta_cols": float(np.median(dc)),
            "n_parts_lost_gt_10pct_of_old_valid": int((frac_old < -0.10).sum()),
            "n_parts_lost_gt_25pct_of_old_valid": int((frac_old < -0.25).sum()),
            "n_parts_gained_gt_10pct_of_old_valid": int((frac_old > 0.10).sum()),
            "worst_part": {
                "fiber": int(worst["fiber"]), "part": int(worst["part"]),
                "n_valid_cols_old": int(worst["n_valid_cols_old"]),
                "n_valid_cols_new": int(worst["n_valid_cols_new"]),
                "d_valid_cols": int(worst["d_valid_cols"]),
                "lost_pct_of_old_valid": float(worst["d_cols_frac_of_old_valid"] * 100),
                "rms_med_old_px": float(worst["rms_med_old_px"]),
                "rms_med_new_px": float(worst["rms_med_new_px"]),
                "rms_change_px": float(worst["rms_change_px"]),
            },
            "spearman_dcols_vs_drms": {"rho": float(rho), "p_value": float(rho_p),
                                       "note": "positive rho = parts that GAINED "
                                               "columns also got WORSE rms"},
            "subgroups": _subgroups(frac_old, d),
        },
        "per_part": per_part,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2))

    r = result["rms"]
    print(f"parts paired: {result['n_parts_paired']}  "
          f"improved={r['n_improved']} worsened={r['n_worsened']} tied={r['n_tied']}")
    print(f"rms p50: old={r['median_old_px']:.3f} -> new={r['median_new_px']:.3f} px")
    print(f"median paired change = {r['median_paired_change_px']:+.3f} px  "
          f"95% CI [{r['bootstrap_median_change_ci95_px'][0]:+.3f}, "
          f"{r['bootstrap_median_change_ci95_px'][1]:+.3f}] "
          f"({r['bootstrap_n_resamples']} resamples, seed {r['bootstrap_seed']})")
    print(f"Wilcoxon W={r['wilcoxon_statistic']:.1f}  p={r['wilcoxon_p_value']:.3e}")
    st = result["per_fiber_sign_test"]
    print(f"per-fiber sign test: {st['n_improved']}/{st['n_fibers']} improved, "
          f"p={st['binomial_p_value_two_sided']:.4f}")
    cc = result["column_counts"]
    print(f"columns: mean |delta| = {cc['mean_abs_delta_cols']:.1f} cols "
          f"({cc['mean_abs_delta_pct_of_old_valid']:.2f}% of old valid, "
          f"{cc['mean_abs_delta_pct_of_new_grid']:.2f}% of new grid); "
          f"lost>10%: {cc['n_parts_lost_gt_10pct_of_old_valid']}, "
          f"lost>25%: {cc['n_parts_lost_gt_25pct_of_old_valid']}")
    w = cc["worst_part"]
    print(f"worst part: C1_{w['fiber']:02d}_p{w['part']}  "
          f"{w['n_valid_cols_old']} -> {w['n_valid_cols_new']} cols "
          f"({w['lost_pct_of_old_valid']:+.1f}%), rms {w['rms_med_old_px']:.2f} -> "
          f"{w['rms_med_new_px']:.2f} px")
    sp = cc["spearman_dcols_vs_drms"]
    print(f"Spearman(d_cols, d_rms) rho={sp['rho']:+.3f} p={sp['p_value']:.3e}")
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
