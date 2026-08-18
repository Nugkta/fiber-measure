# report/

Once my folder has changes, please update me.

**One folder per study**, named `<NN>_<slug>/` with the same stem as the study's
labbook file in `docs/labbooks/`. Each folder holds `report.html` (the report:
intro, method, results, discussion — styled, self-contained CSS, open locally in
a browser) plus that study's evidence figures, referenced relatively. Written
when the study's main experiment has run. (Study 01 predates the HTML
convention and keeps its `report.md`.)

## Folders

- `02_multiangle_xsection/` — current (revised 2026-08-16 post-review) —
  study 02 (multi-angle cross-section, C1): `report.html` + six figures
  (`02_xsec_*.png`). Synthetic recovery table; real-data results: all fibers
  elliptical (ratio med 1.127), anisotropy Wilcoxon p=0.048 after the
  alignment fix, A_mean within ±3% of the circular assumption, φ-transfer
  negative (twist), uniformity 0.57–0.82 with six low-confidence fibers;
  scale adjudication and discussion.
- `01_esf_edge_consistency/` — current — study 01 (erf edge refinement):
  `report.md` + `01_esf_ab_summary.png`. M5 A/B validation on the full MasP2
  set: method, tuning path, per-group results, acceptance verdict (coverage
  PASS, replicate-consistency FAIL) and recommended next steps.
