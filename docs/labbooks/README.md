# labbooks/

Once my folder has changes, please update me.

One markdown file per study, named `<NN>_<slug>.md` (`<NN>` = next unused
two-digit number here). Each file: hypothesis/status/started header, then one
`## YYYY-MM-DD` entry per day, newest first. Indexed by `docs/metalabbook.md`.

## Files

- `02_multiangle_xsection.md` — done — multi-angle cross-section study (C1):
  survey + angle-structure identification, scale adjudication (scale-bar →
  `Scaling/Items` 0.3889 µm/px), full implementation log of pipeline stage 3
  (`run_xsection`), synthetic recovery, real-data validation and close-out;
  report in `docs/report/02_multiangle_xsection/report.html`.
- `01_esf_edge_consistency.md` — active — erf-midpoint edge refinement study: validated the blurred-step (PSF) model on MasP2 images, then A/B-tested replicate consistency of the refined edges over the full 141-image set (null result — see `docs/report/01_esf_edge_consistency/report.md`); kept open pending a better replicate set.
- `01_esf_fits.png` — evidence figure for study 01 (edge-aligned profiles, erf fits, residual clouds on three MasP2 images).
