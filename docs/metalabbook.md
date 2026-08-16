# Metalabbook — study index

Project: `fibrecv` — fibre diameter measurement from optical microscopy images
(per-column edge detection on a desaturation z-map → QC → replicate registration).
This file is the single index of ALL studies; each row points to a
`docs/labbooks/<NN>_<slug>.md` file. Rows sorted by Last Updated, descending.

New studies are only added after the slug + hypothesis have been proposed and
confirmed in chat (see CLAUDE.md, Lab Notebook Protocol).

| Slug | Status | Started | Last Updated | One-line Summary |
|------|--------|---------|--------------|------------------|
| 02_multiangle_xsection | done | 2026-08-15 | 2026-08-16 | multi-angle cross-section (C1: 15 fibers × 6 angles × 5 parts) SHIPPED as pipeline stage 3 (`run_xsection`: w²-space ellipse fit + hexagon QC + A_mean/A_harm/A_min): scale adjudicated to `Scaling/Items` 0.3889 µm/px by scale-bar (camera 0.220 wrong, 1.77×); all fibers elliptical (ratio med ≈1.13, φ persists along fibers, Fisher p=1.4e-4) BUT A_mean within ±3% of π(d̄/2)² — circular area OK for C1, A_min/A_mean 0.71–0.92 is the new information; anisotropy Wilcoxon marginal (p=0.076, focus systematics dominate); v2 angle offsets not needed; report `docs/report/02_multiangle_xsection.md`, code on `worktree-multiangle-xsection` |
| 01_esf_edge_consistency | done | 2026-08-03 | 2026-08-04 | erf-midpoint edge refinement: accurate on synthetic ground truth but **rejected on real MasP2** — A/B (144 images, 46 groups) showed no consistency gain (23/46, p=1.00), and the owner's visual acceptance (2026-08-04) found two structural defects: boundary bites inside the visible edge (50%-midpoint assumes a step edge; a cylinder's profile is a shouldered ramp) and 1.55× median column jitter from unsmoothed 16-col block fits (rougher in 140/144 images); `refine_on=False` stays the default, branch `worktree-esf-edge-refinement` pushed and kept unmerged as the failed-attempt record; full report in `docs/report/01_esf_edge_consistency.md` |
