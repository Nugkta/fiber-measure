# Metalabbook — study index

Project: `fibrecv` — fibre diameter measurement from optical microscopy images
(per-column edge detection on a desaturation z-map → QC → replicate registration).
This file is the single index of ALL studies; each row points to a
`docs/labbooks/<NN>_<slug>.md` file. Rows sorted by Last Updated, descending.

New studies are only added after the slug + hypothesis have been proposed and
confirmed in chat (see CLAUDE.md, Lab Notebook Protocol).

| Slug | Status | Started | Last Updated | One-line Summary |
|------|--------|---------|--------------|------------------|
| 01_esf_edge_consistency | active | 2026-08-03 | 2026-08-04 | erf-midpoint edge refinement: full MasP2 A/B done — coverage 90% after tuning `refine_relmax` 0.08→0.15, but between-replicate std unchanged (23/46 groups improve, p=1.00); hypothesis unsupported because the replicates are not the same fibre segment (median spread 22.9%); owner flipped `refine_on` default to `False` (opt-in) on 2026-08-04 pending a focus-sweep validation; study kept open pending a proper replicate/focus-sweep set |
