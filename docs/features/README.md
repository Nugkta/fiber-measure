# features/

Once my folder has changes, please update me.

One markdown file per substantial new feature (new behaviour or algorithm —
not bugfixes or refactors), named `<NN>_<slug>.md` (`<NN>` = next unused
two-digit number here, independent of the labbook numbering). Template:
What it does / Design choices / Algorithm details / Caveats (see CLAUDE.md).

## Files

- `01_tilt_corrected_diameter.md` — current — perpendicular (tilt-corrected)
  diameter measurement: why the vertical chord over-read, the four
  compensations, design decisions and caveats.
- `03_multiangle_xsection.md` — current — third pipeline stage for
  multi-angle image sets: w²-space per-column ellipse fit, cross-angle
  alignment gates, Zeiss-XML µm scale (adjudicated `Scaling/Items`),
  A_mean/A_harm/A_min export and the anisotropy/φ-transfer validation
  design.
- `02_gui-redesign.md` — current — "clean lab" visual redesign of the
  Streamlit GUI: violet-accent theme, numbered card sections + jump menu,
  shared `_fmt`/`_styled_fig` formatting/plot policy; detection logic and
  exports unchanged.
- `04_anomaly_flags.md` — current — advisory anomaly flagging: edge_jump /
  large_gap / diameter_step per image + replicate_outlier per group, the
  `anomaly_exclude` switch, `"anomaly"` meta sub-dict and
  `summary/per_image_summary.csv`.
- `05_edge_criteria.md` — current — bright-mode edge criteria upgrade:
  median(R,G,B) z-map (rejects single-channel chromatic-aberration fringes)
  + clamped relative edge threshold (`edge_frac`/`edge_cap`, prevents
  inward bite at 50%) + coupled `k_band` 4→6 recalibration; desat mode
  (MasP2) unchanged (bit-identical).
- `06_gui_multiangle.md` — current — GUI multi-angle cross-section mode:
  sidebar analysis-mode radio, per-angle tabs, cross-section panel
  (ellipse fit + split-half uncertainty), condition-scoped batch card;
  manual µm/px scale and numeric `--scale-source` bypassing XML sidecars.
