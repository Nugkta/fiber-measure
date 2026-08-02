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
- `02_gui-redesign.md` — current — "clean lab" visual redesign of the
  Streamlit GUI: violet-accent theme, numbered card sections + jump menu,
  shared `_fmt`/`_styled_fig` formatting/plot policy; detection logic and
  exports unchanged.
