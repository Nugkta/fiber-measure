# fibrecv GUI — local desktop app (Mac / Windows)

A local [Streamlit](https://streamlit.io) web app for tuning the fibre-diameter
detector, previewing detected boundaries, and batch-processing + exporting — all
from your own machine. It runs **entirely locally**: it reads images from a
folder you point it at and writes outputs to a folder you choose. No cluster
access is needed at runtime.

The app is a thin front-end over the validated `fibrecv` pipeline — it changes
**no** detection logic and uses the same calibrated defaults as the CLI.

It uses a "clean lab" light theme (violet accent, slim header with a
one-line status chip, numbered sections with a jump menu on the right)
defined by `.streamlit/config.toml` plus the app's own styling code — purely
visual, see
[docs/features/02_gui-redesign.md](docs/features/02_gui-redesign.md). A
sidebar **Analysis** switch picks between two main-area layouts: **Replicates**
(four sections — the default, described first below) and **Multi-angle
cross-section** (three sections, for C1-style rotation-angle image sets —
described in its own part further down), see
[docs/features/06_gui_multiangle.md](docs/features/06_gui_multiangle.md).

---

## 1. One-time setup

You need **Python 3.12+** and the `fibrecv` project folder plus a folder of
microscopy images on your machine.

Open a terminal **in the `fibrecv` folder** (the one containing `pyproject.toml`):

### Option A — uv (recommended)
```bash
# install uv once: https://docs.astral.sh/uv/getting-started/installation/
uv sync                 # creates .venv and installs everything incl. streamlit
```

### Option B — pip + virtualenv
```bash
python -m venv .venv
# macOS / Linux:
source .venv/bin/activate
# Windows (PowerShell):
.venv\Scripts\Activate.ps1

pip install -e .        # installs fibrecv + streamlit and the fibrecv-gui command
```

---

## 2. Launch

From the `fibrecv` folder:

```bash
# uv:
uv run fibrecv-gui
# or, with the venv activated (pip route):
fibrecv-gui
```

Equivalent, if you prefer to call Streamlit directly:
```bash
uv run streamlit run src/fibrecv/gui_app.py        # uv
streamlit run src/fibrecv/gui_app.py               # activated venv
```

The app opens automatically at **http://localhost:8501**. If it doesn't open,
paste that URL into your browser. Stop the app with `Ctrl+C` in the terminal.

To use a different port: `fibrecv-gui --server.port 8600` (any extra arguments
are passed straight through to `streamlit run`).

---

## 3. Using the app

**Sidebar — Analysis** (mode)
- A radio at the top of the sidebar picks **Replicates** (default) or
  **Multi-angle cross-section**. Replicates mode splits the main area into
  four numbered sections — 01 Replicates, 02 Group panel, 03 Tensile,
  04 Export & batch — described below in order. Multi-angle mode replaces
  all four with three different sections — 01 Angles, 02 Cross-section,
  03 Batch & export — described in their own part at the end of this
  section. A fixed jump menu on the right (hidden on narrow windows) links
  straight to whichever set of sections is showing, once data is loaded.
  Entering multi-angle mode once applies that mode's calibrated `bright`
  image-mode defaults (`edge_frac`/`k_band`); switching back to Replicates
  does not revert them automatically.

The rest of this section (Data / Detection sidebar, 01–04) describes
**Replicates mode**; multi-angle mode is described afterwards under
"Multi-angle cross-section mode".

**Sidebar — Data** (source)
- **Local folder**: type the path to your images folder. The app scans for
  image files (`.jpg/.jpeg/.png/.tif/.tiff/.bmp`) and groups them by the
  numbers at the end of each name: the last number is the replicate, the
  numbers before it form the group (`masp2 3_1_2.jpg` → group `3_1`,
  replicate 2; `3-1-2.png` and `sampleA 3_1_2.tif` work too). Pick a group
  from the **Group** dropdown; files whose names don't end in numbers appear
  under **ungrouped**.
- **Upload**: drop any number of image files. They are grouped by the same
  naming rule (a Group dropdown appears if there is more than one group);
  export needs names ending in numbers to derive a group label.

**Sidebar — Detection** (parameters)
- Four boundary/calibration knobs plus the anomaly-flag section below are
  exposed; everything else uses the validated defaults (the CLI retains full
  parameter control):
  - **`edge_z`** (slider) — where on the fibre wall the boundary is drawn:
    higher = higher up the wall = further inside the fibre = thinner reading;
    lower = thicker. Line cutting into the fibre → lower it; line sitting in
    the shadow outside → raise it.
  - **`edge_frac`** — faint-fibre protection: caps the crossing level at a
    fraction of the wall's own height so a weak wall keeps its boundary. It
    only acts when the wall is weaker than `edge_z`.
  - **`wcol`** — horizontal smoothing width (px): raise it for a smoother,
    more stable line; lower it to preserve fine thickness variation.
    Default 41.
  - **`ppu`** (Calibration) — camera pixels per micron; every µm value is
    px / ppu. Change it if your images come from a different
    microscope/magnification.
- A collapsed **Anomaly flags** expander holds the four advisory-warning thresholds
  (edge-jump px, missing-stretch fraction, diameter-step fraction, replicate
  deviation fraction) and the **exclude** checkbox that makes flagged images
  drop out of the group statistics (off by default — flags are advisory).
  Flagged replicates get a `⚠` tab prefix, the anomaly names in the flags
  badge and an `anomalies` column in the per-image stats table; see
  `docs/features/04_anomaly_flags.md`.
- Edits are **staged**: change what you want, then click **Apply** to re-render
  (a few seconds for a group). **Reset to defaults** restores the calibrated
  values.
- Below Data/Detection, a collapsible **Tensile** expander holds the tensile
  data source (folder or upload, same pattern as the image source), the
  gauge length and the modulus-fit window — see **03 Tensile** below.

**01 Replicates**
- One tab per replicate: the full-resolution **overlay** (cyan top edge / yellow
  bottom edge / dashed centerline / green perpendicular measurement chords —
  each green line is the exact diameter being reported at that column, so its
  far end must touch the bottom edge), the **diameter-vs-position** plot (raw
  points + smoothed line, µm), and scalar metrics (median Ø, coverage, tilt,
  QC flags).
- **Edit boundaries (manual correction)**: when detection fails locally, open
  the expander in a replicate tab, pick the top or bottom line, and click 2+
  points along the true edge in the zoomed strip — the line is redrawn through
  your points (blended at the ends) and shown in **magenta**. Points are
  grouped into **sets**: one set corrects one stretch, and clicking far away
  inside the same set still connects the points — so use **Start new set** to
  fix a different stretch independently (the line between sets is left
  untouched). The set radio picks which set new clicks extend (active set's
  markers are white, others grey); **Undo last point** / **Delete set** manage
  them. Whole-line nudge inputs handle uniform offsets. Corrections re-run QC
  and flow into the profile plot, the group statistics and **Export current
  group**; they survive parameter changes, but **Run batch** recomputes from
  disk and ignores them.

**02 Group panel**
- Every replicate's aligned diameter curve (thin lines) behind
  the registered **mean ± std** curve, a caption listing the applied alignment
  shifts, and a scalar summary (group mean, between-replicate std, CV,
  replicates used, overlap, registration status — hover the ⓘ icons for exact
  definitions).

**03 Tensile**
- Shown once the group panel produces a mean diameter and tensile data was
  matched in the sidebar's Tensile expander: the stress-strain curve for the
  group (fitted Young's modulus segment, toughness as the area under the
  curve, and the detected break point) plus the derived scalar metrics.

**04 Export & batch**
- **Output folder**: where results are written (default `./fibrecv_output`).
- **Export current group**: writes the standard output tree for the loaded group
  (`overlays/`, `per_image/{csv,plots,diagnostics}/`,
  `per_sample/{csv,plots,shifts}/`) at the current parameters.
- **Run batch (whole folder)**: measures every image in the selected folder
  in-process (with a progress bar and a parallel-jobs selector), aggregates all
  groups, writes the full tree including `summary/master_summary.csv` and
  `summary/run_config.json`, then shows `master_summary` with a CSV download
  button.

---

## 4. Multi-angle cross-section mode

For C1-style image sets named `<condition>_<fibre>_a<angle>_part<part>.tiff`
(six rotation-angle images per fibre part). Switch to it with the **Analysis**
radio at the top of the sidebar; see
[docs/features/06_gui_multiangle.md](docs/features/06_gui_multiangle.md) for
the full design writeup.

**Sidebar — source and scale**
- **Local folder** (reuses the same "Image folder" field as Replicates mode)
  or **Upload** (no whole-folder-upload option here — a real C1 directory is
  tens of GB). Filenames are parsed into condition / fibre / part groups;
  scale-bar name twins and non-matching files are counted and skipped.
  Pick a condition (only shown if the folder has more than one), a fibre and
  a part; the up-to-six matching angle images load.
- **Scale (µm per pixel)**: a manual number field, default `0.388924` (the
  C1 microscope's resolved `Scaling/Items` value from study 03) — **not**
  the sidebar's `ppu` calibration field, which is unused in this mode and
  stays visible but inert. Every µm number in this mode, and the batch's
  scale, comes from this field; changing it does not need an Apply.

**01 Angles**
- One tab per angle (`a1`…`a6`), each with its overlay, per-angle metrics
  (median width in µm and in px, coverage, tilt, flags) and, behind an
  "Edit boundaries" checkbox, the same manual boundary editor as Replicates
  mode. A tab gets a `⚠` prefix when its image is QC-excluded from the fit,
  carries an anomaly flag, or has an uncertain/saturated alignment shift.

**02 Cross-section**
- Six status chips (in fit / no image / QC-excluded / uncertain shift), the
  six aligned width curves in µm, and metrics for the fitted per-position
  ellipse: median area with a split-half uncertainty band, axis ratio,
  orientation, semi-axes, and fit residual. Fitting needs all 3 of the 3
  projection directions covered (`a1`/`a4` = 0°, `a2`/`a5` = 60°,
  `a3`/`a6` = 120°); fewer than that shows an error instead of a result.
  With fewer than all six angles present, the fit can still run but the
  split-half uncertainty needs all six, so the area shows `±—` with an
  explanatory caption. A drawn cross-section, an area-vs-position plot, and
  a per-angle QC table (shift, correlation peak, uncertain/saturated,
  residual) round out the section.

**03 Batch & export**
- One button, **Run multi-angle batch (measure + cross-sections)**: measures
  every image for the selected condition (all fibres/parts, not just the one
  loaded above — "roughly 20 minutes for 450 images"), then fits every
  cross-section, and offers `xsection_summary.csv` for download with a
  `status` column flagging low-confidence fibres. Folder source only —
  uploads have no batch button since there is nothing on disk to re-read.
  The batch always recomputes from disk and **ignores** manual boundary
  edits made in 01 Angles; point it at a fresh output folder, since it
  reuses whatever `per_image/*` is already there.

---

## 5. Notes & troubleshooting

- **Performance**: one image is ~1–2 s to compute (HSV on a 5 MP photo), so a
  3-replicate group is ~3–6 s per **Apply**. Results are cached, so switching
  tabs or re-applying unchanged parameters is instant; the JPEG is decoded only
  once.
- **Batch parallelism**: defaults to 4 worker processes. On Windows the app
  automatically falls back to sequential processing if the process pool can't
  start — the progress bar works either way.
- **`fibrecv-gui: command not found`**: make sure you're in the activated venv
  (pip route) or prefix with `uv run` (uv route), from inside the `fibrecv`
  folder.
- **Outputs are identical to the CLI**: the GUI reuses the exact compute and
  aggregation code, so a group exported from the GUI matches
  `python -m fibrecv.run_measure` + `python -m fibrecv.run_aggregate` for the
  same parameters.
