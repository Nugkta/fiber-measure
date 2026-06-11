# fibrecv GUI — local desktop app (Mac / Windows)

A local [Streamlit](https://streamlit.io) web app for tuning the fibre-diameter
detector, previewing detected boundaries, and batch-processing + exporting — all
from your own machine. It runs **entirely locally**: it reads images from a
folder you point it at and writes outputs to a folder you choose. No cluster
access is needed at runtime.

The app is a thin front-end over the validated `fibrecv` pipeline — it changes
**no** detection logic and uses the same calibrated defaults as the CLI.

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

**Sidebar — Data source**
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

**Sidebar — Parameters**
- Four knobs are exposed; everything else uses the validated defaults (the CLI
  retains full parameter control):
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
- Edits are **staged**: change what you want, then click **Apply** to re-render
  (a few seconds for a group). **Reset to defaults** restores the calibrated
  values.

**Main area**
- One tab per replicate: the full-resolution **overlay** (cyan top edge / yellow
  bottom edge / dashed centerline), the **diameter-vs-position** plot (raw points
  + smoothed line, µm), and scalar metrics (median Ø, coverage, tilt, QC flags).
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
- **Group panel**: every replicate's aligned diameter curve (thin lines) behind
  the registered **mean ± std** curve, a caption listing the applied alignment
  shifts, and a scalar summary (group mean, between-replicate std, CV,
  replicates used, overlap, registration status — hover the ⓘ icons for exact
  definitions).

**Export & batch**
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

## 4. Notes & troubleshooting

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
