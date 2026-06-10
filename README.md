# fibrecv

Fibre diameter measurement from optical microscopy images.

`fibrecv` detects the two walls of a fibre lying roughly horizontally in a
micrograph, measures its diameter column by column, applies quality control,
and registers repeated images of the same sample into a mean ± std diameter
profile. It ships as a two-stage command-line pipeline plus a local
[Streamlit](https://streamlit.io) GUI for interactive tuning, preview, and
batch processing.

Detection runs on a per-image *desaturation z-map*: the fibre desaturates the
background colour, so thresholds are self-normalising per image and robust to
illumination changes. Boundaries are placed by per-column edge detection with
wall/shadow discrimination, followed by outlier rejection and smoothing.

## Install

Requires Python 3.12+.

```bash
# with uv (recommended)
uv sync

# or with pip
python -m venv .venv && source .venv/bin/activate
pip install -e .
```

## Input images and naming

Any folder of `.jpg/.jpeg/.png/.tif/.tiff/.bmp` images. Replicates of the
same sample are grouped by the numbers at the end of the filename: the last
number is the replicate, the numbers before it form the group label.

| Filename             | Group  | Replicate |
|----------------------|--------|-----------|
| `masp2 3_1_2.jpg`    | `3_1`  | 2         |
| `3-1-2.png`          | `3_1`  | 2         |
| `sampleA 10_5_1.tif` | `10_5` | 1         |
| `IMG_0123.jpg`       | `IMG`  | 123       |

## GUI

```bash
uv run fibrecv-gui
```

Opens at http://localhost:8501. Point it at an image folder (or drag in any
number of files), pick a group, and tune the three boundary knobs — `edge_z`
(tightness), `edge_frac` (faint-fibre safeguard), `wcol` (anti-jitter
smoothing) — then export results or batch-process the whole folder. See
[GUI_README.md](GUI_README.md) for the full guide.

## CLI pipeline

Stage 1 — measure every image (per-image CSV profiles, overlays, diagnostics):

```bash
uv run python -m fibrecv.run_measure --root "path/to/images" --all
```

Stage 2 — register replicates into per-sample mean ± std curves and a master
summary table:

```bash
uv run python -m fibrecv.run_aggregate --all
```

Both stages write to `./fibrecv_output` by default (`--out` to change).
Select a subset with `--groups 3_1 10_5` or `--glob "*_1_*.jpg"`; parallelise
stage 1 with `--jobs N`. Every detection parameter can be overridden by flag
(`--edge-z`, `--wcol`, ...) — see `--help` and the documented defaults in
`src/fibrecv/config.py`.

## Output tree

```
fibrecv_output/
├── overlays/                 # full-res images with detected boundaries drawn
├── per_image/
│   ├── csv/                  # diameter-vs-position profile per image
│   ├── plots/
│   └── diagnostics/          # per-image meta (coverage, QC flags, ...)
├── per_sample/
│   ├── csv/                  # registered mean ± std curve per group
│   ├── plots/
│   └── shifts/
└── summary/
    ├── master_summary.csv    # one row per group: mean Ø, std, CV, n, QC
    ├── run_config.json       # full parameter + version provenance
    └── run_log.txt
```

Calibration: diameters are converted to microns with `ppu` (pixels per
micron, default 1.3680); pass `--ppu` to match your optics.

## Tests

```bash
uv run pytest
```

## License

[MIT](LICENSE)
