# GUI simplification: unlimited uploads, 3-knob sidebar, robust naming

Date: 2026-06-10
Status: approved

## Goal

Make the local Streamlit GUI (`src/fibrecv/gui_app.py`) simpler to use:

1. The upload source accepts any number of images (currently capped at 3).
2. The sidebar exposes exactly three tuning parameters (currently ~28).
3. Filename group detection works for naming formats beyond `masp2 A_B_C.jpg`,
   in the GUI and everywhere else names are parsed.

The detection pipeline itself is untouched; `CONFIG` keeps every field and both
CLIs keep full parameter control.

## 1. Unlimited uploads

- `st.file_uploader` loses its "1-3" label and the `uploads[:3]` truncation in
  `_load_reps`; any number of files is accepted.
- The folder-fallback multiselect ("Pick up to 3 images") loses
  `max_selections=3`.
- Uploaded files are grouped with the same `parse_name` rule as folder mode.
  If the uploads parse into more than one group, a group selectbox appears
  (same UX as folder mode) and the main area renders the selected group.
  Files whose names do not parse fall into a single "ungrouped" bucket that is
  treated as one ad-hoc group, exactly like today's unnamed uploads.

## 2. Three-parameter sidebar

Visible parameters, with widget type and rationale:

| Param       | Widget | Why it stays                                              |
|-------------|--------|-----------------------------------------------------------|
| `edge_z`    | slider | Primary strictness knob: boundary level above background. |
| `edge_frac` | number | Keeps the boundary on faint/weak walls (relative cap).    |
| `wcol`      | number | Column averaging; suppresses internal iridescence banding.|

- `PARAM_GROUPS` shrinks to these three specs; the form renders them flat
  (no expanders) with the existing Apply / Reset buttons.
- All other `CONFIG` fields are fixed at the validated defaults and are not
  shown or editable in the GUI. Export/batch provenance (`run_config.json`)
  still records the full effective config.
- Reset restores the three knobs to `CONFIG()` defaults.

## 3. Robust name parsing (`io_utils.py`)

`parse_name` is rewritten around a trailing-numbers rule so every consumer
(GUI grouping, `compute.py` meta, `run_aggregate` bucketing, `run_measure`
selection, `export_group`) benefits from one shared parser:

- Strip the extension; find the trailing run of integers at the end of the
  stem, where integers may be separated by spaces, `_`, `-`, `.`, or wrapped
  in parentheses/brackets.
- The last integer is the replicate. The integers before it, joined with `_`,
  form the group.
- If only one trailing integer exists, the normalized text prefix becomes the
  group (e.g. `IMG_0123.jpg` -> group `IMG`, replicate 123).
- No digits at the end -> `ValueError`, as today. The GUI catches this and
  shows the file under an "ungrouped" bucket instead of hiding it; the CLIs
  keep skipping such files.

Compatibility examples:

| Filename            | Group  | Replicate |
|---------------------|--------|-----------|
| `masp2 10_5_2.jpg`  | `10_5` | 2 (unchanged from today) |
| `3-1-2.png`         | `3_1`  | 2 |
| `sampleA 10_5_2.tif`| `10_5` | 2 |
| `fiber_3_1 (2).jpg` | `3_1`  | 2 |
| `IMG_0123.jpg`      | `IMG`  | 123 |

Known trade-off: `sampleA 1_2_3` and `sampleB 1_2_3` in the same folder merge
into group `1_2`. Accepted; the group selector makes it visible.

### Discovery

- `discover_images` drops the hardcoded `"masp2 *_*.jpg"` glob: default glob
  becomes `*`, filtered to suffixes `.jpg/.jpeg/.png/.tif/.tiff/.bmp`
  (case-insensitive). Sidecar `*.jpg_metadata.xml` files remain excluded by
  the suffix check.
- `run_measure --groups` stops building per-group globs and instead discovers
  all images, then filters by `parse_name(p)[0] == group`.
- GUI folder mode no longer has a "no masp2 names found" fallback path; all
  discovered images appear grouped, with unparseable ones under "ungrouped".

## Error handling

- Upload of a corrupt/undecodable image: existing behavior (exception surfaces
  per file) is unchanged.
- `export_group` still refuses groups without a parseable label, but the
  robust parser makes that rare; the error message no longer references the
  `masp2` convention.

## Testing

- Unit tests for `parse_name` covering the table above plus: no-digit names
  (ValueError), multi-digit replicates, mixed separators, and the
  `masp2 10_5_2` backward-compatibility case.
- Unit test for `discover_images` extension filtering and metadata-sidecar
  exclusion.
- Manual GUI smoke run against the existing image folder: groups list,
  >3 uploads, 3-knob form, export, and batch.
