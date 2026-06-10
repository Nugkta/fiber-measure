# GUI Simplification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Unlimited GUI uploads, a 3-knob sidebar (`edge_z`/`edge_frac`/`wcol`), and a robust trailing-numbers filename parser shared by the GUI and both CLIs.

**Architecture:** The naming rule lives in `src/fibrecv/io_utils.py` (`parse_name`, `discover_images`, new `natural_key`) so the GUI, `compute.py`, `run_measure.py` and `run_aggregate.py` all agree. The GUI (`src/fibrecv/gui_app.py`) shrinks its parameter form to three specs and groups both folder images and uploads through the same parser, with an "ungrouped" bucket for unparseable names. `CONFIG` and CLI flags are untouched.

**Tech Stack:** Python 3.12, Streamlit, pytest (new dev dependency), uv.

**Spec:** `docs/superpowers/specs/2026-06-10-gui-simplification-design.md`

**Working directory for all commands:** the project root
`/Users/stan/Library/CloudStorage/OneDrive-TheUniversityofManchester/Y2_onedrive/Projects/spins/fiber-measure`
(quote the path — it contains spaces).

---

### Task 1: Add pytest as a dev dependency

**Files:**
- Modify: `pyproject.toml` (via `uv add`)
- Modify: `uv.lock` (via `uv add`)

- [ ] **Step 1: Add pytest**

```bash
uv add --dev pytest
```

- [ ] **Step 2: Verify pytest runs**

Run: `uv run pytest --version`
Expected: `pytest 8.x.y` (any 8+ version)

- [ ] **Step 3: Commit**

```bash
git add pyproject.toml uv.lock
git commit -m "build: add pytest dev dependency"
```

---

### Task 2: Robust `parse_name` + `natural_key` in io_utils (TDD)

**Files:**
- Create: `tests/test_io_utils.py`
- Modify: `src/fibrecv/io_utils.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_io_utils.py`:

```python
"""Tests for filename parsing and image discovery (io_utils)."""

import pytest

from fibrecv.io_utils import discover_images, natural_key, parse_name


@pytest.mark.parametrize("name,group,rep", [
    ("masp2 10_5_2.jpg", "10_5", 2),           # current convention, unchanged
    ("masp2 3_1_10.jpg", "3_1", 10),           # multi-digit replicate
    ("MASP2 3_1_2.jpg", "3_1", 2),             # prefix case is irrelevant
    ("3-1-2.png", "3_1", 2),                   # dash separators, no prefix
    ("sampleA 10_5_2.tif", "10_5", 2),         # arbitrary text prefix
    ("fiber_3_1 (2).jpg", "3_1", 2),           # parenthesised replicate
    ("3_1_2.jpg", "3_1", 2),                   # bare numbers
    ("IMG_0123.jpg", "IMG", 123),              # single trailing number
    ("scan 7.jpeg", "scan", 7),                # single number, space separator
])
def test_parse_name(name, group, rep):
    assert parse_name(name) == (group, rep)


@pytest.mark.parametrize("name", [
    "background.jpg",        # no digits at all
    "masp2 3_1_2 copy.jpg",  # text after the trailing numbers
    "3.jpg",                 # single number with no prefix -> no group
])
def test_parse_name_rejects(name):
    with pytest.raises(ValueError):
        parse_name(name)


def test_natural_key_orders_numeric_groups():
    groups = ["10_5", "3_3", "3_1", "IMG"]
    assert sorted(groups, key=natural_key) == ["3_1", "3_3", "10_5", "IMG"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_io_utils.py -v`
Expected: FAIL — `ImportError: cannot import name 'natural_key'` (and, once that
exists, parse failures for non-`masp2` names).

- [ ] **Step 3: Implement the parser**

In `src/fibrecv/io_utils.py`, replace the `_NAME_RE` constant and the
`parse_name` and `_sort_key` functions (keep `discover_images` and `load_rgb`
as they are for now) with:

```python
# Separators that may appear between name tokens: space _ - . ( ) [ ]
_SEP_RE = re.compile(r"[\s_\-.()\[\]]+")
_DIGITS_RE = re.compile(r"(\d+)")


def natural_key(text: str):
    """Sort key interleaving numeric and text runs: '3_1' < '3_3' < '10_5'."""
    return tuple((0, int(t)) if t.isdigit() else (1, t.lower())
                 for t in _DIGITS_RE.split(text) if t)


def parse_name(path: str | Path) -> tuple[str, int]:
    """Parse an image name into ``(group, replicate)`` via trailing numbers.

    The stem is split on spaces/underscores/dashes/dots/brackets and the
    trailing run of integer tokens drives the result:

    - two or more trailing integers: the last is the replicate, the rest
      joined with ``_`` form the group ('masp2 10_5_2' -> ('10_5', 2),
      '3-1-2' -> ('3_1', 2), 'sampleA 10_5_2' -> ('10_5', 2));
    - exactly one trailing integer: it is the replicate and the text prefix
      is the group ('IMG_0123' -> ('IMG', 123)).

    Raises ``ValueError`` if the stem does not end in an integer token, or a
    single trailing integer has no text prefix (e.g. '3.jpg').
    """
    stem = Path(path).stem.strip()
    tokens = [t for t in _SEP_RE.split(stem) if t]
    i = len(tokens)
    while i > 0 and tokens[i - 1].isdigit():
        i -= 1
    nums = tokens[i:]
    if not nums:
        raise ValueError(f"unrecognised image name: {path!r}")
    replicate = int(nums[-1])
    if len(nums) >= 2:
        return "_".join(nums[:-1]), replicate
    prefix = " ".join(tokens[:i])
    if not prefix:
        raise ValueError(f"unrecognised image name: {path!r}")
    return prefix, replicate


def _sort_key(path: Path):
    """Sort parseable names by (group, replicate); the rest by filename."""
    try:
        group, rep = parse_name(path)
        return (0, natural_key(group), rep)
    except ValueError:
        return (1, natural_key(Path(path).name), 0)
```

Also update the module docstring's Inputs/Output bullets to describe the new
rule (replace the two `masp2`-specific lines):

```
- A directory of microscopy images whose stems end in numbers (e.g.
  ``masp2 A_B_C.jpg``, ``3-1-2.png``, ``IMG_0123.jpg``) plus sidecar
  ``*.jpg_metadata.xml`` files that must be ignored.
...
- ``parse_name(path)`` -> ``(group, replicate)`` from the trailing run of
  numbers in the stem (last number = replicate, the rest = group).
```

- [ ] **Step 4: Run tests to verify the parse tests pass**

Run: `uv run pytest tests/test_io_utils.py -v`
Expected: PASS (all `test_parse_name*` and `test_natural_key*` cases).

- [ ] **Step 5: Commit**

```bash
git add tests/test_io_utils.py src/fibrecv/io_utils.py
git commit -m "feat: robust trailing-numbers filename parsing in io_utils"
```

---

### Task 3: Extension-tolerant discovery + parse-based group selection (TDD)

**Files:**
- Modify: `tests/test_io_utils.py` (append)
- Modify: `src/fibrecv/io_utils.py:38-46` (`discover_images`)
- Modify: `src/fibrecv/run_measure.py:69-81` (`select_images`) and its io_utils import

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_io_utils.py`:

```python
def test_discover_images_filters_extensions_and_sidecars(tmp_path):
    names = [
        "masp2 3_1_2.jpg", "3-1-2.PNG", "scan 7.tif", "IMG_0123.jpeg",
        "notes.txt", "masp2 3_1_2.jpg_metadata.xml",
    ]
    for n in names:
        (tmp_path / n).write_bytes(b"")
    found = {p.name for p in discover_images(tmp_path)}
    assert found == {"masp2 3_1_2.jpg", "3-1-2.PNG", "scan 7.tif",
                     "IMG_0123.jpeg"}


def test_discover_images_natural_order(tmp_path):
    for n in ["masp2 10_1_1.jpg", "masp2 3_1_10.jpg", "masp2 3_1_2.jpg"]:
        (tmp_path / n).write_bytes(b"")
    assert [p.name for p in discover_images(tmp_path)] == [
        "masp2 3_1_2.jpg", "masp2 3_1_10.jpg", "masp2 10_1_1.jpg"]


def test_select_images_groups_filter(tmp_path):
    import argparse

    from fibrecv.run_measure import select_images

    for n in ["masp2 3_1_1.jpg", "masp2 3_1_2.jpg", "sampleA 10_5_1.png",
              "background.jpg"]:
        (tmp_path / n).write_bytes(b"")
    args = argparse.Namespace(root=str(tmp_path), glob=None, groups=["10_5"])
    assert [p.name for p in select_images(args)] == ["sampleA 10_5_1.png"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_io_utils.py -v`
Expected: the three new tests FAIL (old glob `masp2 *_*.jpg` misses
`3-1-2.PNG`/`scan 7.tif`/`IMG_0123.jpeg`; `select_images` builds
`masp2 10_5_*.jpg` globs so it misses `sampleA 10_5_1.png`).

- [ ] **Step 3: Implement discovery**

In `src/fibrecv/io_utils.py`, add the suffix set near the regexes and replace
`discover_images`:

```python
# Accepted raster formats; the suffix check also excludes the sidecar
# "*.jpg_metadata.xml" files, whose suffix is ".xml".
IMAGE_SUFFIXES = {".jpg", ".jpeg", ".png", ".tif", ".tiff", ".bmp"}


def discover_images(root: str | Path, glob: str = "*") -> list[Path]:
    """Return sorted image files under ``root`` matching ``glob``.

    Keeps only files whose suffix is in ``IMAGE_SUFFIXES`` (case-insensitive).
    """
    root = Path(root)
    paths = [p for p in root.glob(glob)
             if p.is_file() and p.suffix.lower() in IMAGE_SUFFIXES]
    return sorted(paths, key=_sort_key)
```

- [ ] **Step 4: Implement parse-based group selection**

In `src/fibrecv/run_measure.py`, change the io_utils import (line 37) to:

```python
from .io_utils import discover_images, parse_name
```

and replace `select_images` (lines 69-81) with:

```python
def select_images(args: argparse.Namespace) -> list[Path]:
    """Resolve the image selector flags to a sorted, de-duplicated path list."""
    root = Path(args.root)
    if args.glob:
        return discover_images(root, args.glob)
    paths = discover_images(root)
    if args.groups:
        wanted = set(args.groups)
        kept: list[Path] = []
        for p in paths:
            try:
                group, _ = parse_name(p)
            except ValueError:
                continue
            if group in wanted:
                kept.append(p)
        return kept
    # default / --all
    return paths
```

- [ ] **Step 5: Run the full test file**

Run: `uv run pytest tests/test_io_utils.py -v`
Expected: PASS (all tests).

- [ ] **Step 6: Commit**

```bash
git add tests/test_io_utils.py src/fibrecv/io_utils.py src/fibrecv/run_measure.py
git commit -m "feat: extension-tolerant image discovery + parse-based group selection"
```

---

### Task 4: GUI — three-parameter sidebar

**Files:**
- Modify: `src/fibrecv/gui_app.py:69-149` (`PARAM_GROUPS` -> `PARAM_SPECS`, `_INT_FIELDS`)
- Modify: `src/fibrecv/gui_app.py:369-406` (`_param_form`)

- [ ] **Step 1: Replace the parameter metadata**

In `src/fibrecv/gui_app.py`, delete the whole `PARAM_GROUPS` block and the
`_INT_FIELDS` comprehension (lines 66-149, from the `# --- parameter
metadata...` comment through the `_INT_FIELDS` assignment) and replace with:

```python
# --- visible parameters: the three knobs that move the detected boundary.
# Everything else in CONFIG stays at the validated defaults (CLI keeps full
# control). spec = (name, kind, help, step, lo, hi, fmt). ---
PARAM_SPECS: list[tuple] = [
    ("edge_z", "slider",
     "STRICTNESS KNOB: boundary level z-units above the wall-local background. "
     "Higher -> tighter (validated on 3_1: ez4 tracks the true wall; higher gets "
     "dragged inward by internal reflections).", 0.5, 1.0, 12.0, "%.1f"),
    ("edge_frac", "float",
     "Relative cap on the edge level for faint/weak walls (keeps the crossing "
     "on the wall when amplitude A is low).", 0.05, 0.0, 1.0, "%.2f"),
    ("wcol", "int",
     "Column-neighbourhood width averaged before edging (>=15 suppresses "
     "internal iridescence banding).", 1, 1, 201, None),
]

# names of int-typed visible fields (so widgets return int, not float)
_INT_FIELDS = {name for (name, kind, *_rest) in PARAM_SPECS if kind == "int"}
```

Note: `_cfg_from_items` only coerces fields listed in `_INT_FIELDS`; all other
int-typed CONFIG fields now come straight from `DEFAULTS.as_dict()` and are
already `int`, so nothing else needs coercion.

- [ ] **Step 2: Flatten the form and merge Apply with defaults**

Replace `_param_form` with:

```python
def _param_form() -> None:
    """Render the 3-knob parameter form; updates session_state on Apply/Reset."""
    applied = st.session_state.cfg_dict
    ver = st.session_state.form_version

    with st.sidebar.form("params", clear_on_submit=False):
        st.markdown("**Parameters** — edit, then click **Apply** to re-render.")
        new_vals: dict = {}
        for (name, kind, help_txt, step, lo, hi, fmt) in PARAM_SPECS:
            key = f"p_{name}_v{ver}"
            cur = applied[name]
            if kind == "slider":
                new_vals[name] = st.slider(
                    name, min_value=float(lo), max_value=float(hi),
                    value=float(cur), step=float(step), help=help_txt, key=key)
            elif kind == "int":
                new_vals[name] = int(st.number_input(
                    name, min_value=int(lo), max_value=int(hi),
                    value=int(cur), step=int(step), help=help_txt, key=key))
            else:  # float
                kwargs = dict(min_value=float(lo), max_value=float(hi),
                              value=float(cur), step=float(step), help=help_txt, key=key)
                if fmt:
                    kwargs["format"] = fmt
                new_vals[name] = float(st.number_input(name, **kwargs))
        c1, c2 = st.columns(2)
        apply = c1.form_submit_button("Apply", type="primary", width="stretch")
        reset = c2.form_submit_button("Reset to defaults", width="stretch")

    if reset:
        st.session_state.cfg_dict = DEFAULTS.as_dict()
        st.session_state.form_version += 1
        st.rerun()
    if apply:
        # merge: the three knobs override a full defaults dict, so hidden
        # fields always carry the validated values
        st.session_state.cfg_dict = {**DEFAULTS.as_dict(), **new_vals}
        st.rerun()  # re-run top-to-bottom so reps recompute with the new params
```

- [ ] **Step 3: Verify the module still imports and tests pass**

Run: `uv run python -c "import fibrecv.gui_app" && uv run pytest tests/ -q`
Expected: no import error; all tests PASS.

- [ ] **Step 4: Commit**

```bash
git add src/fibrecv/gui_app.py
git commit -m "feat: reduce GUI sidebar to edge_z / edge_frac / wcol"
```

---

### Task 5: GUI — unlimited uploads + grouping via the shared parser

**Files:**
- Modify: `src/fibrecv/gui_app.py:57` (io_utils import)
- Modify: `src/fibrecv/gui_app.py:215-221` (`_group_sort_key`)
- Modify: `src/fibrecv/gui_app.py:236-240` (`export_group` error message)
- Modify: `src/fibrecv/gui_app.py:412-473` (`_load_reps`)
- Modify: `src/fibrecv/gui_app.py:579-582` (batch warning text)
- Modify: `src/fibrecv/gui_app.py:629-630` (empty-state info text)
- Modify: `src/fibrecv/gui_app.py:11-14` (module docstring Inputs)

- [ ] **Step 1: Import `natural_key` and rewrite the group sort key**

Change the io_utils import (line 57) to:

```python
from fibrecv.io_utils import discover_images, natural_key, parse_name  # noqa: E402
```

Replace `_group_sort_key` (and add the `UNGROUPED` constant just above it):

```python
UNGROUPED = "ungrouped"


def _group_sort_key(group: str):
    """Natural sort for group labels ('3_1' < '3_3' < '10_5'); ungrouped last."""
    return (group == UNGROUPED, natural_key(group))
```

- [ ] **Step 2: Add the shared grouping helpers**

Insert directly below `_group_sort_key`:

```python
def _group_by_name(items: list, key) -> dict[str, list]:
    """Bucket items by their parse_name group; unparseable names -> UNGROUPED."""
    groups: dict[str, list] = {}
    for it in items:
        try:
            g, _ = parse_name(key(it))
        except ValueError:
            g = UNGROUPED
        groups.setdefault(g, []).append(it)
    return groups


def _sorted_reps(items: list, key) -> list:
    """Sort one group's items by replicate number; unparseable last, by name."""
    def _k(it):
        try:
            return (0, parse_name(key(it))[1], key(it))
        except ValueError:
            return (1, 0, key(it))
    return sorted(items, key=_k)
```

- [ ] **Step 3: Rewrite `_load_reps`**

Replace the whole function with:

```python
def _load_reps(cfg_items: tuple) -> tuple[list[dict], str | None, str | None]:
    """Resolve the data-source controls to a list of replicate dicts.

    Returns ``(reps, group_label, folder)``. Each rep is
    ``{"name", "rgb", "mr", "idx"}``. ``folder`` is the scanned folder path
    (for batch) or None for the upload source. Both sources group images with
    the shared ``parse_name`` rule; unparseable names land in an "ungrouped"
    bucket instead of being hidden.
    """
    st.sidebar.markdown("### Data source")
    source = st.sidebar.radio("source", ["Local folder", "Upload"],
                              horizontal=True, label_visibility="collapsed")
    reps: list[dict] = []
    group_label: str | None = None
    folder: str | None = None

    if source == "Local folder":
        folder = st.sidebar.text_input("Image folder", value=st.session_state.get(
            "folder", DEFAULT_ROOT))
        st.session_state.folder = folder
        if not folder or not Path(folder).is_dir():
            st.sidebar.warning("Enter a valid local folder path.")
            return reps, None, None
        paths = discover_images(folder)
        if not paths:
            st.sidebar.warning("No image files found in this folder.")
            return reps, None, folder
        groups = _group_by_name(paths, key=lambda p: p.name)
        keys = sorted(groups, key=_group_sort_key)
        group_label = st.sidebar.selectbox("Group", keys)
        st.sidebar.caption(f"{len(paths)} images, {len(keys)} groups")
        for i, p in enumerate(_sorted_reps(groups[group_label],
                                           key=lambda p: p.name)):
            mtime = Path(p).stat().st_mtime
            mr = _cached_compute_path(str(p), mtime, cfg_items)
            rgb = _cached_rgb_from_path(str(p), mtime)
            reps.append({"name": Path(p).stem, "rgb": rgb, "mr": mr, "idx": i})
    else:
        uploads = st.sidebar.file_uploader(
            "Upload images", type=["jpg", "jpeg", "png", "tif", "tiff", "bmp"],
            accept_multiple_files=True)
        if uploads:
            groups = _group_by_name(list(uploads), key=lambda u: u.name)
            keys = sorted(groups, key=_group_sort_key)
            group_label = (st.sidebar.selectbox("Group", keys)
                           if len(keys) > 1 else keys[0])
            st.sidebar.caption(f"{len(uploads)} files, {len(keys)} groups")
            for i, up in enumerate(_sorted_reps(groups[group_label],
                                                key=lambda u: u.name)):
                data = up.getvalue()
                stem = Path(up.name).stem
                mr, rgb = _cached_compute_upload(stem, data, cfg_items)
                reps.append({"name": stem, "rgb": rgb, "mr": mr, "idx": i})
            if group_label == UNGROUPED:
                group_label = None  # header falls back to "Replicates (uploaded)"

    return reps, group_label, folder
```

- [ ] **Step 4: Update the remaining masp2-specific texts**

In `export_group`, replace the `ValueError` message:

```python
        raise ValueError(
            "Cannot export: image names must end in numbers "
            "(e.g. 'name 3_1_2.jpg') to derive a group label."
        )
```

In `_render_export_batch`, replace the empty-folder warning:

```python
                st.warning("No image files found in the folder.")
```

In `main()`, replace the empty-state info:

```python
        st.info("Pick a folder + group, or upload images, to begin.")
```

In the module docstring (lines 12-14), replace the Inputs bullet pair:

```
- Images from a local folder OR any number of uploaded files; both are
  auto-grouped via ``parse_name``'s trailing-numbers rule, with unparseable
  names collected in an "ungrouped" bucket.
```

- [ ] **Step 5: Verify import and tests**

Run: `uv run python -c "import fibrecv.gui_app" && uv run pytest tests/ -q`
Expected: no import error; all tests PASS.

- [ ] **Step 6: Commit**

```bash
git add src/fibrecv/gui_app.py
git commit -m "feat: unlimited uploads with shared trailing-numbers grouping in GUI"
```

---

### Task 6: Docs + final verification

**Files:**
- Modify: `GUI_README.md:67-82`

- [ ] **Step 1: Update GUI_README.md section 3**

Replace the "Sidebar — Data source" and "Sidebar — Parameters" blocks
(lines 67-82) with:

```markdown
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
- Three knobs control the detected boundary; everything else uses the
  validated defaults (the CLI retains full parameter control):
  - **`edge_z`** (slider) — the main strictness knob: higher → tighter
    boundaries.
  - **`edge_frac`** — keeps the boundary on faint/weak walls.
  - **`wcol`** — column averaging; raise it if internal iridescence banding
    disturbs the boundary.
- Edits are **staged**: change what you want, then click **Apply** to re-render
  (a few seconds for a group). **Reset to defaults** restores the calibrated
  values.
```

- [ ] **Step 2: Run the full test suite**

Run: `uv run pytest -q`
Expected: all tests PASS.

- [ ] **Step 3: Manual GUI smoke test**

```bash
uv run streamlit run src/fibrecv/gui_app.py --server.headless true
```

Check in the browser (http://localhost:8501):
1. Sidebar shows exactly edge_z / edge_frac / wcol + Apply / Reset.
2. Local folder mode lists groups for an images folder; an oddly named file
   shows under "ungrouped" rather than disappearing.
3. Upload mode accepts more than 3 files; multiple groups produce a Group
   dropdown; the selected group renders tabs + group panel.
4. Apply with a changed edge_z re-renders the overlays.
5. Export current group still writes the output tree.

Stop with Ctrl+C.

- [ ] **Step 4: Commit**

```bash
git add GUI_README.md
git commit -m "docs: GUI README for 3-knob sidebar, unlimited uploads, generic naming"
```
