# fibrecv GUI redesign — design spec

Date: 2026-08-02
Branch: `feat/gui-redesign`
Status: approved in brainstorming (visual companion session); pending spec review

## Goal

Restyle the Streamlit GUI (`src/fibrecv/gui_app.py`) from the stock Streamlit look
into a clean, scientific-instrument aesthetic — light surfaces, Manchester violet
accent, matched plots, tidy bounded numbers — without changing any behaviour.

## Decisions made (with the user)

- **Aesthetic**: "Clean lab (light)" — near-white surfaces, hairline borders,
  white sidebar, one accent colour. Accent = Manchester deep violet `#660099`.
- **Layout**: one sectioned page (keep the current top-to-bottom flow), numbered
  card sections, small fixed "On this page" jump menu on the right. Rejected:
  top-level tabs (splits the tuning loop), hybrid tabs.
- **Framework**: stay in Streamlit 1.58. Rejected: any rewrite.
- **Implementation**: on `feat/gui-redesign`, executed by subagents.
- **Special emphasis**: numbers are *bounded* — inputs clamped to valid ranges,
  displays at fixed precision, nothing overflows its container.

## Non-goals

- No changes to detection logic, `compute.py`, QC, registration, export formats,
  CLI, or session-state keys/caching.
- No new features; no removal of any existing control or readout.
- No dark mode in this pass (theme structure must not preclude adding one later).

## Design

### 1. Theme foundation

New committed file `.streamlit/config.toml`:

- `base = "light"`, `primaryColor = "#660099"` (violet: buttons, slider, focus,
  active tab), `backgroundColor = "#F8FAFC"`, `secondaryBackgroundColor =
  "#FFFFFF"` (sidebar/cards), `textColor = "#0F172A"`, sans-serif font.
- Radius/border options (`baseRadius`, `borderColor`, …) as supported by 1.58.

One CSS block injected once at the top of `main()` via `st.markdown(...,
unsafe_allow_html=True)`, kept in a module-level constant (or small helper
`_inject_css()`): section cards, header bar, metric cards, jump menu,
sidebar group captions, hide Streamlit chrome (Deploy button, main menu,
footer). All selectors grouped in one place so Streamlit upgrades only require
touching one block.

### 2. Header

Replace `st.title` + caption with a slim header: violet square glyph, app name
**fibrecv**, muted tagline "fibre diameter measurement", and a right-aligned
violet chip showing current state: `group 1_1 · 3 replicates · edge_z 4.0`.
When nothing is loaded the chip reads `no data loaded`.

### 3. Sidebar

Three visually separated groups with small uppercase captions:

- **Data** — source radio, folder path / uploader, group select, count caption.
- **Detection** — the three knobs + calibration + violet full-width Apply,
  quiet Reset button. Plain-language labels with technical names kept:
  "Boundary tightness (edge_z)", "Faint-fibre guard (edge_frac)",
  "Smoothing width (wcol)", "Pixels per micron (ppu)". Existing help text kept
  but tightened.
- **Tensile** — everything tensile (source, gauge length, modulus window)
  inside one expander, collapsed by default.

### 4. Main page sections

Same order as today, each wrapped in a bordered white card with a numbered
header (`01 Replicates`, `02 Group panel`, `03 Tensile`, `04 Export & batch`)
and a muted one-line subtitle. Per-image statistics stay inside the Group card.
Replicate tabs relabelled `Rep 1`, `Rep 2`, … (uploaded/ungrouped files keep
their filename as label). The manual-edit expander stays inside the Replicates
card, retitled "Edit boundaries (manual correction)" with its existing controls
untouched.

Jump menu: fixed-position column on the right with anchor links to the four
section headings (Streamlit headings already emit anchors). Pure CSS; hidden
below ~1200 px viewport width. If anchors prove unreliable it degrades to
nothing — the page works identically without it.

### 5. Metrics

`st.metric` rows become metric cards (CSS on the existing widgets — no layout
rewrite): white card, hairline border, small uppercase label, large
tabular-numeral value. Primary metric per row (mean Ø / group mean Ø) in
violet. Status-like values coloured: green for good (`flags: none`), amber for
warnings (`registration: uncertain`).

### 6. Plots

One shared style helper (e.g. `_plot_style()` context or rcParams dict applied
in `_profile_fig`, `_group_fig`, `_tensile_fig`):

- No top/right spines; light `#E2E8F0` grid; labels/ticks in `#64748B`.
- Primary line / registered mean violet `#660099`; individual replicate
  curves use a fixed muted categorical cycle (`#A78BFA` violet, `#F9A8D4`
  pink, `#93C5FD` blue, `#FCD34D` amber, recycled) so replicates stay
  distinguishable; ±std band violet at low alpha; manual-edit magenta and
  overlay colours (cyan/yellow/green) unchanged — they encode meaning on the
  images.
- Same font family as the app (matplotlib default sans is acceptable — no
  font-file shipping), consistent sizes; figure backgrounds transparent so
  cards show through; `use_container_width=True` everywhere.
- Titles move out of matplotlib into the card headers/captions; axes keep
  their labels (µm, px).

### 7. Bounded numbers (user-flagged as especially important)

- **Inputs**: every numeric widget passes explicit `min_value`/`max_value`.
  Existing bounds kept: `PARAM_SPECS` lo/hi (edge_z 1–12, edge_frac 0–1,
  wcol 1–201, ppu 0.1–10), modulus window 0.02–0.5, parallel jobs 1–16.
  Gap to fix: gauge length currently has `min_value=0.1` but no max — add
  `max_value=1000.0` mm. Audit sweep over all `number_input`/`slider` calls
  (incl. the manual-edit nudge inputs) to confirm both bounds are present.
  Nothing applies out-of-range: clamping happens in the widget, and
  `_cfg_from_items` continues to coerce types.
- **Display precision**: one formatting helper used by all metrics/tables:
  µm → 2 dp; CV → 3 dp; coverage → integer % in metrics, 1 dp in tables;
  px → integer; ppu → 4 dp; tilt slope → 4 dp. Tabular numerals everywhere.
- **Containment**: metric values never overflow their card (CSS `min-width: 0`,
  ellipsis + full value in tooltip for long text like flag lists); long file
  names/paths in tables and captions ellipsized; figures always fit container
  width; sidebar path inputs scroll internally (native input behaviour), never
  widen the sidebar.

### 8. Wording pass

Plain language for all user-facing text; technical identifiers kept in
parentheses where they aid reproducibility. Examples: header caption dropped
("Strictness knob edge_z = 4.0" → state chip); "n reps used" → "reps used";
info/empty states written as sentences ("Set a tensile data folder … to see
stress–strain."). Help text (tooltips) tightened to ≤3 sentences each, keeping
the recommended ranges.

## Implementation notes

- Files touched: `src/fibrecv/gui_app.py` (styling, labels, formatting helper,
  plot style), new `.streamlit/config.toml`, `GUI_README.md` screenshots/refs
  if any, `README.md` one-liner, `docs/features/<NN>_gui-redesign.md` at the end
  (per CLAUDE.md), folder README headers per documentation standards.
- Suggested decomposition for subagents (final plan via writing-plans):
  theme+CSS+header / sidebar+wording / metrics+formatting helper / plots /
  verification (screenshots + smoke test).
- Verification: `uv run pytest tests/test_gui_smoke.py` green; Playwright
  screenshot sweep (same harness used in brainstorming) compared against the
  approved mockup; manual-edit flow exercised once end-to-end.

## Risks / caveats

- Streamlit CSS class names are not a stable API — mitigated by keeping all
  selectors in one block and preferring `data-testid` selectors.
- The jump menu depends on Streamlit heading anchors; acceptable degradation
  is documented above.
- `streamlit-image-coordinates` (manual edit) renders its own component; it
  gets container-width sizing but not deep restyling.
