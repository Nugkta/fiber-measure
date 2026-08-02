# GUI visual redesign
Added: 2026-08-02
Code: `src/fibrecv/gui_app.py` (`_CSS`, `_inject_css`, `_render_header`,
`_render_jump_menu`, `_fmt`, `_styled_fig`, `_ACCENT`/`_REP_CYCLE`, `main`'s
card layout), `.streamlit/config.toml` (Streamlit-engine theme); tests
untouched (`tests/`, 69 passing).

## What it does

Restyles the Streamlit app's look without touching any detection logic,
session-state keys or export/aggregate behaviour. The plain default-Streamlit
page (bare `st.title`, unstyled widgets stacked top to bottom, default
matplotlib figures) becomes a "clean lab" light theme: a slim violet-accent
header with a one-line state chip (group/replicate count/edge_z, or "no data
loaded"), the sidebar regrouped into labelled DATA / DETECTION sections plus a
collapsible Tensile expander with plain-language knob labels, and the main
area split into four numbered card sections — 01 Replicates, 02 Group panel,
03 Tensile, 04 Export & batch — with a fixed jump menu on the right that links
to each (hidden on narrow windows). Metric numbers, plot colours and
NaN/missing-value formatting are now consistent across every card. The
detection pipeline, `CONFIG`, manual-edit flow and output tree are unchanged;
this is a presentation layer over the same data path.

## Design choices

- **"Clean lab" light theme, violet #660099 accent** — chosen to read as a
  measurement tool rather than a generic dashboard; a dark theme and a
  tab-based (rather than single-scroll) layout were both considered and
  rejected — dark mode fights with the bright micrograph overlays it has to
  display, and tabs would hide the group/tensile/export state a user often
  wants to keep visible while tuning replicates.
- **One sectioned page, not multiple pages/tabs.** All four sections
  (replicates/group/tensile/export) stay on one scrollable page inside
  `st.container(key="card_...")` blocks, with a fixed jump menu for
  navigation instead of Streamlit's multipage nav — keeps the existing
  single-`main()` control flow and session-state model untouched (Task 1's
  `_render_group` return-value re-plumb was the only structural change,
  gated by the smoke test's full-flow assertions).
- **CSS is layered onto real Streamlit widgets, not custom HTML replacements.**
  Metric "cards" are ordinary `st.metric` calls wrapped in a keyed
  `st.container`; `_CSS` styles them via `[class*="st-key-metrics_"]
  [data-testid="stMetric"]` selectors. This was required, not just preferred:
  the existing AppTest-based test suite asserts against `at.metric[...]`
  objects, so replacing them with styled `st.markdown` HTML would have broken
  69 passing tests for a purely cosmetic change.
- **A single formatting policy (`_fmt`) and a single plot style
  (`_styled_fig`/`_ACCENT`/`_REP_CYCLE`).** Every metric, caption and the
  header chip previously formatted numbers ad hoc (inconsistent decimal
  places, "nan" leaking into the UI); centralising it in one function/table
  makes "missing value" render identically everywhere (em dash) and every
  matplotlib figure share spine/grid/colour treatment, instead of restyling
  each plot builder separately.
- **Bounded inputs, not free-form numeric fields.** Gauge length is capped at
  a `max_value=1000.0` mm sanity ceiling and manual-edit nudges are clamped to
  ±200 px; both were previously unbounded `st.number_input`s that could accept
  values that silently break downstream geometry (e.g. a nudge larger than the
  image).

## Algorithm details

- **Theme split across two layers.** `.streamlit/config.toml` sets what
  Streamlit's own theming engine controls: `[theme]` base colours
  (`primaryColor = "#660099"`, background/text/border colours), radii
  (`baseRadius`, `buttonRadius`), metric font size/weight, and
  `[client] toolbarMode = "minimal"` to hide the default chrome.
  `gui_app.py`'s `_CSS` constant (module-level, injected once by
  `_inject_css()` at the top of `main()`) covers what the theme engine
  cannot reach: the custom header/chip markup, section-card framing, the
  fixed jump nav, and a chrome backstop (`#MainMenu`/`footer`/deploy-button
  `display: none` in case `toolbarMode` doesn't hide them in a given
  Streamlit version).
- **Hooks, not brittle selectors.** Every custom rule targets either a
  `data-testid="st..."` attribute (Streamlit's own semantic markup) or a
  `.st-key-<name>` class, which Streamlit derives verbatim from
  `st.container(key="...")`. Card sections use `card_replicates` /
  `card_group` / `card_tensile` / `card_export`; per-replicate metric rows use
  `metrics_rep_<safe_name>` (sanitised via `_safe_key`) and the group panel
  uses `metrics_group`; a flags metric additionally gets `status_ok_...` /
  `status_warn_...` for green/amber colouring. Prefix matching
  (`[class*="st-key-metrics_"]`) lets one rule cover every replicate's key
  without enumerating them.
- **`_fmt(value, kind, dash="—")`** looks `kind` up in a module-level
  `_FMT_KINDS` table of `(format_string, unit_suffix)` pairs (e.g.
  `"um" -> ("{:.2f}", " µm")`, `"pct100" -> ("{:.2f}", " %")`,
  `"int" -> ("{:.0f}", "")`); `None` or non-finite input returns `dash`
  unconditionally, so a NaN never reaches the page as the literal string
  `"nan"`.
- **`_styled_fig(figsize=(9, 3))`** builds a `Figure`/`Axes` pair with
  transparent figure/axes backgrounds (so the plot sits on the card instead
  of a white rectangle), top/right spines dropped, remaining spines/ticks/
  labels in muted slate (`_MUTED = "#64748B"`), and a faint slate grid
  (`_GRID = "#E2E8F0"`). `_ACCENT = "#660099"` is the primary series colour
  (e.g. the stress-strain curve); `_REP_CYCLE` is a 4-colour cycle for
  per-replicate overlays in the group panel.
- **Jump menu** (`_render_jump_menu`) is a `position: fixed` nav
  (`.fcv-jump`, `top: 5.5rem; right: 1rem`) with anchor links
  (`#replicates`, `#group-panel`, `#tensile`, `#export`) matching each card's
  `st.subheader(..., anchor=...)`; a `@media (max-width: 1200px)` query hides
  it on narrow windows rather than trying to reflow it. It is only rendered
  once data is loaded (`main()` returns early on the "pick a folder" state
  before calling it), and card 03 (Tensile) is conditionally skipped when the
  group panel produces no mean, so its jump link is simply inert on that page
  — the nav does not special-case a missing target.
- **Folder-upload remount fix** (`_enable_folder_upload`, not new this task
  but part of the same UI layer): sets `webkitdirectory`/`directory` on the
  hidden `<input type=file>` of the matching `st.file_uploader` via a
  zero-height `components.html` script reaching `window.parent.document`. A
  `MutationObserver` on `doc.body` re-applies the attribute whenever a
  matching uploader (re)appears in the DOM — needed because a collapsed
  sidebar expander re-mounts the uploader and drops the one-shot attribute.
  The observer is created once per label (guarded by a flag on
  `window.parent`) and disconnected on `unload`, clearing the flag so a later
  re-creation (e.g. toggling the folder-mode checkbox off then on) gets a
  fresh observer instead of being silently blocked by a stale one.

## Caveats

- **Streamlit CSS class churn.** `.st-key-*` classes and `data-testid`
  attributes are Streamlit-internal, undocumented-as-a-contract naming; a
  future Streamlit upgrade could rename or restructure them, silently
  un-styling (not breaking) the app. All such selectors are kept in the one
  `_CSS` constant specifically so a fix is a single-location diff.
- **Jump-menu anchor dependence.** The nav's `href="#..."` targets rely on
  `st.subheader(..., anchor=...)` continuing to render a matching HTML id;
  if a card's subheader text or anchor argument changes without updating
  `_render_jump_menu`'s hrefs, that link silently does nothing (no error,
  just no scroll).
- **`streamlit-image-coordinates` component not restyled.** The manual-edit
  zoomed-strip click widget (in `_render_edit_expander`) is a third-party
  component rendered in its own iframe; `_CSS` cannot reach into it, so its
  chrome is unchanged from the component's own default look.
- **Folder-upload `MutationObserver` cleanup is unload-based.** The observer
  disconnects on the *component iframe's* `unload` event, not on a Streamlit
  rerun per se — this is believed correct for both expander re-mounts and
  full reruns (both remount the iframe), but was not exhaustively tested
  against every Streamlit rerun path; a regression here would surface as the
  folder-picker silently reverting to single-file mode after some sequence of
  sidebar interactions.
- **Fixed breakpoint, not responsive redesign.** The jump menu's 1200px
  cutoff and the card layout are tuned for a normal desktop browser window;
  very narrow windows lose the jump menu (by design) but the cards themselves
  do not reflow into a mobile-friendly layout.
