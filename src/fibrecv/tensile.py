"""Tensile (stress-strain) ingestion and single-fibre metric computation.

Dependencies
------------
``numpy``, ``pandas`` (CSV/Excel read + rolling median), ``re`` / ``pathlib`` /
``dataclasses`` / ``math``. Reuses ``CONFIG`` (``config``) and ``natural_key``
(``io_utils``). No Streamlit, so it is headlessly testable.

Inputs
------
- A folder of TA Instruments WinTest CSV exports named like
  ``masp2 10_1 06102026 122956_std.CSV`` (the ``_std`` full test; ``_rdr`` quick
  tests are only used when no ``_std`` exists for that fibre). Each file carries
  ~38 metadata lines then a data block headed
  ``"Points","Elapsed Time","Disp","Load 3",`` (the ``Points`` counter resets
  every scan, so it is never used as an index).
- The matched image group's measured **mean** diameter (microns) for the
  cross-sectional area, and a **gauge length** (mm) for strain.

Output
------
- ``parse_tensile_name`` / ``discover_tensile`` -> a ``{group: path}`` mapping
  whose keys match ``io_utils.parse_name``'s group strings (e.g. ``10_1``), so a
  curve pairs with exactly one image group.
- ``read_trace`` -> tidy ``DataFrame[disp_mm, load_n]``.
- ``compute_tensile`` -> ``TensileResult``: breaking force, tensile strength,
  extension/strain at break, Young's modulus (steepest well-fit initial slope)
  and toughness (area under the stress-strain curve).
- ``tensile_row`` / ``build_matrix`` -> one-row-per-fibre metrics for export.

Metrics follow the user's reference table: A = pi*(d/2)^2; strain = dL/L0;
strength = Fmax/A; modulus = slope of the initial linear region; toughness =
area under the stress-strain curve.

Pos
---
Parallel to the diameter pipeline. The GUI feeds it ``register_sample``'s group
mean diameter live; the batch export feeds it ``master_summary.csv``'s
``mean_um`` column.
"""

from __future__ import annotations

import io
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from fibrecv.config import CONFIG
from fibrecv.io_utils import natural_key

# A tensile filename's group token: two-or-more underscore-joined integers, e.g.
# '10_1' -- the exact string parse_name yields for 'masp2 10_1_3.jpg'.
_GROUP_RE = re.compile(r"^\d+(?:_\d+)+$")
_KIND_RE = re.compile(r"_(std|rdr)$", re.IGNORECASE)
TENSILE_SUFFIXES = {".csv", ".xls", ".xlsx"}

# candidate header names (case-insensitive), most-specific first
_DISP_NAMES = ("Disp", "Displacement", "Extension")
_LOAD_NAMES = ("Load 3", "Load 2", "Load 1", "Load")


# --------------------------------------------------------------------------- #
# Filename parsing / discovery                                                 #
# --------------------------------------------------------------------------- #
def parse_tensile_name(path: str | Path) -> tuple[str, str]:
    """Parse a tensile filename into ``(group, kind)``.

    The group is the ``batch_fibre`` token (e.g. ``10_1``); ``kind`` is ``std``
    or ``rdr`` from the trailing ``_std``/``_rdr`` (default ``std``). Raises
    ``ValueError`` if no group token is present.
    """
    stem = Path(path).stem.strip()
    kind_m = _KIND_RE.search(stem)
    kind = kind_m.group(1).lower() if kind_m else "std"
    for token in stem.split():
        if _GROUP_RE.match(token):
            return token, kind
    raise ValueError(f"unrecognised tensile name: {path!r}")


def _select_by_group(triples) -> dict[str, object]:
    """From ``(group, kind, value)`` triples keep one value per group, ``_std``
    winning over ``_rdr`` (an ``_rdr``-only fibre is still kept)."""
    chosen: dict[str, object] = {}
    kinds: dict[str, str] = {}
    for group, kind, value in triples:
        if group not in chosen or (kinds[group] != "std" and kind == "std"):
            chosen[group] = value
            kinds[group] = kind
    return chosen


def discover_tensile(folder: str | Path) -> dict[str, Path]:
    """Map ``{group: path}`` for tensile files in ``folder``, preferring ``_std``.

    Files that do not parse (no ``batch_fibre`` token) are skipped. When a fibre
    has both a ``_std`` and a ``_rdr`` export the ``_std`` wins; an ``_rdr``-only
    fibre is still included.
    """
    folder = Path(folder)
    triples = []
    for p in sorted(folder.iterdir(), key=lambda q: q.name):
        if not p.is_file() or p.suffix.lower() not in TENSILE_SUFFIXES:
            continue
        try:
            group, kind = parse_tensile_name(p)
        except ValueError:
            continue
        triples.append((group, kind, p))
    return _select_by_group(triples)


def discover_tensile_files(files, name_of=None) -> dict[str, object]:
    """Like ``discover_tensile`` but over an arbitrary iterable of file objects
    (e.g. Streamlit uploads) rather than a folder.

    ``name_of`` extracts a filename from each item (default: its ``.name``).
    Items whose suffix is not a tensile type or whose name does not parse are
    skipped; ``_std`` wins over ``_rdr`` per fibre. Returns ``{group: file}``;
    the file objects pair with the polymorphic ``read_trace``.
    """
    if name_of is None:
        def name_of(f):
            return f.name
    triples = []
    for f in sorted(files, key=name_of):
        nm = name_of(f)
        if Path(nm).suffix.lower() not in TENSILE_SUFFIXES:
            continue
        try:
            group, kind = parse_tensile_name(nm)
        except ValueError:
            continue
        triples.append((group, kind, f))
    return _select_by_group(triples)


# --------------------------------------------------------------------------- #
# Trace reading                                                                #
# --------------------------------------------------------------------------- #
def read_trace(source, name: str | None = None) -> pd.DataFrame:
    """Read a tensile export into a tidy ``DataFrame[disp_mm, load_n]``.

    ``source`` may be a path (str/``Path``), raw ``bytes``, or a file-like upload
    (anything with ``getvalue()``/``read()`` and usually a ``.name``), so the GUI
    can feed both a local folder and drag-and-drop uploads. ``name`` overrides the
    filename used to pick the CSV-vs-Excel reader (required for raw bytes).
    """
    data, nm = _as_bytes_name(source, name)
    if nm.lower().endswith((".xls", ".xlsx")):
        df = _read_excel_bytes(data, nm)
    else:
        df = _read_csv_bytes(data, nm)
    return _tidy_trace(df, nm)


def _as_bytes_name(source, name: str | None) -> tuple[bytes, str]:
    """Normalise a path / bytes / file-like upload to ``(bytes, filename)``."""
    if isinstance(source, (str, Path)):
        p = Path(source)
        return p.read_bytes(), (name or p.name)
    if isinstance(source, (bytes, bytearray)):
        if not name:
            raise ValueError("read_trace: 'name' is required for bytes input")
        return bytes(source), name
    # file-like (e.g. Streamlit UploadedFile): getvalue() does not consume it
    if hasattr(source, "getvalue"):
        data = source.getvalue()
    elif hasattr(source, "read"):
        data = source.read()
    else:
        raise TypeError(f"read_trace: unsupported source {type(source)!r}")
    nm = name or getattr(source, "name", None)
    if not nm:
        raise ValueError("read_trace: 'name' is required for this input")
    data = data if isinstance(data, (bytes, bytearray)) else str(data).encode("latin-1")
    return bytes(data), nm


def _read_csv_bytes(data: bytes, name: str) -> pd.DataFrame:
    """Skip the WinTest metadata header and read the data block from bytes."""
    text = data.decode("latin-1")
    lines = text.splitlines()
    hdr = next((i for i, line in enumerate(lines)
                if line.split(",", 1)[0].strip().strip('"') == "Points"), None)
    if hdr is None:
        raise ValueError(f"no data header ('Points' row) found in {name!r}")
    df = pd.read_csv(io.StringIO(text), skiprows=hdr, engine="python")
    return df.iloc[1:]  # drop the units row ("","Sec","mm","N")


def _read_excel_bytes(data: bytes, name: str) -> pd.DataFrame:
    """Optional .xls/.xlsx support; only imports openpyxl/xlrd if actually used."""
    try:
        raw = pd.read_excel(io.BytesIO(data), header=None)
    except ImportError as exc:  # pragma: no cover - depends on optional extra
        raise ImportError(
            "Reading .xls/.xlsx tensile files needs 'openpyxl' (xlsx) or 'xlrd' "
            "(xls); install it, or export the data as .csv."
        ) from exc
    hdr = next((i for i in range(len(raw))
                if str(raw.iloc[i, 0]).strip().strip('"') == "Points"), None)
    if hdr is None:
        raise ValueError(f"no data header ('Points' row) found in {name!r}")
    df = raw.iloc[hdr + 1:].copy()
    df.columns = [str(c).strip().strip('"') for c in raw.iloc[hdr]]
    return df.iloc[1:]


def _pick_column(columns, candidates) -> str | None:
    """Case-insensitive column match: exact candidate first, then a prefix match."""
    low = {str(c).strip().lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    prefix = candidates[0].split()[0].lower()
    for c in columns:
        if str(c).strip().lower().startswith(prefix):
            return c
    return None


def _tidy_trace(df: pd.DataFrame, source: Path) -> pd.DataFrame:
    """Coerce to numeric, pick the Disp/Load columns, drop NaNs."""
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    disp_col = _pick_column(df.columns, _DISP_NAMES)
    load_col = _pick_column(df.columns, _LOAD_NAMES)
    if disp_col is None or load_col is None:
        raise ValueError(
            f"could not find displacement/load columns in {source!r}; "
            f"saw {list(df.columns)}")
    out = pd.DataFrame({
        "disp_mm": df[disp_col].to_numpy(float),
        "load_n": df[load_col].to_numpy(float),
    })
    return out.dropna().reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Metric computation                                                           #
# --------------------------------------------------------------------------- #
@dataclass
class TensileResult:
    """Per-fibre stress-strain analysis (SI internally; converted at export)."""

    group: str | None
    source: Path | None
    disp_mm: np.ndarray
    load_n: np.ndarray
    strain: np.ndarray          # dimensionless (disp_mm / gauge_length_mm)
    stress_pa: np.ndarray       # NaN-filled when no diameter is available
    diameter_um: float | None
    area_m2: float
    gauge_length_mm: float
    fmax_n: float               # breaking force = peak load up to fracture ("the bump")
    tensile_strength_pa: float  # Fmax / area
    extension_at_break_mm: float
    strain_at_break: float      # fraction (multiply by 100 for %)
    youngs_modulus_pa: float    # steepest well-fit initial slope of stress-strain
    toughness_j_m3: float       # area under stress-strain up to fracture
    modulus_fit: dict           # slope/intercept/r2/strain_lo/strain_hi for plotting
    break_index: int
    flags: list[str] = field(default_factory=list)


def compute_tensile(df: pd.DataFrame, diameter_um: float | None,
                    gauge_length_mm: float, cfg: CONFIG | None = None,
                    break_index: int | None = None) -> TensileResult:
    """Compute tensile metrics from a ``DataFrame[disp_mm, load_n]``.

    ``diameter_um`` is the matched image group's mean diameter; pass ``None`` to
    get force-only metrics (Fmax, extension/strain at break) with the
    stress-based metrics left NaN. ``break_index`` forces the fracture to a
    user-chosen sample (the GUI's manual break point), overriding the automatic
    collapse detector; everything downstream (Fmax, strength, strain at break,
    toughness) is then taken up to that sample.
    """
    cfg = cfg or CONFIG()
    disp = np.asarray(df["disp_mm"], dtype=float)
    load = np.asarray(df["load_n"], dtype=float)
    flags: list[str] = []

    L0 = float(gauge_length_mm)
    strain = disp / L0 if L0 > 0 else np.full_like(disp, np.nan)

    has_d = (diameter_um is not None and np.isfinite(diameter_um)
             and float(diameter_um) > 0)
    area = math.pi * (float(diameter_um) * 1e-6 / 2.0) ** 2 if has_d else math.nan
    stress = load / area if has_d else np.full_like(load, np.nan)

    if disp.size < 3 or not np.isfinite(load).any():
        flags.append("too_few_points")
        return TensileResult(
            None, None, disp, load, strain, stress,
            float(diameter_um) if has_d else None, area, L0,
            math.nan, math.nan, math.nan, math.nan, math.nan, math.nan,
            {}, max(disp.size - 1, 0), flags)

    # Fracture first: the rupture is a sudden load collapse (see _find_fracture).
    # Locating it before taking Fmax keeps a post-test recoil spike -- which the
    # raw global max would otherwise latch onto -- out of the breaking force. A
    # manual break_index (from the GUI) overrides the detector outright.
    if break_index is not None:
        break_idx = int(np.clip(break_index, 1, load.size - 1))
        fractured = True
        flags.append("manual_break")
    else:
        break_idx, fractured = _find_fracture(disp, load, cfg)
        if not fractured:
            flags.append("no_fracture_drop")

    # breaking force = the peak of the raw load up to the fracture ("the bump")
    pre = load[:break_idx + 1]
    i_peak = int(np.nanargmax(pre)) if np.isfinite(pre).any() else 0
    fmax = float(load[i_peak])
    # a larger raw spike after the break means we correctly ignored a recoil artifact
    if fractured and np.isfinite(load).any() and float(np.nanmax(load)) > fmax:
        flags.append("artifact_after_break")

    extension_at_break = float(disp[break_idx])
    strain_at_break = float(strain[break_idx])

    strength = fmax / area if has_d else math.nan

    if has_d:
        modulus, fit = _young_modulus(strain, stress, i_peak, cfg)
        sl = slice(0, break_idx + 1)
        toughness = float(np.trapz(stress[sl], strain[sl]))
        if fit and fit.get("r2", 1.0) < cfg.modulus_r2_min:
            flags.append("modulus_fit_weak")
    else:
        modulus, fit, toughness = math.nan, {}, math.nan
        flags.append("no_diameter")

    return TensileResult(
        None, None, disp, load, strain, stress,
        float(diameter_um) if has_d else None, area, L0,
        fmax, strength, extension_at_break, strain_at_break,
        modulus, toughness, fit, break_idx, flags)


def _find_fracture(disp: np.ndarray, load: np.ndarray,
                   cfg: CONFIG) -> tuple[int, bool]:
    """Locate the rupture as the first sudden load collapse; return (index, found).

    A single tensile pull ends when the fibre snaps: the load drops sharply over a
    short extension. Two pitfalls make the naive "global max, then first sample
    below 20% of it" rule fail on real crosshead data:

    * the recorded trace continues past the break, and the slack recoil produces a
      spike *larger* than the true peak -- so ``argmax``/``max`` land in the dead
      tail, and "post-peak" then scans only garbage; and
    * a brittle snap often leaves a high residual load (grip friction holding the
      broken end), so the load never falls below 20% of the peak at all.

    We therefore smooth the load (median, to kill lone spikes), take a *robust*
    peak ``P`` (99.5th percentile, immune to a brief recoil spike), arm only once
    the load has climbed past half of ``P``, then walk forward to the first sample
    where either the load has fallen by ``break_event_frac*P`` within
    ``break_window_mm`` of extension (a sudden collapse, whatever the residual) or
    it has dropped below ``break_drop_frac*P`` outright. Gradual ductile declines
    trip neither test, so they return ``(last, False)`` and are flagged upstream.
    """
    n = load.size
    if n == 0:
        return 0, False
    sm = pd.Series(load).rolling(5, center=True, min_periods=1).median().to_numpy()
    P = float(np.nanpercentile(sm, 99.5))
    if not np.isfinite(P) or P <= 0:
        return n - 1, False

    run = np.maximum.accumulate(np.nan_to_num(sm, nan=-np.inf))
    arm = 0.5 * P
    event = cfg.break_event_frac * P
    floor = cfg.break_drop_frac * P
    win = float(cfg.break_window_mm)

    k = 0  # left edge of the trailing extension window [disp[k], disp[j]]
    for j in range(1, n):
        while k < j and disp[j] - disp[k] > win:
            k += 1
        if run[j] < arm:                      # not yet climbed to a real load
            continue
        local_peak = float(np.nanmax(sm[k:j + 1]))
        if (local_peak - sm[j]) >= event or sm[j] < floor:
            return j, True
    return n - 1, False


def _linfit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Least-squares line fit returning (slope, intercept, R^2)."""
    if x.size < 2 or np.ptp(x) == 0:
        return math.nan, math.nan, 0.0
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return float(slope), float(intercept), float(r2)


def _young_modulus(strain: np.ndarray, stress: np.ndarray, i_peak: int,
                   cfg: CONFIG) -> tuple[float, dict]:
    """Steepest well-fit straight segment in the initial (pre-peak) region.

    Young's modulus is the slope of the stiffest linear stretch of the loading
    curve. We search the whole rise (origin..peak) rather than only the start,
    because some fibres have an initial slack/toe region where the load stays
    near zero and the real elastic rise comes later -- restricting to the first
    part would fit the slack and badly understate the modulus.

    Crosshead displacement is jittery, so a raw steepest-window search latches
    onto a noise-steepened sub-segment and overstates the modulus (taking the max
    over many noisy slope estimates is upward-biased). To counter that we smooth
    both axes with a rolling mean before fitting, and require a window wide enough
    to average over the jitter. A line is fit to each sliding window (with a
    stride) and the steepest slope among windows with R^2 >= ``modulus_r2_min``
    is returned (falling back to the steepest overall if none qualify). Bounded
    work per curve regardless of sample count.
    """
    hi = max(int(i_peak), 4)
    s = np.asarray(strain[:hi + 1], dtype=float)
    t = np.asarray(stress[:hi + 1], dtype=float)
    ok = np.isfinite(s) & np.isfinite(t)
    s, t = s[ok], t[ok]
    n = s.size
    if n < 5:
        return math.nan, {}

    # smooth both axes (~2% of the rising points) to suppress crosshead jitter
    smooth = min(max(7, n // 50), n)
    s = pd.Series(s).rolling(smooth, center=True, min_periods=1).mean().to_numpy()
    t = pd.Series(t).rolling(smooth, center=True, min_periods=1).mean().to_numpy()

    # window must be wide enough to average over noise, not chase a local spike
    win = min(max(12, int(round(cfg.modulus_window * n))), n)
    stride = max(1, win // 4)
    search_end = max(0, n - win)  # search the whole rise (handles slack/toe fibres)

    best_key, best = None, None
    for lo in range(0, search_end + 1, stride):
        hiw = min(lo + win, n)
        slope, intercept, r2 = _linfit(s[lo:hiw], t[lo:hiw])
        if not np.isfinite(slope):
            continue
        key = (r2 >= cfg.modulus_r2_min, slope)  # prefer qualifying, then steepest
        if best_key is None or key > best_key:
            best_key = key
            best = {"slope": slope, "intercept": intercept, "r2": r2,
                    "strain_lo": float(s[lo]), "strain_hi": float(s[hiw - 1])}
    if best is None:
        slope, intercept, r2 = _linfit(s[:win], t[:win])
        return float(slope), {"slope": slope, "intercept": intercept, "r2": r2,
                              "strain_lo": float(s[0]), "strain_hi": float(s[win - 1])}
    return float(best["slope"]), best


# --------------------------------------------------------------------------- #
# Export matrix                                                                #
# --------------------------------------------------------------------------- #
MATRIX_COLUMNS = [
    "group", "diameter_um", "area_um2", "gauge_length_mm", "fmax_N",
    "tensile_strength_MPa", "extension_at_break_mm", "strain_at_break_pct",
    "youngs_modulus_GPa", "toughness_MJ_m3", "flag", "notes",
]


def tensile_row(res: TensileResult, group: str | None = None,
                flag: str = "") -> dict:
    """Flatten a ``TensileResult`` to one export row (GPa / MPa / MJ·m⁻³ / %)."""
    area_um2 = res.area_m2 * 1e12 if np.isfinite(res.area_m2) else math.nan
    return {
        "group": group if group is not None else res.group,
        "diameter_um": res.diameter_um if res.diameter_um is not None else math.nan,
        "area_um2": area_um2,
        "gauge_length_mm": res.gauge_length_mm,
        "fmax_N": res.fmax_n,
        "tensile_strength_MPa": res.tensile_strength_pa / 1e6,
        "extension_at_break_mm": res.extension_at_break_mm,
        "strain_at_break_pct": res.strain_at_break * 100.0,
        "youngs_modulus_GPa": res.youngs_modulus_pa / 1e9,
        "toughness_MJ_m3": res.toughness_j_m3 / 1e6,
        "flag": flag,
        "notes": ";".join(res.flags),
    }


def _empty_row(group: str, diameter_um: float | None, flag: str) -> dict:
    """Row for a fibre with no usable curve (only an image, or a read error)."""
    return {
        "group": group,
        "diameter_um": diameter_um if diameter_um is not None else math.nan,
        "area_um2": math.nan, "gauge_length_mm": math.nan, "fmax_N": math.nan,
        "tensile_strength_MPa": math.nan, "extension_at_break_mm": math.nan,
        "strain_at_break_pct": math.nan, "youngs_modulus_GPa": math.nan,
        "toughness_MJ_m3": math.nan, "flag": flag, "notes": "",
    }


def build_matrix(diameters: dict[str, float], tensile: dict[str, Path],
                 cfg: CONFIG | None = None,
                 breaks: dict[str, int] | None = None) -> pd.DataFrame:
    """One-row-per-fibre metrics matrix over the union of both maps.

    ``diameters`` maps group -> mean diameter (microns); ``tensile`` maps group
    -> trace path. A fibre with only an image is flagged ``unmatched_tensile``;
    one with only a curve is flagged ``unmatched_image`` (force metrics still
    filled, stress metrics NaN). Read failures are flagged ``read_error: ...``.
    ``breaks`` optionally maps group -> manual break sample index (from the GUI),
    so an exported fibre uses the same fracture the user picked on screen.
    """
    cfg = cfg or CONFIG()
    breaks = breaks or {}
    groups = sorted(set(diameters) | set(tensile), key=natural_key)
    rows: list[dict] = []
    for g in groups:
        d = diameters.get(g)
        path = tensile.get(g)
        if path is None:
            rows.append(_empty_row(g, d, "unmatched_tensile"))
            continue
        try:
            res = compute_tensile(read_trace(path), d, cfg.gauge_length_mm, cfg,
                                  break_index=breaks.get(g))
        except Exception as exc:  # noqa: BLE001 - surface read/parse errors as a row
            rows.append(_empty_row(g, d, f"read_error: {exc}"))
            continue
        has_d = d is not None and np.isfinite(d) and float(d) > 0
        rows.append(tensile_row(res, group=g,
                                flag="" if has_d else "unmatched_image"))
    return pd.DataFrame(rows, columns=MATRIX_COLUMNS)
