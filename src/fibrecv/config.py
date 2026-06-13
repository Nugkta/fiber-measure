"""Central configuration for the fibre-diameter pipeline.

Dependencies
------------
Standard library only (``dataclasses``). No third-party imports so this module
can be imported cheaply by every other module and by both CLIs.

Inputs
------
None at import time. ``CONFIG`` instances are built with defaults (calibrated
from 5 representative MasP2 images) and selectively overridden by CLI flags in
``run_measure.py`` / ``run_aggregate.py``.

Output
------
The frozen-ish ``CONFIG`` dataclass carrying every tunable parameter, plus the
``px_to_um`` helper.

Pos
---
Bottom of the dependency graph. Imported by features, band, edges, qc, overlay,
register, measure and both run_* CLIs. Changing a default here changes the whole
pipeline's behaviour; the load-bearing strictness knob is ``edge_frac``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass
class CONFIG:
    """All default parameters for detection, QC and registration.

    The defaults come from the empirical calibration in the plan: background
    HSV saturation is stable (~0.45-0.48) while fibre-core saturation is
    0.10-0.36, so detection thresholds operate on a per-image, self-normalising
    desaturation z-map ``D`` rather than on absolute brightness.
    """

    # --- calibration ---
    ppu: float = 1.3680  # pixels per micron; diameter_um = diameter_px / ppu

    # --- desaturation feature / background estimate ---
    margin: float = 0.12       # top+bottom fraction of rows treated as background
    eps: float = 1e-6          # numerical floor
    mad_scale: float = 1.4826  # MAD -> sigma conversion for the robust z-map

    # --- band localisation ---
    k_band: float = 4.0        # D threshold for the coarse band mask (z units)
    min_width: float = 0.85    # component must span >= this fraction of image width
    close_width: int = 25      # horizontal morphological closing length (px)
    min_object: int = 800      # remove_small_objects min area (px)

    # --- centerline ---
    reject_dev: float = 25.0   # reject columns whose band center deviates > this (px)

    # --- per-column edge detection (the shadow-critical recipe) ---
    wcol: int = 41             # column-neighbourhood width averaged before edging
    #                            (default 41; >=15 needed to suppress internal
    #                             iridescence banding)
    sigma_y: float = 4.0       # vertical Gaussian smoothing sigma (px); 4.0
    #                            merges the plateau-fragmented ramps of
    #                            defocused walls (2.0 left them speckled;
    #                            recalibrated on the full 144-image MasP2 set)
    # STRICTNESS: boundary placed where D crosses level = local_bg + min(edge_z, edge_frac*A).
    # edge_z is the PRIMARY knob (absolute z-units above background): higher -> tighter,
    # and stable because it does not track the per-column specular peak. edge_frac is a
    # relative CAP that keeps the level inside the hump for faint fibres (low amplitude A).
    edge_z: float = 4.0        # STRICTNESS KNOB (z above wall-local bg); higher -> tighter
    #                            (user-validated on pilot 3_1: ez4 tracks the true wall;
    #                             higher values get dragged inward by internal reflections)
    edge_frac: float = 0.65    # relative cap on the level for faint/weak walls
    guard: int = 12            # px just outside the wall for local bg + recovery checks
    # --- wall finding (separates true fibre walls from shadow/vignette ramps) ---
    slope_min: float = 0.05    # absolute slope floor (z/px) for a wall candidate
    slope_rel: float = 0.2     # wall must reach this fraction of the side's max slope
    slope_cap: float = 0.12    # absolute slope that always qualifies (caps the relative
    #                            gate so a defocus-softened true wall is not skipped when
    #                            a sharp internal reflection dominates the side's max)
    rise_min: float = 2.0      # wall run must rise at least this many z-units
    wall_gap_frac: float = 0.12  # plateau-bridging length inside a wall run, as a
    #                              fraction of band thickness (clipped 4..16 px in
    #                              edges.py); heals defocus-fragmented soft walls
    #                              without merging distinct features
    band_ratio_min: float = 0.5  # median diameter / coarse-band thickness below this
    #                              -> low_confidence (detector likely locked onto an
    #                              internal feature inside a defocused blur band)
    band_ratio_max: float = 1.6  # ...and above this -> band_mismatch too (boundary
    #                              likely grabbed a shadow/halo outside the fibre)
    amin: float = 3.0          # minimum band amplitude A (z units) to accept a column

    # --- search window sizing (multiples of detected band thickness) ---
    window_thick_mult: float = 1.5  # half-window = mult*band_half + window_pad + 3*guard
    window_pad: int = 10

    # --- QC / smoothing ---
    min_coverage: float = 0.5  # below this fraction of valid columns -> low_confidence
    roll_window: int = 51      # rolling-MAD outlier window (odd, px)
    roll_k: float = 5.0        # rolling-MAD outlier threshold
    median_k: int = 11         # median-filter kernel for the smoothed profile (odd)
    savgol_window: int = 31    # Savitzky-Golay window (odd)
    savgol_poly: int = 3       # Savitzky-Golay polynomial order

    # --- replicate registration ---
    max_shift: int = 400       # bound on cross-correlation lag (px)
    min_corr: float = 0.3      # normalised corr-peak below this -> registration_uncertain

    # --- tensile (stress-strain) analysis ---
    # The tensile tester records only crosshead displacement (ΔL) and force, so
    # strain = ΔL / gauge_length. Cross-section A = π(d/2)² uses the matched image
    # group's measured mean diameter, tying the two batches together.
    gauge_length_mm: float = 10.0  # grip separation L0; strain = disp_mm / this
    modulus_window: float = 0.03   # sliding-fit width as a fraction of strain-to-peak
    #                                (auto Young's modulus = steepest well-fit segment)
    modulus_r2_min: float = 0.98   # min R² for a window to count as "linear"; the
    #                                steepest qualifying window sets the modulus
    # Fracture detection works on a median-smoothed load and a *robust* peak P
    # (99.5th percentile, so a single post-test spike cannot become Fmax). The
    # rupture is the first sudden collapse: the load falls by >= break_event_frac*P
    # within break_window_mm of extension (catches brittle snaps even when grip
    # friction leaves a high residual), OR drops below break_drop_frac*P outright.
    break_window_mm: float = 0.05  # extension over which a catastrophic drop completes
    break_event_frac: float = 0.50  # sudden fall (fraction of robust peak) = a snap
    break_drop_frac: float = 0.20  # absolute collapse: load < this·P also marks fracture

    def px_to_um(self, px: float) -> float:
        """Convert a pixel length to microns using the calibration."""
        return px / self.ppu

    def as_dict(self) -> dict:
        """Plain-dict snapshot for provenance JSON."""
        return asdict(self)
