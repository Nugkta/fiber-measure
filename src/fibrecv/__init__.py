"""fibrecv -- fibre diameter profiling from optical microscopy images.

Three-stage pipeline:
  * ``run_measure``   : per-image detection -> profile CSV/plot/overlay/meta
  * ``run_aggregate`` : group replicates by filename -> registered mean+/-variance curve
  * ``run_xsection``  : multi-angle (C1) groups -> per-position ellipse cross-sections

See module docstrings for the per-stage contracts.
"""

__version__ = "0.1.0"


def main() -> None:
    print("fibrecv: use `python -m fibrecv.run_measure`, "
          "`python -m fibrecv.run_aggregate` or "
          "`python -m fibrecv.run_xsection`.")
