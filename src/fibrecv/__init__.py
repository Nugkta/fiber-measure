"""fibrecv -- fibre diameter profiling from MasP2 microscopy images.

Two-stage pipeline:
  * ``run_measure``   : per-image detection -> profile CSV/plot/overlay/meta
  * ``run_aggregate`` : group replicates by A_B -> registered mean+/-variance curve

See module docstrings for the per-stage contracts.
"""

__version__ = "0.1.0"


def main() -> None:
    print("fibrecv: use `python -m fibrecv.run_measure` or "
          "`python -m fibrecv.run_aggregate`.")
