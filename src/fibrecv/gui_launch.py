"""Console-script launcher that starts the Streamlit GUI.

Dependencies
------------
``streamlit`` (invoked as a subprocess) and the standard library only.

Inputs
------
None required. Any extra CLI args after ``fibrecv-gui`` are forwarded verbatim to
``streamlit run`` (e.g. ``fibrecv-gui --server.port 8600``).

Output
------
Launches the local Streamlit server running ``gui_app.py`` (opens
``http://localhost:8501`` in the default browser). Returns Streamlit's exit code.

Pos
---
The ``[project.scripts] fibrecv-gui`` entry point. A thin wrapper so users can
type ``fibrecv-gui`` instead of the full ``streamlit run .../gui_app.py`` path;
the GUI logic itself lives in ``gui_app.py``.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    """Shell ``streamlit run gui_app.py``, forwarding any extra CLI args."""
    try:
        import streamlit  # noqa: F401
    except ImportError:
        print(
            "streamlit is not installed. Inside the fibrecv folder run "
            "`pip install -e .` (or `uv sync`) first, then `fibrecv-gui`.",
            file=sys.stderr,
        )
        return 1

    app = Path(__file__).resolve().parent / "gui_app.py"
    cmd = [sys.executable, "-m", "streamlit", "run", str(app), *sys.argv[1:]]
    return subprocess.call(cmd)


if __name__ == "__main__":
    raise SystemExit(main())
