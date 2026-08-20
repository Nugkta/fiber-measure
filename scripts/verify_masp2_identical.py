"""Verify the desat (MasP2) path is bit-identical between the base branch and HEAD.

Input
-----
- ``/Users/stan/Documents/UOM/spins/Images MasP2`` (144 ``*.jpg``).
- The ``fibrecv`` package resolved from the *current working directory's*
  worktree (each worktree's ``.venv`` carries a ``fibrecv.pth`` pointing at its
  own ``src/``), so running this file with ``uv run`` from a given worktree
  measures that worktree's code.

Output
------
A JSON side-record (``--out``) holding, per image: the full ``repr`` of the
median raw diameter (px), the coverage, and SHA-256 digests of
``res.diameter_raw.tobytes()`` and of the ``D`` z-map bytes -- plus the git
rev, the resolved ``fibrecv`` path and the chosen image indices.
``--compare A B`` merges two such side-records into the final verdict JSON.

Pos
---
Study-03 reproducibility artifact. Substantiates (or refutes) the claim that
the median-RGB bright z-map left the desat/MasP2 pipeline untouched. Reads
only; never writes into the package.

Usage
-----
    # in the HEAD worktree
    uv run python scripts/verify_masp2_identical.py --out /tmp/masp2_head.json
    # in the base worktree (worktree-multiangle-xsection)
    uv run python <this file> --out /tmp/masp2_base.json
    # then, from anywhere
    uv run python scripts/verify_masp2_identical.py \
        --compare /tmp/masp2_base.json /tmp/masp2_head.json \
        --out scripts/masp2_identical.json

Reminder: once I am updated, update my header comments and the folder's md.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

MASP2_ROOT = Path("/Users/stan/Documents/UOM/spins/Images MasP2")


def _sha256(arr: np.ndarray) -> str:
    """SHA-256 of an array's raw bytes (C-contiguous copy if needed)."""
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()


def _git(*args: str) -> str:
    try:
        return subprocess.run(["git", *args], capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception as exc:  # pragma: no cover - provenance only
        return f"<unavailable: {exc}>"


def pick_images(root: Path) -> list[tuple[int, Path]]:
    """Deterministic 3-image pick: sorted ``*.jpg``, indices 0, len//2, -1."""
    imgs = sorted(p for p in root.glob("*.jpg"))
    if len(imgs) < 3:
        raise SystemExit(f"need >=3 jpgs under {root}, found {len(imgs)}")
    idx = [0, len(imgs) // 2, len(imgs) - 1]
    return [(i, imgs[i]) for i in idx]


def measure(out_path: Path) -> dict:
    """Run ``compute_measurement(rgb, CONFIG(), name)`` on the 3 picked images."""
    from fibrecv.compute import compute_measurement
    from fibrecv.config import CONFIG
    from fibrecv.io_utils import load_rgb
    import fibrecv

    picks = pick_images(MASP2_ROOT)
    rec = {
        "kind": "side_record",
        "cwd": str(Path.cwd()),
        "git_branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "git_rev": _git("rev-parse", "HEAD"),
        "fibrecv_path": str(Path(fibrecv.__file__).resolve().parent),
        "masp2_root": str(MASP2_ROOT),
        "n_jpg_total": len(sorted(MASP2_ROOT.glob("*.jpg"))),
        "picked_indices": [i for i, _ in picks],
        "config": "CONFIG()  (dataclass defaults, feature_mode='desat')",
        "images": [],
    }
    for i, p in picks:
        rgb = load_rgb(str(p))
        mr = compute_measurement(rgb, CONFIG(), p.stem)
        med = float(np.nanmedian(mr.res.diameter_raw))
        rec["images"].append({
            "index": i,
            "file": p.name,
            "stem": p.stem,
            "image_shape": list(rgb.shape),
            "median_diameter_raw_px_repr": repr(med),
            "median_diameter_raw_px": med,
            "median_diameter_um": mr.meta["median_diameter_um"],
            "coverage": float(mr.meta["coverage"]),
            "n_valid": int(mr.meta["n_valid"]),
            "sha256_diameter_raw": _sha256(mr.res.diameter_raw),
            "sha256_D_zmap": _sha256(mr.D),
            "D_dtype": str(mr.D.dtype),
            "diameter_raw_dtype": str(mr.res.diameter_raw.dtype),
        })
        print(f"[{i:>3}] {p.name}  med={med!r}  cov={mr.meta['coverage']:.6f}\n"
              f"      sha(diam)={rec['images'][-1]['sha256_diameter_raw']}\n"
              f"      sha(D)   ={rec['images'][-1]['sha256_D_zmap']}")
        sys.stdout.flush()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(rec, indent=2))
    print(f"\nside record -> {out_path}")
    return rec


def compare(a_path: Path, b_path: Path, out_path: Path) -> int:
    """Merge two side records and emit the identity verdict."""
    a = json.loads(Path(a_path).read_text())
    b = json.loads(Path(b_path).read_text())
    if a["picked_indices"] != b["picked_indices"]:
        raise SystemExit("side records disagree on picked indices")

    per_image, all_ok = [], True
    for ia, ib in zip(a["images"], b["images"]):
        assert ia["file"] == ib["file"]
        same_d = ia["sha256_diameter_raw"] == ib["sha256_diameter_raw"]
        same_z = ia["sha256_D_zmap"] == ib["sha256_D_zmap"]
        same_med = ia["median_diameter_raw_px_repr"] == ib["median_diameter_raw_px_repr"]
        same_cov = ia["coverage"] == ib["coverage"]
        ok = same_d and same_z and same_med and same_cov
        all_ok &= ok
        per_image.append({
            "index": ia["index"], "file": ia["file"],
            "sha256_diameter_raw_match": same_d,
            "sha256_D_zmap_match": same_z,
            "median_repr_match": same_med,
            "coverage_match": same_cov,
            "bit_identical": ok,
            "base": {k: ia[k] for k in ("median_diameter_raw_px_repr", "coverage",
                                        "sha256_diameter_raw", "sha256_D_zmap")},
            "head": {k: ib[k] for k in ("median_diameter_raw_px_repr", "coverage",
                                        "sha256_diameter_raw", "sha256_D_zmap")},
        })

    verdict = {
        "kind": "comparison",
        "question": "Is the desat (MasP2) path bit-identical between base and HEAD?",
        "base": {k: a[k] for k in ("git_branch", "git_rev", "fibrecv_path", "cwd")},
        "head": {k: b[k] for k in ("git_branch", "git_rev", "fibrecv_path", "cwd")},
        "masp2_root": a["masp2_root"],
        "n_jpg_total": a["n_jpg_total"],
        "picked_indices": a["picked_indices"],
        "picked_files": [im["file"] for im in a["images"]],
        "config": a["config"],
        "per_image": per_image,
        "bit_identical_all": bool(all_ok),
        "verdict": ("MasP2/desat path is BIT-IDENTICAL between base and HEAD"
                    if all_ok else
                    "MasP2/desat path DIFFERS between base and HEAD"),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(verdict, indent=2))
    for r in per_image:
        print(f"{r['file']}: diam_sha={r['sha256_diameter_raw_match']} "
              f"D_sha={r['sha256_D_zmap_match']} med={r['median_repr_match']} "
              f"cov={r['coverage_match']}")
    print(f"\n{verdict['verdict']}\n-> {out_path}")
    return 0 if all_ok else 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--compare", nargs=2, type=Path, default=None,
                    metavar=("BASE_JSON", "HEAD_JSON"))
    args = ap.parse_args(argv)
    if args.compare:
        return compare(args.compare[0], args.compare[1], args.out)
    measure(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
