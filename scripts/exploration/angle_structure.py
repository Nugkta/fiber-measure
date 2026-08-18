"""Infer the angular spacing of the multiangle dataset.

Physics: the projected silhouette of a fiber rotated about its long axis is
180-degree periodic. If the 6 angles are 60-degree steps over 360 deg, angle
pairs (a1,a4), (a2,a5), (a3,a6) are the same projection (vertically mirrored)
and their width profiles should correlate far above other pairs. If the 6
angles are 30-degree steps over 180 deg, no pair repeats.
"""
import re
import sys
from pathlib import Path

import numpy as np
from PIL import Image

DATA = Path("/Users/stan/Documents/UOM/spins/multiangle/C1")
Image.MAX_IMAGE_PIXELS = None

MARGIN = 200          # ignore columns near left/right image border (vignetting)
SMOOTH = 31           # columns, moving-average window for profile smoothing
MAX_SHIFT = 400       # px, alignment search range between angle shots


def width_profile(path):
    """Per-column projected width (px) of the fiber, NaN where unreliable."""
    g = np.asarray(Image.open(path).convert("L"), dtype=np.float32)
    H, W = g.shape
    # global threshold: midpoint between dark background and bright fiber
    lo, hi = np.percentile(g, [20, 90])
    thr = (lo + hi) / 2.0
    mask = g > thr
    w = np.full(W, np.nan)
    top = np.full(W, np.nan)
    bot = np.full(W, np.nan)
    idx = np.arange(H)
    for x in range(MARGIN, W - MARGIN):
        col = mask[:, x]
        if not col.any():
            continue
        # longest True run
        d = np.diff(np.concatenate(([0], col.view(np.int8), [0])))
        starts = np.flatnonzero(d == 1)
        ends = np.flatnonzero(d == -1)
        k = np.argmax(ends - starts)
        t, b = starts[k], ends[k] - 1
        if t <= 2 or b >= H - 3:            # touches frame -> unreliable
            continue
        w[x] = b - t
        top[x] = t
        bot[x] = b
    return w, top, bot


def smooth(a, n=SMOOTH):
    out = np.convolve(np.nan_to_num(a, nan=0.0), np.ones(n) / n, mode="same")
    norm = np.convolve(~np.isnan(a), np.ones(n) / n, mode="same")
    with np.errstate(invalid="ignore"):
        return np.where(norm > 0.5, out / np.maximum(norm, 1e-9), np.nan)


def aligned_corr(a, b, max_shift=MAX_SHIFT):
    """Max Pearson r between profiles a and b over integer shifts."""
    best = -2.0, 0
    for s in range(-max_shift, max_shift + 1, 4):
        if s >= 0:
            x, y = a[s:], b[: len(b) - s]
        else:
            x, y = a[: len(a) + s], b[-s:]
        m = ~np.isnan(x) & ~np.isnan(y)
        if m.sum() < 500:
            continue
        r = np.corrcoef(x[m], y[m])[0, 1]
        if r > best[0]:
            best = (r, s)
    return best


def main(fibers):
    files = sorted(DATA.glob("C1_*_a*_part*.tiff"))
    files = [f for f in files if not f.stem.endswith("s")]
    pat = re.compile(r"C1_(\d+)_a(\d)_part(\d)$")

    # profiles[(fiber, part)][angle] = deviation profile
    profiles = {}
    widths = {}  # mean width per (fiber, angle)
    for f in files:
        m = pat.match(f.stem)
        if not m:
            continue
        fib, ang, part = int(m.group(1)), int(m.group(2)), int(m.group(3))
        if fib not in fibers:
            continue
        w, top, bot = width_profile(f)
        ws = smooth(w)
        dev = ws - np.nanmedian(ws)            # width deviation profile
        profiles.setdefault((fib, part), {})[ang] = dev
        widths.setdefault((fib, ang), []).append(np.nanmedian(w))
        sys.stderr.write(".")
        sys.stderr.flush()
    sys.stderr.write("\n")

    # pairwise angle similarity, aggregated over (fiber, part)
    R = np.zeros((6, 6))
    N = np.zeros((6, 6))
    for (fib, part), d in sorted(profiles.items()):
        for i in range(1, 7):
            for j in range(i + 1, 7):
                if i in d and j in d:
                    r, s = aligned_corr(d[i], d[j])
                    if r > -2:
                        R[i - 1, j - 1] += r
                        N[i - 1, j - 1] += 1

    with np.errstate(invalid="ignore"):
        M = R / np.maximum(N, 1)
    print("mean aligned correlation of width-deviation profiles (upper tri):")
    print("      " + "  ".join(f"a{j}   " for j in range(2, 7)))
    for i in range(5):
        row = ["      " * i]
        for j in range(i + 1, 6):
            row.append(f"{M[i, j]:+.3f}")
        print(f"a{i + 1}  " + "  ".join(row))

    print("\npair means: 180-deg-pair hypothesis vs rest")
    pair_ids = [(0, 3), (1, 4), (2, 5)]
    inpair = [M[i, j] for i, j in pair_ids]
    rest = [M[i, j] for i in range(5) for j in range(i + 1, 6)
            if (i, j) not in pair_ids]
    print(f"  (a1,a4),(a2,a5),(a3,a6): {np.mean(inpair):+.3f}  "
          f"individually {['%+.3f' % v for v in inpair]}")
    print(f"  other 12 pairs:          {np.mean(rest):+.3f}  "
          f"(min {min(rest):+.3f}, max {max(rest):+.3f})")

    print("\nmedian width (px) per fiber x angle:")
    fibs = sorted({f for f, _ in widths})
    print("fib  " + "  ".join(f"  a{a}  " for a in range(1, 7)))
    for fib in fibs:
        row = [f"{np.mean(widths.get((fib, a), [np.nan])):6.1f}"
               for a in range(1, 7)]
        print(f"{fib:02d}   " + "  ".join(row))


if __name__ == "__main__":
    fibers = [int(a) for a in sys.argv[1:]] or [1, 2, 3]
    main(fibers)
