"""Sharper angle-spacing tests: edge-mirror correlation + ellipse projection fit.

Caches per-image (top, bot) edge profiles as npz so reruns are cheap.
"""
import re
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.optimize import least_squares

DATA = Path("/Users/stan/Documents/UOM/spins/multiangle/C1")
CACHE = Path(__file__).parent / "cache"
CACHE.mkdir(exist_ok=True)
Image.MAX_IMAGE_PIXELS = None

MARGIN = 200
SMOOTH = 31
ROUGH_SMOOTH = 301    # long window: subtract slow trend to get roughness
MAX_SHIFT = 400


def edges(path):
    npz = CACHE / (path.stem + ".npz")
    if npz.exists():
        d = np.load(npz)
        return d["top"], d["bot"]
    g = np.asarray(Image.open(path).convert("L"), dtype=np.float32)
    H, W = g.shape
    lo, hi = np.percentile(g, [20, 90])
    thr = (lo + hi) / 2.0
    mask = g > thr
    top = np.full(W, np.nan)
    bot = np.full(W, np.nan)
    for x in range(MARGIN, W - MARGIN):
        col = mask[:, x]
        if not col.any():
            continue
        d = np.diff(np.concatenate(([0], col.view(np.int8), [0])))
        starts = np.flatnonzero(d == 1)
        ends = np.flatnonzero(d == -1)
        k = np.argmax(ends - starts)
        t, b = starts[k], ends[k] - 1
        if t <= 2 or b >= H - 3:
            continue
        top[x] = t
        bot[x] = b
    np.savez_compressed(npz, top=top, bot=bot)
    return top, bot


def smooth(a, n):
    out = np.convolve(np.nan_to_num(a, nan=0.0), np.ones(n) / n, mode="same")
    norm = np.convolve(~np.isnan(a), np.ones(n) / n, mode="same")
    with np.errstate(invalid="ignore"):
        return np.where(norm > 0.5, out / np.maximum(norm, 1e-9), np.nan)


def roughness(a):
    """Short-scale edge detail: smoothed profile minus long-window trend."""
    return smooth(a, SMOOTH) - smooth(a, ROUGH_SMOOTH)


def aligned_corr(a, b, max_shift=MAX_SHIFT, step=4):
    best_r, best_s = -2.0, 0
    for s in range(-max_shift, max_shift + 1, step):
        if s >= 0:
            x, y = a[s:], b[: len(b) - s]
        else:
            x, y = a[: len(a) + s], b[-s:]
        m = ~np.isnan(x) & ~np.isnan(y)
        if m.sum() < 500:
            continue
        r = np.corrcoef(x[m], y[m])[0, 1]
        if r > best_r:
            best_r, best_s = r, s
    return best_r, best_s


def ellipse_resid(widths, step_deg):
    """RMS residual of fitting w(theta)=sqrt(a^2 cos^2 + b^2 sin^2) at the
    hypothesized angular positions k*step_deg."""
    th = np.deg2rad(step_deg) * np.arange(6)
    w = np.asarray(widths, float)

    def model(p):
        a, b, phi = p
        return np.sqrt((a * np.cos(th + phi)) ** 2 + (b * np.sin(th + phi)) ** 2)

    p0 = [w.max(), w.min(), 0.0]
    try:
        res = least_squares(lambda p: model(p) - w, p0, method="lm",
                            max_nfev=2000)
        return np.sqrt(np.mean(res.fun ** 2))
    except Exception:
        return np.nan


def main(fibers):
    files = sorted(DATA.glob("C1_*_a*_part*.tiff"))
    files = [f for f in files if not f.stem.endswith("s")]
    pat = re.compile(r"C1_(\d+)_a(\d)_part(\d)$")

    data = {}  # (fib, part) -> {ang: (top, bot)}
    for f in files:
        m = pat.match(f.stem)
        if not m:
            continue
        fib, ang, part = int(m.group(1)), int(m.group(2)), int(m.group(3))
        if fib not in fibers:
            continue
        data.setdefault((fib, part), {})[ang] = edges(f)
        sys.stderr.write(".")
        sys.stderr.flush()
    sys.stderr.write("\n")

    same = np.zeros((6, 6)); mirr = np.zeros((6, 6)); n = np.zeros((6, 6))
    for (fib, part), d in sorted(data.items()):
        rough = {a: (roughness(t), roughness(b)) for a, (t, b) in d.items()}
        # align on width deviation first, reuse shift for edge tests
        wdev = {a: smooth(b - t, SMOOTH) - np.nanmedian(b - t)
                for a, (t, b) in d.items()}
        for i in range(1, 7):
            for j in range(1, 7):
                if i == j or i not in d or j not in d:
                    continue
                _, s = aligned_corr(wdev[i], wdev[j])
                ti, bi = rough[i]
                tj, bj = rough[j]
                # same-side: top_i vs top_j ; mirror: top_i vs -bot_j
                r_same, _ = aligned_corr(ti, tj, max_shift=abs(s) + 60, step=2)
                r_mirr, _ = aligned_corr(ti, -bj, max_shift=abs(s) + 60, step=2)
                same[i - 1, j - 1] += r_same
                mirr[i - 1, j - 1] += r_mirr
                n[i - 1, j - 1] += 1

    with np.errstate(invalid="ignore"):
        S = same / np.maximum(n, 1)
        Mi = mirr / np.maximum(n, 1)

    def show(mat, name):
        print(f"\n{name} (row=i col=j):")
        print("      " + "  ".join(f"a{j}   " for j in range(1, 7)))
        for i in range(6):
            print(f"a{i+1}  " + "  ".join(
                "  --  " if i == j else f"{mat[i, j]:+.3f}" for j in range(6)))

    show(S, "same-side edge-roughness corr: top(a_i) vs top(a_j)")
    show(Mi, "mirror edge-roughness corr: top(a_i) vs -bottom(a_j)")

    by_lag_s = {}; by_lag_m = {}
    for i in range(6):
        for j in range(6):
            if i == j:
                continue
            lag = min((i - j) % 6, (j - i) % 6)
            by_lag_s.setdefault(lag, []).append(S[i, j])
            by_lag_m.setdefault(lag, []).append(Mi[i, j])
    print("\nmean by angular index lag:")
    print("lag  same-side  mirror")
    for lag in sorted(by_lag_s):
        print(f"{lag}    {np.mean(by_lag_s[lag]):+.3f}     "
              f"{np.mean(by_lag_m[lag]):+.3f}")

    # ellipse projection fit on per-part median widths
    print("\nellipse fit RMS residual (px) per (fiber,part), 30 vs 60 deg steps:")
    r30_all, r60_all = [], []
    for (fib, part), d in sorted(data.items()):
        if len(d) < 6:
            continue
        w = [np.nanmedian(d[a][1] - d[a][0]) for a in range(1, 7)]
        r30 = ellipse_resid(w, 30)
        r60 = ellipse_resid(w, 60)
        r30_all.append(r30); r60_all.append(r60)
    print(f"  mean over {len(r30_all)} parts:  step30 {np.mean(r30_all):.2f} px"
          f"   step60 {np.mean(r60_all):.2f} px")
    print(f"  median:                 step30 {np.median(r30_all):.2f} px"
          f"   step60 {np.median(r60_all):.2f} px")


if __name__ == "__main__":
    fibers = [int(a) for a in sys.argv[1:]] or [1, 2, 3]
    main(fibers)
