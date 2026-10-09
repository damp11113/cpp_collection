#!/usr/bin/env python3
"""Bjontegaard delta rate from bench/compare_codecs.sh results.

    python3 bench/bd_rate.py results.txt [more results ...] [--ref mjpeg] [--psnr 38]

For each codec: the average bitrate difference against the reference codec at the
same PSNR-Y (negative = fewer bits), the bitrate needed for --psnr dB, and the
encode time there. With several result files (clips) the BD-rates are averaged.
Needs numpy.
"""
import argparse
import collections
import math

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('files', nargs='+')
ap.add_argument('--ref', default='mjpeg')
ap.add_argument('--psnr', type=float, default=38.0)
args = ap.parse_args()

order = ['mjpeg', 'h261', 'h263', 'h263p', 'mpeg1', 'mpeg2', 'mpeg4', 'x264', 'tjc2', 'tjc3', 'tjc4']


def load(path):
    d = collections.defaultdict(list)
    for line in open(path):
        f = line.split()
        if len(f) != 7 or f[0].startswith('#') or 'inf' in f[4]:
            continue
        d[f[0]].append((float(f[3]), float(f[4]), float(f[6])))  # kbit/s, PSNR-Y, seconds
    return d


def hull(pts):
    out = []
    for p in sorted(pts):
        if not out or p[1] > out[-1][1]:
            out.append(p)
    return out


def bd(a, b):
    a, b = hull(a), hull(b)
    if len(a) < 3 or len(b) < 3:
        return None
    lo = max(a[0][1], b[0][1])
    hi = min(a[-1][1], b[-1][1])
    if hi - lo < 1.0:
        return None
    fa = np.polyfit([p[1] for p in a], np.log([p[0] for p in a]), min(3, len(a) - 1))
    fb = np.polyfit([p[1] for p in b], np.log([p[0] for p in b]), min(3, len(b) - 1))
    ia, ib = np.polyint(fa), np.polyint(fb)
    avg = ((np.polyval(ib, hi) - np.polyval(ib, lo)) - (np.polyval(ia, hi) - np.polyval(ia, lo))) / (hi - lo)
    return (math.exp(avg) - 1) * 100


def rate_at(pts, target):
    c = hull(pts)
    for p, q in zip(c, c[1:]):
        if p[1] <= target <= q[1]:
            f = (target - p[1]) / (q[1] - p[1])
            return math.exp(math.log(p[0]) + f * (math.log(q[0]) - math.log(p[0]))), p[2] + f * (q[2] - p[2])
    return None


clips = [load(f) for f in args.files]
codecs = [c for c in order if any(c in d for d in clips)]
print(f"BD-rate vs {args.ref} (PSNR-Y; negative = fewer bits for the same quality)")
print(f"{'codec':8} {'BD-rate':>9}  " + ''.join(f"{'kbit/s@%gdB' % args.psnr:>14}{'enc s':>7}" for _ in clips))
for c in codecs:
    vals = [bd(d[args.ref], d[c]) for d in clips if args.ref in d and c in d]
    vals = [v for v in vals if v is not None]
    row = f"{c:8} " + (f"{sum(vals) / len(vals):>8.1f}% " if vals else f"{'n/a':>9} ")
    for d in clips:
        r = rate_at(d.get(c, []), args.psnr)
        row += f"{r[0]:>14.0f}{r[1]:>7.2f}" if r else f"{'n/a':>14}{'':>7}"
    print(row)
