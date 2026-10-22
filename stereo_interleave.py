#!/usr/bin/env python3
"""
Stereo bin-interleave codec.

Encode: transform both stereo channels to the frequency domain, then keep
only every other bin from each channel:

    bin 1 -> channel 1 (L)
    bin 2 -> channel 2 (R)
    bin 3 -> channel 1 (L)
    bin 4 -> channel 2 (R)
    ...    up to Nyquist or --max-freq

so the file holds a single spectrum the size of one channel.

The encoder also stores a small side channel: the true L/R level ratio of
every band of --band-width bins, as one signed byte (0.5 dB steps) per band
per frame.

Decode: every channel is missing every other bin. Rebuild modes:

    avg    - missing bin = mean of the neighbouring bins of the same channel
             (bin 2 of L = mean(bin 1 of L, bin 3 of L))
    blind  - missing bin = the other channel's bin at that frequency, scaled
             by an L/R ratio guessed from the bins each channel kept
    cross  - missing bin = the other channel's bin at that frequency, scaled
             by the stored per-band L/R ratio (default, best quality)

Transforms:
    fft    - STFT, sine window, 50% overlap (complex bins)
    mdct   - MDCT, sine window, critically sampled (real bins, half the data)

Bin width is sample_rate / frame Hz for fft, so --frame 44100 at 44.1 kHz
gives exact 1 Hz bins (with heavy time smearing).

Usage:
    python3 stereo_interleave.py encode in.wav out.npz [--transform fft] [--frame 2048] [--max-freq 16000] [--band-width 4]
    python3 stereo_interleave.py decode in.npz out.wav [--rebuild cross] [--width 8]
    python3 stereo_interleave.py test [in.wav] [encode/decode options]

Requires numpy. Input must be a PCM stereo WAV (8/16/24/32-bit integer).
"""

import argparse
import io
import os
import sys
import wave

import numpy as np


# ---------------------------------------------------------------- WAV I/O

def read_wav(path):
    with wave.open(path, "rb") as w:
        nch, width, sr, n = w.getnchannels(), w.getsampwidth(), w.getframerate(), w.getnframes()
        raw = w.readframes(n)

    if width == 1:
        x = (np.frombuffer(raw, np.uint8).astype(np.float64) - 128) / 128
    elif width == 2:
        x = np.frombuffer(raw, "<i2").astype(np.float64) / 32768
    elif width == 3:
        b = np.frombuffer(raw, np.uint8).reshape(-1, 3).astype(np.int32)
        v = b[:, 0] | (b[:, 1] << 8) | (b[:, 2] << 16)
        v = np.where(v & 0x800000, v - 0x1000000, v)
        x = v.astype(np.float64) / 8388608
    elif width == 4:
        x = np.frombuffer(raw, "<i4").astype(np.float64) / 2147483648
    else:
        raise ValueError(f"unsupported sample width: {width} bytes")

    return x.reshape(-1, nch), sr


def write_wav(path, x, sr):
    pcm = np.clip(np.round(x * 32767), -32768, 32767).astype("<i2")
    with wave.open(path, "wb") as w:
        w.setnchannels(x.shape[1])
        w.setsampwidth(2)
        w.setframerate(sr)
        w.writeframes(pcm.tobytes())


# ---------------------------------------------------------------- transforms

def sine_window(n):
    # w[i]^2 + w[i + n/2]^2 == 1, so it gives perfect reconstruction at 50% overlap
    return np.sin(np.pi * (np.arange(n) + 0.5) / n)


def frame_signal(x, n):
    hop = n // 2
    n_frames = -(-(len(x) + hop) // hop)
    padded = np.zeros((n_frames + 1) * hop)
    padded[hop:hop + len(x)] = x
    idx = np.arange(n)[None, :] + hop * np.arange(n_frames)[:, None]
    return padded[idx]


def overlap_add(frames, n, length):
    hop = n // 2
    out = np.zeros((len(frames) + 1) * hop)
    for i, f in enumerate(frames):
        out[i * hop:i * hop + n] += f
    return out[hop:hop + length]


def mdct_basis(n):
    m = n // 2
    k = np.arange(m)[:, None]
    i = np.arange(n)[None, :]
    return np.cos(np.pi / m * (i + 0.5 + m / 2) * (k + 0.5))


def forward(x, n, transform):
    frames = frame_signal(x, n) * sine_window(n)
    if transform == "fft":
        return np.fft.rfft(frames, axis=1)
    return frames @ mdct_basis(n).T


def inverse(spec, n, transform, length):
    if transform == "fft":
        frames = np.fft.irfft(spec, n, axis=1)
    else:
        frames = (spec @ mdct_basis(n)) * (4 / n)
    return overlap_add(frames * sine_window(n), n, length)


# ---------------------------------------------------------------- codec

def l_bins(k):
    # odd bins (1, 3, 5, ...) belong to channel 1, even bins to channel 2
    return np.arange(k) % 2 == 1


def bin_hz(k, n, sr, transform):
    k = np.asarray(k, dtype=np.float64)
    return k * sr / n if transform == "fft" else (k + 0.5) * sr / n


GAIN_STEP_DB = 0.5


def band_energy(spec, band_width):
    # sum of |bin|^2 over consecutive groups of band_width bins
    n_frames, n_bins = spec.shape
    n_bands = -(-n_bins // band_width)
    e = np.zeros((n_frames, n_bands * band_width))
    e[:, :n_bins] = np.abs(spec) ** 2
    return e.reshape(n_frames, n_bands, band_width).sum(axis=2)


def encode(x, sr, transform="fft", frame=2048, max_freq=None, band_width=4):
    if x.shape[1] != 2:
        raise ValueError(f"need a stereo file, got {x.shape[1]} channel(s)")
    if frame % 2:
        raise ValueError("--frame must be even")

    L = forward(x[:, 0], frame, transform)
    R = forward(x[:, 1], frame, transform)

    n_bins = L.shape[1]
    if max_freq is not None:
        n_bins = int(np.sum(bin_hz(np.arange(n_bins), frame, sr, transform) <= max_freq))

    L, R = L[:, :n_bins], R[:, :n_bins]
    mask = l_bins(n_bins)
    packed = np.where(mask, L, R)

    # per-band L/R level difference, in GAIN_STEP_DB steps
    eps = 1e-20
    ratio_db = 10 * np.log10((band_energy(L, band_width) + eps) / (band_energy(R, band_width) + eps))
    lr_gain = np.clip(np.round(ratio_db / GAIN_STEP_DB), -127, 127).astype(np.int8)

    if transform == "fft":
        packed = np.stack([packed.real, packed.imag], axis=-1)

    return {
        "coeffs": packed.astype(np.float16),
        "lr_gain": lr_gain,
        "band_width": band_width,
        "sr": sr,
        "length": len(x),
        "frame": frame,
        "transform": transform,
    }


def box_mean(a, width):
    # moving average along the bin axis, same length as input
    pad = width // 2
    p = np.pad(a, ((0, 0), (pad, width - pad)), mode="edge")
    c = np.cumsum(p, axis=1)
    c = np.concatenate([np.zeros((a.shape[0], 1)), c], axis=1)
    return (c[:, width:] - c[:, :-width])[:, :a.shape[1]] / width


def rebuild_avg(own, kept):
    out = own.copy()
    n = own.shape[1]
    for i in np.where(~kept)[0]:
        nb = [j for j in (i - 1, i + 1) if 0 <= j < n]
        out[:, i] = np.mean(own[:, nb], axis=1)
    return out


def rebuild_cross(own, kept, other, gain):
    return np.where(kept, own, gain * other)


def rebuild_blind(own, kept, other, width, max_gain=10.0):
    # local energy of each channel, measured only from the bins it kept
    e_own = box_mean(np.abs(own) ** 2 * kept, width) / np.maximum(box_mean(kept[None, :] * 1.0, width), 1e-9)
    e_oth = box_mean(np.abs(other) ** 2 * ~kept, width) / np.maximum(box_mean(~kept[None, :] * 1.0, width), 1e-9)
    eps = 1e-12 * (np.max(e_own) + np.max(e_oth) + 1e-30)
    gain = np.minimum(np.sqrt((e_own + eps) / (e_oth + eps)), max_gain)
    return np.where(kept, own, gain * other)


def decode(enc, rebuild="cross", width=8):
    coeffs = enc["coeffs"].astype(np.float64)
    transform = str(enc["transform"])
    frame = int(enc["frame"])
    length = int(enc["length"])

    if transform == "fft":
        coeffs = coeffs[..., 0] + 1j * coeffs[..., 1]
        full_bins = frame // 2 + 1
    else:
        full_bins = frame // 2

    n_bins = coeffs.shape[1]
    mask = l_bins(n_bins)

    # split the packed spectrum back into two half-empty channels
    L = np.where(mask, coeffs, 0)
    R = np.where(~mask, coeffs, 0)

    if rebuild == "avg":
        L, R = rebuild_avg(L, mask), rebuild_avg(R, ~mask)
    elif rebuild == "blind":
        L, R = rebuild_blind(L, mask, R, width), rebuild_blind(R, ~mask, L, width)
    else:
        bw = int(enc["band_width"])
        g = 10 ** (enc["lr_gain"].astype(np.float64) * GAIN_STEP_DB / 20)
        g = np.repeat(g, bw, axis=1)[:, :n_bins]   # amplitude ratio L/R per bin
        L, R = rebuild_cross(L, mask, R, g), rebuild_cross(R, ~mask, L, 1 / g)

    pad = ((0, 0), (0, full_bins - n_bins))
    L, R = np.pad(L, pad), np.pad(R, pad)

    y = np.stack([inverse(L, frame, transform, length),
                  inverse(R, frame, transform, length)], axis=1)
    return y, int(enc["sr"])


def save_encoded(enc, path_or_file):
    if isinstance(path_or_file, str):
        with open(path_or_file, "wb") as f:
            np.savez_compressed(f, **enc)
    else:
        np.savez_compressed(path_or_file, **enc)


def load_encoded(path):
    with np.load(path) as z:
        return {k: z[k] for k in z.files}


# ---------------------------------------------------------------- test

def test_signal(sr=44100, seconds=5.0, seed=0):
    # panned harmonics + plucks + a bit of uncorrelated noise
    rng = np.random.default_rng(seed)
    t = np.arange(int(sr * seconds)) / sr
    L = np.zeros_like(t)
    R = np.zeros_like(t)
    for f in [110.3, 220.7, 330.1, 441.5, 987.2, 2500.9, 5200.4]:
        pan = rng.uniform(0.2, 0.8)
        s = np.sin(2 * np.pi * f * t + rng.uniform(0, 2 * np.pi)) * rng.uniform(0.05, 0.2)
        L += s * np.sqrt(1 - pan)
        R += s * np.sqrt(pan)
    for start in np.arange(0.5, seconds, 1.0):
        env = np.where(t >= start, np.exp(-(t - start) * 6), 0)
        pluck = env * np.sin(2 * np.pi * 659.3 * t) * 0.3
        L += pluck * 0.8
        R += pluck * 0.4
    L += 0.003 * rng.standard_normal(len(t))
    R += 0.003 * rng.standard_normal(len(t))
    x = np.stack([L, R], axis=1)
    return x / np.max(np.abs(x)) * 0.9, sr


def snr_db(ref, y):
    return 10 * np.log10(np.sum(ref ** 2) / max(np.sum((ref - y) ** 2), 1e-30))


def run_test(args):
    if args.input:
        x, sr = read_wav(args.input)
        src = args.input
    else:
        x, sr = test_signal()
        src = "synthetic test signal"

    enc = encode(x, sr, args.transform, args.frame, args.max_freq, args.band_width)
    buf = io.BytesIO()
    save_encoded(enc, buf)
    pcm_bytes = x.shape[0] * x.shape[1] * 2

    print(f"source      : {src} ({x.shape[0] / sr:.2f} s, {sr} Hz)")
    print(f"transform   : {args.transform}, frame {args.frame}, "
          f"bin width {sr / args.frame:.2f} Hz, {enc['coeffs'].shape[1]} bins kept")
    print(f"size        : PCM16 stereo {pcm_bytes:,} B -> encoded {buf.tell():,} B "
          f"({pcm_bytes / buf.tell():.2f}:1)")
    print()

    # what you would get if the decoder had every bin (float16 + transform limit)
    ref = decode_full(x, sr, args)
    mid = np.repeat(x.mean(axis=1, keepdims=True), 2, axis=1)
    rows = [("all bins (no interleave)", ref),
            ("mono mid (L+R)/2", mid),
            ("rebuild avg", decode(enc, "avg")[0]),
            (f"rebuild blind (width {args.width})", decode(enc, "blind", args.width)[0]),
            (f"rebuild cross (band {args.band_width})", decode(enc, "cross")[0])]
    print(f"{'method':30s} {'L SNR':>9s} {'R SNR':>9s}")
    for name, y in rows:
        print(f"{name:30s} {snr_db(x[:, 0], y[:, 0]):7.2f}dB {snr_db(x[:, 1], y[:, 1]):7.2f}dB")


def decode_full(x, sr, args):
    out = []
    for c in range(2):
        s = forward(x[:, c], args.frame, args.transform)
        if args.max_freq is not None:
            s[:, bin_hz(np.arange(s.shape[1]), args.frame, sr, args.transform) > args.max_freq] = 0
        if args.transform == "fft":
            s = s.real.astype(np.float16) + 1j * s.imag.astype(np.float16)
        else:
            s = s.astype(np.float16).astype(np.float64)
        out.append(inverse(s, args.frame, args.transform, len(x)))
    return np.stack(out, axis=1)


# ---------------------------------------------------------------- CLI

def main():
    ap = argparse.ArgumentParser(description="Stereo FFT/MDCT bin-interleave codec")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def enc_opts(p):
        p.add_argument("--transform", choices=["fft", "mdct"], default="fft")
        p.add_argument("--frame", type=int, default=2048, help="frame size in samples (even)")
        p.add_argument("--max-freq", type=float, default=None, help="drop bins above this frequency (Hz)")
        p.add_argument("--band-width", type=int, default=4, help="bins per stored L/R gain band")

    def dec_opts(p):
        p.add_argument("--rebuild", choices=["cross", "blind", "avg"], default="cross")
        p.add_argument("--width", type=int, default=8, help="bins used to guess L/R ratio (blind)")

    p = sub.add_parser("encode", help="stereo WAV -> .npz")
    p.add_argument("input")
    p.add_argument("output")
    enc_opts(p)

    p = sub.add_parser("decode", help=".npz -> stereo WAV")
    p.add_argument("input")
    p.add_argument("output")
    dec_opts(p)

    p = sub.add_parser("test", help="round-trip a WAV (or a synthetic signal) and print SNR")
    p.add_argument("input", nargs="?")
    enc_opts(p)
    dec_opts(p)

    args = ap.parse_args()

    if args.cmd == "encode":
        x, sr = read_wav(args.input)
        enc = encode(x, sr, args.transform, args.frame, args.max_freq, args.band_width)
        save_encoded(enc, args.output)
        print(f"{args.input} -> {args.output}: {os.path.getsize(args.input):,} B -> "
              f"{os.path.getsize(args.output):,} B")
    elif args.cmd == "decode":
        y, sr = decode(load_encoded(args.input), args.rebuild, args.width)
        peak = np.max(np.abs(y))
        if peak > 1:
            print(f"warning: output peaks at {peak:.2f}, clipping", file=sys.stderr)
        write_wav(args.output, y, sr)
        print(f"{args.input} -> {args.output}")
    else:
        run_test(args)


if __name__ == "__main__":
    main()
