#!/usr/bin/env python3
"""Rate-distortion comparison of TJC against older codecs (and H.264 as a reference).
Runs on Windows, Linux and macOS: needs Python 3, ffmpeg in PATH and the tjc tool.

1. Cut a test clip. This writes clip.yuv (the reference every codec is measured
   against) and clip.avi (the same frames as uncompressed RGB, for encoders that
   only take ordinary video files, such as RAD Video Tools for Bink):

       python bench/compare_codecs.py prepare test.mp4 --seconds 4

2. Optionally encode clip.avi with an outside encoder at several bitrates, e.g.
   Bink: bink_1000.bik, bink_2000.bik, ...

3. Measure everything:

       python bench/compare_codecs.py run --tjc ./tjc --extra bink=bink_*.bik

   One line per encode goes to results.txt (codec setting bytes kbit/s psnr_y
   psnr_all encode_seconds), followed by a BD-rate table (bench/bd_rate.py).

The older codecs get their better encoder options (RD macroblock decisions,
trellis quantization, 4 vectors per macroblock where allowed), no B-frames (TJC
has none) and a keyframe every 4 seconds. H.261 and baseline H.263 only allow
CIF-family sizes and are included only for those. Everything runs on one thread.
"""
import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bd_rate  # noqa: E402


def ffmpeg(*args, capture=False):
    cmd = ['ffmpeg', '-nostdin', '-hide_banner', '-y'] + [str(a) for a in args]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0 and not capture:
        raise RuntimeError('ffmpeg failed: ' + ' '.join(cmd) + '\n' + r.stderr[-2000:])
    return r.stderr


def probe(path):
    out = subprocess.run(['ffprobe', '-v', 'error', '-select_streams', 'v:0', '-show_entries',
                          'stream=width,height,r_frame_rate', '-of', 'json', path],
                         stdout=subprocess.PIPE, text=True, check=True).stdout
    s = json.loads(out)['streams'][0]
    num, den = s['r_frame_rate'].split('/')
    return int(s['width']), int(s['height']), float(num) / float(den)


# ---------------------------------------------------------------------------
# prepare
# ---------------------------------------------------------------------------
def prepare(a):
    w, h, fps = probe(a.input)
    if a.size:
        w, h = (int(v) for v in a.size.lower().split('x'))
    w -= w % 2
    h -= h % 2
    vf = f'fps={fps},scale={w}:{h}:flags=bicubic'
    # RGB first, then the YUV reference made from that RGB, so a codec fed the AVI
    # starts from exactly the pixels the others get (only its own RGB->YUV step differs).
    ffmpeg('-ss', a.start, '-i', a.input, '-t', a.seconds, '-an', '-vf', vf, '-c:v', 'rawvideo',
           '-pix_fmt', 'bgr24', a.avi)
    ffmpeg('-i', a.avi, '-f', 'rawvideo', '-pix_fmt', 'yuv420p', a.yuv)
    frames = os.path.getsize(a.yuv) // (w * h * 3 // 2)
    info = {'width': w, 'height': h, 'fps': fps, 'frames': frames, 'yuv': a.yuv, 'avi': a.avi}
    with open(a.info, 'w') as f:
        json.dump(info, f, indent=1)
    print(f'{a.yuv}: {w}x{h}, {frames} frames at {fps:g} fps (reference)')
    print(f'{a.avi}: the same frames as uncompressed RGB, for outside encoders')
    print(f'settings saved in {a.info}')


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------
class Bench:
    def __init__(self, a, info, tmp):
        self.a = a
        self.w, self.h, self.fps, self.n = info['width'], info['height'], info['fps'], info['frames']
        self.ref = info['yuv']
        self.info = info
        self.tmp = tmp
        self.secs = self.n / self.fps
        self.raw = ['-f', 'rawvideo', '-pix_fmt', 'yuv420p', '-s', f'{self.w}x{self.h}', '-r', f'{self.fps:g}']
        self.lines = []

    def psnr(self, decoded, ref=None):
        log = ffmpeg(*self.raw, '-i', decoded, *self.raw, '-i', ref or self.ref, '-lavfi', '[0:v][1:v]psnr',
                     '-f', 'null', '-', capture=True)
        for line in log.splitlines():
            if 'PSNR y:' in line:
                py = line.split('PSNR y:')[1].split()[0]
                pa = line.split('average:')[1].split()[0]
                return float(py if py != 'inf' else 99), float(pa if pa != 'inf' else 99)
        raise RuntimeError('psnr failed:\n' + log[-1000:])

    # How an outside encoder may have turned the RGB clip into YUV. For its files
    # the reference is rebuilt from the RGB clip the same way, so a different
    # matrix or range costs it nothing.
    COLOUR = {
        'BT.601': ('bt601', 'tv'),
        'BT.601 full range': ('bt601', 'pc'),
        'BT.709': ('bt709', 'tv'),
        'BT.709 full range': ('bt709', 'pc'),
    }

    def reference(self, name):
        if name == 'BT.601':
            return self.info['yuv']
        matrix, rng = self.COLOUR[name]
        out = os.path.join(self.tmp, 'ref_' + matrix + rng + '.yuv')
        if not os.path.exists(out):
            ffmpeg('-i', self.info['avi'], '-frames:v', self.n, '-vf',
                   f'scale=out_color_matrix={matrix}:out_range={rng}', '-f', 'rawvideo', '-pix_fmt', 'yuv420p', out)
        return out

    def decode_ffmpeg(self, path, fmt=None, keep_range=False):
        out = os.path.join(self.tmp, 'dec.yuv')
        args = ['-threads', '1'] + (['-f', fmt] if fmt else []) + ['-i', path, '-frames:v', self.n]
        if keep_range:  # leave full-range (yuvj) output as it is
            args += ['-vf', 'scale=in_range=pc:out_range=pc']
        ffmpeg(*args, '-f', 'rawvideo', '-pix_fmt', 'yuv420p', out)
        return out

    def record(self, codec, setting, path, seconds, py, pa):
        nbytes = os.path.getsize(path)
        line = f'{codec} {setting} {nbytes} {nbytes * 8 / self.secs / 1000:.1f} {py:.4f} {pa:.4f} {seconds:.3f}'
        print(line, flush=True)
        self.lines.append(line)

    def run_ffmpeg(self, codec, setting, args, fmt):
        out = os.path.join(self.tmp, 'out')
        t0 = time.time()
        ffmpeg(*self.raw, '-i', self.ref, *args, '-threads', 1, '-f', fmt, out)
        t = time.time() - t0
        py, pa = self.psnr(self.decode_ffmpeg(out, 'mjpeg' if codec == 'mjpeg' else None))
        self.record(codec, setting, out, t, py, pa)

    def run_tjc(self, codec, q, opts):
        out = os.path.join(self.tmp, 'out.tjc')
        dec = os.path.join(self.tmp, 'dec.yuv')
        t0 = time.time()
        subprocess.run([self.a.tjc, 'encode', '--size', f'{self.w}x{self.h}', '--fps', f'{self.fps:g}',
                        '--quality', str(q), '--threads', '1', *opts, '-i', self.ref, '-o', out, '-q'], check=True)
        t = time.time() - t0
        subprocess.run([self.a.tjc, 'decode', '-i', out, '-o', dec, '-q'], check=True)
        py, pa = self.psnr(dec)
        self.record(codec, q, out, t, py, pa)

    def run_extra(self, label, pattern):
        files = sorted(glob.glob(pattern))
        if not files:
            print(f'# {label}: no files match {pattern}', flush=True)
            return
        # Encoders differ in how they turn RGB into YUV (matrix, range). Try the
        # common choices on the largest file and keep the best match.
        big = max(files, key=os.path.getsize)
        best_name, best = None, None
        try:
            for name, (_, rng) in self.COLOUR.items():
                py = self.psnr(self.decode_ffmpeg(big, keep_range=rng == 'pc'), self.reference(name))[0]
                if best is None or py > best + 0.2:
                    best_name, best = name, py
        except RuntimeError as e:
            print(f'# {label}: {e}', flush=True)
            return
        ref, full = self.reference(best_name), self.COLOUR[best_name][1] == 'pc'
        print(f'# {label}: its YUV matches {best_name}; measured against a reference made that way', flush=True)
        for f in files:
            setting = os.path.splitext(os.path.basename(f))[0]
            py, pa = self.psnr(self.decode_ffmpeg(f, keep_range=full), ref)
            self.record(label, setting, f, 0.0, py, pa)


def run(a):
    with open(a.info) as f:
        info = json.load(f)
    if not shutil.which(a.tjc) and not os.path.exists(a.tjc):
        sys.exit(f'tjc tool not found: {a.tjc} (use --tjc)')
    encoders = subprocess.run(['ffmpeg', '-hide_banner', '-encoders'], stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL, text=True).stdout
    with tempfile.TemporaryDirectory() as tmp:
        b = Bench(a, info, tmp)
        gop = int(round(b.fps * 4))
        old = ['-mbd', 'rd', '-trellis', 1, '-cmp', 2, '-subcmp', 2, '-g', gop, '-bf', 0]
        cif = (b.w, b.h) in [(128, 96), (176, 144), (352, 288), (704, 576), (1408, 1152)]
        want = set(a.codecs.split(',')) if a.codecs else None
        use = lambda c: want is None or c in want  # noqa: E731
        qs = [2, 3, 4, 6, 9, 13, 20, 31]
        if use('mjpeg'):
            for q in qs[:-1]:
                b.run_ffmpeg('mjpeg', q, ['-c:v', 'mjpeg', '-strict', -1, '-q:v', q], 'mjpeg')
        for q in qs:
            if cif and use('h261'):
                b.run_ffmpeg('h261', q, ['-c:v', 'h261', '-q:v', q, '-g', gop, '-mbd', 'rd'], 'h261')
            if cif and use('h263'):
                b.run_ffmpeg('h263', q, ['-c:v', 'h263', '-q:v', q, *old, '-flags', '+mv4'], 'h263')
            if use('h263p'):
                b.run_ffmpeg('h263p', q, ['-c:v', 'h263p', '-q:v', q, *old, '-flags', '+mv4+aic', '-umv', 1,
                                          '-aiv', 1], 'h263')
            if use('mpeg1'):
                b.run_ffmpeg('mpeg1', q, ['-c:v', 'mpeg1video', '-q:v', q, *old], 'mpeg1video')
            if use('mpeg2'):
                b.run_ffmpeg('mpeg2', q, ['-c:v', 'mpeg2video', '-q:v', q, *old], 'mpeg2video')
            if use('mpeg4'):
                b.run_ffmpeg('mpeg4', q, ['-c:v', 'mpeg4', '-q:v', q, *old, '-flags', '+mv4'], 'm4v')
        if use('x264') and 'libx264' in encoders:
            for c in [16, 20, 24, 28, 32, 36, 40]:
                b.run_ffmpeg('x264', c, ['-c:v', 'libx264', '-crf', c, '-preset', 'medium', '-tune', 'psnr'], 'h264')
        for q in [30, 45, 60, 75, 85, 92, 97]:
            if use('tjc2'):
                b.run_tjc('tjc2', q, ['--format', 'tjc2'])
            if use('tjc3'):
                b.run_tjc('tjc3', q, ['--format', 'tjc3'])
            if use('tjc4'):
                b.run_tjc('tjc4', q, ['--quadtree'])
        for e in a.extra or []:
            label, _, pattern = e.partition('=')
            b.run_extra(label, pattern)
    with open(a.out, 'w') as f:
        f.write('\n'.join(b.lines) + '\n')
    print(f'\nresults in {a.out}\n')
    present = {line.split()[0] for line in b.lines}
    refs = [r for r in dict.fromkeys([a.ref, 'h263p', 'h263', 'mpeg4', 'tjc3']) if r in present][:2]
    for ref in refs:
        bd_rate.report([a.out], ref, a.psnr_target)
        print()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('prepare', help='cut a test clip from a video')
    p.add_argument('input')
    p.add_argument('--start', default='0', help='start time (default 0)')
    p.add_argument('--seconds', default='4', help='length (default 4)')
    p.add_argument('--size', help='scale to WxH (default: source size)')
    p.add_argument('--yuv', default='clip.yuv')
    p.add_argument('--avi', default='clip.avi')
    p.add_argument('--info', default='clip.json')
    r = sub.add_parser('run', help='encode, decode and measure')
    r.add_argument('--info', default='clip.json', help='written by prepare')
    r.add_argument('--tjc', default='tjc.exe' if os.name == 'nt' else './tjc')
    r.add_argument('--codecs', help='comma list (default all): mjpeg,h261,h263,h263p,mpeg1,mpeg2,mpeg4,x264,'
                                    'tjc2,tjc3,tjc4')
    r.add_argument('--extra', action='append',
                   help='LABEL=GLOB: files made by another encoder from clip.avi, e.g. bink=bink_*.bik')
    r.add_argument('--out', default='results.txt')
    r.add_argument('--ref', default='mjpeg', help='BD-rate reference codec')
    r.add_argument('--psnr-target', type=float, default=38.0)
    a = ap.parse_args()
    prepare(a) if a.cmd == 'prepare' else run(a)


if __name__ == '__main__':
    main()
