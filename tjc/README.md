# TJC — Tiled-JPEG Codec

A small video codec for low-power decoders (written with a single-core MIPS 24KEc in
mind: integer-only, no FPU needed). Each frame is split into tiles; only tiles that
changed (plus tiles picked by a refresh policy) are sent, each as a tiny
baseline-JPEG-style blob or, in TJC3, as a motion-compensated copy of the previous
frame plus a small correction. TJC4 adds variable tile size (a quadtree in each
tile) and a quality per frame, which enables rate control (average target and a
hard cap). The decoder keeps a persistent framebuffer and patches it. An optional audio track is coded with [QOA](https://qoaformat.org) and
interleaved per video frame, so audio and video stay in sync by construction.

## Compression tools (TJC3)

| Tool | What it does | Format | Decoder cost |
|------|--------------|--------|--------------|
| Motion compensation (`--motion`, default on) | A dirty tile can be "the previous frame at offset (dx, dy), half-pel precision, plus a coded correction" | TJC3 | +1 frame of RAM; a block copy and an add per moving tile |
| Adaptive Huffman tables (`--huffman`, default adaptive) | Entropy tables fitted to the video, sent only when they pay off | TJC3 | builds a table when one arrives (~2 KB) |
| Deadzone quantization (`--deadzone`, default 0.33) | Small coefficients round to zero | any | none |
| Smart rolling refresh (`--smart-refresh`, default on) | Refresh skips tiles whose content is recent | any | none |
| Skip invisible updates (`--skip-invisible T`, default off) | Drops updates nobody could see (grain) | any | none; static areas keep their last grain |

Bytes needed for the same picture quality (BD-rate on PSNR), compared with the
previous release, on 720p test clips:

| Clip | `--format tjc2` | **TJC3** | TJC3 + `--skip-invisible 2` |
|------|-----------------|----------|-----------------------------|
| Screen content | -3.7% | **-65%** | -65% |
| Motion graphics | -5.0% | **-33%** | -33% |
| Noisy camera | -8.5% | **-35%** | **-51%** |
| Zoom | -3.4% | **-40%** | -40% |
| Camera pan | -3.2% | **-83%** | -83% |
| Noisy pan | -6.2% | **-71%** | **-76%** |

`--format tjc2` turns motion compensation and adaptive tables off: the file then
plays on TJC2 decoders. New decoders play TJC1, TJC2 and TJC3.

Decoding TJC3 is about as fast as TJC2 (faster on high-motion clips, since there are
fewer bits to unpack; ~25% slower in the worst case measured). A size-optimized
decoder-only build (`-Os`, `TJC_NO_ENCODER`) is about 21 KB of x86-64 code including
motion compensation and QOA audio.

Encoding TJC3 costs more than TJC2: the motion search. `--effort` picks the trade-off
(it never changes the format); at 720p on 4 threads:

| Effort | Noisy clip | Size vs `best` |
|--------|------------|----------------|
| `fast` | 42 fps | up to +9% |
| `normal` (default) | 37 fps | within ~0.5% (2% on motion graphics) |
| `best` | 23 fps | — |

## TJC4: variable tile size and rate control

### Variable tile size (`--quadtree`)

Each tile, 64x64 by default, is a quadtree. A large area that moves together costs
one vector. A small change inside a big tile sends only a small leaf, down to 8x8,
and the rest of the tile is marked "skip". The encoder decides split versus no split
by rate-distortion cost, so it only splits where splitting pays.

Same picture quality, bytes compared with TJC3 at 16x16 (q75, 1 thread):

| Clip | TJC3 | **TJC4 `--quadtree`** (64 → 8) | 64 → 16 (`--tile-size 64 --split-levels 2`) | Encode time vs TJC3 |
|------|------|-------------------------------|------------------------------|---------------------|
| Screen content 720p | 199 KB | **95 KB (-52%)** | 102 KB (-49%) | 2.5x |
| Camera pan 720p | 1.20 MB | **729 KB (-39%)** | 734 KB (-39%) | 3.0x |
| Noisy camera 1080p | 5.28 MB | **3.76 MB (-29%)** | 3.79 MB (-28%) | 1.7x |
| Noisy camera 720p | 3.65 MB | **2.93 MB (-20%)** | 2.97 MB (-19%) | 1.7x |
| Camera 1080p | 1.45 MB | **1.18 MB (-18%)**, +0.8 dB | 1.21 MB (-16%) | 3.1x |
| Noisy pan 720p | 2.95 MB | 2.41 MB (-18%) | **2.37 MB (-20%)** | 3.2x / 1.8x |
| Motion graphics 720p | 1.07 MB | **889 KB (-17%)** | 914 KB (-14%) | 3.3x |
| Zoom 720p | 1.70 MB | 1.55 MB (-9%) | **1.50 MB (-12%)** | 4.7x / 3.6x |
| Thin moving edges 720p | 335 KB | **238 KB (-29%)** | 243 KB (-28%) | 3.6x |

The decoder does the same work per pixel, and decoding is 3-15% slower than TJC3
because of more, smaller blocks. It needs no extra frame memory, only a 6-byte
entry per smallest leaf (86 KB at 720p with 8x8 leaves, 22 KB with 16x16) that
tells a player how each area was coded. The decoder-only build grows from 21 KB to
24 KB. So TJC4 still fits the MIPS target; the extra cost is on the encoder side.
`--split-levels 0` gives fixed tiles with the TJC4 syntax, which is what you get
when you only want rate control.

### Rate control (`--bitrate`, `--maxrate`)

- `--bitrate 2M`: average target, audio included. The encoder predicts each frame's
  size from a rate model (bits ∝ complexity × quantizer step^-α) and chooses a
  quality between `--min-quality` and `--quality`. When a frame lands far off its
  target, the encoder codes it again at a corrected quality. Overspending and
  underspending are paid back over about a second. Quality rises at most 3 steps per
  frame, so it doesn't pump.
- `--maxrate 3M --bufsize 1000`: a hard cap. No window of `--bufsize` ms (in whole
  frames) carries more than the max rate. It is a strict sliding window, not a
  leaky bucket: a leaky bucket lets the first second burst to twice the cap. The
  encoder plans the next half-window evenly, so a burst cannot starve the frames
  after it. If a frame still doesn't fit, it is re-coded lower, down to quality 1,
  and then the least useful dirty tiles wait for the next frame. Refresh-only tiles
  go first, then the tiles that improve the picture least per bit. Only a full
  refresh frame can break the cap, and only in TJC2/TJC3.
- Both together: target the average, never exceed the max.
- When even `--min-quality` is too big for the target, the encoder raises the
  threshold, deadzone and skip-invisible in steps. With `--format tjc2`/`tjc3` there
  is no per-frame quality, so those steps (plus holding tiles back under a cap) are
  the only levers. TJC3 with `--maxrate` keeps to the cap by sending updates evenly
  over time, which shows up as lower temporal resolution on heavy content.

Measured on 720p30 clips, `--quadtree`:

| Clip | Settings | Average | Busiest 1 s | Quality used |
|------|----------|---------|-------------|--------------|
| Noisy camera | `--bitrate 2000` | 1992 kbit/s | 2011 kbit/s | 10..30 |
| Camera pan | `--bitrate 2000 --maxrate 3000` | 1770 kbit/s | 1920 kbit/s | 11..57 |
| Camera pan | `--maxrate 1500` | 1405 kbit/s | 1476 kbit/s | 10..75 |
| Noisy camera | `--maxrate 2500` | 2443 kbit/s | 2459 kbit/s | 15..75 |
| Noisy camera 1080p | `--bitrate 8000` | 8044 kbit/s | 8105 kbit/s | 32..75 |

### How high can the bitrate go?

The theoretical worst case per 8x8 block is 27 + 63×26 bits. In practice, the
heaviest content possible (full-range random noise, every tile, every frame) at
quality 100 gives:

| | 720p30 | 1080p30 (scaled) |
|---|---|---|
| Random noise, q100 | ~315 Mbit/s | ~710 Mbit/s |
| Random noise, q75 | ~124 Mbit/s | ~280 Mbit/s |
| Real camera clip, q100 | ~89 Mbit/s | ~200 Mbit/s |
| Raw YUV420P, for comparison | 332 Mbit/s | 746 Mbit/s |

So the codec tops out near the raw rate. Use `--maxrate` to keep a link or an SD
card within its budget.

### Ghost specks (`--peak-threshold`, all formats)

A tile is "dirty" when the mean difference in one of its 8x8 blocks exceeds
`--motion-threshold`. A thin edge that moves away leaves 1-2 pixels behind, with
too little mean difference to trigger. Motion compensation then dragged such specks
along with the edge. Now a tile is also re-sent when any single pixel changed by
more than `--peak-threshold` (default 32). In addition, the encoder rejects a motion
compensated tile if it would leave a pixel that far off and intra coding gets it
clearly closer. On a moving orange edge over a dark panel, pixels more than 40 levels
off dropped from 3476 to 93, and the file got 4% smaller. Real video costs about
nothing extra. This is an encoder change only: re-encoded TJC2/TJC3 files play on the
old decoders.

## How TJC compares with older codecs

These are rate-distortion sweeps on six 60-frame clips: camera pan, noisy camera,
zoom, screen content, motion graphics and noisy pan. All codecs run single-threaded.
The older codecs get their better ffmpeg options: RD macroblock decisions,
trellis, 4 vectors per macroblock, no B-frames and a keyframe every 4 s. H.264 runs
`-tune psnr`. H.261 and baseline H.263 only take CIF-family sizes, so they are
measured on a CIF version of the same clips. The scripts are in `bench/`, so you
can run your own clip:

```bash
ffmpeg -i test.mp4 -t 4 -f rawvideo -pix_fmt yuv420p clip.yuv
TJC=./tjc bench/compare_codecs.sh clip.yuv 1280 720 120 30 > results.txt
python3 bench/bd_rate.py results.txt --ref h263p
```

BD-rate (PSNR-Y), the average bitrate difference at equal quality. Negative means
fewer bits:

| Codec (year) | 720p vs MJPEG | 720p vs H.263+ | CIF vs MJPEG | CIF vs H.263 |
|---|---|---|---|---|
| MJPEG (1992) | 0% | +509% | 0% | +559% |
| H.261 (1990) | — | — | -68% | +56% |
| H.263 (1996) | — | — | -80% | 0% |
| H.263+ (1998) | -72% | 0% | — | — |
| MPEG-1 (1993) | -72% | -1% | — | — |
| MPEG-2 (1995) | -71% | +6% | — | — |
| MPEG-4 Part 2 (1999) | -73% | -4% | — | — |
| H.264 / x264 (2003, reference) | -87% | -60% | -94% | -67% |
| TJC2 | +5% | +447% | -27% | +263% |
| TJC3 | -45% | +115% | -60% | +101% |
| **TJC4 `--quadtree`** | **-55%** | **+61%** | **-63%** | **+78%** |

Per clip, TJC4 against H.263+ at 720p: screen content +5%, motion graphics +1%,
zoom +60%, camera pan +49%, noisy camera +77%, noisy pan +173%. At high quality
(above ~43 dB) TJC4 overtakes H.263+/MPEG-4 on the noisy clips and zoom, where the
older codecs' curves flatten. At CIF, TJC4 ends up between H.261 and H.263: 26%
behind H.261 on average, 10% ahead of it on zoom.

What this means:

- **TJC beats MJPEG by a mile** (TJC3/TJC4 need about half the bits), with a
  decoder about as simple.
- **TJC is not yet an H.263 killer.** At low and medium bitrates on camera video,
  1990s inter codecs still need 30-60% fewer bits. It catches up on screen content
  and motion graphics, and at high quality.
- Where the bits go, measured on the camera pan:
  - The per-tile 2-byte length, byte padding and dirty bitmap are roughly a
    quarter of the stream when most tiles change.
  - About 12% goes to motion-compensated residual blocks that are all zero. H.263
    marks those with one coded-block-pattern code per macroblock.
  - The JPEG perceptual quant tables also cost PSNR against the flat matrices of
    H.263/MPEG. That is a measurement artifact as much as a real loss.
- **Decode speed** (x86, single thread, 720p at 38 dB): TJC3/TJC4 decode motion
  graphics at 780-920 fps, faster than ffmpeg's SIMD-optimized MPEG-1/2/4 and H.263+
  decoders (550-630 fps). The camera pan decodes at about 300 fps against 500-600,
  and the noisy camera at about 130 fps against 340-375. There TJC also unpacks 2-3x
  as many bits for the same quality. TJC's decoder is plain portable C++ with no
  SIMD, which is the relevant case for the MIPS target.
- TJC's design advantages are elsewhere:
  - Partial-screen updates with no macroblock overhead for static areas.
  - Exact seeking without keyframes.
  - Per-tile error isolation.
  - Integer-exact output on every platform.
  - Built-in synchronized audio.
  - A 24 KB decoder.

| File          | What                                                              |
|---------------|-------------------------------------------------------------------|
| `tjc.h`       | single-header library: encoder + decoder (stb-style)              |
| `tjc.cpp`     | `tjc` command line tool: `encode`, `decode`, `info`               |
| `tjc_test.cpp`| unit tests (DCT, Huffman, motion, seeking, QOA, fuzzing, ...) + end-to-end |
| `gui/`        | **TJC Studio**: desktop player + encoder for Windows and Linux/X11 |
| `bench/`      | rate-distortion comparison against MJPEG/H.261/H.263/MPEG-1/2/4/H.264 |

## Build

```bash
g++ -std=c++17 -O2 -pthread -o tjc tjc.cpp
g++ -std=c++17 -O2 -pthread -o tjc_test tjc_test.cpp && ./tjc_test
```

or `cmake -B build && cmake --build build && ./build/tjc_test`.

Optional: `ln -s tjc tjc_encode && ln -s tjc tjc_decode` lets you call the tool
without the subcommand.

## TJC Studio (GUI)

A desktop app for playing and encoding TJC streams. One codebase for Windows and
Linux/X11, built with GLFW, Dear ImGui (OpenGL 2) and miniaudio. CMake downloads
all three when you configure the build.

**Player tab**
- Plays video and audio in sync, timed by the audio clock (or the wall clock when
  there is no audio or no sound device).
- Exact seeking: slider, Left/Right = 5 s, `,` `.` = one frame, Home = start.
  An index pass computes for each frame the oldest frame whose tiles are still on
  screen. Seeking decodes from there, so the picture is bit-identical to playing
  from the start without decoding the whole file. That's one refresh cycle for
  rolling streams, back to the last full refresh for periodic ones.
- Tile overlay (`T`): red = tile coded from scratch, green = motion compensated,
  blue = full refresh.
- Info panel: stream settings, per-frame tiles/bytes/audio/quality, the frame each
  seek restarts from, a frame-size graph, and bitrate. Bitrate shows now (last
  second), average and peak (busiest second), plus a bitrate graph and a per-frame
  quality graph for TJC4.
- With TJC4 quadtree streams, the overlay shows the leaves themselves; unchanged
  leaves stay clear.
- Loop, volume, BT.601/BT.709 color matrix, click the picture to play/pause.

**Encoder tab**
- Input: any file ffmpeg can read. ffprobe fills in size, frame rate and audio.
- Output size presets (source/1080p/720p/480p/360p/custom), frame rate, tile size,
  quality, motion threshold, refresh mode, keyframes, threads, audio on/off with
  resampling (rate, stereo/mono).
- Compression: format (auto/TJC2/TJC3/TJC4), motion compensation, adaptive Huffman
  tables, variable tile size (quadtree, smallest leaf), effort, skip-invisible and
  peak threshold.
- Rate control: constant quality, average bitrate, average + max, or max only,
  with buffer and minimum quality. The progress line shows the last second's
  bitrate and the current quality.
- Progress, speed, ETA, bitrate, a live preview of what the decoder will show,
  Cancel (deletes the partial file), and "Play output".
- Produces the same bytes as `tjc encode` with the same settings.

Drop a `.tjc` on the window to play it, or any video to encode it. A file given on
the command line works the same way (`tjc_studio clip.tjc`).

ffmpeg/ffprobe must be in PATH, or set their folder under "ffmpeg location" in the
Encoder tab. The player doesn't need them.

### Build on Windows

```bat
:: Visual Studio (Developer Command Prompt) or MinGW; needs CMake 3.16+
cmake -S gui -B build-studio
cmake --build build-studio --config Release
:: -> build-studio\Release\tjc_studio.exe (VS) or build-studio\tjc_studio.exe (MinGW)
```

It's a normal GUI app with no console window. A MinGW build links its runtime
statically, so the `.exe` needs only system DLLs.

### Build on Linux (X11)

```bash
sudo apt install build-essential cmake libx11-dev libxrandr-dev libxinerama-dev \
                 libxcursor-dev libxi-dev libgl1-mesa-dev
cmake -S gui -B build-studio && cmake --build build-studio -j
./build-studio/tjc_studio
```

It also runs on Wayland desktops through XWayland. Sound goes through
PulseAudio/PipeWire/ALSA, picked at runtime.

Or build everything at once: `cmake -B build -DTJC_BUILD_STUDIO=ON`. Offline builds
can point `FETCHCONTENT_SOURCE_DIR_GLFW`, `..._IMGUI` and `..._MINIAUDIO` at local
copies of GLFW 3.4, Dear ImGui 1.91.9b and miniaudio 0.11.22.

## Windows (command line tool)

Build with any of these (output: `tjc.exe`):

```bat
:: Visual Studio: "x64 Native Tools Command Prompt for VS"
cl /std:c++17 /O2 /EHsc tjc.cpp

:: MinGW-w64 / MSYS2
g++ -std=c++17 -O2 -static -o tjc.exe tjc.cpp

:: CMake (either toolchain)
cmake -B build && cmake --build build --config Release
```

Get ffmpeg with `winget install Gyan.FFmpeg` (or a build from ffmpeg.org).

**Run the pipes from `cmd.exe`, not Windows PowerShell 5.1.** PowerShell 5.1 sends
pipe data through as text and corrupts binary video. It also has no `<` redirect.
PowerShell 7.4+ passes bytes between programs unchanged, but still has no `<`, so
use `-i`/`-o` there.

```bat
:: cmd.exe
ffmpeg -i input.mp4 -f rawvideo -pix_fmt yuv420p -s 320x240 - | tjc encode --size 320x240 > stream.tjc
tjc decode < stream.tjc | ffplay -f rawvideo -pix_fmt yuv420p -video_size 320x240 -
```

```powershell
# PowerShell (7.4+ for the pipe; the -i/-o form works in any version)
ffmpeg -i input.mp4 -f rawvideo -pix_fmt yuv420p -s 320x240 raw.yuv
.\tjc encode --size 320x240 -i raw.yuv -o stream.tjc
.\tjc decode -i stream.tjc -o out.yuv
ffplay -f rawvideo -pix_fmt yuv420p -video_size 320x240 out.yuv
```

With audio (cmd.exe; check the frame rate with `ffprobe input.mp4` first):

```bat
ffmpeg -i input.mp4 -vn -c:a pcm_s16le audio.wav
ffmpeg -i input.mp4 -f rawvideo -pix_fmt yuv420p - | tjc encode --size 1920x1080 --fps 29.97 --audio audio.wav -o stream.tjc
tjc decode -i stream.tjc -o out.yuv --audio-out out.wav
ffmpeg -f rawvideo -pix_fmt yuv420p -s 1920x1080 -framerate 30000/1001 -i out.yuv -i out.wav -c:v libx264 -c:a aac out.mp4
```

`copy tjc.exe tjc_encode.exe` (and `tjc_decode.exe`) gives the subcommand-free names,
since Windows has no handy symlinks.

## Using the tool with ffmpeg

Frames in and out are raw YUV420P, so the tool sits in an ffmpeg pipe.

```bash
# Encode
ffmpeg -i input.mp4 -f rawvideo -pix_fmt yuv420p -s 320x240 - \
  | ./tjc encode --width 320 --height 240 --tile-size 16x16 --motion-threshold 3 \
      --refresh-mode rolling --refresh-param 30 --quality 75 > stream.tjc

# Preview
./tjc decode < stream.tjc | ffplay -f rawvideo -pix_fmt yuv420p -video_size 320x240 -

# Back to a normal file
./tjc decode < stream.tjc \
  | ffmpeg -f rawvideo -pix_fmt yuv420p -s 320x240 -i - -c:v libx264 out.mp4

# Inspect a stream (header + per-frame dirty tiles and sizes)
./tjc info -i stream.tjc
```

The decoder reads width/height/tile size/frame rate from the stream; it needs no
options. Add `-v` for a per-frame log on stderr and `--psnr` (encoder) for quality
numbers.

### Audio

```bash
# 1. Pull the audio out as 16-bit PCM WAV (any rate, 1..8 channels)
ffmpeg -i input.mp4 -vn -c:a pcm_s16le audio.wav

# 2. Encode video + audio. --fps must match the video for the audio to stay in sync.
ffmpeg -i input.mp4 -f rawvideo -pix_fmt yuv420p - \
  | ./tjc encode --size 1920x1080 --fps 29.97 --audio audio.wav > stream.tjc

# 3. Decode both
./tjc decode -i stream.tjc -o out.yuv --audio-out out.wav

# 4. Watch it, or put it back into an mp4
ffmpeg -f rawvideo -pix_fmt yuv420p -s 1920x1080 -framerate 30000/1001 -i out.yuv \
  -i out.wav -c:v libx264 -c:a aac out.mp4
```

- `--fps` accepts `30`, `25`, `29.97`, `23.976`, `59.94` or an exact `30000/1001`.
  NTSC decimals are mapped to their exact `x000/1001` rates.
- Each video frame carries the audio that plays during it. Frame *n* ends at sample
  `floor((n+1) * rate * fps_den / fps_num)`, so 48 kHz at 29.97 fps alternates
  between 1601 and 1602 samples per frame and never drifts.
- If the audio is shorter than the video it is padded with silence; if it is longer
  the rest is dropped. The tool warns in both cases.
- To resample for the target (for example a 22.05 kHz mono DAC), let ffmpeg do it in
  step 1: `ffmpeg -i input.mp4 -vn -ar 22050 -ac 1 -c:a pcm_s16le audio.wav`.

QOA uses 3.2 bits per sample: about 280 kbit/s for 44.1 kHz stereo and 320 kbit/s
for 48 kHz stereo. It's integer-only (a 4-tap sign-sign LMS predictor, about 20
integer ops per sample to decode), so it suits a CPU without an FPU. The QOA code in
`tjc.h` is a port of the MIT-licensed [reference](https://github.com/phoboslab/qoa).
It is bit-compatible with it: the same input gives the same bytes, and the frames
inside a `.tjc` are standard QOA frames.

### Encoder options

| Option                   | Default      | Meaning |
|--------------------------|--------------|---------|
| `--width` / `--height` / `--size WxH` | required | frame size |
| `--fps RATE`             | `30`         | frame rate of the input, stored in the stream |
| `--audio FILE.wav`       | none         | add a QOA audio track from a 16-bit PCM WAV |
| `--tile-size WxH` or `N` | `16x16`      | 1..255 each; TJC2/TJC3: `w*h` at most 12544 (e.g. 64x64, 112x112). See "Choosing a tile size" |
| `--motion-threshold K`   | `3`          | tile is dirty when any 8x8 block's mean abs diff per sample is > K; `0` = any change |
| `--refresh-mode`         | `rolling`    | `none`, `full` (full refresh every N frames), `rolling` (N tiles per frame, round robin) |
| `--refresh-param N`      | `30`         | N for the refresh mode |
| `--quality Q`            | `75`         | 1..100, libjpeg-style scaling of the standard quant tables |
| `--diff-ref`             | `last-coded` | what the new frame is compared to: `last-coded` source pixels or `recon` (decoded pixels) |
| `--keyframe-every N`     | `0`          | additionally force a full refresh every N frames |
| `--threads N`            | `0`          | encoder threads, `0` = all cores. The output is identical for any value |
| `--motion on\|off`        | `on`         | motion compensation (TJC3) |
| `--huffman adaptive\|standard` | `adaptive` | Huffman tables fitted to the video (TJC3) |
| `--effort fast\|normal\|best` | `normal` | encoder speed vs size |
| `--deadzone X`           | `0.33`       | AC rounding point, 0..0.5 (0.5 = plain rounding) |
| `--smart-refresh on\|off` | `on`         | rolling refresh skips recently sent tiles |
| `--skip-invisible T`     | `0` (off)    | drop updates that change no 8x8 block by more than T per pixel; 1-2 for grainy video |
| `--motion-range N`       | `32`         | motion search range in pixels |
| `--sync-window N`        | automatic    | how many frames back a decoder may have to go to seek or join |
| `--peak-threshold P`     | `32`         | also re-send a tile when any pixel changed by more than P (`0` = off); fixes ghost specks |
| `--format tjc2\|tjc3\|tjc4` | automatic | stream format; automatic = the oldest one the options need |
| `--quadtree`             |              | TJC4 variable tile size: `--tile-size 64 --split-levels 3` |
| `--split-levels N`       | `0`          | TJC4 quadtree depth 0..4 (square tiles, smallest leaf even and >= 8) |
| `--bitrate RATE`         | off          | average target, audio included: `800k`, `2.5M`, `1500` (kbit/s); TJC4 by default |
| `--maxrate RATE`         | off          | hard cap over any `--bufsize` window |
| `--bufsize MS`           | `1000`       | cap window |
| `--min-quality Q`        | `10`         | lowest quality `--bitrate` may use (`--maxrate` may go lower) |

### Choosing a tile size

Smaller tiles skip more unchanged area but cost more per tile: 2 bytes of `tile_len`,
one bitmap bit, and at least one 8x8 block per component even when the tile covers
only a few pixels. Below 8x8 the overhead outweighs the savings fast. Here are the
numbers for 1920x1080 `testsrc2`, q75, K=5, 4 threads:

| Tile  | Bytes/frame | Dirty | vs raw | Encode speed |
|-------|-------------|-------|--------|--------------|
| 32x32 | 50.6 KB     | 17.0% | 61:1   | 118 fps |
| 16x16 | 50.7 KB     | 11.9% | 61:1   | 140 fps |
| 8x8   | 79.8 KB     |  9.3% | 39:1   | 111 fps |
| 4x4   | 201 KB      |  7.9% | 15:1   | 52 fps  |
| 2x2   | 468 KB      |  7.1% | 6.6:1  | 22 fps  |
| 1x1   | 909 KB      |  6.4% | 3.4:1  | 13 fps  |

A 1x1 tile turns each pixel into three 8x8 blocks plus a length field, which is
more bytes than the raw pixel. A full refresh at 1x1 is bigger than the raw frame.
Use 16x16 or 8x8 unless you have a specific reason not to.

### Speed

- Build optimized: `/O2` (MSVC) or `-O2`/`-O3` (gcc/clang). An unoptimized build
  is about 3.5x slower.
- The encoder uses every core by default (`--threads`). Tiles are split into
  contiguous runs and joined in order, so the stream doesn't depend on the thread
  count. The decoder stays single-threaded for the embedded target.
- Worst case (every tile changes every frame) at 1080p: about 30 fps on 1 thread,
  about 60 fps on 4 threads.

With `--diff-ref recon` at low quality the quantization error alone can exceed a
small K and keep static tiles dirty; `last-coded` doesn't have that problem and still
catches slow drift, because the difference is measured from what was last sent.

## Using the library

```cpp
#define TJC_IMPLEMENTATION   // in exactly one .cpp
#include "tjc.h"

// Encoder
tjc::Config cfg;
cfg.width = 320; cfg.height = 240;
cfg.tile_w = cfg.tile_h = 16;
cfg.fps_num = 30000; cfg.fps_den = 1001;
cfg.audio_channels = 2; cfg.audio_rate = 48000;  // leave 0 for no audio
tjc::Encoder enc;
if (!enc.init(cfg)) puts(enc.error());
std::vector<uint8_t> out;
enc.write_stream_header(out);
enc.push_audio(pcm, enc.samples_for_next_frame());  // interleaved int16; push ahead is fine
enc.encode_frame(yuv420p_frame, out);               // appends one frame record (video + audio)

// Decoder (pull-style input callback; tjc::MemoryReader for buffers)
size_t rd(void* f, void* dst, size_t n) { return fread(dst, 1, n, (FILE*)f); }
tjc::Decoder dec;
dec.read_header(rd, stdin);
while (dec.decode_frame(rd, stdin) == tjc::Status::Ok) {
    dec.copy_frame(frame);                     // or use dec.plane(c) / dec.stride(c)
    play(dec.audio(), dec.audio_samples());    // this frame's audio, interleaved
}
```

Compile-time switches (define before the implementation include):

- `TJC_NO_ENCODER`: decoder only (video, motion compensation, audio): about 21 KB of
  x86-64 code with `-Os -ffunction-sections -Wl,--gc-sections`, 41 KB at `-O2`.
- `TJC_STREAM_BIG_ENDIAN 1`: big-endian multi-byte stream fields.
- `TJC_MAX_PIXELS n`: largest `width*height` the decoder will accept.
- `TJC_NO_THREADS`: single-threaded encoder with no `<thread>` dependency.

The decode path is integer-only (IJG "islow" DCT) and deterministic: decoder output
is bit-identical to the encoder's internal reconstruction on every platform, so the
two framebuffers never drift apart. Decoder RAM is the padded frame (1.5 bytes per
pixel), a second frame of the same size for motion compensated streams only (1.4 MB
at 720p, 3.1 MB at 1080p), one tile payload buffer, and one frame's worth of audio.

TJC4 settings live in `Config`: `split_levels`, `bitrate_kbps`, `max_bitrate_kbps`,
`buffer_ms`, `min_quality` and `format`. `peak_threshold` works with every format.
`FrameStats` reports each frame's `quality`, quadtree `leaves` and the
`deferred_tiles` held back by the cap.

Seeking: `tjc::SyncTracker` follows what every tile's content depends on and gives,
for any frame, the frame to start decoding from to reproduce it exactly. While
indexing, feed it `Decoder::decode_frame(..., apply=false)` and
`unit_layout()`/`unit_info()`: the quadtree leaves for TJC4, the tiles otherwise.
Before decoding from that frame, load the Huffman tables in force there. Parse the
last frame at or before it that has `FrameStats::tables` set (again with
`apply=false`), or call `reset_tables()` if there is none. TJC Studio does exactly
this.

## Stream format (TJC4)

Multi-byte fields are little-endian (see `TJC_STREAM_BIG_ENDIAN`), except inside
the QOA frames, which keep QOA's own big-endian layout.

```
Stream header, 34 bytes
  magic "TJC4" (4) | width (2) | height (2) | tile_w (1) | tile_h (1)
  chroma_format (1, 0 = 4:2:0) | refresh_mode (1) | refresh_param (2)
  quality (1, the highest any frame uses) | reserved (3, zero)
  fps_num (4) | fps_den (4)
  audio_codec (1, 0 = none, 1 = QOA) | audio_channels (1) | audio_rate (4)
  flags (1: bit0 adaptive Huffman, bit1 motion compensation)
  split_levels (1: quadtree depth 0..4; TJC3: reserved, zero)

Frame
  frame_num (4) | flags (1: bit0 force refresh, bit1 tables follow)
  quality (1, 1..100: this frame's quant tables)  [TJC4 only]
  tables, if flags bit1: 4 x (16 code-length counts + symbols), JPEG DHT layout;
      DC luma, AC luma, DC chroma, AC chroma; in force until replaced
  dirty_bitmap, ceil(tiles/8) bytes, only if force refresh is not set
      tile i (raster order) = bit (i & 7) of byte (i >> 3)
  per dirty tile: tile_len | payload (tile_len)
      tile_len: TJC4 LEB128 varint (1-3 bytes, low 7 bits first); TJC3: 2 bytes
  if audio_codec != 0: audio_len (4) | QOA frames (audio_len bytes)
      standard QOA frames of up to 5120 samples per channel, holding the samples
      that play during this video frame

TJC4 tile payload: a quadtree
  node above the deepest level: split (1 bit); 1 = four children follow
      (top-left, top-right, bottom-left, bottom-right)
  leaf: mode, except in force refresh frames (all intra)
      with motion: '0' inter, '10' intra, '11' skip (keep the picture)
      without motion: '0' intra, '1' skip
  intra / inter leaves are coded like a TJC3 tile (below) over the leaf's area;
  the vector predictor is the last inter leaf's vector in the same tile row,
  (0, 0) at the start of each row. With split_levels 0 a tile is one leaf.

TJC3 tile payload
  [motion streams only] mode (1 bit): 0 intra, 1 inter;
      inter: mvx, mvy (signed Exp-Golomb, half-pels, minus the previous coded tile's
      vector in the same tile row if it was inter), residual flag (1 bit)
  Huffman-coded 8x8 blocks, MSB first, zero-padded, no 0xFF stuffing. Order: Y
  blocks, then Cb, then Cr, raster within the tile. The DC predictor resets per
  component per tile. Inter blocks are residuals added to the prediction (previous
  frame moved by the vector, half-pel averaging; chroma vector = luma >> 1).
```

TJC2 is the same without the header flags, frame tables and tile mode bits (always
standard tables, always intra). TJC1 streams (the first 18 header bytes with magic
`TJC1`, no frame rate, no audio) still decode; they are treated as 30 fps.

Changes from the original plan:

- The first reserved header byte carries `quality`, so the decoder can rebuild the
  same quant tables (the tables themselves stay baked in).
- Tiles can be any size from 1 to 255. Partial 8x8 blocks are padded by edge
  replication, and each 4:2:0 chroma sample belongs to the tile holding its top-left
  luma pixel, so tiles never overlap (see the comment at the top of `tjc.h`).
- Frame sizes don't have to be tile multiples (or even): the codec pads internally
  by edge replication and crops on output.
- The motion threshold applies per 8x8 block instead of averaging over the whole
  tile, so a small object moving inside a large tile is still detected.
- Audio, a v1 non-goal, is now supported (QOA, per-frame interleave), along with a
  frame rate in the stream (TJC2).
- TJC3 adds motion compensation and adaptive Huffman tables; the plan's static
  tables and intra-only tiles remain available as `--format tjc2`.
- TJC4 adds the quadtree (variable tile size), a quality per frame for rate
  control, and varint tile lengths (so tiles up to 255x255 are allowed).
