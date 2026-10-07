# TJC — Tiled-JPEG Codec

A small intra-only video codec for low-power decoders (written with a single-core
MIPS 24KEc in mind). Each frame is split into tiles; only tiles that changed (plus
tiles picked by a refresh policy) are sent, each as a tiny baseline-JPEG-style blob.
The decoder keeps a persistent framebuffer and patches it.

| File          | What                                                              |
|---------------|-------------------------------------------------------------------|
| `tjc.h`       | single-header library: encoder + decoder (stb-style)              |
| `tjc.cpp`     | `tjc` command line tool: `encode`, `decode`, `info`               |
| `tjc_test.cpp`| unit tests (bitstream, DCT, Huffman, bitmap, SAD, refresh) + end-to-end |

## Build

```bash
g++ -std=c++17 -O2 -o tjc tjc.cpp
g++ -std=c++17 -O2 -o tjc_test tjc_test.cpp && ./tjc_test
```

or `cmake -B build && cmake --build build && ./build/tjc_test`.

Optional: `ln -s tjc tjc_encode && ln -s tjc tjc_decode` lets you call the tool
without the subcommand.

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

The decoder reads width/height/tile size from the stream; it needs no options.
Add `-v` for a per-frame log on stderr and `--psnr` (encoder) for quality numbers.

### Encoder options

| Option                   | Default      | Meaning |
|--------------------------|--------------|---------|
| `--width` / `--height` / `--size WxH` | required | frame size |
| `--tile-size WxH` or `N` | `16x16`      | multiples of 16, up to 240; `w*h` at most 12544 (e.g. 64x64, 128x64, 112x112) |
| `--motion-threshold K`   | `3`          | tile is dirty when any 8x8 block's mean abs diff per sample is > K; `0` = any change |
| `--refresh-mode`         | `rolling`    | `none`, `full` (full refresh every N frames), `rolling` (N tiles per frame, round robin) |
| `--refresh-param N`      | `30`         | N for the refresh mode |
| `--quality Q`            | `75`         | 1..100, libjpeg-style scaling of the standard quant tables |
| `--diff-ref`             | `last-coded` | what the new frame is compared to: `last-coded` source pixels or `recon` (decoded pixels) |
| `--keyframe-every N`     | `0`          | additionally force a full refresh every N frames |

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
tjc::Encoder enc;
if (!enc.init(cfg)) puts(enc.error());
std::vector<uint8_t> out;
enc.write_stream_header(out);
enc.encode_frame(yuv420p_frame, out);          // appends one frame record

// Decoder (pull-style input callback; tjc::MemoryReader for buffers)
size_t rd(void* f, void* dst, size_t n) { return fread(dst, 1, n, (FILE*)f); }
tjc::Decoder dec;
dec.read_header(rd, stdin);
while (dec.decode_frame(rd, stdin) == tjc::Status::Ok)
    dec.copy_frame(frame);                     // or use dec.plane(c) / dec.stride(c)
```

Compile-time switches (define before the implementation include):

- `TJC_NO_ENCODER`: decoder only. About 18 KB of code on x86-64 at `-O2`.
- `TJC_STREAM_BIG_ENDIAN 1`: big-endian multi-byte stream fields.
- `TJC_MAX_PIXELS n`: largest `width*height` the decoder will accept.

The decode path is integer-only (IJG "islow" DCT) and deterministic: decoder output
is bit-identical to the encoder's internal reconstruction on every platform, so the
two framebuffers never drift apart. Decoder RAM is the padded frame (1.5 bytes per
pixel) plus one tile payload buffer.

## Stream format (TJC1)

Multi-byte fields are little-endian (see `TJC_STREAM_BIG_ENDIAN`).

```
Stream header, 18 bytes
  magic "TJC1" (4) | width (2) | height (2) | tile_w (1) | tile_h (1)
  chroma_format (1, 0 = 4:2:0) | refresh_mode (1) | refresh_param (2)
  quality (1) | reserved (3, zero)

Frame
  frame_num (4) | force_refresh (1)
  dirty_bitmap, ceil(tiles/8) bytes, only if force_refresh == 0
      tile i (raster order) = bit (i & 7) of byte (i >> 3)
  per dirty tile: tile_len (2) | payload (tile_len)

Tile payload
  Huffman-coded 8x8 blocks (standard JPEG Annex K tables), MSB first, zero-padded,
  no 0xFF stuffing. Order: Y blocks, then Cb, then Cr, raster within the tile.
  The DC predictor resets per component per tile, so each tile decodes on its own.
```

Changes from the original plan:

- The first reserved header byte carries `quality`, so the decoder can rebuild the
  same quant tables (the tables themselves stay baked in).
- Tiles must be multiples of 16 so that 4:2:0 chroma tiles are whole 8x8 blocks.
- Frame sizes don't have to be tile multiples (or even): the codec pads internally
  by edge replication and crops on output.
- The motion threshold applies per 8x8 block instead of averaging over the whole
  tile, so a small object moving inside a large tile is still detected.
