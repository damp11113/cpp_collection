// tjc.h - Tiled-JPEG Codec (stream format "TJC1"), single-header C++17 library.
//
// A frame is split into fixed-size tiles. Each frame only carries the tiles that
// changed (plus tiles picked by the refresh policy); every carried tile is coded as
// a small baseline-JPEG-style blob (integer DCT, static quant tables, static JPEG
// Huffman tables). The decoder keeps a persistent framebuffer and patches it.
//
// Usage: in exactly ONE .cpp file
//     #define TJC_IMPLEMENTATION
//     #include "tjc.h"
// and include "tjc.h" normally everywhere else.
//
// Compile-time options (define before including in the implementation file):
//     TJC_STREAM_BIG_ENDIAN 1  write/read multi-byte stream fields big-endian (default 0 = little)
//     TJC_NO_ENCODER           leave out the encoder (decode-only embedded builds)
//     TJC_MAX_PIXELS n         largest width*height the decoder accepts (default 1<<26)
//
// The decode path is integer-only and deterministic: the same stream gives
// bit-identical output on every platform, and the encoder's internal
// reconstruction is bit-identical to what the decoder produces.
//
// ---------------------------------------------------------------------------
// Stream format (all multi-byte fields little-endian unless TJC_STREAM_BIG_ENDIAN)
//
//   Stream header, 18 bytes, once:
//     magic         4  "TJC1"
//     width         2  pixels (luma)
//     height        2  pixels (luma)
//     tile_w        1  multiple of 16, 16..240
//     tile_h        1  multiple of 16, 16..240
//     chroma_format 1  0 = 4:2:0 (only value accepted in v1)
//     refresh_mode  1  0 = NONE, 1 = FULL_PERIODIC, 2 = ROLLING (informational)
//     refresh_param 2  N for FULL_PERIODIC, tiles per frame for ROLLING
//     quality       1  1..100, selects the quant tables (was reserved byte 0)
//     reserved      3  zero
//
//   Per frame:
//     frame_num     4  uint32
//     force_refresh 1  0 or 1; 1 = every tile follows and no bitmap is sent
//     dirty_bitmap  ceil(tiles/8) bytes, only when force_refresh == 0.
//                   Tile i (raster order) is bit (i & 7) of byte (i >> 3), LSB first.
//     per dirty tile, raster order:
//       tile_len    2  payload size in bytes (1..65535)
//       payload     tile_len bytes
//
//   Tile payload: Huffman-coded 8x8 blocks, MSB-first bit packing, zero-padded to a
//   byte, no marker byte stuffing. Block order: all Y blocks of the tile (raster),
//   then all Cb blocks, then all Cr blocks. The DC predictor resets to 0 at the start
//   of each component of each tile, so every tile decodes on its own.
//
// Frame dimensions need not be multiples of the tile size: the codec works on an
// internal frame padded up to whole tiles (edge pixels replicated) and crops on output.
// ---------------------------------------------------------------------------

#ifndef TJC_H_INCLUDED
#define TJC_H_INCLUDED

#include <cstddef>
#include <cstdint>
#include <vector>

#ifndef TJC_STREAM_BIG_ENDIAN
#define TJC_STREAM_BIG_ENDIAN 0
#endif

#ifndef TJC_MAX_PIXELS
#define TJC_MAX_PIXELS (1u << 26)
#endif

namespace tjc {

constexpr int kStreamHeaderSize = 18;
constexpr int kFrameHeaderSize = 5;  // without the bitmap
constexpr int kMaxTilePayload = 65535;

enum class RefreshMode : uint8_t { None = 0, FullPeriodic = 1, Rolling = 2 };

// What the encoder compares the incoming frame against to decide if a tile is dirty.
//   LastCoded:     the source pixels as they were when the tile was last sent. Slow
//                  drift still accumulates until it crosses the threshold, and static
//                  content is never re-sent just because of quantization noise.
//   Reconstructed: the decoded (lossy) pixels the decoder is showing. Catches every
//                  deviation, but at low quality the quantization error alone can
//                  exceed a small threshold and keep static tiles dirty.
enum class DiffReference : uint8_t { LastCoded = 0, Reconstructed = 1 };

struct Config {
    uint16_t width = 0;
    uint16_t height = 0;
    uint8_t tile_w = 16;
    uint8_t tile_h = 16;

    // A tile is dirty when any of its 8x8 blocks (Y, Cb or Cr) has SAD > 64 * K,
    // i.e. the mean absolute difference per sample in that block exceeds K. Judging
    // per block keeps the meaning of K independent of the tile size: a small moving
    // object inside a big tile is not averaged away. K = 0 means any change.
    uint32_t motion_threshold_k = 3;

    RefreshMode refresh_mode = RefreshMode::Rolling;
    uint32_t refresh_param = 30;  // N for FullPeriodic, tiles per frame for Rolling

    uint8_t quality = 75;  // 1..100, libjpeg-style scaling of the standard tables

    DiffReference diff_reference = DiffReference::LastCoded;
};

struct StreamHeader {
    uint16_t width = 0;
    uint16_t height = 0;
    uint8_t tile_w = 0;
    uint8_t tile_h = 0;
    uint8_t chroma_format = 0;
    RefreshMode refresh_mode = RefreshMode::None;
    uint16_t refresh_param = 0;
    uint8_t quality = 0;
};

struct FrameStats {
    uint32_t frame_num = 0;
    bool force_refresh = false;
    uint32_t dirty_tiles = 0;
    uint32_t total_tiles = 0;
    size_t bytes = 0;  // encoded size of the frame record
};

enum class Status {
    Ok = 0,
    EndOfStream,     // clean end of input at a frame boundary
    Truncated,       // input ended in the middle of a record
    BadHeader,       // stream header invalid or unsupported
    Corrupt,         // malformed frame or tile data
    OutOfMemory,
    NotInitialized,  // decode_frame() before read_header()
};

const char* status_string(Status s);

// Pull-style input. Return the number of bytes copied into dst (0 = end of input).
// Short reads are fine; the decoder keeps calling until it has what it needs.
typedef size_t (*ReadFn)(void* user, void* dst, size_t n);

struct MemoryReader {
    const uint8_t* data = nullptr;
    size_t size = 0;
    size_t pos = 0;
    static size_t read(void* user, void* dst, size_t n);
};

// Geometry derived from the stream header.
struct Layout {
    int width = 0, height = 0;        // visible luma size
    int cwidth = 0, cheight = 0;      // visible chroma size, ceil(w/2) x ceil(h/2)
    int tile_w = 0, tile_h = 0;
    int cols = 0, rows = 0, tiles = 0;
    int pwidth = 0, pheight = 0;      // padded luma size (whole tiles)
    int pcwidth = 0, pcheight = 0;    // padded chroma size
};

// Size in bytes of one tightly packed YUV420P frame.
size_t frame_size(int width, int height);

// Upper bound of one tile payload for the given tile size.
size_t worst_case_tile_bytes(int tile_w, int tile_h);

// Returns nullptr if the geometry is usable, otherwise a reason.
const char* validate_geometry(int width, int height, int tile_w, int tile_h);

Layout make_layout(int width, int height, int tile_w, int tile_h);

void write_stream_header(const StreamHeader& h, uint8_t out[kStreamHeaderSize]);
Status parse_stream_header(const uint8_t in[kStreamHeaderSize], StreamHeader* h);

// Low-level building blocks, exposed for tests and for anyone porting the hot loops.
namespace detail {

extern const uint8_t kZigzag[64];  // zigzag position -> natural (row-major) index

// Forward DCT of an 8x8 block of 8-bit samples (level-shifted internally).
// Output is natural order, scaled up by 8 relative to the orthonormal DCT.
void fdct8x8(const uint8_t* src, int stride, int32_t out[64]);

// Inverse DCT of dequantized natural-order coefficients (each within [-2048, 2047]).
void idct8x8(const int32_t in[64], uint8_t* dst, int stride);

// Fast path for blocks whose AC coefficients are all zero. Bit-identical to idct8x8.
void idct8x8_dc(int32_t dc, uint8_t* dst, int stride);

// Quant tables for a quality (1..100), natural order.
void build_quant_tables(int quality, uint16_t luma[64], uint16_t chroma[64]);

struct HuffEnc {
    uint16_t code[256];
    uint8_t size[256];  // 0 = symbol not in table
};

struct HuffDec {
    uint16_t lut[512];  // 9-bit lookahead: (length << 8) | symbol, 0 = longer code
    int32_t maxcode[17];
    int32_t valoff[17];
    uint8_t vals[256];
};

struct HuffTables {
    HuffEnc enc_dc[2], enc_ac[2];  // [0] = luma, [1] = chroma
    HuffDec dec_dc[2], dec_ac[2];
};

const HuffTables& huff_tables();

class BitWriter {
public:
    explicit BitWriter(std::vector<uint8_t>& out) : out_(out) {}
    void put(uint32_t bits, int len);
    void flush();  // pad the last byte with zeros

private:
    std::vector<uint8_t>& out_;
    uint32_t acc_ = 0;
    int n_ = 0;
};

class BitReader {
public:
    BitReader(const uint8_t* data, size_t size) : p_(data), size_(size) {}
    void refill() {
        while (cnt_ <= 24) {
            uint32_t b = pos_ < size_ ? p_[pos_] : 0u;
            ++pos_;
            buf_ |= b << (24 - cnt_);
            cnt_ += 8;
        }
    }
    uint32_t peek(int n) const { return buf_ >> (32 - n); }  // 1 <= n <= 24, after refill()
    void skip(int n) { buf_ <<= n; cnt_ -= n; }
    uint32_t get(int n) {
        if (n == 0) return 0;
        refill();
        uint32_t v = peek(n);
        skip(n);
        return v;
    }
    // True once more bits were consumed than the buffer holds.
    bool overrun() const { return pos_ * 8 - size_t(cnt_) > size_ * 8; }

private:
    const uint8_t* p_;
    size_t size_;
    size_t pos_ = 0;
    uint32_t buf_ = 0;
    int cnt_ = 0;
};

// Returns the decoded symbol, or -1 for an invalid code.
int huff_decode(BitReader& br, const HuffDec& t);

// Entropy-code one block of quantized coefficients given in zigzag order.
// Coefficients must lie in [-1023, 1023]. dc_pred is updated.
void encode_block(BitWriter& bw, const int16_t zz[64], int& dc_pred, int table);

// Decode one block into natural-order dequantized coefficients (clamped to
// [-2048, 2047]). dq is the quant table in zigzag order. dc_pred is updated.
// Returns -1 on bad data, 0 when only DC is present, 1 when any AC is present.
int decode_block(BitReader& br, int table, const uint16_t dq[64], int& dc_pred, int32_t coef[64]);

void pack_bitmap(const uint8_t* flags, int n, uint8_t* out);
void unpack_bitmap(const uint8_t* in, int n, uint8_t* flags);

// Sum of absolute differences between two w x h areas, stopping early once it
// exceeds `limit` (the returned value is then > limit but not exact).
uint64_t sad(const uint8_t* a, int astride, const uint8_t* b, int bstride, int w, int h, uint64_t limit);

class RefreshPolicy {
public:
    void reset(RefreshMode mode, uint32_t param, int tiles);
    // Whether frame `frame_index` (0-based count of frames since start) is a full refresh.
    bool forced(uint32_t frame_index) const;
    // Rolling mode: mark the next `param` tiles (round-robin) as dirty.
    void apply(uint8_t* dirty);
    int cursor() const { return cursor_; }

private:
    RefreshMode mode_ = RefreshMode::None;
    uint32_t param_ = 0;
    int tiles_ = 0;
    int cursor_ = 0;
};

}  // namespace detail

// Padded planar YUV420P storage used by both encoder and decoder.
struct Planes {
    std::vector<uint8_t> p[3];
    int stride[3] = {0, 0, 0};
    bool alloc(const Layout& l);
};

class Decoder {
public:
    // Reads and validates the 18-byte stream header, allocates the framebuffer
    // (initialized to black). Must be called first.
    Status read_header(ReadFn fn, void* user);

    // Reads and applies one frame record. Returns EndOfStream at clean EOF.
    Status decode_frame(ReadFn fn, void* user, FrameStats* stats = nullptr);

    const StreamHeader& header() const { return header_; }
    const Layout& layout() const { return layout_; }

    // Copies the current visible frame as tightly packed YUV420P (frame_size() bytes).
    void copy_frame(uint8_t* dst) const;

    // Direct access to the padded planes (0 = Y, 1 = Cb, 2 = Cr).
    const uint8_t* plane(int c) const { return fb_.p[c].data(); }
    int stride(int c) const { return fb_.stride[c]; }

private:
    bool decode_tile(const uint8_t* data, size_t size, int tile);

    StreamHeader header_;
    Layout layout_;
    Planes fb_;
    uint16_t dq_[2][64] = {};  // dequant tables, zigzag order
    std::vector<uint8_t> dirty_, bitmap_, payload_;
    bool ready_ = false;
};

#ifndef TJC_NO_ENCODER
class Encoder {
public:
    // Validates the config and allocates buffers. On failure returns false and
    // error() says why.
    bool init(const Config& cfg);
    const char* error() const { return error_; }

    const Config& config() const { return cfg_; }
    const Layout& layout() const { return layout_; }
    StreamHeader stream_header() const;

    // Appends the 18-byte stream header.
    void write_stream_header(std::vector<uint8_t>& out) const;

    // Encodes one tightly packed YUV420P frame of frame_size(width, height) bytes
    // and appends the frame record to out.
    void encode_frame(const uint8_t* yuv420p, std::vector<uint8_t>& out, FrameStats* stats = nullptr);

    // Same, with separate planes and strides.
    void encode_frame(const uint8_t* y, int ystride, const uint8_t* u, int ustride,
                      const uint8_t* v, int vstride, std::vector<uint8_t>& out,
                      FrameStats* stats = nullptr);

    // Make the next frame a full refresh (e.g. a new receiver joined).
    void force_keyframe() { keyframe_request_ = true; }

    // The encoder's shadow copy of the decoder framebuffer, as tightly packed YUV420P.
    void copy_recon(uint8_t* dst) const;

private:
    void load(const uint8_t* y, int ys, const uint8_t* u, int us, const uint8_t* v, int vs);
    bool tile_changed(int tile) const;
    void encode_tile(int tile, std::vector<uint8_t>& out);
    void copy_tile(const Planes& from, Planes& to, int tile);

    Config cfg_;
    Layout layout_;
    Planes cur_, last_coded_, recon_;
    uint16_t q_[2][64] = {};  // natural order
    detail::RefreshPolicy policy_;
    std::vector<uint8_t> dirty_;
    uint32_t frame_num_ = 0;
    uint32_t frame_index_ = 0;
    bool keyframe_request_ = false;
    bool ready_ = false;
    const char* error_ = "not initialized";
};
#endif  // TJC_NO_ENCODER

}  // namespace tjc

#endif  // TJC_H_INCLUDED

// ===========================================================================
// Implementation
// ===========================================================================
#ifdef TJC_IMPLEMENTATION
#ifndef TJC_IMPLEMENTATION_DONE
#define TJC_IMPLEMENTATION_DONE

#include <algorithm>
#include <cstring>

static_assert((-5 >> 1) == -3, "tjc needs arithmetic right shift of negative integers");

namespace tjc {

// ---------------------------------------------------------------------------
// Byte order helpers
// ---------------------------------------------------------------------------
namespace {

inline void put_u16(uint8_t* p, uint32_t v) {
#if TJC_STREAM_BIG_ENDIAN
    p[0] = uint8_t(v >> 8); p[1] = uint8_t(v);
#else
    p[0] = uint8_t(v); p[1] = uint8_t(v >> 8);
#endif
}

inline void put_u32(uint8_t* p, uint32_t v) {
#if TJC_STREAM_BIG_ENDIAN
    p[0] = uint8_t(v >> 24); p[1] = uint8_t(v >> 16); p[2] = uint8_t(v >> 8); p[3] = uint8_t(v);
#else
    p[0] = uint8_t(v); p[1] = uint8_t(v >> 8); p[2] = uint8_t(v >> 16); p[3] = uint8_t(v >> 24);
#endif
}

inline uint32_t get_u16(const uint8_t* p) {
#if TJC_STREAM_BIG_ENDIAN
    return (uint32_t(p[0]) << 8) | p[1];
#else
    return uint32_t(p[0]) | (uint32_t(p[1]) << 8);
#endif
}

inline uint32_t get_u32(const uint8_t* p) {
#if TJC_STREAM_BIG_ENDIAN
    return (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | p[3];
#else
    return uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24);
#endif
}

inline uint8_t clamp_u8(int32_t v) { return uint8_t(v < 0 ? 0 : (v > 255 ? 255 : v)); }

// Reads exactly n bytes unless the input ends; returns the count read.
size_t read_full(ReadFn fn, void* user, void* dst, size_t n) {
    size_t got = 0;
    uint8_t* d = static_cast<uint8_t*>(dst);
    while (got < n) {
        size_t r = fn(user, d + got, n - got);
        if (r == 0) break;
        got += r;
    }
    return got;
}

const uint8_t kStdLumaQ[64] = {
    16, 11, 10, 16, 24, 40, 51, 61,
    12, 12, 14, 19, 26, 58, 60, 55,
    14, 13, 16, 24, 40, 57, 69, 56,
    14, 17, 22, 29, 51, 87, 80, 62,
    18, 22, 37, 56, 68, 109, 103, 77,
    24, 35, 55, 64, 81, 104, 113, 92,
    49, 64, 78, 87, 103, 121, 120, 101,
    72, 92, 95, 98, 112, 100, 103, 99,
};

const uint8_t kStdChromaQ[64] = {
    17, 18, 24, 47, 99, 99, 99, 99,
    18, 21, 26, 66, 99, 99, 99, 99,
    24, 26, 56, 99, 99, 99, 99, 99,
    47, 66, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
};

// Standard JPEG Huffman tables (ITU T.81 Annex K.3).
const uint8_t kDcLumaBits[16] = {0, 1, 5, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0};
const uint8_t kDcLumaVals[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
const uint8_t kDcChromaBits[16] = {0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0};
const uint8_t kDcChromaVals[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

const uint8_t kAcLumaBits[16] = {0, 2, 1, 3, 3, 2, 4, 3, 5, 5, 4, 4, 0, 0, 1, 0x7d};
const uint8_t kAcLumaVals[162] = {
    0x01, 0x02, 0x03, 0x00, 0x04, 0x11, 0x05, 0x12, 0x21, 0x31, 0x41, 0x06, 0x13, 0x51, 0x61, 0x07,
    0x22, 0x71, 0x14, 0x32, 0x81, 0x91, 0xa1, 0x08, 0x23, 0x42, 0xb1, 0xc1, 0x15, 0x52, 0xd1, 0xf0,
    0x24, 0x33, 0x62, 0x72, 0x82, 0x09, 0x0a, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x25, 0x26, 0x27, 0x28,
    0x29, 0x2a, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3a, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49,
    0x4a, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5a, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69,
    0x6a, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79, 0x7a, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89,
    0x8a, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9a, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7,
    0xa8, 0xa9, 0xaa, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xc2, 0xc3, 0xc4, 0xc5,
    0xc6, 0xc7, 0xc8, 0xc9, 0xca, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7, 0xd8, 0xd9, 0xda, 0xe1, 0xe2,
    0xe3, 0xe4, 0xe5, 0xe6, 0xe7, 0xe8, 0xe9, 0xea, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8,
    0xf9, 0xfa,
};

const uint8_t kAcChromaBits[16] = {0, 2, 1, 2, 4, 4, 3, 4, 7, 5, 4, 4, 0, 1, 2, 0x77};
const uint8_t kAcChromaVals[162] = {
    0x00, 0x01, 0x02, 0x03, 0x11, 0x04, 0x05, 0x21, 0x31, 0x06, 0x12, 0x41, 0x51, 0x07, 0x61, 0x71,
    0x13, 0x22, 0x32, 0x81, 0x08, 0x14, 0x42, 0x91, 0xa1, 0xb1, 0xc1, 0x09, 0x23, 0x33, 0x52, 0xf0,
    0x15, 0x62, 0x72, 0xd1, 0x0a, 0x16, 0x24, 0x34, 0xe1, 0x25, 0xf1, 0x17, 0x18, 0x19, 0x1a, 0x26,
    0x27, 0x28, 0x29, 0x2a, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3a, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48,
    0x49, 0x4a, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5a, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68,
    0x69, 0x6a, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79, 0x7a, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
    0x88, 0x89, 0x8a, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9a, 0xa2, 0xa3, 0xa4, 0xa5,
    0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xc2, 0xc3,
    0xc4, 0xc5, 0xc6, 0xc7, 0xc8, 0xc9, 0xca, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7, 0xd8, 0xd9, 0xda,
    0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7, 0xe8, 0xe9, 0xea, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8,
    0xf9, 0xfa,
};

void build_huff(const uint8_t bits[16], const uint8_t* vals, detail::HuffEnc& enc, detail::HuffDec& dec) {
    std::memset(&enc, 0, sizeof(enc));
    std::memset(&dec, 0, sizeof(dec));
    uint32_t code = 0;
    int k = 0;
    for (int len = 1; len <= 16; ++len) {
        int n = bits[len - 1];
        dec.maxcode[len] = -1;
        if (n) {
            dec.valoff[len] = k - int32_t(code);
            for (int i = 0; i < n; ++i, ++k, ++code) {
                uint8_t sym = vals[k];
                dec.vals[k] = sym;
                enc.code[sym] = uint16_t(code);
                enc.size[sym] = uint8_t(len);
                if (len <= 9) {
                    uint32_t first = code << (9 - len);
                    uint32_t count = 1u << (9 - len);
                    for (uint32_t j = 0; j < count; ++j) dec.lut[first + j] = uint16_t((len << 8) | sym);
                }
            }
            dec.maxcode[len] = int32_t(code) - 1;
        }
        code <<= 1;
    }
}

inline int bit_count(int v) {
    int n = 0;
    while (v) { ++n; v >>= 1; }
    return n;
}

}  // namespace

// ---------------------------------------------------------------------------
// Public helpers
// ---------------------------------------------------------------------------
const char* status_string(Status s) {
    switch (s) {
        case Status::Ok: return "ok";
        case Status::EndOfStream: return "end of stream";
        case Status::Truncated: return "truncated input";
        case Status::BadHeader: return "bad or unsupported stream header";
        case Status::Corrupt: return "corrupt frame data";
        case Status::OutOfMemory: return "out of memory";
        case Status::NotInitialized: return "stream header not read yet";
    }
    return "unknown";
}

size_t MemoryReader::read(void* user, void* dst, size_t n) {
    MemoryReader* r = static_cast<MemoryReader*>(user);
    size_t left = r->size - r->pos;
    if (n > left) n = left;
    std::memcpy(dst, r->data + r->pos, n);
    r->pos += n;
    return n;
}

size_t frame_size(int width, int height) {
    size_t cw = size_t(width + 1) / 2, ch = size_t(height + 1) / 2;
    return size_t(width) * size_t(height) + 2 * cw * ch;
}

size_t worst_case_tile_bytes(int tile_w, int tile_h) {
    // 4:2:0 -> 1.5 samples per pixel. Per block: DC <= 16+11 bits, 63 AC <= 16+10 bits.
    size_t blocks = size_t(tile_w) * size_t(tile_h) * 3 / 128;
    return (blocks * (27 + 63 * 26) + 7) / 8;
}

const char* validate_geometry(int width, int height, int tile_w, int tile_h) {
    if (width <= 0 || height <= 0 || width > 65535 || height > 65535) return "width/height must be 1..65535";
    if (tile_w < 16 || tile_h < 16 || tile_w > 240 || tile_h > 240 || tile_w % 16 || tile_h % 16)
        return "tile width/height must be a multiple of 16 in 16..240";
    if (worst_case_tile_bytes(tile_w, tile_h) > size_t(kMaxTilePayload))
        return "tile too large: worst-case payload would not fit the 16-bit tile_len (keep w*h <= 12544, e.g. 112x112 or 128x64)";
    return nullptr;
}

Layout make_layout(int width, int height, int tile_w, int tile_h) {
    Layout l;
    l.width = width;
    l.height = height;
    l.cwidth = (width + 1) / 2;
    l.cheight = (height + 1) / 2;
    l.tile_w = tile_w;
    l.tile_h = tile_h;
    l.cols = (width + tile_w - 1) / tile_w;
    l.rows = (height + tile_h - 1) / tile_h;
    l.tiles = l.cols * l.rows;
    l.pwidth = l.cols * tile_w;
    l.pheight = l.rows * tile_h;
    l.pcwidth = l.pwidth / 2;
    l.pcheight = l.pheight / 2;
    return l;
}

void write_stream_header(const StreamHeader& h, uint8_t out[kStreamHeaderSize]) {
    std::memset(out, 0, kStreamHeaderSize);
    out[0] = 'T'; out[1] = 'J'; out[2] = 'C'; out[3] = '1';
    put_u16(out + 4, h.width);
    put_u16(out + 6, h.height);
    out[8] = h.tile_w;
    out[9] = h.tile_h;
    out[10] = h.chroma_format;
    out[11] = uint8_t(h.refresh_mode);
    put_u16(out + 12, h.refresh_param);
    out[14] = h.quality;
}

Status parse_stream_header(const uint8_t in[kStreamHeaderSize], StreamHeader* h) {
    if (in[0] != 'T' || in[1] != 'J' || in[2] != 'C' || in[3] != '1') return Status::BadHeader;
    StreamHeader r;
    r.width = uint16_t(get_u16(in + 4));
    r.height = uint16_t(get_u16(in + 6));
    r.tile_w = in[8];
    r.tile_h = in[9];
    r.chroma_format = in[10];
    if (in[11] > 2) return Status::BadHeader;
    r.refresh_mode = RefreshMode(in[11]);
    r.refresh_param = uint16_t(get_u16(in + 12));
    r.quality = in[14];
    if (r.chroma_format != 0 || r.quality < 1 || r.quality > 100) return Status::BadHeader;
    if (validate_geometry(r.width, r.height, r.tile_w, r.tile_h)) return Status::BadHeader;
    *h = r;
    return Status::Ok;
}

bool Planes::alloc(const Layout& l) {
    size_t ls = size_t(l.pwidth) * size_t(l.pheight);
    size_t cs = size_t(l.pcwidth) * size_t(l.pcheight);
    try {
        p[0].assign(ls, 16);
        p[1].assign(cs, 128);
        p[2].assign(cs, 128);
    } catch (...) {
        return false;
    }
    stride[0] = l.pwidth;
    stride[1] = stride[2] = l.pcwidth;
    return true;
}

namespace {

// Visible plane geometry for component c.
inline void visible_size(const Layout& l, int c, int& w, int& h) {
    w = c ? l.cwidth : l.width;
    h = c ? l.cheight : l.height;
}

// Tile rectangle of component c within its padded plane.
inline void tile_rect(const Layout& l, int tile, int c, int& x, int& y, int& w, int& h) {
    int tx = tile % l.cols, ty = tile / l.cols;
    w = c ? l.tile_w / 2 : l.tile_w;
    h = c ? l.tile_h / 2 : l.tile_h;
    x = tx * w;
    y = ty * h;
}

void copy_cropped(const Planes& src, const Layout& l, uint8_t* dst) {
    for (int c = 0; c < 3; ++c) {
        int w, h;
        visible_size(l, c, w, h);
        const uint8_t* s = src.p[c].data();
        for (int y = 0; y < h; ++y, dst += w) std::memcpy(dst, s + size_t(y) * src.stride[c], size_t(w));
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// detail: DCT, quantization, Huffman, bitstream, bitmap, SAD, refresh policy
// ---------------------------------------------------------------------------
namespace detail {

const uint8_t kZigzag[64] = {
    0, 1, 8, 16, 9, 2, 3, 10, 17, 24, 32, 25, 18, 11, 4, 5,
    12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13, 6, 7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63,
};

// Loeffler-Ligtenberg-Moschytz integer DCT, as in the IJG "islow" code:
// 13-bit fixed-point constants, 2 extra bits of precision between passes.
namespace {
constexpr int kConstBits = 13;
constexpr int kPass1Bits = 2;
constexpr int32_t FIX_0_298631336 = 2446;
constexpr int32_t FIX_0_390180644 = 3196;
constexpr int32_t FIX_0_541196100 = 4433;
constexpr int32_t FIX_0_765366865 = 6270;
constexpr int32_t FIX_0_899976223 = 7373;
constexpr int32_t FIX_1_175875602 = 9633;
constexpr int32_t FIX_1_501321110 = 12299;
constexpr int32_t FIX_1_847759065 = 15137;
constexpr int32_t FIX_1_961570560 = 16069;
constexpr int32_t FIX_2_053119869 = 16819;
constexpr int32_t FIX_2_562915447 = 20995;
constexpr int32_t FIX_3_072711026 = 25172;

inline int32_t descale(int32_t x, int n) { return (x + (int32_t(1) << (n - 1))) >> n; }
}  // namespace

void fdct8x8(const uint8_t* src, int stride, int32_t out[64]) {
    int32_t ws[64];
    for (int y = 0; y < 8; ++y) {
        const uint8_t* s = src + size_t(y) * stride;
        int32_t* o = ws + y * 8;
        int32_t tmp0 = s[0] + s[7] - 256, tmp7 = s[0] - s[7];
        int32_t tmp1 = s[1] + s[6] - 256, tmp6 = s[1] - s[6];
        int32_t tmp2 = s[2] + s[5] - 256, tmp5 = s[2] - s[5];
        int32_t tmp3 = s[3] + s[4] - 256, tmp4 = s[3] - s[4];

        int32_t tmp10 = tmp0 + tmp3, tmp13 = tmp0 - tmp3;
        int32_t tmp11 = tmp1 + tmp2, tmp12 = tmp1 - tmp2;
        o[0] = (tmp10 + tmp11) * (1 << kPass1Bits);
        o[4] = (tmp10 - tmp11) * (1 << kPass1Bits);
        int32_t z1 = (tmp12 + tmp13) * FIX_0_541196100;
        o[2] = descale(z1 + tmp13 * FIX_0_765366865, kConstBits - kPass1Bits);
        o[6] = descale(z1 - tmp12 * FIX_1_847759065, kConstBits - kPass1Bits);

        z1 = tmp4 + tmp7;
        int32_t z2 = tmp5 + tmp6, z3 = tmp4 + tmp6, z4 = tmp5 + tmp7;
        int32_t z5 = (z3 + z4) * FIX_1_175875602;
        tmp4 *= FIX_0_298631336; tmp5 *= FIX_2_053119869;
        tmp6 *= FIX_3_072711026; tmp7 *= FIX_1_501321110;
        z1 *= -FIX_0_899976223; z2 *= -FIX_2_562915447;
        z3 *= -FIX_1_961570560; z4 *= -FIX_0_390180644;
        z3 += z5; z4 += z5;
        o[7] = descale(tmp4 + z1 + z3, kConstBits - kPass1Bits);
        o[5] = descale(tmp5 + z2 + z4, kConstBits - kPass1Bits);
        o[3] = descale(tmp6 + z2 + z3, kConstBits - kPass1Bits);
        o[1] = descale(tmp7 + z1 + z4, kConstBits - kPass1Bits);
    }
    for (int x = 0; x < 8; ++x) {
        const int32_t* c = ws + x;
        int32_t* o = out + x;
        int32_t tmp0 = c[0] + c[56], tmp7 = c[0] - c[56];
        int32_t tmp1 = c[8] + c[48], tmp6 = c[8] - c[48];
        int32_t tmp2 = c[16] + c[40], tmp5 = c[16] - c[40];
        int32_t tmp3 = c[24] + c[32], tmp4 = c[24] - c[32];

        int32_t tmp10 = tmp0 + tmp3, tmp13 = tmp0 - tmp3;
        int32_t tmp11 = tmp1 + tmp2, tmp12 = tmp1 - tmp2;
        o[0] = descale(tmp10 + tmp11, kPass1Bits);
        o[32] = descale(tmp10 - tmp11, kPass1Bits);
        int32_t z1 = (tmp12 + tmp13) * FIX_0_541196100;
        o[16] = descale(z1 + tmp13 * FIX_0_765366865, kConstBits + kPass1Bits);
        o[48] = descale(z1 - tmp12 * FIX_1_847759065, kConstBits + kPass1Bits);

        z1 = tmp4 + tmp7;
        int32_t z2 = tmp5 + tmp6, z3 = tmp4 + tmp6, z4 = tmp5 + tmp7;
        int32_t z5 = (z3 + z4) * FIX_1_175875602;
        tmp4 *= FIX_0_298631336; tmp5 *= FIX_2_053119869;
        tmp6 *= FIX_3_072711026; tmp7 *= FIX_1_501321110;
        z1 *= -FIX_0_899976223; z2 *= -FIX_2_562915447;
        z3 *= -FIX_1_961570560; z4 *= -FIX_0_390180644;
        z3 += z5; z4 += z5;
        o[56] = descale(tmp4 + z1 + z3, kConstBits + kPass1Bits);
        o[40] = descale(tmp5 + z2 + z4, kConstBits + kPass1Bits);
        o[24] = descale(tmp6 + z2 + z3, kConstBits + kPass1Bits);
        o[8] = descale(tmp7 + z1 + z4, kConstBits + kPass1Bits);
    }
}

void idct8x8(const int32_t in[64], uint8_t* dst, int stride) {
    int32_t ws[64];
    // Pass 1: columns. Inputs within 11 bits + sign keep every product in 32 bits.
    for (int x = 0; x < 8; ++x) {
        const int32_t* c = in + x;
        int32_t* w = ws + x;
        if ((c[8] | c[16] | c[24] | c[32] | c[40] | c[48] | c[56]) == 0) {
            int32_t dc = c[0] * (1 << kPass1Bits);
            w[0] = w[8] = w[16] = w[24] = w[32] = w[40] = w[48] = w[56] = dc;
            continue;
        }
        int32_t z2 = c[16], z3 = c[48];
        int32_t z1 = (z2 + z3) * FIX_0_541196100;
        int32_t tmp2 = z1 - z3 * FIX_1_847759065;
        int32_t tmp3 = z1 + z2 * FIX_0_765366865;
        z2 = c[0];
        z3 = c[32];
        int32_t tmp0 = (z2 + z3) * (1 << kConstBits);
        int32_t tmp1 = (z2 - z3) * (1 << kConstBits);
        int32_t tmp10 = tmp0 + tmp3, tmp13 = tmp0 - tmp3;
        int32_t tmp11 = tmp1 + tmp2, tmp12 = tmp1 - tmp2;

        tmp0 = c[56]; tmp1 = c[40]; tmp2 = c[24]; tmp3 = c[8];
        z1 = tmp0 + tmp3;
        z2 = tmp1 + tmp2;
        z3 = tmp0 + tmp2;
        int32_t z4 = tmp1 + tmp3;
        int32_t z5 = (z3 + z4) * FIX_1_175875602;
        tmp0 *= FIX_0_298631336; tmp1 *= FIX_2_053119869;
        tmp2 *= FIX_3_072711026; tmp3 *= FIX_1_501321110;
        z1 *= -FIX_0_899976223; z2 *= -FIX_2_562915447;
        z3 *= -FIX_1_961570560; z4 *= -FIX_0_390180644;
        z3 += z5; z4 += z5;
        tmp0 += z1 + z3; tmp1 += z2 + z4; tmp2 += z2 + z3; tmp3 += z1 + z4;

        const int sh = kConstBits - kPass1Bits;
        w[0] = descale(tmp10 + tmp3, sh);
        w[56] = descale(tmp10 - tmp3, sh);
        w[8] = descale(tmp11 + tmp2, sh);
        w[48] = descale(tmp11 - tmp2, sh);
        w[16] = descale(tmp12 + tmp1, sh);
        w[40] = descale(tmp12 - tmp1, sh);
        w[24] = descale(tmp13 + tmp0, sh);
        w[32] = descale(tmp13 - tmp0, sh);
    }
    // Pass 2: rows, removing the remaining scale (incl. the factor of 8) and level shift.
    const int sh = kConstBits + kPass1Bits + 3;
    for (int y = 0; y < 8; ++y) {
        const int32_t* r = ws + y * 8;
        uint8_t* o = dst + size_t(y) * stride;
        if ((r[1] | r[2] | r[3] | r[4] | r[5] | r[6] | r[7]) == 0) {
            uint8_t v = clamp_u8(descale(r[0] * (1 << kConstBits), sh) + 128);
            std::memset(o, v, 8);
            continue;
        }
        int32_t z2 = r[2], z3 = r[6];
        int32_t z1 = (z2 + z3) * FIX_0_541196100;
        int32_t tmp2 = z1 - z3 * FIX_1_847759065;
        int32_t tmp3 = z1 + z2 * FIX_0_765366865;
        int32_t tmp0 = (r[0] + r[4]) * (1 << kConstBits);
        int32_t tmp1 = (r[0] - r[4]) * (1 << kConstBits);
        int32_t tmp10 = tmp0 + tmp3, tmp13 = tmp0 - tmp3;
        int32_t tmp11 = tmp1 + tmp2, tmp12 = tmp1 - tmp2;

        tmp0 = r[7]; tmp1 = r[5]; tmp2 = r[3]; tmp3 = r[1];
        z1 = tmp0 + tmp3;
        z2 = tmp1 + tmp2;
        z3 = tmp0 + tmp2;
        int32_t z4 = tmp1 + tmp3;
        int32_t z5 = (z3 + z4) * FIX_1_175875602;
        tmp0 *= FIX_0_298631336; tmp1 *= FIX_2_053119869;
        tmp2 *= FIX_3_072711026; tmp3 *= FIX_1_501321110;
        z1 *= -FIX_0_899976223; z2 *= -FIX_2_562915447;
        z3 *= -FIX_1_961570560; z4 *= -FIX_0_390180644;
        z3 += z5; z4 += z5;
        tmp0 += z1 + z3; tmp1 += z2 + z4; tmp2 += z2 + z3; tmp3 += z1 + z4;

        o[0] = clamp_u8(descale(tmp10 + tmp3, sh) + 128);
        o[7] = clamp_u8(descale(tmp10 - tmp3, sh) + 128);
        o[1] = clamp_u8(descale(tmp11 + tmp2, sh) + 128);
        o[6] = clamp_u8(descale(tmp11 - tmp2, sh) + 128);
        o[2] = clamp_u8(descale(tmp12 + tmp1, sh) + 128);
        o[5] = clamp_u8(descale(tmp12 - tmp1, sh) + 128);
        o[3] = clamp_u8(descale(tmp13 + tmp0, sh) + 128);
        o[4] = clamp_u8(descale(tmp13 - tmp0, sh) + 128);
    }
}

void idct8x8_dc(int32_t dc, uint8_t* dst, int stride) {
    // Same arithmetic idct8x8 performs when only coefficient 0 is non-zero.
    int32_t ws = dc * (1 << kPass1Bits);
    uint8_t v = clamp_u8(descale(ws * (1 << kConstBits), kConstBits + kPass1Bits + 3) + 128);
    for (int y = 0; y < 8; ++y) std::memset(dst + size_t(y) * stride, v, 8);
}

void build_quant_tables(int quality, uint16_t luma[64], uint16_t chroma[64]) {
    if (quality < 1) quality = 1;
    if (quality > 100) quality = 100;
    int scale = quality < 50 ? 5000 / quality : 200 - quality * 2;
    for (int i = 0; i < 64; ++i) {
        int l = (kStdLumaQ[i] * scale + 50) / 100;
        int c = (kStdChromaQ[i] * scale + 50) / 100;
        luma[i] = uint16_t(std::min(255, std::max(1, l)));
        chroma[i] = uint16_t(std::min(255, std::max(1, c)));
    }
}

const HuffTables& huff_tables() {
    static const HuffTables tables = [] {
        HuffTables t;
        build_huff(kDcLumaBits, kDcLumaVals, t.enc_dc[0], t.dec_dc[0]);
        build_huff(kDcChromaBits, kDcChromaVals, t.enc_dc[1], t.dec_dc[1]);
        build_huff(kAcLumaBits, kAcLumaVals, t.enc_ac[0], t.dec_ac[0]);
        build_huff(kAcChromaBits, kAcChromaVals, t.enc_ac[1], t.dec_ac[1]);
        return t;
    }();
    return tables;
}

void BitWriter::put(uint32_t bits, int len) {
    if (len == 0) return;
    acc_ = (acc_ << len) | (bits & ((1u << len) - 1));
    n_ += len;
    while (n_ >= 8) {
        n_ -= 8;
        out_.push_back(uint8_t(acc_ >> n_));
    }
}

void BitWriter::flush() {
    if (n_ > 0) out_.push_back(uint8_t(acc_ << (8 - n_)));
    acc_ = 0;
    n_ = 0;
}

int huff_decode(BitReader& br, const HuffDec& t) {
    br.refill();
    uint32_t e = t.lut[br.peek(9)];
    if (e) {
        br.skip(int(e >> 8));
        return int(e & 0xff);
    }
    for (int len = 10; len <= 16; ++len) {
        int32_t code = int32_t(br.peek(len));
        if (code <= t.maxcode[len]) {
            br.skip(len);
            return t.vals[code + t.valoff[len]];
        }
    }
    return -1;
}

void encode_block(BitWriter& bw, const int16_t zz[64], int& dc_pred, int table) {
    const HuffTables& ht = huff_tables();
    const HuffEnc& dc = ht.enc_dc[table];
    const HuffEnc& ac = ht.enc_ac[table];

    int diff = zz[0] - dc_pred;
    dc_pred = zz[0];
    int s = bit_count(diff < 0 ? -diff : diff);
    bw.put(dc.code[s], dc.size[s]);
    bw.put(uint32_t(diff < 0 ? diff - 1 : diff), s);

    int run = 0;
    for (int k = 1; k < 64; ++k) {
        int v = zz[k];
        if (v == 0) {
            ++run;
            continue;
        }
        while (run > 15) {
            bw.put(ac.code[0xf0], ac.size[0xf0]);
            run -= 16;
        }
        s = bit_count(v < 0 ? -v : v);
        int sym = (run << 4) | s;
        bw.put(ac.code[sym], ac.size[sym]);
        bw.put(uint32_t(v < 0 ? v - 1 : v), s);
        run = 0;
    }
    if (run) bw.put(ac.code[0x00], ac.size[0x00]);
}

int decode_block(BitReader& br, int table, const uint16_t dq[64], int& dc_pred, int32_t coef[64]) {
    const HuffTables& ht = huff_tables();
    std::memset(coef, 0, 64 * sizeof(int32_t));

    int s = huff_decode(br, ht.dec_dc[table]);
    if (s < 0 || s > 11) return -1;
    int diff = 0;
    if (s) {
        diff = int(br.get(s));
        if (diff < (1 << (s - 1))) diff -= (1 << s) - 1;
    }
    dc_pred += diff;
    if (dc_pred < -2047 || dc_pred > 2047) return -1;
    int32_t v = dc_pred * int32_t(dq[0]);
    coef[0] = std::min<int32_t>(2047, std::max<int32_t>(-2048, v));

    const HuffDec& ac = ht.dec_ac[table];
    int has_ac = 0;
    for (int k = 1; k < 64;) {
        int rs = huff_decode(br, ac);
        if (rs < 0) return -1;
        int r = rs >> 4;
        s = rs & 15;
        if (s == 0) {
            if (r != 15) break;  // EOB
            k += 16;             // ZRL
            continue;
        }
        k += r;
        if (k > 63) return -1;
        int a = int(br.get(s));
        if (a < (1 << (s - 1))) a -= (1 << s) - 1;
        v = a * int32_t(dq[k]);
        coef[kZigzag[k]] = std::min<int32_t>(2047, std::max<int32_t>(-2048, v));
        has_ac = 1;
        ++k;
    }
    return has_ac;
}

void pack_bitmap(const uint8_t* flags, int n, uint8_t* out) {
    std::memset(out, 0, size_t(n + 7) / 8);
    for (int i = 0; i < n; ++i)
        if (flags[i]) out[i >> 3] |= uint8_t(1u << (i & 7));
}

void unpack_bitmap(const uint8_t* in, int n, uint8_t* flags) {
    for (int i = 0; i < n; ++i) flags[i] = (in[i >> 3] >> (i & 7)) & 1;
}

uint64_t sad(const uint8_t* a, int astride, const uint8_t* b, int bstride, int w, int h, uint64_t limit) {
    uint64_t total = 0;
    for (int y = 0; y < h; ++y) {
        const uint8_t* pa = a + size_t(y) * astride;
        const uint8_t* pb = b + size_t(y) * bstride;
        uint32_t row = 0;
        for (int x = 0; x < w; ++x) {
            int d = int(pa[x]) - int(pb[x]);
            row += uint32_t(d < 0 ? -d : d);
        }
        total += row;
        if (total > limit) return total;
    }
    return total;
}

void RefreshPolicy::reset(RefreshMode mode, uint32_t param, int tiles) {
    mode_ = mode;
    param_ = param;
    tiles_ = tiles;
    cursor_ = 0;
}

bool RefreshPolicy::forced(uint32_t frame_index) const {
    if (frame_index == 0) return true;
    return mode_ == RefreshMode::FullPeriodic && param_ && frame_index % param_ == 0;
}

void RefreshPolicy::apply(uint8_t* dirty) {
    if (mode_ != RefreshMode::Rolling || param_ == 0 || tiles_ == 0) return;
    uint32_t n = std::min<uint32_t>(param_, uint32_t(tiles_));
    for (uint32_t i = 0; i < n; ++i) {
        dirty[cursor_] = 1;
        if (++cursor_ == tiles_) cursor_ = 0;
    }
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Decoder
// ---------------------------------------------------------------------------
Status Decoder::read_header(ReadFn fn, void* user) {
    ready_ = false;
    uint8_t buf[kStreamHeaderSize];
    size_t n = read_full(fn, user, buf, sizeof(buf));
    if (n == 0) return Status::EndOfStream;
    if (n < sizeof(buf)) return Status::Truncated;
    Status st = parse_stream_header(buf, &header_);
    if (st != Status::Ok) return st;
    if (uint64_t(header_.width) * header_.height > uint64_t(TJC_MAX_PIXELS)) return Status::BadHeader;

    layout_ = make_layout(header_.width, header_.height, header_.tile_w, header_.tile_h);
    uint16_t ql[64], qc[64];
    detail::build_quant_tables(header_.quality, ql, qc);
    for (int k = 0; k < 64; ++k) {
        dq_[0][k] = ql[detail::kZigzag[k]];
        dq_[1][k] = qc[detail::kZigzag[k]];
    }
    if (!fb_.alloc(layout_)) return Status::OutOfMemory;
    try {
        dirty_.assign(size_t(layout_.tiles), 0);
        bitmap_.assign(size_t(layout_.tiles + 7) / 8, 0);
    } catch (...) {
        return Status::OutOfMemory;
    }
    ready_ = true;
    return Status::Ok;
}

bool Decoder::decode_tile(const uint8_t* data, size_t size, int tile) {
    detail::BitReader br(data, size);
    int32_t coef[64];
    for (int c = 0; c < 3; ++c) {
        int x0, y0, w, h;
        tile_rect(layout_, tile, c, x0, y0, w, h);
        int stride = fb_.stride[c];
        uint8_t* base = fb_.p[c].data() + size_t(y0) * stride + x0;
        int table = c ? 1 : 0;
        int pred = 0;
        for (int by = 0; by < h; by += 8) {
            for (int bx = 0; bx < w; bx += 8) {
                int r = detail::decode_block(br, table, dq_[table], pred, coef);
                if (r < 0) return false;
                uint8_t* dst = base + size_t(by) * stride + bx;
                if (r) detail::idct8x8(coef, dst, stride);
                else detail::idct8x8_dc(coef[0], dst, stride);
            }
        }
    }
    return !br.overrun();
}

Status Decoder::decode_frame(ReadFn fn, void* user, FrameStats* stats) {
    if (!ready_) return Status::NotInitialized;
    uint8_t hdr[kFrameHeaderSize];
    size_t n = read_full(fn, user, hdr, sizeof(hdr));
    if (n == 0) return Status::EndOfStream;
    if (n < sizeof(hdr)) return Status::Truncated;

    uint32_t frame_num = get_u32(hdr);
    if (hdr[4] > 1) return Status::Corrupt;
    bool force = hdr[4] != 0;
    size_t bytes = kFrameHeaderSize;
    if (force) {
        std::fill(dirty_.begin(), dirty_.end(), uint8_t(1));
    } else {
        if (read_full(fn, user, bitmap_.data(), bitmap_.size()) < bitmap_.size()) return Status::Truncated;
        detail::unpack_bitmap(bitmap_.data(), layout_.tiles, dirty_.data());
        bytes += bitmap_.size();
    }

    uint32_t dirty_count = 0;
    for (int t = 0; t < layout_.tiles; ++t) {
        if (!dirty_[size_t(t)]) continue;
        uint8_t lb[2];
        if (read_full(fn, user, lb, 2) < 2) return Status::Truncated;
        size_t len = get_u16(lb);
        if (len == 0) return Status::Corrupt;
        if (payload_.size() < len) {
            try {
                payload_.resize(len);
            } catch (...) {
                return Status::OutOfMemory;
            }
        }
        if (read_full(fn, user, payload_.data(), len) < len) return Status::Truncated;
        if (!decode_tile(payload_.data(), len, t)) return Status::Corrupt;
        bytes += 2 + len;
        ++dirty_count;
    }

    if (stats) {
        stats->frame_num = frame_num;
        stats->force_refresh = force;
        stats->dirty_tiles = dirty_count;
        stats->total_tiles = uint32_t(layout_.tiles);
        stats->bytes = bytes;
    }
    return Status::Ok;
}

void Decoder::copy_frame(uint8_t* dst) const { copy_cropped(fb_, layout_, dst); }

// ---------------------------------------------------------------------------
// Encoder
// ---------------------------------------------------------------------------
#ifndef TJC_NO_ENCODER
bool Encoder::init(const Config& cfg) {
    ready_ = false;
    if (const char* e = validate_geometry(cfg.width, cfg.height, cfg.tile_w, cfg.tile_h)) {
        error_ = e;
        return false;
    }
    if (cfg.quality < 1 || cfg.quality > 100) {
        error_ = "quality must be 1..100";
        return false;
    }
    if (cfg.refresh_param > 65535) {
        error_ = "refresh_param must be 0..65535";
        return false;
    }
    if (uint8_t(cfg.refresh_mode) > 2 || uint8_t(cfg.diff_reference) > 1) {
        error_ = "invalid refresh mode or diff reference";
        return false;
    }
    cfg_ = cfg;
    layout_ = make_layout(cfg.width, cfg.height, cfg.tile_w, cfg.tile_h);
    if (!cur_.alloc(layout_) || !last_coded_.alloc(layout_) || !recon_.alloc(layout_)) {
        error_ = "out of memory";
        return false;
    }
    dirty_.assign(size_t(layout_.tiles), 0);
    detail::build_quant_tables(cfg.quality, q_[0], q_[1]);
    policy_.reset(cfg.refresh_mode, cfg.refresh_param, layout_.tiles);
    frame_num_ = 0;
    frame_index_ = 0;
    keyframe_request_ = false;
    error_ = nullptr;
    ready_ = true;
    return true;
}

StreamHeader Encoder::stream_header() const {
    StreamHeader h;
    h.width = cfg_.width;
    h.height = cfg_.height;
    h.tile_w = cfg_.tile_w;
    h.tile_h = cfg_.tile_h;
    h.chroma_format = 0;
    h.refresh_mode = cfg_.refresh_mode;
    h.refresh_param = uint16_t(cfg_.refresh_param);
    h.quality = cfg_.quality;
    return h;
}

void Encoder::write_stream_header(std::vector<uint8_t>& out) const {
    uint8_t buf[kStreamHeaderSize];
    tjc::write_stream_header(stream_header(), buf);
    out.insert(out.end(), buf, buf + kStreamHeaderSize);
}

void Encoder::load(const uint8_t* y, int ys, const uint8_t* u, int us, const uint8_t* v, int vs) {
    const uint8_t* src[3] = {y, u, v};
    int sstride[3] = {ys, us, vs};
    for (int c = 0; c < 3; ++c) {
        int w, h;
        visible_size(layout_, c, w, h);
        int pw = cur_.stride[c];
        int ph = c ? layout_.pcheight : layout_.pheight;
        uint8_t* d = cur_.p[c].data();
        for (int r = 0; r < h; ++r) {
            uint8_t* row = d + size_t(r) * pw;
            std::memcpy(row, src[c] + size_t(r) * sstride[c], size_t(w));
            std::memset(row + w, row[w - 1], size_t(pw - w));
        }
        for (int r = h; r < ph; ++r) std::memcpy(d + size_t(r) * pw, d + size_t(h - 1) * pw, size_t(pw));
    }
}

bool Encoder::tile_changed(int tile) const {
    // Judged per 8x8 block so a small moving object is not averaged away in a big tile.
    const Planes& ref = cfg_.diff_reference == DiffReference::Reconstructed ? recon_ : last_coded_;
    uint64_t limit = uint64_t(cfg_.motion_threshold_k) * 64;
    for (int c = 0; c < 3; ++c) {
        int x, y, w, h;
        tile_rect(layout_, tile, c, x, y, w, h);
        int s = cur_.stride[c];
        for (int by = 0; by < h; by += 8) {
            for (int bx = 0; bx < w; bx += 8) {
                size_t off = size_t(y + by) * s + x + bx;
                if (detail::sad(cur_.p[c].data() + off, s, ref.p[c].data() + off, s, 8, 8, limit) > limit)
                    return true;
            }
        }
    }
    return false;
}

void Encoder::copy_tile(const Planes& from, Planes& to, int tile) {
    for (int c = 0; c < 3; ++c) {
        int x, y, w, h;
        tile_rect(layout_, tile, c, x, y, w, h);
        int s = from.stride[c];
        for (int r = 0; r < h; ++r) {
            size_t off = size_t(y + r) * s + x;
            std::memcpy(to.p[c].data() + off, from.p[c].data() + off, size_t(w));
        }
    }
}

void Encoder::encode_tile(int tile, std::vector<uint8_t>& out) {
    size_t len_pos = out.size();
    out.push_back(0);
    out.push_back(0);
    detail::BitWriter bw(out);
    int32_t dct[64], deq[64];
    int16_t zz[64];
    for (int c = 0; c < 3; ++c) {
        int x0, y0, w, h;
        tile_rect(layout_, tile, c, x0, y0, w, h);
        int stride = cur_.stride[c];
        size_t base = size_t(y0) * stride + x0;
        int table = c ? 1 : 0;
        const uint16_t* q = q_[table];
        int pred = 0;
        for (int by = 0; by < h; by += 8) {
            for (int bx = 0; bx < w; bx += 8) {
                size_t off = base + size_t(by) * stride + bx;
                detail::fdct8x8(cur_.p[c].data() + off, stride, dct);
                bool has_ac = false;
                for (int k = 0; k < 64; ++k) {
                    int nat = detail::kZigzag[k];
                    int32_t d = int32_t(q[nat]) * 8;
                    int32_t v = dct[nat];
                    int32_t qv = v < 0 ? -((-v + d / 2) / d) : (v + d / 2) / d;
                    qv = std::min<int32_t>(1023, std::max<int32_t>(-1023, qv));
                    zz[k] = int16_t(qv);
                    // Shadow decode: same dequantization and clamp as decode_block().
                    deq[nat] = std::min<int32_t>(2047, std::max<int32_t>(-2048, qv * int32_t(q[nat])));
                    if (k && qv) has_ac = true;
                }
                detail::encode_block(bw, zz, pred, table);
                uint8_t* dst = recon_.p[c].data() + off;
                if (has_ac) detail::idct8x8(deq, dst, stride);
                else detail::idct8x8_dc(deq[0], dst, stride);
            }
        }
    }
    bw.flush();
    put_u16(out.data() + len_pos, uint32_t(out.size() - len_pos - 2));
}

void Encoder::encode_frame(const uint8_t* yuv420p, std::vector<uint8_t>& out, FrameStats* stats) {
    size_t ys = size_t(cfg_.width) * cfg_.height;
    size_t cs = size_t(layout_.cwidth) * layout_.cheight;
    encode_frame(yuv420p, cfg_.width, yuv420p + ys, layout_.cwidth, yuv420p + ys + cs, layout_.cwidth, out,
                 stats);
}

void Encoder::encode_frame(const uint8_t* y, int ystride, const uint8_t* u, int ustride, const uint8_t* v,
                           int vstride, std::vector<uint8_t>& out, FrameStats* stats) {
    if (!ready_) return;
    size_t start = out.size();
    load(y, ystride, u, ustride, v, vstride);

    bool force = keyframe_request_ || policy_.forced(frame_index_);
    keyframe_request_ = false;
    int tiles = layout_.tiles;
    if (force) {
        std::fill(dirty_.begin(), dirty_.end(), uint8_t(1));
    } else {
        for (int t = 0; t < tiles; ++t) dirty_[size_t(t)] = tile_changed(t) ? 1 : 0;
        policy_.apply(dirty_.data());
    }

    uint8_t hdr[kFrameHeaderSize];
    put_u32(hdr, frame_num_);
    hdr[4] = force ? 1 : 0;
    out.insert(out.end(), hdr, hdr + kFrameHeaderSize);
    if (!force) {
        size_t pos = out.size();
        out.resize(pos + size_t(tiles + 7) / 8);
        detail::pack_bitmap(dirty_.data(), tiles, out.data() + pos);
    }

    uint32_t dirty_count = 0;
    for (int t = 0; t < tiles; ++t) {
        if (!dirty_[size_t(t)]) continue;
        encode_tile(t, out);
        copy_tile(cur_, last_coded_, t);
        ++dirty_count;
    }

    if (stats) {
        stats->frame_num = frame_num_;
        stats->force_refresh = force;
        stats->dirty_tiles = dirty_count;
        stats->total_tiles = uint32_t(tiles);
        stats->bytes = out.size() - start;
    }
    ++frame_num_;
    ++frame_index_;
}

void Encoder::copy_recon(uint8_t* dst) const { copy_cropped(recon_, layout_, dst); }
#endif  // TJC_NO_ENCODER

}  // namespace tjc

#endif  // TJC_IMPLEMENTATION_DONE
#endif  // TJC_IMPLEMENTATION
