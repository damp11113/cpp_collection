// Unit and end-to-end tests for tjc.h.
// Build: g++ -std=c++17 -O2 -o tjc_test tjc_test.cpp && ./tjc_test

#define TJC_IMPLEMENTATION
#include "tjc.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <set>
#include <string>
#include <vector>

static int g_failed = 0;
static int g_checks = 0;

#define CHECK(cond)                                                              \
    do {                                                                         \
        ++g_checks;                                                              \
        if (!(cond)) {                                                           \
            ++g_failed;                                                          \
            std::printf("  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);        \
        }                                                                        \
    } while (0)

using namespace tjc;

static const double kPi = 3.14159265358979323846;

static uint32_t g_rng = 12345;
static uint32_t rnd() {
    g_rng = g_rng * 1664525u + 1013904223u;
    return g_rng >> 8;
}

static void ref_dct(const uint8_t* src, double out[64]) {
    for (int v = 0; v < 8; ++v)
        for (int u = 0; u < 8; ++u) {
            double s = 0;
            for (int y = 0; y < 8; ++y)
                for (int x = 0; x < 8; ++x)
                    s += (src[y * 8 + x] - 128.0) * std::cos((2 * x + 1) * u * kPi / 16) *
                         std::cos((2 * y + 1) * v * kPi / 16);
            double cu = u ? 1 : std::sqrt(0.5), cv = v ? 1 : std::sqrt(0.5);
            out[v * 8 + u] = 0.25 * cu * cv * s;
        }
}

static void test_bitstream() {
    std::printf("bitstream\n");
    std::vector<uint8_t> buf;
    detail::BitWriter bw(buf);
    std::vector<std::pair<uint32_t, int>> items;
    for (int i = 0; i < 5000; ++i) {
        int len = int(rnd() % 17);
        uint32_t v = len ? rnd() & ((1u << len) - 1) : 0;
        items.push_back({v, len});
        bw.put(v, len);
    }
    bw.flush();
    detail::BitReader br(buf.data(), buf.size());
    bool ok = true;
    for (auto& it : items) ok &= br.get(it.second) == it.first;
    CHECK(ok);
    CHECK(!br.overrun());
    br.get(16);
    br.get(16);
    CHECK(br.overrun());
}

static void test_dct() {
    std::printf("dct\n");
    int max_fwd_err = 0, max_rt_err = 0;
    bool fast_path_exact = true;
    for (int iter = 0; iter < 2000; ++iter) {
        uint8_t blk[64];
        int kind = iter % 4;
        for (int i = 0; i < 64; ++i) {
            if (kind == 0) blk[i] = uint8_t(rnd());
            else if (kind == 1) blk[i] = uint8_t((i % 8) * 30 + (i / 8) * 2);
            else if (kind == 2) blk[i] = (i + iter) % 2 ? 255 : 0;
            else blk[i] = uint8_t(iter & 255);
        }
        int32_t f[64];
        double r[64];
        detail::fdct8x8(blk, 8, f);
        ref_dct(blk, r);
        for (int i = 0; i < 64; ++i) {
            int e = int(std::lround(std::fabs(f[i] / 8.0 - r[i])));
            if (e > max_fwd_err) max_fwd_err = e;
        }
        // Unquantized round trip (q = 1).
        int32_t c[64];
        for (int i = 0; i < 64; ++i) {
            int32_t v = f[i];
            c[i] = v < 0 ? -((-v + 4) / 8) : (v + 4) / 8;
        }
        uint8_t out[64];
        detail::idct8x8(c, out, 8);
        for (int i = 0; i < 64; ++i) {
            int e = std::abs(int(out[i]) - int(blk[i]));
            if (e > max_rt_err) max_rt_err = e;
        }
        // DC-only fast path must match the full IDCT bit for bit.
        int32_t dc[64] = {};
        dc[0] = int32_t(rnd() % 4096) - 2048;
        uint8_t a[64], b[64];
        detail::idct8x8(dc, a, 8);
        detail::idct8x8_dc(dc[0], b, 8);
        fast_path_exact &= std::memcmp(a, b, 64) == 0;
    }
    std::printf("  max |fdct - float dct| = %d, max round-trip pixel error = %d\n", max_fwd_err, max_rt_err);
    CHECK(max_fwd_err <= 1);
    CHECK(max_rt_err <= 1);
    CHECK(fast_path_exact);
}

static void test_huffman_tables() {
    std::printf("huffman tables\n");
    const detail::HuffTables& ht = detail::huff_tables();
    for (int t = 0; t < 2; ++t) {
        for (int s = 0; s <= 11; ++s) CHECK(ht.enc_dc[t].size[s] > 0);
        // Every (run, size) symbol baseline coding can emit must be present.
        int missing = 0;
        for (int run = 0; run < 16; ++run)
            for (int size = 1; size <= 10; ++size) missing += ht.enc_ac[t].size[(run << 4) | size] == 0;
        missing += ht.enc_ac[t].size[0x00] == 0;
        missing += ht.enc_ac[t].size[0xf0] == 0;
        CHECK(missing == 0);
        // Codes must be prefix-free: decode each code on its own.
        for (int sym = 0; sym < 256; ++sym) {
            for (int dc = 0; dc < 2; ++dc) {
                const detail::HuffEnc& e = dc ? ht.enc_dc[t] : ht.enc_ac[t];
                if (!e.size[sym]) continue;
                std::vector<uint8_t> buf;
                detail::BitWriter bw(buf);
                bw.put(e.code[sym], e.size[sym]);
                bw.flush();
                detail::BitReader br(buf.data(), buf.size());
                CHECK(detail::huff_decode(br, dc ? ht.dec_dc[t] : ht.dec_ac[t]) == sym);
            }
        }
    }
}

static void test_block_roundtrip() {
    std::printf("block entropy round trip\n");
    uint16_t ones[64];
    for (auto& v : ones) v = 1;
    for (int table = 0; table < 2; ++table) {
        std::vector<uint8_t> buf;
        detail::BitWriter bw(buf);
        std::vector<std::vector<int16_t>> blocks;
        int pred = 0;
        for (int b = 0; b < 500; ++b) {
            std::vector<int16_t> zz(64, 0);
            int density = int(rnd() % 4);
            for (int k = 0; k < 64; ++k) {
                bool nz = density == 3 || (rnd() % (k + 1 + density * 8)) < 2;
                if (k == 0 || nz) zz[size_t(k)] = int16_t(int(rnd() % 2047) - 1023);
            }
            if (b % 50 == 0) std::fill(zz.begin() + 1, zz.end(), int16_t(0));
            if (b % 77 == 0) { std::fill(zz.begin() + 1, zz.end(), int16_t(0)); zz[63] = -1; }  // long ZRL run
            blocks.push_back(zz);
            detail::encode_block(bw, zz.data(), pred, table);
        }
        bw.flush();
        detail::BitReader br(buf.data(), buf.size());
        pred = 0;
        bool ok = true;
        for (auto& zz : blocks) {
            int32_t coef[64];
            int r = detail::decode_block(br, table, ones, pred, coef);
            ok &= r >= 0;
            bool any_ac = false;
            for (int k = 0; k < 64; ++k) {
                int32_t expect = zz[size_t(k)];
                ok &= coef[detail::kZigzag[k]] == expect;
                if (k && expect) any_ac = true;
            }
            ok &= (r == 1) == any_ac;
        }
        CHECK(ok);
        CHECK(!br.overrun());
    }
}

static void test_quant_div() {
    std::printf("reciprocal quantization\n");
    bool ok = true;
    for (uint32_t q = 1; q <= 255 && ok; ++q) {
        uint32_t d = q * 8;
        uint64_t r = detail::quant_recip(d);
        for (int32_t v = -(1 << 16); v <= (1 << 16); ++v) {
            int32_t expect = v < 0 ? -((-v + int32_t(d) / 2) / int32_t(d)) : (v + int32_t(d) / 2) / int32_t(d);
            if (detail::quant_div(v, d, r) != expect) {
                std::printf("  mismatch v=%d d=%u\n", v, d);
                ok = false;
                break;
            }
        }
    }
    CHECK(ok);
}

static void test_bitmap() {
    std::printf("bitmap packing\n");
    for (int n : {1, 7, 8, 9, 300, 1021}) {
        std::vector<uint8_t> flags(static_cast<size_t>(n)), back(static_cast<size_t>(n));
        for (auto& f : flags) f = rnd() & 1;
        std::vector<uint8_t> packed(size_t(n + 7) / 8);
        detail::pack_bitmap(flags.data(), n, packed.data());
        detail::unpack_bitmap(packed.data(), n, back.data());
        CHECK(flags == back);
    }
    // Bit order is fixed by the spec: tile i -> bit (i & 7) of byte i >> 3.
    uint8_t flags[10] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    uint8_t packed[2];
    detail::pack_bitmap(flags, 10, packed);
    CHECK(packed[0] == 0x01 && packed[1] == 0x02);
}

static void test_sad() {
    std::printf("sad / threshold\n");
    std::vector<uint8_t> a(32 * 32), b(32 * 32);
    for (auto& v : a) v = uint8_t(rnd());
    b = a;
    CHECK(detail::sad(a.data(), 32, b.data(), 32, 32, 32, ~0ull) == 0);
    std::vector<uint8_t> z(32 * 32, 0), f(32 * 32, 255);
    CHECK(detail::sad(z.data(), 32, f.data(), 32, 32, 32, ~0ull) == 255ull * 32 * 32);
    CHECK(detail::sad(z.data(), 32, f.data(), 32, 32, 32, 100) > 100);

    // Threshold semantics through the encoder: mean |diff| == K stays clean, K+1 is dirty.
    Config cfg;
    cfg.width = 32;
    cfg.height = 32;
    cfg.tile_w = cfg.tile_h = 16;
    cfg.motion_threshold_k = 3;
    cfg.refresh_mode = RefreshMode::None;
    Encoder enc;
    CHECK(enc.init(cfg));
    std::vector<uint8_t> frame(frame_size(32, 32), 100), out;
    FrameStats st;
    enc.encode_frame(frame.data(), out, &st);
    CHECK(st.force_refresh && st.dirty_tiles == 4);
    enc.encode_frame(frame.data(), out, &st);
    CHECK(!st.force_refresh && st.dirty_tiles == 0);
    std::vector<uint8_t> f3(frame.size(), 103), f4(frame.size(), 104);
    enc.encode_frame(f3.data(), out, &st);
    CHECK(st.dirty_tiles == 0);
    enc.encode_frame(f4.data(), out, &st);
    CHECK(st.dirty_tiles == 4);
    // Only the top-left tile changes.
    std::vector<uint8_t> f5 = f4;
    for (int y = 0; y < 16; ++y)
        for (int x = 0; x < 16; ++x) f5[size_t(y) * 32 + x] = 200;
    enc.encode_frame(f5.data(), out, &st);
    CHECK(st.dirty_tiles == 1);
}

static void test_refresh_policy() {
    std::printf("refresh policy\n");
    detail::RefreshPolicy p;
    p.reset(RefreshMode::None, 30, 100);
    CHECK(p.forced(0));
    CHECK(!p.forced(1) && !p.forced(30));

    p.reset(RefreshMode::FullPeriodic, 10, 100);
    CHECK(p.forced(0) && p.forced(10) && p.forced(20));
    CHECK(!p.forced(1) && !p.forced(9) && !p.forced(11));

    // Rolling: every tile refreshed at least once within ceil(tiles / param) frames.
    for (int tiles : {1, 7, 100, 301}) {
        for (uint32_t param : {1u, 3u, 30u, 500u}) {
            p.reset(RefreshMode::Rolling, param, tiles);
            CHECK(!p.forced(1));
            int bound = int((uint32_t(tiles) + param - 1) / param);
            std::vector<int> last(size_t(tiles), -1);
            bool ok = true;
            for (int frame = 0; frame < bound * 4; ++frame) {
                std::vector<uint8_t> dirty(size_t(tiles), 0);
                p.apply(dirty.data());
                int count = 0;
                for (int t = 0; t < tiles; ++t)
                    if (dirty[size_t(t)]) {
                        ++count;
                        last[size_t(t)] = frame;
                    }
                ok &= count == int(std::min<uint32_t>(param, uint32_t(tiles)));
                if (frame >= bound - 1)
                    for (int t = 0; t < tiles; ++t) ok &= frame - last[size_t(t)] < bound;
            }
            CHECK(ok);
        }
    }
}

static void test_header() {
    std::printf("stream header\n");
    StreamHeader h;
    h.width = 1918;
    h.height = 1080;
    h.tile_w = 32;
    h.tile_h = 48;
    h.refresh_mode = RefreshMode::FullPeriodic;
    h.refresh_param = 600;
    h.quality = 42;
    h.fps_num = 30000;
    h.fps_den = 1001;
    h.audio_codec = 1;
    h.audio_channels = 2;
    h.audio_rate = 48000;
    h.version = 2;
    uint8_t buf[kStreamHeaderSize];
    write_stream_header(h, buf);
    CHECK(std::memcmp(buf, "TJC2", 4) == 0);
    {
        StreamHeader h3 = h;
        h3.version = 3;
        h3.flags = kFlagMotion | kFlagAdaptiveHuffman;
        uint8_t b3[kStreamHeaderSize];
        write_stream_header(h3, b3);
        StreamHeader r3;
        CHECK(std::memcmp(b3, "TJC3", 4) == 0);
        CHECK(parse_stream_header(b3, kStreamHeaderSize, &r3) == Status::Ok);
        CHECK(r3.version == 3 && r3.flags == (kFlagMotion | kFlagAdaptiveHuffman) && r3.fps_num == 30000);
        b3[32] = 0x80;  // unknown flag
        CHECK(parse_stream_header(b3, kStreamHeaderSize, &r3) == Status::BadHeader);
    }
#if !TJC_STREAM_BIG_ENDIAN
    CHECK(buf[4] == (1918 & 255) && buf[5] == (1918 >> 8));
#endif
    StreamHeader r;
    CHECK(parse_stream_header(buf, kStreamHeaderSize, &r) == Status::Ok);
    CHECK(r.width == h.width && r.height == h.height && r.tile_w == 32 && r.tile_h == 48 &&
          r.refresh_mode == h.refresh_mode && r.refresh_param == 600 && r.quality == 42);
    CHECK(r.version == 2 && r.fps_num == 30000 && r.fps_den == 1001 && r.audio_codec == 1 &&
          r.audio_channels == 2 && r.audio_rate == 48000);
    CHECK(parse_stream_header(buf, kStreamHeaderSizeV1, &r) == Status::Truncated);
    buf[8] = 0;  // tile width 0
    CHECK(parse_stream_header(buf, kStreamHeaderSize, &r) == Status::BadHeader);
    buf[8] = 1;
    CHECK(parse_stream_header(buf, kStreamHeaderSize, &r) == Status::Ok);
    buf[27] = 9;  // too many audio channels
    CHECK(parse_stream_header(buf, kStreamHeaderSize, &r) == Status::BadHeader);
    buf[27] = 2;
    buf[18] = buf[19] = buf[20] = buf[21] = 0;  // fps 0
    CHECK(parse_stream_header(buf, kStreamHeaderSize, &r) == Status::BadHeader);
    // A TJC1 header is the first 18 bytes with the old magic: 30 fps, no audio.
    buf[3] = '1';
    CHECK(parse_stream_header(buf, kStreamHeaderSizeV1, &r) == Status::Ok);
    CHECK(r.version == 1 && r.fps_num == 30 && r.fps_den == 1 && r.audio_codec == 0);
    CHECK(validate_geometry(64, 64, 128, 128) != nullptr);  // payload could overflow tile_len
    CHECK(validate_geometry(64, 64, 112, 112) == nullptr);
}

// Synthetic moving test pattern.
static void make_frame(std::vector<uint8_t>& f, int w, int h, int t) {
    int cw = (w + 1) / 2, ch = (h + 1) / 2;
    f.resize(frame_size(w, h));
    uint8_t* y = f.data();
    uint8_t* u = y + size_t(w) * h;
    uint8_t* v = u + size_t(cw) * ch;
    for (int j = 0; j < h; ++j)
        for (int i = 0; i < w; ++i) {
            int val = (i * 3 + j * 2) & 255;  // static gradient
            int bx = 10 + t * 3, by = 20 + t;  // moving box
            if (i >= bx && i < bx + 24 && j >= by && j < by + 24) val = 230 - ((i ^ j) & 31);
            y[size_t(j) * w + i] = uint8_t(val);
        }
    for (int j = 0; j < ch; ++j)
        for (int i = 0; i < cw; ++i) {
            u[size_t(j) * cw + i] = uint8_t(128 + ((i - cw / 2) / 2));
            v[size_t(j) * cw + i] = uint8_t(128 + ((j - ch / 2) / 2));
        }
}

static double psnr(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b) {
    double sse = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        double d = double(a[i]) - double(b[i]);
        sse += d * d;
    }
    return sse == 0 ? 99.0 : 10 * std::log10(255.0 * 255.0 * double(a.size()) / sse);
}

static void run_pipeline(int w, int h, int tile_w, int tile_h, RefreshMode mode, uint32_t param, int quality,
                         DiffReference ref) {
    Config cfg;
    cfg.width = uint16_t(w);
    cfg.height = uint16_t(h);
    cfg.tile_w = uint8_t(tile_w);
    cfg.tile_h = uint8_t(tile_h);
    cfg.refresh_mode = mode;
    cfg.refresh_param = param;
    cfg.quality = uint8_t(quality);
    cfg.diff_reference = ref;
    Encoder enc;
    CHECK(enc.init(cfg));

    std::vector<uint8_t> stream;
    enc.write_stream_header(stream);
    const int frames = 24;
    std::vector<std::vector<uint8_t>> src(frames), recon(frames);
    uint64_t dirty = 0, total = 0;
    for (int t = 0; t < frames; ++t) {
        make_frame(src[size_t(t)], w, h, t);
        FrameStats st;
        enc.encode_frame(src[size_t(t)].data(), stream, &st);
        dirty += st.dirty_tiles;
        total += st.total_tiles;
        recon[size_t(t)].resize(src[size_t(t)].size());
        enc.copy_recon(recon[size_t(t)].data());
    }

    MemoryReader mr;
    mr.data = stream.data();
    mr.size = stream.size();
    Decoder dec;
    CHECK(dec.read_header(MemoryReader::read, &mr) == Status::Ok);
    bool exact = true;
    double min_psnr = 99;
    std::vector<uint8_t> out(frame_size(w, h));
    int decoded = 0;
    for (;;) {
        Status st = dec.decode_frame(MemoryReader::read, &mr);
        if (st == Status::EndOfStream) break;
        CHECK(st == Status::Ok);
        if (st != Status::Ok) break;
        dec.copy_frame(out.data());
        exact &= out == recon[size_t(decoded)];
        min_psnr = std::min(min_psnr, psnr(out, src[size_t(decoded)]));
        ++decoded;
    }
    CHECK(decoded == frames);
    CHECK(exact);  // decoder output == encoder shadow reconstruction, bit for bit
    double raw = double(frame_size(w, h)) * frames;
    std::printf("  %4dx%-4d tile %3dx%-3d q%-3d %-7s ref=%-5s  %7zu B (%5.1f:1)  dirty %5.1f%%  min PSNR %.2f dB\n", w,
                h, tile_w, tile_h, quality,
                mode == RefreshMode::None ? "none" : mode == RefreshMode::Rolling ? "rolling" : "full",
                ref == DiffReference::LastCoded ? "last" : "recon", stream.size(), raw / double(stream.size()),
                100.0 * double(dirty) / double(total), min_psnr);
    if (quality >= 75) CHECK(min_psnr > 30.0);
    CHECK(dirty < total);  // skip tiles are actually being skipped
}

static void test_pipeline() {
    std::printf("full pipeline\n");
    run_pipeline(320, 240, 16, 16, RefreshMode::Rolling, 30, 75, DiffReference::LastCoded);
    run_pipeline(320, 240, 32, 32, RefreshMode::FullPeriodic, 10, 90, DiffReference::Reconstructed);
    run_pipeline(321, 179, 48, 16, RefreshMode::None, 0, 50, DiffReference::LastCoded);  // odd size, padding
    run_pipeline(100, 60, 64, 64, RefreshMode::Rolling, 1, 100, DiffReference::LastCoded);
    run_pipeline(64, 64, 16, 16, RefreshMode::Rolling, 5, 5, DiffReference::LastCoded);
    // Tiles that are not multiples of 16: partial blocks and shared chroma ownership.
    run_pipeline(320, 240, 8, 8, RefreshMode::Rolling, 30, 75, DiffReference::LastCoded);
    run_pipeline(97, 61, 1, 1, RefreshMode::Rolling, 200, 75, DiffReference::LastCoded);
    run_pipeline(97, 61, 3, 5, RefreshMode::FullPeriodic, 7, 75, DiffReference::Reconstructed);
    run_pipeline(131, 77, 13, 9, RefreshMode::None, 0, 90, DiffReference::LastCoded);
    run_pipeline(64, 48, 2, 1, RefreshMode::Rolling, 50, 60, DiffReference::LastCoded);
}

static void test_threads_deterministic() {
    std::printf("thread count does not change the stream\n");
    for (int tile : {1, 7, 16, 64}) {
        std::vector<uint8_t> ref;
        for (int threads : {1, 2, 3, 8}) {
            Config cfg;
            cfg.width = 200;
            cfg.height = 120;
            cfg.tile_w = cfg.tile_h = uint8_t(tile);
            cfg.threads = threads;
            Encoder enc;
            CHECK(enc.init(cfg));
            std::vector<uint8_t> stream, f;
            enc.write_stream_header(stream);
            for (int t = 0; t < 6; ++t) {
                make_frame(f, 200, 120, t);
                enc.encode_frame(f.data(), stream);
            }
            if (ref.empty()) ref = stream;
            CHECK(stream == ref);
        }
    }
}

static void make_audio(std::vector<int16_t>& a, size_t start, size_t n, int ch, uint32_t rate) {
    a.resize(n * size_t(ch));
    for (size_t i = 0; i < n; ++i)
        for (int c = 0; c < ch; ++c) {
            double t = double(start + i) / rate;
            double v = 0.4 * std::sin(2 * kPi * (330 + 110 * c) * t) + 0.2 * std::sin(2 * kPi * 2093 * t) *
                       (std::fmod(t, 0.25) < 0.05 ? 1.0 : 0.0);
            a[i * size_t(ch) + size_t(c)] = int16_t(std::lround(v * 32767));
        }
}

static double snr16(const int16_t* a, const int16_t* b, size_t n) {
    double sig = 0, err = 0;
    for (size_t i = 0; i < n; ++i) {
        sig += double(a[i]) * a[i];
        err += double(a[i] - b[i]) * (a[i] - b[i]);
    }
    return err == 0 ? 99.0 : 10 * std::log10(sig / err);
}

static void test_qoa() {
    std::printf("qoa audio frames\n");
    for (int ch : {1, 2, 8}) {
        for (int len : {1, 19, 20, 21, 1000, detail::qoa::kMaxFrameLen}) {
            std::vector<int16_t> src;
            make_audio(src, 0, size_t(len) * 3, ch, 44100);
            detail::qoa::Lms lms[8];
            for (auto& l : lms) detail::qoa::init_lms(l);
            std::vector<uint8_t> enc;
            for (int f = 0; f < 3; ++f)
                detail::qoa::encode_frame(src.data() + size_t(f) * size_t(len) * size_t(ch), ch, 44100, len, lms, enc);
            CHECK(enc.size() == 3 * detail::qoa::frame_size(ch, len));
            std::vector<int16_t> out(size_t(detail::qoa::kMaxFrameLen) * size_t(ch) * 3);
            size_t pos = 0, got = 0;
            bool ok = true;
            for (int f = 0; f < 3; ++f) {
                size_t used = 0;
                int n = detail::qoa::decode_frame(enc.data() + pos, enc.size() - pos, ch, 44100,
                                                  out.data() + got * size_t(ch), &used);
                ok &= n == len;
                if (n < 0) break;
                pos += used;
                got += size_t(n);
            }
            CHECK(ok && pos == enc.size());
            // Very short frames are dominated by the predictor warming up.
            if (len >= 1000) CHECK(snr16(src.data(), out.data(), src.size()) > 35.0);
            // Wrong channel count or rate must be refused.
            size_t used = 0;
            CHECK(detail::qoa::decode_frame(enc.data(), enc.size(), ch == 1 ? 2 : 1, 44100, out.data(), &used) < 0);
            CHECK(detail::qoa::decode_frame(enc.data(), enc.size(), ch, 48000, out.data(), &used) < 0);
        }
    }
}

// Video + audio through encoder and decoder: exact per-frame sample counts,
// audio quality, and video still bit-exact.
static void test_av(uint32_t fps_num, uint32_t fps_den, uint32_t rate, int ch, int frames, bool short_audio) {
    Config cfg;
    cfg.width = 96;
    cfg.height = 64;
    cfg.fps_num = fps_num;
    cfg.fps_den = fps_den;
    cfg.audio_channels = uint8_t(ch);
    cfg.audio_rate = rate;
    Encoder enc;
    CHECK(enc.init(cfg));
    const size_t total = size_t(uint64_t(frames) * rate * fps_den / fps_num);
    const size_t provided = short_audio ? total / 2 : total + 777;
    std::vector<int16_t> src;
    make_audio(src, 0, provided, ch, rate);

    std::vector<uint8_t> stream, f;
    enc.write_stream_header(stream);
    std::vector<std::vector<uint8_t>> recon(static_cast<size_t>(frames));
    size_t pushed = 0;
    uint64_t padded = 0;
    for (int t = 0; t < frames; ++t) {
        // Push in uneven chunks, sometimes ahead of what the frame needs.
        size_t want = std::min(provided - pushed, enc.samples_for_next_frame() + (t % 3 == 0 ? 37 : 0));
        enc.push_audio(src.data() + pushed * size_t(ch), want);
        pushed += want;
        make_frame(f, 96, 64, t);
        FrameStats st;
        enc.encode_frame(f.data(), stream, &st);
        padded += st.audio_padded;
        recon[size_t(t)].resize(f.size());
        enc.copy_recon(recon[size_t(t)].data());
    }

    MemoryReader mr;
    mr.data = stream.data();
    mr.size = stream.size();
    Decoder dec;
    CHECK(dec.read_header(MemoryReader::read, &mr) == Status::Ok);
    CHECK(dec.header().audio_codec == 1 && dec.header().audio_rate == rate);
    std::vector<int16_t> audio;
    std::vector<uint8_t> out(frame_size(96, 64));
    bool counts_ok = true, video_ok = true;
    int n = 0;
    FrameStats st;
    while (dec.decode_frame(MemoryReader::read, &mr, &st) == Status::Ok) {
        size_t expect = size_t(uint64_t(n + 1) * rate * fps_den / fps_num - uint64_t(n) * rate * fps_den / fps_num);
        counts_ok &= dec.audio_samples() == expect && st.audio_samples == expect;
        audio.insert(audio.end(), dec.audio(), dec.audio() + dec.audio_samples() * size_t(ch));
        dec.copy_frame(out.data());
        video_ok &= out == recon[size_t(n)];
        ++n;
    }
    CHECK(n == frames);
    CHECK(counts_ok);
    CHECK(video_ok);
    CHECK(audio.size() == total * size_t(ch));
    size_t real = std::min(total, provided) * size_t(ch);
    double q = snr16(src.data(), audio.data(), real);
    CHECK(q > 35.0);
    // Padding is digital silence; QOA is lossy, so after the predictor settles it
    // decodes to within a few LSB of zero.
    bool silence = true;
    for (size_t i = real + 200 * size_t(ch); i < audio.size(); ++i) silence &= std::abs(audio[i]) <= 8;
    CHECK(silence);
    CHECK(padded == (short_audio ? total - provided : 0));
    std::printf("  %u/%u fps, %u Hz x%d, %d frames: %zu samples, SNR %.1f dB%s\n", fps_num, fps_den, rate, ch, frames,
                total, q, short_audio ? ", short audio padded with silence" : "");
}

static void test_audio_video() {
    std::printf("audio + video\n");
    test_av(30000, 1001, 48000, 2, 90, false);
    test_av(25, 1, 44100, 1, 50, false);
    test_av(24000, 1001, 22050, 6, 40, true);
    test_av(1, 2, 8000, 1, 3, false);  // 0.5 fps: 16000 samples per frame, several QOA frames each
}

static void test_sync_drift() {
    std::printf("no audio drift over a long run\n");
    Config cfg;
    cfg.width = cfg.height = 16;
    cfg.fps_num = 30000;
    cfg.fps_den = 1001;
    cfg.audio_channels = 1;
    cfg.audio_rate = 1000;
    cfg.threads = 1;
    Encoder enc;
    CHECK(enc.init(cfg));
    std::vector<uint8_t> frame(frame_size(16, 16), 128), out;
    std::vector<int16_t> zeros(64, 0);
    uint64_t sum = 0;
    bool ok = true;
    const int frames = 20000;
    for (int t = 0; t < frames; ++t) {
        enc.push_audio(zeros.data(), enc.samples_for_next_frame());
        FrameStats st;
        out.clear();
        enc.encode_frame(frame.data(), out, &st);
        sum += st.audio_samples;
        ok &= sum == uint64_t(t + 1) * 1000 * 1001 / 30000;
    }
    CHECK(ok);
}

static void test_tjc1_compat() {
    std::printf("TJC1 streams still decode\n");
    Config cfg;
    cfg.width = 80;
    cfg.height = 48;
    cfg.motion = false;
    cfg.adaptive_huffman = false;
    Encoder enc;
    CHECK(enc.init(cfg));
    std::vector<uint8_t> v2, f;
    enc.write_stream_header(v2);
    for (int t = 0; t < 5; ++t) {
        make_frame(f, 80, 48, t);
        enc.encode_frame(f.data(), v2);
    }
    // Without audio, a TJC1 stream is the TJC2 one minus the 16 extra header bytes.
    std::vector<uint8_t> v1(v2.begin(), v2.begin() + kStreamHeaderSizeV1);
    v1[3] = '1';
    v1.insert(v1.end(), v2.begin() + kStreamHeaderSize, v2.end());
    std::vector<uint8_t> a(frame_size(80, 48)), b(a.size());
    MemoryReader r1{v1.data(), v1.size(), 0}, r2{v2.data(), v2.size(), 0};
    Decoder d1, d2;
    CHECK(d1.read_header(MemoryReader::read, &r1) == Status::Ok);
    CHECK(d2.read_header(MemoryReader::read, &r2) == Status::Ok);
    CHECK(d1.header().version == 1 && d1.header().audio_codec == 0);
    int frames = 0;
    bool same = true;
    while (d1.decode_frame(MemoryReader::read, &r1) == Status::Ok) {
        CHECK(d2.decode_frame(MemoryReader::read, &r2) == Status::Ok);
        d1.copy_frame(a.data());
        d2.copy_frame(b.data());
        same &= a == b && d1.audio_samples() == 0;
        ++frames;
    }
    CHECK(frames == 5 && same);
}

// Exact seeking as TJC Studio does it: a parse-only pass (apply = false) records
// per frame the oldest "last sent" frame over all tiles; decoding from there
// must reproduce frame T exactly.
static void test_seek_sync(RefreshMode mode, uint32_t param, uint32_t keyframe_every, int tile, bool motion = true) {
    Config cfg;
    cfg.width = 96;
    cfg.height = 64;
    cfg.motion = motion;
    cfg.tile_w = cfg.tile_h = uint8_t(tile);
    cfg.refresh_mode = mode;
    cfg.refresh_param = param;
    cfg.audio_channels = 1;
    cfg.audio_rate = 8000;
    Encoder enc;
    CHECK(enc.init(cfg));
    std::vector<uint8_t> stream, f;
    std::vector<int16_t> pcm;
    enc.write_stream_header(stream);
    const int frames = 60;
    for (int t = 0; t < frames; ++t) {
        if (keyframe_every && t && t % int(keyframe_every) == 0) enc.force_keyframe();
        make_frame(f, 96, 64, t);
        make_audio(pcm, 0, enc.samples_for_next_frame(), 1, 8000);
        enc.push_audio(pcm.data(), pcm.size());
        enc.encode_frame(f.data(), stream);
    }

    // Reference: every frame decoded from the start.
    std::vector<std::vector<uint8_t>> ref;
    std::vector<FrameStats> ref_stats;
    {
        MemoryReader mr{stream.data(), stream.size(), 0};
        Decoder d;
        CHECK(d.read_header(MemoryReader::read, &mr) == Status::Ok);
        FrameStats st;
        while (d.decode_frame(MemoryReader::read, &mr, &st) == Status::Ok) {
            ref.emplace_back(frame_size(96, 64));
            d.copy_frame(ref.back().data());
            ref_stats.push_back(st);
        }
    }
    CHECK(int(ref.size()) == frames);

    // Index pass without decoding.
    std::vector<size_t> offsets;
    std::vector<uint32_t> sync;
    std::vector<int> tables_at;  // last frame <= t that carried Huffman tables
    {
        MemoryReader mr{stream.data(), stream.size(), 0};
        Decoder d;
        CHECK(d.read_header(MemoryReader::read, &mr) == Status::Ok);
        const size_t tiles = size_t(d.layout().tiles);
        (void)tiles;
        SyncTracker tracker;
        tracker.reset(d.layout());
        int last_tables = -1;
        std::vector<uint8_t> before(frame_size(96, 64)), after(before.size());
        d.copy_frame(before.data());
        FrameStats st;
        bool stats_match = true;
        for (uint32_t t = 0;; ++t) {
            size_t off = mr.pos;
            if (d.decode_frame(MemoryReader::read, &mr, &st, false) != Status::Ok) break;
            stats_match &= st.bytes == ref_stats[t].bytes && st.dirty_tiles == ref_stats[t].dirty_tiles &&
                           st.audio_bytes == ref_stats[t].audio_bytes && d.audio_samples() == 0;
            if (st.tables) last_tables = int(t);
            offsets.push_back(off);
            sync.push_back(tracker.add_frame(t, d.tile_info()));
            tables_at.push_back(last_tables);
        }
        d.copy_frame(after.data());
        CHECK(stats_match);
        CHECK(before == after);  // apply = false leaves the picture alone
    }

    bool exact = true;
    uint64_t decoded = 0;
    for (int t = 0; t < frames; ++t) {
        MemoryReader mr{stream.data(), stream.size(), 0};
        Decoder d;
        d.read_header(MemoryReader::read, &mr);
        uint32_t s0 = sync[size_t(t)];
        int tf = tables_at[size_t(s0)];
        if (tf >= 0 && uint32_t(tf) < s0) {  // load the tables in force at the start frame
            mr.pos = offsets[size_t(tf)];
            d.decode_frame(MemoryReader::read, &mr, nullptr, false);
        }
        mr.pos = offsets[size_t(s0)];
        for (int k = int(sync[size_t(t)]); k <= t; ++k, ++decoded) d.decode_frame(MemoryReader::read, &mr);
        std::vector<uint8_t> out(frame_size(96, 64));
        d.copy_frame(out.data());
        exact &= out == ref[size_t(t)];
    }
    CHECK(exact);
    std::printf("  %-8s param %-3u keyframe %-3u tile %-2d motion %d: avg %.1f frames decoded per seek\n",
                mode == RefreshMode::None ? "none" : mode == RefreshMode::Rolling ? "rolling" : "periodic", param,
                keyframe_every, tile, int(motion), double(decoded) / frames);
}

static void test_seeking() {
    std::printf("exact seeking from sync points\n");
    test_seek_sync(RefreshMode::Rolling, 4, 0, 16);
    test_seek_sync(RefreshMode::Rolling, 1, 0, 8);
    test_seek_sync(RefreshMode::FullPeriodic, 10, 0, 16);
    test_seek_sync(RefreshMode::None, 0, 0, 16);
    test_seek_sync(RefreshMode::None, 0, 15, 7);
    test_seek_sync(RefreshMode::Rolling, 4, 0, 16, false);
    test_seek_sync(RefreshMode::Rolling, 6, 0, 13);
}

// --- TJC3 -------------------------------------------------------------------

static void test_optimal_huffman() {
    std::printf("optimal Huffman tables\n");
    bool ok = true;
    for (int iter = 0; iter < 400; ++iter) {
        uint32_t f[256] = {};
        int kind = iter % 4;
        for (int i = 0; i < 256; ++i) {
            if (kind == 0) f[i] = rnd() % 1000;                        // random
            else if (kind == 1) f[i] = (rnd() % 8 == 0) ? rnd() % 50 : 0; // sparse
            else if (kind == 2) f[i] = i < 40 ? (1u << (i % 30)) : 0;   // forces the 16-bit limit
            else f[i] = i == int(rnd() % 256) ? 1000 : 0;               // a single symbol
        }
        detail::HuffSet hs;
        for (int t = 0; t < 4; ++t) hs.spec[t] = detail::optimal_huff(f);
        bool any = false;
        for (uint32_t v : f) any |= v != 0;
        if (!any) continue;
        ok &= hs.build();
        int maxlen = 0;
        for (int sym = 0; sym < 256; ++sym) {
            ok &= (f[sym] != 0) == (hs.enc[0].size[sym] != 0);
            maxlen = std::max(maxlen, int(hs.enc[0].size[sym]));
        }
        ok &= maxlen <= 16;
        // Every code decodes back to its symbol.
        for (int sym = 0; sym < 256 && ok; ++sym) {
            if (!f[sym]) continue;
            std::vector<uint8_t> buf;
            detail::BitWriter bw(buf);
            bw.put(hs.enc[0].code[sym], hs.enc[0].size[sym]);
            bw.put(0xffff, 16);
            bw.flush();
            detail::BitReader br(buf.data(), buf.size());
            ok &= detail::huff_decode(br, hs.dec[0]) == sym;
        }
    }
    CHECK(ok);
    // Malformed tables are refused.
    detail::HuffSet bad = detail::standard_huff();
    bad.spec[1].bits[0] = 3;  // three 1-bit codes
    CHECK(!bad.build());
    bad = detail::standard_huff();
    bad.spec[0].vals[1] = bad.spec[0].vals[0];  // duplicate symbol
    CHECK(!bad.build());
}

static void test_exp_golomb() {
    std::printf("exp-golomb\n");
    std::vector<uint8_t> buf;
    detail::BitWriter bw(buf);
    std::vector<int> vals;
    for (int v = -3000; v <= 3000; v += 7) vals.push_back(v);
    vals.push_back(0);
    vals.push_back(65535);
    vals.push_back(-65536);
    size_t bits = 0;
    for (int v : vals) {
        detail::put_se(bw, v);
        bits += size_t(detail::se_bits(v));
    }
    bw.flush();
    CHECK(buf.size() == (bits + 7) / 8);
    detail::BitReader br(buf.data(), buf.size());
    bool ok = true;
    for (int v : vals) {
        int got = 0;
        ok &= detail::get_se(br, &got) && got == v;
    }
    CHECK(ok);
    std::vector<uint8_t> zeros(16, 0);
    detail::BitReader bz(zeros.data(), zeros.size());
    int dummy;
    CHECK(!detail::get_se(bz, &dummy));  // runaway prefix is rejected
}

static void test_predict_block() {
    std::printf("half-pel prediction\n");
    std::vector<uint8_t> plane(40 * 40);
    for (auto& p : plane) p = uint8_t(rnd());
    bool ok = true;
    for (int mvy = -6; mvy <= 6; ++mvy)
        for (int mvx = -6; mvx <= 6; ++mvx) {
            uint8_t out[64];
            detail::predict_block(plane.data(), 40, 16, 16, mvx, mvy, 8, 8, out, 8);
            for (int y = 0; y < 8; ++y)
                for (int x = 0; x < 8; ++x) {
                    int sx = 16 + x + (mvx >> 1), sy = 16 + y + (mvy >> 1), fx = mvx & 1, fy = mvy & 1;
                    auto P = [&](int a, int b) { return int(plane[size_t(b) * 40 + size_t(a)]); };
                    int e = !fx && !fy ? P(sx, sy)
                            : fx && !fy ? (P(sx, sy) + P(sx + 1, sy) + 1) >> 1
                            : !fx       ? (P(sx, sy) + P(sx, sy + 1) + 1) >> 1
                                        : (P(sx, sy) + P(sx + 1, sy) + P(sx, sy + 1) + P(sx + 1, sy + 1) + 2) >> 2;
                    ok &= out[y * 8 + x] == e;
                }
        }
    CHECK(ok);
}

// Smooth texture panning by `speed` half-pels per frame (sub-pixel motion).
static void make_pan(std::vector<uint8_t>& f, int w, int h, int t, int speed, bool noise) {
    int cw = (w + 1) / 2, ch = (h + 1) / 2;
    f.resize(frame_size(w, h));
    double ox = t * speed * 0.5;
    for (int j = 0; j < h; ++j)
        for (int i = 0; i < w; ++i) {
            double x = i + ox;
            double v = 128 + 55 * std::sin(x * 0.21) * std::cos(j * 0.13) + 40 * std::sin((x * 0.05 + j * 0.07) * 3);
            if (noise) v += double(int(rnd() % 9) - 4);
            f[size_t(j) * w + i] = uint8_t(std::max(0.0, std::min(255.0, v)));
        }
    uint8_t* u = f.data() + size_t(w) * h;
    uint8_t* v = u + size_t(cw) * ch;
    for (int j = 0; j < ch; ++j)
        for (int i = 0; i < cw; ++i) {
            double x = 2 * i + ox;
            u[size_t(j) * cw + i] = uint8_t(128 + 40 * std::sin(x * 0.03));
            v[size_t(j) * cw + i] = uint8_t(128 + 40 * std::cos(x * 0.04 + j * 0.05));
        }
}

struct RunResult {
    std::vector<uint8_t> stream;
    bool exact = true;
    double min_psnr = 99;
    uint64_t inter = 0, dropped = 0, tables = 0;
};

static RunResult run_tjc3(Config cfg, int frames, int speed, bool noise, bool box = false) {
    RunResult r;
    Encoder enc;
    CHECK(enc.init(cfg));
    enc.write_stream_header(r.stream);
    std::vector<std::vector<uint8_t>> src, rec;
    g_rng = 777;  // same noise for every configuration
    for (int t = 0; t < frames; ++t) {
        std::vector<uint8_t> f;
        if (box) make_frame(f, cfg.width, cfg.height, t);
        else make_pan(f, cfg.width, cfg.height, t, speed, noise);
        FrameStats st;
        enc.encode_frame(f.data(), r.stream, &st);
        r.inter += st.inter_tiles;
        r.dropped += st.dropped_tiles;
        r.tables += st.tables;
        rec.emplace_back(f.size());
        enc.copy_recon(rec.back().data());
        src.push_back(std::move(f));
    }
    MemoryReader mr{r.stream.data(), r.stream.size(), 0};
    Decoder d;
    CHECK(d.read_header(MemoryReader::read, &mr) == Status::Ok);
    std::vector<uint8_t> out(frame_size(cfg.width, cfg.height));
    int n = 0;
    FrameStats st;
    while (d.decode_frame(MemoryReader::read, &mr, &st) == Status::Ok) {
        d.copy_frame(out.data());
        r.exact &= out == rec[size_t(n)];
        r.min_psnr = std::min(r.min_psnr, psnr(out, src[size_t(n)]));
        ++n;
    }
    CHECK(n == frames);
    return r;
}

static void test_tjc3_pipeline() {
    std::printf("TJC3: motion compensation, adaptive Huffman, skip-invisible\n");
    // Small frames have few tiles, so the rolling refresh is kept slow (refresh_param)
    // or it would force most tiles to intra and hide what motion compensation does.
    struct Case { const char* name; int w, h, tile, speed; bool noise, motion, huff; double skip; uint32_t refresh; };
    const Case cases[] = {
        {"pan 1.5px/frame", 160, 96, 16, 3, false, true, true, 0, 2},
        {"pan, no motion comp.", 160, 96, 16, 3, false, false, true, 0, 2},
        {"pan, standard tables", 160, 96, 16, 3, false, true, false, 0, 2},
        {"TJC2 (both off)", 160, 96, 16, 3, false, false, false, 0, 2},
        {"static noise + skip", 160, 96, 16, 0, true, true, true, 2.0, 2},
        {"static noise, no skip", 160, 96, 16, 0, true, true, true, 0, 2},
        {"odd tiles 13x9", 131, 77, 13, 5, false, true, true, 0, 2},
        {"tiles 7x5 + skip", 97, 61, 7, 0, true, true, true, 1.0, 4},
        {"fast pan 9.5px", 160, 96, 32, 19, false, true, true, 0, 1},
        {"default refresh", 160, 96, 16, 3, false, true, true, 0, 30},
    };
    size_t with_mc = 0, without_mc = 0, std_tables = 0, skip_on = 0, skip_off = 0;
    for (const Case& c : cases) {
        Config cfg;
        cfg.width = uint16_t(c.w);
        cfg.height = uint16_t(c.h);
        cfg.tile_w = cfg.tile_h = uint8_t(c.tile);
        cfg.motion = c.motion;
        cfg.adaptive_huffman = c.huff;
        cfg.skip_invisible = c.skip;
        cfg.refresh_param = c.refresh;
        cfg.motion_threshold_k = c.noise ? 1 : 3;
        RunResult r = run_tjc3(cfg, 24, c.speed, c.noise);
        CHECK(r.exact);
        CHECK(r.min_psnr > 30.0);
        bool v2 = !c.motion && !c.huff;
        CHECK(std::memcmp(r.stream.data(), v2 ? "TJC2" : "TJC3", 4) == 0);
        if (!c.motion) CHECK(r.inter == 0);
        if (c.motion && c.speed) CHECK(r.inter > 0);
        if (!c.huff) CHECK(r.tables == 0);
        if (c.skip > 0) CHECK(r.dropped > 0);
        if (std::string(c.name) == "pan 1.5px/frame") with_mc = r.stream.size();
        if (std::string(c.name) == "pan, no motion comp.") without_mc = r.stream.size();
        if (std::string(c.name) == "pan, standard tables") std_tables = r.stream.size();
        if (std::string(c.name) == "static noise + skip") skip_on = r.stream.size();
        if (std::string(c.name) == "static noise, no skip") skip_off = r.stream.size();
        std::printf("  %-22s %7zu B  min PSNR %5.2f  inter tiles %5llu  table updates %3llu  dropped %llu\n", c.name,
                    r.stream.size(), r.min_psnr, (unsigned long long)r.inter, (unsigned long long)r.tables,
                    (unsigned long long)r.dropped);
    }
    CHECK(with_mc * 10 < without_mc * 6);  // motion compensation must pay off on a pan
    CHECK(with_mc <= std_tables);           // adaptive tables never lose
    CHECK(skip_on * 10 < skip_off * 8);     // skipping invisible updates saves on grain

    // Output does not depend on the thread count.
    std::vector<uint8_t> ref;
    for (int threads : {1, 2, 3, 8}) {
        Config cfg;
        cfg.width = 200;
        cfg.height = 120;
        cfg.threads = threads;
        cfg.skip_invisible = 1.0;
        RunResult r = run_tjc3(cfg, 10, 3, true);
        if (ref.empty()) ref = r.stream;
        CHECK(r.stream == ref);
    }
    // The moving-box pattern with every tool on, against the decoder.
    Config cfg;
    cfg.width = 96;
    cfg.height = 64;
    RunResult r = run_tjc3(cfg, 30, 0, false, true);
    CHECK(r.exact);
}

static void test_corrupt_input() {
    std::printf("corrupt / truncated input\n");
    Config cfg;
    cfg.width = 64;
    cfg.height = 48;
    cfg.audio_channels = 2;
    cfg.audio_rate = 44100;
    Encoder enc;
    CHECK(enc.init(cfg));
    std::vector<uint8_t> stream, f;
    std::vector<int16_t> pcm;
    enc.write_stream_header(stream);
    for (int t = 0; t < 4; ++t) {
        make_frame(f, 64, 48, t);
        make_audio(pcm, size_t(t) * 1470, enc.samples_for_next_frame(), 2, 44100);
        enc.push_audio(pcm.data(), pcm.size() / 2);
        enc.encode_frame(f.data(), stream);
    }
    // Truncation at every possible length must never crash and never report Ok past the end.
    for (size_t cut = 0; cut < stream.size(); cut += 7) {
        MemoryReader mr;
        mr.data = stream.data();
        mr.size = cut;
        Decoder dec;
        if (dec.read_header(MemoryReader::read, &mr) != Status::Ok) continue;
        while (dec.decode_frame(MemoryReader::read, &mr) == Status::Ok) {
        }
    }
    // Random byte flips must be handled gracefully as well.
    int errors = 0;
    for (int iter = 0; iter < 300; ++iter) {
        std::vector<uint8_t> bad = stream;
        for (int k = 0; k < 3; ++k) bad[kStreamHeaderSize + rnd() % (bad.size() - kStreamHeaderSize)] ^= uint8_t(1 + rnd() % 255);
        MemoryReader mr;
        mr.data = bad.data();
        mr.size = bad.size();
        Decoder dec;
        CHECK(dec.read_header(MemoryReader::read, &mr) == Status::Ok);
        Status st;
        while ((st = dec.decode_frame(MemoryReader::read, &mr)) == Status::Ok) {
        }
        errors += st != Status::EndOfStream;
    }
    std::printf("  %d/300 damaged streams reported an error (rest decoded to garbage pixels/audio)\n", errors);
    CHECK(errors > 0);
}

int main() {
    test_bitstream();
    test_dct();
    test_huffman_tables();
    test_block_roundtrip();
    test_bitmap();
    test_sad();
    test_refresh_policy();
    test_header();
    test_quant_div();
    test_pipeline();
    test_threads_deterministic();
    test_qoa();
    test_audio_video();
    test_sync_drift();
    test_tjc1_compat();
    test_seeking();
    test_optimal_huffman();
    test_exp_golomb();
    test_predict_block();
    test_tjc3_pipeline();
    test_corrupt_input();
    std::printf("\n%d checks, %d failed\n", g_checks, g_failed);
    return g_failed ? 1 : 0;
}
