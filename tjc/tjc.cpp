// tjc - command line encoder/decoder for the Tiled-JPEG Codec (TJC1).
//
//   tjc encode --width W --height H [options] [-i in.yuv] [-o out.tjc]
//   tjc decode [-i in.tjc] [-o out.yuv]
//   tjc info   [-i in.tjc]
//
// Input/output default to stdin/stdout, so it sits in an ffmpeg pipe:
//   ffmpeg -i in.mp4 -f rawvideo -pix_fmt yuv420p - | tjc encode --width 320 --height 240 > s.tjc
//   tjc decode < s.tjc | ffplay -f rawvideo -pix_fmt yuv420p -video_size 320x240 -
//
// When the binary is named (or symlinked as) tjc_encode / tjc_decode the
// subcommand can be left out.
//
// Build: g++ -std=c++17 -O2 -o tjc tjc.cpp

#define TJC_IMPLEMENTATION
#include "tjc.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#ifdef _WIN32
    #include <io.h>
    #include <fcntl.h>
#endif

namespace {

void usage() {
    std::fprintf(stderr,
        "Tiled-JPEG Codec (TJC1)\n"
        "\n"
        "usage:\n"
        "  tjc encode --width W --height H [options] [-i in.yuv] [-o out.tjc]\n"
        "  tjc decode [-i in.tjc] [-o out.yuv] [-v] [-q]\n"
        "  tjc info   [-i in.tjc] [-v]\n"
        "\n"
        "encode options:\n"
        "  --width, -W N              frame width in pixels (required)\n"
        "  --height, -H N             frame height in pixels (required)\n"
        "  --size WxH                 shorthand for --width/--height\n"
        "  --tile-size WxH | N        tile size 1..255 (default 16x16; below 8x8 costs more bits, see README)\n"
        "  --motion-threshold K       dirty if mean |diff| per sample > K (default 3, 0 = any change)\n"
        "  --refresh-mode MODE        none | full | rolling (default rolling)\n"
        "  --refresh-param N          full: refresh every N frames, rolling: N tiles per frame (default 30)\n"
        "  --quality Q                1..100 (default 75)\n"
        "  --diff-ref REF             last-coded | recon: what new frames are compared to (default last-coded)\n"
        "  --keyframe-every N         also force a full refresh every N frames (default 0 = off)\n"
        "  --psnr                     report PSNR of the reconstruction vs the source\n"
        "  --threads N                encoder threads, 0 = all cores (default 0)\n"
        "\n"
        "common options:\n"
        "  -i FILE                    input (default stdin, '-' = stdin)\n"
        "  -o FILE                    output (default stdout, '-' = stdout)\n"
        "  --fps N                    frame rate used for the bitrate in the summary (default 30)\n"
        "  -v, --verbose              per-frame log on stderr\n"
        "  -q, --quiet                no summary\n"
        "\n"
        "Input/output frames are raw YUV420P (ffmpeg: -f rawvideo -pix_fmt yuv420p).\n");
}

struct Options {
    std::string mode;
    std::string in = "-", out = "-";
    tjc::Config cfg;
    bool have_w = false, have_h = false;
    uint32_t keyframe_every = 0;
    double fps = 30.0;
    bool verbose = false, quiet = false, psnr = false;
};

[[noreturn]] void die(const char* fmt, const char* arg = "") {
    std::fprintf(stderr, "tjc: ");
    std::fprintf(stderr, fmt, arg);
    std::fprintf(stderr, "\n");
    std::exit(1);
}

long parse_int(const char* s, const char* what, long lo, long hi) {
    char* end = nullptr;
    long v = std::strtol(s, &end, 10);
    if (!*s || *end || v < lo || v > hi) die("invalid value for %s", what);
    return v;
}

void parse_pair(const char* s, const char* what, long& a, long& b) {
    std::string str(s);
    size_t x = str.find_first_of("xX");
    if (x == std::string::npos) {
        a = b = parse_int(s, what, 1, 65535);
        return;
    }
    a = parse_int(str.substr(0, x).c_str(), what, 1, 65535);
    b = parse_int(str.substr(x + 1).c_str(), what, 1, 65535);
}

Options parse_args(int argc, char** argv) {
    Options o;
    int i = 1;
    std::string prog = argv[0];
    size_t slash = prog.find_last_of("/\\");
    if (slash != std::string::npos) prog = prog.substr(slash + 1);
    if (prog.rfind("tjc_encode", 0) == 0) o.mode = "encode";
    else if (prog.rfind("tjc_decode", 0) == 0) o.mode = "decode";
    else if (prog.rfind("tjc_info", 0) == 0) o.mode = "info";

    if (o.mode.empty()) {
        if (argc < 2) {
            usage();
            std::exit(1);
        }
        o.mode = argv[1];
        i = 2;
        if (o.mode == "-h" || o.mode == "--help" || o.mode == "help") {
            usage();
            std::exit(0);
        }
        if (o.mode != "encode" && o.mode != "decode" && o.mode != "info") die("unknown command '%s'", argv[1]);
    }

    for (; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> const char* {
            if (i + 1 >= argc) die("missing value for %s", argv[i]);
            return argv[++i];
        };
        if (a == "-h" || a == "--help") {
            usage();
            std::exit(0);
        } else if (a == "-i" || a == "--input") {
            o.in = next();
        } else if (a == "-o" || a == "--output") {
            o.out = next();
        } else if (a == "-v" || a == "--verbose") {
            o.verbose = true;
        } else if (a == "-q" || a == "--quiet") {
            o.quiet = true;
        } else if (a == "--fps") {
            o.fps = std::atof(next());
            if (!(o.fps > 0)) die("invalid value for --fps");
        } else if (a == "--width" || a == "-W") {
            o.cfg.width = uint16_t(parse_int(next(), "--width", 1, 65535));
            o.have_w = true;
        } else if (a == "--height" || a == "-H") {
            o.cfg.height = uint16_t(parse_int(next(), "--height", 1, 65535));
            o.have_h = true;
        } else if (a == "--size" || a == "-s" || a == "--video-size") {
            long w, h;
            parse_pair(next(), "--size", w, h);
            o.cfg.width = uint16_t(w);
            o.cfg.height = uint16_t(h);
            o.have_w = o.have_h = true;
        } else if (a == "--tile-size") {
            long w, h;
            parse_pair(next(), "--tile-size", w, h);
            if (w > 255 || h > 255) die("--tile-size must be at most 255x255");
            o.cfg.tile_w = uint8_t(w);
            o.cfg.tile_h = uint8_t(h);
        } else if (a == "--motion-threshold") {
            o.cfg.motion_threshold_k = uint32_t(parse_int(next(), "--motion-threshold", 0, 255));
        } else if (a == "--refresh-mode") {
            std::string m = next();
            if (m == "none") o.cfg.refresh_mode = tjc::RefreshMode::None;
            else if (m == "full" || m == "full-periodic" || m == "periodic") o.cfg.refresh_mode = tjc::RefreshMode::FullPeriodic;
            else if (m == "rolling") o.cfg.refresh_mode = tjc::RefreshMode::Rolling;
            else die("unknown refresh mode '%s' (none | full | rolling)", m.c_str());
        } else if (a == "--refresh-param") {
            o.cfg.refresh_param = uint32_t(parse_int(next(), "--refresh-param", 0, 65535));
        } else if (a == "--quality") {
            o.cfg.quality = uint8_t(parse_int(next(), "--quality", 1, 100));
        } else if (a == "--diff-ref") {
            std::string m = next();
            if (m == "last-coded" || m == "source") o.cfg.diff_reference = tjc::DiffReference::LastCoded;
            else if (m == "recon" || m == "reconstructed") o.cfg.diff_reference = tjc::DiffReference::Reconstructed;
            else die("unknown diff reference '%s' (last-coded | recon)", m.c_str());
        } else if (a == "--keyframe-every") {
            o.keyframe_every = uint32_t(parse_int(next(), "--keyframe-every", 0, 1L << 30));
        } else if (a == "--threads") {
            o.cfg.threads = int(parse_int(next(), "--threads", 0, 256));
        } else if (a == "--psnr") {
            o.psnr = true;
        } else {
            die("unknown option '%s' (see --help)", a.c_str());
        }
    }
    return o;
}

FILE* open_in(const std::string& path) {
    if (path == "-") {
#ifdef _WIN32
        _setmode(_fileno(stdin), _O_BINARY);
#endif
        return stdin;
    }
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) die("cannot open input '%s'", path.c_str());
    return f;
}

FILE* open_out(const std::string& path) {
    if (path == "-") {
#ifdef _WIN32
        _setmode(_fileno(stdout), _O_BINARY);
#endif
        return stdout;
    }
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) die("cannot open output '%s'", path.c_str());
    return f;
}

void write_all(FILE* f, const void* data, size_t n) {
    if (n && std::fwrite(data, 1, n, f) != n) die("write failed");
}

size_t file_read(void* user, void* dst, size_t n) { return std::fread(dst, 1, n, static_cast<FILE*>(user)); }

double psnr(double sse, double count) {
    if (sse <= 0) return 99.0;
    return 10.0 * std::log10(255.0 * 255.0 * count / sse);
}

const char* refresh_name(tjc::RefreshMode m) {
    switch (m) {
        case tjc::RefreshMode::None: return "none";
        case tjc::RefreshMode::FullPeriodic: return "full-periodic";
        case tjc::RefreshMode::Rolling: return "rolling";
    }
    return "?";
}

void print_header(const tjc::StreamHeader& h, const tjc::Layout& l) {
    std::fprintf(stderr, "stream: %ux%u, tiles %ux%u (%dx%d = %d), quality %u, refresh %s/%u\n", h.width,
                 h.height, h.tile_w, h.tile_h, l.cols, l.rows, l.tiles, h.quality, refresh_name(h.refresh_mode),
                 h.refresh_param);
}

void log_frame(const tjc::FrameStats& s) {
    std::fprintf(stderr, "frame %6u %c dirty %5u/%-5u (%5.1f%%) %8zu B\n", s.frame_num, s.force_refresh ? 'I' : 'P',
                 s.dirty_tiles, s.total_tiles, 100.0 * s.dirty_tiles / (s.total_tiles ? s.total_tiles : 1), s.bytes);
}

struct Totals {
    uint64_t frames = 0, iframes = 0, bytes = 0, dirty = 0, tiles = 0;
    void add(const tjc::FrameStats& s) {
        ++frames;
        iframes += s.force_refresh;
        bytes += s.bytes;
        dirty += s.dirty_tiles;
        tiles += s.total_tiles;
    }
    void print(const char* what, size_t raw_frame, double fps, uint64_t stream_bytes) const {
        double raw = double(raw_frame) * double(frames);
        std::fprintf(stderr,
                     "%s: %llu frames (%llu full refresh), %llu bytes, %.1f B/frame avg, %.1f kbit/s @ %.4g fps\n"
                     "      dirty tiles %.1f%%, compression %.1f:1 vs raw yuv420p\n",
                     what, (unsigned long long)frames, (unsigned long long)iframes, (unsigned long long)stream_bytes,
                     frames ? double(bytes) / frames : 0.0, frames ? double(stream_bytes) * 8.0 * fps / frames / 1000.0 : 0.0,
                     fps, tiles ? 100.0 * dirty / tiles : 0.0, stream_bytes ? raw / double(stream_bytes) : 0.0);
    }
};

int run_encode(const Options& o) {
    if (!o.have_w || !o.have_h) die("encode needs --width and --height (or --size WxH)");
    tjc::Encoder enc;
    if (!enc.init(o.cfg)) die("%s", enc.error());

    FILE* in = open_in(o.in);
    FILE* out = open_out(o.out);
    const size_t fsize = tjc::frame_size(o.cfg.width, o.cfg.height);
    std::vector<uint8_t> frame(fsize), recon(o.psnr ? fsize : 0), buf;
    buf.reserve(fsize);

    enc.write_stream_header(buf);
    write_all(out, buf.data(), buf.size());
    uint64_t stream_bytes = buf.size();
    if (o.verbose) print_header(enc.stream_header(), enc.layout());

    Totals tot;
    double sse[3] = {0, 0, 0};
    const tjc::Layout& l = enc.layout();
    const size_t plane_off[3] = {0, size_t(l.width) * l.height, size_t(l.width) * l.height + size_t(l.cwidth) * l.cheight};
    const size_t plane_len[3] = {size_t(l.width) * l.height, size_t(l.cwidth) * l.cheight, size_t(l.cwidth) * l.cheight};

    for (;;) {
        size_t got = std::fread(frame.data(), 1, fsize, in);
        if (got == 0) break;
        if (got < fsize) {
            std::fprintf(stderr, "tjc: warning: dropped trailing partial frame (%zu of %zu bytes)\n", got, fsize);
            break;
        }
        if (o.keyframe_every && tot.frames && tot.frames % o.keyframe_every == 0) enc.force_keyframe();
        tjc::FrameStats st;
        buf.clear();
        enc.encode_frame(frame.data(), buf, &st);
        write_all(out, buf.data(), buf.size());
        stream_bytes += buf.size();
        tot.add(st);
        if (o.verbose) log_frame(st);
        if (o.psnr) {
            enc.copy_recon(recon.data());
            for (int c = 0; c < 3; ++c) {
                const uint8_t* a = frame.data() + plane_off[c];
                const uint8_t* b = recon.data() + plane_off[c];
                double e = 0;
                for (size_t k = 0; k < plane_len[c]; ++k) {
                    double d = double(a[k]) - double(b[k]);
                    e += d * d;
                }
                sse[c] += e;
            }
        }
    }
    std::fflush(out);
    if (!o.quiet) {
        tot.print("encoded", fsize, o.fps, stream_bytes);
        if (o.psnr && tot.frames) {
            double n[3] = {double(plane_len[0]) * tot.frames, double(plane_len[1]) * tot.frames,
                           double(plane_len[2]) * tot.frames};
            std::fprintf(stderr, "      PSNR Y %.2f dB, Cb %.2f dB, Cr %.2f dB, all %.2f dB\n", psnr(sse[0], n[0]),
                         psnr(sse[1], n[1]), psnr(sse[2], n[2]), psnr(sse[0] + sse[1] + sse[2], n[0] + n[1] + n[2]));
        }
    }
    if (in != stdin) std::fclose(in);
    if (out != stdout) std::fclose(out);
    return 0;
}

int run_decode(const Options& o, bool write_frames) {
    FILE* in = open_in(o.in);
    FILE* out = write_frames ? open_out(o.out) : nullptr;
    tjc::Decoder dec;
    tjc::Status st = dec.read_header(file_read, in);
    if (st == tjc::Status::EndOfStream) die("empty input");
    if (st != tjc::Status::Ok) die("%s", tjc::status_string(st));
    const tjc::StreamHeader& h = dec.header();
    if (o.verbose || !write_frames) print_header(h, dec.layout());
    if (write_frames && !o.quiet)
        std::fprintf(stderr, "output: rawvideo yuv420p %ux%u  (ffplay -f rawvideo -pix_fmt yuv420p -video_size %ux%u -)\n",
                     h.width, h.height, h.width, h.height);

    const size_t fsize = tjc::frame_size(h.width, h.height);
    std::vector<uint8_t> frame(write_frames ? fsize : 0);
    Totals tot;
    uint64_t stream_bytes = tjc::kStreamHeaderSize;
    int rc = 0;
    for (;;) {
        tjc::FrameStats fs;
        st = dec.decode_frame(file_read, in, &fs);
        if (st == tjc::Status::EndOfStream) break;
        if (st != tjc::Status::Ok) {
            std::fprintf(stderr, "tjc: frame %llu: %s\n", (unsigned long long)tot.frames, tjc::status_string(st));
            rc = 1;
            break;
        }
        tot.add(fs);
        stream_bytes += fs.bytes;
        if (o.verbose || !write_frames) log_frame(fs);
        if (write_frames) {
            dec.copy_frame(frame.data());
            write_all(out, frame.data(), fsize);
        }
    }
    if (out) std::fflush(out);
    if (!o.quiet) tot.print(write_frames ? "decoded" : "stream", fsize, o.fps, stream_bytes);
    if (in != stdin) std::fclose(in);
    if (out && out != stdout) std::fclose(out);
    return rc;
}

}  // namespace

int main(int argc, char** argv) {
    Options o = parse_args(argc, argv);
    if (o.mode == "encode") return run_encode(o);
    if (o.mode == "decode") return run_decode(o, true);
    return run_decode(o, false);
}
