#include "encode_job.h"

#include "platform.h"
#include "yuv.h"

#include <chrono>
#include <cstdlib>
#include <map>
#include <sstream>

namespace gui {

namespace {

bool parse_rational(const std::string& s, uint32_t& num, uint32_t& den) {
    size_t slash = s.find('/');
    if (slash == std::string::npos) return false;
    unsigned long long n = std::strtoull(s.c_str(), nullptr, 10);
    unsigned long long d = std::strtoull(s.c_str() + slash + 1, nullptr, 10);
    if (n == 0 || d == 0 || n > 0xffffffffull || d > 0xffffffffull) return false;
    num = uint32_t(n);
    den = uint32_t(d);
    return true;
}

std::string last_lines(const std::string& text, size_t lines) {
    size_t end = text.size();
    while (end && (text[end - 1] == '\n' || text[end - 1] == '\r')) --end;
    size_t pos = end;
    for (size_t n = 0; pos && n < lines;) {
        --pos;
        if (text[pos] == '\n') ++n;
    }
    if (pos && text[pos] == '\n') ++pos;
    return text.substr(pos, end - pos);
}

std::string fps_string(uint32_t num, uint32_t den) { return std::to_string(num) + "/" + std::to_string(den); }

}  // namespace

MediaInfo probe_media(const std::string& ffprobe, const std::string& path) {
    MediaInfo info;
    std::string log = temp_path("probe.log");
    Process p;
    std::string err;
    if (!p.start({ffprobe, "-v", "error", "-show_entries",
                  "stream=codec_type,codec_name,width,height,r_frame_rate,avg_frame_rate,sample_rate,channels:"
                  "format=duration",
                  "-of", "compact=p=0", path},
                 log, &err)) {
        info.error = err;
        return info;
    }
    std::string out;
    char buf[4096];
    while (size_t n = p.read(buf, sizeof(buf))) out.append(buf, n);
    int code = p.wait();
    std::string errtext = read_text_file(log);
    remove_utf8(log);
    if (code != 0) {
        info.error = errtext.empty() ? "ffprobe could not read the file" : last_lines(errtext, 3);
        return info;
    }

    std::istringstream lines(out);
    std::string line;
    bool have_video = false;
    while (std::getline(lines, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::map<std::string, std::string> kv;
        std::istringstream fields(line);
        std::string field;
        while (std::getline(fields, field, '|')) {
            size_t eq = field.find('=');
            if (eq != std::string::npos) kv[field.substr(0, eq)] = field.substr(eq + 1);
        }
        if (kv.count("duration") && !kv.count("codec_type")) {
            info.duration = std::atof(kv["duration"].c_str());
        } else if (kv["codec_type"] == "video" && !have_video) {
            have_video = true;
            info.width = std::atoi(kv["width"].c_str());
            info.height = std::atoi(kv["height"].c_str());
            info.video_codec = kv["codec_name"];
            if (!parse_rational(kv["avg_frame_rate"], info.fps_num, info.fps_den) &&
                !parse_rational(kv["r_frame_rate"], info.fps_num, info.fps_den)) {
                info.fps_num = 30;
                info.fps_den = 1;
            }
        } else if (kv["codec_type"] == "audio" && !info.has_audio) {
            info.has_audio = true;
            info.audio_rate = std::atoi(kv["sample_rate"].c_str());
            info.audio_channels = std::atoi(kv["channels"].c_str());
            info.audio_codec = kv["codec_name"];
        }
    }
    if (!have_video || info.width <= 0 || info.height <= 0) {
        info.error = "no video stream found";
        return info;
    }
    if (info.audio_rate <= 0 || info.audio_channels <= 0) info.has_audio = false;
    info.ok = true;
    return info;
}

EncodeJob::~EncodeJob() {
    cancel();
    if (thread_.joinable()) thread_.join();
}

bool EncodeJob::start(const EncodeSettings& s, std::string* err) {
    if (busy()) {
        if (err) *err = "an encode is already running";
        return false;
    }
    if (thread_.joinable()) thread_.join();
    tjc::Encoder probe;
    if (!probe.init(s.cfg)) {
        if (err) *err = probe.error();
        return false;
    }
    s_ = s;
    cancel_ = false;
    finished_ = false;
    {
        std::lock_guard<std::mutex> lk(mu_);
        p_ = EncodeProgress();
        p_.running = true;
        p_.message = "starting ffmpeg...";
        double fps = double(s.cfg.fps_num) / s.cfg.fps_den;
        p_.total_estimate = s.duration > 0 ? uint64_t(s.duration * fps + 0.5) : 0;
        preview_new_ = false;
    }
    thread_ = std::thread(&EncodeJob::run, this);
    return true;
}

void EncodeJob::cancel() { cancel_ = true; }

EncodeProgress EncodeJob::progress() const {
    std::lock_guard<std::mutex> lk(mu_);
    return p_;
}

bool EncodeJob::take_preview(std::vector<uint8_t>& rgba, int& w, int& h) {
    std::lock_guard<std::mutex> lk(mu_);
    if (!preview_new_) return false;
    rgba.swap(preview_);
    w = s_.cfg.width;
    h = s_.cfg.height;
    preview_new_ = false;
    return true;
}

void EncodeJob::set_message(const std::string& m) {
    std::lock_guard<std::mutex> lk(mu_);
    p_.message = m;
}

void EncodeJob::run() {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    const tjc::Config& cfg = s_.cfg;
    const std::string vlog = temp_path("video.log"), alog = temp_path("audio.log");
    auto finish = [&](bool failed, bool cancelled, const std::string& msg) {
        std::lock_guard<std::mutex> lk(mu_);
        p_.running = false;
        p_.done = !failed && !cancelled;
        p_.failed = failed;
        p_.cancelled = cancelled;
        p_.message = msg;
        p_.elapsed = std::chrono::duration<double>(clock::now() - t0).count();
        finished_ = true;
    };

    std::string vf = "fps=" + fps_string(cfg.fps_num, cfg.fps_den);
    if (s_.scale) vf += ",scale=" + std::to_string(cfg.width) + ":" + std::to_string(cfg.height) + ":flags=bicubic";
    Process video;
    std::string err;
    if (!video.start({s_.ffmpeg, "-nostdin", "-v", "error", "-i", s_.input, "-map", "0:v:0", "-vf", vf, "-pix_fmt",
                      "yuv420p", "-f", "rawvideo", "-"},
                     vlog, &err)) {
        finish(true, false, err);
        return;
    }
    Process audio;
    if (cfg.audio_channels &&
        !audio.start({s_.ffmpeg, "-nostdin", "-v", "error", "-i", s_.input, "-map", "0:a:0", "-ac",
                      std::to_string(cfg.audio_channels), "-ar", std::to_string(cfg.audio_rate), "-f", "s16le",
                      "-acodec", "pcm_s16le", "-"},
                     alog, &err)) {
        video.kill();
        video.wait();
        finish(true, false, err);
        return;
    }

    FILE* out = fopen_utf8(s_.output, "wb");
    if (!out) {
        video.kill();
        audio.kill();
        finish(true, false, "cannot create " + s_.output);
        return;
    }
    tjc::Encoder enc;
    enc.init(cfg);
    std::vector<uint8_t> buf, frame(tjc::frame_size(cfg.width, cfg.height)), raw_audio;
    std::vector<int16_t> pcm;
    enc.write_stream_header(buf);
    bool write_ok = std::fwrite(buf.data(), 1, buf.size(), out) == buf.size();
    uint64_t bytes = buf.size(), audio_bytes = 0, dirty = 0, tiles = 0, padded = 0;
    bool audio_eof = false;
    auto last_preview = clock::now() - std::chrono::seconds(1);
    set_message("encoding");

    uint64_t frames = 0;
    while (write_ok && !cancel_) {
        if (video.read_full(frame.data(), frame.size()) < frame.size()) break;
        if (cfg.audio_channels && !audio_eof) {
            size_t need = enc.samples_for_next_frame(), queued = enc.queued_audio();
            if (need > queued) {
                size_t ch = cfg.audio_channels, want = (need - queued) * ch * 2;
                raw_audio.resize(want);
                size_t got = audio.read_full(raw_audio.data(), want);
                if (got < want) audio_eof = true;
                size_t samples = got / (ch * 2);
                pcm.resize(samples * ch);
                for (size_t i = 0; i < pcm.size(); ++i)
                    pcm[i] = int16_t(uint16_t(raw_audio[2 * i] | (raw_audio[2 * i + 1] << 8)));
                if (samples) enc.push_audio(pcm.data(), samples);
            }
        }
        if (s_.keyframe_every && frames && frames % s_.keyframe_every == 0) enc.force_keyframe();
        tjc::FrameStats st;
        buf.clear();
        enc.encode_frame(frame.data(), buf, &st);
        write_ok = std::fwrite(buf.data(), 1, buf.size(), out) == buf.size();
        ++frames;
        bytes += buf.size();
        audio_bytes += st.audio_bytes;
        dirty += st.dirty_tiles;
        tiles += st.total_tiles;
        padded += st.audio_padded;

        auto now = clock::now();
        bool want_preview = now - last_preview > std::chrono::milliseconds(250);
        std::vector<uint8_t> rgba;
        if (want_preview) {
            last_preview = now;
            std::vector<uint8_t> yuv(frame.size());
            enc.copy_recon(yuv.data());
            rgba.resize(size_t(cfg.width) * cfg.height * 4);
            yuv420_to_rgba(yuv.data(), cfg.width, cfg.height, ColorMatrix::Auto, rgba.data());
        }
        std::lock_guard<std::mutex> lk(mu_);
        p_.frames = frames;
        p_.bytes = bytes;
        p_.audio_bytes = audio_bytes;
        p_.dirty_percent = tiles ? 100.0 * double(dirty) / double(tiles) : 0;
        p_.audio_padded = padded;
        p_.elapsed = std::chrono::duration<double>(now - t0).count();
        if (want_preview) {
            preview_.swap(rgba);
            preview_new_ = true;
        }
    }

    // Anything left in the audio pipe is longer than the video.
    uint64_t extra = 0;
    if (cfg.audio_channels && !cancel_ && !audio_eof) {
        while (size_t n = audio.read(raw_audio.data(), raw_audio.size() ? raw_audio.size() : 1)) extra += n;
        extra /= uint64_t(cfg.audio_channels) * 2;
    }
    extra += enc.queued_audio();
    if (cancel_) {
        video.kill();
        audio.kill();
    }
    int vcode = video.wait();
    if (cfg.audio_channels) audio.wait();
    bool close_ok = std::fclose(out) == 0;
    std::string vtext = read_text_file(vlog), atext = read_text_file(alog);
    remove_utf8(vlog);
    remove_utf8(alog);

    if (cancel_) {
        remove_utf8(s_.output);
        finish(false, true, "cancelled; partial output deleted");
        return;
    }
    if (!write_ok || !close_ok) {
        finish(true, false, "writing " + s_.output + " failed (disk full?)");
        return;
    }
    if (frames == 0) {
        std::string why = last_lines(vtext, 3);
        finish(true, false, why.empty() ? "ffmpeg produced no frames" : why);
        return;
    }
    {
        std::lock_guard<std::mutex> lk(mu_);
        p_.audio_extra = extra;
    }
    std::string msg = "done: " + std::to_string(frames) + " frames";
    if (vcode != 0 && !vtext.empty()) msg += " (ffmpeg reported: " + last_lines(vtext, 1) + ")";
    if (cfg.audio_channels && padded)
        msg += "; audio ended early, " + std::to_string(double(padded) / cfg.audio_rate).substr(0, 5) + " s of silence added";
    if (cfg.audio_channels && extra > uint64_t(cfg.audio_rate) / 10)
        msg += "; " + std::to_string(double(extra) / cfg.audio_rate).substr(0, 5) + " s of extra audio dropped";
    if (!atext.empty() && cfg.audio_channels) msg += " (audio: " + last_lines(atext, 1) + ")";
    finish(false, false, msg);
}

}  // namespace gui
