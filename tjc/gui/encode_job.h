// Background encode: ffmpeg decodes any input to raw YUV420P + PCM, the frames go
// straight into tjc::Encoder. Also: probing an input with ffprobe.
#pragma once

#include "tjc.h"

#include <atomic>
#include <cstdint>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace gui {

struct MediaInfo {
    bool ok = false;
    std::string error;
    int width = 0, height = 0;
    uint32_t fps_num = 0, fps_den = 1;
    double duration = 0;  // seconds, 0 if unknown
    std::string video_codec;
    bool has_audio = false;
    int audio_rate = 0, audio_channels = 0;
    std::string audio_codec;
};

// Runs ffprobe (blocking; call from a worker thread).
MediaInfo probe_media(const std::string& ffprobe, const std::string& path);

struct EncodeSettings {
    std::string input, output;
    std::string ffmpeg = "ffmpeg";
    tjc::Config cfg;         // width/height/fps/audio filled in by the caller
    bool scale = false;      // resize to cfg.width x cfg.height
    uint32_t keyframe_every = 0;
    double duration = 0;     // from the probe, for progress
};

struct EncodeProgress {
    bool running = false, done = false, failed = false, cancelled = false;
    uint64_t frames = 0, total_estimate = 0;
    uint64_t bytes = 0, audio_bytes = 0;
    double elapsed = 0;
    double dirty_percent = 0;
    double inter_percent = 0;  // of the sent tiles, motion compensated
    uint64_t audio_padded = 0, audio_extra = 0;
    std::string message;
};

class EncodeJob {
public:
    ~EncodeJob();
    bool start(const EncodeSettings& s, std::string* err);
    void cancel();
    EncodeProgress progress() const;
    // Latest reconstruction preview (RGBA). Returns false if nothing new.
    bool take_preview(std::vector<uint8_t>& rgba, int& w, int& h);
    bool busy() const { return thread_.joinable() && !finished_; }

private:
    void run();
    void set_message(const std::string& m);

    EncodeSettings s_;
    std::thread thread_;
    std::atomic<bool> cancel_{false};
    std::atomic<bool> finished_{false};
    mutable std::mutex mu_;
    EncodeProgress p_;
    std::vector<uint8_t> preview_;
    bool preview_new_ = false;
};

}  // namespace gui
