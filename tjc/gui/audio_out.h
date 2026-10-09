// Audio output through miniaudio, fed from a single-producer ring buffer. The
// number of samples the device has consumed doubles as the playback clock.
#pragma once

#include <atomic>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace gui {

class AudioOut {
public:
    AudioOut();
    ~AudioOut();
    AudioOut(const AudioOut&) = delete;
    AudioOut& operator=(const AudioOut&) = delete;

    // Opens the default device. buffer_seconds sizes the ring buffer.
    bool open(int channels, uint32_t rate, double buffer_seconds, std::string* err);
    void close();
    bool is_open() const { return open_; }

    // Producer side (one thread): copies up to `frames` frames, returns how many fit.
    size_t write(const int16_t* interleaved, size_t frames);
    size_t free_frames() const;
    size_t buffered_frames() const;

    // Drops everything buffered and resets the clock to 0.
    void clear();
    void set_paused(bool p) { paused_ = p; }
    void set_volume(float v) { volume_ = v; }

    // Frames played since the last clear().
    uint64_t played() const { return played_.load(); }
    uint32_t rate() const { return rate_; }
    int channels() const { return channels_; }

    // Called from the device thread.
    void render(int16_t* out, uint32_t frames);

private:
    struct Device;
    std::unique_ptr<Device> dev_;
    bool open_ = false;
    int channels_ = 0;
    uint32_t rate_ = 0;
    std::vector<int16_t> ring_;
    size_t capacity_ = 0;  // frames
    std::atomic<uint64_t> read_{0}, write_{0}, played_{0};
    std::atomic<bool> paused_{true};
    std::atomic<float> volume_{1.0f};
    std::mutex render_lock_;  // held by render() (try_lock) and clear()
};

}  // namespace gui
