#define MINIAUDIO_IMPLEMENTATION
#define MA_NO_DECODING
#define MA_NO_ENCODING
#define MA_NO_GENERATION
#define MA_NO_RESOURCE_MANAGER
#define MA_NO_NODE_GRAPH
#define MA_NO_ENGINE
#include "miniaudio.h"

#include "audio_out.h"

#include <algorithm>
#include <cstring>

namespace gui {

struct AudioOut::Device {
    ma_device device;
};

namespace {
void data_callback(ma_device* d, void* out, const void*, ma_uint32 frames) {
    static_cast<AudioOut*>(d->pUserData)->render(static_cast<int16_t*>(out), frames);
}
}  // namespace

AudioOut::AudioOut() = default;

AudioOut::~AudioOut() { close(); }

bool AudioOut::open(int channels, uint32_t rate, double buffer_seconds, std::string* err) {
    close();
    channels_ = channels;
    rate_ = rate;
    capacity_ = std::max<size_t>(4096, size_t(buffer_seconds * rate));
    ring_.assign(capacity_ * size_t(channels), 0);
    read_ = write_ = played_ = 0;
    paused_ = true;

    dev_ = std::make_unique<Device>();
    ma_device_config cfg = ma_device_config_init(ma_device_type_playback);
    cfg.playback.format = ma_format_s16;
    cfg.playback.channels = ma_uint32(channels);
    cfg.sampleRate = rate;
    cfg.dataCallback = data_callback;
    cfg.pUserData = this;
    if (ma_device_init(nullptr, &cfg, &dev_->device) != MA_SUCCESS) {
        dev_.reset();
        if (err) *err = "no audio output device available";
        return false;
    }
    if (ma_device_start(&dev_->device) != MA_SUCCESS) {
        ma_device_uninit(&dev_->device);
        dev_.reset();
        if (err) *err = "could not start the audio device";
        return false;
    }
    open_ = true;
    return true;
}

void AudioOut::close() {
    if (dev_) {
        ma_device_uninit(&dev_->device);
        dev_.reset();
    }
    open_ = false;
}

size_t AudioOut::free_frames() const { return capacity_ - size_t(write_.load() - read_.load()); }

size_t AudioOut::buffered_frames() const { return size_t(write_.load() - read_.load()); }

size_t AudioOut::write(const int16_t* in, size_t frames) {
    if (!capacity_) return 0;
    uint64_t w = write_.load();
    size_t n = std::min(frames, capacity_ - size_t(w - read_.load()));
    const size_t ch = size_t(channels_);
    for (size_t i = 0; i < n;) {
        size_t pos = size_t((w + i) % capacity_);
        size_t run = std::min(n - i, capacity_ - pos);
        std::memcpy(&ring_[pos * ch], in + i * ch, run * ch * sizeof(int16_t));
        i += run;
    }
    write_.store(w + n);
    return n;
}

void AudioOut::clear() {
    std::lock_guard<std::mutex> lock(render_lock_);
    read_.store(write_.load());
    played_ = 0;
}

void AudioOut::render(int16_t* out, uint32_t frames) {
    const size_t ch = size_t(channels_);
    std::memset(out, 0, size_t(frames) * ch * sizeof(int16_t));
    if (paused_) return;
    std::unique_lock<std::mutex> lock(render_lock_, std::try_to_lock);
    if (!lock.owns_lock()) return;
    uint64_t r = read_.load();
    size_t n = std::min<size_t>(frames, size_t(write_.load() - r));
    float vol = volume_.load();
    for (size_t i = 0; i < n; ++i) {
        const int16_t* s = &ring_[size_t((r + i) % capacity_) * ch];
        for (size_t c = 0; c < ch; ++c) {
            float v = float(s[c]) * vol;
            out[i * ch + c] = int16_t(v > 32767.f ? 32767 : (v < -32768.f ? -32768 : v));
        }
    }
    read_.store(r + n);
    played_ += n;
}

}  // namespace gui
