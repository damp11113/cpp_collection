// TJC playback engine: index thread + decode thread + audio clock.
//
// Seeking is exact. Every tile is self-contained, so the picture at frame T is
// the last copy of each tile sent at or before T. The index records, for every
// T, the oldest of those "last sent" frames (sync[T]); decoding from there to T
// reproduces frame T bit for bit without starting from the beginning.
#pragma once

#include "audio_out.h"
#include "yuv.h"
#include "tjc.h"

#include <atomic>
#include <condition_variable>
#include <cstdio>
#include <deque>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace gui {

struct VideoFrame {
    uint32_t index = 0;
    double pts = 0;
    uint64_t generation = 0;
    std::vector<uint8_t> rgba;   // width * height * 4
    std::vector<uint8_t> dirty;  // one flag per tile
    tjc::FrameStats stats;
};

class Player {
public:
    Player() = default;
    ~Player();
    Player(const Player&) = delete;
    Player& operator=(const Player&) = delete;

    bool open(const std::string& path, std::string* err);
    void close();
    bool is_open() const { return open_; }
    const std::string& path() const { return path_; }
    const tjc::StreamHeader& header() const { return header_; }
    const tjc::Layout& layout() const { return layout_; }
    double fps() const { return double(header_.fps_num) / double(header_.fps_den); }
    double pts_of(uint32_t frame) const { return double(frame) * header_.fps_den / header_.fps_num; }

    void play();
    void pause();
    bool playing() const { return playing_; }
    void seek(uint32_t frame);
    void step(int delta);
    void set_volume(float v) { audio_.set_volume(v); }
    void set_loop(bool loop) { loop_ = loop; }
    void set_matrix(ColorMatrix m);

    // Call once per UI redraw. Returns true when `out` holds a new frame to show
    // (the previous contents of `out` are recycled).
    bool update(VideoFrame& out);

    uint32_t current_frame() const { return current_; }
    uint32_t frames_indexed() const;
    bool index_complete() const;
    double index_fraction() const;
    bool index_entry(uint32_t frame, uint32_t* bytes, uint32_t* sync) const;
    bool seeking() const { return pending_show_; }
    bool has_audio() const { return audio_.is_open(); }
    const std::string& audio_status() const { return audio_status_; }
    std::string error() const;

private:
    void index_main();
    void decode_main();
    void do_seek(uint32_t target, uint64_t gen);
    bool decode_one();
    void recycle(VideoFrame&& f);

    static constexpr size_t kQueue = 8;

    std::string path_;
    bool open_ = false;
    tjc::StreamHeader header_;
    tjc::Layout layout_;
    int64_t file_size_ = 0;

    // Decode thread state.
    FILE* dfile_ = nullptr;
    tjc::Decoder dec_;
    uint32_t next_frame_ = 0;
    std::vector<uint8_t> yuv_;
    std::atomic<int> matrix_{0};

    // Index (index thread writes, others read under index_mu_).
    mutable std::mutex index_mu_;
    std::vector<uint64_t> offsets_;
    std::vector<uint32_t> sync_, sizes_;
    uint64_t index_bytes_ = 0;
    bool index_done_ = false;
    std::string index_error_;

    // Frame queue shared by decode thread and UI thread.
    mutable std::mutex mu_;
    std::condition_variable cv_;
    std::deque<VideoFrame> queue_;
    std::vector<VideoFrame> pool_;
    bool eof_ = false;
    std::string decode_error_;

    std::atomic<bool> stop_{false};
    std::atomic<bool> playing_{false};
    std::atomic<bool> loop_{false};
    std::atomic<int64_t> seek_req_{-1};
    std::atomic<uint64_t> gen_{0};

    // UI-thread clock state.
    bool pending_show_ = true;
    double base_pts_ = 0;
    double wall_ = 0;
    double last_update_ = -1;
    uint32_t current_ = 0;

    AudioOut audio_;
    std::string audio_status_;

    std::thread index_thread_, decode_thread_;
};

}  // namespace gui
