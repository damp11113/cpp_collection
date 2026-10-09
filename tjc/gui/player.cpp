#include "player.h"

#include "platform.h"

#include <algorithm>
#include <chrono>
#include <cmath>

namespace gui {

namespace {

size_t file_read(void* user, void* dst, size_t n) { return std::fread(dst, 1, n, static_cast<FILE*>(user)); }

bool seek64(FILE* f, uint64_t off) {
#ifdef _WIN32
    return _fseeki64(f, (long long)off, SEEK_SET) == 0;
#else
    return fseeko(f, off_t(off), SEEK_SET) == 0;
#endif
}

double now_seconds() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

}  // namespace

Player::~Player() { close(); }

bool Player::open(const std::string& path, std::string* err) {
    close();
    dfile_ = fopen_utf8(path, "rb");
    if (!dfile_) {
        if (err) *err = "cannot open " + path;
        return false;
    }
    tjc::Status st = dec_.read_header(file_read, dfile_);
    if (st != tjc::Status::Ok) {
        if (err) *err = std::string("not a playable TJC stream: ") + tjc::status_string(st);
        std::fclose(dfile_);
        dfile_ = nullptr;
        return false;
    }
    path_ = path;
    header_ = dec_.header();
    layout_ = dec_.layout();
    units_ = dec_.unit_layout();
    file_size_ = file_size_utf8(path);
    yuv_.resize(tjc::frame_size(layout_.width, layout_.height));

    audio_status_.clear();
    if (header_.audio_codec) {
        std::string aerr;
        if (audio_.open(header_.audio_channels, header_.audio_rate, 1.0, &aerr))
            audio_status_ = "QOA " + std::to_string(header_.audio_rate) + " Hz, " +
                            std::to_string(header_.audio_channels) + " ch";
        else
            audio_status_ = "audio muted: " + aerr;
    } else {
        audio_status_ = "no audio track";
    }

    stop_ = false;
    playing_ = false;
    seek_req_ = -1;
    gen_ = 0;
    eof_ = false;
    next_frame_ = 0;
    pending_show_ = true;
    base_pts_ = 0;
    wall_ = 0;
    last_update_ = -1;
    current_ = 0;
    {
        std::lock_guard<std::mutex> lk(index_mu_);
        offsets_.clear();
        sync_.clear();
        tables_.clear();
        sizes_.clear();
        qualities_.clear();
        cum_.assign(1, 0);
        peak_kbps_ = 0;
        index_bytes_ = 0;
        index_done_ = false;
        index_error_.clear();
    }
    open_ = true;
    index_thread_ = std::thread(&Player::index_main, this);
    decode_thread_ = std::thread(&Player::decode_main, this);
    return true;
}

void Player::close() {
    if (!open_) return;
    stop_ = true;
    cv_.notify_all();
    if (index_thread_.joinable()) index_thread_.join();
    if (decode_thread_.joinable()) decode_thread_.join();
    audio_.close();
    if (dfile_) std::fclose(dfile_);
    dfile_ = nullptr;
    queue_.clear();
    pool_.clear();
    open_ = false;
    playing_ = false;
}

void Player::play() {
    if (!open_) return;
    std::lock_guard<std::mutex> lk(mu_);
    // Pressing play at the end starts over.
    if (eof_ && queue_.empty() && current_ + 1 >= frames_indexed() && index_complete()) {
        playing_ = true;
        audio_.set_paused(false);
        seek_req_ = 0;
        ++gen_;
        pending_show_ = true;
        base_pts_ = 0;
        wall_ = 0;
        cv_.notify_all();
        return;
    }
    playing_ = true;
    audio_.set_paused(false);
}

void Player::pause() {
    playing_ = false;
    audio_.set_paused(true);
}

void Player::seek(uint32_t frame) {
    if (!open_) return;
    uint32_t known = frames_indexed();
    if (known == 0) return;
    frame = std::min(frame, known - 1);
    std::lock_guard<std::mutex> lk(mu_);
    ++gen_;
    seek_req_ = int64_t(frame);
    pending_show_ = true;
    base_pts_ = pts_of(frame);
    wall_ = 0;
    current_ = frame;
    audio_.clear();
    cv_.notify_all();
}

void Player::step(int delta) {
    int64_t t = int64_t(current_) + delta;
    seek(uint32_t(std::max<int64_t>(0, t)));
}

void Player::set_matrix(ColorMatrix m) {
    if (matrix_.exchange(int(m)) != int(m) && open_ && !playing_) seek(current_);
}

uint32_t Player::frames_indexed() const {
    std::lock_guard<std::mutex> lk(index_mu_);
    return uint32_t(offsets_.size());
}

bool Player::index_complete() const {
    std::lock_guard<std::mutex> lk(index_mu_);
    return index_done_;
}

double Player::index_fraction() const {
    std::lock_guard<std::mutex> lk(index_mu_);
    if (index_done_) return 1.0;
    return file_size_ > 0 ? double(index_bytes_) / double(file_size_) : 0.0;
}

bool Player::index_entry(uint32_t frame, uint32_t* bytes, uint32_t* sync, uint8_t* quality) const {
    std::lock_guard<std::mutex> lk(index_mu_);
    if (frame >= offsets_.size()) return false;
    if (bytes) *bytes = sizes_[frame];
    if (sync) *sync = sync_[frame];
    if (quality) *quality = qualities_[frame];
    return true;
}

double Player::bitrate_kbps(uint32_t frame, uint32_t window) const {
    std::lock_guard<std::mutex> lk(index_mu_);
    if (frame >= offsets_.size() || window == 0) return 0;
    uint32_t first = frame + 1 >= window ? frame + 1 - window : 0;
    double secs = double(frame + 1 - first) * header_.fps_den / header_.fps_num;
    return double(cum_[frame + 1] - cum_[first]) * 8.0 / secs / 1000.0;
}

double Player::peak_kbps() const {
    std::lock_guard<std::mutex> lk(index_mu_);
    return peak_kbps_;
}

std::string Player::error() const {
    {
        std::lock_guard<std::mutex> lk(mu_);
        if (!decode_error_.empty()) return decode_error_;
    }
    std::lock_guard<std::mutex> lk(index_mu_);
    return index_error_;
}

void Player::index_main() {
    FILE* f = fopen_utf8(path_, "rb");
    tjc::Decoder dec;
    if (!f || dec.read_header(file_read, f) != tjc::Status::Ok) {
        std::lock_guard<std::mutex> lk(index_mu_);
        index_error_ = "cannot index the stream";
        index_done_ = true;
        if (f) std::fclose(f);
        return;
    }
    uint64_t off = header_.version == 1 ? tjc::kStreamHeaderSizeV1 : tjc::kStreamHeaderSize;
    // The tracker follows what each tile's content depends on (motion vectors
    // included) and gives the frame a seek has to start decoding from.
    tjc::SyncTracker tracker;
    tracker.reset(dec.unit_layout());
    int64_t last_tables = -1;
    const uint32_t second = uint32_t(std::max(1.0, std::round(double(header_.fps_num) / header_.fps_den)));
    for (uint32_t frame = 0; !stop_; ++frame) {
        tjc::FrameStats fs;
        tjc::Status st = dec.decode_frame(file_read, f, &fs, false);
        if (st == tjc::Status::EndOfStream) break;
        if (st != tjc::Status::Ok) {
            std::lock_guard<std::mutex> lk(index_mu_);
            index_error_ = std::string("stream damaged after frame ") + std::to_string(frame) + ": " +
                           tjc::status_string(st);
            break;
        }
        uint32_t sync = tracker.add_frame(frame, dec.unit_info());
        if (fs.tables) last_tables = frame;
        {
            std::lock_guard<std::mutex> lk(index_mu_);
            offsets_.push_back(off);
            sync_.push_back(sync);
            tables_.push_back(int32_t(last_tables));
            sizes_.push_back(uint32_t(fs.bytes));
            qualities_.push_back(fs.quality);
            cum_.push_back(cum_.back() + fs.bytes);
            if (frame + 1 >= second)  // busiest full one-second window
                peak_kbps_ = std::max(peak_kbps_, double(cum_[frame + 1] - cum_[frame + 1 - second]) * 8.0 /
                                                      (double(second) * header_.fps_den / header_.fps_num) / 1000.0);
            index_bytes_ = off + fs.bytes;
        }
        off += fs.bytes;
    }
    std::fclose(f);
    std::lock_guard<std::mutex> lk(index_mu_);
    index_done_ = true;
    if (peak_kbps_ == 0 && !offsets_.empty())  // shorter than a second: the average
        peak_kbps_ = double(cum_.back()) * 8.0 / (double(offsets_.size()) * header_.fps_den / header_.fps_num) / 1000.0;
}

void Player::recycle(VideoFrame&& f) {
    if (pool_.size() < kQueue + 2) pool_.push_back(std::move(f));
}

bool Player::decode_one() {
    VideoFrame f;
    {
        std::lock_guard<std::mutex> lk(mu_);
        if (!pool_.empty()) {
            f = std::move(pool_.back());
            pool_.pop_back();
        }
    }
    uint64_t gen = gen_.load();
    tjc::Status st = dec_.decode_frame(file_read, dfile_, &f.stats, true);
    if (st != tjc::Status::Ok) {
        std::lock_guard<std::mutex> lk(mu_);
        eof_ = true;
        if (st != tjc::Status::EndOfStream)
            decode_error_ = std::string("frame ") + std::to_string(next_frame_) + ": " + tjc::status_string(st);
        return false;
    }
    f.index = next_frame_++;
    f.pts = pts_of(f.index);
    f.generation = gen;
    dec_.copy_frame(yuv_.data());
    f.rgba.resize(size_t(layout_.width) * size_t(layout_.height) * 4);
    yuv420_to_rgba(yuv_.data(), layout_.width, layout_.height, ColorMatrix(matrix_.load()), f.rgba.data());
    // Per unit: 0 = unchanged, 1 = intra, 2 = motion compensated, 3 = full refresh.
    f.dirty.resize(size_t(units_.tiles));
    const tjc::TileInfo* info = dec_.unit_info();
    for (int t = 0; t < units_.tiles; ++t)
        f.dirty[size_t(t)] = f.stats.force_refresh && info[t].mode ? 3 : info[t].mode;

    if (audio_.is_open() && dec_.audio_samples()) {
        const int16_t* a = dec_.audio();
        size_t left = dec_.audio_samples();
        while (left && !stop_ && seek_req_.load() < 0) {
            size_t n = audio_.write(a, left);
            a += n * size_t(header_.audio_channels);
            left -= n;
            if (left) std::this_thread::sleep_for(std::chrono::milliseconds(4));
        }
        if (left) return false;  // interrupted by a seek or close
    }
    std::lock_guard<std::mutex> lk(mu_);
    if (gen != gen_.load()) {
        recycle(std::move(f));
        return false;
    }
    queue_.push_back(std::move(f));
    return true;
}

void Player::do_seek(uint32_t target, uint64_t gen) {
    {
        std::lock_guard<std::mutex> lk(mu_);
        while (!queue_.empty()) {
            recycle(std::move(queue_.front()));
            queue_.pop_front();
        }
        eof_ = false;
        decode_error_.clear();
    }
    audio_.clear();
    uint64_t off, tables_off = 0;
    uint32_t start;
    bool need_tables = false;
    {
        std::lock_guard<std::mutex> lk(index_mu_);
        if (offsets_.empty()) return;
        target = std::min<uint32_t>(target, uint32_t(offsets_.size() - 1));
        start = sync_[target];
        off = offsets_[start];
        int32_t tf = tables_[start];
        if (tf >= 0 && uint32_t(tf) < start) {  // Huffman tables in force at the start frame
            need_tables = true;
            tables_off = offsets_[size_t(tf)];
        }
    }
    if (need_tables) {
        if (!seek64(dfile_, tables_off)) return;
        dec_.decode_frame(file_read, dfile_, nullptr, false);
    } else {
        dec_.reset_tables();  // the start frame still uses the standard tables
    }
    if (!seek64(dfile_, off)) return;
    // Bring every tile up to date; only the target frame is shown.
    for (uint32_t f = start; f < target; ++f) {
        if (stop_ || gen != gen_.load()) return;
        if (dec_.decode_frame(file_read, dfile_, nullptr, true) != tjc::Status::Ok) break;
    }
    next_frame_ = target;
    audio_.clear();
}

void Player::decode_main() {
    while (!stop_) {
        int64_t req = seek_req_.exchange(-1);
        if (req >= 0) {
            do_seek(uint32_t(req), gen_.load());
            continue;
        }
        {
            std::unique_lock<std::mutex> lk(mu_);
            cv_.wait_for(lk, std::chrono::milliseconds(20), [&] {
                return stop_ || seek_req_.load() >= 0 || (!eof_ && queue_.size() < kQueue);
            });
            if (stop_ || seek_req_.load() >= 0 || eof_ || queue_.size() >= kQueue) continue;
        }
        decode_one();
    }
}

bool Player::update(VideoFrame& out) {
    if (!open_) return false;
    double now = now_seconds();
    double dt = last_update_ < 0 ? 0 : now - last_update_;
    last_update_ = now;

    std::unique_lock<std::mutex> lk(mu_);
    uint64_t gen = gen_.load();
    while (!queue_.empty() && queue_.front().generation != gen) {
        recycle(std::move(queue_.front()));
        queue_.pop_front();
    }

    bool got = false;
    auto take_front = [&] {
        std::swap(out, queue_.front());
        recycle(std::move(queue_.front()));
        queue_.pop_front();
        current_ = out.index;
        got = true;
    };

    if (pending_show_) {
        // First frame after open or a seek: show it right away, then run the clock from it.
        if (!queue_.empty()) {
            take_front();
            pending_show_ = false;
            base_pts_ = out.pts;
            wall_ = 0;
        }
    } else if (playing_) {
        if (!audio_.is_open()) wall_ += dt;
        double clock = base_pts_ + (audio_.is_open() ? double(audio_.played()) / audio_.rate() : wall_);
        while (!queue_.empty() && queue_.front().pts <= clock) take_front();
        bool drained = !audio_.is_open() || audio_.buffered_frames() == 0;
        if (queue_.empty() && eof_ && drained) {
            if (loop_ && frames_indexed() > 1) {
                lk.unlock();
                seek(0);
                return got;
            }
            playing_ = false;
            audio_.set_paused(true);
        }
    }
    cv_.notify_all();
    return got;
}

}  // namespace gui
