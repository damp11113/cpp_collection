// TJC Studio: play and encode TJC streams. GLFW + Dear ImGui (OpenGL 2) +
// miniaudio; runs on Windows and Linux/X11.
//
//   tjc_studio [file]   .tjc files open in the player, anything else in the encoder

#include "encode_job.h"
#include "file_browser.h"
#include "platform.h"
#include "player.h"
#include "yuv.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"

#include <GLFW/glfw3.h>
#ifdef _WIN32
    #define GLFW_EXPOSE_NATIVE_WIN32
    #include <GLFW/glfw3native.h>
#endif

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <future>
#include <string>
#include <vector>

#ifndef GL_CLAMP_TO_EDGE
    #define GL_CLAMP_TO_EDGE 0x812F
#endif

namespace {

using namespace gui;

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

struct Texture {
    GLuint id = 0;
    int w = 0, h = 0;

    void upload(const uint8_t* rgba, int width, int height, bool nearest) {
        if (!id) glGenTextures(1, &id);
        glBindTexture(GL_TEXTURE_2D, id);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        GLint filter = nearest ? GL_NEAREST : GL_LINEAR;
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filter);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filter);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        if (width != w || height != h) {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgba);
            w = width;
            h = height;
        } else {
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, rgba);
        }
    }
    ImTextureID im() const { return ImTextureID(uintptr_t(id)); }
    void release() {
        if (id) glDeleteTextures(1, &id);
        id = 0;
        w = h = 0;
    }
};

std::string format_time(double s) {
    if (s < 0) s = 0;
    int total = int(s);
    char buf[32];
    if (total >= 3600)
        std::snprintf(buf, sizeof(buf), "%d:%02d:%05.2f", total / 3600, total / 60 % 60, s - (total / 60) * 60);
    else
        std::snprintf(buf, sizeof(buf), "%02d:%05.2f", total / 60, s - (total / 60) * 60);
    return buf;
}

std::string format_bytes(double b) {
    char buf[32];
    if (b >= 1e9) std::snprintf(buf, sizeof(buf), "%.2f GB", b / 1e9);
    else if (b >= 1e6) std::snprintf(buf, sizeof(buf), "%.2f MB", b / 1e6);
    else if (b >= 1e3) std::snprintf(buf, sizeof(buf), "%.1f KB", b / 1e3);
    else std::snprintf(buf, sizeof(buf), "%.0f B", b);
    return buf;
}

std::string lower_ext(const std::string& path) {
    std::string e = std::filesystem::u8path(path).extension().u8string();
    for (char& c : e) c = char(std::tolower(static_cast<unsigned char>(c)));
    return e;
}

std::string file_name(const std::string& path) { return std::filesystem::u8path(path).filename().u8string(); }

// "30", "25", "29.97", "30000/1001" -> rational. NTSC decimals map to x000/1001.
bool parse_fps(const char* s, uint32_t& num, uint32_t& den) {
    std::string str(s);
    size_t slash = str.find('/');
    if (slash != std::string::npos) {
        unsigned long n = std::strtoul(str.c_str(), nullptr, 10), d = std::strtoul(str.c_str() + slash + 1, nullptr, 10);
        if (!n || !d) return false;
        num = uint32_t(n);
        den = uint32_t(d);
        return true;
    }
    char* end = nullptr;
    double v = std::strtod(s, &end);
    if (end == s || *end || !(v > 0) || v > 100000) return false;
    for (double ntsc : {24.0, 30.0, 48.0, 60.0, 120.0}) {
        if (std::fabs(v - ntsc * 1000.0 / 1001.0) < 0.006) {
            num = uint32_t(ntsc * 1000);
            den = 1001;
            return true;
        }
    }
    if (v == std::floor(v)) {
        num = uint32_t(v);
        den = 1;
    } else {
        num = uint32_t(std::lround(v * 1000));
        den = 1000;
    }
    return true;
}

void help_marker(const char* text) {
    ImGui::SameLine();
    ImGui::TextDisabled("(?)");
    if (ImGui::BeginItemTooltip()) {
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 30.0f);
        ImGui::TextUnformatted(text);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

// Fits a w x h picture into the available region, centered. Returns top-left and size.
void fit_rect(int w, int h, ImVec2 avail, ImVec2 origin, ImVec2& p0, ImVec2& size) {
    float scale = std::min(avail.x / float(w), avail.y / float(h));
    size = ImVec2(float(w) * scale, float(h) * scale);
    p0 = ImVec2(origin.x + (avail.x - size.x) * 0.5f, origin.y + (avail.y - size.y) * 0.5f);
}

const std::vector<std::pair<std::string, std::string>> kTjcFilter = {{"TJC streams (*.tjc)", "*.tjc"},
                                                                     {"All files", "*.*"}};
const std::vector<std::pair<std::string, std::string>> kVideoFilter = {
    {"Video files", "*.mp4;*.mkv;*.mov;*.avi;*.webm;*.m4v;*.ts;*.flv;*.wmv;*.mpg;*.mpeg;*.y4m;*.gif"},
    {"All files", "*.*"}};

// ---------------------------------------------------------------------------
// Application state
// ---------------------------------------------------------------------------

enum class Browse { None, OpenTjc, EncodeInput, EncodeOutput };

struct App {
    GLFWwindow* window = nullptr;
    FileBrowser browser;
    Browse browse_target = Browse::None;
    int select_tab = -1;  // 0 = player, 1 = encoder

    // Player
    Player player;
    VideoFrame shown;
    bool have_frame = false;
    Texture video_tex, overlay_tex;
    std::vector<uint8_t> overlay_rgba;
    bool show_overlay = false, show_info = true, loop = false;
    int matrix = 0;
    float volume = 1.0f;
    std::string player_msg;

    // Encoder
    std::string enc_input, enc_output, ffmpeg_dir;
    MediaInfo info;
    std::future<MediaInfo> probe_future;
    std::string probed_path;
    int res_mode = 0, custom_w = 1280, custom_h = 720;
    int fps_mode = 0;
    char fps_text[32] = "30";
    int tile_mode = 1, tile_w = 16, tile_h = 16;
    int quality = 75, motion = 3, refresh_mode = 2, refresh_param = 30, diff_ref = 0;
    int keyframe_every = 0, threads = 0;
    bool use_motion = true, adaptive_huff = true;
    int effort = 1;
    float skip_invisible = 0.0f;
    bool audio_on = true;
    int audio_rate_mode = 0, audio_ch_mode = 0;
    EncodeJob job;
    Texture preview_tex;
    std::vector<uint8_t> preview_rgba;
    std::string enc_msg;
    std::string last_output;

    std::string tool(const char* name) const {
        if (ffmpeg_dir.empty()) return name;
        return (std::filesystem::u8path(ffmpeg_dir) / name).u8string();
    }
};

App* g_app = nullptr;

// ---------------------------------------------------------------------------
// File handling
// ---------------------------------------------------------------------------

void open_in_player(App& app, const std::string& path) {
    std::string err;
    app.have_frame = false;
    app.player_msg.clear();
    if (!app.player.open(path, &err)) {
        app.player_msg = err;
        return;
    }
    app.player.set_volume(app.volume);
    app.player.set_loop(app.loop);
    app.player.set_matrix(ColorMatrix(app.matrix));
    std::string title = "TJC Studio - " + file_name(path);
    glfwSetWindowTitle(app.window, title.c_str());
    app.select_tab = 0;
}

void set_encode_input(App& app, const std::string& path) {
    app.enc_input = path;
    std::filesystem::path p = std::filesystem::u8path(path);
    app.enc_output = (p.parent_path() / p.stem()).u8string() + ".tjc";
    app.select_tab = 1;
}

void open_path(App& app, const std::string& path) {
    if (lower_ext(path) == ".tjc") open_in_player(app, path);
    else set_encode_input(app, path);
}

void browse(App& app, Browse what) {
    const bool save = what == Browse::EncodeOutput;
    const auto& filters = what == Browse::EncodeInput ? kVideoFilter : kTjcFilter;
    std::string start = what == Browse::EncodeOutput ? app.enc_output
                        : what == Browse::EncodeInput ? app.enc_input
                                                      : app.player.path();
    if (has_native_dialogs()) {
        void* owner = nullptr;
#ifdef _WIN32
        owner = glfwGetWin32Window(app.window);
#endif
        std::string path = start;
        bool ok = save ? native_save_dialog(owner, "Save TJC stream as", filters, "tjc", path)
                       : native_open_dialog(owner, what == Browse::OpenTjc ? "Open TJC stream" : "Choose a video",
                                            filters, path);
        if (!ok) return;
        if (what == Browse::OpenTjc) open_in_player(app, path);
        else if (what == Browse::EncodeInput) set_encode_input(app, path);
        else app.enc_output = path;
        return;
    }
    std::vector<std::string> exts;
    if (what != Browse::EncodeInput) exts = {".tjc"};
    else exts = {".mp4", ".mkv", ".mov", ".avi", ".webm", ".m4v", ".ts", ".flv", ".wmv", ".mpg", ".mpeg", ".y4m", ".gif"};
    app.browser.open(save ? "Save TJC stream as" : (what == Browse::OpenTjc ? "Open TJC stream" : "Choose a video"),
                     save, exts, start);
    app.browse_target = what;
}

void drop_callback(GLFWwindow*, int count, const char** paths) {
    if (g_app && count > 0) open_path(*g_app, paths[0]);
}

// ---------------------------------------------------------------------------
// Player tab
// ---------------------------------------------------------------------------

void update_overlay(App& app) {
    const tjc::Layout& l = app.player.layout();
    app.overlay_rgba.assign(size_t(l.tiles) * 4, 0);
    for (int t = 0; t < l.tiles && size_t(t) < app.shown.dirty.size(); ++t) {
        uint8_t* p = &app.overlay_rgba[size_t(t) * 4];
        switch (app.shown.dirty[size_t(t)]) {
            case 1: p[0] = 255; p[1] = 60; p[2] = 40; p[3] = 110; break;   // intra
            case 2: p[0] = 60; p[1] = 230; p[2] = 90; p[3] = 90; break;    // motion compensated
            case 3: p[0] = 60; p[1] = 130; p[2] = 255; p[3] = 70; break;   // full refresh
            default: break;
        }
    }
    app.overlay_tex.upload(app.overlay_rgba.data(), l.cols, l.rows, true);
}

void player_keys(App& app) {
    if (!app.player.is_open() || ImGui::GetIO().WantTextInput) return;
    Player& p = app.player;
    double fps = p.fps();
    if (ImGui::IsKeyPressed(ImGuiKey_Space, false)) p.playing() ? p.pause() : p.play();
    if (ImGui::IsKeyPressed(ImGuiKey_RightArrow)) p.seek(p.current_frame() + uint32_t(std::lround(5 * fps)));
    if (ImGui::IsKeyPressed(ImGuiKey_LeftArrow)) p.step(-int(std::lround(5 * fps)));
    if (ImGui::IsKeyPressed(ImGuiKey_Period)) { p.pause(); p.step(1); }
    if (ImGui::IsKeyPressed(ImGuiKey_Comma)) { p.pause(); p.step(-1); }
    if (ImGui::IsKeyPressed(ImGuiKey_Home)) p.seek(0);
    if (ImGui::IsKeyPressed(ImGuiKey_T, false)) app.show_overlay = !app.show_overlay;
}

void draw_info_panel(App& app) {
    Player& p = app.player;
    const tjc::StreamHeader& h = p.header();
    const tjc::Layout& l = p.layout();
    ImGui::SeparatorText("Stream");
    ImGui::Text("TJC%u  %ux%u", h.version, h.width, h.height);
    ImGui::Text("%.3f fps (%u/%u)", p.fps(), h.fps_num, h.fps_den);
    ImGui::Text("Tiles %ux%u  (%d x %d = %d)", h.tile_w, h.tile_h, l.cols, l.rows, l.tiles);
    ImGui::Text("Quality %u", h.quality);
    const char* modes[] = {"none", "full periodic", "rolling"};
    ImGui::Text("Refresh %s / %u", modes[int(h.refresh_mode) % 3], h.refresh_param);
    ImGui::TextWrapped("Audio: %s", p.audio_status().c_str());
    if (h.version >= 3)
        ImGui::Text("Motion comp. %s, adaptive tables %s", (h.flags & tjc::kFlagMotion) ? "on" : "off",
                    (h.flags & tjc::kFlagAdaptiveHuffman) ? "on" : "off");

    uint32_t frames = p.frames_indexed();
    ImGui::SeparatorText("File");
    int64_t size = file_size_utf8(p.path());
    ImGui::Text("%s", format_bytes(double(size)).c_str());
    if (p.index_complete() && frames) {
        double secs = frames / p.fps();
        ImGui::Text("%u frames, %s", frames, format_time(secs).c_str());
        ImGui::Text("Average %.0f kbit/s", double(size) * 8 / secs / 1000);
    } else {
        ImGui::Text("Indexing... %u frames", frames);
    }

    if (app.have_frame) {
        const tjc::FrameStats& s = app.shown.stats;
        ImGui::SeparatorText("This frame");
        ImGui::Text("#%u  %s", app.shown.index, s.force_refresh ? "full refresh" : "update");
        ImGui::Text("Tiles %u / %u (%.1f%%)", s.dirty_tiles, s.total_tiles,
                    s.total_tiles ? 100.0 * s.dirty_tiles / s.total_tiles : 0.0);
        if (h.flags & tjc::kFlagMotion) ImGui::Text("Motion compensated %u", s.inter_tiles);
        if (s.tables) ImGui::TextUnformatted("New Huffman tables");
        ImGui::Text("Size %s", format_bytes(double(s.bytes)).c_str());
        if (h.audio_codec) ImGui::Text("Audio %u samples, %s", s.audio_samples, format_bytes(double(s.audio_bytes)).c_str());
        uint32_t sync = 0;
        if (p.index_entry(app.shown.index, nullptr, &sync))
            ImGui::Text("Seek restarts at #%u", sync);
    }

    // Frame sizes around the current position.
    if (frames) {
        ImGui::SeparatorText("Frame sizes (KB)");
        const int span = 150;
        uint32_t cur = p.current_frame();
        uint32_t first = cur > uint32_t(span / 2) ? cur - uint32_t(span / 2) : 0;
        static std::vector<float> vals;
        vals.clear();
        float vmax = 1;
        for (uint32_t f = first; f < first + span && f < frames; ++f) {
            uint32_t b = 0;
            p.index_entry(f, &b, nullptr);
            vals.push_back(float(b) / 1000.0f);
            vmax = std::max(vmax, vals.back());
        }
        char overlay[48];
        std::snprintf(overlay, sizeof(overlay), "peak %.0f KB", vmax);
        ImGui::PlotHistogram("##sizes", vals.data(), int(vals.size()), 0, overlay, 0, vmax * 1.1f,
                             ImVec2(-1, ImGui::GetFontSize() * 6));
    }
}

void draw_player(App& app) {
    Player& p = app.player;
    if (ImGui::Button("Open...")) browse(app, Browse::OpenTjc);
    ImGui::SameLine();
    if (p.is_open()) ImGui::TextUnformatted(file_name(p.path()).c_str());
    else ImGui::TextDisabled("Open a .tjc file, or drop one on the window");

    float right = ImGui::GetContentRegionMax().x;
    ImGui::SameLine(right - ImGui::GetFontSize() * 25);
    if (ImGui::Checkbox("Tile overlay", &app.show_overlay)) {}
    help_marker("Tints the tiles the current frame carried (T). Red = coded from scratch (intra), "
                "green = motion compensated, blue = full refresh frame.");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5);
    const char* matrices[] = {"Auto", "BT.601", "BT.709"};
    if (ImGui::Combo("##matrix", &app.matrix, matrices, 3)) p.set_matrix(ColorMatrix(app.matrix));
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("YUV to RGB matrix (Auto = BT.709 for 720p and up)");
    ImGui::SameLine();
    ImGui::Checkbox("Info", &app.show_info);

    // New frame from the engine?
    if (p.is_open() && p.update(app.shown)) {
        app.have_frame = true;
        app.video_tex.upload(app.shown.rgba.data(), p.layout().width, p.layout().height, false);
        update_overlay(app);
    }

    const float controls_h = ImGui::GetFrameHeightWithSpacing() * 2 + ImGui::GetStyle().ItemSpacing.y;
    const float info_w = app.show_info && p.is_open() ? ImGui::GetFontSize() * 19 : 0;
    ImVec2 avail = ImGui::GetContentRegionAvail();
    ImVec2 area(avail.x - info_w - (info_w > 0 ? ImGui::GetStyle().ItemSpacing.x : 0), avail.y - controls_h);
    if (ImGui::BeginChild("##video", area, ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar)) {
        ImVec2 o = ImGui::GetCursorScreenPos(), a = ImGui::GetContentRegionAvail();
        ImDrawList* dl = ImGui::GetWindowDrawList();
        dl->AddRectFilled(o, ImVec2(o.x + a.x, o.y + a.y), IM_COL32(0, 0, 0, 255));
        if (app.have_frame && a.x > 1 && a.y > 1) {
            const tjc::Layout& l = p.layout();
            ImVec2 p0, sz;
            fit_rect(l.width, l.height, a, o, p0, sz);
            ImVec2 p1(p0.x + sz.x, p0.y + sz.y);
            dl->AddImage(app.video_tex.im(), p0, p1);
            if (app.show_overlay) {
                ImVec2 uv1(float(l.width) / float(l.cols * l.tile_w), float(l.height) / float(l.rows * l.tile_h));
                dl->AddImage(app.overlay_tex.im(), p0, p1, ImVec2(0, 0), uv1);
            }
        } else if (!p.is_open()) {
            const char* t = "Drop a .tjc file here to play it, or a video to encode it";
            ImVec2 ts = ImGui::CalcTextSize(t);
            dl->AddText(ImVec2(o.x + (a.x - ts.x) * 0.5f, o.y + (a.y - ts.y) * 0.5f), IM_COL32(140, 140, 140, 255), t);
        }
        if (p.is_open() && p.seeking() && p.frames_indexed() > 0) {
            dl->AddText(ImVec2(o.x + 10, o.y + 10), IM_COL32(255, 220, 80, 255), "seeking...");
        }
        // Click on the picture toggles play.
        ImGui::InvisibleButton("##videoclick", a);
        if (ImGui::IsItemClicked() && p.is_open()) p.playing() ? p.pause() : p.play();
    }
    ImGui::EndChild();

    if (info_w > 0) {
        ImGui::SameLine();
        if (ImGui::BeginChild("##info", ImVec2(info_w, area.y), ImGuiChildFlags_Borders)) draw_info_panel(app);
        ImGui::EndChild();
    }

    // Transport controls.
    bool open = p.is_open();
    if (!open) ImGui::BeginDisabled();
    float bw = ImGui::GetFontSize() * 4.5f;
    if (ImGui::Button(p.playing() ? "Pause" : "Play", ImVec2(bw, 0))) p.playing() ? p.pause() : p.play();
    ImGui::SameLine();
    if (ImGui::Button("|<")) p.seek(0);
    ImGui::SameLine();
    if (ImGui::ArrowButton("##back", ImGuiDir_Left)) { p.pause(); p.step(-1); }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Previous frame (,)");
    ImGui::SameLine();
    if (ImGui::ArrowButton("##fwd", ImGuiDir_Right)) { p.pause(); p.step(1); }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Next frame (.)");
    ImGui::SameLine();

    uint32_t frames = open ? p.frames_indexed() : 0;
    int cur = int(p.current_frame());
    double total_s = open && frames ? frames / p.fps() : 0;
    std::string label = format_time(open ? p.pts_of(uint32_t(cur)) : 0) + " / " + format_time(total_s) +
                        (open && !p.index_complete() ? "+" : "");
    float tail = ImGui::GetFontSize() * 15;
    ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x - tail);
    int maxf = frames ? int(frames) - 1 : 0;
    if (ImGui::SliderInt("##seek", &cur, 0, maxf, label.c_str(), ImGuiSliderFlags_NoInput)) p.seek(uint32_t(cur));
    ImGui::SameLine();
    if (ImGui::Checkbox("Loop", &app.loop)) p.set_loop(app.loop);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x);
    if (ImGui::SliderFloat("##vol", &app.volume, 0.0f, 1.5f, "vol %.2f")) p.set_volume(app.volume);
    if (!open) ImGui::EndDisabled();

    // Status line.
    std::string err = open ? p.error() : std::string();
    if (!app.player_msg.empty()) ImGui::TextColored(ImVec4(1, 0.45f, 0.4f, 1), "%s", app.player_msg.c_str());
    else if (!err.empty()) ImGui::TextColored(ImVec4(1, 0.45f, 0.4f, 1), "%s", err.c_str());
    else if (open && !p.index_complete())
        ImGui::Text("Indexing for seeking: %.0f%%", p.index_fraction() * 100);
    else if (open)
        ImGui::TextDisabled("Space play/pause   Left/Right 5 s   , . frame step   Home start   T tile overlay   "
                            "click picture to play/pause");
    player_keys(app);
}

// ---------------------------------------------------------------------------
// Encoder tab
// ---------------------------------------------------------------------------

void poll_probe(App& app) {
    if (app.enc_input != app.probed_path && !app.enc_input.empty() && !app.probe_future.valid()) {
        app.probed_path = app.enc_input;
        app.info = MediaInfo();
        app.probe_future = std::async(std::launch::async, probe_media, app.tool("ffprobe"), app.enc_input);
    }
    if (app.probe_future.valid() &&
        app.probe_future.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
        app.info = app.probe_future.get();
        if (app.info.ok) {
            std::snprintf(app.fps_text, sizeof(app.fps_text), "%u/%u", app.info.fps_num, app.info.fps_den);
            app.audio_on = app.info.has_audio;
        }
    }
}

// Resolution the encoder will produce for the chosen preset.
void target_size(const App& app, int& w, int& h) {
    w = app.info.width;
    h = app.info.height;
    static const int heights[] = {0, 1080, 720, 480, 360};
    if (app.res_mode >= 1 && app.res_mode <= 4 && app.info.height > 0) {
        h = heights[app.res_mode];
        w = int(std::lround(double(app.info.width) * h / app.info.height / 2.0)) * 2;
    } else if (app.res_mode == 5) {
        w = app.custom_w;
        h = app.custom_h;
    }
}

bool build_settings(App& app, EncodeSettings& s, std::string& err) {
    if (app.enc_input.empty()) { err = "choose an input video"; return false; }
    if (app.enc_output.empty()) { err = "choose an output file"; return false; }
    if (!app.info.ok) { err = app.info.error.empty() ? "the input has not been probed yet" : app.info.error; return false; }
    s = EncodeSettings();
    s.input = app.enc_input;
    s.output = app.enc_output;
    s.ffmpeg = app.tool("ffmpeg");
    s.duration = app.info.duration;
    s.keyframe_every = uint32_t(std::max(0, app.keyframe_every));
    tjc::Config& c = s.cfg;
    int w, h;
    target_size(app, w, h);
    if (w < 1 || h < 1 || w > 65535 || h > 65535) { err = "invalid output size"; return false; }
    c.width = uint16_t(w);
    c.height = uint16_t(h);
    s.scale = w != app.info.width || h != app.info.height;
    if (app.fps_mode == 0) {
        c.fps_num = app.info.fps_num;
        c.fps_den = app.info.fps_den;
    } else if (!parse_fps(app.fps_text, c.fps_num, c.fps_den)) {
        err = "invalid frame rate (use 30, 29.97 or 30000/1001)";
        return false;
    }
    static const int tiles[] = {8, 16, 32, 48, 64};
    if (app.tile_mode < 5) c.tile_w = c.tile_h = uint8_t(tiles[app.tile_mode]);
    else { c.tile_w = uint8_t(std::clamp(app.tile_w, 1, 255)); c.tile_h = uint8_t(std::clamp(app.tile_h, 1, 255)); }
    c.quality = uint8_t(std::clamp(app.quality, 1, 100));
    c.motion_threshold_k = uint32_t(std::max(0, app.motion));
    c.refresh_mode = tjc::RefreshMode(app.refresh_mode);
    c.refresh_param = uint32_t(std::clamp(app.refresh_param, 0, 65535));
    c.diff_reference = tjc::DiffReference(app.diff_ref);
    c.threads = std::max(0, app.threads);
    c.motion = app.use_motion;
    c.adaptive_huffman = app.adaptive_huff;
    c.effort = std::clamp(app.effort, 0, 2);
    c.skip_invisible = double(std::max(0.0f, app.skip_invisible));
    if (app.audio_on && app.info.has_audio) {
        static const int rates[] = {0, 48000, 44100, 32000, 22050, 16000, 8000};
        c.audio_rate = uint32_t(app.audio_rate_mode == 0 ? app.info.audio_rate : rates[app.audio_rate_mode]);
        int ch = app.audio_ch_mode == 0 ? app.info.audio_channels : (app.audio_ch_mode == 1 ? 2 : 1);
        c.audio_channels = uint8_t(std::clamp(ch, 1, 8));
    }
    return true;
}

void draw_encoder(App& app) {
    poll_probe(app);
    EncodeProgress pr = app.job.progress();
    const bool busy = app.job.busy();
    const float label_w = ImGui::GetFontSize() * 7;
    const float btn_w = ImGui::GetFontSize() * 6;

    // Input / output
    if (busy) ImGui::BeginDisabled();
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("Input video");
    ImGui::SameLine(label_w);
    char buf[4096];
    std::snprintf(buf, sizeof(buf), "%s", app.enc_input.c_str());
    ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x - btn_w);
    if (ImGui::InputTextWithHint("##in", "any file ffmpeg can read (or drop it on the window)", buf, sizeof(buf),
                                 ImGuiInputTextFlags_EnterReturnsTrue))
        set_encode_input(app, buf);
    ImGui::SameLine();
    if (ImGui::Button("Browse##in", ImVec2(-1, 0))) browse(app, Browse::EncodeInput);

    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("Output .tjc");
    ImGui::SameLine(label_w);
    std::snprintf(buf, sizeof(buf), "%s", app.enc_output.c_str());
    ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x - btn_w);
    if (ImGui::InputText("##out", buf, sizeof(buf))) app.enc_output = buf;
    ImGui::SameLine();
    if (ImGui::Button("Browse##out", ImVec2(-1, 0))) browse(app, Browse::EncodeOutput);

    // Probe result
    ImGui::SetCursorPosX(label_w);
    if (app.probe_future.valid()) ImGui::TextDisabled("reading file info...");
    else if (app.info.ok) {
        std::string a = app.info.has_audio ? app.info.audio_codec + " " + std::to_string(app.info.audio_rate) + " Hz " +
                                                 std::to_string(app.info.audio_channels) + " ch"
                                           : "no audio";
        ImGui::TextDisabled("%s %dx%d @ %.3f fps, %s, %s", app.info.video_codec.c_str(), app.info.width,
                            app.info.height, double(app.info.fps_num) / app.info.fps_den,
                            format_time(app.info.duration).c_str(), a.c_str());
    } else if (!app.info.error.empty()) {
        ImGui::TextColored(ImVec4(1, 0.45f, 0.4f, 1), "%s", app.info.error.c_str());
    } else {
        ImGui::TextDisabled(" ");
    }

    // Settings
    if (ImGui::BeginTable("##settings", 2, ImGuiTableFlags_SizingStretchSame)) {
        ImGui::TableNextColumn();
        ImGui::SeparatorText("Video");
        const char* res[] = {"Source", "1080p", "720p", "480p", "360p", "Custom"};
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::Combo("Resolution", &app.res_mode, res, 6);
        if (app.res_mode == 5) {
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
            ImGui::InputInt("Width", &app.custom_w);
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
            ImGui::InputInt("Height", &app.custom_h);
        } else if (app.info.ok) {
            int w, h;
            target_size(app, w, h);
            ImGui::SameLine();
            ImGui::TextDisabled("%dx%d", w, h);
        }
        const char* fpsm[] = {"Source", "Custom"};
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::Combo("Frame rate", &app.fps_mode, fpsm, 2);
        if (app.fps_mode == 1) {
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
            ImGui::InputText("fps (30, 29.97, 30000/1001)", app.fps_text, sizeof(app.fps_text));
        }
        const char* tiles[] = {"8x8", "16x16", "32x32", "48x48", "64x64", "Custom"};
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::Combo("Tile size", &app.tile_mode, tiles, 6);
        help_marker("16x16 is the best default. Sizes that are not multiples of 16 waste bits on chroma padding; "
                    "below 8x8 the per-tile overhead grows fast.");
        if (app.tile_mode == 5) {
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
            ImGui::InputInt("Tile width", &app.tile_w);
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
            ImGui::InputInt("Tile height", &app.tile_h);
        }
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 12);
        ImGui::SliderInt("Quality", &app.quality, 1, 100);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 12);
        ImGui::SliderInt("Motion threshold", &app.motion, 0, 30);
        help_marker("A tile is re-sent when the mean difference in any of its 8x8 blocks exceeds this. "
                    "Raise it for noisy camera footage; 0 = any change.");

        ImGui::TableNextColumn();
        ImGui::SeparatorText("Refresh");
        const char* rm[] = {"None", "Full periodic", "Rolling"};
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::Combo("Mode", &app.refresh_mode, rm, 3);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::InputInt(app.refresh_mode == 1 ? "Every N frames" : "Tiles per frame", &app.refresh_param);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::InputInt("Keyframe every", &app.keyframe_every);
        help_marker("Also force a full refresh every N frames (0 = off).");
        const char* dr[] = {"Last coded", "Reconstructed"};
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::Combo("Diff against", &app.diff_ref, dr, 2);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::InputInt("Threads (0=all)", &app.threads);

        ImGui::SeparatorText("Compression (TJC3)");
        ImGui::Checkbox("Motion compensation", &app.use_motion);
        help_marker("Tiles can be predicted from the previous frame: big savings on pans and moving objects. "
                    "The player/decoder then keeps one extra frame in memory.");
        ImGui::Checkbox("Adaptive Huffman tables", &app.adaptive_huff);
        help_marker("Entropy tables fitted to the video, resent only when it pays off.");
        bool tjc2 = !app.use_motion && !app.adaptive_huff;
        if (ImGui::Checkbox("TJC2 compatible", &tjc2)) {
            app.use_motion = !tjc2;
            app.adaptive_huff = !tjc2;
        }
        help_marker("Turns both off: the file then plays on older TJC2 decoders (still ~5% smaller than before).");
        const char* efforts[] = {"Fast", "Normal", "Best"};
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::Combo("Effort", &app.effort, efforts, 3);
        help_marker("Encoder speed vs file size. Normal is within ~0.5% of Best at 1.5x its speed. "
                    "Never changes the file format.");
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::SliderFloat("Skip invisible", &app.skip_invisible, 0.0f, 4.0f, app.skip_invisible > 0 ? "%.1f" : "off");
        help_marker("Leave out tile updates that change no 8x8 block by more than this per pixel. "
                    "1-2 shrinks grainy video a lot; static areas keep their last grain pattern.");

        ImGui::SeparatorText("Audio");
        bool has_audio = app.info.ok && app.info.has_audio;
        if (!has_audio) ImGui::BeginDisabled();
        ImGui::Checkbox("Include audio (QOA)", &app.audio_on);
        const char* rates[] = {"Source", "48000", "44100", "32000", "22050", "16000", "8000"};
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::Combo("Sample rate", &app.audio_rate_mode, rates, 7);
        const char* chs[] = {"Source", "Stereo", "Mono"};
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        ImGui::Combo("Channels", &app.audio_ch_mode, chs, 3);
        if (!has_audio) ImGui::EndDisabled();
        ImGui::EndTable();
    }
    if (ImGui::CollapsingHeader("ffmpeg location")) {
        std::snprintf(buf, sizeof(buf), "%s", app.ffmpeg_dir.c_str());
        ImGui::SetNextItemWidth(-1);
        if (ImGui::InputTextWithHint("##ffdir", "folder with ffmpeg and ffprobe (empty = use PATH)", buf, sizeof(buf))) {
            app.ffmpeg_dir = buf;
            app.probed_path.clear();
        }
    }
    if (busy) ImGui::EndDisabled();

    // Start / cancel
    ImGui::Separator();
    if (!busy) {
        if (ImGui::Button("Start encoding", ImVec2(btn_w * 2, 0))) {
            EncodeSettings s;
            std::string err;
            if (build_settings(app, s, err) && app.job.start(s, &err)) {
                app.enc_msg.clear();
                app.last_output = s.output;
            } else {
                app.enc_msg = err;
            }
        }
    } else if (ImGui::Button("Cancel", ImVec2(btn_w * 2, 0))) {
        app.job.cancel();
    }
    ImGui::SameLine();
    if (!pr.running && pr.done && !app.last_output.empty() && ImGui::Button("Play output", ImVec2(btn_w * 2, 0)))
        open_in_player(app, app.last_output);

    // Progress
    double fps = app.info.ok ? double(app.info.fps_num) / app.info.fps_den : 30.0;
    if (pr.running || pr.done || pr.failed || pr.cancelled) {
        float frac = pr.total_estimate ? std::min(1.0f, float(pr.frames) / float(pr.total_estimate)) : 0.0f;
        if (pr.done) frac = 1.0f;
        char ov[64];
        std::snprintf(ov, sizeof(ov), "%llu / %llu frames", (unsigned long long)pr.frames,
                      (unsigned long long)pr.total_estimate);
        ImGui::ProgressBar(frac, ImVec2(-1, 0), ov);
        double speed = pr.elapsed > 0 ? double(pr.frames) / pr.elapsed : 0;
        double secs = pr.frames / fps;
        double eta = speed > 0 && pr.total_estimate > pr.frames ? double(pr.total_estimate - pr.frames) / speed : 0;
        ImGui::Text("%.1f fps (%.2fx realtime)   elapsed %s   ETA %s", speed, speed / fps,
                    format_time(pr.elapsed).c_str(), pr.running ? format_time(eta).c_str() : "-");
        ImGui::Text("Output %s   %.0f kbit/s   dirty tiles %.1f%% (%.0f%% motion comp.)   audio %.0f kbit/s",
                    format_bytes(double(pr.bytes)).c_str(), secs > 0 ? double(pr.bytes) * 8 / secs / 1000 : 0.0,
                    pr.dirty_percent, pr.inter_percent, secs > 0 ? double(pr.audio_bytes) * 8 / secs / 1000 : 0.0);
        ImVec4 col = pr.failed ? ImVec4(1, 0.45f, 0.4f, 1) : (pr.done ? ImVec4(0.5f, 0.9f, 0.5f, 1) : ImVec4(1, 1, 1, 1));
        ImGui::PushTextWrapPos(0);
        ImGui::TextColored(col, "%s", pr.message.c_str());
        ImGui::PopTextWrapPos();
    }
    if (!app.enc_msg.empty()) ImGui::TextColored(ImVec4(1, 0.45f, 0.4f, 1), "%s", app.enc_msg.c_str());

    // Preview of what the decoder will show.
    int pw, ph;
    if (app.job.take_preview(app.preview_rgba, pw, ph)) app.preview_tex.upload(app.preview_rgba.data(), pw, ph, false);
    if (app.preview_tex.id && (pr.running || pr.done)) {
        ImVec2 a = ImGui::GetContentRegionAvail();
        if (a.y > 40) {
            ImVec2 p0, sz;
            fit_rect(app.preview_tex.w, app.preview_tex.h, a, ImGui::GetCursorScreenPos(), p0, sz);
            ImDrawList* dl = ImGui::GetWindowDrawList();
            dl->AddImage(app.preview_tex.im(), p0, ImVec2(p0.x + sz.x, p0.y + sz.y));
            const char* label = "decoder preview";
            ImVec2 ts = ImGui::CalcTextSize(label);
            dl->AddRectFilled(ImVec2(p0.x, p0.y + sz.y - ts.y - 6), ImVec2(p0.x + ts.x + 12, p0.y + sz.y),
                              IM_COL32(0, 0, 0, 170));
            dl->AddText(ImVec2(p0.x + 6, p0.y + sz.y - ts.y - 3), IM_COL32(255, 255, 255, 220), label);
        }
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

void glfw_error(int code, const char* desc) { std::fprintf(stderr, "GLFW error %d: %s\n", code, desc); }

}  // namespace

int main(int argc, char** argv) {
    glfwSetErrorCallback(glfw_error);
    if (!glfwInit()) {
        std::fprintf(stderr, "tjc_studio: cannot initialize the windowing system (no display?)\n");
        return 1;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    GLFWwindow* window = glfwCreateWindow(1280, 800, "TJC Studio", nullptr, nullptr);
    if (!window) {
        std::fprintf(stderr, "tjc_studio: cannot create an OpenGL window\n");
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.IniFilename = nullptr;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    float xs = 1, ys = 1;
    glfwGetWindowContentScale(window, &xs, &ys);
    float scale = std::max(1.0f, std::max(xs, ys));
    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 0;
    style.FrameRounding = 3;
    style.GrabRounding = 3;
    style.ScaleAllSizes(scale);
    ImFontConfig font_cfg;
    font_cfg.SizePixels = std::round(13.0f * scale);
    io.Fonts->AddFontDefault(&font_cfg);
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();

    App app;
    app.window = window;
    g_app = &app;
    glfwSetDropCallback(window, drop_callback);
    std::vector<std::string> args = utf8_args(argc, argv);
    if (args.size() > 1) open_path(app, args[1]);

    while (!glfwWindowShouldClose(window)) {
        bool active = app.player.playing() || app.job.busy() || app.probe_future.valid() ||
                      (app.player.is_open() && (app.player.seeking() || !app.player.index_complete()));
        if (active) glfwPollEvents();
        else glfwWaitEventsTimeout(0.25);

        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        const ImGuiViewport* vp = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(vp->WorkPos);
        ImGui::SetNextWindowSize(vp->WorkSize);
        ImGui::Begin("##main", nullptr,
                     ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings |
                         ImGuiWindowFlags_NoBringToFrontOnFocus);
        if (ImGui::BeginTabBar("##tabs")) {
            ImGuiTabItemFlags f0 = app.select_tab == 0 ? ImGuiTabItemFlags_SetSelected : 0;
            ImGuiTabItemFlags f1 = app.select_tab == 1 ? ImGuiTabItemFlags_SetSelected : 0;
            app.select_tab = -1;
            if (ImGui::BeginTabItem("Player", nullptr, f0)) {
                draw_player(app);
                ImGui::EndTabItem();
            } else if (app.player.is_open()) {
                // Keep the clock and queue moving while the other tab is shown.
                if (app.player.update(app.shown)) {
                    app.have_frame = true;
                    app.video_tex.upload(app.shown.rgba.data(), app.player.layout().width,
                                         app.player.layout().height, false);
                    update_overlay(app);
                }
            }
            if (ImGui::BeginTabItem("Encoder", nullptr, f1)) {
                draw_encoder(app);
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }
        std::string picked;
        if (app.browser.draw(picked)) {
            if (app.browse_target == Browse::OpenTjc) open_in_player(app, picked);
            else if (app.browse_target == Browse::EncodeInput) set_encode_input(app, picked);
            else if (app.browse_target == Browse::EncodeOutput) app.enc_output = picked;
        }
        ImGui::End();

        ImGui::Render();
        int fw, fh;
        glfwGetFramebufferSize(window, &fw, &fh);
        glViewport(0, 0, fw, fh);
        glClearColor(0.08f, 0.08f, 0.09f, 1);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    app.job.cancel();
    app.player.close();
    app.video_tex.release();
    app.overlay_tex.release();
    app.preview_tex.release();
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
