#include "file_browser.h"

#include "imgui.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <filesystem>

namespace fs = std::filesystem;

namespace gui {

namespace {

std::string lower(std::string s) {
    for (char& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string human_size(unsigned long long b) {
    char buf[32];
    if (b >= (1ull << 30)) std::snprintf(buf, sizeof(buf), "%.1f GB", double(b) / double(1ull << 30));
    else if (b >= (1ull << 20)) std::snprintf(buf, sizeof(buf), "%.1f MB", double(b) / double(1ull << 20));
    else if (b >= 1024) std::snprintf(buf, sizeof(buf), "%.1f KB", double(b) / 1024.0);
    else std::snprintf(buf, sizeof(buf), "%llu B", b);
    return buf;
}

}  // namespace

void FileBrowser::open(const std::string& title, bool save, const std::vector<std::string>& exts,
                       const std::string& start_path) {
    title_ = title;
    save_ = save;
    exts_ = exts;
    show_all_ = exts.empty();
    name_.clear();
    std::error_code ec;
    fs::path p = fs::u8path(start_path);
    if (!start_path.empty() && fs::is_directory(p, ec)) {
        dir_ = p.u8string();
    } else if (!start_path.empty() && p.has_parent_path() && fs::is_directory(p.parent_path(), ec)) {
        dir_ = p.parent_path().u8string();
        if (save) name_ = p.filename().u8string();
    } else {
        dir_ = fs::current_path(ec).u8string();
    }
    refresh();
    visible_ = true;
    request_open_ = true;
}

void FileBrowser::refresh() {
    entries_.clear();
    error_.clear();
    selected_ = -1;
    std::error_code ec;
    for (fs::directory_iterator it(fs::u8path(dir_), fs::directory_options::skip_permission_denied, ec), end;
         !ec && it != end; it.increment(ec)) {
        std::string name = it->path().filename().u8string();
        if (name.empty() || name[0] == '.') continue;
        std::error_code e2;
        bool dir = it->is_directory(e2);
        if (!dir) {
            if (!show_all_ && !exts_.empty()) {
                std::string ext = lower(it->path().extension().u8string());
                if (std::find(exts_.begin(), exts_.end(), ext) == exts_.end()) continue;
            }
        }
        unsigned long long size = dir ? 0 : (unsigned long long)it->file_size(e2);
        entries_.push_back({name, dir, size});
    }
    if (ec) error_ = ec.message();
    std::sort(entries_.begin(), entries_.end(), [](const Entry& a, const Entry& b) {
        if (a.dir != b.dir) return a.dir;
        return lower(a.name) < lower(b.name);
    });
}

bool FileBrowser::draw(std::string& result) {
    if (!visible_) return false;
    if (request_open_) {
        ImGui::OpenPopup(title_.c_str());
        request_open_ = false;
    }
    bool done = false;
    ImVec2 vp = ImGui::GetMainViewport()->Size;
    ImGui::SetNextWindowSize(ImVec2(vp.x * 0.7f, vp.y * 0.75f), ImGuiCond_Appearing);
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    if (ImGui::BeginPopupModal(title_.c_str(), &visible_)) {
        if (ImGui::Button("Up")) {
            fs::path p = fs::u8path(dir_);
            if (p.has_parent_path() && p.parent_path() != p) {
                dir_ = p.parent_path().u8string();
                refresh();
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Home")) {
            const char* home = std::getenv("HOME");
            if (!home) home = std::getenv("USERPROFILE");
            if (home) {
                dir_ = home;
                refresh();
            }
        }
        ImGui::SameLine();
        char path_buf[4096];
        std::snprintf(path_buf, sizeof(path_buf), "%s", dir_.c_str());
        ImGui::SetNextItemWidth(-1);
        if (ImGui::InputText("##dir", path_buf, sizeof(path_buf), ImGuiInputTextFlags_EnterReturnsTrue)) {
            std::error_code ec;
            if (fs::is_directory(fs::u8path(path_buf), ec)) {
                dir_ = path_buf;
                refresh();
            }
        }

        float footer = ImGui::GetFrameHeightWithSpacing() * (save_ ? 2.2f : 1.2f);
        if (ImGui::BeginChild("##list", ImVec2(0, -footer), ImGuiChildFlags_Borders)) {
            if (!error_.empty()) ImGui::TextColored(ImVec4(1, 0.4f, 0.4f, 1), "%s", error_.c_str());
            for (int i = 0; i < int(entries_.size()); ++i) {
                const Entry& e = entries_[size_t(i)];
                std::string label = e.dir ? "[" + e.name + "]" : e.name;
                if (ImGui::Selectable(label.c_str(), selected_ == i, ImGuiSelectableFlags_AllowDoubleClick)) {
                    selected_ = i;
                    if (!e.dir) name_ = e.name;
                    if (ImGui::IsMouseDoubleClicked(0)) {
                        if (e.dir) {
                            dir_ = (fs::u8path(dir_) / fs::u8path(e.name)).u8string();
                            refresh();
                            break;
                        }
                        result = (fs::u8path(dir_) / fs::u8path(e.name)).u8string();
                        done = true;
                    }
                }
                if (!e.dir) {
                    ImGui::SameLine(ImGui::GetContentRegionAvail().x - 80);
                    ImGui::TextDisabled("%s", human_size(e.size).c_str());
                }
            }
        }
        ImGui::EndChild();

        if (save_) {
            char name_buf[1024];
            std::snprintf(name_buf, sizeof(name_buf), "%s", name_.c_str());
            ImGui::SetNextItemWidth(-1);
            if (ImGui::InputTextWithHint("##name", "file name", name_buf, sizeof(name_buf))) name_ = name_buf;
        }
        if (!exts_.empty() && ImGui::Checkbox("Show all files", &show_all_)) refresh();
        ImGui::SameLine();
        float bw = ImGui::CalcTextSize("Cancel").x + ImGui::GetStyle().FramePadding.x * 2;
        ImGui::SetCursorPosX(ImGui::GetWindowContentRegionMax().x - bw * 2.6f);
        bool can_ok = !name_.empty();
        if (!can_ok) ImGui::BeginDisabled();
        if (ImGui::Button(save_ ? "Save" : "Open", ImVec2(bw * 1.4f, 0))) {
            result = (fs::u8path(dir_) / fs::u8path(name_)).u8string();
            if (save_ && !exts_.empty() && fs::u8path(result).extension().empty()) result += exts_[0];
            done = true;
        }
        if (!can_ok) ImGui::EndDisabled();
        ImGui::SameLine();
        if (ImGui::Button("Cancel")) visible_ = false;
        if (done) visible_ = false;
        if (!visible_) ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    } else {
        visible_ = false;
    }
    return done;
}

}  // namespace gui
