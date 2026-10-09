// Minimal ImGui file browser, used where no native dialog is available (Linux/X11).
#pragma once

#include <string>
#include <vector>

namespace gui {

class FileBrowser {
public:
    // exts: lower-case extensions including the dot (".tjc"); empty = all files.
    void open(const std::string& title, bool save, const std::vector<std::string>& exts,
              const std::string& start_path);
    // Draw every frame. Returns true once when the user confirmed; `result` then holds the path.
    bool draw(std::string& result);
    bool visible() const { return visible_; }

private:
    struct Entry {
        std::string name;
        bool dir;
        unsigned long long size;
    };
    void refresh();

    std::string title_, dir_, name_, error_;
    std::vector<std::string> exts_;
    std::vector<Entry> entries_;
    bool save_ = false, visible_ = false, request_open_ = false, show_all_ = false;
    int selected_ = -1;
};

}  // namespace gui
