// Small platform layer for the GUI: UTF-8 file access, child processes with a
// stdout pipe (used to run ffmpeg/ffprobe), and native file dialogs.
#pragma once

#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

namespace gui {

// fopen() that takes a UTF-8 path on every platform.
FILE* fopen_utf8(const std::string& path, const char* mode);
bool remove_utf8(const std::string& path);
int64_t file_size_utf8(const std::string& path);

// Child process whose stdout we read. stdin is the null device, stderr goes to
// a file (or the null device), and on Windows no console window is shown.
class Process {
public:
    Process() = default;
    ~Process();
    Process(const Process&) = delete;
    Process& operator=(const Process&) = delete;

    // args[0] is looked up in PATH. Returns false and sets *err on failure.
    bool start(const std::vector<std::string>& args, const std::string& stderr_path, std::string* err);
    // Reads up to n bytes; returns 0 at end of output.
    size_t read(void* buf, size_t n);
    // Reads exactly n bytes unless the output ends first; returns bytes read.
    size_t read_full(void* buf, size_t n);
    void kill();
    // Waits for exit and returns the exit code (-1 if unknown).
    int wait();
    bool running() const { return started_; }

private:
    bool started_ = false;
#ifdef _WIN32
    void* process_ = nullptr;
    void* out_ = nullptr;
#else
    int pid_ = -1;
    int out_ = -1;
#endif
};

// Path of a fresh file in the temp directory (not created).
std::string temp_path(const std::string& name);
std::string read_text_file(const std::string& path, size_t max_bytes = 64 * 1024);

// Native file dialogs (Windows only; elsewhere these return false and the
// caller uses the built-in ImGui file browser). filters: {"description", "*.ext;*.ext2"}.
bool has_native_dialogs();
bool native_open_dialog(void* owner, const char* title,
                        const std::vector<std::pair<std::string, std::string>>& filters, std::string& path);
bool native_save_dialog(void* owner, const char* title,
                        const std::vector<std::pair<std::string, std::string>>& filters, const char* default_ext,
                        std::string& path);

// Command line as UTF-8 strings (on Windows taken from GetCommandLineW so that
// non-ASCII paths survive).
std::vector<std::string> utf8_args(int argc, char** argv);

}  // namespace gui
