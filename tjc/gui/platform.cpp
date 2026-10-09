#include "platform.h"

#include <cstring>
#include <filesystem>

#ifdef _WIN32
    #ifndef NOMINMAX
        #define NOMINMAX
    #endif
    #include <windows.h>
    #include <commdlg.h>
    #include <shellapi.h>
#else
    #include <cerrno>
    #include <csignal>
    #include <fcntl.h>
    #include <spawn.h>
    #include <sys/stat.h>
    #include <sys/wait.h>
    #include <unistd.h>
extern char** environ;
#endif

namespace gui {

#ifdef _WIN32

namespace {

std::wstring widen(const std::string& s) {
    if (s.empty()) return std::wstring();
    int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), nullptr, 0);
    std::wstring w(size_t(n), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), w.data(), n);
    return w;
}

std::string narrow(const wchar_t* w) {
    int n = WideCharToMultiByte(CP_UTF8, 0, w, -1, nullptr, 0, nullptr, nullptr);
    if (n <= 1) return std::string();
    std::string s(size_t(n - 1), '\0');
    WideCharToMultiByte(CP_UTF8, 0, w, -1, s.data(), n, nullptr, nullptr);
    return s;
}

// Quotes one argument following the rules of CommandLineToArgvW / the MSVC CRT.
void append_quoted(std::wstring& cmd, const std::wstring& arg) {
    if (!arg.empty() && arg.find_first_of(L" \t\n\v\"") == std::wstring::npos) {
        cmd += arg;
        return;
    }
    cmd += L'"';
    for (size_t i = 0;; ++i) {
        size_t backslashes = 0;
        while (i < arg.size() && arg[i] == L'\\') {
            ++i;
            ++backslashes;
        }
        if (i == arg.size()) {
            cmd.append(backslashes * 2, L'\\');
            break;
        }
        if (arg[i] == L'"') {
            cmd.append(backslashes * 2 + 1, L'\\');
            cmd += L'"';
        } else {
            cmd.append(backslashes, L'\\');
            cmd += arg[i];
        }
    }
    cmd += L'"';
}

std::wstring filter_string(const std::vector<std::pair<std::string, std::string>>& filters) {
    std::wstring f;
    for (const auto& p : filters) {
        f += widen(p.first);
        f += L'\0';
        f += widen(p.second);
        f += L'\0';
    }
    f += L'\0';
    return f;
}

}  // namespace

FILE* fopen_utf8(const std::string& path, const char* mode) {
    return _wfopen(widen(path).c_str(), widen(mode).c_str());
}

bool remove_utf8(const std::string& path) { return _wremove(widen(path).c_str()) == 0; }

int64_t file_size_utf8(const std::string& path) {
    WIN32_FILE_ATTRIBUTE_DATA d;
    if (!GetFileAttributesExW(widen(path).c_str(), GetFileExInfoStandard, &d)) return -1;
    return (int64_t(d.nFileSizeHigh) << 32) | d.nFileSizeLow;
}

bool Process::start(const std::vector<std::string>& args, const std::string& stderr_path, std::string* err) {
    if (started_ || args.empty()) return false;
    SECURITY_ATTRIBUTES sa = {sizeof(sa), nullptr, TRUE};
    HANDLE rd = nullptr, wr = nullptr;
    if (!CreatePipe(&rd, &wr, &sa, 1 << 20)) {
        if (err) *err = "CreatePipe failed";
        return false;
    }
    SetHandleInformation(rd, HANDLE_FLAG_INHERIT, 0);
    HANDLE in = CreateFileW(L"NUL", GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, &sa, OPEN_EXISTING, 0, nullptr);
    HANDLE errh = stderr_path.empty()
                      ? CreateFileW(L"NUL", GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, &sa, OPEN_EXISTING, 0, nullptr)
                      : CreateFileW(widen(stderr_path).c_str(), GENERIC_WRITE, FILE_SHARE_READ, &sa, CREATE_ALWAYS,
                                    FILE_ATTRIBUTE_NORMAL, nullptr);

    std::wstring cmd;
    for (size_t i = 0; i < args.size(); ++i) {
        if (i) cmd += L' ';
        append_quoted(cmd, widen(args[i]));
    }

    // Only these three handles go to the child, so a process started from another
    // thread at the same moment can't keep our pipe open.
    HANDLE inherit[3] = {in, wr, errh};
    SIZE_T attr_size = 0;
    InitializeProcThreadAttributeList(nullptr, 1, 0, &attr_size);
    std::vector<char> attr_buf(attr_size);
    auto attrs = reinterpret_cast<LPPROC_THREAD_ATTRIBUTE_LIST>(attr_buf.data());
    InitializeProcThreadAttributeList(attrs, 1, 0, &attr_size);
    UpdateProcThreadAttribute(attrs, 0, PROC_THREAD_ATTRIBUTE_HANDLE_LIST, inherit, sizeof(inherit), nullptr, nullptr);

    STARTUPINFOEXW si = {};
    si.StartupInfo.cb = sizeof(si);
    si.StartupInfo.dwFlags = STARTF_USESTDHANDLES;
    si.StartupInfo.hStdInput = in;
    si.StartupInfo.hStdOutput = wr;
    si.StartupInfo.hStdError = errh;
    si.lpAttributeList = attrs;
    PROCESS_INFORMATION pi = {};
    BOOL ok = CreateProcessW(nullptr, cmd.data(), nullptr, nullptr, TRUE,
                             CREATE_NO_WINDOW | EXTENDED_STARTUPINFO_PRESENT, nullptr, nullptr, &si.StartupInfo, &pi);
    DWORD code = GetLastError();
    DeleteProcThreadAttributeList(attrs);
    CloseHandle(wr);
    CloseHandle(in);
    CloseHandle(errh);
    if (!ok) {
        CloseHandle(rd);
        if (err) {
            *err = code == ERROR_FILE_NOT_FOUND ? "'" + args[0] + "' was not found (is it installed and in PATH?)"
                                                : "could not start '" + args[0] + "' (error " + std::to_string(code) + ")";
        }
        return false;
    }
    CloseHandle(pi.hThread);
    process_ = pi.hProcess;
    out_ = rd;
    started_ = true;
    return true;
}

size_t Process::read(void* buf, size_t n) {
    if (!out_) return 0;
    DWORD got = 0;
    DWORD want = DWORD(n > 0x40000000 ? 0x40000000 : n);
    if (!ReadFile(out_, buf, want, &got, nullptr)) return 0;
    return got;
}

void Process::kill() {
    if (process_) TerminateProcess(process_, 1);
}

int Process::wait() {
    if (!started_) return -1;
    if (out_) {
        CloseHandle(out_);
        out_ = nullptr;
    }
    DWORD code = DWORD(-1);
    WaitForSingleObject(process_, INFINITE);
    GetExitCodeProcess(process_, &code);
    CloseHandle(process_);
    process_ = nullptr;
    started_ = false;
    return int(code);
}

std::string temp_path(const std::string& name) {
    wchar_t dir[MAX_PATH + 1];
    DWORD n = GetTempPathW(MAX_PATH + 1, dir);
    std::string d = n ? narrow(dir) : std::string(".\\");
    return d + "tjc_gui_" + std::to_string(GetCurrentProcessId()) + "_" + name;
}

bool has_native_dialogs() { return true; }

bool native_open_dialog(void* owner, const char* title,
                        const std::vector<std::pair<std::string, std::string>>& filters, std::string& path) {
    wchar_t buf[32768] = {};
    std::wstring f = filter_string(filters), t = widen(title);
    OPENFILENAMEW ofn = {};
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = static_cast<HWND>(owner);
    ofn.lpstrFilter = f.c_str();
    ofn.lpstrFile = buf;
    ofn.nMaxFile = 32768;
    ofn.lpstrTitle = t.c_str();
    ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST | OFN_NOCHANGEDIR;
    if (!GetOpenFileNameW(&ofn)) return false;
    path = narrow(buf);
    return true;
}

bool native_save_dialog(void* owner, const char* title,
                        const std::vector<std::pair<std::string, std::string>>& filters, const char* default_ext,
                        std::string& path) {
    wchar_t buf[32768] = {};
    std::wstring init = widen(path);
    if (init.size() < 32767) std::memcpy(buf, init.c_str(), (init.size() + 1) * sizeof(wchar_t));
    std::wstring f = filter_string(filters), t = widen(title), ext = widen(default_ext);
    OPENFILENAMEW ofn = {};
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = static_cast<HWND>(owner);
    ofn.lpstrFilter = f.c_str();
    ofn.lpstrFile = buf;
    ofn.nMaxFile = 32768;
    ofn.lpstrTitle = t.c_str();
    ofn.lpstrDefExt = ext.c_str();
    ofn.Flags = OFN_EXPLORER | OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST | OFN_NOCHANGEDIR;
    if (!GetSaveFileNameW(&ofn)) return false;
    path = narrow(buf);
    return true;
}

std::vector<std::string> utf8_args(int argc, char** argv) {
    (void)argc;
    (void)argv;
    std::vector<std::string> out;
    int n = 0;
    LPWSTR* w = CommandLineToArgvW(GetCommandLineW(), &n);
    if (!w) return out;
    for (int i = 0; i < n; ++i) out.push_back(narrow(w[i]));
    LocalFree(w);
    return out;
}

#else  // POSIX

FILE* fopen_utf8(const std::string& path, const char* mode) { return std::fopen(path.c_str(), mode); }

bool remove_utf8(const std::string& path) { return std::remove(path.c_str()) == 0; }

int64_t file_size_utf8(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return -1;
    return int64_t(st.st_size);
}

bool Process::start(const std::vector<std::string>& args, const std::string& stderr_path, std::string* err) {
    if (started_ || args.empty()) return false;
    int fds[2];
    if (pipe2(fds, O_CLOEXEC) != 0) {
        if (err) *err = std::string("pipe failed: ") + std::strerror(errno);
        return false;
    }
    posix_spawn_file_actions_t fa;
    posix_spawn_file_actions_init(&fa);
    posix_spawn_file_actions_addopen(&fa, 0, "/dev/null", O_RDONLY, 0);
    posix_spawn_file_actions_adddup2(&fa, fds[1], 1);
    posix_spawn_file_actions_addopen(&fa, 2, stderr_path.empty() ? "/dev/null" : stderr_path.c_str(),
                                     O_WRONLY | O_CREAT | O_TRUNC, 0644);
    std::vector<char*> argv;
    for (const auto& a : args) argv.push_back(const_cast<char*>(a.c_str()));
    argv.push_back(nullptr);
    pid_t pid = -1;
    int rc = posix_spawnp(&pid, argv[0], &fa, nullptr, argv.data(), environ);
    posix_spawn_file_actions_destroy(&fa);
    close(fds[1]);
    if (rc != 0) {
        close(fds[0]);
        if (err) {
            *err = rc == ENOENT ? "'" + args[0] + "' was not found (is it installed and in PATH?)"
                                : "could not start '" + args[0] + "': " + std::strerror(rc);
        }
        return false;
    }
    pid_ = pid;
    out_ = fds[0];
    started_ = true;
    return true;
}

size_t Process::read(void* buf, size_t n) {
    if (out_ < 0) return 0;
    for (;;) {
        ssize_t r = ::read(out_, buf, n);
        if (r >= 0) return size_t(r);
        if (errno != EINTR) return 0;
    }
}

void Process::kill() {
    if (pid_ > 0) ::kill(pid_, SIGKILL);
}

int Process::wait() {
    if (!started_) return -1;
    if (out_ >= 0) {
        close(out_);
        out_ = -1;
    }
    int status = 0;
    while (waitpid(pid_, &status, 0) < 0 && errno == EINTR) {
    }
    pid_ = -1;
    started_ = false;
    return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
}

std::string temp_path(const std::string& name) {
    std::error_code ec;
    std::filesystem::path dir = std::filesystem::temp_directory_path(ec);
    if (ec) dir = "/tmp";
    return (dir / ("tjc_gui_" + std::to_string(getpid()) + "_" + name)).string();
}

bool has_native_dialogs() { return false; }

bool native_open_dialog(void*, const char*, const std::vector<std::pair<std::string, std::string>>&, std::string&) {
    return false;
}

bool native_save_dialog(void*, const char*, const std::vector<std::pair<std::string, std::string>>&, const char*,
                        std::string&) {
    return false;
}

std::vector<std::string> utf8_args(int argc, char** argv) { return std::vector<std::string>(argv, argv + argc); }

#endif

Process::~Process() {
    if (started_) {
        kill();
        wait();
    }
}

size_t Process::read_full(void* buf, size_t n) {
    size_t got = 0;
    auto* p = static_cast<uint8_t*>(buf);
    while (got < n) {
        size_t r = read(p + got, n - got);
        if (r == 0) break;
        got += r;
    }
    return got;
}

std::string read_text_file(const std::string& path, size_t max_bytes) {
    std::string s;
    FILE* f = fopen_utf8(path, "rb");
    if (!f) return s;
    s.resize(max_bytes);
    s.resize(std::fread(&s[0], 1, max_bytes, f));
    std::fclose(f);
    return s;
}

}  // namespace gui
