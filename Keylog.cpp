#include <windows.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <thread>

class KeyboardTracker {
private:
    std::map<int, int> keyCount;
    std::map<std::string, int> keySequences;
    std::string lastKey;
    HHOOK hKeyboardHook;
    bool isRunning;
    std::chrono::steady_clock::time_point lastSaveTime;
    int saveIntervalMinutes;
    
    // Virtual key code to readable string mapping
    std::map<int, std::string> vkToString = {
        {VK_SPACE, "SPACE"}, {VK_RETURN, "ENTER"}, {VK_BACK, "BACKSPACE"},
        {VK_TAB, "TAB"}, {VK_SHIFT, "SHIFT"}, {VK_CONTROL, "CTRL"},
        {VK_MENU, "ALT"}, {VK_ESCAPE, "ESC"}, {VK_DELETE, "DELETE"},
        {VK_INSERT, "INSERT"}, {VK_HOME, "HOME"}, {VK_END, "END"},
        {VK_PRIOR, "PAGEUP"}, {VK_NEXT, "PAGEDOWN"}, {VK_UP, "UP"},
        {VK_DOWN, "DOWN"}, {VK_LEFT, "LEFT"}, {VK_RIGHT, "RIGHT"},
        {VK_F1, "F1"}, {VK_F2, "F2"}, {VK_F3, "F3"}, {VK_F4, "F4"},
        {VK_F5, "F5"}, {VK_F6, "F6"}, {VK_F7, "F7"}, {VK_F8, "F8"},
        {VK_F9, "F9"}, {VK_F10, "F10"}, {VK_F11, "F11"}, {VK_F12, "F12"}
    };

public:
    KeyboardTracker(int saveInterval = 5) : saveIntervalMinutes(saveInterval), isRunning(false) {
        lastSaveTime = std::chrono::steady_clock::now();
    }

    std::string getCurrentDateTime() {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&time_t), "%Y%m%d_%H%M%S");
        return ss.str();
    }

    std::string vkToReadableString(int vk) {
        if (vkToString.find(vk) != vkToString.end()) {
            return vkToString[vk];
        }
        
        // Handle alphanumeric keys
        if (vk >= 'A' && vk <= 'Z') {
            return std::string(1, (char)vk);
        }
        if (vk >= '0' && vk <= '9') {
            return std::string(1, (char)vk);
        }
        
        // Handle some special characters
        switch (vk) {
            case VK_OEM_1: return "SEMICOLON";      // ;:
            case VK_OEM_PLUS: return "PLUS";        // =+
            case VK_OEM_COMMA: return "COMMA";      // ,<
            case VK_OEM_MINUS: return "MINUS";      // -_
            case VK_OEM_PERIOD: return "PERIOD";    // .>
            case VK_OEM_2: return "SLASH";          // /?
            case VK_OEM_3: return "GRAVE";          // `~
            case VK_OEM_4: return "LBRACKET";       // [{
            case VK_OEM_5: return "BACKSLASH";      // \|
            case VK_OEM_6: return "RBRACKET";       // ]}
            case VK_OEM_7: return "QUOTE";          // '"
            default: return "KEY_" + std::to_string(vk);
        }
    }

    static LRESULT CALLBACK KeyboardProc(int nCode, WPARAM wParam, LPARAM lParam) {
        if (nCode >= 0 && wParam == WM_KEYDOWN) {
            KBDLLHOOKSTRUCT* pKeyboard = (KBDLLHOOKSTRUCT*)lParam;
            KeyboardTracker* tracker = reinterpret_cast<KeyboardTracker*>(
                GetWindowLongPtr(GetConsoleWindow(), GWLP_USERDATA));
            
            if (tracker) {
                tracker->recordKeyPress(pKeyboard->vkCode);
            }
        }
        return CallNextHookEx(NULL, nCode, wParam, lParam);
    }

    void recordKeyPress(int vk) {
        std::string keyName = vkToReadableString(vk);
        
        // Count individual keys
        keyCount[vk]++;
        
        // Track key sequences (bigrams)
        if (!lastKey.empty()) {
            std::string sequence = lastKey + "->" + keyName;
            keySequences[sequence]++;
        }
        lastKey = keyName;
        
        // Check if it's time to save
        auto now = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::minutes>(now - lastSaveTime);
        if (duration.count() >= saveIntervalMinutes) {
            saveData();
            lastSaveTime = now;
        }
    }

    void saveData() {
        std::string filename = std::string("logs/") + "keyboard_data_" + getCurrentDateTime() + ".json";
        std::ofstream file(filename);
        
        if (!file.is_open()) {
            return;
        }
        
        file << "{\n";
        file << "  \"timestamp\": \"" << getCurrentDateTime() << "\",\n";
        file << "  \"save_interval_minutes\": " << saveIntervalMinutes << ",\n";
        
        // Save key counts
        file << "  \"key_counts\": {\n";
        bool first = true;
        for (const auto& pair : keyCount) {
            if (!first) file << ",\n";
            file << "    \"" << vkToReadableString(pair.first) << "\": " << pair.second;
            first = false;
        }
        file << "\n  },\n";
        
        // Save key sequences
        file << "  \"key_sequences\": {\n";
        first = true;
        for (const auto& pair : keySequences) {
            if (!first) file << ",\n";
            file << "    \"" << pair.first << "\": " << pair.second;
            first = false;
        }
        file << "\n  }\n";
        
        file << "}\n";
        file.close();
        
        // Clear data after saving to avoid double counting
        keyCount.clear();
        keySequences.clear();
        lastKey.clear();
    }

    bool startTracking() {
        if (isRunning) return false;
        
        // Hide console window
        HWND consoleWindow = GetConsoleWindow();
        ShowWindow(consoleWindow, SW_HIDE);
        
        // Store pointer to this instance for the hook procedure
        SetWindowLongPtr(consoleWindow, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
        
        // Install keyboard hook
        hKeyboardHook = SetWindowsHookEx(WH_KEYBOARD_LL, KeyboardProc, GetModuleHandle(NULL), 0);
        
        if (hKeyboardHook == NULL) {
            return false;
        }
        
        isRunning = true;
        return true;
    }

    void stopTracking() {
        if (isRunning) {
            UnhookWindowsHookEx(hKeyboardHook);
            isRunning = false;
            saveData(); // Save final data before stopping
        }
    }

    void messageLoop() {
        MSG msg;
        while (isRunning && GetMessage(&msg, NULL, 0, 0)) {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }
};

int main() {
    std::cout << "Keyboard Usage Tracker\n";
    std::cout << "=====================\n";
    std::cout << "This tool will track keyboard usage patterns to help create custom layouts.\n";
    std::cout << "Data will be saved every 1 minutes in JSON format.\n\n";
    
    int interval = 1;
    
    std::cout << "Starting keyboard tracking...\n";
    std::cout << "The application will now hide and run in the background.\n";
    std::cout << "Data files will be saved as 'keyboard_data_YYYYMMDD_HHMMSS.json'\n\n";
    
	std::cout << "This window will close in 3 sec.\n\nPress Ctrl+C to stop";
	
    // Wait a moment before hiding
    std::this_thread::sleep_for(std::chrono::seconds(3));
    
    KeyboardTracker tracker(interval);
    
    if (!tracker.startTracking()) {
        std::cout << "Failed to start keyboard tracking. Please run as administrator.\n";
        system("pause");
        return 1;
    }
    
    // Set up Ctrl+C handler
    SetConsoleCtrlHandler([](DWORD ctrlType) -> BOOL {
        if (ctrlType == CTRL_C_EVENT) {
            PostQuitMessage(0);
            return TRUE;
        }
        return FALSE;
    }, TRUE);
    
    // Run message loop
    tracker.messageLoop();
    
    tracker.stopTracking();
    std::cout << "Keyboard tracking stopped. Final data saved.\n";
    
    return 0;
}