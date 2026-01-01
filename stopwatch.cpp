#include <iostream>
#include <chrono>
#include <thread>
#include <atomic>
#include <functional>
#include <iomanip>
#include <sstream>

class MultiIntervalStopwatch {
private:
    std::atomic<bool> running{false};
    std::atomic<long long> milliseconds{0};
    std::chrono::steady_clock::time_point start_time;
    std::thread timer_thread;
    
    // Counter for each interval
    long long ms1_counter = 0;
    long long ms10_counter = 0;
    long long ms100_counter = 0;
    long long s1_counter = 0;
    long long s10_counter = 0;
    long long m1_counter = 0;
    long long m10_counter = 0;
    long long h1_counter = 0;
    long long h10_counter = 0;
    long long d_counter = 0;
    
public:
    MultiIntervalStopwatch() = default;
    
    ~MultiIntervalStopwatch() {
        stop();
    }
    
    // Function callbacks - can be customized
    std::function<void(long long)> ms1 = [](long long count) {
        // Trigger every 1ms
        // std::cout << "1ms trigger - Count: " << count << std::endl;
    };
    
    std::function<void(long long)> ms10 = [](long long count) {
        std::cout << "10ms trigger - Count: " << count << " - Time: " << getCurrentTimeString() << std::endl;
    };
    
    std::function<void(long long)> ms100 = [](long long count) {
        std::cout << "100ms trigger - Count: " << count << " - Time: " << getCurrentTimeString() << std::endl;
    };
    
    std::function<void(long long)> s1 = [](long long count) {
        std::cout << "1s trigger - Count: " << count << " - Time: " << getCurrentTimeString() << std::endl;
    };
    
    std::function<void(long long)> s10 = [](long long count) {
        std::cout << "10s trigger - Count: " << count << " - Time: " << getCurrentTimeString() << std::endl;
    };
    
    std::function<void(long long)> m1 = [](long long count) {
        std::cout << "1min trigger - Count: " << count << " - Time: " << getCurrentTimeString() << std::endl;
    };
    
    std::function<void(long long)> m10 = [](long long count) {
        std::cout << "10min trigger - Count: " << count << " - Time: " << getCurrentTimeString() << std::endl;
    };
    
    std::function<void(long long)> h1 = [](long long count) {
        std::cout << "1hour trigger - Count: " << count << " - Time: " << getCurrentTimeString() << std::endl;
    };
    
    std::function<void(long long)> h10 = [](long long count) {
        std::cout << "10hour trigger - Count: " << count << " - Time: " << getCurrentTimeString() << std::endl;
    };
    
    std::function<void(long long)> d = [](long long count) {
        std::cout << "1day trigger - Count: " << count << " - Time: " << getCurrentTimeString() << std::endl;
    };
    
    void start() {
        if (running.load()) {
            return; // Already running
        }
        
        running = true;
        start_time = std::chrono::steady_clock::now();
        
        timer_thread = std::thread([this]() {
            auto next_tick = std::chrono::steady_clock::now();
            
            while (running.load()) {
                next_tick += std::chrono::milliseconds(1);
                std::this_thread::sleep_until(next_tick);
                
                if (!running.load()) break;
                
                milliseconds++;
                long long current_ms = milliseconds.load();
                
                // Trigger functions at their respective intervals
                ms1_counter++;
                ms1(ms1_counter);
                
                if (current_ms % 10 == 0) {
                    ms10_counter++;
                    ms10(ms10_counter);
                }
                
                if (current_ms % 100 == 0) {
                    ms100_counter++;
                    ms100(ms100_counter);
                }
                
                if (current_ms % 1000 == 0) {
                    s1_counter++;
                    s1(s1_counter);
                }
                
                if (current_ms % 10000 == 0) {
                    s10_counter++;
                    s10(s10_counter);
                }
                
                if (current_ms % 60000 == 0) {
                    m1_counter++;
                    m1(m1_counter);
                }
                
                if (current_ms % 600000 == 0) {
                    m10_counter++;
                    m10(m10_counter);
                }
                
                if (current_ms % 3600000 == 0) {
                    h1_counter++;
                    h1(h1_counter);
                }
                
                if (current_ms % 36000000 == 0) {
                    h10_counter++;
                    h10(h10_counter);
                }
                
                if (current_ms % 86400000 == 0) {
                    d_counter++;
                    d(d_counter);
                }
            }
        });
    }
    
    void stop() {
        if (!running.load()) {
            return; // Already stopped
        }
        
        running = false;
        if (timer_thread.joinable()) {
            timer_thread.join();
        }
    }
    
    void reset() {
        bool was_running = running.load();
        stop();
        
        milliseconds = 0;
        ms1_counter = 0;
        ms10_counter = 0;
        ms100_counter = 0;
        s1_counter = 0;
        s10_counter = 0;
        m1_counter = 0;
        m10_counter = 0;
        h1_counter = 0;
        h10_counter = 0;
        d_counter = 0;
        
        if (was_running) {
            start();
        }
    }
    
    // Get current time - triggers every 1ms access
    long long getCurrentTime() const {
        return milliseconds.load();
    }
    
    // Get formatted current time string
    std::string getCurrentTimeFormatted() const {
        long long total_ms = milliseconds.load();
        
        long long days = total_ms / 86400000;
        total_ms %= 86400000;
        
        long long hours = total_ms / 3600000;
        total_ms %= 3600000;
        
        long long minutes = total_ms / 60000;
        total_ms %= 60000;
        
        long long seconds = total_ms / 1000;
        long long ms = total_ms % 1000;
        
        std::ostringstream oss;
        if (days > 0) {
            oss << days << "d ";
        }
        oss << std::setfill('0') << std::setw(2) << hours << ":"
            << std::setfill('0') << std::setw(2) << minutes << ":"
            << std::setfill('0') << std::setw(2) << seconds << "."
            << std::setfill('0') << std::setw(3) << ms;
        
        return oss.str();
    }
    
    bool isRunning() const {
        return running.load();
    }
    
    // Static function to get current system time as string
    static std::string getCurrentTimeString() {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()) % 1000;
        
        std::ostringstream oss;
        oss << std::put_time(std::localtime(&time_t), "%H:%M:%S");
        oss << "." << std::setfill('0') << std::setw(3) << ms.count();
        return oss.str();
    }
};

// Example usage and test
int main() {
    MultiIntervalStopwatch stopwatch;
    
    // Customize callback functions if needed
    stopwatch.ms1 = [](long long count) {
        std::cout << "Custom 10ms - Count: " << count << std::endl;
    };
    
    stopwatch.s1 = [](long long count) {
        std::cout << "Second passed! Count: " << count 
                  << " - Elapsed: " << MultiIntervalStopwatch::getCurrentTimeString() << std::endl;
    };
    
    std::cout << "Starting stopwatch..." << std::endl;
    stopwatch.start();
    
    // Let it run for a demo period
    std::this_thread::sleep_for(std::chrono::seconds(5));
    
    std::cout << "\nCurrent time: " << stopwatch.getCurrentTime() << "ms" << std::endl;
    std::cout << "Formatted time: " << stopwatch.getCurrentTimeFormatted() << std::endl;
    
    std::cout << "\nStopping stopwatch..." << std::endl;
    stopwatch.stop();
    
    std::cout << "Final time: " << stopwatch.getCurrentTime() << "ms" << std::endl;
    
    return 0;
}