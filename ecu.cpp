#include <iostream>
#include <thread>
#include <vector>
#include <chrono>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <algorithm>

// Configuration constants
const int DEFAULT_CYLINDERS = 4;
const int MIN_RPM = 800;
const int MAX_RPM = 8000;
const int SAFE_LOOP_THRESHOLD_US = 500;  // Minimum safe loop time in microseconds
const int MAX_QUEUED_CYCLES = 10;        // Maximum allowed queued cycles before failure
const int RPM_RAMP_STEP = 50;            // RPM increase per step during auto-ramp
const int RAMP_INTERVAL_MS = 100;        // Time between RPM increases during auto-ramp

// Thread-safe console output
std::mutex console_mutex;

// Global simulation state
std::atomic<bool> simulation_running{true};
std::atomic<int> current_rpm{MIN_RPM};
std::atomic<int> total_failures{0};
std::atomic<int> total_overloads{0};

// Cylinder thread statistics
struct CylinderStats {
    std::atomic<uint64_t> cycle_count{0};
    std::atomic<uint64_t> missed_cycles{0};
    std::atomic<uint64_t> overload_count{0};
    std::atomic<double> avg_loop_time_us{0.0};
    std::atomic<bool> is_failing{false};
};

// Thread-safe logging function
void log_message(const std::string& message) {
    std::lock_guard<std::mutex> lock(console_mutex);
    auto now = std::chrono::steady_clock::now();
    auto time_since_epoch = now.time_since_epoch();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_since_epoch).count();
    
    std::cout << "[" << std::setw(8) << ms % 100000 << "ms] " << message << std::endl;
}

// Calculate firing interval based on RPM
// For a 4-stroke engine: each cylinder fires once every 2 crankshaft rotations
// Firing interval = (60 seconds / RPM) * 2 / number_of_cylinders
double calculate_firing_interval_ms(int rpm, int num_cylinders) {
    return (60.0 * 1000.0 * 2.0) / (rpm * num_cylinders);
}

// Simulate ECU processing load (sensor reading, calculations, etc.)
void simulate_ecu_processing(int cylinder_id, double processing_time_us) {
    // Simulate variable processing time based on RPM (higher RPM = more processing)
    double load_factor = std::min(2.0, current_rpm.load() / 4000.0);
    auto sleep_time = std::chrono::microseconds(static_cast<int>(processing_time_us * load_factor));
    std::this_thread::sleep_for(sleep_time);
}

// Main cylinder thread function
void cylinder_thread(int cylinder_id, int num_cylinders, CylinderStats& stats) {
    std::queue<std::chrono::steady_clock::time_point> cycle_queue;
    auto last_cycle_time = std::chrono::steady_clock::now();
    double total_loop_time = 0.0;
    uint64_t local_cycle_count = 0;
    
    // Calculate phase offset for this cylinder (evenly distribute firing times)
    double phase_offset_ms = (cylinder_id * calculate_firing_interval_ms(current_rpm.load(), num_cylinders)) / num_cylinders;
    
    log_message("Cylinder " + std::to_string(cylinder_id) + " thread started with phase offset: " + 
                std::to_string(phase_offset_ms) + "ms");
    
    // Initial delay to stagger cylinder firing
    std::this_thread::sleep_for(std::chrono::microseconds(static_cast<int>(phase_offset_ms * 1000)));
    
    while (simulation_running.load()) {
        auto cycle_start = std::chrono::steady_clock::now();
        
        // Calculate current firing interval
        int rpm = current_rpm.load();
        double firing_interval_ms = calculate_firing_interval_ms(rpm, num_cylinders);
        double firing_interval_us = firing_interval_ms * 1000.0;
        
        // Check if we're running too fast (safety threshold)
        if (firing_interval_us < SAFE_LOOP_THRESHOLD_US) {
            stats.is_failing.store(true);
            total_failures.fetch_add(1);
            
            std::stringstream ss;
            ss << "⚠️  CRITICAL: Cylinder " << cylinder_id 
               << " firing interval (" << std::fixed << std::setprecision(1) 
               << firing_interval_us << "μs) below safety threshold (" 
               << SAFE_LOOP_THRESHOLD_US << "μs) at " << rpm << " RPM";
            log_message(ss.str());
        } else {
            stats.is_failing.store(false);
        }
        
        // Check for cycle queue overflow (simulating missed cycles)
        cycle_queue.push(cycle_start);
        if (cycle_queue.size() > MAX_QUEUED_CYCLES) {
            cycle_queue.pop();
            stats.missed_cycles.fetch_add(1);
            
            std::stringstream ss;
            ss << "🔥 OVERLOAD: Cylinder " << cylinder_id 
               << " dropping cycles (queue size: " << cycle_queue.size() 
               << ", missed: " << stats.missed_cycles.load() << ")";
            log_message(ss.str());
            
            stats.overload_count.fetch_add(1);
            total_overloads.fetch_add(1);
        }
        
        // Simulate ignition event
        {
            std::stringstream ss;
            ss << "🔥 Cylinder " << cylinder_id << " IGNITION at " << rpm 
               << " RPM (interval: " << std::fixed << std::setprecision(1) 
               << firing_interval_ms << "ms)";
            log_message(ss.str());
        }
        
        // Simulate ECU processing (sensor reading, fuel injection calc, etc.)
        simulate_ecu_processing(cylinder_id, 50 + (rpm / 100)); // Base 50μs + RPM-dependent load
        
        // Simulate sensor events (knock sensor, O2 sensor, etc.)
        {
            std::stringstream ss;
            ss << "📊 Cylinder " << cylinder_id << " sensors: Knock=OK, Temp=" 
               << (80 + (rpm / 100)) << "°C, Pressure=" << (1.0 + rpm / 10000.0) << "bar";
            log_message(ss.str());
        }
        
        // Update statistics
        auto cycle_end = std::chrono::steady_clock::now();
        auto loop_time = std::chrono::duration_cast<std::chrono::microseconds>(cycle_end - cycle_start);
        total_loop_time += loop_time.count();
        local_cycle_count++;
        
        stats.cycle_count.store(local_cycle_count);
        stats.avg_loop_time_us.store(total_loop_time / local_cycle_count);
        
        // Check for soft lockup (loop taking too long)
        if (loop_time.count() > firing_interval_us * 0.8) {
            std::stringstream ss;
            ss << "⚠️  SOFT LOCKUP: Cylinder " << cylinder_id 
               << " loop took " << loop_time.count() << "μs (limit: " 
               << std::fixed << std::setprecision(1) << firing_interval_us * 0.8 << "μs)";
            log_message(ss.str());
        }
        
        // Wait for next firing cycle
        auto next_cycle = cycle_start + std::chrono::microseconds(static_cast<int>(firing_interval_us));
        auto now = std::chrono::steady_clock::now();
        
        if (next_cycle > now) {
            std::this_thread::sleep_until(next_cycle);
        } else {
            // We're already behind schedule
            stats.missed_cycles.fetch_add(1);
            std::stringstream ss;
            ss << "⏰ TIMING MISS: Cylinder " << cylinder_id 
               << " behind by " << std::chrono::duration_cast<std::chrono::microseconds>(now - next_cycle).count() << "μs";
            log_message(ss.str());
        }
        
        // Clean up old cycles from queue
        while (!cycle_queue.empty() && 
               std::chrono::duration_cast<std::chrono::milliseconds>(cycle_start - cycle_queue.front()).count() > 1000) {
            cycle_queue.pop();
        }
    }
    
    log_message("Cylinder " + std::to_string(cylinder_id) + " thread terminated");
}

// Statistics display thread
void stats_display_thread(const std::vector<CylinderStats>& cylinder_stats, int num_cylinders) {
    while (simulation_running.load()) {
        std::this_thread::sleep_for(std::chrono::seconds(2));
        
        {
            std::lock_guard<std::mutex> lock(console_mutex);
            std::cout << "\n" << std::string(80, '=') << std::endl;
            std::cout << "📊 ECU SIMULATION STATISTICS" << std::endl;
            std::cout << "Current RPM: " << current_rpm.load() << std::endl;
            std::cout << "Total System Failures: " << total_failures.load() << std::endl;
            std::cout << "Total Overloads: " << total_overloads.load() << std::endl;
            std::cout << std::string(80, '-') << std::endl;
            
            for (int i = 0; i < num_cylinders; ++i) {
                const auto& stats = cylinder_stats[i];
                std::cout << "Cylinder " << i << ": "
                         << "Cycles=" << stats.cycle_count.load()
                         << ", Missed=" << stats.missed_cycles.load()
                         << ", Overloads=" << stats.overload_count.load()
                         << ", AvgLoop=" << std::fixed << std::setprecision(1) 
                         << stats.avg_loop_time_us.load() << "μs"
                         << (stats.is_failing.load() ? " [FAILING]" : " [OK]")
                         << std::endl;
            }
            std::cout << std::string(80, '=') << "\n" << std::endl;
        }
    }
}

// RPM control thread (for automatic ramp-up)
void rpm_ramp_thread() {
    log_message("🚀 Starting automatic RPM ramp from " + std::to_string(MIN_RPM) + " to " + std::to_string(MAX_RPM));
    
    int rpm = MIN_RPM;
    while (simulation_running.load() && rpm <= MAX_RPM) {
        current_rpm.store(rpm);
        
        if (rpm % 500 == 0) {  // Log every 500 RPM
            log_message("🏁 RPM ramp: " + std::to_string(rpm));
        }
        
        std::this_thread::sleep_for(std::chrono::milliseconds(RAMP_INTERVAL_MS));
        rpm += RPM_RAMP_STEP;
    }
    
    log_message("🏁 Maximum RPM reached: " + std::to_string(MAX_RPM));
    
    // Hold at max RPM for observation
    std::this_thread::sleep_for(std::chrono::seconds(5));
    
    log_message("🛑 Simulation complete - shutting down");
    simulation_running.store(false);
}

int main(int argc, char* argv[]) {
    int num_cylinders = DEFAULT_CYLINDERS;
    bool auto_ramp = false;
    
    // Parse command line arguments
    if (argc > 1) {
        num_cylinders = std::atoi(argv[1]);
        if (num_cylinders < 1 || num_cylinders > 16) {
            std::cerr << "Error: Number of cylinders must be between 1 and 16" << std::endl;
            return 1;
        }
    }
    
    if (argc > 2 && std::string(argv[2]) == "auto") {
        auto_ramp = true;
    }
    
    std::cout << "🏭 Engine ECU Real-Time Simulator" << std::endl;
    std::cout << "Configuration:" << std::endl;
    std::cout << "  - Cylinders: " << num_cylinders << std::endl;
    std::cout << "  - Safety threshold: " << SAFE_LOOP_THRESHOLD_US << "μs" << std::endl;
    std::cout << "  - Max queued cycles: " << MAX_QUEUED_CYCLES << std::endl;
    std::cout << "  - RPM range: " << MIN_RPM << " - " << MAX_RPM << std::endl;
    
    if (auto_ramp) {
        std::cout << "  - Mode: Automatic RPM ramp" << std::endl;
    } else {
        std::cout << "  - Mode: Manual RPM control" << std::endl;
        std::cout << "Usage: Enter RPM values (800-8000) or 'q' to quit" << std::endl;
    }
    
    std::cout << std::string(80, '=') << std::endl;
    
    // Initialize cylinder statistics
    std::vector<CylinderStats> cylinder_stats(num_cylinders);
    
    // Start cylinder threads
    std::vector<std::thread> cylinder_threads;
    for (int i = 0; i < num_cylinders; ++i) {
        cylinder_threads.emplace_back(cylinder_thread, i, num_cylinders, std::ref(cylinder_stats[i]));
    }
    
    // Start statistics display thread
    std::thread stats_thread(stats_display_thread, std::ref(cylinder_stats), num_cylinders);
    
    // Start RPM control
    std::thread rpm_thread;
    if (auto_ramp) {
        rpm_thread = std::thread(rpm_ramp_thread);
    }
    
    // Manual RPM control or wait for auto-ramp completion
    if (!auto_ramp) {
        std::string input;
        while (simulation_running.load()) {
            std::cout << "Enter RPM (800-8000) or 'q' to quit: ";
            std::getline(std::cin, input);
            
            if (input == "q" || input == "quit") {
                break;
            }
            
            try {
                int rpm = std::stoi(input);
                if (rpm >= MIN_RPM && rpm <= MAX_RPM) {
                    current_rpm.store(rpm);
                    log_message("🎛️  RPM set to: " + std::to_string(rpm));
                } else {
                    std::cout << "RPM must be between " << MIN_RPM << " and " << MAX_RPM << std::endl;
                }
            } catch (const std::exception& e) {
                std::cout << "Invalid input. Please enter a number between " << MIN_RPM << " and " << MAX_RPM << std::endl;
            }
        }
    } else {
        // Wait for auto-ramp to complete
        rpm_thread.join();
    }
    
    // Shutdown sequence
    log_message("🛑 Initiating shutdown sequence...");
    simulation_running.store(false);
    
    // Wait for all threads to complete
    for (auto& thread : cylinder_threads) {
        thread.join();
    }
    
    stats_thread.join();
    
    // Final statistics
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "📊 FINAL SIMULATION RESULTS" << std::endl;
    std::cout << "Total System Failures: " << total_failures.load() << std::endl;
    std::cout << "Total Overloads: " << total_overloads.load() << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    uint64_t total_cycles = 0;
    uint64_t total_missed = 0;
    for (int i = 0; i < num_cylinders; ++i) {
        const auto& stats = cylinder_stats[i];
        total_cycles += stats.cycle_count.load();
        total_missed += stats.missed_cycles.load();
        
        std::cout << "Cylinder " << i << " Final Stats:" << std::endl;
        std::cout << "  Total Cycles: " << stats.cycle_count.load() << std::endl;
        std::cout << "  Missed Cycles: " << stats.missed_cycles.load() << std::endl;
        std::cout << "  Overloads: " << stats.overload_count.load() << std::endl;
        std::cout << "  Avg Loop Time: " << std::fixed << std::setprecision(1) 
                 << stats.avg_loop_time_us.load() << "μs" << std::endl;
        std::cout << "  Success Rate: " << std::fixed << std::setprecision(2)
                 << (100.0 * (stats.cycle_count.load() - stats.missed_cycles.load()) / stats.cycle_count.load()) 
                 << "%" << std::endl;
    }
    
    std::cout << std::string(80, '-') << std::endl;
    std::cout << "Overall Success Rate: " << std::fixed << std::setprecision(2)
             << (100.0 * (total_cycles - total_missed) / total_cycles) << "%" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    
    return 0;
}