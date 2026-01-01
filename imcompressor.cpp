#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <opencv2/opencv.hpp>

#ifdef _WIN32
    #include <io.h>
    #include <fcntl.h>
#else
    #include <unistd.h>
#endif

struct ColorSegment {
    uint32_t start_x;
    uint32_t end_x;
    uint32_t y;
    uint8_t r;
    uint8_t g;
    uint8_t b;
};

class ImageProcessor {
private:
    std::string detect_algorithm;
    
public:
    ImageProcessor(const std::string& algorithm = "euclidean") : detect_algorithm(algorithm) {}
    
    double compareRgbColors(const cv::Vec3b& color1, const cv::Vec3b& color2) {
        double similarity = 0.0;
        
        if (detect_algorithm == "euclidean") {
            double distance = std::sqrt(
                std::pow(color1[0] - color2[0], 2) + 
                std::pow(color1[1] - color2[1], 2) + 
                std::pow(color1[2] - color2[2], 2)
            );
            double max_distance = std::sqrt(3 * 255 * 255);
            similarity = (1.0 - distance / max_distance) * 100.0;
        }
        else if (detect_algorithm == "manhattan") {
            double distance = std::abs(color1[0] - color2[0]) + 
                             std::abs(color1[1] - color2[1]) + 
                             std::abs(color1[2] - color2[2]);
            double max_distance = 3 * 255;
            similarity = (1.0 - distance / max_distance) * 100.0;
        }
        else if (detect_algorithm == "cosine") {
            double dot_product = color1[0] * color2[0] + color1[1] * color2[1] + color1[2] * color2[2];
            double norm1 = std::sqrt(color1[0] * color1[0] + color1[1] * color1[1] + color1[2] * color1[2]);
            double norm2 = std::sqrt(color2[0] * color2[0] + color2[1] * color2[1] + color2[2] * color2[2]);
            
            if (norm1 == 0 || norm2 == 0) {
                similarity = 0.0;
            } else {
                double cosine_sim = dot_product / (norm1 * norm2);
                similarity = cosine_sim * 100.0;
            }
        }
        
        return std::max(0.0, std::min(100.0, similarity));
    }
    
    std::vector<uint8_t> packDataToBytes(const std::vector<ColorSegment>& segments) {
        std::vector<uint8_t> packed_data;
        packed_data.reserve(segments.size() * 15); // 15 bytes per segment
        
        for (const auto& segment : segments) {
            // Pack as little-endian (matching Python struct format)
            // start_x (4 bytes)
            packed_data.push_back(segment.start_x & 0xFF);
            packed_data.push_back((segment.start_x >> 8) & 0xFF);
            packed_data.push_back((segment.start_x >> 16) & 0xFF);
            packed_data.push_back((segment.start_x >> 24) & 0xFF);
            
            // end_x (4 bytes)
            packed_data.push_back(segment.end_x & 0xFF);
            packed_data.push_back((segment.end_x >> 8) & 0xFF);
            packed_data.push_back((segment.end_x >> 16) & 0xFF);
            packed_data.push_back((segment.end_x >> 24) & 0xFF);
            
            // y (4 bytes)
            packed_data.push_back(segment.y & 0xFF);
            packed_data.push_back((segment.y >> 8) & 0xFF);
            packed_data.push_back((segment.y >> 16) & 0xFF);
            packed_data.push_back((segment.y >> 24) & 0xFF);
            
            // RGB (3 bytes)
            packed_data.push_back(segment.r);
            packed_data.push_back(segment.g);
            packed_data.push_back(segment.b);
        }
        
        return packed_data;
    }
    
    std::vector<ColorSegment> unpackDataFromBytes(const std::vector<uint8_t>& packed_data) {
        std::vector<ColorSegment> segments;
        const size_t segment_size = 15; // 4+4+4+1+1+1 bytes
        
        for (size_t i = 0; i < packed_data.size(); i += segment_size) {
            if (i + segment_size > packed_data.size()) break;
            
            ColorSegment segment;
            
            // Unpack little-endian
            segment.start_x = packed_data[i] | 
                             (packed_data[i+1] << 8) | 
                             (packed_data[i+2] << 16) | 
                             (packed_data[i+3] << 24);
            
            segment.end_x = packed_data[i+4] | 
                           (packed_data[i+5] << 8) | 
                           (packed_data[i+6] << 16) | 
                           (packed_data[i+7] << 24);
            
            segment.y = packed_data[i+8] | 
                       (packed_data[i+9] << 8) | 
                       (packed_data[i+10] << 16) | 
                       (packed_data[i+11] << 24);
            
            segment.r = packed_data[i+12];
            segment.g = packed_data[i+13];
            segment.b = packed_data[i+14];
            
            segments.push_back(segment);
        }
        
        return segments;
    }
    
    std::vector<ColorSegment> processImage(const cv::Mat& image, double similarity_threshold = 80.0) {
        std::vector<ColorSegment> segments;
        int height = image.rows;
        int width = image.cols;
        
        for (int y = 0; y < height; y++) {
            int start_x = 0;
            cv::Vec3b start_color = image.at<cv::Vec3b>(y, 0);
            
            for (int x = 1; x < width; x++) {
                cv::Vec3b current_color = image.at<cv::Vec3b>(y, x);
                double similarity = compareRgbColors(start_color, current_color);
                
                if (similarity < similarity_threshold) {
                    ColorSegment segment;
                    segment.start_x = start_x;
                    segment.end_x = x - 1;
                    segment.y = y;
                    segment.r = start_color[2]; // OpenCV uses BGR, convert to RGB
                    segment.g = start_color[1];
                    segment.b = start_color[0];
                    
                    segments.push_back(segment);
                    
                    start_x = x;
                    start_color = current_color;
                }
            }
            
            // Add the last segment of the row
            if (start_x < width) {
                ColorSegment segment;
                segment.start_x = start_x;
                segment.end_x = width - 1;
                segment.y = y;
                segment.r = start_color[2]; // OpenCV uses BGR, convert to RGB
                segment.g = start_color[1];
                segment.b = start_color[0];
                
                segments.push_back(segment);
            }
        }
        
        return segments;
    }
    
    cv::Mat reconstructImage(const std::vector<ColorSegment>& segments, int height, int width) {
        cv::Mat reconstructed = cv::Mat::zeros(height, width, CV_8UC3);
        
        for (const auto& segment : segments) {
            cv::Scalar color(segment.b, segment.g, segment.r); // Convert RGB back to BGR for OpenCV
            cv::line(reconstructed, 
                    cv::Point(segment.start_x, segment.y), 
                    cv::Point(segment.end_x, segment.y), 
                    color, 1);
        }
        
        return reconstructed;
    }
};

// Function to read binary data from stdin
std::vector<uint8_t> readBinaryFromStdin() {
    std::vector<uint8_t> data;
    
#ifdef _WIN32
    _setmode(_fileno(stdin), _O_BINARY);
#endif
    
    // First, read the size of the data (4 bytes)
    uint32_t data_size;
    if (std::cin.read(reinterpret_cast<char*>(&data_size), 4)) {
        data.resize(data_size);
        std::cin.read(reinterpret_cast<char*>(data.data()), data_size);
    }
    
    return data;
}

// Function to write binary data to stdout
void writeBinaryToStdout(const std::vector<uint8_t>& data) {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    
    // First, write the size of the data (4 bytes)
    uint32_t data_size = data.size();
    std::cout.write(reinterpret_cast<const char*>(&data_size), 4);
    
    // Then write the actual data
    std::cout.write(reinterpret_cast<const char*>(data.data()), data.size());
    std::cout.flush();
}

int main(int argc, char* argv[]) {
    try {
        std::string mode = "process"; // Default mode
        std::string algorithm = "euclidean";
        std::string image_path;
        double threshold = 80.0;
        
        // Parse command line arguments
        for (int i = 1; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == "--mode" && i + 1 < argc) {
                mode = argv[++i];
            } else if (arg == "--algorithm" && i + 1 < argc) {
                algorithm = argv[++i];
            } else if (arg == "--image" && i + 1 < argc) {
                image_path = argv[++i];
            } else if (arg == "--threshold" && i + 1 < argc) {
                threshold = std::stod(argv[++i]);
            }
        }
        
        ImageProcessor processor(algorithm);
        
        if (mode == "process") {
            // Mode 1: Process image and return packed segments
            if (image_path.empty()) {
                std::cerr << "Error: Image path required for process mode" << std::endl;
                return 1;
            }
            
            cv::Mat image = cv::imread(image_path);
            if (image.empty()) {
                std::cerr << "Error: Could not load image: " << image_path << std::endl;
                return 1;
            }
            
            auto segments = processor.processImage(image, threshold);
            auto packed_data = processor.packDataToBytes(segments);
            
            writeBinaryToStdout(packed_data);
            
        } else if (mode == "reconstruct") {
            // Mode 2: Read packed data from stdin and reconstruct image
            auto packed_data = readBinaryFromStdin();
            if (packed_data.empty()) {
                std::cerr << "Error: No data received from stdin" << std::endl;
                return 1;
            }
            
            auto segments = processor.unpackDataFromBytes(packed_data);
            
            // For reconstruction, we need image dimensions
            // We'll find the maximum dimensions from the segments
            int max_width = 0, max_height = 0;
            for (const auto& segment : segments) {
                max_width = std::max(max_width, static_cast<int>(segment.end_x) + 1);
                max_height = std::max(max_height, static_cast<int>(segment.y) + 1);
            }
            
            cv::Mat reconstructed = processor.reconstructImage(segments, max_height, max_width);
            
            // Encode image as PNG and send back
            std::vector<uint8_t> encoded_image;
            cv::imencode(".png", reconstructed, encoded_image);
            
            writeBinaryToStdout(encoded_image);
            
        } else if (mode == "roundtrip") {
            // Mode 3: Read packed data, unpack, repack, and send back
            auto packed_data = readBinaryFromStdin();
            auto segments = processor.unpackDataFromBytes(packed_data);
            auto repacked_data = processor.packDataToBytes(segments);
            
            writeBinaryToStdout(repacked_data);
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}