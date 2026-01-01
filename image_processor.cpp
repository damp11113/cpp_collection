#include <opencv2/opencv.hpp>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <set>

struct Color {
    int r, g, b;
    Color(int r = 0, int g = 0, int b = 0) : r(r), g(g), b(b) {}
};

struct Segment {
    int start_x, end_x, y;
    Color color;
    Segment(int sx, int ex, int y, Color c) : start_x(sx), end_x(ex), y(y), color(c) {}
};

struct Box {
    int start_x, end_x, start_y, end_y;
    Color color;
    Box(int sx, int ex, int sy, int ey, Color c) : start_x(sx), end_x(ex), start_y(sy), end_y(ey), color(c) {}
};

class ImageBoxProcessor {
private:
    cv::Mat image;
    int height, width;
    double colorThreshold;
    std::string algorithm;

public:
    ImageBoxProcessor(const std::string& imagePath, double threshold = 90.0, const std::string& alg = "manhattan") 
        : colorThreshold(threshold), algorithm(alg) {
        image = cv::imread(imagePath);
        if (image.empty()) {
            throw std::runtime_error("Could not load image: " + imagePath);
        }
        height = image.rows;
        width = image.cols;
    }

    double compareRGBColors(const Color& c1, const Color& c2, const std::string& method = "euclidean") {
        if (method == "euclidean") {
            double distance = sqrt(pow(c1.r - c2.r, 2) + pow(c1.g - c2.g, 2) + pow(c1.b - c2.b, 2));
            double maxDistance = sqrt(3 * pow(255, 2));
            return std::max(0.0, std::min(100.0, (1 - distance / maxDistance) * 100));
        }
        else if (method == "manhattan") {
            double distance = abs(c1.r - c2.r) + abs(c1.g - c2.g) + abs(c1.b - c2.b);
            double maxDistance = 3 * 255;
            return std::max(0.0, std::min(100.0, (1 - distance / maxDistance) * 100));
        }
        else if (method == "cosine") {
            double dotProduct = c1.r * c2.r + c1.g * c2.g + c1.b * c2.b;
            double norm1 = sqrt(c1.r * c1.r + c1.g * c1.g + c1.b * c1.b);
            double norm2 = sqrt(c2.r * c2.r + c2.g * c2.g + c2.b * c2.b);
            
            if (norm1 == 0 || norm2 == 0) {
                return 0.0;
            }
            
            double cosineSim = dotProduct / (norm1 * norm2);
            return std::max(0.0, std::min(100.0, cosineSim * 100));
        }
        return 0.0;
    }

    bool segmentsAreSimilar(const Segment& seg1, const Segment& seg2, double threshold = 80.0) {
        // Check x coordinate tolerance
        int xTolerance = 5;
        if (abs(seg1.start_x - seg2.start_x) > xTolerance || abs(seg1.end_x - seg2.end_x) > xTolerance) {
            return false;
        }

        // Check color similarity
        double similarity = compareRGBColors(seg1.color, seg2.color, "euclidean");
        return similarity >= threshold;
    }

    std::vector<char> packBoxData(const std::vector<Box>& boxes) {
        std::vector<char> packed;
        for (const auto& box : boxes) {
            // Pack as: start_x(4), end_x(4), start_y(4), end_y(4), R(1), G(1), B(1) = 19 bytes per box
            packed.insert(packed.end(), reinterpret_cast<const char*>(&box.start_x), 
                         reinterpret_cast<const char*>(&box.start_x) + sizeof(int));
            packed.insert(packed.end(), reinterpret_cast<const char*>(&box.end_x), 
                         reinterpret_cast<const char*>(&box.end_x) + sizeof(int));
            packed.insert(packed.end(), reinterpret_cast<const char*>(&box.start_y), 
                         reinterpret_cast<const char*>(&box.start_y) + sizeof(int));
            packed.insert(packed.end(), reinterpret_cast<const char*>(&box.end_y), 
                         reinterpret_cast<const char*>(&box.end_y) + sizeof(int));
            
            char r = static_cast<char>(box.color.r);
            char g = static_cast<char>(box.color.g);
            char b = static_cast<char>(box.color.b);
            packed.push_back(r);
            packed.push_back(g);
            packed.push_back(b);
        }
        return packed;
    }

    std::vector<Box> unpackBoxData(const std::vector<char>& packed) {
        std::vector<Box> boxes;
        const int boxSize = 19; // 4*4 + 3*1 bytes per box
        
        for (size_t i = 0; i < packed.size(); i += boxSize) {
            int start_x = *reinterpret_cast<const int*>(&packed[i]);
            int end_x = *reinterpret_cast<const int*>(&packed[i + 4]);
            int start_y = *reinterpret_cast<const int*>(&packed[i + 8]);
            int end_y = *reinterpret_cast<const int*>(&packed[i + 12]);
            
            int r = static_cast<unsigned char>(packed[i + 16]);
            int g = static_cast<unsigned char>(packed[i + 17]);
            int b = static_cast<unsigned char>(packed[i + 18]);
            
            boxes.emplace_back(start_x, end_x, start_y, end_y, Color(r, g, b));
        }
        return boxes;
    }

    std::vector<Segment> createHorizontalSegments() {
        std::vector<Segment> segments;
        
        std::cout << "Phase 1: Creating horizontal segments..." << std::endl;
        
        for (int y = 0; y < height; y++) {
            int startX = 0;
            cv::Vec3b startColor = image.at<cv::Vec3b>(y, 0);
            Color currentStartColor(startColor[2], startColor[1], startColor[0]); // BGR to RGB
            
            for (int x = 1; x < width; x++) {
                cv::Vec3b pixelColor = image.at<cv::Vec3b>(y, x);
                Color currentColor(pixelColor[2], pixelColor[1], pixelColor[0]); // BGR to RGB
                
                double similarity = compareRGBColors(currentStartColor, currentColor, algorithm);
                
                if (similarity < colorThreshold) {
                    segments.emplace_back(startX, x - 1, y, currentStartColor);
                    startX = x;
                    currentStartColor = currentColor;
                }
            }
            
            // Add the last segment of the row
            if (startX < width) {
                segments.emplace_back(startX, width - 1, y, currentStartColor);
            }
            
            // Simple progress indicator
            if (y % (height / 20) == 0) {
                std::cout << "Progress: " << (y * 100) / height << "%" << std::endl;
            }
        }
        
        return segments;
    }

    std::vector<Box> groupSegmentsIntoBoxes(std::vector<Segment>& segments) {
        std::cout << "Phase 2: Grouping segments into boxes..." << std::endl;
        
        // Sort segments by y coordinate, then by x coordinate
        std::sort(segments.begin(), segments.end(), 
                 [](const Segment& a, const Segment& b) {
                     return a.y < b.y || (a.y == b.y && a.start_x < b.start_x);
                 });

        std::vector<Box> boxes;
        std::set<int> usedSegments;

        for (size_t i = 0; i < segments.size(); i++) {
            if (usedSegments.count(i)) continue;

            const Segment& currentSeg = segments[i];
            
            // Start a new box with the current segment
            int boxStartX = currentSeg.start_x;
            int boxEndX = currentSeg.end_x;
            int boxStartY = currentSeg.y;
            int boxEndY = currentSeg.y;
            Color boxColor = currentSeg.color;

            usedSegments.insert(i);

            // Look for vertically adjacent segments that can be merged
            int currentY = currentSeg.y;

            while (true) {
                bool foundAdjacent = false;

                // Look for segments in the next row that match
                for (size_t j = 0; j < segments.size(); j++) {
                    if (usedSegments.count(j) || segments[j].y != currentY + 1) {
                        continue;
                    }

                    // Create a temporary segment for comparison
                    Segment tempSeg(boxStartX, boxEndX, segments[j].y, boxColor);

                    if (segmentsAreSimilar(tempSeg, segments[j], colorThreshold)) {
                        // Extend the box
                        boxEndY = segments[j].y;
                        usedSegments.insert(j);
                        foundAdjacent = true;
                        currentY = segments[j].y;
                        break;
                    }
                }

                if (!foundAdjacent) break;
            }

            // Add the completed box
            boxes.emplace_back(boxStartX, boxEndX, boxStartY, boxEndY, boxColor);
            
            // Progress indicator
            if (i % (segments.size() / 20) == 0) {
                std::cout << "Grouping progress: " << (i * 100) / segments.size() << "%" << std::endl;
            }
        }

        return boxes;
    }

    cv::Mat reconstructImage(const std::vector<Box>& boxes) {
        cv::Mat reconstructed = cv::Mat::zeros(height, width, CV_8UC3);

        for (const auto& box : boxes) {
            cv::Rect rect(box.start_x, box.start_y, 
                         box.end_x - box.start_x + 1, 
                         box.end_y - box.start_y + 1);
            
            cv::Scalar color(box.color.b, box.color.g, box.color.r); // RGB to BGR for OpenCV
            reconstructed(rect) = color;
        }

        return reconstructed;
    }

    void processImage(const std::string& outputPath = "reconstructed_image_boxes.png") {
        std::cout << "Starting image processing..." << std::endl;
        std::cout << "Image dimensions: " << width << "x" << height << std::endl;
        std::cout << "Color threshold: " << colorThreshold << std::endl;
        std::cout << "Algorithm: " << algorithm << std::endl;

        // Phase 1: Create horizontal segments
        std::vector<Segment> segments = createHorizontalSegments();
        std::cout << "Created " << segments.size() << " horizontal segments" << std::endl;

        // Phase 2: Group segments into boxes
        std::vector<Box> boxes = groupSegmentsIntoBoxes(segments);
        std::cout << "Grouped into " << boxes.size() << " boxes" << std::endl;

        // Calculate compression ratio
        double compressionRatio = static_cast<double>(segments.size()) / boxes.size();
        std::cout << "Compression ratio: " << compressionRatio << "x" << std::endl;

        // Pack and unpack data for verification
        std::vector<char> packedData = packBoxData(boxes);
        std::vector<Box> unpackedBoxes = unpackBoxData(packedData);

        // Size comparison
        size_t originalSize = segments.size() * (4 * sizeof(int) + 3 * sizeof(char)); // Approximate
        size_t compressedSize = packedData.size();
        double sizeReduction = (1.0 - static_cast<double>(compressedSize) / originalSize) * 100;
        
        std::cout << "Size comparison:" << std::endl;
        std::cout << "Original segments (estimated): " << originalSize << " bytes" << std::endl;
        std::cout << "Box data packed: " << compressedSize << " bytes" << std::endl;
        std::cout << "Size reduction: " << sizeReduction << "%" << std::endl;

        // Verify packing/unpacking
        bool verificationPassed = (boxes.size() == unpackedBoxes.size());
        std::cout << "Verification passed: " << (verificationPassed ? "Yes" : "No") << std::endl;

        // Reconstruct and save image
        cv::Mat reconstructed = reconstructImage(unpackedBoxes);
        cv::imwrite(outputPath, reconstructed);
        std::cout << "Reconstructed image saved as: " << outputPath << std::endl;

        // Display images
        cv::imshow("Original Image", image);
        cv::imshow("Reconstructed Image", reconstructed);
        
        std::cout << "Press any key to close windows..." << std::endl;
        cv::waitKey(0);
        cv::destroyAllWindows();
    }
};

int main() {
    try {
        // Update this path to your image file
        std::string imagePath = "input_image.png";
        
        ImageBoxProcessor processor(imagePath, 90.0, "manhattan");
        processor.processImage("reconstructed_image_boxes.png");
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}