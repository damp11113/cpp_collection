#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>
#include <filesystem>

int main(int argc, char** argv) {
    // === Argument Check ===
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_video_path>" << std::endl;
        return -1;
    }

    std::string input_path = argv[1];

    // === Open Input Video ===
    cv::VideoCapture cap(input_path);
    if (!cap.isOpened()) {
        std::cerr << "Error: Could not open video file " << input_path << std::endl;
        return -1;
    }

    // === Get Video Info ===
    double fps = cap.get(cv::CAP_PROP_FPS);
    int width = static_cast<int>(cap.get(cv::CAP_PROP_FRAME_WIDTH));
    int height = static_cast<int>(cap.get(cv::CAP_PROP_FRAME_HEIGHT));
    int frame_count = static_cast<int>(cap.get(cv::CAP_PROP_FRAME_COUNT));

    std::cout << "Video info: " << width << "x" << height << " @ " << fps << " fps, "
              << frame_count << " frames" << std::endl;

    // === Generate Output Paths Based on Input Filename ===
    std::filesystem::path input_file(input_path);
    std::string stem = input_file.stem().string();  // filename without extension
    std::string gray_path = stem + "_gray.avi";
    std::string color_path = stem + "_color.avi";

    // === Define Codec and Writers ===
    int fourcc = cv::VideoWriter::fourcc('M', 'J', 'P', 'G');
    cv::VideoWriter out_gray(gray_path, fourcc, fps, cv::Size(width, height), false);
    cv::VideoWriter out_color(color_path, fourcc, fps, cv::Size(width, height), true);

    if (!out_gray.isOpened() || !out_color.isOpened()) {
        std::cerr << "Error: Could not open output video writers" << std::endl;
        return -1;
    }

    // === Process Frames ===
    cv::Mat frame, gray, hsv, pure_color;
    int processed_frames = 0;

    std::cout << "Processing video..." << std::endl;

    while (cap.read(frame)) {
        // Grayscale for LCD1
        cv::cvtColor(frame, gray, cv::COLOR_BGR2GRAY);

        // Color Boost for LCD2
        cv::cvtColor(frame, hsv, cv::COLOR_BGR2HSV);
        std::vector<cv::Mat> hsv_channels;
        cv::split(hsv, hsv_channels);
        hsv_channels[2] = cv::Scalar(255); // Max V channel
        cv::merge(hsv_channels, hsv);
        cv::cvtColor(hsv, pure_color, cv::COLOR_HSV2BGR);

        // Write both
        out_gray.write(gray);
        out_color.write(pure_color);

        processed_frames++;
        if (processed_frames % 100 == 0 || processed_frames == frame_count) {
            std::cout << "Progress: " << processed_frames << "/" << frame_count 
                      << " (" << (100.0 * processed_frames / frame_count) << "%)" << std::endl;
        }
    }

    // === Cleanup ===
    cap.release();
    out_gray.release();
    out_color.release();

    std::cout << "✅ Done: Output saved as '" << gray_path << "' and '" << color_path << "'" << std::endl;
    std::cout << "Processed " << processed_frames << " frames total." << std::endl;

    return 0;
}
