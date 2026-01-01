#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>

using namespace cv;
using namespace std;

int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " <input_video_path>" << endl;
        return -1;
    }
    
    string inputPath = argv[1];
    
    // Open input video
    VideoCapture cap(inputPath);
    if (!cap.isOpened()) {
        cerr << "Error: Cannot open video file " << inputPath << endl;
        return -1;
    }
    
    // Get video properties
    int frameWidth = static_cast<int>(cap.get(CAP_PROP_FRAME_WIDTH));
    int frameHeight = static_cast<int>(cap.get(CAP_PROP_FRAME_HEIGHT));
    double fps = cap.get(CAP_PROP_FPS);
    int totalFrames = static_cast<int>(cap.get(CAP_PROP_FRAME_COUNT));
    
    cout << "Video Properties:" << endl;
    cout << "Resolution: " << frameWidth << "x" << frameHeight << endl;
    cout << "FPS: " << fps << endl;
    cout << "Total frames: " << totalFrames << endl;
    
    // Define codec (you can change this based on your preference)
    int fourcc = VideoWriter::fourcc('M', 'P', '4', 'V');
    
    // Create output video writers for each channel
    VideoWriter redWriter("red_channel.mp4", fourcc, fps, Size(frameWidth, frameHeight), true);
    VideoWriter greenWriter("green_channel.mp4", fourcc, fps, Size(frameWidth, frameHeight), true);
    VideoWriter blueWriter("blue_channel.mp4", fourcc, fps, Size(frameWidth, frameHeight), true);
    
    // Check if writers are initialized
    if (!redWriter.isOpened() || !greenWriter.isOpened() || !blueWriter.isOpened()) {
        cerr << "Error: Cannot create output video files" << endl;
        return -1;
    }
    
    cout << "Processing video..." << endl;
    
    Mat frame, bgrChannels[3];
    Mat redFrame, greenFrame, blueFrame;
    Mat zeros;
    
    int frameCount = 0;
    
    while (true) {
        // Read frame
        cap >> frame;
        if (frame.empty()) {
            break;
        }
        
        // Create a zero matrix for unused channels
        zeros = Mat::zeros(frame.size(), CV_8UC1);
        
        // Split BGR channels (OpenCV uses BGR, not RGB)
        split(frame, bgrChannels);
        
        // Create channel-specific frames
        // For red channel (B=0, G=0, R=original)
        vector<Mat> redChannelMats = {zeros, zeros, bgrChannels[2]};
        merge(redChannelMats, redFrame);
        
        // For green channel (B=0, G=original, R=0)
        vector<Mat> greenChannelMats = {zeros, bgrChannels[1], zeros};
        merge(greenChannelMats, greenFrame);
        
        // For blue channel (B=original, G=0, R=0)
        vector<Mat> blueChannelMats = {bgrChannels[0], zeros, zeros};
        merge(blueChannelMats, blueFrame);
        
        // Write frames to respective output videos
        redWriter.write(redFrame);
        greenWriter.write(greenFrame);
        blueWriter.write(blueFrame);
        
        frameCount++;
        
        // Show progress
        if (frameCount % 30 == 0) {
            cout << "Processed " << frameCount << "/" << totalFrames << " frames" << endl;
        }
        
        // Optional: Display frames (comment out for faster processing)
        /*
        imshow("Original", frame);
        imshow("Red Channel", redFrame);
        imshow("Green Channel", greenFrame);
        imshow("Blue Channel", blueFrame);
        
        if (waitKey(1) == 'q') {
            break;
        }
        */
    }
    
    // Release resources
    cap.release();
    redWriter.release();
    greenWriter.release();
    blueWriter.release();
    destroyAllWindows();
    
    cout << "Processing complete!" << endl;
    cout << "Output files created:" << endl;
    cout << "- red_channel.mp4" << endl;
    cout << "- green_channel.mp4" << endl;
    cout << "- blue_channel.mp4" << endl;
    
    return 0;
}

// Alternative version: Extract to grayscale channel videos
/*
void extractToGrayscaleChannels(const string& inputPath) {
    VideoCapture cap(inputPath);
    if (!cap.isOpened()) {
        cerr << "Error: Cannot open video file" << endl;
        return;
    }
    
    int frameWidth = static_cast<int>(cap.get(CAP_PROP_FRAME_WIDTH));
    int frameHeight = static_cast<int>(cap.get(CAP_PROP_FRAME_HEIGHT));
    double fps = cap.get(CAP_PROP_FPS);
    int fourcc = VideoWriter::fourcc('M', 'P', '4', 'V');
    
    VideoWriter redWriter("red_gray.mp4", fourcc, fps, Size(frameWidth, frameHeight), false);
    VideoWriter greenWriter("green_gray.mp4", fourcc, fps, Size(frameWidth, frameHeight), false);
    VideoWriter blueWriter("blue_gray.mp4", fourcc, fps, Size(frameWidth, frameHeight), false);
    
    Mat frame, bgrChannels[3];
    
    while (true) {
        cap >> frame;
        if (frame.empty()) break;
        
        split(frame, bgrChannels);
        
        redWriter.write(bgrChannels[2]);    // Red channel
        greenWriter.write(bgrChannels[1]);  // Green channel
        blueWriter.write(bgrChannels[0]);   // Blue channel
    }
    
    cap.release();
    redWriter.release();
    greenWriter.release();
    blueWriter.release();
}
*/