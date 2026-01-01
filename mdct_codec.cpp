#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <omp.h>  // OpenMP for multi-threading

const double PI = 3.14159265358979323846;

// WAV header structure (same as before)
struct WAVHeader {
    char riff[4];
    int32_t chunkSize;
    char wave[4];
    char fmt[4];
    int32_t subchunk1Size;
    int16_t audioFormat;
    int16_t numChannels;
    int32_t sampleRate;
    int32_t byteRate;
    int16_t blockAlign;
    int16_t bitsPerSample;
    char data[4];
    int32_t dataSize;
};

class MDCTCodec {
private:
    static const int N = 2048; // Window size
    std::vector<double> window;
    std::vector<double> cosTable;
    std::vector<double> sinTable;

    void initWindow() {
        window.resize(N);
        for (int i = 0; i < N; i++) {
            window[i] = sin(PI * (i + 0.5) / N);
        }
    }

    void initTables() {
        cosTable.resize(N);
        sinTable.resize(N);
        for (int i = 0; i < N; i++) {
            cosTable[i] = cos(PI * (i + 0.5) / N);
            sinTable[i] = sin(PI * (i + 0.5) / N);
        }
    }

public:
    MDCTCodec() {
        initWindow();
        initTables();
    }

    std::vector<double> encode(const std::vector<double>& input) {
        std::vector<double> output(N/2);

        // Parallelize the encoding process
        #pragma omp parallel for
        for (int k = 0; k < N/2; k++) {
            double sum = 0.0;
            for (int n = 0; n < N; n++) {
                sum += input[n] * window[n] * cosTable[(n + 0.5 + N/4) * (k + 0.5)];
            }
            output[k] = sum;
        }
        return output;
    }

    std::vector<double> decode(const std::vector<double>& input) {
        std::vector<double> output(N);

        // Parallelize the decoding process
        #pragma omp parallel for
        for (int n = 0; n < N; n++) {
            double sum = 0.0;
            for (int k = 0; k < N/2; k++) {
                sum += input[k] * cosTable[(n + 0.5 + N/4) * (k + 0.5)];
            }
            output[n] = sum * window[n] * 4.0 / N;
        }
        return output;
    }

    int getWindowSize() const {
        return N;
    }
};

class AudioProcessor {
private:
    WAVHeader header;
    std::vector<double> audioData;
    MDCTCodec codec;

    void readWAV(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Cannot open input file");
        }

        // Read WAV header
        file.read(reinterpret_cast<char*>(&header), sizeof(WAVHeader));

        // Read audio data
        std::vector<int16_t> rawData(header.dataSize / 2);
        file.read(reinterpret_cast<char*>(rawData.data()), header.dataSize);

        // Convert to double
        audioData.resize(rawData.size());
        for (size_t i = 0; i < rawData.size(); i++) {
            audioData[i] = rawData[i] / 32768.0;
        }

        file.close();
    }

    void writeWAV(const std::string& filename, const std::vector<double>& data) {
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Cannot open output file");
        }

        // Write WAV header
        file.write(reinterpret_cast<char*>(&header), sizeof(WAVHeader));

        // Convert back to int16_t and write
        std::vector<int16_t> rawData(data.size());
        for (size_t i = 0; i < data.size(); i++) {
            double sample = std::max(-1.0, std::min(1.0, data[i]));
            rawData[i] = static_cast<int16_t>(sample * 32767);
        }
        file.write(reinterpret_cast<char*>(rawData.data()), header.dataSize);

        file.close();
    }

public:
    void processMDCT(const std::string& inputFile, const std::string& mdctFile, 
                 const std::string& outputFile) {
		// Read input WAV
		readWAV(inputFile);

		// Process audio in overlapping windows
		std::vector<std::vector<double>> mdctCoefficients;
		const int hopSize = codec.getWindowSize() / 2;
		const int windowSize = codec.getWindowSize();  // Store window size outside the loop

		// Show progress for encoding
		int totalFrames = (audioData.size() - windowSize) / hopSize + 1;
		int frameCount = 0;

		std::cout << "Encoding progress: " << std::endl;
		#pragma omp parallel for
		for (size_t i = 0; i + windowSize <= audioData.size(); i += hopSize) {
			std::vector<double> window(audioData.begin() + i, audioData.begin() + i + windowSize);
			std::vector<double> encodedFrame = codec.encode(window);
			mdctCoefficients.push_back(encodedFrame);

			// Update progress
			frameCount++;
			if (frameCount % 10 == 0) {
				int progress = (frameCount * 100) / totalFrames;
				std::cout << "\rEncoding: " << progress << "%   ";
				std::flush(std::cout);
			}
		}
		std::cout << std::endl;

		// Save MDCT coefficients
		std::ofstream mdctOut(mdctFile, std::ios::binary);
		for (const auto& frame : mdctCoefficients) {
			mdctOut.write(reinterpret_cast<const char*>(frame.data()), frame.size() * sizeof(double));
		}
		mdctOut.close();

		// Reconstruct audio
		std::vector<double> reconstructedAudio(audioData.size(), 0.0);
		frameCount = 0;
		std::cout << "Decoding progress: " << std::endl;
		for (size_t i = 0; i < mdctCoefficients.size(); i++) {
			std::vector<double> reconstructedWindow = codec.decode(mdctCoefficients[i]);

			// Overlap-add
			for (size_t j = 0; j < windowSize; j++) {
				if (i * hopSize + j < reconstructedAudio.size()) {
					reconstructedAudio[i * hopSize + j] += reconstructedWindow[j];
				}
			}

			// Update progress
			frameCount++;
			if (frameCount % 10 == 0) {
				int progress = (frameCount * 100) / mdctCoefficients.size();
				std::cout << "\rDecoding: " << progress << "%   ";
				std::flush(std::cout);
			}
		}
		std::cout << std::endl;

		// Write output WAV
		writeWAV(outputFile, reconstructedAudio);
	}

};

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " input.wav mdct.bin output.wav" << std::endl;
        return 1;
    }

    try {
        AudioProcessor processor;
        processor.processMDCT(argv[1], argv[2], argv[3]);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}
