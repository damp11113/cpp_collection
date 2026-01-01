#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdint> 

// WAV header structure
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

    // Constructor to create a default 8kHz, 16-bit mono WAV header
    WAVHeader(int32_t numSamples = 0) {
        std::memcpy(riff, "RIFF", 4);
        std::memcpy(wave, "WAVE", 4);
        std::memcpy(fmt, "fmt ", 4);
        std::memcpy(data, "data", 4);
        
        subchunk1Size = 16;
        audioFormat = 1; // PCM
        numChannels = 1; // Mono
        sampleRate = 8000;
        bitsPerSample = 16;
        blockAlign = numChannels * bitsPerSample / 8;
        byteRate = sampleRate * blockAlign;
        
        dataSize = numSamples * blockAlign;
        chunkSize = 36 + dataSize;
    }
};

// LPC-10 frame structure
struct LPC10Frame {
    std::vector<float> coefficients; // 10 LPC coefficients
    float pitch;                     // Pitch period
    float gain;                      // Frame gain
};

class LPC10Codec {
public:
    static constexpr int FRAME_SIZE = 180;       // 22.5ms at 8kHz
    static const int NUM_COEFFICIENTS = 10;  // LPC order

    // Autocorrelation function
    std::vector<float> computeAutocorrelation(const std::vector<float>& frame) {
        std::vector<float> r(NUM_COEFFICIENTS + 1, 0.0f);
        for (int lag = 0; lag <= NUM_COEFFICIENTS; lag++) {
            for (int n = 0; n < FRAME_SIZE - lag; n++) {
                r[lag] += frame[n] * frame[n + lag];
            }
        }
        return r;
    }

    // Levinson-Durbin algorithm for LPC analysis
    std::vector<float> levinsonDurbin(const std::vector<float>& r) {
        std::vector<float> a(NUM_COEFFICIENTS + 1, 0.0f);
        std::vector<float> e(NUM_COEFFICIENTS + 1, 0.0f);
        std::vector<float> k(NUM_COEFFICIENTS + 1, 0.0f);
        
        e[0] = r[0];
        a[0] = 1.0f;

        for (int i = 1; i <= NUM_COEFFICIENTS; i++) {
            float sum = 0.0f;
            for (int j = 1; j < i; j++) {
                sum += a[j] * r[i - j];
            }
            k[i] = (r[i] - sum) / e[i - 1];
            a[i] = k[i];

            for (int j = 1; j < i; j++) {
                float aj = a[j];
                a[j] = aj - k[i] * a[i - j];
            }

            e[i] = (1.0f - k[i] * k[i]) * e[i - 1];
        }

        return std::vector<float>(a.begin() + 1, a.end());
    }

    // Pitch detection using autocorrelation
    float detectPitch(const std::vector<float>& frame) {
        std::vector<float> r = computeAutocorrelation(frame);
        float maxCorr = 0.0f;
        int maxLag = 0;
        
        // Search for pitch in the range 20Hz - 500Hz
        int minLag = FRAME_SIZE / 25;  // 500Hz
        int maxLagSearch = FRAME_SIZE / 2;  // 20Hz
        
        for (int lag = minLag; lag < maxLagSearch; lag++) {
            if (r[lag] > maxCorr) {
                maxCorr = r[lag];
                maxLag = lag;
            }
        }
        
        return maxLag > 0 ? FRAME_SIZE / static_cast<float>(maxLag) : 0.0f;
    }

    // Compute frame gain
    float computeGain(const std::vector<float>& frame) {
        float sum = 0.0f;
        for (const float& sample : frame) {
            sum += sample * sample;
        }
        return std::sqrt(sum / frame.size());
    }

    // Encode a frame of audio
    LPC10Frame encode(const std::vector<float>& frame) {
        LPC10Frame encoded;
        
        // Compute autocorrelation
        std::vector<float> r = computeAutocorrelation(frame);
        
        // Compute LPC coefficients
        encoded.coefficients = levinsonDurbin(r);
        
        // Detect pitch
        encoded.pitch = detectPitch(frame);
        
        // Compute gain
        encoded.gain = computeGain(frame);
        
        return encoded;
    }

    // Synthesize speech from LPC parameters
    std::vector<float> decode(const LPC10Frame& frame) {
        std::vector<float> output(FRAME_SIZE, 0.0f);
        std::vector<float> excitation(FRAME_SIZE, 0.0f);
        
        // Generate excitation signal
        if (frame.pitch > 0) {
            // Voiced speech - use pulse train
            int period = static_cast<int>(FRAME_SIZE / frame.pitch);
            for (int i = 0; i < FRAME_SIZE; i += period) {
                if (i < FRAME_SIZE) {
                    excitation[i] = 1.0f;
                }
            }
        } else {
            // Unvoiced speech - use white noise
            for (int i = 0; i < FRAME_SIZE; i++) {
                excitation[i] = (static_cast<float>(rand()) / RAND_MAX) * 2.0f - 1.0f;
            }
        }
        
        // Apply gain
        for (float& sample : excitation) {
            sample *= frame.gain;
        }
        
        // Apply LPC synthesis filter
        for (int n = 0; n < FRAME_SIZE; n++) {
            output[n] = excitation[n];
            for (size_t i = 0; i < frame.coefficients.size(); i++) {
                if (n > i) {
                    output[n] += frame.coefficients[i] * output[n - i - 1];
                }
            }
        }
        
        return output;
    }
};

class WAVFile {
public:
    static bool readWAV(const std::string& filename, std::vector<float>& samples, WAVHeader& header) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error opening input file" << std::endl;
            return false;
        }

        // Read WAV header
        file.read(reinterpret_cast<char*>(&header), sizeof(WAVHeader));

        // Verify format
        if (header.audioFormat != 1 || header.numChannels != 1 || 
            header.sampleRate != 8000 || header.bitsPerSample != 16) {
            std::cerr << "Unsupported WAV format. Please use 8kHz, 16-bit, mono." << std::endl;
            return false;
        }

        // Read samples
        int numSamples = header.dataSize / sizeof(int16_t);
        std::vector<int16_t> rawSamples(numSamples);
        file.read(reinterpret_cast<char*>(rawSamples.data()), header.dataSize);

        // Convert to float
        samples.resize(numSamples);
        for (int i = 0; i < numSamples; i++) {
            samples[i] = rawSamples[i] / 32768.0f;
        }

        return true;
    }

    static bool writeWAV(const std::string& filename, const std::vector<float>& samples, const WAVHeader& header) {
        std::ofstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error opening output file" << std::endl;
            return false;
        }

        // Write header
        file.write(reinterpret_cast<const char*>(&header), sizeof(WAVHeader));

        // Convert samples to int16_t and write
        std::vector<int16_t> rawSamples(samples.size());
        for (size_t i = 0; i < samples.size(); i++) {
            float clipped = std::max(-1.0f, std::min(1.0f, samples[i]));
            rawSamples[i] = static_cast<int16_t>(clipped * 32767.0f);
        }
        file.write(reinterpret_cast<const char*>(rawSamples.data()), rawSamples.size() * sizeof(int16_t));

        return true;
    }
};

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <mode> <input.wav> <output.wav>" << std::endl;
        std::cerr << "Mode: encode (wav to lpc) or decode (lpc to wav)" << std::endl;
        return 1;
    }

    std::string mode = argv[1];
    std::string inputFile = argv[2];
    std::string outputFile = argv[3];

    LPC10Codec codec;
    std::vector<float> inputSamples;
    WAVHeader header;

    if (mode == "encode") {
        // Read input WAV file
        if (!WAVFile::readWAV(inputFile, inputSamples, header)) {
            return 1;
        }

        // Process frames
        std::vector<LPC10Frame> frames;
        for (size_t i = 0; i < inputSamples.size(); i += LPC10Codec::FRAME_SIZE) {
            std::vector<float> frame(LPC10Codec::FRAME_SIZE);
            size_t remainingSamples = std::min(LPC10Codec::FRAME_SIZE, 
                                             static_cast<int>(inputSamples.size() - i));
            
            // Copy frame samples
            for (size_t j = 0; j < remainingSamples; j++) {
                frame[j] = inputSamples[i + j];
            }
            
            // Zero-pad if necessary
            for (size_t j = remainingSamples; j < LPC10Codec::FRAME_SIZE; j++) {
                frame[j] = 0.0f;
            }
            
            frames.push_back(codec.encode(frame));
        }

        // Write encoded frames to binary file
        std::ofstream lpcFile(outputFile, std::ios::binary);
        for (const auto& frame : frames) {
            lpcFile.write(reinterpret_cast<const char*>(frame.coefficients.data()), 
                         frame.coefficients.size() * sizeof(float));
            lpcFile.write(reinterpret_cast<const char*>(&frame.pitch), sizeof(float));
            lpcFile.write(reinterpret_cast<const char*>(&frame.gain), sizeof(float));
        }
    } else if (mode == "decode") {
        // Read LPC frames
        std::vector<LPC10Frame> frames;
        std::ifstream lpcFile(inputFile, std::ios::binary);
        
        while (true) {
            LPC10Frame frame;
            frame.coefficients.resize(LPC10Codec::NUM_COEFFICIENTS);
            
            lpcFile.read(reinterpret_cast<char*>(frame.coefficients.data()), 
                        frame.coefficients.size() * sizeof(float));
            lpcFile.read(reinterpret_cast<char*>(&frame.pitch), sizeof(float));
            lpcFile.read(reinterpret_cast<char*>(&frame.gain), sizeof(float));
            
            if (!lpcFile) break;
            frames.push_back(frame);
        }

        // Decode frames
        std::vector<float> outputSamples;
        for (const auto& frame : frames) {
            std::vector<float> decodedFrame = codec.decode(frame);
            outputSamples.insert(outputSamples.end(), decodedFrame.begin(), decodedFrame.end());
        }

        // Create new WAV header for output file
        WAVHeader outputHeader(outputSamples.size());
        
        // Write output WAV file
        if (!WAVFile::writeWAV(outputFile, outputSamples, outputHeader)) {
            return 1;
        }
    } else {
        std::cerr << "Invalid mode. Use 'encode' or 'decode'" << std::endl;
        return 1;
    }

    return 0;
}