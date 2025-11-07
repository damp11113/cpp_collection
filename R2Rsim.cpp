#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <cstdint>
#include <cstdlib>

// WAV file header structure
struct WAVHeader {
    char riff[4] = {'R', 'I', 'F', 'F'};
    uint32_t fileSize;
    char wave[4] = {'W', 'A', 'V', 'E'};
    char fmt[4] = {'f', 'm', 't', ' '};
    uint32_t fmtSize = 16;
    uint16_t audioFormat = 1; // PCM
    uint16_t numChannels;
    uint32_t sampleRate;
    uint32_t byteRate;
    uint16_t blockAlign;
    uint16_t bitsPerSample;
    char data[4] = {'d', 'a', 't', 'a'};
    uint32_t dataSize;
};

class R2RDAC {
private:
    int bitDepth;
    double r1; // R value (lower resistor)
    double r2; // 2R value (upper resistor)
    double tolerance; // Resistor tolerance (e.g., 0.01 for 1%)
    
    // Convert sample to digital code
    uint32_t sampleToDigital(double sample, int inputBits) {
        // Normalize sample from [-1.0, 1.0] to [0, maxValue]
        double normalized = (sample + 1.0) / 2.0;
        uint32_t maxValue = (1 << inputBits) - 1;
        uint32_t digital = static_cast<uint32_t>(normalized * maxValue + 0.5);
        return std::min(digital, maxValue);
    }
    
    // Simulate R2R DAC ladder network
    double digitalToAnalog(uint32_t digitalValue) {
        double vref = 1.0; // Reference voltage
        double output = 0.0;
        
        // Apply resistor tolerance variations
        std::vector<double> r1Values(bitDepth);
        std::vector<double> r2Values(bitDepth);
        
        for (int i = 0; i < bitDepth; i++) {
            // Add tolerance variation (simplified - same for all for now)
            r1Values[i] = r1 * (1.0 + tolerance * (rand() / (double)RAND_MAX - 0.5) * 2.0);
            r2Values[i] = r2 * (1.0 + tolerance * (rand() / (double)RAND_MAX - 0.5) * 2.0);
        }
        
        // Simulate R2R ladder from MSB to LSB
        for (int bit = bitDepth - 1; bit >= 0; bit--) {
            bool bitValue = (digitalValue >> bit) & 1;
            
            // Calculate contribution of this bit using voltage division
            // In R2R ladder, each bit contributes Vref * (1/2)^(bitDepth-bit)
            double bitWeight = vref / pow(2.0, bitDepth - bit);
            
            // Apply non-idealities from resistor tolerance
            double actualWeight = bitWeight * (r1Values[bit] / r1);
            
            if (bitValue) {
                output += actualWeight;
            }
        }
        
        return output;
    }
    
public:
    R2RDAC(int bits = 16, double r1Val = 10000.0, double r2Val = 20000.0, double tol = 0.01)
        : bitDepth(bits), r1(r1Val), r2(r2Val), tolerance(tol) {
        if (bitDepth < 1 || bitDepth > 32) {
            throw std::runtime_error("Bit depth must be between 1 and 32");
        }
    }
    
    // Process a single sample through the DAC
    double process(double sample, int inputBits) {
        // Convert to digital
        uint32_t digital = sampleToDigital(sample, inputBits);
        
        // Quantize to DAC bit depth if needed
        if (bitDepth < inputBits) {
            digital >>= (inputBits - bitDepth);
        }
        
        // Convert back to analog through R2R ladder
        double analog = digitalToAnalog(digital);
        
        // Normalize output back to [-1.0, 1.0]
        double normalized = (analog / 1.0) * 2.0 - 1.0;
        
        // Clamp output
        return std::max(-1.0, std::min(1.0, normalized));
    }
    
    void setBitDepth(int bits) { bitDepth = bits; }
    void setR1(double val) { r1 = val; }
    void setR2(double val) { r2 = val; }
    void setTolerance(double tol) { tolerance = tol; }
};

class AudioProcessor {
private:
    R2RDAC dac;
    
public:
    AudioProcessor(const R2RDAC& dacConfig) : dac(dacConfig) {}
    
    // Read WAV header
    bool readWAVHeader(std::istream& input, WAVHeader& header) {
        input.read(reinterpret_cast<char*>(&header), sizeof(WAVHeader));
        
        if (std::strncmp(header.riff, "RIFF", 4) != 0 ||
            std::strncmp(header.wave, "WAVE", 4) != 0) {
            return false;
        }
        
        return true;
    }
    
    // Write WAV header
    void writeWAVHeader(std::ostream& output, const WAVHeader& header) {
        output.write(reinterpret_cast<const char*>(&header), sizeof(WAVHeader));
    }
    
    // Process stereo audio stream
    void processStream(std::istream& input, std::ostream& output, size_t bufferSize = 4096) {
        WAVHeader inHeader, outHeader;
        
        if (!readWAVHeader(input, inHeader)) {
            throw std::runtime_error("Invalid WAV file");
        }
        
        // Setup output header
        outHeader = inHeader;
        writeWAVHeader(output, outHeader);
        
        int bytesPerSample = inHeader.bitsPerSample / 8;
        int inputBits = inHeader.bitsPerSample;
        size_t samplesPerBuffer = bufferSize / (bytesPerSample * inHeader.numChannels);
        
        std::vector<char> inputBuffer(bufferSize);
        std::vector<char> outputBuffer(bufferSize);
        
        size_t totalSamples = inHeader.dataSize / (bytesPerSample * inHeader.numChannels);
        size_t processedSamples = 0;
        
        while (processedSamples < totalSamples) {
            size_t samplesToRead = std::min(samplesPerBuffer, totalSamples - processedSamples);
            size_t bytesToRead = samplesToRead * bytesPerSample * inHeader.numChannels;
            
            input.read(inputBuffer.data(), bytesToRead);
            size_t bytesRead = input.gcount();
            
            if (bytesRead == 0) break;
            
            // Process samples
            for (size_t i = 0; i < bytesRead; i += bytesPerSample * inHeader.numChannels) {
                for (int ch = 0; ch < inHeader.numChannels; ch++) {
                    // Read sample based on bit depth
                    double sample = 0.0;
                    int offset = i + ch * bytesPerSample;
                    
                    if (inputBits == 16) {
                        int16_t val = *reinterpret_cast<int16_t*>(&inputBuffer[offset]);
                        sample = val / 32768.0;
                    } else if (inputBits == 24) {
                        int32_t val = 0;
                        std::memcpy(&val, &inputBuffer[offset], 3);
                        if (val & 0x800000) val |= 0xFF000000; // Sign extend
                        sample = val / 8388608.0;
                    } else if (inputBits == 32) {
                        int32_t val = *reinterpret_cast<int32_t*>(&inputBuffer[offset]);
                        sample = val / 2147483648.0;
                    }
                    
                    // Process through R2R DAC
                    double processed = dac.process(sample, inputBits);
                    
                    // Write output sample
                    if (inputBits == 16) {
                        int16_t outVal = static_cast<int16_t>(processed * 32767.0);
                        std::memcpy(&outputBuffer[offset], &outVal, sizeof(int16_t));
                    } else if (inputBits == 24) {
                        int32_t outVal = static_cast<int32_t>(processed * 8388607.0);
                        std::memcpy(&outputBuffer[offset], &outVal, 3);
                    } else if (inputBits == 32) {
                        int32_t outVal = static_cast<int32_t>(processed * 2147483647.0);
                        std::memcpy(&outputBuffer[offset], &outVal, sizeof(int32_t));
                    }
                }
            }
            
            output.write(outputBuffer.data(), bytesRead);
            processedSamples += samplesToRead;
        }
    }
    
    R2RDAC& getDAC() { return dac; }
};

void printUsage(const char* progName) {
    std::cout << "R2R DAC Simulator - Stereo WAV Processor\n\n";
    std::cout << "Usage: " << progName << " [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -i <file>     Input WAV file (or - for stdin)\n";
    std::cout << "  -o <file>     Output WAV file (or - for stdout)\n";
    std::cout << "  -b <bits>     DAC bit depth (1-32, default: 16)\n";
    std::cout << "  -r1 <ohms>    R1 resistor value (default: 10000)\n";
    std::cout << "  -r2 <ohms>    R2 resistor value (default: 20000)\n";
    std::cout << "  -t <percent>  Resistor tolerance (default: 0.01 = 1%)\n";
    std::cout << "  -h            Show this help\n\n";
    std::cout << "Examples:\n";
    std::cout << "  " << progName << " -i input.wav -o output.wav -b 12\n";
    std::cout << "  ffmpeg -i input.mp3 -f wav - | " << progName << " -i - -o - -b 8 | ffplay -\n";
}

int main(int argc, char* argv[]) {
    std::string inputFile, outputFile;
    int bitDepth = 16;
    double r1 = 10000.0;
    double r2 = 20000.0;
    double tolerance = 0.01;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-i" && i + 1 < argc) {
            inputFile = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            outputFile = argv[++i];
        } else if (arg == "-b" && i + 1 < argc) {
            bitDepth = std::atoi(argv[++i]);
        } else if (arg == "-r1" && i + 1 < argc) {
            r1 = std::atof(argv[++i]);
        } else if (arg == "-r2" && i + 1 < argc) {
            r2 = std::atof(argv[++i]);
        } else if (arg == "-t" && i + 1 < argc) {
            tolerance = std::atof(argv[++i]);
        }
    }
    
    if (inputFile.empty() || outputFile.empty()) {
        printUsage(argv[0]);
        return 1;
    }
    
    try {
        // Setup DAC with specified parameters
        R2RDAC dac(bitDepth, r1, r2, tolerance);
        AudioProcessor processor(dac);
        
        // Setup input stream
        std::istream* input;
        std::ifstream inputFileStream;
        if (inputFile == "-") {
            input = &std::cin;
            std::cin.sync_with_stdio(false);
        } else {
            inputFileStream.open(inputFile, std::ios::binary);
            if (!inputFileStream) {
                std::cerr << "Error: Cannot open input file: " << inputFile << std::endl;
                return 1;
            }
            input = &inputFileStream;
        }
        
        // Setup output stream
        std::ostream* output;
        std::ofstream outputFileStream;
        if (outputFile == "-") {
            output = &std::cout;
            std::cout.sync_with_stdio(false);
        } else {
            outputFileStream.open(outputFile, std::ios::binary);
            if (!outputFileStream) {
                std::cerr << "Error: Cannot open output file: " << outputFile << std::endl;
                return 1;
            }
            output = &outputFileStream;
        }
        
        // Process audio
        std::cerr << "Processing with R2R DAC: " << bitDepth << " bits, "
                  << "R1=" << r1 << "Ω, R2=" << r2 << "Ω, "
                  << "Tolerance=" << (tolerance * 100) << "%\n";
        
        processor.processStream(*input, *output);
        
        std::cerr << "Processing complete!\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}