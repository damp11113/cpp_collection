#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <limits>
#include <cstdint>
#include <numeric>
#include <thread>
#include <mutex>
#include <atomic>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Simple WAV file handling
struct WavHeader {
    char riff[4] = {'R', 'I', 'F', 'F'};
    uint32_t fileSize;
    char wave[4] = {'W', 'A', 'V', 'E'};
    char fmt[4] = {'f', 'm', 't', ' '};
    uint32_t fmtSize = 16;
    uint16_t audioFormat = 1;
    uint16_t numChannels;
    uint32_t sampleRate;
    uint32_t byteRate;
    uint16_t blockAlign;
    uint16_t bitsPerSample;
    char data[4] = {'d', 'a', 't', 'a'};
    uint32_t dataSize;
};

class WavFile {
public:
    static bool read(const std::string& filename, std::vector<float>& audio, int& sampleRate) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) return false;

        WavHeader header;
        file.read(reinterpret_cast<char*>(&header), sizeof(WavHeader));

        sampleRate = header.sampleRate;
        int numSamples = header.dataSize / (header.bitsPerSample / 8);

        if (header.bitsPerSample == 16) {
            std::vector<int16_t> samples(numSamples);
            file.read(reinterpret_cast<char*>(samples.data()), header.dataSize);
            
            audio.resize(numSamples / header.numChannels);
            for (size_t i = 0; i < audio.size(); i++) {
                if (header.numChannels == 1) {
                    audio[i] = samples[i] / 32768.0f;
                } else {
                    audio[i] = (samples[i * 2] + samples[i * 2 + 1]) / 65536.0f;
                }
            }
        }
        return true;
    }

    static bool write(const std::string& filename, const std::vector<float>& audio, int sampleRate) {
        std::ofstream file(filename, std::ios::binary);
        if (!file) return false;

        WavHeader header;
        header.numChannels = 1;
        header.sampleRate = sampleRate;
        header.bitsPerSample = 16;
        header.byteRate = sampleRate * header.numChannels * header.bitsPerSample / 8;
        header.blockAlign = header.numChannels * header.bitsPerSample / 8;
        header.dataSize = audio.size() * header.bitsPerSample / 8;
        header.fileSize = 36 + header.dataSize;

        file.write(reinterpret_cast<const char*>(&header), sizeof(WavHeader));

        std::vector<int16_t> samples(audio.size());
        for (size_t i = 0; i < audio.size(); i++) {
            samples[i] = static_cast<int16_t>(std::clamp(audio[i], -1.0f, 1.0f) * 32767.0f);
        }
        file.write(reinterpret_cast<const char*>(samples.data()), samples.size() * sizeof(int16_t));

        return true;
    }
};

// FFT implementation (simple Cooley-Tukey)
class FFT {
public:
    static void fft(std::vector<std::complex<float>>& x) {
        int N = x.size();
        if (N <= 1) return;

        std::vector<std::complex<float>> even(N / 2), odd(N / 2);
        for (int i = 0; i < N / 2; i++) {
            even[i] = x[i * 2];
            odd[i] = x[i * 2 + 1];
        }

        fft(even);
        fft(odd);

        for (int k = 0; k < N / 2; k++) {
            float angle = -2.0f * static_cast<float>(M_PI) * k / N;
            std::complex<float> t = std::polar(1.0f, angle) * odd[k];
            x[k] = even[k] + t;
            x[k + N / 2] = even[k] - t;
        }
    }

    static std::vector<float> hannWindow(int size) {
        std::vector<float> window(size);
        for (int i = 0; i < size; i++) {
            window[i] = 0.5f * (1.0f - std::cos(2.0f * static_cast<float>(M_PI) * i / (size - 1)));
        }
        return window;
    }
};

// Harmonic structures
struct Harmonic {
    float amp;
    float phase;
};

struct HarmonicObject {
    float freq;
    int nHarmonic;
    float mainAmp;
    std::vector<Harmonic> harmonics;
};

// Peak detection
std::vector<int> findPeaks(const std::vector<float>& data, float threshold) {
    std::vector<int> peaks;
    for (size_t i = 1; i < data.size() - 1; i++) {
        if (data[i] > data[i-1] && data[i] > data[i+1] && data[i] > threshold) {
            peaks.push_back(i);
        }
    }
    return peaks;
}

// Harmonic Extractor
class HarmonicExtractor {
private:
    int sampleRate;
    int windowSize;
    int hopSize;
    float minF0Freq;
    float maxF0Freq;
    float peakThreshold;
    int maxHarmonicsPerF0;
    float maxHarmonicFreqOutput;
    int maxHarmonicFreqObject;
    std::vector<float> window;
    float freqResolution;
    int minF0Bin;
    int maxF0Bin;

public:
    HarmonicExtractor(int sr, int ws, int hs, float minF0, float maxF0,
                     float pt, int maxHarm, float maxHarmFreq, int maxHarmObj)
        : sampleRate(sr), windowSize(ws), hopSize(hs), minF0Freq(minF0),
          maxF0Freq(maxF0), peakThreshold(pt), maxHarmonicsPerF0(maxHarm),
          maxHarmonicFreqOutput(maxHarmFreq), maxHarmonicFreqObject(maxHarmObj) {
        
        window = FFT::hannWindow(windowSize);
        freqResolution = static_cast<float>(sampleRate) / windowSize;
        minF0Bin = static_cast<int>(minF0Freq / freqResolution);
        maxF0Bin = std::min(static_cast<int>(maxF0Freq / freqResolution), windowSize / 2);
    }

    std::vector<HarmonicObject> processChunk(const std::vector<float>& chunk) {
        std::vector<float> paddedChunk(windowSize, 0.0f);
        std::copy(chunk.begin(), chunk.end(), paddedChunk.begin());

        // Apply window
        for (int i = 0; i < windowSize; i++) {
            paddedChunk[i] *= window[i];
        }

        // FFT
        std::vector<std::complex<float>> fftInput(windowSize);
        for (int i = 0; i < windowSize; i++) {
            fftInput[i] = std::complex<float>(paddedChunk[i], 0.0f);
        }
        FFT::fft(fftInput);

        // Get magnitude and phase
        std::vector<float> magnitude(windowSize / 2);
        std::vector<float> phase(windowSize / 2);
        float maxMag = 0.0f;

        for (int i = 0; i < windowSize / 2; i++) {
            magnitude[i] = std::abs(fftInput[i]);
            phase[i] = std::arg(fftInput[i]);
            maxMag = std::max(maxMag, magnitude[i]);
        }

        if (maxMag < 1e-6f) return {};

        // Find peaks
        std::vector<int> peaks = findPeaks(magnitude, peakThreshold * maxMag);
        
        // Filter peaks in frequency range
        std::vector<int> validPeaks;
        std::vector<float> peakMags;
        for (int peak : peaks) {
            if (peak >= minF0Bin && peak < maxF0Bin) {
                validPeaks.push_back(peak);
                peakMags.push_back(magnitude[peak]);
            }
        }

        if (validPeaks.empty()) return {};

        // Sort by magnitude and keep top N
        std::vector<size_t> indices(validPeaks.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), 
                 [&peakMags](size_t a, size_t b) { return peakMags[a] > peakMags[b]; });

        int numPeaks = std::min(static_cast<int>(validPeaks.size()), maxHarmonicFreqObject);
        
        // Extract harmonics
        std::vector<HarmonicObject> harmonicObjects;
        for (int i = 0; i < numPeaks; i++) {
            int peakBin = validPeaks[indices[i]];
            float peakMag = magnitude[peakBin];
            float p0Freq = peakBin * freqResolution;
            float p0Phase = phase[peakBin];

            HarmonicObject obj;
            obj.freq = p0Freq;
            obj.mainAmp = peakMag;
            obj.harmonics.push_back({peakMag, p0Phase});

            // Extract harmonics
            for (int n = 2; n <= maxHarmonicsPerF0; n++) {
                float targetFreq = n * p0Freq;
                if (targetFreq > maxHarmonicFreqOutput || targetFreq >= sampleRate / 2.0f) break;

                int harmBin = static_cast<int>(targetFreq / freqResolution);
                if (harmBin < windowSize / 2) {
                    obj.harmonics.push_back({magnitude[harmBin], phase[harmBin]});
                }
            }

            obj.nHarmonic = obj.harmonics.size();
            harmonicObjects.push_back(obj);
        }

        return harmonicObjects;
    }
};

// Harmonic Packer
class HarmonicPacker {
private:
    int sampleRate;
    int windowSize;
    float binSize;
    float minAmpLinear;

public:
    HarmonicPacker(int sr, int ws) 
        : sampleRate(sr), windowSize(ws) {
        binSize = static_cast<float>(sampleRate) / windowSize;
        minAmpLinear = std::pow(10.0f, -80.0f / 20.0f);
    }

    std::vector<uint8_t> packChunk(const std::vector<HarmonicObject>& objects) {
        std::vector<uint8_t> data;
        data.push_back(0); // Prediction flag (no prediction)

        for (const auto& obj : objects) {
            uint16_t freqBin = static_cast<uint16_t>(obj.freq / binSize);
            uint8_t hcount = static_cast<uint8_t>(obj.harmonics.size());

            float maxAmp = 0.0f;
            for (const auto& h : obj.harmonics) {
                maxAmp = std::max(maxAmp, std::abs(h.amp));
            }
            maxAmp = std::max(maxAmp, minAmpLinear);

            // Write header
            data.push_back(freqBin & 0xFF);
            data.push_back((freqBin >> 8) & 0xFF);
            data.push_back(hcount);

            // Write max amplitude (float32)
            uint8_t* maxAmpBytes = reinterpret_cast<uint8_t*>(&maxAmp);
            for (int i = 0; i < 4; i++) {
                data.push_back(maxAmpBytes[i]);
            }

            // Write harmonics
            for (const auto& h : obj.harmonics) {
                // Encode amplitude (uint8, log scale)
                float ampSafe = std::max(std::abs(h.amp), minAmpLinear);
                float norm = std::log(ampSafe / minAmpLinear) / std::log(maxAmp / minAmpLinear);
                norm = std::clamp(norm, 0.0f, 1.0f);
                uint8_t encodedAmp = static_cast<uint8_t>(norm * 255);

                // Encode phase (uint8)
                float phaseNorm = (h.phase + static_cast<float>(M_PI)) / (2.0f * static_cast<float>(M_PI));
                uint8_t encodedPhase = static_cast<uint8_t>(phaseNorm * 255);

                data.push_back(encodedAmp);
                data.push_back(encodedPhase);
            }
        }

        return data;
    }
};

// Harmonic Unpacker
class HarmonicUnpacker {
private:
    int sampleRate;
    int windowSize;
    float binSize;
    float minAmpLinear;

public:
    HarmonicUnpacker(int sr, int ws)
        : sampleRate(sr), windowSize(ws) {
        binSize = static_cast<float>(sampleRate) / windowSize;
        minAmpLinear = std::pow(10.0f, -80.0f / 20.0f);
    }

    std::vector<HarmonicObject> unpackChunk(const std::vector<uint8_t>& data) {
        std::vector<HarmonicObject> objects;
        size_t pos = 1; // Skip prediction flag

        while (pos < data.size()) {
            if (pos + 7 > data.size()) break;

            uint16_t freqBin = data[pos] | (data[pos + 1] << 8);
            uint8_t hcount = data[pos + 2];
            pos += 3;

            float maxAmp;
            std::memcpy(&maxAmp, &data[pos], sizeof(float));
            pos += 4;

            HarmonicObject obj;
            obj.freq = freqBin * binSize;

            for (int i = 0; i < hcount; i++) {
                if (pos + 2 > data.size()) break;

                uint8_t encodedAmp = data[pos++];
                uint8_t encodedPhase = data[pos++];

                // Decode amplitude
                float norm = encodedAmp / 255.0f;
                float logRatio = std::log(maxAmp / minAmpLinear);
                float amp = minAmpLinear * std::exp(norm * logRatio);

                // Decode phase
                float phaseNorm = encodedPhase / 255.0f;
                float phase = phaseNorm * 2.0f * static_cast<float>(M_PI) - static_cast<float>(M_PI);

                obj.harmonics.push_back({amp, phase});
            }

            obj.nHarmonic = obj.harmonics.size();
            obj.mainAmp = obj.harmonics.empty() ? 0.0f : obj.harmonics[0].amp;
            objects.push_back(obj);
        }

        return objects;
    }
};

// Harmonic Generator
class HarmonicGenerator {
private:
    int sampleRate;
    int windowSize;
    int hopSize;
    std::vector<float> window;

public:
    HarmonicGenerator(int sr, int ws, int hs)
        : sampleRate(sr), windowSize(ws), hopSize(hs) {
        window = FFT::hannWindow(windowSize);
    }

    std::vector<float> processChunk(const std::vector<HarmonicObject>& objects) {
        std::vector<float> output(windowSize, 0.0f);

        for (const auto& obj : objects) {
            float p0 = obj.freq;
            for (size_t n = 0; n < obj.harmonics.size(); n++) {
                float freq = p0 * (n + 1);
                float amp = obj.harmonics[n].amp;
                float phase = obj.harmonics[n].phase;

                for (int t = 0; t < windowSize; t++) {
                    float time = static_cast<float>(t) / sampleRate;
                    output[t] += amp * std::sin(2.0f * static_cast<float>(M_PI) * freq * time + phase);
                }
            }
        }

        // Apply window
        for (int i = 0; i < windowSize; i++) {
            output[i] *= window[i];
        }

        return output;
    }
};

// Main processing function
void processAudio(const std::string& inputFile, const std::string& outputFile, int numThreads = 0) {
    if (numThreads <= 0) {
        numThreads = std::thread::hardware_concurrency();
        if (numThreads == 0) numThreads = 4;
    }
    
    std::cout << "Reading input file: " << inputFile << std::endl;

    std::vector<float> audio;
    int sampleRate;
    if (!WavFile::read(inputFile, audio, sampleRate)) {
        std::cerr << "Failed to read input file!" << std::endl;
        return;
    }

    std::cout << "Sample rate: " << sampleRate << " Hz" << std::endl;
    std::cout << "Audio length: " << audio.size() << " samples" << std::endl;
    std::cout << "Using " << numThreads << " threads" << std::endl;

    // Parameters
    const int WINDOW_SIZE = 4096;
    const int HOP_SIZE = 512;

    // Calculate total chunks
    int totalChunks = 0;
    for (size_t i = 0; i < audio.size(); i += HOP_SIZE) {
        totalChunks++;
    }
    
    std::cout << "Total chunks to process: " << totalChunks << std::endl;

    // Storage for encoded data
    std::vector<std::vector<uint8_t>> encodedChunks(totalChunks);
    std::atomic<int> processedChunks(0);
    std::atomic<long long> totalBits(0);
    std::mutex progressMutex;

    // --- ANALYSIS PHASE (Multi-threaded) ---
    std::cout << "Extracting harmonics..." << std::endl;
    
    auto analysisWorker = [&](int threadId, int startIdx, int endIdx) {
        HarmonicExtractor extractor(sampleRate, WINDOW_SIZE, HOP_SIZE, 20.0f, 12000.0f, 
                                    0.0f, 20, 12000.0f, 20);
        HarmonicPacker packer(sampleRate, WINDOW_SIZE);
        
        for (int chunkIdx = startIdx; chunkIdx < endIdx; chunkIdx++) {
            size_t i = chunkIdx * HOP_SIZE;
            if (i >= audio.size()) break;
            
            size_t chunkEnd = std::min(i + WINDOW_SIZE, audio.size());
            std::vector<float> chunk(audio.begin() + i, audio.begin() + chunkEnd);

            // Extract harmonics
            auto harmonics = extractor.processChunk(chunk);

            // Pack
            encodedChunks[chunkIdx] = packer.packChunk(harmonics);
            totalBits += encodedChunks[chunkIdx].size() * 8;

            int processed = ++processedChunks;
            if (processed % 100 == 0) {
                std::lock_guard<std::mutex> lock(progressMutex);
                std::cout << "Processed " << processed << "/" << totalChunks << " chunks" << std::endl;
            }
        }
    };

    // Launch analysis threads
    std::vector<std::thread> threads;
    int chunksPerThread = (totalChunks + numThreads - 1) / numThreads;
    
    for (int t = 0; t < numThreads; t++) {
        int startIdx = t * chunksPerThread;
        int endIdx = std::min(startIdx + chunksPerThread, totalChunks);
        if (startIdx < totalChunks) {
            threads.emplace_back(analysisWorker, t, startIdx, endIdx);
        }
    }

    // Wait for all analysis threads
    for (auto& thread : threads) {
        thread.join();
    }
    threads.clear();

    std::cout << "Extraction complete. Total Chunks: " << totalChunks << std::endl;
    
    long long avgBits = totalBits / totalChunks;
    float avgBitrate = (avgBits * (sampleRate / static_cast<float>(WINDOW_SIZE))) / 1000.0f;
    std::cout << "Average bitrate: " << avgBitrate << " kbps" << std::endl;

    // --- SYNTHESIS PHASE (Multi-threaded) ---
    std::cout << "Generating audio..." << std::endl;
    
    std::vector<float> outputSignal(audio.size() + WINDOW_SIZE, 0.0f);
    std::mutex outputMutex;
    std::atomic<int> synthesizedChunks(0);

    auto synthesisWorker = [&](int threadId, int startIdx, int endIdx) {
        HarmonicUnpacker unpacker(sampleRate, WINDOW_SIZE);
        HarmonicGenerator generator(sampleRate, WINDOW_SIZE, HOP_SIZE);
        
        // Local buffer for this thread
        std::vector<float> localOutput(outputSignal.size(), 0.0f);
        
        for (int chunkIdx = startIdx; chunkIdx < endIdx; chunkIdx++) {
            // Unpack
            auto decoded = unpacker.unpackChunk(encodedChunks[chunkIdx]);

            // Generate
            auto synthesized = generator.processChunk(decoded);

            // Add to local buffer
            size_t startIndex = chunkIdx * HOP_SIZE;
            for (size_t t = 0; t < synthesized.size() && (startIndex + t) < localOutput.size(); t++) {
                localOutput[startIndex + t] += synthesized[t];
            }

            int synth = ++synthesizedChunks;
            if (synth % 100 == 0) {
                std::lock_guard<std::mutex> lock(progressMutex);
                std::cout << "Synthesized " << synth << "/" << totalChunks << " chunks" << std::endl;
            }
        }
        
        // Merge local buffer into output (needs synchronization)
        std::lock_guard<std::mutex> lock(outputMutex);
        for (size_t i = 0; i < outputSignal.size(); i++) {
            outputSignal[i] += localOutput[i];
        }
    };

    // Launch synthesis threads
    for (int t = 0; t < numThreads; t++) {
        int startIdx = t * chunksPerThread;
        int endIdx = std::min(startIdx + chunksPerThread, totalChunks);
        if (startIdx < totalChunks) {
            threads.emplace_back(synthesisWorker, t, startIdx, endIdx);
        }
    }

    // Wait for all synthesis threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Normalize and save
    std::cout << "Normalizing and saving output..." << std::endl;
    float maxVal = 0.0f;
    for (float sample : outputSignal) {
        maxVal = std::max(maxVal, std::abs(sample));
    }
    if (maxVal > 1e-6f) {
        for (float& sample : outputSignal) {
            sample /= maxVal;
        }
    }

    // Write output
    std::cout << "Writing output file: " << outputFile << std::endl;
    WavFile::write(outputFile, outputSignal, sampleRate);
    std::cout << "Done!" << std::endl;
}

int main(int argc, char* argv[]) {
    std::string inputFile = "sample.wav";
    std::string outputFile = "output.cpp.wav";
    int numThreads = 0; // 0 = auto-detect

    if (argc >= 2) inputFile = argv[1];
    if (argc >= 3) outputFile = argv[2];
    if (argc >= 4) numThreads = std::atoi(argv[3]);

    processAudio(inputFile, outputFile, numThreads);
    return 0;
}