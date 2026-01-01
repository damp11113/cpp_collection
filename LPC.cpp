#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdint>

class LPCCodec {
private:
    static constexpr int ORDER = 10;  // LPC order
    static constexpr float MIN_ENERGY = 1e-10f;

    // Autocorrelation calculation
    std::vector<float> computeAutocorrelation(const std::vector<float>& frame) {
        std::vector<float> r(ORDER + 1, 0.0f);
        for (int lag = 0; lag <= ORDER; lag++) {
            for (size_t i = lag; i < frame.size(); i++) {
                r[lag] += frame[i] * frame[i - lag];
            }
        }
        return r;
    }

    // Levinson-Durbin recursion for LPC coefficients
    std::vector<float> levinsonDurbin(const std::vector<float>& r) {
        std::vector<float> lpc(ORDER + 1, 0.0f);
        std::vector<float> temp(ORDER + 1, 0.0f);
        
        if (std::abs(r[0]) < MIN_ENERGY) {
            return lpc;
        }

        float k;
        float error = r[0];

        lpc[0] = 1.0f;
        
        for (int i = 1; i <= ORDER; i++) {
            k = r[i];
            for (int j = 1; j < i; j++) {
                k += lpc[j] * r[i - j];
            }
            k = -k / error;

            for (int j = 1; j < i; j++) {
                temp[j] = lpc[j] + k * lpc[i - j];
            }
            for (int j = 1; j < i; j++) {
                lpc[j] = temp[j];
            }

            lpc[i] = k;
            error *= (1.0f - k * k);
        }

        return lpc;
    }

public:
    struct EncodedFrame {
        std::vector<float> coefficients;
        std::vector<float> residual;
        float gain;
    };

    // Encode a frame of audio
    EncodedFrame encode(const std::vector<float>& frame) {
        EncodedFrame encoded;
        
        // Compute autocorrelation
        std::vector<float> r = computeAutocorrelation(frame);
        
        // Get LPC coefficients
        encoded.coefficients = levinsonDurbin(r);
        
        // Calculate residual
        encoded.residual.resize(frame.size());
        for (size_t i = 0; i < frame.size(); i++) {
            float predicted = 0.0f;
            for (int j = 1; j <= std::min(static_cast<int>(i), ORDER); j++) {
                predicted += encoded.coefficients[j] * frame[i - j];
            }
            encoded.residual[i] = frame[i] - predicted;
        }

        // Calculate gain
        float sum = 0.0f;
        for (float sample : encoded.residual) {
            sum += sample * sample;
        }
        encoded.gain = std::sqrt(sum / frame.size());

        // Normalize residual
        if (encoded.gain > MIN_ENERGY) {
            for (float& sample : encoded.residual) {
                sample /= encoded.gain;
            }
        }

        return encoded;
    }

    // Decode a frame of audio
    std::vector<float> decode(const EncodedFrame& encoded) {
        std::vector<float> output(encoded.residual.size(), 0.0f);
        
        // Denormalize residual
        std::vector<float> scaled_residual = encoded.residual;
        for (float& sample : scaled_residual) {
            sample *= encoded.gain;
        }

        // Synthesis
        for (size_t i = 0; i < output.size(); i++) {
            output[i] = scaled_residual[i];
            for (int j = 1; j <= std::min(static_cast<int>(i), ORDER); j++) {
                output[i] += encoded.coefficients[j] * output[i - j];
            }
        }

        return output;
    }
};

class WAVProcessor {
private:
    static constexpr int FRAME_SIZE = 480;  // 30ms at 16kHz
    static constexpr int FRAME_OVERLAP = 40;  // 2.5ms overlap

public:
	struct WAVHeader {
        // RIFF chunk
        char riff_header[4];    // Contains "RIFF"
        uint32_t wav_size;      // Size of WAV file
        char wave_header[4];    // Contains "WAVE"
        
        // Format chunk
        char fmt_header[4];     // Contains "fmt "
        uint32_t fmt_chunk_size;
        uint16_t audio_format;
        uint16_t num_channels;
        uint32_t sample_rate;
        uint32_t byte_rate;
        uint16_t sample_alignment;
        uint16_t bit_depth;
        
        // Data chunk
        char data_header[4];    // Contains "data"
        uint32_t data_bytes;
    };

    static std::vector<float> readWAV(const std::string& filename, WAVHeader& header) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open WAV file for reading");
        }

        // Read WAV header
        file.read(reinterpret_cast<char*>(&header), sizeof(WAVHeader));

        // Verify format
        if (header.audio_format != 1) { // PCM = 1
            throw std::runtime_error("Only PCM WAV files are supported");
        }

        // Read audio data
        std::vector<float> audio_data;
        const size_t num_samples = header.data_bytes / (header.bit_depth / 8);
        audio_data.reserve(num_samples);

        if (header.bit_depth == 16) {
            int16_t sample;
            while (file.read(reinterpret_cast<char*>(&sample), sizeof(int16_t))) {
                audio_data.push_back(sample / 32768.0f);
            }
        } else {
            throw std::runtime_error("Only 16-bit WAV files are supported");
        }

        file.close();
        return audio_data;
    }

    static void writeWAV(const std::string& filename, const std::vector<float>& audio_data, const WAVHeader& header) {
        std::ofstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open WAV file for writing");
        }

        // Write header
        WAVHeader new_header = header;
        new_header.data_bytes = audio_data.size() * (header.bit_depth / 8);
        new_header.wav_size = new_header.data_bytes + sizeof(WAVHeader) - 8;
        file.write(reinterpret_cast<const char*>(&new_header), sizeof(WAVHeader));

        // Write audio data
        for (float sample : audio_data) {
            int16_t pcm_sample = static_cast<int16_t>(std::clamp(sample * 32768.0f, -32768.0f, 32767.0f));
            file.write(reinterpret_cast<const char*>(&pcm_sample), sizeof(int16_t));
        }

        file.close();
    }

    static std::vector<std::vector<float>> splitIntoFrames(const std::vector<float>& audio_data) {
        std::vector<std::vector<float>> frames;
        size_t pos = 0;
        
        while (pos < audio_data.size()) {
            std::vector<float> frame;
            frame.reserve(FRAME_SIZE);
            
            for (int i = 0; i < FRAME_SIZE && (pos + i) < audio_data.size(); i++) {
                frame.push_back(audio_data[pos + i]);
            }
            
            // Pad last frame if necessary
            while (frame.size() < FRAME_SIZE) {
                frame.push_back(0.0f);
            }
            
            frames.push_back(frame);
            pos += (FRAME_SIZE - FRAME_OVERLAP);
        }
        
        return frames;
    }

    static std::vector<float> combineFrames(const std::vector<std::vector<float>>& frames) {
        if (frames.empty()) return std::vector<float>();
        
        size_t total_size = (frames.size() - 1) * (FRAME_SIZE - FRAME_OVERLAP) + FRAME_SIZE;
        std::vector<float> audio_data(total_size, 0.0f);
        std::vector<float> overlap_count(total_size, 0.0f);
        
        size_t pos = 0;
        for (const auto& frame : frames) {
            for (size_t i = 0; i < frame.size(); i++) {
                audio_data[pos + i] += frame[i];
                overlap_count[pos + i] += 1.0f;
            }
            pos += (FRAME_SIZE - FRAME_OVERLAP);
        }
        
        // Average overlapping regions
        for (size_t i = 0; i < audio_data.size(); i++) {
            if (overlap_count[i] > 0) {
                audio_data[i] /= overlap_count[i];
            }
        }
        
        return audio_data;
    }
};

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <input.wav> <output.wav>\n";
        return 1;
    }

    try {
        // Create codec instance
        LPCCodec codec;
        WAVProcessor::WAVHeader header;

        // Read input WAV file
        std::cout << "Reading WAV file...\n";
        auto audio_data = WAVProcessor::readWAV(argv[1], header);

        // Split into frames
        std::cout << "Processing frames...\n";
        auto frames = WAVProcessor::splitIntoFrames(audio_data);

        // Process each frame
        std::vector<std::vector<float>> processed_frames;
        for (const auto& frame : frames) {
            // Encode
            auto encoded = codec.encode(frame);
            
            // Decode
            auto decoded = codec.decode(encoded);
            processed_frames.push_back(decoded);
        }

        // Combine frames
        std::cout << "Combining frames...\n";
        auto output_audio = WAVProcessor::combineFrames(processed_frames);

        // Write output WAV file
        std::cout << "Writing WAV file...\n";
        WAVProcessor::writeWAV(argv[2], output_audio, header);

        std::cout << "Processing complete!\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}