// CVSD (Continuously Variable Slope Delta modulation) codec
//
// 1-bit-per-sample delta modulator with syllabic step-size adaptation, using the
// Bluetooth (HFP/SCO) parameter set by default:
//   step min 10, step max 1280, run length J = 4,
//   accumulator decay h = 1 - 1/32, step decay beta = 1 - 1/1024.
//
// The core runs at the bit rate. PCM input is upsampled (linear interpolation)
// from the WAV rate to the bit rate before encoding, and the decoder's output is
// low-pass filtered and decimated back to the WAV rate.
//
// Build: g++ -O2 -std=c++17 cvsd.cpp -o cvsd

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// CVSD core
// ---------------------------------------------------------------------------

struct CVSDParams {
    int32_t stepMin = 10;
    int32_t stepMax = 1280;
    int runLength = 4;       // J: equal bits in a row that trigger step growth
    int accumShift = 5;      // h = 1 - 2^-accumShift
    int stepShift = 10;      // beta = 1 - 2^-stepShift
};

// Shared state machine: encoder and decoder run the exact same integer math so
// they stay bit-exact in sync.
class CVSDState {
public:
    explicit CVSDState(const CVSDParams& p = CVSDParams()) : params(p) { reset(); }

    void reset() {
        accum = 0;
        step = params.stepMin;
        history = 0;
    }

    // Feed one bit, return the new reconstructed sample.
    int16_t update(int bit) {
        const uint32_t mask = (1u << params.runLength) - 1;
        history = ((history << 1) | (bit & 1)) & mask;

        // Coincidence: last J bits all ones or all zeros -> slope overload, grow step.
        if (history == 0 || history == mask) {
            step = std::min(step + params.stepMin, params.stepMax);
        } else {
            step = std::max(step - (step >> params.stepShift), params.stepMin);
        }

        // Leaky integrator.
        accum -= accum >> params.accumShift;
        accum += bit ? step : -step;
        accum = std::clamp<int32_t>(accum, -32768, 32767);
        return static_cast<int16_t>(accum);
    }

    int16_t estimate() const { return static_cast<int16_t>(accum); }

private:
    CVSDParams params;
    int32_t accum;
    int32_t step;
    uint32_t history;
};

class CVSDEncoder {
public:
    explicit CVSDEncoder(const CVSDParams& p = CVSDParams()) : state(p) {}

    void reset() { state.reset(); }

    int encodeSample(int16_t sample) {
        int bit = sample >= state.estimate() ? 1 : 0;
        state.update(bit);
        return bit;
    }

    std::vector<uint8_t> encode(const std::vector<int16_t>& samples) {
        std::vector<uint8_t> bits(samples.size());
        for (size_t i = 0; i < samples.size(); i++) {
            bits[i] = static_cast<uint8_t>(encodeSample(samples[i]));
        }
        return bits;
    }

private:
    CVSDState state;
};

class CVSDDecoder {
public:
    explicit CVSDDecoder(const CVSDParams& p = CVSDParams()) : state(p) {}

    void reset() { state.reset(); }

    int16_t decodeBit(int bit) { return state.update(bit); }

    std::vector<int16_t> decode(const std::vector<uint8_t>& bits) {
        std::vector<int16_t> samples(bits.size());
        for (size_t i = 0; i < bits.size(); i++) {
            samples[i] = decodeBit(bits[i]);
        }
        return samples;
    }

private:
    CVSDState state;
};

// ---------------------------------------------------------------------------
// Rate conversion helpers
// ---------------------------------------------------------------------------

// 2nd-order Butterworth low-pass (RBJ cookbook biquad).
class Biquad {
public:
    Biquad(double sampleRate, double cutoff) {
        const double pi = 3.14159265358979323846;
        double w0 = 2.0 * pi * cutoff / sampleRate;
        double alpha = std::sin(w0) / (2.0 * std::sqrt(0.5));
        double c = std::cos(w0);
        double a0 = 1.0 + alpha;
        b0 = (1.0 - c) / 2.0 / a0;
        b1 = (1.0 - c) / a0;
        b2 = b0;
        a1 = -2.0 * c / a0;
        a2 = (1.0 - alpha) / a0;
    }

    double process(double x) {
        double y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        x2 = x1; x1 = x;
        y2 = y1; y1 = y;
        return y;
    }

private:
    double b0, b1, b2, a1, a2;
    double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
};

std::vector<int16_t> upsample(const std::vector<int16_t>& in, int factor) {
    std::vector<int16_t> out;
    out.reserve(in.size() * factor);
    for (size_t i = 0; i < in.size(); i++) {
        double a = in[i];
        double b = (i + 1 < in.size()) ? in[i + 1] : in[i];
        for (int k = 0; k < factor; k++) {
            out.push_back(static_cast<int16_t>(std::lround(a + (b - a) * k / factor)));
        }
    }
    return out;
}

std::vector<int16_t> downsample(const std::vector<int16_t>& in, int factor, double outRate) {
    // Two cascaded biquads -> 4th-order Butterworth anti-alias filter.
    const double pi = 3.14159265358979323846;
    const double cutoff = outRate * 0.45;
    double inRate = outRate * factor;
    Biquad f1(inRate, cutoff), f2(inRate, cutoff);

    // Compensate the filter's group delay (2.6131 / wc for a 4th-order Butterworth)
    // so the output lines up with the original PCM.
    size_t delay = static_cast<size_t>(std::lround(2.6131 / (2.0 * pi * cutoff) * inRate));

    std::vector<int16_t> out;
    out.reserve(in.size() / factor + 1);
    for (size_t i = 0; i < in.size() + delay; i++) {
        double y = f2.process(f1.process(in[std::min(i, in.size() - 1)]));
        if (i >= delay && (i - delay) % factor == 0) {
            out.push_back(static_cast<int16_t>(std::clamp(std::lround(y), -32768L, 32767L)));
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// WAV I/O (16-bit PCM, any channel count in, mono out)
// ---------------------------------------------------------------------------

struct Audio {
    int sampleRate = 8000;
    std::vector<int16_t> samples;
};

template <typename T>
static T readLE(std::istream& in) {
    T v{};
    in.read(reinterpret_cast<char*>(&v), sizeof(T));
    return v;
}

template <typename T>
static void writeLE(std::ostream& out, T v) {
    out.write(reinterpret_cast<const char*>(&v), sizeof(T));
}

Audio readWav(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + path);

    char id[4];
    in.read(id, 4);
    readLE<uint32_t>(in);
    char wave[4];
    in.read(wave, 4);
    if (std::memcmp(id, "RIFF", 4) != 0 || std::memcmp(wave, "WAVE", 4) != 0) {
        throw std::runtime_error(path + " is not a RIFF/WAVE file");
    }

    Audio audio;
    int channels = 0, bits = 0;
    bool haveFmt = false;
    while (in.read(id, 4)) {
        uint32_t size = readLE<uint32_t>(in);
        if (std::memcmp(id, "fmt ", 4) == 0) {
            uint16_t format = readLE<uint16_t>(in);
            channels = readLE<uint16_t>(in);
            audio.sampleRate = static_cast<int>(readLE<uint32_t>(in));
            readLE<uint32_t>(in);
            readLE<uint16_t>(in);
            bits = readLE<uint16_t>(in);
            in.seekg(size - 16 + (size & 1), std::ios::cur);
            if (format != 1 || bits != 16) {
                throw std::runtime_error("only 16-bit PCM WAV is supported");
            }
            haveFmt = true;
        } else if (std::memcmp(id, "data", 4) == 0) {
            if (!haveFmt) throw std::runtime_error("data chunk before fmt chunk");
            size_t frames = size / (2 * channels);
            audio.samples.resize(frames);
            for (size_t i = 0; i < frames; i++) {
                int32_t sum = 0;
                for (int c = 0; c < channels; c++) sum += readLE<int16_t>(in);
                audio.samples[i] = static_cast<int16_t>(sum / channels);
            }
            return audio;
        } else {
            in.seekg(size + (size & 1), std::ios::cur);
        }
    }
    throw std::runtime_error("no data chunk in " + path);
}

void writeWav(const std::string& path, const Audio& audio) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write " + path);

    uint32_t dataSize = static_cast<uint32_t>(audio.samples.size() * 2);
    out.write("RIFF", 4);
    writeLE<uint32_t>(out, 36 + dataSize);
    out.write("WAVE", 4);
    out.write("fmt ", 4);
    writeLE<uint32_t>(out, 16);
    writeLE<uint16_t>(out, 1);
    writeLE<uint16_t>(out, 1);
    writeLE<uint32_t>(out, audio.sampleRate);
    writeLE<uint32_t>(out, audio.sampleRate * 2);
    writeLE<uint16_t>(out, 2);
    writeLE<uint16_t>(out, 16);
    out.write("data", 4);
    writeLE<uint32_t>(out, dataSize);
    out.write(reinterpret_cast<const char*>(audio.samples.data()), dataSize);
}

// ---------------------------------------------------------------------------
// .cvsd bitstream file
//   "CVSD" | u32 pcm sample rate | u32 bit rate | u32 bit count | packed bits (MSB first)
// ---------------------------------------------------------------------------

struct CVSDStream {
    uint32_t sampleRate = 8000;
    uint32_t bitRate = 64000;
    std::vector<uint8_t> bits;
};

void writeStream(const std::string& path, const CVSDStream& s) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write " + path);

    out.write("CVSD", 4);
    writeLE<uint32_t>(out, s.sampleRate);
    writeLE<uint32_t>(out, s.bitRate);
    writeLE<uint32_t>(out, static_cast<uint32_t>(s.bits.size()));

    std::vector<uint8_t> packed((s.bits.size() + 7) / 8, 0);
    for (size_t i = 0; i < s.bits.size(); i++) {
        if (s.bits[i]) packed[i / 8] |= static_cast<uint8_t>(0x80 >> (i % 8));
    }
    out.write(reinterpret_cast<const char*>(packed.data()), packed.size());
}

CVSDStream readStream(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + path);

    char magic[4];
    in.read(magic, 4);
    if (std::memcmp(magic, "CVSD", 4) != 0) throw std::runtime_error(path + " is not a CVSD stream");

    CVSDStream s;
    s.sampleRate = readLE<uint32_t>(in);
    s.bitRate = readLE<uint32_t>(in);
    uint32_t count = readLE<uint32_t>(in);

    std::vector<uint8_t> packed((count + 7) / 8);
    in.read(reinterpret_cast<char*>(packed.data()), packed.size());
    if (!in) throw std::runtime_error(path + " is truncated");

    s.bits.resize(count);
    for (uint32_t i = 0; i < count; i++) {
        s.bits[i] = (packed[i / 8] >> (7 - i % 8)) & 1;
    }
    return s;
}

// ---------------------------------------------------------------------------
// High-level encode/decode
// ---------------------------------------------------------------------------

CVSDStream encodeAudio(const Audio& audio, uint32_t bitRate) {
    if (audio.sampleRate <= 0 || bitRate == 0 || bitRate % audio.sampleRate != 0) {
        throw std::runtime_error("bit rate must be an integer multiple of the sample rate");
    }
    int factor = static_cast<int>(bitRate / audio.sampleRate);

    CVSDEncoder encoder;
    CVSDStream s;
    s.sampleRate = audio.sampleRate;
    s.bitRate = bitRate;
    s.bits = encoder.encode(upsample(audio.samples, factor));
    return s;
}

Audio decodeStream(const CVSDStream& s) {
    if (s.sampleRate == 0 || s.bitRate == 0 || s.bitRate % s.sampleRate != 0) {
        throw std::runtime_error("invalid rates in CVSD stream header");
    }
    int factor = static_cast<int>(s.bitRate / s.sampleRate);

    CVSDDecoder decoder;
    Audio audio;
    audio.sampleRate = static_cast<int>(s.sampleRate);
    audio.samples = downsample(decoder.decode(s.bits), factor, s.sampleRate);
    return audio;
}

double snrDb(const std::vector<int16_t>& ref, const std::vector<int16_t>& test) {
    size_t n = std::min(ref.size(), test.size());
    double sig = 0, noise = 0;
    for (size_t i = 0; i < n; i++) {
        double d = static_cast<double>(ref[i]) - test[i];
        sig += static_cast<double>(ref[i]) * ref[i];
        noise += d * d;
    }
    return noise > 0 ? 10.0 * std::log10(sig / noise) : INFINITY;
}

// ---------------------------------------------------------------------------

static void usage(const char* prog) {
    std::cout << "Usage:\n"
              << "  " << prog << " encode input.wav output.cvsd [bitrate]\n"
              << "  " << prog << " decode input.cvsd output.wav\n"
              << "  " << prog << " roundtrip input.wav output.wav [bitrate]\n"
              << "bitrate defaults to 64000 and must be a multiple of the WAV sample rate.\n";
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        usage(argv[0]);
        return 1;
    }

    std::string mode = argv[1];
    try {
        if (mode == "encode") {
            uint32_t bitRate = argc > 4 ? static_cast<uint32_t>(std::stoul(argv[4])) : 64000;
            Audio audio = readWav(argv[2]);
            CVSDStream s = encodeAudio(audio, bitRate);
            writeStream(argv[3], s);
            std::cout << "Encoded " << audio.samples.size() << " samples @ " << audio.sampleRate
                      << " Hz -> " << s.bits.size() << " bits @ " << bitRate << " bps" << std::endl;
        } else if (mode == "decode") {
            CVSDStream s = readStream(argv[2]);
            Audio audio = decodeStream(s);
            writeWav(argv[3], audio);
            std::cout << "Decoded " << s.bits.size() << " bits -> " << audio.samples.size()
                      << " samples @ " << audio.sampleRate << " Hz" << std::endl;
        } else if (mode == "roundtrip") {
            uint32_t bitRate = argc > 4 ? static_cast<uint32_t>(std::stoul(argv[4])) : 64000;
            Audio audio = readWav(argv[2]);
            Audio decoded = decodeStream(encodeAudio(audio, bitRate));
            writeWav(argv[3], decoded);
            std::cout << "Round trip @ " << bitRate << " bps, SNR: "
                      << snrDb(audio.samples, decoded.samples) << " dB" << std::endl;
        } else {
            usage(argv[0]);
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
