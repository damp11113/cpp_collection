// Simple I/Q AM/FM/SSB modulator and demodulator (complex baseband)
//
// Build: g++ -O2 -std=c++17 iq_modem.cpp -o iq_modem
//
// Usage:
//   iq_modem test
//       Run a self-test: tone -> modulate -> noise -> demodulate, prints SNR.
//   iq_modem mod   <am|fm|usb|lsb> <in.wav> <out.iq> [param]
//       Modulate a 16-bit PCM WAV into interleaved float32 I/Q (cf32, GNU Radio / SDR# compatible).
//       param = AM modulation index (default 0.8) or FM deviation in Hz (default 5000).
//   iq_modem demod <am|fm|usb|lsb> <in.iq> <out.wav> <sample_rate> [param]
//       Demodulate interleaved float32 I/Q into a 16-bit PCM mono WAV.
//       param = FM deviation in Hz (default 5000), ignored for AM/SSB.

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <random>
#include <algorithm>

using namespace std;

typedef complex<float> iq_t;

const double PI = 3.14159265358979323846;

// ---------------------------------------------------------------------------
// AM
// ---------------------------------------------------------------------------

// AM (DSB with carrier): s(t) = (1 + m * x(t)), carrier at 0 Hz -> purely real I, Q = 0.
vector<iq_t> amModulate(const vector<float>& audio, float modIndex) {
    vector<iq_t> iq(audio.size());
    for (size_t i = 0; i < audio.size(); ++i) {
        iq[i] = iq_t(1.0f + modIndex * audio[i], 0.0f);
    }
    return iq;
}

// AM envelope detector: |I + jQ|, then a DC blocker to remove the carrier.
// Using the magnitude makes it immune to carrier phase offset.
vector<float> amDemodulate(const vector<iq_t>& iq) {
    vector<float> audio(iq.size());
    // Start the blocker at the first envelope value so the carrier isn't seen as a step
    float prevIn = iq.empty() ? 0.0f : abs(iq[0]), prevOut = 0.0f;
    const float R = 0.999f; // DC blocker pole
    for (size_t i = 0; i < iq.size(); ++i) {
        float env = abs(iq[i]);
        float out = env - prevIn + R * prevOut;
        prevIn = env;
        prevOut = out;
        audio[i] = out;
    }
    return audio;
}

// ---------------------------------------------------------------------------
// FM
// ---------------------------------------------------------------------------

// FM: phase is the running integral of the message, phi[n] = phi[n-1] + 2*pi*dev*x[n]/fs.
vector<iq_t> fmModulate(const vector<float>& audio, double sampleRate, double deviation) {
    vector<iq_t> iq(audio.size());
    double phase = 0.0;
    const double k = 2.0 * PI * deviation / sampleRate;
    for (size_t i = 0; i < audio.size(); ++i) {
        phase += k * audio[i];
        phase = fmod(phase, 2.0 * PI); // keep precision over long signals
        iq[i] = iq_t((float)cos(phase), (float)sin(phase));
    }
    return iq;
}

// FM quadrature discriminator: angle of z[n] * conj(z[n-1]) is the phase step,
// which is proportional to instantaneous frequency.
vector<float> fmDemodulate(const vector<iq_t>& iq, double sampleRate, double deviation) {
    vector<float> audio(iq.size());
    const double gain = sampleRate / (2.0 * PI * deviation);
    iq_t prev(1.0f, 0.0f);
    for (size_t i = 0; i < iq.size(); ++i) {
        iq_t d = iq[i] * conj(prev);
        audio[i] = (float)(gain * atan2(d.imag(), d.real()));
        prev = iq[i];
    }
    return audio;
}

// ---------------------------------------------------------------------------
// SSB
// ---------------------------------------------------------------------------

const int HILBERT_RADIUS = 127; // 255-tap FIR

// Windowed FIR Hilbert transformer (90 degree phase shift).
// Ideal response is h[n] = 2/(pi*n) for odd n, 0 for even n; a Blackman window
// tames the ripple. Applied centered (non-causal), so output lines up with input.
vector<float> hilbert(const vector<float>& in) {
    static vector<float> h;
    if (h.empty()) {
        h.resize(2 * HILBERT_RADIUS + 1, 0.0f);
        const int N = 2 * HILBERT_RADIUS;
        for (int k = -HILBERT_RADIUS; k <= HILBERT_RADIUS; ++k) {
            if (k % 2 == 0) continue;
            double m = k + HILBERT_RADIUS;
            double w = 0.42 - 0.5 * cos(2.0 * PI * m / N) + 0.08 * cos(4.0 * PI * m / N);
            h[k + HILBERT_RADIUS] = (float)(w * 2.0 / (PI * k));
        }
    }

    const long n = (long)in.size();
    vector<float> out(n, 0.0f);
    for (long i = 0; i < n; ++i) {
        double acc = 0.0;
        for (int k = -HILBERT_RADIUS; k <= HILBERT_RADIUS; k += 2) { // even taps are zero
            if (k == 0) continue;
            long j = i - k;
            if (j < 0 || j >= n) continue;
            acc += h[k + HILBERT_RADIUS] * in[j];
        }
        out[i] = (float)acc;
    }
    return out;
}

// SSB via the phasing method: the analytic signal x + j*H{x} only has positive
// frequencies (USB); x - j*H{x} only has negative frequencies (LSB).
vector<iq_t> ssbModulate(const vector<float>& audio, bool upper) {
    vector<float> q = hilbert(audio);
    vector<iq_t> iq(audio.size());
    for (size_t i = 0; i < audio.size(); ++i) {
        iq[i] = iq_t(audio[i], upper ? q[i] : -q[i]);
    }
    return iq;
}

// SSB demod: I - H{Q} keeps positive frequencies (USB), I + H{Q} keeps negative
// ones (LSB). The opposite sideband cancels, so it is rejected.
// Unlike the AM envelope detector this needs the carrier phase to be correct;
// a phase error only phase-shifts the audio, which the ear barely notices.
vector<float> ssbDemodulate(const vector<iq_t>& iq, bool upper) {
    vector<float> i(iq.size()), q(iq.size());
    for (size_t n = 0; n < iq.size(); ++n) {
        i[n] = iq[n].real();
        q[n] = iq[n].imag();
    }
    vector<float> hq = hilbert(q);

    vector<float> audio(iq.size());
    for (size_t n = 0; n < iq.size(); ++n) {
        audio[n] = 0.5f * (upper ? i[n] - hq[n] : i[n] + hq[n]);
    }
    return audio;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Centered moving-average low-pass to clean up demodulated audio.
// Averages `radius` samples on each side, so it adds no delay.
vector<float> movingAverage(const vector<float>& in, int radius) {
    if (radius <= 0 || in.empty()) return in;
    const long n = (long)in.size();
    vector<double> prefix(n + 1, 0.0);
    for (long i = 0; i < n; ++i) prefix[i + 1] = prefix[i] + in[i];

    vector<float> out(n);
    for (long i = 0; i < n; ++i) {
        long lo = max(0L, i - radius);
        long hi = min(n - 1, i + radius);
        out[i] = (float)((prefix[hi + 1] - prefix[lo]) / (hi - lo + 1));
    }
    return out;
}

void addNoise(vector<iq_t>& iq, float sigma, unsigned seed = 1234) {
    mt19937 rng(seed);
    normal_distribution<float> n(0.0f, sigma);
    for (auto& s : iq) s += iq_t(n(rng), n(rng));
}

// Rotate the whole signal by a fixed phase (simulates unknown carrier phase).
void rotatePhase(vector<iq_t>& iq, float radians) {
    iq_t r = polar(1.0f, radians);
    for (auto& s : iq) s *= r;
}

// Best-fit SNR of `test` against `ref` (removes gain and DC before comparing).
double measureSnrDb(const vector<float>& ref, const vector<float>& test, size_t skip) {
    size_t n = min(ref.size(), test.size());
    double mr = 0, mt = 0;
    for (size_t i = skip; i < n; ++i) { mr += ref[i]; mt += test[i]; }
    mr /= (n - skip); mt /= (n - skip);

    double rr = 0, rt = 0;
    for (size_t i = skip; i < n; ++i) {
        rr += (ref[i] - mr) * (ref[i] - mr);
        rt += (ref[i] - mr) * (test[i] - mt);
    }
    double g = rt / rr;

    double sig = 0, err = 0;
    for (size_t i = skip; i < n; ++i) {
        double r = g * (ref[i] - mr);
        double e = (test[i] - mt) - r;
        sig += r * r;
        err += e * e;
    }
    return 10.0 * log10(sig / max(err, 1e-20));
}

// ---------------------------------------------------------------------------
// File I/O
// ---------------------------------------------------------------------------

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

// Reads 16-bit PCM WAV, mixes to mono, returns samples in [-1, 1].
bool readWav(const string& path, vector<float>& audio, int& sampleRate) {
    ifstream f(path, ios::binary);
    if (!f) return false;

    char riff[12];
    f.read(riff, 12);
    if (!f || strncmp(riff, "RIFF", 4) != 0 || strncmp(riff + 8, "WAVE", 4) != 0) return false;

    int16_t channels = 0, bits = 0, format = 0;
    int32_t rate = 0;
    // Walk chunks so files with LIST/fact chunks still load
    while (f) {
        char id[4];
        int32_t size;
        f.read(id, 4);
        f.read((char*)&size, 4);
        if (!f) return false;

        if (strncmp(id, "fmt ", 4) == 0) {
            vector<char> buf(size);
            f.read(buf.data(), size);
            memcpy(&format, &buf[0], 2);
            memcpy(&channels, &buf[2], 2);
            memcpy(&rate, &buf[4], 4);
            memcpy(&bits, &buf[14], 2);
        } else if (strncmp(id, "data", 4) == 0) {
            if (format != 1 || bits != 16 || channels < 1) {
                cerr << "Only 16-bit PCM WAV is supported" << endl;
                return false;
            }
            size_t frames = size / (2 * channels);
            vector<int16_t> pcm(frames * channels);
            f.read((char*)pcm.data(), pcm.size() * 2);
            frames = f.gcount() / (2 * channels);

            audio.resize(frames);
            for (size_t i = 0; i < frames; ++i) {
                float sum = 0;
                for (int c = 0; c < channels; ++c) sum += pcm[i * channels + c];
                audio[i] = sum / (channels * 32768.0f);
            }
            sampleRate = rate;
            return true;
        } else {
            f.seekg(size + (size & 1), ios::cur);
        }
    }
    return false;
}

bool writeWav(const string& path, const vector<float>& audio, int sampleRate) {
    ofstream f(path, ios::binary);
    if (!f) return false;

    WAVHeader h;
    memcpy(h.riff, "RIFF", 4);
    memcpy(h.wave, "WAVE", 4);
    memcpy(h.fmt, "fmt ", 4);
    memcpy(h.data, "data", 4);
    h.subchunk1Size = 16;
    h.audioFormat = 1;
    h.numChannels = 1;
    h.sampleRate = sampleRate;
    h.bitsPerSample = 16;
    h.blockAlign = 2;
    h.byteRate = sampleRate * 2;
    h.dataSize = (int32_t)(audio.size() * 2);
    h.chunkSize = 36 + h.dataSize;
    f.write((char*)&h, sizeof(h));

    for (float s : audio) {
        float c = max(-1.0f, min(1.0f, s));
        int16_t v = (int16_t)lrintf(c * 32767.0f);
        f.write((char*)&v, 2);
    }
    return true;
}

// Interleaved float32 I, Q, I, Q, ... (cf32)
bool writeIq(const string& path, const vector<iq_t>& iq) {
    ofstream f(path, ios::binary);
    if (!f) return false;
    f.write((const char*)iq.data(), iq.size() * sizeof(iq_t));
    return true;
}

bool readIq(const string& path, vector<iq_t>& iq) {
    ifstream f(path, ios::binary | ios::ate);
    if (!f) return false;
    size_t bytes = f.tellg();
    f.seekg(0);
    iq.resize(bytes / sizeof(iq_t));
    f.read((char*)iq.data(), iq.size() * sizeof(iq_t));
    return true;
}

// Scale so the loudest sample sits at ~0.9 full scale
void normalize(vector<float>& audio) {
    float peak = 0;
    for (float s : audio) peak = max(peak, fabs(s));
    if (peak > 0) for (float& s : audio) s *= 0.9f / peak;
}

// ---------------------------------------------------------------------------
// Commands
// ---------------------------------------------------------------------------

int runTest() {
    const double fs = 48000.0;
    const size_t n = (size_t)fs; // 1 second
    const double toneHz = 1000.0;

    vector<float> tone(n);
    for (size_t i = 0; i < n; ++i) {
        tone[i] = 0.8f * (float)sin(2.0 * PI * toneHz * i / fs);
    }

    const size_t skip = 4800; // ignore DC blocker settling
    const float noise = 0.02f;

    // AM with unknown carrier phase - envelope detector shouldn't care
    vector<iq_t> am = amModulate(tone, 0.8f);
    rotatePhase(am, 1.1f);
    addNoise(am, noise);
    vector<float> amOut = movingAverage(amDemodulate(am), 4);
    double amSnr = measureSnrDb(tone, amOut, skip);

    // FM
    const double dev = 5000.0;
    vector<iq_t> fm = fmModulate(tone, fs, dev);
    rotatePhase(fm, 1.1f);
    addNoise(fm, noise);
    vector<float> fmOut = movingAverage(fmDemodulate(fm, fs, dev), 4);
    double fmSnr = measureSnrDb(tone, fmOut, skip);

    // SSB: no phase rotation here, the product detector needs a locked carrier
    vector<iq_t> usb = ssbModulate(tone, true);
    addNoise(usb, noise);
    vector<float> usbOut = ssbDemodulate(usb, true);
    double usbSnr = measureSnrDb(tone, usbOut, skip);

    vector<iq_t> lsb = ssbModulate(tone, false);
    addNoise(lsb, noise);
    vector<float> lsbOut = ssbDemodulate(lsb, false);
    double lsbSnr = measureSnrDb(tone, lsbOut, skip);

    // Opposite sideband rejection: demodulate the USB signal as LSB
    vector<float> wrongOut = ssbDemodulate(ssbModulate(tone, true), false);
    double rightPow = 0, wrongPow = 0;
    for (size_t i = skip; i < n - skip; ++i) {
        rightPow += (double)tone[i] * tone[i];
        wrongPow += (double)wrongOut[i] * wrongOut[i];
    }
    double rejectDb = 10.0 * log10(rightPow / max(wrongPow, 1e-20));

    cout << "Self-test: 1 kHz tone, fs = " << fs << " Hz, noise sigma = " << noise << endl;
    cout << "  First I/Q samples (AM): ";
    for (int i = 0; i < 3; ++i) cout << am[i] << " ";
    cout << endl << "  First I/Q samples (FM): ";
    for (int i = 0; i < 3; ++i) cout << fm[i] << " ";
    cout << endl;
    cout << "  AM demod SNR: " << amSnr << " dB" << endl;
    cout << "  FM demod SNR: " << fmSnr << " dB" << endl;
    cout << "  USB demod SNR: " << usbSnr << " dB" << endl;
    cout << "  LSB demod SNR: " << lsbSnr << " dB" << endl;
    cout << "  Opposite sideband rejection: " << rejectDb << " dB" << endl;

    bool ok = amSnr > 20.0 && fmSnr > 20.0 && usbSnr > 20.0 && lsbSnr > 20.0 && rejectDb > 30.0;
    cout << (ok ? "PASS" : "FAIL") << endl;
    return ok ? 0 : 1;
}

int runMod(const string& mode, const string& inPath, const string& outPath, double param) {
    vector<float> audio;
    int fs;
    if (!readWav(inPath, audio, fs)) {
        cerr << "Failed to read " << inPath << endl;
        return 1;
    }

    vector<iq_t> iq;
    if (mode == "am") {
        iq = amModulate(audio, (float)(param > 0 ? param : 0.8));
    } else if (mode == "fm") {
        iq = fmModulate(audio, fs, param > 0 ? param : 5000.0);
    } else if (mode == "usb" || mode == "lsb") {
        iq = ssbModulate(audio, mode == "usb");
    } else {
        cerr << "Unknown mode: " << mode << endl;
        return 1;
    }

    if (!writeIq(outPath, iq)) {
        cerr << "Failed to write " << outPath << endl;
        return 1;
    }
    cout << "Wrote " << iq.size() << " I/Q samples @ " << fs << " Hz to " << outPath << endl;
    return 0;
}

int runDemod(const string& mode, const string& inPath, const string& outPath, int fs, double param) {
    vector<iq_t> iq;
    if (!readIq(inPath, iq)) {
        cerr << "Failed to read " << inPath << endl;
        return 1;
    }

    vector<float> audio;
    if (mode == "am") {
        audio = amDemodulate(iq);
    } else if (mode == "fm") {
        audio = fmDemodulate(iq, fs, param > 0 ? param : 5000.0);
    } else if (mode == "usb" || mode == "lsb") {
        audio = ssbDemodulate(iq, mode == "usb");
    } else {
        cerr << "Unknown mode: " << mode << endl;
        return 1;
    }

    audio = movingAverage(audio, 2);
    normalize(audio);

    if (!writeWav(outPath, audio, fs)) {
        cerr << "Failed to write " << outPath << endl;
        return 1;
    }
    cout << "Wrote " << audio.size() << " audio samples @ " << fs << " Hz to " << outPath << endl;
    return 0;
}

void usage(const char* prog) {
    cout << "Usage:" << endl
         << "  " << prog << " test" << endl
         << "  " << prog << " mod   <am|fm|usb|lsb> <in.wav> <out.iq> [am_index|fm_deviation_hz]" << endl
         << "  " << prog << " demod <am|fm|usb|lsb> <in.iq> <out.wav> <sample_rate> [fm_deviation_hz]" << endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    string cmd = argv[1];
    if (cmd == "test") {
        return runTest();
    }
    if (cmd == "mod" && argc >= 5) {
        double param = argc >= 6 ? atof(argv[5]) : 0.0;
        return runMod(argv[2], argv[3], argv[4], param);
    }
    if (cmd == "demod" && argc >= 6) {
        double param = argc >= 7 ? atof(argv[6]) : 0.0;
        return runDemod(argv[2], argv[3], argv[4], atoi(argv[5]), param);
    }

    usage(argv[0]);
    return 1;
}
