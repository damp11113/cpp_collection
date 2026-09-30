// Simple I/Q AM/FM/SSB modulator and demodulator (complex baseband)
//
// Build: g++ -O2 -std=c++17 iq_modem.cpp -o iq_modem
//
// Usage:
//   iq_modem test
//       Run a self-test: tone -> modulate -> noise -> demodulate, prints SNR.
//   iq_modem mod   <am|fm|usb|lsb> <in.wav> <out.wav|out.iq> [options]
//       Modulate a WAV (16-bit PCM or float, mixed to mono) into I/Q.
//       out.wav -> stereo WAV, I = left, Q = right (float32, or 16-bit with --pcm16)
//       out.iq  -> raw interleaved float32 (cf32, GNU Radio / SDR# compatible)
//   iq_modem demod <am|fm|usb|lsb> <in.wav|in.iq> <out.wav> [options]
//       Demodulate I/Q into a 16-bit PCM mono WAV.
//
// Options:
//   --audio-bw <hz>    low-pass the audio (mod: before modulating, demod: after)
//   --channel-bw <hz>  RF channel filter, total occupied bandwidth
//                      (AM/FM: centered on the carrier, USB/LSB: one-sided)
//   --iq-rate <hz>     mod: I/Q output sample rate (default: audio rate)
//   --rate <hz>        demod: sample rate of a raw .iq input (WAV uses its header)
//   --audio-rate <hz>  demod: output audio rate (default: min(I/Q rate, 48000))
//   --dev <hz>         FM peak deviation (default 5000)
//   --am-index <m>     AM modulation index (default 0.8)
//   --pcm16            write I/Q WAV as 16-bit PCM instead of float32
//   --noise-floor <db> add white noise across the whole I/Q band at this total
//                      power in dBFS (0 dBFS = full-scale carrier, |I + jQ| = 1)
//   --cnr <db>         add white noise so the carrier-to-noise ratio inside the
//                      channel bandwidth (or the whole band) is this value
//   --seed <n>         noise seed for repeatable output (default: random)
//                      Noise is added after modulation (mod) or before the channel
//                      filter (demod), so both can simulate a noisy receiver.
//
// Example (narrowband FM, 3 kHz audio, 16 kHz channel, 192 kHz I/Q WAV):
//   iq_modem mod fm voice.wav nbfm_iq.wav --audio-bw 3000 --channel-bw 16000 --iq-rate 192000
//   iq_modem demod fm nbfm_iq.wav out.wav --channel-bw 16000 --audio-bw 3000

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <string>
#include <random>
#include <sstream>
#include <iomanip>
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
// FFT and fast convolution
// ---------------------------------------------------------------------------

typedef complex<double> cd;

// In-place iterative radix-2 FFT; size must be a power of two.
void fft(vector<cd>& a, bool inverse) {
    const size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * PI / len * (inverse ? 1 : -1);
        cd wlen(cos(ang), sin(ang));
        for (size_t i = 0; i < n; i += len) {
            cd w(1.0, 0.0);
            for (size_t k = 0; k < len / 2; ++k) {
                cd u = a[i + k];
                cd v = a[i + k + len / 2] * w;
                a[i + k] = u + v;
                a[i + k + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (inverse) for (auto& x : a) x /= (double)n;
}

inline cd toCd(float x) { return cd(x, 0.0); }
inline cd toCd(const iq_t& x) { return cd(x.real(), x.imag()); }
inline void addTo(float& y, const cd& v) { y += (float)v.real(); }
inline void addTo(iq_t& y, const cd& v) { y += iq_t((float)v.real(), (float)v.imag()); }

// Centered (zero-delay) FIR convolution. Works for real or complex data/taps
// (the output has the input's type). Long filters use FFT overlap-add, which
// costs O(log taps) per sample instead of O(taps).
template <typename T, typename C>
vector<T> firFilter(const vector<T>& in, const vector<C>& h) {
    const long n = (long)in.size();
    const long taps = (long)h.size();
    const long r = (taps - 1) / 2;
    vector<T> out(n);

    if (taps <= 64) {
        for (long i = 0; i < n; ++i) {
            cd acc = 0.0;
            for (long k = max(-r, i - (n - 1)); k <= min(r, i); ++k) {
                acc += toCd(in[i - k]) * toCd(h[k + r]);
            }
            addTo(out[i], acc);
        }
        return out;
    }

    size_t fftSize = 1;
    while (fftSize < 4 * (size_t)taps) fftSize <<= 1;
    const long block = (long)fftSize - taps + 1; // new input samples per FFT

    vector<cd> H(fftSize);
    for (long k = 0; k < taps; ++k) H[k] = toCd(h[k]);
    fft(H, false);

    vector<cd> buf(fftSize);
    for (long start = 0; start < n; start += block) {
        long len = min(block, n - start);
        fill(buf.begin(), buf.end(), cd());
        for (long i = 0; i < len; ++i) buf[i] = toCd(in[start + i]);

        fft(buf, false);
        for (size_t i = 0; i < fftSize; ++i) buf[i] *= H[i];
        fft(buf, true);

        // buf[j] is the causal output at start + j; the centered output is r samples earlier
        for (long j = 0; j < len + taps - 1; ++j) {
            long o = start + j - r;
            if (o >= 0 && o < n) addTo(out[o], buf[j]);
        }
    }
    return out;
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
    return firFilter(in, h);
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
// Filters and resampling
// ---------------------------------------------------------------------------

const int MAX_FIR_RADIUS = 512;

// Blackman window over x in [-1, 1]
double blackman(double x) {
    return 0.42 + 0.5 * cos(PI * x) + 0.08 * cos(2.0 * PI * x);
}

double sinc(double x) {
    return fabs(x) < 1e-12 ? 1.0 : sin(PI * x) / (PI * x);
}

// Windowed-sinc low-pass, 2*radius+1 taps, unity gain at DC.
// Transition band is ~20% of the cutoff (longer filters for sharper edges, capped).
vector<float> designLowpass(double cutoff, double fs) {
    double transition = max(cutoff * 0.2, fs * 0.002);
    int radius = min((int)ceil(2.75 * fs / transition), MAX_FIR_RADIUS);
    double fc = cutoff / fs;

    vector<float> h(2 * radius + 1);
    double sum = 0.0;
    for (int k = -radius; k <= radius; ++k) {
        double v = 2.0 * fc * sinc(2.0 * fc * k) * blackman((double)k / (radius + 1));
        h[k + radius] = (float)v;
        sum += v;
    }
    for (float& v : h) v = (float)(v / sum);
    return h;
}

// Arbitrary-ratio resampler (windowed-sinc interpolation). Band-limits to the
// lower of the two Nyquist rates so downsampling doesn't alias.
template <typename T>
vector<T> resample(const vector<T>& in, double fsIn, double fsOut) {
    if (fsIn == fsOut || in.empty()) return in;

    const double ratio = fsIn / fsOut;                          // input samples per output sample
    const double fc = 0.5 * min(1.0, fsOut / fsIn) * 0.95;      // cutoff in cycles/input sample, 5% guard
    const int half = (int)ceil(16.0 / (2.0 * fc));              // 16 zero crossings each side
    const long n = (long)in.size();

    // Kernel is symmetric, so tabulate it over |d| at 512 points per input sample
    // and interpolate linearly instead of calling sin/cos for every tap.
    const int P = 512;
    vector<float> table((size_t)half * P + 2);
    for (size_t i = 0; i < table.size(); ++i) {
        double d = (double)i / P;
        table[i] = d >= half ? 0.0f : (float)(2.0 * fc * sinc(2.0 * fc * d) * blackman(d / half));
    }

    vector<T> out((size_t)floor(n / ratio));
    for (size_t m = 0; m < out.size(); ++m) {
        double t = m * ratio;
        long c = (long)floor(t);
        T acc = T();
        for (long k = max(0L, c - half + 1); k <= min(n - 1, c + half); ++k) {
            double x = fabs(t - k) * P;
            size_t i0 = (size_t)x;
            float frac = (float)(x - i0);
            acc += in[k] * (table[i0] + frac * (table[i0 + 1] - table[i0]));
        }
        out[m] = acc;
    }
    return out;
}

// Complex channel filter. AM/FM occupy -bw/2..+bw/2; USB occupies 0..bw and
// LSB -bw..0, so the low-pass prototype is shifted to the middle of the sideband.
vector<iq_t> channelFilter(const vector<iq_t>& iq, const string& mode, double bw, double fs) {
    double center = mode == "usb" ? bw / 2 : mode == "lsb" ? -bw / 2 : 0.0;
    vector<float> lp = designLowpass(bw / 2, fs);
    const int r = (int)(lp.size() - 1) / 2;

    vector<iq_t> h(lp.size());
    for (int k = -r; k <= r; ++k) {
        h[k + r] = lp[k + r] * polar(1.0f, (float)(2.0 * PI * center * k / fs));
    }
    return firFilter(iq, h);
}

// ---------------------------------------------------------------------------
// Modulation / demodulation pipelines
// ---------------------------------------------------------------------------

struct ModemOptions {
    string mode;
    double amIndex = 0.8;
    double deviation = 5000.0; // FM peak deviation, Hz
    double audioBw = 0.0;      // audio low-pass cutoff, Hz (0 = off)
    double channelBw = 0.0;    // total occupied RF bandwidth, Hz (0 = off)
    double iqRate = 0.0;       // mod: I/Q output rate (0 = same as audio). demod: rate of a raw .iq input
    double audioRate = 0.0;    // demod: audio output rate (0 = min(I/Q rate, 48000))
    bool pcm16 = false;        // write I/Q WAV as 16-bit PCM instead of float32
    double noiseFloor = NAN;   // total noise power over the I/Q band, dBFS (NAN = off)
    double cnr = NAN;          // carrier-to-noise ratio in the channel, dB (NAN = off)
    unsigned seed = 0;         // noise seed (0 = random)
};

bool isValidMode(const string& mode) {
    return mode == "am" || mode == "fm" || mode == "usb" || mode == "lsb";
}

// Warn about settings that can't work at the given I/Q rate
void checkSettings(const ModemOptions& o, double iqRate) {
    if (o.channelBw > iqRate) {
        cerr << "Warning: channel bandwidth " << o.channelBw << " Hz is wider than the I/Q rate "
             << iqRate << " Hz, raise --iq-rate" << endl;
    }
    if (o.mode == "fm") {
        double audioBw = o.audioBw > 0 ? o.audioBw : 15000.0;
        double carson = 2.0 * (o.deviation + audioBw);
        if (carson > iqRate) {
            cerr << "Warning: FM needs about " << carson << " Hz (Carson's rule) but the I/Q rate is only "
                 << iqRate << " Hz, the signal will alias. Raise --iq-rate or lower --dev/--audio-bw" << endl;
        }
    }
}

// Fixed-point formatting for printouts (avoids "1e+06" and "-8.6e-08")
string fmt(double v, int decimals = 1) {
    ostringstream os;
    os << fixed << setprecision(decimals) << (fabs(v) < 0.5 * pow(10.0, -decimals) ? 0.0 : v);
    return os.str();
}

double meanPower(const vector<iq_t>& iq) {
    double p = 0.0;
    for (const auto& s : iq) p += norm(s);
    return iq.empty() ? 0.0 : p / iq.size();
}

double toDb(double power) {
    return 10.0 * log10(max(power, 1e-30));
}

// Adds a white noise floor across the whole I/Q band, set either as an absolute
// level (--noise-floor) or as a carrier-to-noise ratio inside the channel (--cnr).
void applyNoise(vector<iq_t>& iq, double iqRate, const ModemOptions& o) {
    bool useFloor = !isnan(o.noiseFloor);
    if (!useFloor && isnan(o.cnr)) return;

    double sigPow = meanPower(iq);
    double bw = o.channelBw > 0 && o.channelBw < iqRate ? o.channelBw : iqRate;
    double inBand = bw / iqRate; // share of the white noise that falls inside the channel
    double noisePow = useFloor ? pow(10.0, o.noiseFloor / 10.0)
                               : sigPow / pow(10.0, o.cnr / 10.0) / inBand;

    unsigned seed = o.seed ? o.seed : random_device{}();
    addNoise(iq, (float)sqrt(noisePow / 2.0), seed);

    cout << "Noise floor: " << fmt(toDb(noisePow)) << " dBFS total (" << fmt(toDb(noisePow / iqRate))
         << " dBFS/Hz), signal " << fmt(toDb(sigPow)) << " dBFS, CNR in " << fmt(bw, 0) << " Hz: "
         << fmt(toDb(sigPow / (noisePow * inBand))) << " dB, seed " << seed << endl;
}

vector<iq_t> modulatePipeline(vector<float> audio, double audioRate, const ModemOptions& o, double& iqRate) {
    if (o.audioBw > 0 && o.audioBw < audioRate / 2) {
        audio = firFilter(audio, designLowpass(o.audioBw, audioRate));
    }

    iqRate = o.iqRate > 0 ? o.iqRate : audioRate;
    audio = resample(audio, audioRate, iqRate);

    vector<iq_t> iq;
    if (o.mode == "am") {
        iq = amModulate(audio, (float)o.amIndex);
    } else if (o.mode == "fm") {
        iq = fmModulate(audio, iqRate, o.deviation);
    } else {
        iq = ssbModulate(audio, o.mode == "usb");
    }

    if (o.channelBw > 0 && o.channelBw < iqRate) {
        iq = channelFilter(iq, o.mode, o.channelBw, iqRate);
    }
    return iq;
}

vector<float> demodulatePipeline(vector<iq_t> iq, double iqRate, const ModemOptions& o, double& audioRate) {
    if (o.channelBw > 0 && o.channelBw < iqRate) {
        iq = channelFilter(iq, o.mode, o.channelBw, iqRate);
    }

    vector<float> audio;
    if (o.mode == "am") {
        audio = amDemodulate(iq);
    } else if (o.mode == "fm") {
        audio = fmDemodulate(iq, iqRate, o.deviation);
    } else {
        audio = ssbDemodulate(iq, o.mode == "usb");
    }

    if (o.audioBw > 0 && o.audioBw < iqRate / 2) {
        audio = firFilter(audio, designLowpass(o.audioBw, iqRate));
    }

    audioRate = o.audioRate > 0 ? o.audioRate : min(iqRate, 48000.0);
    return resample(audio, iqRate, audioRate);
}

// ---------------------------------------------------------------------------
// File I/O
// ---------------------------------------------------------------------------

bool hasWavExtension(const string& path) {
    if (path.size() < 4) return false;
    string ext = path.substr(path.size() - 4);
    transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext == ".wav";
}

// Reads a 16-bit PCM or 32-bit float WAV. Returns one vector per channel in [-1, 1].
bool readWav(const string& path, vector<vector<float>>& channels, int& sampleRate) {
    ifstream f(path, ios::binary);
    if (!f) return false;

    char riff[12];
    f.read(riff, 12);
    if (!f || strncmp(riff, "RIFF", 4) != 0 || strncmp(riff + 8, "WAVE", 4) != 0) return false;

    uint16_t numChannels = 0, bits = 0, format = 0;
    uint32_t rate = 0;
    // Walk chunks so files with LIST/fact chunks still load
    while (f) {
        char id[4];
        uint32_t size;
        f.read(id, 4);
        f.read((char*)&size, 4);
        if (!f) return false;

        if (strncmp(id, "fmt ", 4) == 0) {
            vector<char> buf(max(size, 16u));
            f.read(buf.data(), size);
            if (size & 1) f.seekg(1, ios::cur);
            memcpy(&format, &buf[0], 2);
            memcpy(&numChannels, &buf[2], 2);
            memcpy(&rate, &buf[4], 4);
            memcpy(&bits, &buf[14], 2);
            if (format == 0xFFFE && size >= 26) memcpy(&format, &buf[24], 2); // WAVE_FORMAT_EXTENSIBLE
        } else if (strncmp(id, "data", 4) == 0) {
            bool isPcm16 = format == 1 && bits == 16;
            bool isFloat = format == 3 && bits == 32;
            if (!(isPcm16 || isFloat) || numChannels < 1) {
                cerr << "Only 16-bit PCM or 32-bit float WAV is supported" << endl;
                return false;
            }
            const size_t bytesPerFrame = numChannels * (bits / 8);
            vector<char> raw(size);
            f.read(raw.data(), size);
            size_t frames = f.gcount() / bytesPerFrame;

            channels.assign(numChannels, vector<float>(frames));
            for (size_t i = 0; i < frames; ++i) {
                for (int c = 0; c < numChannels; ++c) {
                    const char* p = &raw[i * bytesPerFrame + c * (bits / 8)];
                    if (isPcm16) {
                        int16_t v;
                        memcpy(&v, p, 2);
                        channels[c][i] = v / 32768.0f;
                    } else {
                        memcpy(&channels[c][i], p, 4);
                    }
                }
            }
            sampleRate = (int)rate;
            return true;
        } else {
            f.seekg(size + (size & 1), ios::cur);
        }
    }
    return false;
}

// Writes channels as 16-bit PCM or 32-bit IEEE float WAV
bool writeWav(const string& path, const vector<vector<float>>& channels, int sampleRate, bool asFloat) {
    ofstream f(path, ios::binary);
    if (!f || channels.empty()) return false;

    const uint16_t numChannels = (uint16_t)channels.size();
    const uint16_t bits = asFloat ? 32 : 16;
    const uint16_t blockAlign = numChannels * bits / 8;
    const size_t frames = channels[0].size();
    const uint64_t dataSize = (uint64_t)frames * blockAlign;
    // Non-PCM formats carry a cbSize field and a fact chunk
    const uint32_t fmtSize = asFloat ? 18 : 16;
    const uint64_t riffSize = 4 + (8 + fmtSize) + (asFloat ? 12 : 0) + 8 + dataSize;
    if (riffSize > 0xFFFFFFFFull) {
        cerr << "Output too large for WAV (4 GB limit)" << endl;
        return false;
    }

    auto put16 = [&](uint16_t v) { f.write((const char*)&v, 2); };
    auto put32 = [&](uint32_t v) { f.write((const char*)&v, 4); };

    f.write("RIFF", 4);
    put32((uint32_t)riffSize);
    f.write("WAVE", 4);
    f.write("fmt ", 4);
    put32(fmtSize);
    put16(asFloat ? 3 : 1);
    put16(numChannels);
    put32((uint32_t)sampleRate);
    put32((uint32_t)sampleRate * blockAlign);
    put16(blockAlign);
    put16(bits);
    if (asFloat) {
        put16(0); // cbSize
        f.write("fact", 4);
        put32(4);
        put32((uint32_t)frames);
    }
    f.write("data", 4);
    put32((uint32_t)dataSize);

    for (size_t i = 0; i < frames; ++i) {
        for (const auto& ch : channels) {
            if (asFloat) {
                f.write((const char*)&ch[i], 4);
            } else {
                float c = max(-1.0f, min(1.0f, ch[i]));
                int16_t v = (int16_t)lrintf(c * 32767.0f);
                f.write((const char*)&v, 2);
            }
        }
    }
    return (bool)f;
}

// Reads an audio WAV and mixes it down to mono
bool readAudio(const string& path, vector<float>& audio, int& sampleRate) {
    vector<vector<float>> channels;
    if (!readWav(path, channels, sampleRate)) return false;

    audio.assign(channels[0].size(), 0.0f);
    for (const auto& ch : channels) {
        for (size_t i = 0; i < audio.size(); ++i) audio[i] += ch[i] / channels.size();
    }
    return true;
}

// Scale so the loudest sample sits at ~0.9 full scale
void normalize(vector<float>& audio) {
    float peak = 0;
    for (float s : audio) peak = max(peak, fabs(s));
    if (peak > 0) for (float& s : audio) s *= 0.9f / peak;
}

// Loads I/Q from a stereo WAV (I = left, Q = right, rate from header) or a raw
// interleaved float32 file (cf32, rate must be supplied).
bool loadIq(const string& path, vector<iq_t>& iq, double& rate) {
    if (hasWavExtension(path)) {
        vector<vector<float>> channels;
        int wavRate;
        if (!readWav(path, channels, wavRate)) return false;
        if (channels.size() != 2) {
            cerr << "I/Q WAV must have 2 channels (I = left, Q = right)" << endl;
            return false;
        }
        if (rate > 0 && rate != wavRate) {
            cerr << "Note: using the WAV header rate " << wavRate << " Hz" << endl;
        }
        rate = wavRate;
        iq.resize(channels[0].size());
        for (size_t i = 0; i < iq.size(); ++i) iq[i] = iq_t(channels[0][i], channels[1][i]);
        return true;
    }

    if (rate <= 0) {
        cerr << "Raw .iq input needs its sample rate (--rate)" << endl;
        return false;
    }
    ifstream f(path, ios::binary | ios::ate);
    if (!f) return false;
    size_t bytes = f.tellg();
    f.seekg(0);
    iq.resize(bytes / sizeof(iq_t));
    f.read((char*)iq.data(), iq.size() * sizeof(iq_t));
    return true;
}

// Saves I/Q as a stereo WAV (by extension) or raw interleaved float32
bool saveIq(const string& path, const vector<iq_t>& iq, double rate, bool pcm16) {
    if (hasWavExtension(path)) {
        vector<vector<float>> channels(2, vector<float>(iq.size()));
        float peak = 0;
        for (size_t i = 0; i < iq.size(); ++i) {
            channels[0][i] = iq[i].real();
            channels[1][i] = iq[i].imag();
            peak = max(peak, max(fabs(iq[i].real()), fabs(iq[i].imag())));
        }
        // 16-bit can't go past full scale (AM peaks at 1 + index), so scale down if needed
        if (pcm16 && peak > 0.99f) {
            for (auto& ch : channels) for (float& s : ch) s *= 0.99f / peak;
        }
        return writeWav(path, channels, (int)lround(rate), !pcm16);
    }

    ofstream f(path, ios::binary);
    if (!f) return false;
    f.write((const char*)iq.data(), iq.size() * sizeof(iq_t));
    return (bool)f;
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

    const size_t skip = 4800; // ignore filter/DC blocker settling
    const float noise = 0.02f;
    bool ok = true;

    cout << "Self-test: 1 kHz tone, fs = " << fs << " Hz, noise sigma = " << noise << endl;
    auto report = [&](const string& name, double value, double minValue) {
        bool pass = value > minValue;
        ok = ok && pass;
        cout << "  " << name << ": " << value << " dB" << (pass ? "" : "  <-- FAIL") << endl;
    };

    // AM with unknown carrier phase - envelope detector shouldn't care
    vector<iq_t> am = amModulate(tone, 0.8f);
    rotatePhase(am, 1.1f);
    addNoise(am, noise);
    report("AM demod SNR", measureSnrDb(tone, movingAverage(amDemodulate(am), 4), skip), 20.0);

    // FM
    const double dev = 5000.0;
    vector<iq_t> fm = fmModulate(tone, fs, dev);
    rotatePhase(fm, 1.1f);
    addNoise(fm, noise);
    report("FM demod SNR", measureSnrDb(tone, movingAverage(fmDemodulate(fm, fs, dev), 4), skip), 20.0);

    // SSB: no phase rotation here, the product detector needs a locked carrier
    vector<iq_t> usb = ssbModulate(tone, true);
    addNoise(usb, noise);
    report("USB demod SNR", measureSnrDb(tone, ssbDemodulate(usb, true), skip), 20.0);

    vector<iq_t> lsb = ssbModulate(tone, false);
    addNoise(lsb, noise);
    report("LSB demod SNR", measureSnrDb(tone, ssbDemodulate(lsb, false), skip), 20.0);

    // Opposite sideband rejection: demodulate the USB signal as LSB
    vector<float> wrongOut = ssbDemodulate(ssbModulate(tone, true), false);
    double rightPow = 0, wrongPow = 0;
    for (size_t i = skip; i < n - skip; ++i) {
        rightPow += (double)tone[i] * tone[i];
        wrongPow += (double)wrongOut[i] * wrongOut[i];
    }
    report("Opposite sideband rejection", 10.0 * log10(rightPow / max(wrongPow, 1e-20)), 30.0);

    // Full pipeline: 1 kHz wanted tone + 12 kHz tone that the audio filter must remove,
    // resampled to 192 kHz I/Q, with an interferer 40 kHz away that is stronger than the
    // wanted signal. Without a channel filter it takes over the AM envelope detector and
    // captures the FM discriminator.
    cout << "Pipeline (audio BW 3 kHz, I/Q rate 192 kHz, 1.5x stronger interferer at +40 kHz):" << endl;
    vector<float> mixed(n);
    for (size_t i = 0; i < n; ++i) {
        mixed[i] = tone[i] * 0.6f + 0.3f * (float)sin(2.0 * PI * 12000.0 * i / fs);
    }
    const double iqRate = 192000.0;

    struct Case { string mode; double channelBw; };
    for (const Case& c : { Case{"am", 8000}, Case{"fm", 18000}, Case{"usb", 3000}, Case{"lsb", 3000} }) {
        ModemOptions o;
        o.mode = c.mode;
        o.audioBw = 3000;
        o.iqRate = iqRate;
        o.deviation = 5000;

        double rateOut;
        vector<iq_t> sig = modulatePipeline(mixed, fs, o, rateOut);
        for (size_t i = 0; i < sig.size(); ++i) {
            sig[i] += polar(1.5f, (float)(2.0 * PI * 40000.0 * i / iqRate));
        }
        addNoise(sig, noise);

        double audioRate;
        vector<float> noFilter = demodulatePipeline(sig, iqRate, o, audioRate);
        o.channelBw = c.channelBw;
        vector<float> filtered = demodulatePipeline(sig, iqRate, o, audioRate);

        string name = c.mode;
        transform(name.begin(), name.end(), name.begin(), ::toupper);
        cout << "  " << name << " without channel filter: " << measureSnrDb(tone, noFilter, skip) << " dB" << endl;
        report(name + " with " + to_string((int)c.channelBw) + " Hz channel filter",
               measureSnrDb(tone, filtered, skip), 20.0);
    }

    // Noise floor calibration: level in dBFS and CNR inside the channel
    cout << "Noise floor:" << endl;
    auto checkNear = [&](const string& name, double value, double expected, double tol) {
        bool pass = fabs(value - expected) <= tol;
        ok = ok && pass;
        cout << "  " << name << ": " << value << " dB (expected " << expected << ")"
             << (pass ? "" : "  <-- FAIL") << endl;
    };

    ModemOptions no;
    no.mode = "fm";
    no.noiseFloor = -30.0;
    no.seed = 42;
    vector<iq_t> silence(200000);
    applyNoise(silence, iqRate, no);
    checkNear("Measured noise floor", toDb(meanPower(silence)), -30.0, 0.1);

    no.noiseFloor = NAN;
    no.cnr = 15.0;
    no.channelBw = 18000;
    vector<iq_t> clean = fmModulate(tone, iqRate, 5000);
    vector<iq_t> noisy = clean;
    applyNoise(noisy, iqRate, no);
    vector<iq_t> noiseOnly(noisy.size());
    for (size_t i = 0; i < noisy.size(); ++i) noiseOnly[i] = noisy[i] - clean[i];
    noiseOnly = channelFilter(noiseOnly, "fm", no.channelBw, iqRate);
    checkNear("Measured CNR in 18 kHz", toDb(meanPower(clean) / meanPower(noiseOnly)), 15.0, 0.5);

    cout << (ok ? "PASS" : "FAIL") << endl;
    return ok ? 0 : 1;
}

int runMod(const string& inPath, const string& outPath, const ModemOptions& o) {
    vector<float> audio;
    int fs;
    if (!readAudio(inPath, audio, fs)) {
        cerr << "Failed to read " << inPath << endl;
        return 1;
    }

    checkSettings(o, o.iqRate > 0 ? o.iqRate : fs);
    double iqRate;
    vector<iq_t> iq = modulatePipeline(audio, fs, o, iqRate);
    applyNoise(iq, iqRate, o);

    if (!saveIq(outPath, iq, iqRate, o.pcm16)) {
        cerr << "Failed to write " << outPath << endl;
        return 1;
    }
    cout << "Wrote " << iq.size() << " I/Q samples @ " << fmt(iqRate, 0) << " Hz to " << outPath
         << (hasWavExtension(outPath) ? (o.pcm16 ? " (stereo 16-bit WAV)" : " (stereo float32 WAV)") : " (raw cf32)")
         << endl;
    return 0;
}

int runDemod(const string& inPath, const string& outPath, const ModemOptions& o) {
    vector<iq_t> iq;
    double iqRate = o.iqRate;
    if (!loadIq(inPath, iq, iqRate)) {
        cerr << "Failed to read " << inPath << endl;
        return 1;
    }

    checkSettings(o, iqRate);
    applyNoise(iq, iqRate, o);
    double audioRate;
    vector<float> audio = demodulatePipeline(iq, iqRate, o, audioRate);
    normalize(audio);

    if (!writeWav(outPath, { audio }, (int)lround(audioRate), false)) {
        cerr << "Failed to write " << outPath << endl;
        return 1;
    }
    cout << "Wrote " << audio.size() << " audio samples @ " << fmt(audioRate, 0) << " Hz to " << outPath << endl;
    return 0;
}

void usage(const char* prog) {
    cout << "Usage:" << endl
         << "  " << prog << " test" << endl
         << "  " << prog << " mod   <am|fm|usb|lsb> <in.wav> <out.wav|out.iq> [options]" << endl
         << "  " << prog << " demod <am|fm|usb|lsb> <in.wav|in.iq> <out.wav> [options]" << endl
         << endl
         << "Options:" << endl
         << "  --audio-bw <hz>    low-pass the audio (mod: before modulating, demod: after)" << endl
         << "  --channel-bw <hz>  RF channel filter, total occupied bandwidth" << endl
         << "  --iq-rate <hz>     mod: I/Q output sample rate (default: audio rate)" << endl
         << "  --rate <hz>        demod: sample rate of a raw .iq input (WAV uses its header)" << endl
         << "  --audio-rate <hz>  demod: output audio rate (default: min(I/Q rate, 48000))" << endl
         << "  --dev <hz>         FM peak deviation (default 5000)" << endl
         << "  --am-index <m>     AM modulation index (default 0.8)" << endl
         << "  --pcm16            write I/Q WAV as 16-bit PCM instead of float32" << endl
         << "  --noise-floor <db> add white noise at this total power in dBFS over the I/Q band" << endl
         << "  --cnr <db>         add white noise for this carrier-to-noise ratio in the channel" << endl
         << "  --seed <n>         noise seed for repeatable output (default: random)" << endl
         << endl
         << "I/Q WAV files are stereo with I = left, Q = right." << endl;
}

// Parses options after the three positional arguments. Also accepts the older
// positional forms: "mod ... [am_index|fm_dev]" and "demod ... <rate> [fm_dev]".
bool parseOptions(int argc, char* argv[], int start, bool isMod, ModemOptions& o) {
    int positional = 0;
    for (int i = start; i < argc; ++i) {
        string a = argv[i];
        if (a == "--pcm16") {
            o.pcm16 = true;
            continue;
        }
        if (a.rfind("--", 0) == 0) {
            if (i + 1 >= argc) {
                cerr << "Missing value for " << a << endl;
                return false;
            }
            double v = atof(argv[++i]);
            if (a == "--audio-bw") o.audioBw = v;
            else if (a == "--channel-bw") o.channelBw = v;
            else if (a == "--iq-rate" || a == "--rate") o.iqRate = v;
            else if (a == "--audio-rate") o.audioRate = v;
            else if (a == "--dev") o.deviation = v;
            else if (a == "--am-index") o.amIndex = v;
            else if (a == "--noise-floor") o.noiseFloor = v;
            else if (a == "--cnr") o.cnr = v;
            else if (a == "--seed") o.seed = (unsigned)v;
            else {
                cerr << "Unknown option " << a << endl;
                return false;
            }
            continue;
        }

        double v = atof(a.c_str());
        if (isMod && positional == 0) {
            if (o.mode == "am") o.amIndex = v;
            else if (o.mode == "fm") o.deviation = v;
        } else if (!isMod && positional == 0) {
            o.iqRate = v;
        } else if (!isMod && positional == 1) {
            o.deviation = v;
        } else {
            cerr << "Unexpected argument " << a << endl;
            return false;
        }
        ++positional;
    }
    return true;
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
    if ((cmd == "mod" || cmd == "demod") && argc >= 5) {
        ModemOptions o;
        o.mode = argv[2];
        if (!isValidMode(o.mode)) {
            cerr << "Unknown mode: " << o.mode << endl;
            return 1;
        }
        if (!parseOptions(argc, argv, 5, cmd == "mod", o)) return 1;
        if (!isnan(o.noiseFloor) && !isnan(o.cnr)) {
            cerr << "Use either --noise-floor or --cnr, not both" << endl;
            return 1;
        }
        return cmd == "mod" ? runMod(argv[3], argv[4], o) : runDemod(argv[3], argv[4], o);
    }

    usage(argv[0]);
    return 1;
}
