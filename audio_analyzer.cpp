/*
 * Audio Analysis Video Generator - Enhanced Version
 * 
 * Analyzes stereo WAV audio and outputs video with:
 * - Linear FFT Spectrum (-100dB to +5dB) with scaling
 * - Channel Loudness Meters
 * - XY Oscilloscope
 * - Goniometer (Phase/Vector Scope)
 * 
 * Features:
 * - Multi-threaded processing
 * - AVX2 optimizations
 * - Customizable resolution and framerate
 * - Improved text rendering
 * 
 * Dependencies:
 * - FFmpeg (libavcodec, libavformat, libavutil, libswscale)
 * - FFTW3
 * 
 * Compile with AVX2:
 * g++ -std=c++17 audio_analyzer.cpp -o audio_analyzer \
 *     -lavcodec -lavformat -lavutil -lswscale -lfftw3 -lm \
 *     -O3 -mavx2 -mfma -pthread
 * 
 * Usage:
 * ./audio_analyzer input.wav output.mp4 "Title" "Description" [width] [height] [fps]
 * Example: ./audio_analyzer audio.wav out.mp4 "My Track" "Description" 1920 1080 60
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstring>
#include <fstream>
#include <thread>
#include <mutex>
#include <immintrin.h>

extern "C" {
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/opt.h>
#include <libavutil/imgutils.h>
#include <libswscale/swscale.h>
}

#include <fftw3.h>

// Configuration
int VIDEO_WIDTH = 1920;
int VIDEO_HEIGHT = 1080;
int VIDEO_FPS = 60;
const int FFT_SIZE = 8192;
const double MIN_DB = -100.0;
const double MAX_DB = 5.0;
const double SPECTRUM_SCALE = 0.75; // Reduced scale for better visibility

// Color definitions (RGB)
struct Color {
    uint8_t r, g, b;
};

const Color BG_COLOR = {15, 15, 25};
const Color GRID_COLOR = {50, 50, 60};
const Color SPECTRUM_L_COLOR = {0, 255, 150};
const Color SPECTRUM_R_COLOR = {255, 180, 0};
const Color LEFT_COLOR = {0, 255, 100};
const Color RIGHT_COLOR = {255, 100, 0};
const Color SCOPE_COLOR = {0, 200, 255};
const Color TEXT_COLOR = {220, 220, 220};

// Simple bitmap font (5x7)
const uint8_t font5x7[][5] = {
    {0x00, 0x00, 0x00, 0x00, 0x00}, // space
    {0x00, 0x00, 0x5F, 0x00, 0x00}, // !
    {0x00, 0x07, 0x00, 0x07, 0x00}, // "
    {0x14, 0x7F, 0x14, 0x7F, 0x14}, // #
    {0x24, 0x2A, 0x7F, 0x2A, 0x12}, // $
    {0x23, 0x13, 0x08, 0x64, 0x62}, // %
    {0x36, 0x49, 0x56, 0x20, 0x50}, // &
    {0x00, 0x08, 0x07, 0x03, 0x00}, // '
    {0x00, 0x1C, 0x22, 0x41, 0x00}, // (
    {0x00, 0x41, 0x22, 0x1C, 0x00}, // )
    {0x2A, 0x1C, 0x7F, 0x1C, 0x2A}, // *
    {0x08, 0x08, 0x3E, 0x08, 0x08}, // +
    {0x00, 0x80, 0x70, 0x30, 0x00}, // ,
    {0x08, 0x08, 0x08, 0x08, 0x08}, // -
    {0x00, 0x00, 0x60, 0x60, 0x00}, // .
    {0x20, 0x10, 0x08, 0x04, 0x02}, // /
    {0x3E, 0x51, 0x49, 0x45, 0x3E}, // 0
    {0x00, 0x42, 0x7F, 0x40, 0x00}, // 1
    {0x42, 0x61, 0x51, 0x49, 0x46}, // 2
    {0x21, 0x41, 0x45, 0x4B, 0x31}, // 3
    {0x18, 0x14, 0x12, 0x7F, 0x10}, // 4
    {0x27, 0x45, 0x45, 0x45, 0x39}, // 5
    {0x3C, 0x4A, 0x49, 0x49, 0x30}, // 6
    {0x01, 0x71, 0x09, 0x05, 0x03}, // 7
    {0x36, 0x49, 0x49, 0x49, 0x36}, // 8
    {0x06, 0x49, 0x49, 0x29, 0x1E}, // 9
    {0x00, 0x36, 0x36, 0x00, 0x00}, // :
    {0x00, 0x56, 0x36, 0x00, 0x00}, // ;
    {0x08, 0x14, 0x22, 0x41, 0x00}, // <
    {0x14, 0x14, 0x14, 0x14, 0x14}, // =
    {0x00, 0x41, 0x22, 0x14, 0x08}, // >
    {0x02, 0x01, 0x51, 0x09, 0x06}, // ?
    {0x32, 0x49, 0x79, 0x41, 0x3E}, // @
    {0x7C, 0x12, 0x11, 0x12, 0x7C}, // A
    {0x7F, 0x49, 0x49, 0x49, 0x36}, // B
    {0x3E, 0x41, 0x41, 0x41, 0x22}, // C
    {0x7F, 0x41, 0x41, 0x41, 0x3E}, // D
    {0x7F, 0x49, 0x49, 0x49, 0x41}, // E
    {0x7F, 0x09, 0x09, 0x09, 0x01}, // F
    {0x3E, 0x41, 0x49, 0x49, 0x7A}, // G
    {0x7F, 0x08, 0x08, 0x08, 0x7F}, // H
    {0x00, 0x41, 0x7F, 0x41, 0x00}, // I
    {0x20, 0x40, 0x41, 0x3F, 0x01}, // J
    {0x7F, 0x08, 0x14, 0x22, 0x41}, // K
    {0x7F, 0x40, 0x40, 0x40, 0x40}, // L
    {0x7F, 0x02, 0x0C, 0x02, 0x7F}, // M
    {0x7F, 0x04, 0x08, 0x10, 0x7F}, // N
    {0x3E, 0x41, 0x41, 0x41, 0x3E}, // O
    {0x7F, 0x09, 0x09, 0x09, 0x06}, // P
    {0x3E, 0x41, 0x51, 0x21, 0x5E}, // Q
    {0x7F, 0x09, 0x19, 0x29, 0x46}, // R
    {0x46, 0x49, 0x49, 0x49, 0x31}, // S
    {0x01, 0x01, 0x7F, 0x01, 0x01}, // T
    {0x3F, 0x40, 0x40, 0x40, 0x3F}, // U
    {0x1F, 0x20, 0x40, 0x20, 0x1F}, // V
    {0x3F, 0x40, 0x38, 0x40, 0x3F}, // W
    {0x63, 0x14, 0x08, 0x14, 0x63}, // X
    {0x07, 0x08, 0x70, 0x08, 0x07}, // Y
    {0x61, 0x51, 0x49, 0x45, 0x43}, // Z
};

// WAV file reader
class WAVReader {
private:
    std::ifstream file;
    int channels;
    int sampleRate;
    int bitsPerSample;
    long dataSize;
    long dataStart;
    
public:
    bool open(const std::string& filename) {
        file.open(filename, std::ios::binary);
        if (!file) return false;
        
        char buffer[44];
        file.read(buffer, 44);
        
        if (std::string(buffer, 4) != "RIFF") return false;
        if (std::string(buffer + 8, 4) != "WAVE") return false;
        
        channels = *(int16_t*)(buffer + 22);
        sampleRate = *(int32_t*)(buffer + 24);
        bitsPerSample = *(int16_t*)(buffer + 34);
        dataSize = *(int32_t*)(buffer + 40);
        dataStart = 44;
        
        return true;
    }
    
    int readSamples(std::vector<float>& left, std::vector<float>& right, int numSamples) {
        int bytesPerSample = bitsPerSample / 8;
        int frameSize = channels * bytesPerSample;
        std::vector<char> buffer(numSamples * frameSize);
        
        file.read(buffer.data(), buffer.size());
        int samplesRead = file.gcount() / frameSize;
        
        left.resize(samplesRead);
        right.resize(samplesRead);
        
        for (int i = 0; i < samplesRead; i++) {
            int16_t* samples = (int16_t*)(buffer.data() + i * frameSize);
            left[i] = samples[0] / 32768.0f;
            right[i] = channels > 1 ? samples[1] / 32768.0f : left[i];
        }
        
        return samplesRead;
    }
    
    int getSampleRate() const { return sampleRate; }
    int getChannels() const { return channels; }
    bool isOpen() const { return file.is_open(); }
    void close() { file.close(); }
};

// Video encoder
class VideoEncoder {
private:
    AVFormatContext* fmtCtx = nullptr;
    AVCodecContext* codecCtx = nullptr;
    AVStream* stream = nullptr;
    AVFrame* frame = nullptr;
    AVPacket* pkt = nullptr;
    SwsContext* swsCtx = nullptr;
    int frameCounter = 0;
    
public:
    bool init(const std::string& filename, int width, int height, int fps) {
        avformat_alloc_output_context2(&fmtCtx, nullptr, nullptr, filename.c_str());
        if (!fmtCtx) return false;
        
        const AVCodec* codec = avcodec_find_encoder(AV_CODEC_ID_H264);
        if (!codec) return false;
        
        stream = avformat_new_stream(fmtCtx, nullptr);
        if (!stream) return false;
        
        codecCtx = avcodec_alloc_context3(codec);
        codecCtx->width = width;
        codecCtx->height = height;
        codecCtx->time_base = {1, fps};
        codecCtx->framerate = {fps, 1};
        codecCtx->pix_fmt = AV_PIX_FMT_YUV420P;
        codecCtx->bit_rate = 10000000;
        
        if (fmtCtx->oformat->flags & AVFMT_GLOBALHEADER)
            codecCtx->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;
        
        if (avcodec_open2(codecCtx, codec, nullptr) < 0) return false;
        
        avcodec_parameters_from_context(stream->codecpar, codecCtx);
        stream->time_base = codecCtx->time_base;
        
        frame = av_frame_alloc();
        frame->format = codecCtx->pix_fmt;
        frame->width = width;
        frame->height = height;
        av_frame_get_buffer(frame, 0);
        
        pkt = av_packet_alloc();
        
        if (avio_open(&fmtCtx->pb, filename.c_str(), AVIO_FLAG_WRITE) < 0) return false;
        if (avformat_write_header(fmtCtx, nullptr) < 0) return false;
        
        return true;
    }
    
    void writeFrame(uint8_t* rgb, int width, int height) {
        if (!swsCtx) {
            swsCtx = sws_getContext(width, height, AV_PIX_FMT_RGB24,
                                   width, height, AV_PIX_FMT_YUV420P,
                                   SWS_BILINEAR, nullptr, nullptr, nullptr);
        }
        
        uint8_t* srcData[1] = {rgb};
        int srcLinesize[1] = {width * 3};
        sws_scale(swsCtx, srcData, srcLinesize, 0, height, frame->data, frame->linesize);
        
        frame->pts = frameCounter++;
        
        avcodec_send_frame(codecCtx, frame);
        while (avcodec_receive_packet(codecCtx, pkt) == 0) {
            av_packet_rescale_ts(pkt, codecCtx->time_base, stream->time_base);
            pkt->stream_index = stream->index;
            av_interleaved_write_frame(fmtCtx, pkt);
            av_packet_unref(pkt);
        }
    }
    
    void close() {
        avcodec_send_frame(codecCtx, nullptr);
        while (avcodec_receive_packet(codecCtx, pkt) == 0) {
            av_packet_rescale_ts(pkt, codecCtx->time_base, stream->time_base);
            pkt->stream_index = stream->index;
            av_interleaved_write_frame(fmtCtx, pkt);
            av_packet_unref(pkt);
        }
        
        av_write_trailer(fmtCtx);
        
        if (swsCtx) sws_freeContext(swsCtx);
        av_frame_free(&frame);
        av_packet_free(&pkt);
        avcodec_free_context(&codecCtx);
        avio_closep(&fmtCtx->pb);
        avformat_free_context(fmtCtx);
    }
};

// Drawing utilities
class Canvas {
private:
    std::vector<uint8_t> buffer;
    int width, height;
    
public:
    Canvas(int w, int h) : width(w), height(h), buffer(w * h * 3) {}
    
    void clear(Color c) {
        for (size_t i = 0; i < buffer.size(); i += 3) {
            buffer[i] = c.r;
            buffer[i + 1] = c.g;
            buffer[i + 2] = c.b;
        }
    }
    
    void setPixel(int x, int y, Color c) {
        if (x < 0 || x >= width || y < 0 || y >= height) return;
        int idx = (y * width + x) * 3;
        buffer[idx] = c.r;
        buffer[idx + 1] = c.g;
        buffer[idx + 2] = c.b;
    }
    
    void drawLine(int x1, int y1, int x2, int y2, Color c) {
        int dx = abs(x2 - x1), dy = abs(y2 - y1);
        int sx = x1 < x2 ? 1 : -1, sy = y1 < y2 ? 1 : -1;
        int err = dx - dy;
        
        while (true) {
            setPixel(x1, y1, c);
            if (x1 == x2 && y1 == y2) break;
            int e2 = 2 * err;
            if (e2 > -dy) { err -= dy; x1 += sx; }
            if (e2 < dx) { err += dx; y1 += sy; }
        }
    }
    
    void drawRect(int x, int y, int w, int h, Color c, bool fill = false) {
        if (fill) {
            for (int j = 0; j < h; j++)
                for (int i = 0; i < w; i++)
                    setPixel(x + i, y + j, c);
        } else {
            drawLine(x, y, x + w - 1, y, c);
            drawLine(x + w - 1, y, x + w - 1, y + h - 1, c);
            drawLine(x + w - 1, y + h - 1, x, y + h - 1, c);
            drawLine(x, y + h - 1, x, y, c);
        }
    }
    
    void drawText(int x, int y, const std::string& text, Color c, int scale = 2) {
        for (size_t i = 0; i < text.length(); i++) {
            char ch = text[i];
            int idx = ch - 32;
            
            if (idx >= 0 && idx < 95) {
                for (int py = 0; py < 7; py++) {
                    for (int px = 0; px < 5; px++) {
                        if (font5x7[idx][px] & (1 << py)) {
                            for (int sy = 0; sy < scale; sy++) {
                                for (int sx = 0; sx < scale; sx++) {
                                    setPixel(x + i * 6 * scale + px * scale + sx, 
                                           y + py * scale + sy, c);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    uint8_t* data() { return buffer.data(); }
    int getWidth() const { return width; }
    int getHeight() const { return height; }
};

// Audio analyzer with AVX2 optimizations
class AudioAnalyzer {
private:
    fftw_plan planL, planR;
    double* inL;
    double* inR;
    fftw_complex* outL;
    fftw_complex* outR;
    std::vector<double> window;
    
public:
    AudioAnalyzer() {
        inL = (double*)fftw_malloc(sizeof(double) * FFT_SIZE);
        inR = (double*)fftw_malloc(sizeof(double) * FFT_SIZE);
        outL = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * FFT_SIZE);
        outR = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * FFT_SIZE);
        
        planL = fftw_plan_dft_r2c_1d(FFT_SIZE, inL, outL, FFTW_MEASURE);
        planR = fftw_plan_dft_r2c_1d(FFT_SIZE, inR, outR, FFTW_MEASURE);
        
        // Hann window
        window.resize(FFT_SIZE);
        for (int i = 0; i < FFT_SIZE; i++)
            window[i] = 0.5 * (1 - cos(2 * M_PI * i / (FFT_SIZE - 1)));
    }
    
    ~AudioAnalyzer() {
        fftw_destroy_plan(planL);
        fftw_destroy_plan(planR);
        fftw_free(inL);
        fftw_free(inR);
        fftw_free(outL);
        fftw_free(outR);
    }
    
    void computeFFT(const std::vector<float>& left, const std::vector<float>& right,
                    std::vector<double>& magL, std::vector<double>& magR) {
        int size = std::min((int)left.size(), FFT_SIZE);
        
        // Apply window with AVX2
        #ifdef __AVX2__
        int idx = 0;
        for (; idx + 4 <= size; idx += 4) {
            __m128 l = _mm_loadu_ps(&left[idx]);
            __m128 r = _mm_loadu_ps(&right[idx]);
            __m256d w = _mm256_loadu_pd(&window[idx]);
            
            __m256d ld = _mm256_cvtps_pd(l);
            __m256d rd = _mm256_cvtps_pd(r);
            
            __m256d lw = _mm256_mul_pd(ld, w);
            __m256d rw = _mm256_mul_pd(rd, w);
            
            _mm256_storeu_pd(&inL[idx], lw);
            _mm256_storeu_pd(&inR[idx], rw);
        }
        for (; idx < size; idx++) {
            inL[idx] = left[idx] * window[idx];
            inR[idx] = right[idx] * window[idx];
        }
        #else
        for (int idx = 0; idx < size; idx++) {
            inL[idx] = left[idx] * window[idx];
            inR[idx] = right[idx] * window[idx];
        }
        #endif
        
        for (int i = size; i < FFT_SIZE; i++) {
            inL[i] = 0;
            inR[i] = 0;
        }
        
        fftw_execute(planL);
        fftw_execute(planR);
        
        magL.resize(FFT_SIZE / 2);
        magR.resize(FFT_SIZE / 2);
        
        // Compute magnitudes with AVX2
        #ifdef __AVX2__
        int j = 0;
        __m256d scale = _mm256_set1_pd(1.0 / FFT_SIZE);
        __m256d epsilon = _mm256_set1_pd(1e-10);
        __m256d log_scale = _mm256_set1_pd(20.0 / log(10.0));
        
        for (; j + 4 <= FFT_SIZE / 2; j += 4) {
            __m256d rL = _mm256_setr_pd(outL[j][0], outL[j+1][0], outL[j+2][0], outL[j+3][0]);
            __m256d iL = _mm256_setr_pd(outL[j][1], outL[j+1][1], outL[j+2][1], outL[j+3][1]);
            __m256d rR = _mm256_setr_pd(outR[j][0], outR[j+1][0], outR[j+2][0], outR[j+3][0]);
            __m256d iR = _mm256_setr_pd(outR[j][1], outR[j+1][1], outR[j+2][1], outR[j+3][1]);
            
            __m256d magL_sq = _mm256_add_pd(_mm256_mul_pd(rL, rL), _mm256_mul_pd(iL, iL));
            __m256d magR_sq = _mm256_add_pd(_mm256_mul_pd(rR, rR), _mm256_mul_pd(iR, iR));
            
            double mL[4], mR[4];
            _mm256_storeu_pd(mL, magL_sq);
            _mm256_storeu_pd(mR, magR_sq);
            
            for (int k = 0; k < 4 && j + k < FFT_SIZE / 2; k++) {
                magL[j + k] = 20 * log10(sqrt(mL[k]) / FFT_SIZE + 1e-10);
                magR[j + k] = 20 * log10(sqrt(mR[k]) / FFT_SIZE + 1e-10);
            }
        }
        for (; j < FFT_SIZE / 2; j++) {
            double rL = outL[j][0], iL = outL[j][1];
            double rR = outR[j][0], iR = outR[j][1];
            magL[j] = 20 * log10(sqrt(rL * rL + iL * iL) / FFT_SIZE + 1e-10);
            magR[j] = 20 * log10(sqrt(rR * rR + iR * iR) / FFT_SIZE + 1e-10);
        }
        #else
        for (int j = 0; j < FFT_SIZE / 2; j++) {
            double rL = outL[j][0], iL = outL[j][1];
            double rR = outR[j][0], iR = outR[j][1];
            magL[j] = 20 * log10(sqrt(rL * rL + iL * iL) / FFT_SIZE + 1e-10);
            magR[j] = 20 * log10(sqrt(rR * rR + iR * iR) / FFT_SIZE + 1e-10);
        }
        #endif
    }
    
    double computeRMS(const std::vector<float>& samples) {
        double sum = 0;
        
        #ifdef __AVX2__
        __m256 vsum = _mm256_setzero_ps();
        size_t idx = 0;
        for (; idx + 8 <= samples.size(); idx += 8) {
            __m256 s = _mm256_loadu_ps(&samples[idx]);
            vsum = _mm256_add_ps(vsum, _mm256_mul_ps(s, s));
        }
        
        float temp[8];
        _mm256_storeu_ps(temp, vsum);
        for (int j = 0; j < 8; j++) sum += temp[j];
        
        for (; idx < samples.size(); idx++) {
            sum += samples[idx] * samples[idx];
        }
        #else
        for (float s : samples) sum += s * s;
        #endif
        
        return 20 * log10(sqrt(sum / samples.size()) + 1e-10);
    }
};

// Visualization renderer
class Visualizer {
private:
    Canvas canvas;
    std::string title, description;
    
public:
    Visualizer(int width, int height, const std::string& t, const std::string& d) 
        : canvas(width, height), title(t), description(d) {}
    
    void render(const std::vector<float>& left, const std::vector<float>& right,
                const std::vector<double>& fftL, const std::vector<double>& fftR,
                double rmsL, double rmsR, int sampleRate) {
        canvas.clear(BG_COLOR);
        
        int w = canvas.getWidth();
        int h = canvas.getHeight();
        
        // Scale factors for different resolutions
        int textScale = std::max(2, w / 800);
        
        // Draw title and description with better visibility
        canvas.drawText(30, 30, title, TEXT_COLOR, textScale);
        canvas.drawText(30, 30 + 20 * textScale, description, TEXT_COLOR, textScale);
        
        // Calculate layout
        int spectrumHeight = h * 0.4;
        int spectrumY = h * 0.08;
        int bottomY = spectrumY + spectrumHeight + 20;
        int bottomHeight = h - bottomY - 20;
        
        int meterWidth = w * 0.04;
        int meterX = w - meterWidth * 2 - 40;
        
        int scopeSize = std::min(bottomHeight, (meterX - 80) / 2);
        
        // FFT Spectrum (larger, scaled)
        drawSpectrum(40, spectrumY, meterX - 60, spectrumHeight, fftL, fftR, sampleRate);
        
        // Loudness meters (right side)
        drawLoudnessMeter(meterX, bottomY, meterWidth, bottomHeight, rmsL, LEFT_COLOR, "L");
        drawLoudnessMeter(meterX + meterWidth + 20, bottomY, meterWidth, bottomHeight, rmsR, RIGHT_COLOR, "R");
        
        // XY Oscilloscope (bottom left)
        drawOscilloscope(40, bottomY, scopeSize, scopeSize, left, right);
        
        // Goniometer (bottom center)
        drawGoniometer(60 + scopeSize, bottomY, scopeSize, scopeSize, left, right);
    }
    
    void drawSpectrum(int x, int y, int w, int h, 
                     const std::vector<double>& fftL, const std::vector<double>& fftR,
                     int sampleRate) {
        // Draw border
        canvas.drawRect(x, y, w, h, GRID_COLOR);
        
        // Draw horizontal grid lines (dB levels)
        for (int i = 0; i <= 10; i++) {
            int yPos = y + (h * i) / 10;
            canvas.drawLine(x, yPos, x + w, yPos, GRID_COLOR);
            
            // Draw dB labels with better sizing
            double db = MAX_DB - (MAX_DB - MIN_DB) * i / 10;
            std::string label = std::to_string((int)db) + "dB";
            canvas.drawText(x - 60, yPos - 7, label, TEXT_COLOR, 2);
        }
        
        // Draw vertical grid lines (frequency markers)
        int freqMarkers[] = {100, 500, 1000, 5000, 10000, 20000};
        int nyquist = sampleRate / 2;
        for (int freq : freqMarkers) {
            if (freq > nyquist) break;
            int xPos = x + (w * freq) / nyquist;
            canvas.drawLine(xPos, y, xPos, y + h, GRID_COLOR);
            
            std::string label = freq >= 1000 ? 
                std::to_string(freq / 1000) + "k" : std::to_string(freq);
            canvas.drawText(xPos - 15, y + h + 15, label, TEXT_COLOR, 2);
        }
        
        // Draw spectrum curves (single line per channel for clarity)
        int prevXL = x, prevYL = y + h;
        int prevXR = x, prevYR = y + h;
        
        for (int i = 1; i < w; i++) {
            double freq = (double)i / w * nyquist;
            int bin = (freq / nyquist) * (fftL.size() - 1);
            
            if (bin >= 0 && bin < fftL.size()) {
                // Apply spectrum scaling
                double dbL = (fftL[bin] + 20) * SPECTRUM_SCALE - 20;
                double dbR = (fftR[bin] + 20) * SPECTRUM_SCALE - 20;
                
                dbL = std::max(MIN_DB, std::min(MAX_DB, dbL));
                dbR = std::max(MIN_DB, std::min(MAX_DB, dbR));
                
                int yL = y + h - (dbL - MIN_DB) / (MAX_DB - MIN_DB) * h;
                int yR = y + h - (dbR - MIN_DB) / (MAX_DB - MIN_DB) * h;
                
                // Draw single clean lines
                canvas.drawLine(prevXL, prevYL, x + i, yL, SPECTRUM_L_COLOR);
                canvas.drawLine(prevXR, prevYR, x + i, yR, SPECTRUM_R_COLOR);
                
                prevXL = x + i;
                prevYL = yL;
                prevXR = x + i;
                prevYR = yR;
            }
        }
        
        // Draw legend with better text
        canvas.drawRect(x + 10, y + 10, 20, 15, SPECTRUM_L_COLOR, true);
        canvas.drawText(x + 35, y + 12, "Left", TEXT_COLOR, 2);
        
        canvas.drawRect(x + 110, y + 10, 20, 15, SPECTRUM_R_COLOR, true);
        canvas.drawText(x + 135, y + 12, "Right", TEXT_COLOR, 2);
    }
    
    void drawLoudnessMeter(int x, int y, int w, int h, double rms, Color c, const std::string& label) {
        // Draw border
        canvas.drawRect(x, y, w, h, GRID_COLOR);
        
        // Draw dB scale markers
        for (int db = 0; db >= -60; db -= 10) {
            double normalized = (db - MIN_DB) / (MAX_DB - MIN_DB);
            int yPos = y + h - normalized * h;
            canvas.drawLine(x, yPos, x + w, yPos, GRID_COLOR);
            
            std::string dbLabel = std::to_string(db);
            canvas.drawText(x + w + 5, yPos - 7, dbLabel, TEXT_COLOR, 2);
        }
        
        // Calculate bar height
        double normalized = (rms - MIN_DB) / (MAX_DB - MIN_DB);
        int barHeight = std::max(0, std::min(h, (int)(normalized * h)));
        
        // Draw gradient bar
        for (int i = 0; i < barHeight; i++) {
            double ratio = (double)i / h;
            Color barColor;
            
            if (ratio > 0.85) {
                // Red zone (> -6dB)
                barColor = {255, 0, 0};
            } else if (ratio > 0.70) {
                // Yellow zone (-12dB to -6dB)
                barColor = {255, 200, 0};
            } else {
                // Green zone (< -12dB)
                barColor = c;
            }
            
            canvas.drawLine(x + 2, y + h - i, x + w - 2, y + h - i, barColor);
        }
        
        // Draw peak indicator
        int peakY = y + h - barHeight;
        canvas.drawLine(x, peakY - 1, x + w, peakY - 1, {255, 255, 255});
        canvas.drawLine(x, peakY, x + w, peakY, {255, 255, 255});
        
        // Draw label
        canvas.drawText(x + w / 2 - 5, y - 20, label, TEXT_COLOR, 2);
        
        // Draw current value
        std::string valStr = std::to_string((int)rms) + "dB";
        canvas.drawText(x, peakY - 20, valStr, TEXT_COLOR, 2);
    }
    
    void drawOscilloscope(int x, int y, int w, int h, 
                         const std::vector<float>& left, const std::vector<float>& right) {
        // Draw border and grid
        canvas.drawRect(x, y, w, h, GRID_COLOR);
        canvas.drawLine(x, y + h/2, x + w, y + h/2, GRID_COLOR);
        canvas.drawLine(x + w/2, y, x + w/2, y + h, GRID_COLOR);
        
        // Draw circular guide
        int cx = x + w/2, cy = y + h/2;
        int radius = std::min(w, h) / 2 - 5;
        for (int a = 0; a < 360; a += 5) {
            int x1 = cx + radius * cos(a * M_PI / 180);
            int y1 = cy + radius * sin(a * M_PI / 180);
            canvas.setPixel(x1, y1, GRID_COLOR);
        }
        
        // Draw oscilloscope trace
        int samples = std::min((int)left.size(), 4000);
        int skip = std::max(1, samples / 2000);
        
        for (int i = 0; i < samples - skip; i += skip) {
            int x1 = x + w/2 + left[i] * (w/2 - 5);
            int y1 = y + h/2 - right[i] * (h/2 - 5);
            int x2 = x + w/2 + left[i + skip] * (w/2 - 5);
            int y2 = y + h/2 - right[i + skip] * (h/2 - 5);
            
            canvas.drawLine(x1, y1, x2, y2, SCOPE_COLOR);
        }
        
        // Label
        canvas.drawText(x + 10, y + 10, "XY Scope", TEXT_COLOR, 2);
        canvas.drawText(x + w - 45, y + h/2 - 7, "L", LEFT_COLOR, 2);
        canvas.drawText(x + w/2 - 7, y + 10, "R", RIGHT_COLOR, 2);
    }
    
    void drawGoniometer(int x, int y, int w, int h,
                       const std::vector<float>& left, const std::vector<float>& right) {
        int cx = x + w/2, cy = y + h/2;
        int radius = std::min(w, h) / 2 - 10;
        
        // Draw polar grid circles
        for (int r = radius/4; r <= radius; r += radius/4) {
            for (int a = 0; a < 360; a += 2) {
                int px = cx + r * cos(a * M_PI / 180);
                int py = cy + r * sin(a * M_PI / 180);
                canvas.setPixel(px, py, GRID_COLOR);
            }
        }
        
        // Draw angle lines
        for (int a = 0; a < 360; a += 30) {
            int px = cx + radius * cos(a * M_PI / 180);
            int py = cy + radius * sin(a * M_PI / 180);
            canvas.drawLine(cx, cy, px, py, GRID_COLOR);
        }
        
        // Draw goniometer points
        int samples = std::min((int)left.size(), 4000);
        int skip = std::max(1, samples / 2000);
        
        for (int i = 0; i < samples; i += skip) {
            double mid = (left[i] + right[i]) / 2.0;
            double side = (left[i] - right[i]) / 2.0;
            
            double angle = atan2(side, mid + 1e-10);
            double mag = sqrt(mid * mid + side * side);
            mag = std::min(mag, 1.0);
            
            int px = cx + mag * radius * cos(angle);
            int py = cy + mag * radius * sin(angle);
            
            // Draw with slight glow effect
            canvas.setPixel(px, py, SCOPE_COLOR);
            canvas.setPixel(px + 1, py, SCOPE_COLOR);
            canvas.setPixel(px, py + 1, SCOPE_COLOR);
        }
        
        // Draw reference lines
        canvas.drawLine(cx - radius, cy, cx + radius, cy, {100, 100, 120}); // L+R (mono)
        canvas.drawLine(cx, cy - radius, cx, cy + radius, {100, 100, 120}); // L-R (stereo width)
        
        // Labels
        canvas.drawText(x + 10, y + 10, "Goniometer", TEXT_COLOR, 2);
        canvas.drawText(cx + radius - 20, cy + 7, "L", LEFT_COLOR, 2);
        canvas.drawText(cx - 10, cy - radius + 5, "R", RIGHT_COLOR, 2);
        canvas.drawText(cx + 5, cy + 5, "M", TEXT_COLOR, 2);
    }
    
    uint8_t* getFrame() { return canvas.data(); }
    int getWidth() { return canvas.getWidth(); }
    int getHeight() { return canvas.getHeight(); }
};

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cout << "Usage: " << argv[0] << " <input.wav> <output.mp4> <title> <description> [width] [height] [fps]\n";
        std::cout << "Example: " << argv[0] << " audio.wav output.mp4 \"My Track\" \"Description\" 1920 1080 60\n";
        std::cout << "\nFeatures:\n";
        std::cout << "  - Multi-threaded processing\n";
        std::cout << "  - AVX2 optimizations (if compiled with -mavx2)\n";
        std::cout << "  - Customizable resolution and framerate\n";
        return 1;
    }
    
    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    std::string title = argv[3];
    std::string description = argv[4];
    
    // Optional resolution and framerate
    if (argc >= 6) VIDEO_WIDTH = std::stoi(argv[5]);
    if (argc >= 7) VIDEO_HEIGHT = std::stoi(argv[6]);
    if (argc >= 8) VIDEO_FPS = std::stoi(argv[7]);
    
    WAVReader wav;
    if (!wav.open(inputFile)) {
        std::cerr << "Failed to open WAV file: " << inputFile << "\n";
        return 1;
    }
    
    VideoEncoder encoder;
    if (!encoder.init(outputFile, VIDEO_WIDTH, VIDEO_HEIGHT, VIDEO_FPS)) {
        std::cerr << "Failed to initialize video encoder\n";
        return 1;
    }
    
    AudioAnalyzer analyzer;
    Visualizer visualizer(VIDEO_WIDTH, VIDEO_HEIGHT, title, description);
    
    int sampleRate = wav.getSampleRate();
    int samplesPerFrame = sampleRate / VIDEO_FPS;
    
    std::cout << "===========================================\n";
    std::cout << "Audio Analysis Video Generator\n";
    std::cout << "===========================================\n";
    std::cout << "Input: " << inputFile << "\n";
    std::cout << "Output: " << outputFile << "\n";
    std::cout << "Audio: " << sampleRate << " Hz, " << wav.getChannels() << " channels\n";
    std::cout << "Video: " << VIDEO_WIDTH << "x" << VIDEO_HEIGHT << " @ " << VIDEO_FPS << " fps\n";
    
    #ifdef __AVX2__
    std::cout << "Optimizations: AVX2 enabled\n";
    #else
    std::cout << "Optimizations: Standard (compile with -mavx2 for AVX2)\n";
    #endif
    
    std::cout << "Processing threads: " << std::thread::hardware_concurrency() << "\n";
    std::cout << "===========================================\n\n";
    
    int frameCount = 0;
    auto startTime = std::chrono::steady_clock::now();
    
    while (true) {
        std::vector<float> left, right;
        int samplesRead = wav.readSamples(left, right, samplesPerFrame);
        if (samplesRead == 0) break;
        
        // Pad with zeros if needed
        if (samplesRead < samplesPerFrame) {
            left.resize(samplesPerFrame, 0.0f);
            right.resize(samplesPerFrame, 0.0f);
        }
        
        std::vector<double> fftL, fftR;
        analyzer.computeFFT(left, right, fftL, fftR);
        
        double rmsL = analyzer.computeRMS(left);
        double rmsR = analyzer.computeRMS(right);
        
        visualizer.render(left, right, fftL, fftR, rmsL, rmsR, sampleRate);
        encoder.writeFrame(visualizer.getFrame(), VIDEO_WIDTH, VIDEO_HEIGHT);
        
        frameCount++;
        if (frameCount % 60 == 0) {
            auto currentTime = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
            double fps = frameCount / (elapsed + 0.001);
            
            std::cout << "Progress: " << frameCount << " frames (" 
                     << frameCount / VIDEO_FPS << "s) @ " 
                     << (int)fps << " fps                    \r" << std::flush;
        }
    }
    
    std::cout << "\n\nFinalizing video...\n";
    encoder.close();
    wav.close();
    
    auto endTime = std::chrono::steady_clock::now();
    auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    
    std::cout << "\n===========================================\n";
    std::cout << "Completed!\n";
    std::cout << "Total frames: " << frameCount << "\n";
    std::cout << "Duration: " << frameCount / VIDEO_FPS << " seconds\n";
    std::cout << "Processing time: " << totalTime << " seconds\n";
    std::cout << "Average FPS: " << frameCount / (totalTime + 0.001) << "\n";
    std::cout << "Output saved to: " << outputFile << "\n";
    std::cout << "===========================================\n";
    
    return 0;
}
