// YUV420P (limited range) to RGBA conversion, integer only.
#pragma once

#include <cstddef>
#include <cstdint>

namespace gui {

enum class ColorMatrix { Auto = 0, BT601 = 1, BT709 = 2 };

// Auto picks BT.709 for HD (height >= 720) and BT.601 below, like most players do
// for video without color metadata.
inline ColorMatrix resolve_matrix(ColorMatrix m, int height) {
    if (m != ColorMatrix::Auto) return m;
    return height >= 720 ? ColorMatrix::BT709 : ColorMatrix::BT601;
}

// yuv: tightly packed planes Y (w*h), U and V ((w+1)/2 * (h+1)/2). rgba: w*h*4.
inline void yuv420_to_rgba(const uint8_t* yuv, int w, int h, ColorMatrix m, uint8_t* rgba) {
    // 16.16 fixed point coefficients for limited-range input.
    const int cy = 76309;  // 255/219
    int crv, cgu, cgv, cbu;
    if (resolve_matrix(m, h) == ColorMatrix::BT709) {
        crv = 117489; cgu = 13975; cgv = 34925; cbu = 138438;
    } else {
        crv = 104597; cgu = 25675; cgv = 53279; cbu = 132201;
    }
    const int cw = (w + 1) / 2, ch = (h + 1) / 2;
    const uint8_t* Y = yuv;
    const uint8_t* U = Y + size_t(w) * size_t(h);
    const uint8_t* V = U + size_t(cw) * size_t(ch);
    auto clamp = [](int v) -> uint8_t { return uint8_t(v < 0 ? 0 : (v > 255 ? 255 : v)); };
    for (int y = 0; y < h; ++y) {
        const uint8_t* yr = Y + size_t(y) * size_t(w);
        const uint8_t* ur = U + size_t(y / 2) * size_t(cw);
        const uint8_t* vr = V + size_t(y / 2) * size_t(cw);
        uint8_t* o = rgba + size_t(y) * size_t(w) * 4;
        for (int x = 0; x < w; ++x, o += 4) {
            int u = ur[x >> 1] - 128, v = vr[x >> 1] - 128;
            int l = (yr[x] - 16) * cy + 32768;
            o[0] = clamp((l + crv * v) >> 16);
            o[1] = clamp((l - cgu * u - cgv * v) >> 16);
            o[2] = clamp((l + cbu * u) >> 16);
            o[3] = 255;
        }
    }
}

}  // namespace gui
