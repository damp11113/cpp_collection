#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <bitset>

using namespace std;

const int M = 64; // 64-QAM
const int SQRT_M = 8; // Assuming 8x8 for 64-QAM
const double PI = 3.14159265358979323846;

// Convert binary to Gray code
int binaryToGray(int num) {
    return num ^ (num >> 1);
}

// Convert Gray code to binary
int grayToBinary(int gray) {
    int num = gray;
    while (gray >>= 1) {
        num ^= gray;
    }
    return num;
}

// Modulate a sequence of bits into 64-QAM symbols
vector<complex<double>> modulate(const vector<int>& bits) {
    vector<complex<double>> symbols;
    int num_symbols = bits.size() / 6; // 6 bits per symbol

    for (int i = 0; i < num_symbols; ++i) {
        int symbol_bits[6];
        for (int j = 0; j < 6; ++j) {
            symbol_bits[j] = bits[i * 6 + j];
        }

        // Convert to Gray code
        int gray_x = binaryToGray(symbol_bits[0] * 8 + symbol_bits[1] * 4 + symbol_bits[2] * 2 + symbol_bits[3]);
        int gray_y = binaryToGray(symbol_bits[4] * 8 + symbol_bits[5] * 4);

        // Map Gray code to 64-QAM constellation
        int x = grayToBinary(gray_x);
        int y = grayToBinary(gray_y);

        // Normalize to constellation range
        double real = (2 * x - SQRT_M + 1) * 0.5;
        double imag = (2 * y - SQRT_M + 1) * 0.5;

        symbols.push_back(complex<double>(real, imag));
    }

    return symbols;
}

// Demodulate 64-QAM symbols to a sequence of bits
vector<int> demodulate(const vector<complex<double>>& symbols) {
    vector<int> bits;

    for (const auto& symbol : symbols) {
        double real = symbol.real();
        double imag = symbol.imag();

        // Quantize the real and imaginary parts
        int x = round((real + SQRT_M - 1) * 2); // Adjusted to match normalization
        int y = round((imag + SQRT_M - 1) * 2);

        // Ensure x and y are within valid range
        x = max(0, min(x, static_cast<int>(SQRT_M - 1)));
        y = max(0, min(y, static_cast<int>(SQRT_M - 1)));

        // Convert to Gray code
        int gray_x = binaryToGray(x);
        int gray_y = binaryToGray(y);

        int symbol_bits[6];
        symbol_bits[0] = (gray_x & 8) >> 3;
        symbol_bits[1] = (gray_x & 4) >> 2;
        symbol_bits[2] = (gray_x & 2) >> 1;
        symbol_bits[3] = (gray_x & 1);
        symbol_bits[4] = (gray_y & 8) >> 3;
        symbol_bits[5] = (gray_y & 4) >> 2;

        for (int i = 0; i < 6; ++i) {
            bits.push_back(symbol_bits[i]);
        }
    }

    return bits;
}

// Convert string to binary representation
vector<int> stringToBinary(const string& str) {
    vector<int> binary;

    for (char c : str) {
        bitset<8> bits(c);
        for (int i = 7; i >= 0; --i) {
            binary.push_back(bits[i]);
        }
    }

    return binary;
}

// Convert binary representation back to string
string binaryToString(const vector<int>& binary) {
    string str;
    for (size_t i = 0; i < binary.size(); i += 8) {
        bitset<8> bits;
        for (size_t j = 0; j < 8; ++j) {
            bits[j] = binary[i + j];
        }
        str += static_cast<char>(bits.to_ulong());
    }
    return str;
}

// Simple parity check for error detection
bool checkParity(const vector<int>& bits) {
    int count = 0;
    for (int bit : bits) {
        count += bit;
    }
    return count % 2 == 0;
}

int main() {
    // Original message
    string message = "Hello World";

    // Convert message to binary
    vector<int> binary_message = stringToBinary(message);

    // Pad binary message to be divisible by 6 (for 64-QAM)
    while (binary_message.size() % 6 != 0) {
        binary_message.push_back(0); // Padding with zeros
    }

    // Modulate binary message
    vector<complex<double>> symbols = modulate(binary_message);
    cout << "Modulated symbols:" << endl;
    for (const auto& symbol : symbols) {
        cout << symbol << endl;
    }

    // Demodulate symbols back to binary
    vector<int> demod_binary = demodulate(symbols);

    // Verify parity (error detection)
    if (checkParity(demod_binary)) {
        cout << "Parity check passed. No errors detected." << endl;
    } else {
        cout << "Parity check failed. Errors detected." << endl;
    }

    // Convert binary back to string
    string decoded_message = binaryToString(demod_binary);
    cout << "Decoded message: " << decoded_message << endl;

    return 0;
}
