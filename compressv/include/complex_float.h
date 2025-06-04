// complex_float.h - Improved complex number operations using float
#ifndef COMPLEX_FLOAT_H
#define COMPLEX_FLOAT_H

#include <cmath>
#include <cstdint>

struct ComplexFloat {
    float real;
    float imag;
    
    ComplexFloat() : real(0.0f), imag(0.0f) {}
    ComplexFloat(float r, float i) : real(r), imag(i) {}
    
    // Copy constructor and assignment
    ComplexFloat(const ComplexFloat& other) : real(other.real), imag(other.imag) {}
    ComplexFloat& operator=(const ComplexFloat& other) {
        real = other.real;
        imag = other.imag;
        return *this;
    }
};

// Convert from uint8 pixel to normalized complex
inline ComplexFloat pixel_to_complex(uint8_t pixel) {
    // Normalize to [-1, 1] range for better FFT behavior
    float normalized = (static_cast<float>(pixel) - 127.5f) / 127.5f;
    return ComplexFloat(normalized, 0.0f);
}

// Convert complex back to uint8 pixel
inline uint8_t complex_to_pixel(const ComplexFloat& c) {
    // Take only real part and denormalize
    float denormalized = (c.real * 127.5f) + 127.5f;
    
    // Clamp to valid range
    int pixel_int = static_cast<int>(std::round(denormalized));
    return static_cast<uint8_t>(std::max(0, std::min(255, pixel_int)));
}

// Complex arithmetic operations
inline ComplexFloat complex_add(const ComplexFloat& a, const ComplexFloat& b) {
    return ComplexFloat(a.real + b.real, a.imag + b.imag);
}

inline ComplexFloat complex_sub(const ComplexFloat& a, const ComplexFloat& b) {
    return ComplexFloat(a.real - b.real, a.imag - b.imag);
}

inline ComplexFloat complex_mul(const ComplexFloat& a, const ComplexFloat& b) {
    return ComplexFloat(
        a.real * b.real - a.imag * b.imag,
        a.real * b.imag + a.imag * b.real
    );
}

inline float complex_magnitude_squared(const ComplexFloat& a) {
    return a.real * a.real + a.imag * a.imag;
}

inline float complex_magnitude(const ComplexFloat& a) {
    return std::sqrt(complex_magnitude_squared(a));
}

// Twiddle factor generation
inline ComplexFloat twiddle_factor(int k, int N, bool inverse = false) {
    double angle = -2.0 * M_PI * k / static_cast<double>(N);
    if (inverse) angle = -angle;
    return ComplexFloat(std::cos(angle), std::sin(angle));
}

// Conjugate for inverse FFT
inline ComplexFloat complex_conjugate(const ComplexFloat& a) {
    return ComplexFloat(a.real, -a.imag);
}

#endif // COMPLEX_FLOAT_H