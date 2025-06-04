#ifndef FIXED_POINT_H
#define FIXED_POINT_H

#include <stdint.h>
#include <cmath> // For std::round

// Define the number of fractional bits for our fixed-point numbers.
// Q=15 for int32_t means: 1 sign bit, 16 integer bits, 15 fractional bits.
// This allows a range of approx +/- 32767 for the integer part and a fractional precision of 1/32768.
// For FFT, values can grow, so you might need a larger Q for intermediate results
// or larger fixed_point_t (e.g., int64_t for intermediate products).
#define FIXED_POINT_Q_BITS 15
#define FIXED_POINT_SCALE (1 << FIXED_POINT_Q_BITS) // Represents 1.0 in fixed-point

typedef int32_t fixed_point_t; // Our fixed-point number type

// --- Conversion Functions ---

// Convert float to fixed-point_t
inline fixed_point_t float_to_fp(float val) {
    return static_cast<fixed_point_t>(std::round(val * FIXED_POINT_SCALE));
}

// Convert fixed-point_t to float
inline float fp_to_float(fixed_point_t val) {
    return static_cast<float>(val) / FIXED_POINT_SCALE;
}

// --- Arithmetic Operations ---

// Fixed-point multiplication: (a * b) / 2^Q_BITS
// Uses int64_t for intermediate product to prevent overflow before scaling.
inline fixed_point_t fp_mul(fixed_point_t a, fixed_point_t b) {
    return static_cast<fixed_point_t>((static_cast<int64_t>(a) * b) >> FIXED_POINT_Q_BITS);
}

// Fixed-point addition
inline fixed_point_t fp_add(fixed_point_t a, fixed_point_t b) {
    return a + b;
}

// Fixed-point subtraction
inline fixed_point_t fp_sub(fixed_point_t a, fixed_point_t b) {
    return a - b;
}

#endif // FIXED_POINT_H
