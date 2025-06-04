// fixed_point.h - Corrected fixed point arithmetic
#ifndef FIXED_POINT_H
#define FIXED_POINT_H

#include <cstdint>
#include <cmath>

// Use Q15 format: 1 sign bit + 16 integer bits + 15 fractional bits
typedef int32_t fixed_point_t;
const int FIXED_POINT_FRACTIONAL_BITS = 15;
const fixed_point_t FIXED_POINT_ONE = 1 << FIXED_POINT_FRACTIONAL_BITS;

// Convert float to fixed point
inline fixed_point_t float_to_fp(float f) {
    return static_cast<fixed_point_t>(f * FIXED_POINT_ONE);
}

// Convert fixed point to float
inline float fp_to_float(fixed_point_t fp) {
    return static_cast<float>(fp) / FIXED_POINT_ONE;
}

// Fixed point multiplication
inline fixed_point_t fp_mul(fixed_point_t a, fixed_point_t b) {
    int64_t result = static_cast<int64_t>(a) * static_cast<int64_t>(b);
    return static_cast<fixed_point_t>(result >> FIXED_POINT_FRACTIONAL_BITS);
}

// Fixed point addition
inline fixed_point_t fp_add(fixed_point_t a, fixed_point_t b) {
    return a + b;
}

// Fixed point subtraction
inline fixed_point_t fp_sub(fixed_point_t a, fixed_point_t b) {
    return a - b;
}

#endif // FIXED_POINT_H
