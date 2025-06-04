#!/bin/bash

# Define directories
INCLUDE_DIR="include"
SRC_DIR="src"
CROW_INCLUDE_PATH="thirdparty/crow/include" # Adjust this if your Crow path is different

echo "Creating project directories..."
mkdir -p "$INCLUDE_DIR"
mkdir -p "$SRC_DIR"

echo "Creating header files in $INCLUDE_DIR/..."

# fixed_point.h
cat <<EOF > "$INCLUDE_DIR/fixed_point.h"
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
EOF

# complex_int.h
cat <<EOF > "$INCLUDE_DIR/complex_int.h"
#ifndef COMPLEX_INT_H
#define COMPLEX_INT_H

#include "fixed_point.h"

// Structure for a fixed-point complex number
struct ComplexInt {
    fixed_point_t real; // Real part
    fixed_point_t imag; // Imaginary part
};

// --- Complex Arithmetic Operations (Scalar Fixed-Point) ---

// Complex addition: (a.re + b.re) + i(a.im + b.im)
inline ComplexInt cint_add(ComplexInt a, ComplexInt b) {
    return {fp_add(a.real, b.real), fp_add(a.imag, b.imag)};
}

// Complex subtraction: (a.re - b.re) + i(a.im - b.im)
inline ComplexInt cint_sub(ComplexInt a, ComplexInt b) {
    return {fp_sub(a.real, b.real), fp_sub(a.imag, b.imag)};
}

// Complex multiplication: (a+bi)(c+di) = (ac-bd) + (ad+bc)i
// All intermediate products (ac, bd, ad, bc) are fixed-point multiplications.
inline ComplexInt cint_mul(ComplexInt a, ComplexInt b) {
    fixed_point_t res_real = fp_sub(fp_mul(a.real, b.real), fp_mul(a.imag, b.imag));
    fixed_point_t res_imag = fp_add(fp_mul(a.real, b.imag), fp_mul(a.imag, b.real));
    return {res_real, res_imag};
}

// Convert float real/imaginary parts to FixedComplex
inline ComplexInt float_to_cint(float real_val, float imag_val) {
    return {float_to_fp(real_val), float_to_fp(imag_val)};
}

#endif // COMPLEX_INT_H
EOF

# fft_riscv_vec.h
cat <<EOF > "$INCLUDE_DIR/fft_riscv_vec.h"
#ifndef FFT_RISCV_VEC_H
#define FFT_RISCV_VEC_H

#include "complex_int.h" // For ComplexInt struct
#include <vector>        // For std::vector (useful for internal buffers or copies)

// Function to initialize precomputed twiddle factors.
// Should be called once at application startup.
void initialize_twiddle_factors(int max_fft_size_for_twiddles);

// Performs a 1D FFT on a given data array (length N must be a power of 2).
// This function contains the conceptual integration points for RISC-V integer vector intrinsics.
// inverse: true for Inverse FFT, false for Forward FFT.
void fft_1d_riscv_vec(ComplexInt* data, int N, bool inverse);

// Performs a 2D FFT on a square image data (M rows, N columns).
// M and N must be powers of 2.
// The image_data array is assumed to be M * N elements in row-major order.
// inverse: true for Inverse FFT, false for Forward FFT.
void fft_2d_riscv_vec(ComplexInt* image_data, int M, int N, bool inverse);

// Helper for bit-reversal permutation (in-place)
void bit_reversal_permutation(ComplexInt* data, int N);

#endif // FFT_RISCV_VEC_H
EOF

# image_compressor.h
cat <<EOF > "$INCLUDE_DIR/image_compressor.h"
#ifndef IMAGE_COMPRESSOR_H
#define IMAGE_COMPRESSOR_H

#include <vector>
#include <cstdint> // For uint8_t
#include "complex_int.h" // For ComplexInt

class ImageCompressor {
public:
    // Constructor initializes with expected image dimensions.
    ImageCompressor(int width, int height);

    // Processes raw 8-bit grayscale image data:
    // 1. Converts to fixed-point complex numbers.
    // 2. Performs 2D Forward FFT.
    // 3. Applies a simple compression filter (zeroes out high frequencies).
    // 4. Performs 2D Inverse FFT.
    // 5. Converts back to 8-bit grayscale pixels.
    // retention_ratio: A float (0.0 to 1.0) indicating how much of the frequency spectrum to retain.
    // A higher ratio means more detail, less compression.
    std::vector<uint8_t> compressImage(const std::vector<uint8_t>& raw_image_data, float retention_ratio);

private:
    int width_;
    int height_;

    // Helper to apply the compression filter in the frequency domain.
    // Zeroes out coefficients beyond a certain radius.
    void apply_frequency_filter(ComplexInt* fft_coeffs, float retention_ratio);
};

#endif // IMAGE_COMPRESSOR_H
EOF

# vector_test.h
cat <<EOF > "$INCLUDE_DIR/vector_test.h"
#ifndef VECTOR_TEST_H
#define VECTOR_TEST_H

#include <stdio.h> // For printf
#include <riscv_vector.h> // Header for RISC-V Vector Intrinsics

// Function to run the vector extension test
void run_vector_test();

#endif // VECTOR_TEST_H
EOF

echo "Creating source files in $SRC_DIR/..."

# fixed_point.cpp (empty as it's inline)
cat <<EOF > "$SRC_DIR/fixed_point.cpp"
// fixed_point.h is primarily inline functions.
// No implementation needed here.
EOF

# complex_int.cpp (empty as it's inline)
cat <<EOF > "$SRC_DIR/complex_int.cpp"
// complex_int.h is primarily inline functions.
// No implementation needed here.
EOF


# fft_riscv_vec.cpp
cat <<EOF > "$SRC_DIR/fft_riscv_vec.cpp"
#include "fft_riscv_vec.h"
#include <vector>
#include <cmath>     // For std::sin, std::cos during twiddle factor generation
#include <algorithm> // For std::swap, std::max
#include <stdexcept> // For std::runtime_error
#include <iostream>  // For debugging output

// --- RISC-V Vector Intrinsics Inclusion and Conceptual Definitions ---
#ifdef __riscv_vector
#include <riscv_vector.h>
#else
// Fallback/conceptual definitions for compilation on non-RISC-V systems
// These won't actually perform vector operations but allow code to compile.

// Basic vector types for 32-bit elements, m1 (1 vector register group)
typedef void* vint32m1_t;
typedef void* vint64m2_t; // For widening multiplies (int32*int32 -> int64)

// Conceptual intrinsic functions
inline size_t __riscv_vsetvl_e32m1(size_t avl) { return avl; } // Set vector length
inline vint32m1_t __riscv_vle32_v_i32m1(const int32_t* base, size_t vl) { return nullptr; } // Load 32-bit vector
inline void __riscv_vse32_v_i32m1(int32_t* base, vint32m1_t value, size_t vl) {} // Store 32-bit vector
inline vint32m1_t __riscv_vadd_vv_i32m1(vint32m1_t vs1, vint32m1_t vs2, size_t vl) { return nullptr; } // Vector-vector add
inline vint32m1_t __riscv_vsub_vv_i32m1(vint32m1_t vs1, vint32m1_t vs2, size_t vl) { return nullptr; } // Vector-vector sub
inline vint32m1_t __riscv_vmv_v_x_i32m1(int32_t scalar_val, size_t vl) { return nullptr; } // Broadcast scalar to vector

// --- Conceptual Fixed-Point Vector Multiplication ---
// This is a simplified representation. Actual implementation needs careful handling
// of vector register groups (m-factors) and narrowing.
// In a real scenario, this would involve:
// 1. Widening multiply (e.g., __riscv_vwmul_vv_i64m2 for signed multiply)
//    Takes two SEW-bit vectors, produces a 2*SEW-bit result.
// 2. Arithmetic Right Shift (e.g., __riscv_vsra_vi_i64m2)
// 3. Narrowing (e.g., __riscv_vnsra_wi_i32m1)
inline vint32m1_t fp_mul_vv_i32m1(vint32m1_t va, vint32m1_t vb, size_t vl) {
    // This function needs to be replaced with actual RVV intrinsics for fixed-point multiplication.
    // For now, it's a placeholder.
    return __riscv_vmv_v_x_i32m1(0, vl); // Placeholder for actual RVV fixed-point multiply
}

#endif // __riscv_vector


// Global storage for precomputed twiddle factors
// Indexed by (MAX_FFT_SIZE / len) where 'len' is current FFT stage length.
// Stores W_len^0, W_len^1, ... W_len^(len/2 - 1)
static std::vector<ComplexInt> g_twiddle_factors_for_all_lengths;
static bool g_twiddle_factors_initialized = false;
static int g_max_fft_size = 0;


void initialize_twiddle_factors(int max_fft_size_for_twiddles) {
    if (g_twiddle_factors_initialized && g_max_fft_size == max_fft_size_for_twiddles) {
        return; // Already initialized for this max size
    }

    g_max_fft_size = max_fft_size_for_twiddles;
    g_twiddle_factors_for_all_lengths.resize(max_fft_size_for_twiddles / 2);

    for (int k = 0; k < max_fft_size_for_twiddles / 2; ++k) {
        float angle = -2.0f * M_PI * k / static_cast<float>(max_fft_size_for_twiddles);
        // Using std::cos and std::sin from <cmath>
        g_twiddle_factors_for_all_lengths[k] = float_to_cint(std::cos(angle), std::sin(angle));
    }
    g_twiddle_factors_initialized = true;
    std::cout << "Twiddle factors initialized for max_fft_size: " << max_fft_size_for_twiddles << std::endl;
}


// --- Bit-Reversal Permutation ---
void bit_reversal_permutation(ComplexInt* data, int N) {
    int j = 0;
    for (int i = 0; i < N; ++i) {
        if (j > i) {
            std::swap(data[i], data[j]);
        }
        int m = N >> 1;
        while (m >= 1 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
}

// --- 1D FFT Implementation with Conceptual RVV Integration ---
void fft_1d_riscv_vec(ComplexInt* data, int N, bool inverse) {
    if ((N & (N - 1)) != 0 || N == 0) {
        throw std::runtime_error("FFT size N must be a power of 2 and greater than 0.");
    }
    if (!g_twiddle_factors_initialized || g_max_fft_size < N) {
        initialize_twiddle_factors(std::max(N, g_max_fft_size > 0 ? g_max_fft_size : 256));
    }

    bit_reversal_permutation(data, N);

    for (int len = 2; len <= N; len <<= 1) { // Current stage's sub-FFT length
        int half_len = len >> 1;
        
        ComplexInt W_base_fp = g_twiddle_factors_for_all_lengths[g_max_fft_size / len];
        if (inverse) {
            W_base_fp.imag = fp_sub(float_to_fp(0.0f), W_base_fp.imag); // Conjugate for IFFT
        }

        for (int i = 0; i < N; i += len) {
            ComplexInt W_current_fp = float_to_cint(1.0f, 0.0f); // W_len^0 = 1 + 0i

            for (int j = 0; j < half_len; ) {
                size_t vl = __riscv_vsetvl_e32m1(half_len - j);

                int32_t* u_real_ptr = &data[i + j].real;
                int32_t* u_imag_ptr = &data[i + j].imag;
                int32_t* v_real_ptr = &data[i + j + half_len].real;
                int32_t* v_imag_ptr = &data[i + j + half_len].imag;

                vint32m1_t u_real_vec = __riscv_vle32_v_i32m1(u_real_ptr, vl);
                vint32m1_t u_imag_vec = __riscv_vle32_v_i32m1(u_imag_ptr, vl);

                vint32m1_t v_real_vec = __riscv_vle32_v_i32m1(v_real_ptr, vl);
                vint32m1_t v_imag_vec = __riscv_vle32_v_i32m1(v_imag_ptr, vl);

                vint32m1_t w_real_vec = __riscv_vmv_v_x_i32m1(W_current_fp.real, vl);
                vint32m1_t w_imag_vec = __riscv_vmv_v_x_i32m1(W_current_fp.imag, vl);

                // Complex Multiplication: (V_real*W_real - V_imag*W_imag) + i(V_real*W_imag + V_imag*W_real)
                vint32m1_t vw_real_prod1 = fp_mul_vv_i32m1(v_real_vec, w_real_vec, vl);
                vint32m1_t vw_real_prod2 = fp_mul_vv_i32m1(v_imag_vec, w_imag_vec, vl);
                vint32m1_t vw_real_vec = __riscv_vsub_vv_i32m1(vw_real_prod1, vw_real_prod2, vl);

                vint32m1_t vw_imag_prod1 = fp_mul_vv_i32m1(v_real_vec, w_imag_vec, vl);
                vint32m1_t vw_imag_prod2 = fp_mul_vv_i32m1(v_imag_vec, w_real_vec, vl);
                vint32m1_t vw_imag_vec = __riscv_vadd_vv_i32m1(vw_imag_prod1, vw_imag_prod2, vl);

                // Butterfly Additions/Subtractions
                vint32m1_t res_u_real_vec = __riscv_vadd_vv_i32m1(u_real_vec, vw_real_vec, vl);
                vint32m1_t res_u_imag_vec = __riscv_vadd_vv_i32m1(u_imag_vec, vw_imag_vec, vl);

                vint32m1_t res_v_real_vec = __riscv_vsub_vv_i32m1(u_real_vec, vw_real_vec, vl);
                vint32m1_t res_v_imag_vec = __riscv_vsub_vv_i32m1(u_imag_vec, vw_imag_vec, vl);

                __riscv_vse32_v_i32m1(u_real_ptr, res_u_real_vec, vl);
                __riscv_vse32_v_i32m1(u_imag_ptr, res_u_imag_vec, vl);
                __riscv_vse32_v_i32m1(v_real_ptr, res_v_real_vec, vl);
                __riscv_vse32_v_i32m1(v_imag_ptr, res_v_imag_vec, vl);
                
                W_current_fp = cint_mul(W_current_fp, W_base_fp);

                j += vl;
            }
        }
    }

    // Normalization for IFFT (divide by N)
    if (inverse) {
        fixed_point_t N_inv_fp = float_to_fp(1.0f / static_cast<float>(N));

        for (int i = 0; i < N; ) {
            size_t vl = __riscv_vsetvl_e32m1(N - i);

            vint32m1_t v_real = __riscv_vle32_v_i32m1(&data[i].real, vl);
            vint32m1_t v_imag = __riscv_vle32_v_i32m1(&data[i].imag, vl);
            
            vint32m1_t N_inv_vec = __riscv_vmv_v_x_i32m1(N_inv_fp, vl);
            
            v_real = fp_mul_vv_i32m1(v_real, N_inv_vec, vl);
            v_imag = fp_mul_vv_i32m1(v_imag, N_inv_vec, vl);

            __riscv_vse32_v_i32m1(&data[i].real, v_real, vl);
            __riscv_vse32_v_i32m1(&data[i].imag, v_imag, vl);
            i += vl;
        }
    }
}


// --- 2D FFT Implementation ---
void fft_2d_riscv_vec(ComplexInt* image_data, int M, int N, bool inverse) {
    if ((M & (M - 1)) != 0 || (N & (N - 1)) != 0 || M == 0 || N == 0) {
        throw std::runtime_error("Image dimensions (M, N) must be powers of 2 and greater than 0.");
    }
    
    // 1. Perform 1D FFT on each row
    std::cout << "  2D FFT: Processing rows (" << M << " rows of " << N << " elements)..." << std::endl;
    for (int r = 0; r < M; ++r) {
        fft_1d_riscv_vec(&image_data[r * N], N, inverse);
    }

    // 2. Transpose the matrix
    std::vector<ComplexInt> temp_transposed_data(M * N);
    std::cout << "  2D FFT: Transposing matrix (1st transpose)..." << std::endl;
    for (int r = 0; r < M; ++r) {
        for (int c = 0; c < N; ++c) {
            temp_transposed_data[c * M + r] = image_data[r * N + c];
        }
    }
    std::copy(temp_transposed_data.begin(), temp_transposed_data.end(), image_data);


    // 3. Perform 1D FFT on each "column" (which are now rows after transpose)
    std::cout << "  2D FFT: Processing columns (now rows) (" << N << " rows of " << M << " elements)..." << std::endl;
    for (int c_original = 0; c_original < N; ++c_original) {
        fft_1d_riscv_vec(&image_data[c_original * M], M, inverse);
    }

    // 4. Transpose back to original orientation
    std::cout << "  2D FFT: Transposing matrix (2nd transpose)..." << std::endl;
    for (int r_transposed = 0; r_transposed < N; ++r_transposed) {
        for (int c_transposed = 0; c_transposed < M; ++c_transposed) {
             temp_transposed_data[c_transposed * N + r_transposed] = image_data[r_transposed * M + c_transposed];
        }
    }
    std::copy(temp_transposed_data.begin(), temp_transposed_data.end(), image_data);
    std::cout << "  2D FFT: Transpose completed." << std::endl;
}
EOF

# image_compressor.cpp
cat <<EOF > "$SRC_DIR/image_compressor.cpp"
#include "image_compressor.h"
#include "fft_riscv_vec.h" // For 2D FFT functions
#include <cmath>           // For std::sqrt, std::min, std::max
#include <stdexcept>       // For std::runtime_error
#include <iostream>        // For debugging

ImageCompressor::ImageCompressor(int width, int height)
    : width_(width), height_(height)
{
    if ((width & (width - 1)) != 0 || (height & (height - 1)) != 0 || width == 0 || height == 0) {
        throw std::runtime_error("Image dimensions must be powers of 2 and greater than 0 for FFT.");
    }
    initialize_twiddle_factors(std::max(width, height));
}

std::vector<uint8_t> ImageCompressor::compressImage(const std::vector<uint8_t>& raw_image_data, float retention_ratio)
{
    if (raw_image_data.size() != static_cast<size_t>(width_ * height_)) {
        throw std::runtime_error("Raw image data size does not match specified dimensions.");
    }

    std::vector<ComplexInt> image_complex_fp(width_ * height_);
    for (int i = 0; i < width_ * height_; ++i) {
        image_complex_fp[i] = float_to_cint(static_cast<float>(raw_image_data[i]), 0.0f);
    }

    std::cout << "Performing 2D Forward FFT..." << std::endl;
    fft_2d_riscv_vec(image_complex_fp.data(), height_, width_, false); // false for forward FFT
    std::cout << "2D Forward FFT completed." << std::endl;

    std::cout << "Applying frequency filter for compression (retention ratio: " << retention_ratio << ")..." << std::endl;
    apply_frequency_filter(image_complex_fp.data(), retention_ratio);
    std::cout << "Frequency filter applied." << std::endl;

    std::cout << "Performing 2D Inverse FFT..." << std::endl;
    fft_2d_riscv_vec(image_complex_fp.data(), height_, width_, true); // true for inverse FFT
    std::cout << "2D Inverse FFT completed." << std::endl;

    std::vector<uint8_t> compressed_image_data(width_ * height_);
    for (int i = 0; i < width_ * height_; ++i) {
        float pixel_float = fp_to_float(image_complex_fp[i].real);
        
        compressed_image_data[i] = static_cast<uint8_t>(std::min(255, std::max(0, static_cast<int>(std::round(pixel_float)))));
    }
    std::cout << "Image converted back to 8-bit pixels." << std::endl;

    return compressed_image_data;
}

void ImageCompressor::apply_frequency_filter(ComplexInt* fft_coeffs, float retention_ratio) {
    float max_freq_x = static_cast<float>(width_) / 2.0f;
    float max_freq_y = static_cast<float>(height_) / 2.0f;
    float max_total_radius_sq = max_freq_x * max_freq_x + max_freq_y * max_freq_y;
    float cutoff_radius_sq = max_total_radius_sq * (retention_ratio * retention_ratio);

    for (int r = 0; r < height_; ++r) {
        for (int c = 0; c < width_; ++c) {
            float freq_x = static_cast<float>(c);
            if (c > width_ / 2) {
                freq_x = static_cast<float>(c - width_);
            }

            float freq_y = static_cast<float>(r);
            if (r > height_ / 2) {
                freq_y = static_cast<float>(r - height_);
            }

            float current_radius_sq = freq_x * freq_x + freq_y * freq_y;

            if (current_radius_sq > cutoff_radius_sq && (r != 0 || c != 0)) {
                fft_coeffs[r * width_ + c].real = float_to_fp(0.0f);
                fft_coeffs[r * width_ + c].imag = float_to_fp(0.0f);
            }
        }
    }
}
EOF

# vector_test.cpp
cat <<EOF > "$SRC_DIR/vector_test.cpp"
#include "vector_test.h"
#include <stdio.h> // For printf
#include <riscv_vector.h> // Header for RISC-V Vector Intrinsics

// Define the array size
#define ARRAY_SIZE 32

void run_vector_test() {
    printf("--- Running RISC-V Vector Test ---\n");

    // Input arrays - using 64-bit integers
    int64_t a[ARRAY_SIZE];
    int64_t b[ARRAY_SIZE];
    // Output array
    int64_t c[ARRAY_SIZE];

    // Initialize input arrays
    for (int i = 0; i < ARRAY_SIZE; i++) {
        a[i] = (int64_t)i;
        b[i] = (int64_t)(ARRAY_SIZE - 1 - i);
    }

    printf("Input A:\n");
    for (int i = 0; i < ARRAY_SIZE; i++) {
        printf("%ld ", a[i]);
    }
    printf("\n\nInput B:\n");
    for (int i = 0; i < ARRAY_SIZE; i++) {
        printf("%ld ", b[i]);
    }
    printf("\n\nPerforming 64-bit vector addition A + B...\n\n");

    // RISC-V Vector Loop
    for (size_t i = 0; i < ARRAY_SIZE; ) {
        // Set vector length for 64-bit integer elements with m1 multiplier
        size_t vl = __riscv_vsetvl_e64m1(ARRAY_SIZE - i);

        // Load vector elements (64-bit integers)
        vint64m1_t va = __riscv_vle64_v_i64m1(&a[i], vl);
        vint64m1_t vb = __riscv_vle64_v_i64m1(&b[i], vl);

        // Vector integer addition
        vint64m1_t vc = __riscv_vadd_vv_i64m1(va, vb, vl);

        // Store vector elements back to memory
        __riscv_vse64_v_i64m1(&c[i], vc, vl);

        // Increment the loop counter by the actual vector length processed
        i += vl;
    }

    printf("Output C (A + B):\n");
    for (int i = 0; i < ARRAY_SIZE; i++) {
        printf("%ld ", c[i]);
    }
    printf("\n--- RISC-V Vector Test Completed ---\n\n");
}
EOF

# server.cpp
cat <<EOF > "server.cpp"
#include "crow.h"
#include "image_compressor.h" // Our image compression class
#include "fft_riscv_vec.h"    // To ensure twiddle factors are initialized
#include "vector_test.h"      // For the vector test

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

// Define fixed image dimensions for this simple example.
// Must be powers of 2 for FFT.
const int IMAGE_WIDTH = 64; // e.g., 64x64 or 128x128
const int IMAGE_HEIGHT = 64;
const int IMAGE_PIXEL_COUNT = IMAGE_WIDTH * IMAGE_HEIGHT;

int main()
{
    // Run the RISC-V Vector test on startup
    run_vector_test();

    // Initialize twiddle factors for FFT once at application startup.
    initialize_twiddle_factors(std::max(IMAGE_WIDTH, IMAGE_HEIGHT));

    crow::SimpleApp app;

    // Default route
    CROW_ROUTE(app, "/")([]() {
        return "Welcome to the RISC-V Vector FFT Image Compressor.\n"
               "POST raw grayscale " + std::to_string(IMAGE_WIDTH) + "x" + std::to_string(IMAGE_HEIGHT) + " image data (" + std::to_string(IMAGE_PIXEL_COUNT) + " bytes) to /compress_image";
    });

    // POST endpoint for image compression
    CROW_ROUTE(app, "/compress_image").methods(crow::HTTPMethod::POST)(
        [&](const crow::request& req) {
            std::cout << "Received POST request to /compress_image." << std::endl;

            if (req.body.length() != IMAGE_PIXEL_COUNT) {
                return crow::response(400, "Bad Request: Expected exactly " + std::to_string(IMAGE_PIXEL_COUNT) + " bytes for a " + std::to_string(IMAGE_WIDTH) + "x" + std::to_string(IMAGE_HEIGHT) + " grayscale image.");
            }

            try {
                std::vector<uint8_t> input_image(req.body.begin(), req.body.end());

                ImageCompressor compressor(IMAGE_WIDTH, IMAGE_HEIGHT);
                
                float retention_ratio = 0.2f; // Configurable ratio
                std::vector<uint8_t> compressed_image_data = compressor.compressImage(input_image, retention_ratio); 

                std::cout << "Image compression successful. Output image size: " << compressed_image_data.size() << " bytes." << std::endl;

                crow::json::wvalue response;
                response["status"] = "success";
                response["message"] = "Image processed and compressed.";
                response["original_dimensions"] = std::to_string(IMAGE_WIDTH) + "x" + std::to_string(IMAGE_HEIGHT);
                response["processed_pixel_count"] = compressed_image_data.size();
                response["retention_ratio"] = retention_ratio;
                
                return crow::response(200, response);

            } catch (const std::exception& e) {
                std::cerr << "Error processing image: " << e.what() << std::endl;
                return crow::response(500, "Internal Server Error: " + std::string(e.what()));
            }
        }
    );

    CROW_ROUTE(app, "/test")([]() { 
        return "Test's Completed"; 
    });


    app.port(18080).multithreaded().run();

    return 0;
}
EOF

echo "All source and header files created."
echo "Remember to place your Crow header-only library in 'thirdparty/crow/include' relative to your project root."
echo "Now, create a 'build.sh' script to compile your project:"
echo ""

# Create build.sh
cat <<EOF > "build.sh"
#!/bin/bash

# Define compiler and flags
CXX="riscv64-linux-gnu-g++"
CFLAGS="-O3 -march=rv64gcv -mabi=lp64d -static -pthread -lrt -lm"
INCLUDE_FLAGS="-I./include -I${CROW_INCLUDE_PATH}" # Adjust CROW_INCLUDE_PATH if necessary
OUT_DIR="bin"
EXECUTABLE_NAME="image_compressor_server"

echo "Creating output directory..."
mkdir -p "\$OUT_DIR"

echo "Compiling source files..."
# Compile each .cpp file into an object file
"\$CXX" \$CFLAGS \$INCLUDE_FLAGS -c src/fixed_point.cpp -o "\$OUT_DIR"/fixed_point.o || { echo "Fixed-point compilation failed"; exit 1; }
"\$CXX" \$CFLAGS \$INCLUDE_FLAGS -c src/complex_int.cpp -o "\$OUT_DIR"/complex_int.o || { echo "Complex-int compilation failed"; exit 1; }
"\$CXX" \$CFLAGS \$INCLUDE_FLAGS -c src/fft_riscv_vec.cpp -o "\$OUT_DIR"/fft_riscv_vec.o || { echo "FFT compilation failed"; exit 1; }
"\$CXX" \$CFLAGS \$INCLUDE_FLAGS -c src/image_compressor.cpp -o "\$OUT_DIR"/image_compressor.o || { echo "Image compressor compilation failed"; exit 1; }
"\$CXX" \$CFLAGS \$INCLUDE_FLAGS -c src/vector_test.cpp -o "\$OUT_DIR"/vector_test.o || { echo "Vector test compilation failed"; exit 1; }
"\$CXX" \$CFLAGS \$INCLUDE_FLAGS -c server.cpp -o "\$OUT_DIR"/server.o || { echo "Server compilation failed"; exit 1; }

echo "Linking executable..."
# Link all object files to create the final executable
"\$CXX" \$CFLAGS "\$OUT_DIR"/*.o -o "\$OUT_DIR"/"\$EXECUTABLE_NAME" || { echo "Linking failed"; exit 1; }

echo "Build complete! Executable is at \$OUT_DIR/\$EXECUTABLE_NAME"
echo "To run with QEMU:"
echo "qemu-riscv64 -cpu rv64,v=true,vlen=128,elen=64 \$OUT_DIR/\$EXECUTABLE_NAME"

EOF

chmod +x build.sh
echo "build.sh script created and made executable."