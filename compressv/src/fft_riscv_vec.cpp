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
inline vint32m1_t fp_mul_vv_i32m1(vint32m1_t va, vint32m1_t vb, size_t vl) {
    // This function needs to be replaced with actual RVV intrinsics for fixed-point multiplication.
    // For now, it's a placeholder.
    return __riscv_vmv_v_x_i32m1(0, vl); // Placeholder for actual RVV fixed-point multiply
}


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
// Function symbol imported from raw assembly file
extern "C" void fft_butterfly_vec(int32_t* u_real, int32_t* u_imag,
    int32_t* v_real, int32_t* v_imag,
    int32_t W_real, int32_t W_imag,
    size_t vl);

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
