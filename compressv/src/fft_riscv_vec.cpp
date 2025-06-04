// Fixed fft_riscv_vec.cpp
#include "fft_riscv_vec.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cstring>

#ifdef __riscv_vector
#include <riscv_vector.h>

// Proper fixed-point vector multiplication using RVV intrinsics
inline vint32m1_t fp_mul_vv_i32m1_fixed(vint32m1_t va, vint32m1_t vb, size_t vl) {
    // Widen multiply: int32 * int32 -> int64
    vint64m2_t va_wide = __riscv_vsext_vf2_i64m2(va, vl);
    vint64m2_t vb_wide = __riscv_vsext_vf2_i64m2(vb, vl);
    vint64m2_t result_wide = __riscv_vmul_vv_i64m2(va_wide, vb_wide, vl);
    
    // Shift right by 15 bits (Q15 format) and narrow back to int32
    vint64m2_t shifted = __riscv_vsra_vx_i64m2(result_wide, 15, vl);
    return __riscv_vnsra_wx_i32m1(shifted, 0, vl);
}

#else
// Fallback definitions for non-RISC-V compilation

typedef void* vint32m1_t;
typedef void* vint64m2_t;

inline size_t __riscv_vsetvl_e32m1(size_t avl) { return std::min(avl, (size_t)8); }
inline vint32m1_t __riscv_vle32_v_i32m1(const int32_t* base, size_t vl) { return nullptr; }
inline void __riscv_vse32_v_i32m1(int32_t* base, vint32m1_t value, size_t vl) {}
inline vint32m1_t __riscv_vadd_vv_i32m1(vint32m1_t vs1, vint32m1_t vs2, size_t vl) { return nullptr; }
inline vint32m1_t __riscv_vsub_vv_i32m1(vint32m1_t vs1, vint32m1_t vs2, size_t vl) { return nullptr; }
inline vint32m1_t __riscv_vmv_v_x_i32m1(int32_t scalar_val, size_t vl) { return nullptr; }

// Fallback scalar fixed-point multiply
inline vint32m1_t fp_mul_vv_i32m1_fixed(vint32m1_t va, vint32m1_t vb, size_t vl) {
    return __riscv_vmv_v_x_i32m1(0, vl);
}
#endif

// Global twiddle factors
static std::vector<ComplexInt> g_twiddle_factors;
static bool g_twiddle_factors_initialized = false;
static int g_max_fft_size = 0;

void initialize_twiddle_factors(int max_fft_size_for_twiddles) {
    if (g_twiddle_factors_initialized && g_max_fft_size >= max_fft_size_for_twiddles) {
        return;
    }

    g_max_fft_size = max_fft_size_for_twiddles;
    g_twiddle_factors.clear();
    g_twiddle_factors.reserve(max_fft_size_for_twiddles);

    // Generate all possible twiddle factors
    for (int N = 2; N <= max_fft_size_for_twiddles; N <<= 1) {
        for (int k = 0; k < N/2; ++k) {
            float angle = -2.0f * M_PI * k / static_cast<float>(N);
            ComplexInt twiddle = float_to_cint(std::cos(angle), std::sin(angle));
            g_twiddle_factors.push_back(twiddle);
        }
    }
    
    g_twiddle_factors_initialized = true;
    std::cout << "Twiddle factors initialized for max_fft_size: " << max_fft_size_for_twiddles << std::endl;
}

// Get twiddle factor for specific N and k
ComplexInt get_twiddle_factor(int N, int k) {
    // Simple lookup - in practice you'd want a more efficient indexing scheme
    float angle = -2.0f * M_PI * k / static_cast<float>(N);
    return float_to_cint(std::cos(angle), std::sin(angle));
}

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

// Corrected 1D FFT with proper memory handling
void fft_1d_riscv_vec(ComplexInt* data, int N, bool inverse) {
    if ((N & (N - 1)) != 0 || N == 0) {
        throw std::runtime_error("FFT size N must be a power of 2 and greater than 0.");
    }

    bit_reversal_permutation(data, N);

    // FFT butterfly stages
    for (int len = 2; len <= N; len <<= 1) {
        int half_len = len >> 1;
        
        for (int i = 0; i < N; i += len) {
            for (int j = 0; j < half_len; j += 4) { // Process 4 elements at a time (or adjust based on vector length)
                int remaining = std::min(4, half_len - j);
                
                // Get pointers to real and imaginary parts
                // Note: This assumes ComplexInt has real and imag as consecutive members
                int32_t* u_real_ptr = &data[i + j].real;
                int32_t* u_imag_ptr = &data[i + j].imag;
                int32_t* v_real_ptr = &data[i + j + half_len].real;
                int32_t* v_imag_ptr = &data[i + j + half_len].imag;

                // For now, do scalar operations due to stride issues
                for (int elem = 0; elem < remaining; ++elem) {
                    ComplexInt u = {u_real_ptr[elem * 2], u_imag_ptr[elem * 2]};
                    ComplexInt v = {v_real_ptr[elem * 2], v_imag_ptr[elem * 2]};
                    
                    ComplexInt W = get_twiddle_factor(len, j + elem);
                    if (inverse) {
                        W.imag = fp_sub(float_to_fp(0.0f), W.imag); // Conjugate
                    }
                    
                    ComplexInt v_twiddle = cint_mul(v, W);
                    
                    ComplexInt u_new = cint_add(u, v_twiddle);
                    ComplexInt v_new = cint_sub(u, v_twiddle);
                    
                    u_real_ptr[elem * 2] = u_new.real;
                    u_imag_ptr[elem * 2] = u_new.imag;
                    v_real_ptr[elem * 2] = v_new.real;
                    v_imag_ptr[elem * 2] = v_new.imag;
                }
            }
        }
    }

    // Normalization for IFFT
    if (inverse) {
        fixed_point_t N_inv_fp = float_to_fp(1.0f / static_cast<float>(N));
        for (int i = 0; i < N; ++i) {
            data[i].real = fp_mul(data[i].real, N_inv_fp);
            data[i].imag = fp_mul(data[i].imag, N_inv_fp);
        }
    }
}

// Simplified and corrected 2D FFT
void fft_2d_riscv_vec(ComplexInt* image_data, int M, int N, bool inverse) {
    if ((M & (M - 1)) != 0 || (N & (N - 1)) != 0 || M == 0 || N == 0) {
        throw std::runtime_error("Image dimensions (M, N) must be powers of 2 and greater than 0.");
    }
    
    std::cout << "  2D FFT: Processing " << M << "x" << N << " image..." << std::endl;
    
    // Process rows
    std::cout << "  Processing rows..." << std::endl;
    for (int r = 0; r < M; ++r) {
        fft_1d_riscv_vec(&image_data[r * N], N, inverse);
    }

    // Transpose
    std::cout << "  Transposing matrix..." << std::endl;
    std::vector<ComplexInt> temp(M * N);
    for (int r = 0; r < M; ++r) {
        for (int c = 0; c < N; ++c) {
            temp[c * M + r] = image_data[r * N + c];
        }
    }
    std::memcpy(image_data, temp.data(), M * N * sizeof(ComplexInt));

    // Process columns (now rows after transpose)
    std::cout << "  Processing columns..." << std::endl;
    for (int c = 0; c < N; ++c) {
        fft_1d_riscv_vec(&image_data[c * M], M, inverse);
    }

    // Transpose back
    std::cout << "  Final transpose..." << std::endl;
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < M; ++c) {
            temp[c * N + r] = image_data[r * M + c];
        }
    }
    std::memcpy(image_data, temp.data(), M * N * sizeof(ComplexInt));
    
    std::cout << "  2D FFT completed." << std::endl;
}