// debug_fft_riscv_vec.cpp - Debug version with validation
#include "fft_riscv_vec.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cstring>
#include <iomanip>

#ifdef __riscv_vector
#include <riscv_vector.h>
#else
// Fallback definitions
typedef void* vint32m1_t;
typedef void* vint64m2_t;
inline size_t __riscv_vsetvl_e32m1(size_t avl) { return std::min(avl, (size_t)8); }
inline vint32m1_t __riscv_vle32_v_i32m1(const int32_t* base, size_t vl) { return nullptr; }
inline void __riscv_vse32_v_i32m1(int32_t* base, vint32m1_t value, size_t vl) {}
inline vint32m1_t __riscv_vadd_vv_i32m1(vint32m1_t vs1, vint32m1_t vs2, size_t vl) { return nullptr; }
inline vint32m1_t __riscv_vsub_vv_i32m1(vint32m1_t vs1, vint32m1_t vs2, size_t vl) { return nullptr; }
inline vint32m1_t __riscv_vmv_v_x_i32m1(int32_t scalar_val, size_t vl) { return nullptr; }
#endif

// External assembly functions (fallback if not available)
extern "C" {
    void fft_butterfly_vec(int32_t* u_real, int32_t* u_imag,
                          int32_t* v_real, int32_t* v_imag,
                          int32_t W_real, int32_t W_imag,
                          size_t vl) __attribute__((weak));
}

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
            double angle = -2.0 * M_PI * k / static_cast<double>(N);
            ComplexInt twiddle = float_to_cint(std::cos(angle), std::sin(angle));
            g_twiddle_factors.push_back(twiddle);
        }
    }
    
    g_twiddle_factors_initialized = true;
    std::cout << "Debug: Twiddle factors initialized for max_fft_size: " << max_fft_size_for_twiddles << std::endl;
}

// Get twiddle factor for specific N and k
ComplexInt get_twiddle_factor(int N, int k) {
    double angle = -2.0 * M_PI * k / static_cast<double>(N);
    return float_to_cint(std::cos(angle), std::sin(angle));
}

// Debug function to validate FFT
void validate_fft_result(ComplexInt* original, ComplexInt* result, int N, const std::string& operation) {
    double max_error = 0.0;
    double avg_error = 0.0;
    
    for (int i = 0; i < N; ++i) {
        double real_error = std::abs(fp_to_float(original[i].real) - fp_to_float(result[i].real));
        double imag_error = std::abs(fp_to_float(original[i].imag) - fp_to_float(result[i].imag));
        double total_error = real_error + imag_error;
        
        max_error = std::max(max_error, total_error);
        avg_error += total_error;
    }
    
    avg_error /= N;
    std::cout << "Debug " << operation << " - Max error: " << max_error 
              << ", Avg error: " << avg_error << std::endl;
}

void bit_reversal_permutation(ComplexInt* data, int N) {
    std::cout << "Debug: Starting bit-reversal for N=" << N << std::endl;
    
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
    
    std::cout << "Debug: Bit-reversal completed" << std::endl;
}

void debug_butterfly(ComplexInt& u, ComplexInt& v, const ComplexInt& W) {
    // Store original values for debugging
    ComplexInt orig_u = u, orig_v = v;
    
    // Complex multiplication: v * W
    ComplexInt v_twiddle = cint_mul(v, W);
    
    // Butterfly operations
    ComplexInt u_new = cint_add(u, v_twiddle);
    ComplexInt v_new = cint_sub(u, v_twiddle);
    
    u = u_new;
    v = v_new;
}

// Debug 1D FFT with extensive validation
void fft_1d_riscv_vec(ComplexInt* data, int N, bool inverse) {
    if ((N & (N - 1)) != 0 || N == 0) {
        throw std::runtime_error("FFT size N must be a power of 2 and greater than 0.");
    }

    std::cout << "Debug: Starting " << (inverse ? "IFFT" : "FFT") << " for N=" << N << std::endl;
    
    // Save original data for validation
    std::vector<ComplexInt> original_data(data, data + N);

    bit_reversal_permutation(data, N);
    
    std::cout << "Debug: First 4 values after bit-reversal:" << std::endl;
    for (int i = 0; i < std::min(4, N); ++i) {
        std::cout << "  [" << i << "] = " << fp_to_float(data[i].real) 
                  << " + " << fp_to_float(data[i].imag) << "i" << std::endl;
    }

    // FFT butterfly stages
    for (int len = 2; len <= N; len <<= 1) {
        int half_len = len >> 1;
        std::cout << "Debug: Processing stage with len=" << len << std::endl;
        
        for (int i = 0; i < N; i += len) {
            for (int j = 0; j < half_len; ++j) {
                ComplexInt W = get_twiddle_factor(len, j);
                if (inverse) {
                    W.imag = fp_sub(float_to_fp(0.0f), W.imag); // Conjugate
                }
                
                if (fft_butterfly_vec != nullptr && j + 1 < half_len) {
                    int32_t u_real = data[i + j].real;
                    int32_t u_imag = data[i + j].imag;
                    int32_t v_real = data[i + j + half_len].real;
                    int32_t v_imag = data[i + j + half_len].imag;
                    
                    fft_butterfly_vec(&u_real, &u_imag, &v_real, &v_imag, 
                                     W.real, W.imag, 1);
                    
                    data[i + j].real = u_real;
                    data[i + j].imag = u_imag;
                    data[i + j + half_len].real = v_real;
                    data[i + j + half_len].imag = v_imag;
                } else {
                    throw std::runtime_error("fft_butterfly_vec function not implemented. KABEER!!!!");
                }
            }
        }
        
        // Print some values after each stage
        if (len <= 8) {
            std::cout << "Debug: After stage len=" << len << ", first 4 values:" << std::endl;
            for (int i = 0; i < std::min(4, N); ++i) {
                std::cout << "  [" << i << "] = " << fp_to_float(data[i].real) 
                          << " + " << fp_to_float(data[i].imag) << "i" << std::endl;
            }
        }
    }

    // Normalization for IFFT
    if (inverse) {
        std::cout << "Debug: Applying IFFT normalization" << std::endl;
        fixed_point_t N_inv_fp = float_to_fp(1.0f / static_cast<float>(N));
        for (int i = 0; i < N; ++i) {
            data[i].real = fp_mul(data[i].real, N_inv_fp);
            data[i].imag = fp_mul(data[i].imag, N_inv_fp);
        }
    }
    
    std::cout << "Debug: FFT completed. Final first 4 values:" << std::endl;
    for (int i = 0; i < std::min(4, N); ++i) {
        std::cout << "  [" << i << "] = " << fp_to_float(data[i].real) 
                  << " + " << fp_to_float(data[i].imag) << "i" << std::endl;
    }
}

// Debug 2D FFT with validation
void fft_2d_riscv_vec(ComplexInt* image_data, int M, int N, bool inverse) {
    if ((M & (M - 1)) != 0 || (N & (N - 1)) != 0 || M == 0 || N == 0) {
        throw std::runtime_error("Image dimensions (M, N) must be powers of 2 and greater than 0.");
    }
    
    std::cout << "Debug: Starting 2D " << (inverse ? "IFFT" : "FFT") 
              << " for " << M << "x" << N << " image" << std::endl;
    
    // Print some input values
    std::cout << "Debug: Input corner values:" << std::endl;
    std::cout << "  [0,0] = " << fp_to_float(image_data[0].real) 
              << " + " << fp_to_float(image_data[0].imag) << "i" << std::endl;
    std::cout << "  [0,1] = " << fp_to_float(image_data[1].real) 
              << " + " << fp_to_float(image_data[1].imag) << "i" << std::endl;
    
    // Process rows
    std::cout << "Debug: Processing rows..." << std::endl;
    for (int r = 0; r < M; ++r) {
        if (r < 2) std::cout << "Debug: Processing row " << r << std::endl;
        fft_1d_riscv_vec(&image_data[r * N], N, inverse);
    }

    // Print values after row processing
    std::cout << "Debug: After row processing, corner values:" << std::endl;
    std::cout << "  [0,0] = " << fp_to_float(image_data[0].real) 
              << " + " << fp_to_float(image_data[0].imag) << "i" << std::endl;

    // Transpose
    std::cout << "Debug: Transposing matrix..." << std::endl;
    std::vector<ComplexInt> temp(M * N);
    for (int r = 0; r < M; ++r) {
        for (int c = 0; c < N; ++c) {
            temp[c * M + r] = image_data[r * N + c];
        }
    }
    std::memcpy(image_data, temp.data(), M * N * sizeof(ComplexInt));

    // Process columns (now rows after transpose)
    std::cout << "Debug: Processing columns..." << std::endl;
    for (int c = 0; c < N; ++c) {
        if (c < 2) std::cout << "Debug: Processing column " << c << std::endl;
        fft_1d_riscv_vec(&image_data[c * M], M, inverse);
    }

    // Transpose back
    std::cout << "Debug: Final transpose..." << std::endl;
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < M; ++c) {
            temp[c * N + r] = image_data[r * M + c];
        }
    }
    std::memcpy(image_data, temp.data(), M * N * sizeof(ComplexInt));
    
    // Print final values
    std::cout << "Debug: Final corner values:" << std::endl;
    std::cout << "  [0,0] = " << fp_to_float(image_data[0].real) 
              << " + " << fp_to_float(image_data[0].imag) << "i" << std::endl;
    
    std::cout << "Debug: 2D FFT completed." << std::endl;
}