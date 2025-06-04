// fft_riscv_vec_float.cpp - Enhanced FFT with floating point precision
#include "fft_riscv_vec_float.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cstring>

#ifdef __riscv_vector
#include <riscv_vector.h>
#endif

// External assembly functions
extern "C" {
    void fft_butterfly_vec_float(float* u_real, float* u_imag,
                                 float* v_real, float* v_imag,
                                 float W_real, float W_imag,
                                 size_t count) __attribute__((weak));
    
    void bit_reverse_float_vec(float* real_data, float* imag_data, 
                              size_t N) __attribute__((weak));
    
    void matrix_transpose_float_vec(float* src_real, float* src_imag,
                                   float* dst_real, float* dst_imag,
                                   size_t rows, size_t cols) __attribute__((weak));
}

// Global twiddle factor cache
static std::vector<ComplexFloat> g_twiddle_cache;
static bool g_twiddle_initialized = false;
static int g_max_cached_size = 0;

void initialize_twiddle_factors_float(int max_fft_size) {
    if (g_twiddle_initialized && g_max_cached_size >= max_fft_size) {
        return;
    }

    g_max_cached_size = max_fft_size;
    g_twiddle_cache.clear();
    g_twiddle_cache.reserve(max_fft_size);

    std::cout << "Initializing twiddle factors for max size: " << max_fft_size << std::endl;

    // Pre-compute all twiddle factors for all possible FFT sizes up to max_fft_size
    for (int N = 2; N <= max_fft_size; N <<= 1) {
        for (int k = 0; k < N/2; ++k) {
            double angle = -2.0 * M_PI * k / static_cast<double>(N);
            g_twiddle_cache.push_back(ComplexFloat(std::cos(angle), std::sin(angle)));
        }
    }
    
    g_twiddle_initialized = true;
    std::cout << "Twiddle factors initialized. Cache size: " << g_twiddle_cache.size() << std::endl;
}

// Fast twiddle factor lookup
ComplexFloat get_twiddle_factor_cached(int N, int k) {
    // For now, compute on-demand. In production, use cache lookup
    double angle = -2.0 * M_PI * k / static_cast<double>(N);
    return ComplexFloat(std::cos(angle), std::sin(angle));
}

// Software fallback for bit reversal
void bit_reversal_permutation_float(ComplexFloat* data, int N) {
    std::cout << "Performing bit-reversal for N=" << N << std::endl;
    
    if (bit_reverse_float_vec != nullptr) {
        // Use optimized assembly version
        std::vector<float> real_data(N), imag_data(N);
        
        // Extract real and imaginary parts
        for (int i = 0; i < N; ++i) {
            real_data[i] = data[i].real;
            imag_data[i] = data[i].imag;
        }
        
        bit_reverse_float_vec(real_data.data(), imag_data.data(), N);
        
        // Repack into ComplexFloat array
        for (int i = 0; i < N; ++i) {
            data[i].real = real_data[i];
            data[i].imag = imag_data[i];
        }
    } else {
        // Software fallback
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
}

// Optimized 1D FFT with better precision and error handling
void fft_1d_float(ComplexFloat* data, int N, bool inverse) {
    if ((N & (N - 1)) != 0 || N == 0) {
        throw std::runtime_error("FFT size N must be a power of 2 and greater than 0.");
    }

    std::cout << "Starting " << (inverse ? "IFFT" : "FFT") << " for N=" << N << std::endl;

    // Step 1: Bit-reversal permutation
    bit_reversal_permutation_float(data, N);
    
    // Print first few values after bit-reversal for debugging
    if (N <= 16) {
        std::cout << "After bit-reversal:" << std::endl;
        for (int i = 0; i < std::min(8, N); ++i) {
            std::cout << "  [" << i << "] = " << std::fixed << std::setprecision(4) 
                      << data[i].real << " + " << data[i].imag << "i" << std::endl;
        }
    }

    // Step 2: FFT butterfly computation
    for (int len = 2; len <= N; len <<= 1) {
        int half_len = len >> 1;
        
        for (int i = 0; i < N; i += len) {
            if (fft_butterfly_vec_float != nullptr && half_len >= 4) {
                // Use vectorized assembly for larger blocks
                std::vector<float> u_real(half_len), u_imag(half_len);
                std::vector<float> v_real(half_len), v_imag(half_len);
                std::vector<float> W_real_vec(half_len), W_imag_vec(half_len);
                
                // Extract data and compute twiddle factors
                for (int j = 0; j < half_len; ++j) {
                    u_real[j] = data[i + j].real;
                    u_imag[j] = data[i + j].imag;
                    v_real[j] = data[i + j + half_len].real;
                    v_imag[j] = data[i + j + half_len].imag;
                    
                    ComplexFloat W = get_twiddle_factor_cached(len, j);
                    if (inverse) W = complex_conjugate(W);
                    W_real_vec[j] = W.real;
                    W_imag_vec[j] = W.imag;
                }
                
                // Process in chunks that fit vector registers
                for (int j = 0; j < half_len; j += 8) { // Process 8 elements at a time
                    int chunk_size = std::min(8, half_len - j);
                    
                    fft_butterfly_vec_float(&u_real[j], &u_imag[j],
                                           &v_real[j], &v_imag[j],
                                           W_real_vec[j], W_imag_vec[j],
                                           chunk_size);
                }
                
                // Store results back
                for (int j = 0; j < half_len; ++j) {
                    data[i + j].real = u_real[j];
                    data[i + j].imag = u_imag[j];
                    data[i + j + half_len].real = v_real[j];
                    data[i + j + half_len].imag = v_imag[j];
                }
            } else {
                // Software fallback for small blocks or when assembly not available
                for (int j = 0; j < half_len; ++j) {
                    ComplexFloat W = get_twiddle_factor_cached(len, j);
                    if (inverse) W = complex_conjugate(W);
                    
                    ComplexFloat u = data[i + j];
                    ComplexFloat v = data[i + j + half_len];
                    ComplexFloat v_twiddle = complex_mul(v, W);
                    
                    data[i + j] = complex_add(u, v_twiddle);
                    data[i + j + half_len] = complex_sub(u, v_twiddle);
                }
            }
        }
        
        // Debug output for small FFTs
        if (len <= 8 && N <= 16) {
            std::cout << "After stage len=" << len << ":" << std::endl;
            for (int i = 0; i < std::min(8, N); ++i) {
                std::cout << "  [" << i << "] = " << std::fixed << std::setprecision(4)
                          << data[i].real << " + " << data[i].imag << "i" << std::endl;
            }
        }
    }

    // Step 3: Normalization for IFFT
    if (inverse) {
        float norm_factor = 1.0f / static_cast<float>(N);
        for (int i = 0; i < N; ++i) {
            data[i].real *= norm_factor;
            data[i].imag *= norm_factor;
        }
        std::cout << "Applied IFFT normalization factor: " << norm_factor << std::endl;
    }
    
    std::cout << "FFT completed successfully" << std::endl;
}

// Enhanced 2D FFT with optimized transpose
void fft_2d_float(ComplexFloat* image_data, int M, int N, bool inverse) {
    if ((M & (M - 1)) != 0 || (N & (N - 1)) != 0 || M == 0 || N == 0) {
        throw std::runtime_error("Image dimensions (M, N) must be powers of 2 and greater than 0.");
    }
    
    std::cout << "Starting 2D " << (inverse ? "IFFT" : "FFT") 
              << " for " << M << "x" << N << " image" << std::endl;

    // Step 1: Process all rows
    std::cout << "Processing rows..." << std::endl;
    for (int r = 0; r < M; ++r) {
        fft_1d_float(&image_data[r * N], N, inverse);
    }

    // Step 2: Transpose matrix for column processing
    std::cout << "Transposing matrix..." << std::endl;
    std::vector<ComplexFloat> temp(M * N);
    
    if (matrix_transpose_float_vec != nullptr) {
        // Use optimized assembly transpose
        std::vector<float> src_real(M * N), src_imag(M * N);
        std::vector<float> dst_real(M * N), dst_imag(M * N);
        
        // Extract real/imaginary parts
        for (int i = 0; i < M * N; ++i) {
            src_real[i] = image_data[i].real;
            src_imag[i] = image_data[i].imag;
        }
        
        matrix_transpose_float_vec(src_real.data(), src_imag.data(),
                                  dst_real.data(), dst_imag.data(),
                                  M, N);
        
        // Repack transposed data
        for (int i = 0; i < M * N; ++i) {
            image_data[i].real = dst_real[i];
            image_data[i].imag = dst_imag[i];
        }
    } else {
        // Software transpose fallback
        for (int r = 0; r < M; ++r) {
            for (int c = 0; c < N; ++c) {
                temp[c * M + r] = image_data[r * N + c];
            }
        }
        std::memcpy(image_data, temp.data(), M * N * sizeof(ComplexFloat));
    }

    // Step 3: Process all columns (now rows after transpose)
    std::cout << "Processing columns..." << std::endl;
    for (int c = 0; c < N; ++c) {
        fft_1d_float(&image_data[c * M], M, inverse);
    }

    // Step 4: Transpose back to original orientation
    std::cout << "Final transpose..." << std::endl;
    if (matrix_transpose_float_vec != nullptr) {
        // Use optimized assembly transpose
        std::vector<float> src_real(M * N), src_imag(M * N);
        std::vector<float> dst_real(M * N), dst_imag(M * N);
        
        for (int i = 0; i < M * N; ++i) {
            src_real[i] = image_data[i].real;
            src_imag[i] = image_data[i].imag;
        }
        
        matrix_transpose_float_vec(src_real.data(), src_imag.data(),
                                  dst_real.data(), dst_imag.data(),
                                  N, M); // Note: swapped dimensions
        
        for (int i = 0; i < M * N; ++i) {
            image_data[i].real = dst_real[i];
            image_data[i].imag = dst_imag[i];
        }
    } else {
        // Software transpose back
        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < M; ++c) {
                temp[c * N + r] = image_data[r * M + c];
            }
        }
        std::memcpy(image_data, temp.data(), M * N * sizeof(ComplexFloat));
    }
    
    std::cout << "2D FFT completed successfully" << std::endl;
}