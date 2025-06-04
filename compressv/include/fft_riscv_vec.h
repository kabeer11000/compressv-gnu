#ifndef FFT_RISCV_VEC_H
#define FFT_RISCV_VEC_H

#include "complex_int.h" // For ComplexInt struct
#include <vector>        // For std::vector (useful for internal buffers or copies)

#ifdef __cplusplus
extern "C" {
#endif

// Assembly functions for vectorized butterfly operations
// Basic butterfly operation on vector elements
void fft_butterfly_vec(int32_t* u_real, int32_t* u_imag,
                       int32_t* v_real, int32_t* v_imag,
                       int32_t W_real, int32_t W_imag,
                       size_t vl);

// Strided butterfly operation for processing multiple batches
void fft_butterfly_vec_stride(int32_t* data_real, int32_t* data_imag,
                             size_t stride, int32_t W_real, int32_t W_imag,
                             size_t count, size_t vl);

#ifdef __cplusplus
}
#endif

// Function to initialize precomputed twiddle factors.
// Should be called once at application startup.
void initialize_twiddle_factors(int max_fft_size_for_twiddles);

// Performs a 1D FFT on a given data array (length N must be a power of 2).
// This function uses RISC-V vector intrinsics and assembly optimizations.
// inverse: true for Inverse FFT, false for Forward FFT.
void fft_1d_riscv_vec(ComplexInt* data, int N, bool inverse);

// Performs a 2D FFT on a square image data (M rows, N columns).
// M and N must be powers of 2.
// The image_data array is assumed to be M * N elements in row-major order.
// inverse: true for Inverse FFT, false for Forward FFT.
void fft_2d_riscv_vec(ComplexInt* image_data, int M, int N, bool inverse);

// Helper for vectorized bit-reversal permutation (in-place)
void bit_reversal_permutation(ComplexInt* data, int N);

// Utility functions for vector operations
#ifdef __riscv_vector
#include <riscv_vector.h>

// Vector fixed-point multiplication helper
inline vint32m1_t fp_mul_vv_i32m1(vint32m1_t va, vint32m1_t vb, size_t vl);

// Vector twiddle factor computation
void compute_twiddle_vectors(int32_t* twiddle_real, int32_t* twiddle_imag, 
                           int len, int count);

// Vectorized memory operations for complex data
void vector_complex_load(ComplexInt* src, int32_t* real_dest, int32_t* imag_dest, 
                        size_t count);
void vector_complex_store(int32_t* real_src, int32_t* imag_src, ComplexInt* dest, 
                         size_t count);

// Vectorized transpose for 2D FFT
void vector_transpose_complex(ComplexInt* src, ComplexInt* dest, int M, int N);

#endif // __riscv_vector

#endif // FFT_RISCV_VEC_H