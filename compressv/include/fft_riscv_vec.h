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
