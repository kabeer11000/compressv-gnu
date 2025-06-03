#include <riscv_vector.h>
#include <math.h>
#include <stdint.h>
#include "./fft.h"

complex_t input[N][N];
complex_t output[N][N];

void fft_2d_vector(complex_t in[N][N], complex_t out[N][N]) {
    // Apply 1D FFT to rows
    for (int i = 0; i < N; i++) {
        fft_1d_vector(in[i], N);
    }

    // Transpose matrix
    complex_t temp[N][N];
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            temp[j][i] = in[i][j];

    // Apply 1D FFT to columns (now rows of transposed matrix)
    for (int i = 0; i < N; i++) {
        fft_1d_vector(temp[i], N);
    }

    // Transpose back to get final output
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            out[i][j] = temp[j][i];
}
