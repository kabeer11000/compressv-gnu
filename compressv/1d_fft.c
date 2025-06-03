#include <riscv_vector.h>
#include <math.h>
#include <stdint.h>
#include "./fft.h"

#include <riscv_vector.h>
#include <math.h>

void fft_1d_vector(float* real, float* imag, int n) {
    for (int len = 2; len <= n; len *= 2) {
        int half = len / 2;
        float theta = -2.0f * M_PI / len;

        for (int i = 0; i < n; i += len) {
            for (int j = 0; j < half; j += 1) { // Replace with vectorized loop after fixing access
                float angle = theta * j;
                float cos_theta = cosi(angle);
                float sin_theta = sini(angle);

                float tre = real[i + j + half] * cos_theta - imag[i + j + half] * sin_theta;
                float tim = real[i + j + half] * sin_theta + imag[i + j + half] * cos_theta;

                float u_real = real[i + j];
                float u_imag = imag[i + j];

                real[i + j] = u_real + tre;
                imag[i + j] = u_imag + tim;
                real[i + j + half] = u_real - tre;
                imag[i + j + half] = u_imag - tim;
            }
        }
    }
}
