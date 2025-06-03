#include <riscv_vector.h>
#include <math.h>
#include <stdint.h>

#define N 64  // example size, must be power of 2

typedef struct {
    float real;
    float imag;
} complex_t;
