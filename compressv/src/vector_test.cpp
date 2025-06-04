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
