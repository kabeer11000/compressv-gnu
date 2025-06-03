#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> // For memcpy, etc.

// Include the RISC-V Vector intrinsics header
// You might need to adjust the path depending on your toolchain installation
#include <riscv_vector.h>

// --- Configuration and Global Defines ---
#define IMAGE_WIDTH  256
#define IMAGE_HEIGHT 256
#define PI           3.14159265358979323846f

// --- Data Structures ---
// Structure to hold complex numbers for FFT
typedef struct {
    float real;
    float imag;
} complex_float;

// --- Function Prototypes ---
// Basic image I/O (placeholder functions)
void load_grayscale_image(const char* filename, unsigned char* image_data, int width, int height);
void save_grayscale_image(const char* filename, const unsigned char* image_data, int width, int height);

// Core FFT functions
void fft_1d(complex_float* data, int N);
void fft_1d_rvv(complex_float* data, int N); // RVV-optimized 1D FFT
void fft_2d(float* input_pixels, int width, int height, complex_float* fft_output);
void ifft_1d(complex_float* data, int N);
void ifft_1d_rvv(complex_float* data, int N); // RVV-optimized 1D IFFT
void ifft_2d(complex_float* fft_input, int width, int height, float* output_pixels);

// Image processing
void apply_low_pass_filter(complex_float* fft_data, int width, int height, float cutoff_ratio);
void normalize_pixels(float* pixels, int size);
void convert_float_to_uchar(const float* float_pixels, unsigned char* uchar_pixels, int size);

// --- Main Program ---
int main() {
    // 1. Allocate memory for image data
    unsigned char* original_uchar_pixels = (unsigned char*)malloc(IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(unsigned char));
    float* float_pixels = (float*)malloc(IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(float));
    complex_float* fft_data = (complex_float*)malloc(IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(complex_float));
    float* compressed_float_pixels = (float*)malloc(IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(float));
    unsigned char* compressed_uchar_pixels = (unsigned char*)malloc(IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(unsigned char));

    if (!original_uchar_pixels || !float_pixels || !fft_data || !compressed_float_pixels || !compressed_uchar_pixels) {
        fprintf(stderr, "Memory allocation failed!\n");
        return 1;
    }

    // 2. Load Grayscale Image (replace with your actual loading mechanism)
    // For demonstration, let's fill with a simple pattern
    for (int i = 0; i < IMAGE_WIDTH * IMAGE_HEIGHT; ++i) {
        original_uchar_pixels[i] = (unsigned char)(i % 256); // Simple gradient
    }
    printf("Image loaded (dummy data).\n");

    // Convert UCHAR to FLOAT
    for (int i = 0; i < IMAGE_WIDTH * IMAGE_HEIGHT; ++i) {
        float_pixels[i] = (float)original_uchar_pixels[i];
    }

    // 3. Perform 2D FFT
    printf("Performing 2D FFT...\n");
    fft_2d(float_pixels, IMAGE_WIDTH, IMAGE_HEIGHT, fft_data);
    printf("2D FFT complete.\n");

    // 4. Apply Compression (Low-Pass Filter)
    float cutoff_ratio = 0.2f; // Keep 20% of the lowest frequencies
    printf("Applying low-pass filter with cutoff ratio: %.2f\n", cutoff_ratio);
    apply_low_pass_filter(fft_data, IMAGE_WIDTH, IMAGE_HEIGHT, cutoff_ratio);
    printf("Compression applied.\n");

    // 5. Perform Inverse 2D FFT
    printf("Performing Inverse 2D FFT...\n");
    ifft_2d(fft_data, IMAGE_WIDTH, IMAGE_HEIGHT, compressed_float_pixels);
    printf("Inverse 2D FFT complete.\n");

    // 6. Normalize and Convert to 8-bit for saving
    printf("Normalizing and converting to 8-bit...\n");
    normalize_pixels(compressed_float_pixels, IMAGE_WIDTH * IMAGE_HEIGHT); // Ensure values are in 0-255 range
    convert_float_to_uchar(compressed_float_pixels, compressed_uchar_pixels, IMAGE_WIDTH * IMAGE_HEIGHT);
    printf("Conversion complete.\n");

    // 7. Save Compressed Image (replace with your actual saving mechanism)
    // For now, you can print a few pixel values to check
    printf("First 10 compressed pixel values: ");
    for (int i = 0; i < 10; ++i) {
        printf("%d ", compressed_uchar_pixels[i]);
    }
    printf("\n");

    // Example of saving to a simple raw file (you'd need a header for a real image viewer)
    // save_grayscale_image("compressed_image.raw", compressed_uchar_pixels, IMAGE_WIDTH, IMAGE_HEIGHT);

    // Clean up memory
    free(original_uchar_pixels);
    free(float_pixels);
    free(fft_data);
    free(compressed_float_pixels);
    free(compressed_uchar_pixels);

    printf("Image compression process finished.\n");
    return 0;
}

// --- Image I/O Placeholder Functions ---
// You would replace these with actual file reading/writing for your specific raw format.
// For example, if you're loading a raw 256x256 grayscale image, it's just 256*256 bytes.
void load_grayscale_image(const char* filename, unsigned char* image_data, int width, int height) {
    // FILE* fp = fopen(filename, "rb");
    // if (!fp) {
    //     perror("Error loading image");
    //     exit(1);
    // }
    // fread(image_data, 1, width * height, fp);
    // fclose(fp);
    printf("Placeholder: Load image from %s\n", filename);
    // Fill with dummy data for testing
    for (int i = 0; i < width * height; ++i) {
        image_data[i] = (unsigned char)(rand() % 256);
    }
}

void save_grayscale_image(const char* filename, const unsigned char* image_data, int width, int height) {
    // FILE* fp = fopen(filename, "wb");
    // if (!fp) {
    //     perror("Error saving image");
    //     exit(1);
    // }
    // fwrite(image_data, 1, width * height, fp);
    // fclose(fp);
    printf("Placeholder: Image saved to %s\n", filename);
}


// --- 2D FFT and IFFT Implementations ---

// Reordering function for FFT (Bit-Reversal Permutation)
// This is typically done before the main butterfly stages.
void bit_reverse_reorder(complex_float* data, int N) {
    int i, j, k;
    for (i = 0, j = 0; i < N; i++) {
        if (j > i) {
            complex_float temp = data[j];
            data[j] = data[i];
            data[i] = temp;
        }
        k = N / 2;
        while (k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}

// Full 1D FFT implementation (Scalar version for comparison/simplicity)
// This is a Cooley-Tukey FFT algorithm
void fft_1d(complex_float* data, int N) {
    bit_reverse_reorder(data, N);

    for (int len = 2; len <= N; len <<= 1) { // len = current sub-transform size
        float angle_step = -2.0f * PI / len;
        complex_float W_len = {cosf(angle_step), sinf(angle_step)}; // Twiddle factor for this length

        for (int i = 0; i < N; i += len) {
            complex_float W = {1.0f, 0.0f}; // Current twiddle factor for this step
            for (int j = 0; j < len / 2; j++) {
                complex_float t = {
                    data[i + j + len / 2].real * W.real - data[i + j + len / 2].imag * W.imag,
                    data[i + j + len / 2].real * W.imag + data[i + j + len / 2].imag * W.real
                };
                complex_float u = data[i + j];

                data[i + j].real = u.real + t.real;
                data[i + j].imag = u.imag + t.imag;
                data[i + j + len / 2].real = u.real - t.real;
                data[i + j + len / 2].imag = u.imag - t.imag;

                // Update W for next iteration: W = W * W_len
                complex_float W_next = {
                    W.real * W_len.real - W.imag * W_len.imag,
                    W.real * W_len.imag + W.imag * W_len.real
                };
                W = W_next;
            }
        }
    }
}

// --- RVV-Optimized 1D FFT (Conceptual) ---
// This is a highly conceptual sketch. A full, correct RVV FFT is complex.
// It demonstrates how intrinsics might be used.
void fft_1d_rvv(complex_float* data, int N) {
    bit_reverse_reorder(data, N); // Bit-reversal might also be vectorized

    for (int len = 2; len <= N; len <<= 1) {
        float angle_step = -2.0f * PI / len;
        // Pre-calculate twiddle factors for the entire stage or on the fly
        // For RVV, it's often better to pre-calculate if memory allows.
        // Or calculate for `vl` elements and reuse.

        for (int i = 0; i < N; i += len) {
            // Inner loop for butterfly operations
            size_t vl = vsetvl_e32m1(len / 2); // Set vector length for f32, m1 group

            // Load Twiddle Factors (W) - these would be pre-calculated or dynamically generated
            // For simplicity, imagine 'precomputed_twiddle_real' and 'precomputed_twiddle_imag'
            // vfloat32m1_t vw_real = vle32_v_f32m1(precomputed_twiddle_reals + j_start_idx);
            // vfloat32m1_t vw_imag = vle32_v_f32m1(precomputed_twiddle_imags + j_start_idx);

            // This loop processes `vl` elements at a time
            for (size_t j = 0; j < len / 2; j += vl) {
                // Load complex data from data[i + j] and data[i + j + len/2]
                vfloat32m1_t u_real = vle32_v_f32m1(&data[i + j].real);
                vfloat32m1_t u_imag = vle32_v_f32m1(&data[i + j].imag);
                vfloat32m1_t t_real = vle32_v_f32m1(&data[i + j + len / 2].real);
                vfloat32m1_t t_imag = vle32_v_f32m1(&data[i + j + len / 2].imag);

                // For simplified demonstration, let's assume `W` is always 1+0i for now.
                // In a real FFT, you'd multiply `t` by the twiddle factor `W`.
                // vfloat32m1_t new_t_real = vfsub_vv_f32m1(vfmul_vv_f32m1(t_real, vw_real), vfmul_vv_f32m1(t_imag, vw_imag));
                // vfloat32m1_t new_t_imag = vfadd_vv_f32m1(vfmul_vv_f32m1(t_real, vw_imag), vfmul_vv_f32m1(t_imag, vw_real));
                // t_real = new_t_real;
                // t_imag = new_t_imag;

                // Butterfly operations:
                // data[i+j] = u + t
                vfloat32m1_t out_r_u_plus_t = vfadd_vv_f32m1(u_real, t_real);
                vfloat32m1_t out_i_u_plus_t = vfadd_vv_f32m1(u_imag, t_imag);

                // data[i+j+len/2] = u - t
                vfloat32m1_t out_r_u_minus_t = vfsub_vv_f32m1(u_real, t_real);
                vfloat32m1_t out_i_u_minus_t = vfsub_vv_f32m1(u_imag, t_imag);

                // Store results
                vse32_v_f32m1(&data[i + j].real, out_r_u_plus_t);
                vse32_v_f32m1(&data[i + j].imag, out_i_u_plus_t);
                vse32_v_f32m1(&data[i + j + len / 2].real, out_r_u_minus_t);
                vse32_v_f32m1(&data[i + j + len / 2].imag, out_i_u_minus_t);
            }
        }
    }
}


// 2D FFT (Row-Column approach)
void fft_2d(float* input_pixels, int width, int height, complex_float* fft_output) {
    // Copy input pixels to complex_float array for initial real parts
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fft_output[y * width + x].real = input_pixels[y * width + x];
            fft_output[y * width + x].imag = 0.0f; // Initial imaginary parts are zero
        }
    }

    // Perform 1D FFT on each row
    for (int y = 0; y < height; ++y) {
        fft_1d_rvv(&fft_output[y * width], width); // Use RVV version
    }

    // Transpose the matrix for column FFTs
    complex_float* transposed_data = (complex_float*)malloc(width * height * sizeof(complex_float));
    if (!transposed_data) {
        fprintf(stderr, "Transpose memory allocation failed!\n");
        exit(1);
    }
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            transposed_data[x * height + y] = fft_output[y * width + x];
        }
    }
    memcpy(fft_output, transposed_data, width * height * sizeof(complex_float));
    free(transposed_data);


    // Perform 1D FFT on each column (which are now rows in the transposed data)
    for (int x = 0; x < width; ++x) { // Iterate through original columns
        fft_1d_rvv(&fft_output[x * height], height); // Use RVV version
    }

    // Transpose back (optional if you don't need the standard orientation of frequency domain)
    // For proper interpretation, you might want to transpose back.
    // However, the inverse FFT will just operate on the current `fft_output` arrangement.
}

// Full 1D Inverse FFT implementation (Scalar version)
void ifft_1d(complex_float* data, int N) {
    // Conjugate the complex numbers
    for (int i = 0; i < N; ++i) {
        data[i].imag = -data[i].imag;
    }

    // Perform forward FFT
    fft_1d(data, N);

    // Conjugate again and divide by N
    for (int i = 0; i < N; ++i) {
        data[i].real = (data[i].real / N);
        data[i].imag = -(data[i].imag / N);
    }
}

// --- RVV-Optimized 1D Inverse FFT (Conceptual) ---
void ifft_1d_rvv(complex_float* data, int N) {
    // Conjugate using RVV
    size_t vl = vsetvl_e32m1(N);
    for (size_t i = 0; i < N; i += vl) {
        vfloat32m1_t imag_vec = vle32_v_f32m1(&data[i].imag);
        imag_vec = vfneg_v_f32m1(imag_vec); // Negate imaginary parts
        vse32_v_f32m1(&data[i].imag, imag_vec);
    }

    // Perform forward FFT (using the RVV version)
    fft_1d_rvv(data, N);

    // Conjugate again and divide by N using RVV
    for (size_t i = 0; i < N; i += vl) {
        vfloat32m1_t real_vec = vle32_v_f32m1(&data[i].real);
        vfloat32m1_t imag_vec = vle32_v_f32m1(&data[i].imag);

        real_vec = vfdiv_vf_f32m1(real_vec, (float)N);
        imag_vec = vfneg_v_f32m1(imag_vec);
        imag_vec = vfdiv_vf_f32m1(imag_vec, (float)N);

        vse32_v_f32m1(&data[i].real, real_vec);
        vse32_v_f32m1(&data[i].imag, imag_vec);
    }
}

// 2D Inverse FFT
void ifft_2d(complex_float* fft_input, int width, int height, float* output_pixels) {
    // Transpose the matrix for column IFFTs (if previously transposed)
    // Assuming fft_input is already in the "transposed after 2D FFT" state
    // for simplicity, we directly apply IFFT to current rows (which are original columns)
    complex_float* current_fft_data = (complex_float*)malloc(width * height * sizeof(complex_float));
    if (!current_fft_data) {
        fprintf(stderr, "IFFT memory allocation failed!\n");
        exit(1);
    }
    memcpy(current_fft_data, fft_input, width * height * sizeof(complex_float));

    // Perform 1D IFFT on each "column" (rows in current orientation)
    for (int x = 0; x < width; ++x) { // Iterate through original columns
        ifft_1d_rvv(&current_fft_data[x * height], height); // Use RVV version
    }

    // Transpose back to original image orientation
    complex_float* temp_data = (complex_float*)malloc(width * height * sizeof(complex_float));
    if (!temp_data) {
        fprintf(stderr, "Transpose back memory allocation failed!\n");
        exit(1);
    }
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            temp_data[y * width + x] = current_fft_data[x * height + y];
        }
    }
    memcpy(current_fft_data, temp_data, width * height * sizeof(complex_float));
    free(temp_data);

    // Perform 1D IFFT on each row
    for (int y = 0; y < height; ++y) {
        ifft_1d_rvv(&current_fft_data[y * width], width); // Use RVV version
    }

    // Extract real parts to output pixels
    for (int i = 0; i < width * height; ++i) {
        output_pixels[i] = current_fft_data[i].real;
    }
    free(current_fft_data);
}


// --- Compression and Utility Functions ---

// Apply a simple low-pass filter in the frequency domain
// Zeroes out frequencies beyond a certain cutoff ratio from the center (DC component)
void apply_low_pass_filter(complex_float* fft_data, int width, int height, float cutoff_ratio) {
    float max_dist_sq = (width / 2.0f) * (width / 2.0f) + (height / 2.0f) * (height / 2.0f);
    float cutoff_dist_sq = max_dist_sq * (cutoff_ratio * cutoff_ratio);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // Calculate distance from the center (DC component)
            // Handle frequency domain shift (DC component is at (0,0) after IFFT)
            int fx = (x > width / 2) ? (x - width) : x;
            int fy = (y > height / 2) ? (y - height) : y;

            float dist_sq = (float)(fx * fx + fy * fy);

            if (dist_sq > cutoff_dist_sq) {
                fft_data[y * width + x].real = 0.0f;
                fft_data[y * width + x].imag = 0.0f;
            }
        }
    }
}

// Normalize float pixels to 0-255 range and clamp
void normalize_pixels(float* pixels, int size) {
    float min_val = pixels[0];
    float max_val = pixels[0];
    for (int i = 0; i < size; ++i) {
        if (pixels[i] < min_val) min_val = pixels[i];
        if (pixels[i] > max_val) max_val = pixels[i];
    }

    // Avoid division by zero if all pixels are the same
    if (fabs(max_val - min_val) < 1e-6) {
        for (int i = 0; i < size; ++i) {
            pixels[i] = 0.0f; // Or 128.0f for grayscale mid-value
        }
        return;
    }

    // Scale to 0-255
    for (int i = 0; i < size; ++i) {
        pixels[i] = 255.0f * (pixels[i] - min_val) / (max_val - min_val);
    }
}

// Convert float pixel values to unsigned char (0-255)
void convert_float_to_uchar(const float* float_pixels, unsigned char* uchar_pixels, int size) {
    for (int i = 0; i < size; ++i) {
        // Clamp to 0-255 and cast
        uchar_pixels[i] = (unsigned char)fmaxf(0.0f, fminf(255.0f, float_pixels[i]));
    }
}
