#include "image_compressor.h"
#include "fft_riscv_vec.h" // For 2D FFT functions
#include <cmath>           // For std::sqrt, std::min, std::max
#include <stdexcept>       // For std::runtime_error
#include <iostream>        // For debugging

ImageCompressor::ImageCompressor(int width, int height)
    : width_(width), height_(height)
{
    if ((width & (width - 1)) != 0 || (height & (height - 1)) != 0 || width == 0 || height == 0) {
        throw std::runtime_error("Image dimensions must be powers of 2 and greater than 0 for FFT.");
    }
    initialize_twiddle_factors(std::max(width, height));
}

std::vector<uint8_t> ImageCompressor::compressImage(const std::vector<uint8_t>& raw_image_data, float retention_ratio)
{
    if (raw_image_data.size() != static_cast<size_t>(width_ * height_)) {
        throw std::runtime_error("Raw image data size does not match specified dimensions.");
    }

    std::vector<ComplexInt> image_complex_fp(width_ * height_);
    for (int i = 0; i < width_ * height_; ++i) {
        image_complex_fp[i] = float_to_cint(static_cast<float>(raw_image_data[i]), 0.0f);
    }

    std::cout << "Performing 2D Forward FFT..." << std::endl;
    fft_2d_riscv_vec(image_complex_fp.data(), height_, width_, false); // false for forward FFT
    std::cout << "2D Forward FFT completed." << std::endl;

    std::cout << "Applying frequency filter for compression (retention ratio: " << retention_ratio << ")..." << std::endl;
    apply_frequency_filter(image_complex_fp.data(), retention_ratio);
    std::cout << "Frequency filter applied." << std::endl;

    std::cout << "Performing 2D Inverse FFT..." << std::endl;
    fft_2d_riscv_vec(image_complex_fp.data(), height_, width_, true); // true for inverse FFT
    std::cout << "2D Inverse FFT completed." << std::endl;

    std::vector<uint8_t> compressed_image_data(width_ * height_);
    for (int i = 0; i < width_ * height_; ++i) {
        float pixel_float = fp_to_float(image_complex_fp[i].real);
        
        compressed_image_data[i] = static_cast<uint8_t>(std::min(255, std::max(0, static_cast<int>(std::round(pixel_float)))));
    }
    std::cout << "Image converted back to 8-bit pixels." << std::endl;

    return compressed_image_data;
}

void ImageCompressor::apply_frequency_filter(ComplexInt* fft_coeffs, float retention_ratio) {
    float max_freq_x = static_cast<float>(width_) / 2.0f;
    float max_freq_y = static_cast<float>(height_) / 2.0f;
    float max_total_radius_sq = max_freq_x * max_freq_x + max_freq_y * max_freq_y;
    float cutoff_radius_sq = max_total_radius_sq * (retention_ratio * retention_ratio);

    for (int r = 0; r < height_; ++r) {
        for (int c = 0; c < width_; ++c) {
            float freq_x = static_cast<float>(c);
            if (c > width_ / 2) {
                freq_x = static_cast<float>(c - width_);
            }

            float freq_y = static_cast<float>(r);
            if (r > height_ / 2) {
                freq_y = static_cast<float>(r - height_);
            }

            float current_radius_sq = freq_x * freq_x + freq_y * freq_y;

            if (current_radius_sq > cutoff_radius_sq && (r != 0 || c != 0)) {
                fft_coeffs[r * width_ + c].real = float_to_fp(0.0f);
                fft_coeffs[r * width_ + c].imag = float_to_fp(0.0f);
            }
        }
    }
}
