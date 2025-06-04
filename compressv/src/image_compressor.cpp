// Fixed image_compressor.cpp
#include "image_compressor.h"
#include "fft_riscv_vec.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>

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

    std::cout << "Converting image to complex format..." << std::endl;
    std::vector<ComplexInt> image_complex_fp(width_ * height_);
    for (int i = 0; i < width_ * height_; ++i) {
        // Convert uint8 pixel to fixed-point complex number
        float pixel_value = static_cast<float>(raw_image_data[i]) / 255.0f; // Normalize to [0,1]
        image_complex_fp[i] = float_to_cint(pixel_value, 0.0f);
    }

    std::cout << "Performing 2D Forward FFT..." << std::endl;
    fft_2d_riscv_vec(image_complex_fp.data(), height_, width_, false);
    std::cout << "2D Forward FFT completed." << std::endl;

    std::cout << "Applying frequency filter (retention ratio: " << retention_ratio << ")..." << std::endl;
    apply_frequency_filter(image_complex_fp.data(), retention_ratio);
    std::cout << "Frequency filter applied." << std::endl;

    std::cout << "Performing 2D Inverse FFT..." << std::endl;
    fft_2d_riscv_vec(image_complex_fp.data(), height_, width_, true);
    std::cout << "2D Inverse FFT completed." << std::endl;

    std::cout << "Converting back to 8-bit pixels..." << std::endl;
    std::vector<uint8_t> compressed_image_data(width_ * height_);
    for (int i = 0; i < width_ * height_; ++i) {
        // Convert fixed-point back to float and then to uint8
        float pixel_float = fp_to_float(image_complex_fp[i].real) * 255.0f; // Denormalize
        
        // Clamp to valid pixel range
        int pixel_int = static_cast<int>(std::round(pixel_float));
        pixel_int = std::max(0, std::min(255, pixel_int));
        compressed_image_data[i] = static_cast<uint8_t>(pixel_int);
    }
    std::cout << "Image conversion completed." << std::endl;

    return compressed_image_data;
}

void ImageCompressor::apply_frequency_filter(ComplexInt* fft_coeffs, float retention_ratio) {
    std::cout << "  Calculating frequency cutoff..." << std::endl;
    
    // Calculate the maximum possible frequency distance from center
    float center_x = static_cast<float>(width_) / 2.0f;
    float center_y = static_cast<float>(height_) / 2.0f;
    float max_radius = std::sqrt(center_x * center_x + center_y * center_y);
    float cutoff_radius = max_radius * retention_ratio;
    
    int coeffs_kept = 0;
    int coeffs_zeroed = 0;

    for (int r = 0; r < height_; ++r) {
        for (int c = 0; c < width_; ++c) {
            // Calculate frequency coordinates relative to center
            float freq_x = static_cast<float>(c) - center_x;
            float freq_y = static_cast<float>(r) - center_y;
            
            // Handle wrapping for negative frequencies
            if (c > width_ / 2) {
                freq_x = static_cast<float>(c - width_);
            }
            if (r > height_ / 2) {
                freq_y = static_cast<float>(r - height_);
            }

            float current_radius = std::sqrt(freq_x * freq_x + freq_y * freq_y);

            // Keep DC component and low frequencies, zero out high frequencies
            if (current_radius > cutoff_radius && !(r == 0 && c == 0)) {
                fft_coeffs[r * width_ + c].real = float_to_fp(0.0f);
                fft_coeffs[r * width_ + c].imag = float_to_fp(0.0f);
                coeffs_zeroed++;
            } else {
                coeffs_kept++;
            }
        }
    }
    
    std::cout << "  Kept " << coeffs_kept << " coefficients, zeroed " << coeffs_zeroed << " coefficients." << std::endl;
    std::cout << "  Retention ratio achieved: " << (static_cast<float>(coeffs_kept) / (coeffs_kept + coeffs_zeroed)) << std::endl;
}