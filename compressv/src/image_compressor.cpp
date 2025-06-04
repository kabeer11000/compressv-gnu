// Fixed image_compressor.cpp - Corrected frequency domain processing
#include "image_compressor.h"
#include "fft_riscv_vec.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cstring>

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
    
    // Convert to complex with proper normalization and centering
    for (int i = 0; i < width_ * height_; ++i) {
        // Convert uint8 to normalized float, then center around zero
        float pixel_value = (static_cast<float>(raw_image_data[i]) - 127.5f) / 127.5f; // Range [-1, 1]
        image_complex_fp[i] = float_to_cint(pixel_value, 0.0f);
    }

    std::cout << "Performing 2D Forward FFT..." << std::endl;
    fft_2d_riscv_vec(image_complex_fp.data(), height_, width_, false);
    std::cout << "2D Forward FFT completed." << std::endl;

    // Debug: Print some FFT coefficients
    std::cout << "DC component: real=" << fp_to_float(image_complex_fp[0].real) 
              << ", imag=" << fp_to_float(image_complex_fp[0].imag) << std::endl;

    std::cout << "Applying frequency filter (retention ratio: " << retention_ratio << ")..." << std::endl;
    apply_frequency_filter(image_complex_fp.data(), retention_ratio);
    std::cout << "Frequency filter applied." << std::endl;

    std::cout << "Performing 2D Inverse FFT..." << std::endl;
    fft_2d_riscv_vec(image_complex_fp.data(), height_, width_, true);
    std::cout << "2D Inverse FFT completed." << std::endl;

    std::cout << "Converting back to 8-bit pixels..." << std::endl;
    std::vector<uint8_t> compressed_image_data(width_ * height_);
    
    // Find min/max for proper scaling
    float min_val = fp_to_float(image_complex_fp[0].real);
    float max_val = fp_to_float(image_complex_fp[0].real);
    
    for (int i = 0; i < width_ * height_; ++i) {
        float pixel_float = fp_to_float(image_complex_fp[i].real);
        min_val = std::min(min_val, pixel_float);
        max_val = std::max(max_val, pixel_float);
    }
    
    std::cout << "Pixel range after IFFT: [" << min_val << ", " << max_val << "]" << std::endl;
    
    // Convert back with proper scaling
    for (int i = 0; i < width_ * height_; ++i) {
        float pixel_float = fp_to_float(image_complex_fp[i].real);
        
        // Rescale from [-1,1] back to [0,255] with proper clamping
        pixel_float = (pixel_float + 1.0f) * 127.5f; // Convert from [-1,1] to [0,255]
        
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
    
    // Use a more conservative frequency filtering approach
    // Calculate magnitude of all coefficients first
    std::vector<float> magnitudes(width_ * height_);
    std::vector<std::pair<float, int>> mag_index_pairs;
    
    for (int i = 0; i < width_ * height_; ++i) {
        float real_val = fp_to_float(fft_coeffs[i].real);
        float imag_val = fp_to_float(fft_coeffs[i].imag);
        float magnitude = std::sqrt(real_val * real_val + imag_val * imag_val);
        magnitudes[i] = magnitude;
        
        // Skip DC component (index 0) - always keep it
        if (i != 0) {
            mag_index_pairs.push_back({magnitude, i});
        }
    }
    
    // Sort by magnitude (descending)
    std::sort(mag_index_pairs.begin(), mag_index_pairs.end(), 
              [](const std::pair<float, int>& a, const std::pair<float, int>& b) {
                  return a.first > b.first;
              });
    
    // Keep only the top percentage of coefficients by magnitude
    int coeffs_to_keep = static_cast<int>((mag_index_pairs.size() * retention_ratio));
    int coeffs_kept = 1; // Count DC component
    int coeffs_zeroed = 0;
    
    // Zero out the smallest coefficients
    for (size_t i = coeffs_to_keep; i < mag_index_pairs.size(); ++i) {
        int idx = mag_index_pairs[i].second;
        fft_coeffs[idx].real = float_to_fp(0.0f);
        fft_coeffs[idx].imag = float_to_fp(0.0f);
        coeffs_zeroed++;
    }
    
    coeffs_kept += coeffs_to_keep;
    
    std::cout << "  Kept " << coeffs_kept << " coefficients (including DC)" << std::endl;
    std::cout << "  Zeroed " << coeffs_zeroed << " coefficients" << std::endl;
    std::cout << "  Actual retention ratio: " << (static_cast<float>(coeffs_kept) / (width_ * height_)) << std::endl;
    
    // Alternative: Use radial frequency filtering (commented out)
    /*
    float center_x = static_cast<float>(width_) / 2.0f;
    float center_y = static_cast<float>(height_) / 2.0f;
    float max_radius = std::sqrt(center_x * center_x + center_y * center_y);
    float cutoff_radius = max_radius * retention_ratio;
    
    int coeffs_kept = 0;
    int coeffs_zeroed = 0;

    for (int r = 0; r < height_; ++r) {
        for (int c = 0; c < width_; ++c) {
            // Calculate frequency coordinates (handling FFT shift)
            float freq_x = (c <= width_/2) ? static_cast<float>(c) : static_cast<float>(c - width_);
            float freq_y = (r <= height_/2) ? static_cast<float>(r) : static_cast<float>(r - height_);
            
            float current_radius = std::sqrt(freq_x * freq_x + freq_y * freq_y);

            // Always keep DC component, filter based on radius for others
            if (current_radius > cutoff_radius && !(r == 0 && c == 0)) {
                fft_coeffs[r * width_ + c].real = float_to_fp(0.0f);
                fft_coeffs[r * width_ + c].imag = float_to_fp(0.0f);
                coeffs_zeroed++;
            } else {
                coeffs_kept++;
            }
        }
    }
    */
}