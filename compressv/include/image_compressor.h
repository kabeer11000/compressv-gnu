#ifndef IMAGE_COMPRESSOR_H
#define IMAGE_COMPRESSOR_H

#include <vector>
#include <cstdint> // For uint8_t
#include "complex_int.h" // For ComplexInt

class ImageCompressor {
public:
    // Constructor initializes with expected image dimensions.
    ImageCompressor(int width, int height);

    // Processes raw 8-bit grayscale image data:
    // 1. Converts to fixed-point complex numbers.
    // 2. Performs 2D Forward FFT.
    // 3. Applies a simple compression filter (zeroes out high frequencies).
    // 4. Performs 2D Inverse FFT.
    // 5. Converts back to 8-bit grayscale pixels.
    // retention_ratio: A float (0.0 to 1.0) indicating how much of the frequency spectrum to retain.
    // A higher ratio means more detail, less compression.
    std::vector<uint8_t> compressImage(const std::vector<uint8_t>& raw_image_data, float retention_ratio);

private:
    int width_;
    int height_;

    // Helper to apply the compression filter in the frequency domain.
    // Zeroes out coefficients beyond a certain radius.
    void apply_frequency_filter(ComplexInt* fft_coeffs, float retention_ratio);
};

#endif // IMAGE_COMPRESSOR_H
