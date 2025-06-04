#include "crow.h"
#include "image_compressor.h" // Our image compression class
#include "fft_riscv_vec.h"    // To ensure twiddle factors are initialized
#include "vector_test.h"      // For the vector test

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

#include <sstream>
#include <iomanip>
// Define fixed image dimensions for this simple example.
// Must be powers of 2 for FFT.
const int IMAGE_WIDTH = 64; // e.g., 64x64 or 128x128
const int IMAGE_HEIGHT = 64;
const int IMAGE_PIXEL_COUNT = IMAGE_WIDTH * IMAGE_HEIGHT;
std::string toHexString(const std::vector<uint8_t>& data) {
    std::ostringstream oss;
    for (uint8_t byte : data) {
        oss << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(byte);
    }
    return oss.str();
}
int main()
{
    // Run the RISC-V Vector test on startup
    run_vector_test();

    // Initialize twiddle factors for FFT once at application startup.
    initialize_twiddle_factors(std::max(IMAGE_WIDTH, IMAGE_HEIGHT));

    crow::SimpleApp app;

    // Default route
    CROW_ROUTE(app, "/")([]() {
        return "Welcome to the RISC-V Vector FFT Image Compressor.\n"
               "POST raw grayscale " + std::to_string(IMAGE_WIDTH) + "x" + std::to_string(IMAGE_HEIGHT) + " image data (" + std::to_string(IMAGE_PIXEL_COUNT) + " bytes) to /compress_image";
    });

    // POST endpoint for image compression
    CROW_ROUTE(app, "/compress_image").methods(crow::HTTPMethod::POST)(
        [&](const crow::request& req) {
            std::cout << "Received POST request to /compress_image." << std::endl;

            if (req.body.length() != IMAGE_PIXEL_COUNT) {
                return crow::response(400, "Bad Request: Expected exactly " + std::to_string(IMAGE_PIXEL_COUNT) + " bytes for a " + std::to_string(IMAGE_WIDTH) + "x" + std::to_string(IMAGE_HEIGHT) + " grayscale image.");
            }

            try {
                std::vector<uint8_t> input_image(req.body.begin(), req.body.end());

                ImageCompressor compressor(IMAGE_WIDTH, IMAGE_HEIGHT);
                
                float retention_ratio = 0.9f; // Configurable ratio
                std::vector<uint8_t> compressed_image_data = compressor.compressImage(input_image, retention_ratio); 

                std::cout << "Image compression successful. Output image size: " << compressed_image_data.size() << " bytes." << std::endl;

                crow::json::wvalue response;
                response["status"] = "success";
                response["message"] = "Image processed and compressed.";
                response["original_dimensions"] = std::to_string(IMAGE_WIDTH) + "x" + std::to_string(IMAGE_HEIGHT);
                response["processed_pixel_count"] = compressed_image_data.size();
                response["retention_ratio"] = retention_ratio;
                // Send raw bytes as hex string (you could also return it directly as octet-stream in real API)
                response["raw_hex"] = toHexString(compressed_image_data);
                
                return crow::response(200, response);

            } catch (const std::exception& e) {
                std::cerr << "Error processing image: " << e.what() << std::endl;
                return crow::response(500, "Internal Server Error: " + std::string(e.what()));
            }
        }
    );

    CROW_ROUTE(app, "/test")([]() { 
        return "Test's Completed"; 
    });


    app.port(18080).multithreaded().run();

    return 0;
}
