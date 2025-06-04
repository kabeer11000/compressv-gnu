#!/bin/bash

# Define compiler and flags
CXX="riscv64-linux-gnu-g++"
CFLAGS="-O3 -march=rv64gcv_zbb_zbs  -mabi=lp64d -static -pthread -lrt -lm"
INCLUDE_FLAGS="-I./include -Ithirdparty/crow/include" # Adjust CROW_INCLUDE_PATH if necessary
OUT_DIR="bin"
EXECUTABLE_NAME="image_compressor_server"

echo "Creating output directory..."
mkdir -p "$OUT_DIR"

echo "Compiling source files..."
# Compile each .cpp file into an object file
"$CXX" $CFLAGS $INCLUDE_FLAGS -c src/fft_butterfly_vec.s -o "$OUT_DIR"/fft_butterfly_vec.o || { echo "fft_butterfly_vec compilation failed"; exit 1; }
"$CXX" $CFLAGS $INCLUDE_FLAGS -c src/fixed_point.cpp -o "$OUT_DIR"/fixed_point.o || { echo "Fixed-point compilation failed"; exit 1; }
"$CXX" $CFLAGS $INCLUDE_FLAGS -c src/complex_int.cpp -o "$OUT_DIR"/complex_int.o || { echo "Complex-int compilation failed"; exit 1; }
"$CXX" $CFLAGS $INCLUDE_FLAGS -c src/fft_riscv_vec.cpp -o "$OUT_DIR"/fft_riscv_vec.o || { echo "FFT compilation failed"; exit 1; }
"$CXX" $CFLAGS $INCLUDE_FLAGS -c src/image_compressor.cpp -o "$OUT_DIR"/image_compressor.o || { echo "Image compressor compilation failed"; exit 1; }
"$CXX" $CFLAGS $INCLUDE_FLAGS -c src/vector_test.cpp -o "$OUT_DIR"/vector_test.o || { echo "Vector test compilation failed"; exit 1; }
"$CXX" $CFLAGS $INCLUDE_FLAGS -c server.cpp -o "$OUT_DIR"/server.o || { echo "Server compilation failed"; exit 1; }

echo "Linking executable..."
# Link all object files to create the final executable
"$CXX" $CFLAGS "$OUT_DIR"/*.o -o "$OUT_DIR"/"$EXECUTABLE_NAME" || { echo "Linking failed"; exit 1; }

echo "Build complete! Executable is at $OUT_DIR/$EXECUTABLE_NAME"
echo "To run with QEMU:"
echo "qemu-riscv64 -cpu rv64,v=true,vlen=128,elen=64 $OUT_DIR/$EXECUTABLE_NAME"

