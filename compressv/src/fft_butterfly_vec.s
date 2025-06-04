.section .text
.global fft_butterfly_vec_float
.type fft_butterfly_vec_float, @function

# void fft_butterfly_vec_float(float* u_real, float* u_imag,
#                               float* v_real, float* v_imag,
#                               float W_real, float W_imag,
#                               size_t count)
#
# Performs vectorized FFT butterfly operation using single-precision floating point:
# temp = v * W (complex multiplication)
# u_new = u + temp
# v_new = u - temp

fft_butterfly_vec_float:
    # Arguments:
    # a0 = u_real pointer, a1 = u_imag pointer
    # a2 = v_real pointer, a3 = v_imag pointer  
    # fa0 = W_real scalar, fa1 = W_imag scalar
    # a4 = element count

    # Main processing loop
butterfly_loop:
    beqz a4, butterfly_done         # Exit if count == 0
    
    # Set vector length for 32-bit float elements
    vsetvli t0, a4, e32, m1, ta, ma

    # Load u and v vectors (single precision float)
    vle32.v v0, (a0)               # v0 = u_real vector
    vle32.v v1, (a1)               # v1 = u_imag vector  
    vle32.v v2, (a2)               # v2 = v_real vector
    vle32.v v3, (a3)               # v3 = v_imag vector

    # Broadcast twiddle factor components to vectors
    vfmv.v.f v4, fa0               # v4 = W_real (broadcast)
    vfmv.v.f v5, fa1               # v5 = W_imag (broadcast)

    # Complex multiplication: v * W = (v_real + i*v_imag) * (W_real + i*W_imag)
    # Real part: v_real*W_real - v_imag*W_imag
    # Imag part: v_real*W_imag + v_imag*W_real

    # Compute the four multiplication terms
    vfmul.vv v8, v2, v4            # v8  = v_real * W_real
    vfmul.vv v9, v3, v5            # v9  = v_imag * W_imag
    vfmul.vv v10, v2, v5           # v10 = v_real * W_imag
    vfmul.vv v11, v3, v4           # v11 = v_imag * W_real

    # Complex multiplication results
    vfsub.vv v6, v8, v9            # v6 = temp_real = v_real*W_real - v_imag*W_imag
    vfadd.vv v7, v10, v11          # v7 = temp_imag = v_real*W_imag + v_imag*W_real

    # Butterfly computation
    # u_new = u + temp
    # v_new = u - temp
    vfadd.vv v12, v0, v6           # v12 = u_real + temp_real
    vfadd.vv v13, v1, v7           # v13 = u_imag + temp_imag
    vfsub.vv v14, v0, v6           # v14 = u_real - temp_real  
    vfsub.vv v15, v1, v7           # v15 = u_imag - temp_imag

    # Store results back to memory
    vse32.v v12, (a0)              # Store new u_real
    vse32.v v13, (a1)              # Store new u_imag
    vse32.v v14, (a2)              # Store new v_real
    vse32.v v15, (a3)              # Store new v_imag

    # Update pointers and counter
    slli t1, t0, 2                 # t1 = vl * 4 (bytes per float)
    add a0, a0, t1                 # Advance u_real pointer
    add a1, a1, t1                 # Advance u_imag pointer
    add a2, a2, t1                 # Advance v_real pointer
    add a3, a3, t1                 # Advance v_imag pointer
    sub a4, a4, t0                 # Decrement count
    j butterfly_loop

butterfly_done:
    ret

.size fft_butterfly_vec_float, .-fft_butterfly_vec_float

# Optimized bit reversal permutation for float arrays
.global bit_reverse_float_vec
.type bit_reverse_float_vec, @function

# void bit_reverse_float_vec(float* real_data, float* imag_data, size_t N)
# Performs in-place bit-reversal permutation on complex float arrays

bit_reverse_float_vec:
    # a0 = real_data pointer
    # a1 = imag_data pointer  
    # a2 = N (array size, must be power of 2)
    
    li t0, 0                       # i = 0
    li t1, 0                       # j = 0
    
bit_reverse_loop:
    bge t0, a2, bit_reverse_done   # if i >= N, exit
    
    # if j > i, swap elements
    ble t1, t0, skip_swap
    
    # Load elements at positions i and j
    slli t2, t0, 2                 # t2 = i * 4 (byte offset)
    slli t3, t1, 2                 # t3 = j * 4 (byte offset)
    
    add t4, a0, t2                 # t4 = &real_data[i]
    add t5, a0, t3                 # t5 = &real_data[j]
    flw ft0, 0(t4)                 # ft0 = real_data[i]
    flw ft1, 0(t5)                 # ft1 = real_data[j]
    fsw ft1, 0(t4)                 # real_data[i] = real_data[j]
    fsw ft0, 0(t5)                 # real_data[j] = real_data[i]
    
    add t4, a1, t2                 # t4 = &imag_data[i]
    add t5, a1, t3                 # t5 = &imag_data[j]
    flw ft0, 0(t4)                 # ft0 = imag_data[i]
    flw ft1, 0(t5)                 # ft1 = imag_data[j]
    fsw ft1, 0(t4)                 # imag_data[i] = imag_data[j]
    fsw ft0, 0(t5)                 # imag_data[j] = imag_data[i]

skip_swap:
    # Update j using bit-reversal algorithm
    srli t2, a2, 1                 # t2 = N >> 1 (m = N/2)
    
bit_reverse_update_j:
    blt t2, 1, bit_reverse_add_m   # if m < 1, break
    blt t1, t2, bit_reverse_add_m  # if j < m, break
    
    sub t1, t1, t2                 # j -= m
    srli t2, t2, 1                 # m >>= 1
    j bit_reverse_update_j

bit_reverse_add_m:
    add t1, t1, t2                 # j += m
    
    addi t0, t0, 1                 # i++
    j bit_reverse_loop

bit_reverse_done:
    ret

.size bit_reverse_float_vec, .-bit_reverse_float_vec

# Vectorized memory transpose for 2D FFT
.global matrix_transpose_float_vec
.type matrix_transpose_float_vec, @function

# void matrix_transpose_float_vec(float* src_real, float* src_imag,
#                                 float* dst_real, float* dst_imag,
#                                 size_t rows, size_t cols)

matrix_transpose_float_vec:
    # a0 = src_real, a1 = src_imag
    # a2 = dst_real, a3 = dst_imag  
    # a4 = rows, a5 = cols
    
    li t0, 0                       # r = 0 (row index)
    
transpose_row_loop:
    bge t0, a4, transpose_done     # if r >= rows, exit
    
    li t1, 0                       # c = 0 (column index)
    
transpose_col_loop:
    bge t1, a5, transpose_next_row # if c >= cols, next row
    
    # Calculate source index: r * cols + c
    mul t2, t0, a5                 # t2 = r * cols
    add t2, t2, t1                 # t2 = r * cols + c
    slli t2, t2, 2                 # t2 = (r * cols + c) * 4
    
    # Calculate destination index: c * rows + r  
    mul t3, t1, a4                 # t3 = c * rows
    add t3, t3, t0                 # t3 = c * rows + r
    slli t3, t3, 2                 # t3 = (c * rows + r) * 4
    
    # Load and store real part
    add t4, a0, t2                 # t4 = &src_real[r*cols + c]
    add t5, a2, t3                 # t5 = &dst_real[c*rows + r]
    flw ft0, 0(t4)
    fsw ft0, 0(t5)
    
    # Load and store imaginary part
    add t4, a1, t2                 # t4 = &src_imag[r*cols + c]
    add t5, a3, t3                 # t5 = &dst_imag[c*rows + r]
    flw ft0, 0(t4)
    fsw ft0, 0(t5)
    
    addi t1, t1, 1                 # c++
    j transpose_col_loop

transpose_next_row:
    addi t0, t0, 1                 # r++
    j transpose_row_loop

transpose_done:
    ret

.size matrix_transpose_float_vec, .-matrix_transpose_float_vec