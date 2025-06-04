.section .text
.global fft_butterfly_vec
.type fft_butterfly_vec, @function

# void fft_butterfly_vec(int32_t* u_real, int32_t* u_imag,
#                        int32_t* v_real, int32_t* v_imag,
#                        int32_t W_real, int32_t W_imag,
#                        size_t vl)
#
# Performs vectorized FFT butterfly operation:
# temp = v * W (complex multiplication)
# u_new = u + temp
# v_new = u - temp
#
# Where W = W_real + i*W_imag is the twiddle factor
# and complex multiplication: (a+bi)*(c+di) = (ac-bd) + i(ad+bc)

fft_butterfly_vec:
    # Arguments:
    # a0 = u_real pointer, a1 = u_imag pointer
    # a2 = v_real pointer, a3 = v_imag pointer  
    # a4 = W_real scalar, a5 = W_imag scalar
    # a6 = vector length

    # Set vector length for 32-bit elements with m1 multiplier
    vsetvli t0, a6, e32, m1, ta, ma

    # Load u and v vectors
    vle32.v v0, (a0)    # v0 = u_real vector
    vle32.v v1, (a1)    # v1 = u_imag vector  
    vle32.v v2, (a2)    # v2 = v_real vector
    vle32.v v3, (a3)    # v3 = v_imag vector

    # Broadcast twiddle factor components to vectors
    vmv.v.x v4, a4      # v4 = W_real (broadcast)
    vmv.v.x v5, a5      # v5 = W_imag (broadcast)

    # Complex multiplication: v * W = (v_real + i*v_imag) * (W_real + i*W_imag)
    # Real part: v_real*W_real - v_imag*W_imag
    # Imag part: v_real*W_imag + v_imag*W_real

    # First, compute the widening multiplications for precision
    # Set vector length for 64-bit operations (m2 for double width)
    vsetvli t0, a6, e32, m1, ta, ma
    
    # Convert to 64-bit for multiplication to avoid overflow
    vsext.vf2 v8, v2    # v8  = v_real extended to 64-bit (m2)
    vsext.vf2 v10, v3   # v10 = v_imag extended to 64-bit (m2)
    vsext.vf2 v12, v4   # v12 = W_real extended to 64-bit (m2)  
    vsext.vf2 v14, v5   # v14 = W_imag extended to 64-bit (m2)

    # Set vector length for 64-bit operations
    vsetvli t0, a6, e64, m2, ta, ma

    # 64-bit multiplications
    vmul.vv v16, v8, v12   # v16 = v_real * W_real (64-bit)
    vmul.vv v18, v10, v14  # v18 = v_imag * W_imag (64-bit)
    vmul.vv v20, v8, v14   # v20 = v_real * W_imag (64-bit)
    vmul.vv v22, v10, v12  # v22 = v_imag * W_real (64-bit)

    # Compute complex multiplication results
    vsub.vv v24, v16, v18  # v24 = real part (64-bit)
    vadd.vv v26, v20, v22  # v26 = imag part (64-bit)

    # Shift right by 15 bits for Q15 fixed-point format
    vsra.vi v24, v24, 15   # Scale down real part
    vsra.vi v26, v26, 15   # Scale down imag part

    # Convert back to 32-bit
    vsetvli t0, a6, e32, m1, ta, ma
    vnsra.wi v6, v24, 0    # v6 = temp_real (32-bit)
    vnsra.wi v7, v26, 0    # v7 = temp_imag (32-bit)

    # Butterfly computation
    # u_new = u + temp
    # v_new = u - temp
    vadd.vv v28, v0, v6    # v28 = u_real + temp_real
    vadd.vv v29, v1, v7    # v29 = u_imag + temp_imag
    vsub.vv v30, v0, v6    # v30 = u_real - temp_real  
    vsub.vv v31, v1, v7    # v31 = u_imag - temp_imag

    # Store results back to memory
    vse32.v v28, (a0)      # Store new u_real
    vse32.v v29, (a1)      # Store new u_imag
    vse32.v v30, (a2)      # Store new v_real
    vse32.v v31, (a3)      # Store new v_imag

    ret

.size fft_butterfly_vec, .-fft_butterfly_vec

# Additional optimized function for multiple butterfly operations
.global fft_butterfly_vec_stride
.type fft_butterfly_vec_stride, @function

# void fft_butterfly_vec_stride(int32_t* data_real, int32_t* data_imag,
#                               size_t stride, int32_t W_real, int32_t W_imag,
#                               size_t count, size_t vl)
# 
# Processes multiple butterflies with strided access pattern
fft_butterfly_vec_stride:
    # a0 = data_real base, a1 = data_imag base
    # a2 = stride between u and v elements  
    # a3 = W_real, a4 = W_imag
    # a5 = count of butterfly pairs
    # a6 = vector length per operation

butterfly_loop:
    beqz a5, butterfly_done    # Exit if count == 0

    # Set vector length
    vsetvli t0, a6, e32, m1, ta, ma

    # Load strided data
    vlse32.v v0, (a0), a2      # u_real with stride
    vlse32.v v1, (a1), a2      # u_imag with stride
    
    # Calculate v addresses (u + stride/2)
    srai t1, a2, 1             # t1 = stride/2
    add t2, a0, t1             # t2 = v_real address
    add t3, a1, t1             # t3 = v_imag address
    
    vlse32.v v2, (t2), a2      # v_real with stride
    vlse32.v v3, (t3), a2      # v_imag with stride

    # Broadcast twiddle factors
    vmv.v.x v4, a3             # v4 = W_real
    vmv.v.x v5, a4             # v5 = W_imag

    # Complex multiplication with 64-bit precision
    vsext.vf2 v8, v2           # Extend v_real to 64-bit
    vsext.vf2 v10, v3          # Extend v_imag to 64-bit  
    vsext.vf2 v12, v4          # Extend W_real to 64-bit
    vsext.vf2 v14, v5          # Extend W_imag to 64-bit

    # Switch to 64-bit vector operations
    vsetvli t0, a6, e64, m2, ta, ma

    # 64-bit complex multiply: v * W
    vmul.vv v16, v8, v12       # v_real * W_real
    vmul.vv v18, v10, v14      # v_imag * W_imag  
    vmul.vv v20, v8, v14       # v_real * W_imag
    vmul.vv v22, v10, v12      # v_imag * W_real

    vsub.vv v24, v16, v18      # real = v_real*W_real - v_imag*W_imag
    vadd.vv v26, v20, v22      # imag = v_real*W_imag + v_imag*W_real

    # Q15 scaling
    vsra.vi v24, v24, 15
    vsra.vi v26, v26, 15

    # Convert back to 32-bit
    vsetvli t0, a6, e32, m1, ta, ma
    vnsra.wi v6, v24, 0        # temp_real
    vnsra.wi v7, v26, 0        # temp_imag

    # Butterfly operations
    vadd.vv v28, v0, v6        # u_new_real = u_real + temp_real
    vadd.vv v29, v1, v7        # u_new_imag = u_imag + temp_imag
    vsub.vv v30, v0, v6        # v_new_real = u_real - temp_real
    vsub.vv v31, v1, v7        # v_new_imag = u_imag - temp_imag

    # Store results with stride
    vsse32.v v28, (a0), a2     # Store u_new_real
    vsse32.v v29, (a1), a2     # Store u_new_imag
    vsse32.v v30, (t2), a2     # Store v_new_real  
    vsse32.v v31, (t3), a2     # Store v_new_imag

    # Update pointers and counter
    add a0, a0, a6             # Move to next batch
    add a1, a1, a6
    addi a5, a5, -1            # Decrement count
    j butterfly_loop

butterfly_done:
    ret

.size fft_butterfly_vec_stride, .-fft_butterfly_vec_stride