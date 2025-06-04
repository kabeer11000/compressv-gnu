.section .text
.global fft_butterfly_vec
.type fft_butterfly_vec, @function

# void fft_butterfly_vec(int32_t* u_real, int32_t* u_imag,
#                        int32_t* v_real, int32_t* v_imag,
#                        int32_t W_real, int32_t W_imag,
#                        size_t vl);
fft_butterfly_vec:
    # Assumes: a0 = u_real, a1 = u_imag, a2 = v_real, a3 = v_imag
    #          a4 = W_real, a5 = W_imag, a6 = vl

    # Load vectors
    vsetvli t0, a6, e32, m1
    vle32.v v0, (a0)        # u_real
    vle32.v v1, (a1)        # u_imag
    vle32.v v2, (a2)        # v_real
    vle32.v v3, (a3)        # v_imag

    # Broadcast W_real and W_imag
    vmv.v.x v4, a4          # w_real
    vmv.v.x v5, a5          # w_imag

    # Complex Multiply v * W
    vmul.vv v6, v2, v4      # v_real * w_real
    vmul.vv v7, v3, v5      # v_imag * w_imag
    vsub.vv v8, v6, v7      # real result

    vmul.vv v6, v2, v5      # v_real * w_imag
    vmul.vv v7, v3, v4      # v_imag * w_real
    vadd.vv v9, v6, v7      # imag result

    # Butterfly
    vadd.vv v10, v0, v8     # u_real + v_twiddled_real
    vadd.vv v11, v1, v9     # u_imag + v_twiddled_imag
    vsub.vv v12, v0, v8     # u_real - v_twiddled_real
    vsub.vv v13, v1, v9     # u_imag - v_twiddled_imag

    # Store results
    vse32.v v10, (a0)
    vse32.v v11, (a1)
    vse32.v v12, (a2)
    vse32.v v13, (a3)

    ret
