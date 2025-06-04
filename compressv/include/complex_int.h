#ifndef COMPLEX_INT_H
#define COMPLEX_INT_H

#include "fixed_point.h"

// Structure for a fixed-point complex number
struct ComplexInt {
    fixed_point_t real; // Real part
    fixed_point_t imag; // Imaginary part
};

// --- Complex Arithmetic Operations (Scalar Fixed-Point) ---

// Complex addition: (a.re + b.re) + i(a.im + b.im)
inline ComplexInt cint_add(ComplexInt a, ComplexInt b) {
    return {fp_add(a.real, b.real), fp_add(a.imag, b.imag)};
}

// Complex subtraction: (a.re - b.re) + i(a.im - b.im)
inline ComplexInt cint_sub(ComplexInt a, ComplexInt b) {
    return {fp_sub(a.real, b.real), fp_sub(a.imag, b.imag)};
}

// Complex multiplication: (a+bi)(c+di) = (ac-bd) + (ad+bc)i
// All intermediate products (ac, bd, ad, bc) are fixed-point multiplications.
inline ComplexInt cint_mul(ComplexInt a, ComplexInt b) {
    fixed_point_t res_real = fp_sub(fp_mul(a.real, b.real), fp_mul(a.imag, b.imag));
    fixed_point_t res_imag = fp_add(fp_mul(a.real, b.imag), fp_mul(a.imag, b.real));
    return {res_real, res_imag};
}

// Convert float real/imaginary parts to FixedComplex
inline ComplexInt float_to_cint(float real_val, float imag_val) {
    return {float_to_fp(real_val), float_to_fp(imag_val)};
}

#endif // COMPLEX_INT_H
