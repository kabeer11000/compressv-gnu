
// complex_int.h - Corrected complex number operations
#ifndef COMPLEX_INT_H
#define COMPLEX_INT_H

#include "fixed_point.h"

struct ComplexInt {
    fixed_point_t real;
    fixed_point_t imag;
    
    ComplexInt() : real(0), imag(0) {}
    ComplexInt(fixed_point_t r, fixed_point_t i) : real(r), imag(i) {}
};

// Convert float pair to complex fixed point
inline ComplexInt float_to_cint(float real_f, float imag_f) {
    return ComplexInt(float_to_fp(real_f), float_to_fp(imag_f));
}

// Complex addition
inline ComplexInt cint_add(const ComplexInt& a, const ComplexInt& b) {
    return ComplexInt(fp_add(a.real, b.real), fp_add(a.imag, b.imag));
}

// Complex subtraction
inline ComplexInt cint_sub(const ComplexInt& a, const ComplexInt& b) {
    return ComplexInt(fp_sub(a.real, b.real), fp_sub(a.imag, b.imag));
}

// Complex multiplication: (a + bi)(c + di) = (ac - bd) + (ad + bc)i
inline ComplexInt cint_mul(const ComplexInt& a, const ComplexInt& b) {
    fixed_point_t real_part = fp_sub(fp_mul(a.real, b.real), fp_mul(a.imag, b.imag));
    fixed_point_t imag_part = fp_add(fp_mul(a.real, b.imag), fp_mul(a.imag, b.real));
    return ComplexInt(real_part, imag_part);
}

// Complex magnitude squared (for filtering)
inline fixed_point_t cint_mag_squared(const ComplexInt& a) {
    return fp_add(fp_mul(a.real, a.real), fp_mul(a.imag, a.imag));
}

#endif // COMPLEX_INT_H