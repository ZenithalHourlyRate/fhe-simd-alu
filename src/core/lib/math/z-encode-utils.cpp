#include "math/z-encode-utils.h"

namespace lbcrypto {

CSlots ZPolynomial::toCSlots() const {
    return ZLinearTransform::MultZU(zN, coefficients);
}

CSlots RPolynomial::toCSlots() const {
    return multU(getRU(), coefficients);
}

ZPolynomial CSlots::toZPolynomial() const {
    //return ZLinearTransform::MultZUInverse(zBits, slots);
    return multUInverse(getZUInverse(), slots);
}

RPolynomial CSlots::toRPolynomial() const {
    return multUInverse(getRUInverse(), slots);
}

ZPolynomial RPolynomial::toZPolynomial() const {
    return toCSlots().toZPolynomial();
}

RPolynomial ZPolynomial::toRPolynomial() const {
    return toCSlots().toRPolynomial();
}

// This is for interpreting ZPolynomial coeffs as RPolynomial directly
RPolynomial ZPolynomial::interpretAsRPolynomial() const {
    return RPolynomial(coefficients);
}

// This is for interpreting RPolynomial coeffs as ZPolynomial directly
ZPolynomial RPolynomial::interpretAsZPolynomial() const {
    return ZPolynomial(coefficients);
}

}  // namespace lbcrypto