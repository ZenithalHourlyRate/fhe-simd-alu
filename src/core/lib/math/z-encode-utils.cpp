#include "math/z-encode-utils.h"

CSlots ZPolynomial::toCSlots() const {
    return multU(getZU(), coefficients);
}

CSlots RPolynomial::toCSlots() const {
    return multU(getRU(), coefficients);
}

ZPolynomial CSlots::toZPolynomial() const {
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