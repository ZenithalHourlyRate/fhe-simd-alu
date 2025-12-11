#include "math/z-polynomial.h"
#include "math/dftransform-bigcomplex.h"

namespace lbcrypto {

CSlots ZPolynomial::toCSlots() const {
    return CSlots(zN, 1, ZLinearTransform::MultZU(zN, coefficients));
}

ZPolynomial CSlots::getZPolynomial(size_t slotIndex) const {
    if (slotIndex >= zSlots) {
        OPENFHE_THROW("CSlots::getZPolynomial: slotIndex out of range");
    }
    std::vector<BigComplex> cSlotsForIndex(zN / 2);
    for (size_t i = 0; i != zN / 2; ++i) {
        cSlotsForIndex[i] = slots[slotIndex * (zN / 2) + i % (zN / 2)];
    }
    return ZPolynomial(ZLinearTransform::MultZUInverse(zN, cSlotsForIndex));
}

CSlots RPolynomial::toCSlots() const {
    auto m      = coefficients.size() * 2;
    auto cSlots = coefficients.size() / 2;
    BigCVector forward(cSlots);
    for (size_t i = 0; i != cSlots; ++i) {
        forward[i] = BigComplex(coefficients[i], coefficients[i + cSlots]);
    }
    DiscreteFourierTransformBigComplex::FFTSpecial(forward, m);
    return CSlots(zN, zSlots, forward);
}

RPolynomial CSlots::toRPolynomial() const {
    auto m             = slots.size() * 4;
    BigCVector inverse = slots;
    DiscreteFourierTransformBigComplex::FFTSpecialInv(inverse, m);

    std::vector<BigFixedPoint> rValues(2 * slots.size());
    for (size_t i = 0; i != inverse.size(); ++i) {
        rValues[i]                = inverse[i].getReal();
        rValues[i + slots.size()] = inverse[i].getImag();
    }
    return RPolynomial(zN, zSlots, rValues);
}

}  // namespace lbcrypto