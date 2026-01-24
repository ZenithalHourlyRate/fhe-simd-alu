#include "math/z-polynomial.h"
#include "math/dftransform-bigcomplex.h"

namespace lbcrypto {

CSlots ZPolynomial::toCSlots() const {
    return CSlots(ZEncodingParams(ZMode, zN, 1), ZLinearTransform::MultZU(zN, coefficients));
}

ZPolynomial CSlots::getZPolynomial(size_t slotIndex) const {
    if (!params.isZMode()) {
        OPENFHE_THROW("CSlots::getZPolynomial: not in ZMode");
    }
    if (slotIndex >= params.getZSlots()) {
        OPENFHE_THROW("CSlots::getZPolynomial: slotIndex out of range");
    }
    auto zN = params.getZN();
    std::vector<BigComplex> cSlotsForIndex(zN / 2, BigFixedPoint::zero());
    for (size_t i = 0; i != zN / 2; ++i) {
        cSlotsForIndex[i] = slots[(slotIndex * (zN / 2)) + i];
        for (size_t j = 0; j != params.getZDeg(); ++j) {
            // scale down it by zSlots as we did scale up during ZEncodingImpl::encodeZ
            cSlotsForIndex[i] /= BigFixedPoint::positive(params.getZSlots());
        }
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
    return CSlots(params, forward);
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
    return RPolynomial(params, rValues);
}

// TODO: deprecate into BMode...
std::pair<uint64_t, double> CSlots::getIntegerAndErrorAtBooleanMode(size_t slotIndex) const {
    if (!params.isZMode()) {
        OPENFHE_THROW("CSlots::getIntegerAtBooleanMode: not in ZMode");
    }
    if (slotIndex >= params.getZSlots()) {
        OPENFHE_THROW("CSlots::getIntegerAtBooleanMode: slotIndex out of range");
    }
    auto zN                = params.getZN();
    auto zSlots            = params.getZSlots();
    uint64_t reconstructed = 0;
    double log2MaxError    = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i != zN; ++i) {
        auto index = slotIndex * (zN / 2) + i;
        if (i >= zN / 2) {
            index += (zSlots - 1) * (zN / 2);
        }
        auto bit         = slots[index].getReal().convertToDouble();
        auto integerPart = std::round(bit);
        auto fracPart    = bit - integerPart;
        reconstructed += (static_cast<uint64_t>(integerPart) << (i));
        log2MaxError = std::max(log2MaxError, std::log2(std::abs(fracPart)));
    }
    return {reconstructed, log2MaxError};
}

uint64_t CSlots::getIntegerAtBooleanMode(size_t slotIndex) const {
    return getIntegerAndErrorAtBooleanMode(slotIndex).first;
}

double CSlots::getIntegerErrorAtBooleanMode(size_t slotIndex) const {
    return getIntegerAndErrorAtBooleanMode(slotIndex).second;
}

}  // namespace lbcrypto