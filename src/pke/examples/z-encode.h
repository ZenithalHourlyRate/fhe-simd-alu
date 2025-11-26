#ifndef SRC_PKE_EXAMPLES_Z_ENCODE_H_
#define SRC_PKE_EXAMPLES_Z_ENCODE_H_

#include <cassert>
#include "openfhe.h"
#include "high-prec-complex.h"

struct RPolynomial;
struct CSlots;

struct ZPolynomial {
public:
    ZPolynomial() : coefficients(zN) {}
    ZPolynomial(const std::vector<BigFixedPoint>& coeffs) : coefficients(coeffs) {}
    std::vector<BigFixedPoint> getCoefficients() const {
        return coefficients;
    }

    double getLog2Norm() {
        double maxNorm = -1e30;
        for (const auto& coeff : coefficients) {
            double coeffNorm = coeff.log2Norm();
            if (coeffNorm > maxNorm) {
                maxNorm = coeffNorm;
            }
        }
        return maxNorm;
    }

    BigFixedPoint& operator[](size_t index) {
        return coefficients[index];
    }

    const BigFixedPoint& operator[](size_t index) const {
        return coefficients[index];
    }

    static ZPolynomial getT() {
        std::vector<BigFixedPoint> t(zN, 0);
        auto two = BigFixedPoint(2, 0, false).scaleTo(z_upper_roots_scale);
        auto one = BigFixedPoint(1, 0, false).scaleTo(z_upper_roots_scale);
        t[0]     = -two;
        t[1]     = one;
        return ZPolynomial(t);
    }

    static ZPolynomial getTInv() {
        std::vector<BigFixedPoint> tInv;
        // denominator
        auto one = BigFixedPoint(1).scaleTo(z_upper_roots_scale);
        auto d   = BigFixedPoint(BigInteger(1) << zN, 0, false).scaleTo(z_upper_roots_scale);
        for (size_t i = 0; i != zN; ++i) {
            auto powerOf2 = BigFixedPoint(BigInteger(1) << (zN - 1 - i), 0, false).scaleTo(z_upper_roots_scale);
            if (i == 0) {
                tInv.push_back((one - powerOf2) / d);
            }
            else {
                tInv.push_back(-powerOf2 / d);
            }
        }
        return tInv;
    }

    static ZPolynomial addRaw(ZPolynomial a, ZPolynomial b) {
        std::vector<BigFixedPoint> result(zN, BigFixedPoint(0));
        for (size_t i = 0; i != zN; ++i) {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    static ZPolynomial multiplyRaw(ZPolynomial a, ZPolynomial b) {
        std::vector<BigFixedPoint> result(2 * zN - 1, BigFixedPoint(0));
        for (size_t i = 0; i != zN; ++i) {
            for (size_t j = 0; j != zN; ++j) {
                result[i + j] += a[i] * b[j];
            }
        }
        auto two = BigFixedPoint(2, 0, false).scaleTo(z_upper_roots_scale);
        // now euclidean reduction mod X^zN - X + 2
        for (size_t i = result.size() - 1; i >= zN; --i) {
            result[i - zN + 1] += result[i];
            result[i - zN] += -two * result[i];
            result.pop_back();
        }
        assert(result.size() == zN);
        return result;
    }

    static ZPolynomial encode(uint32_t input) {
        std::vector<BigFixedPoint> bits;
        // get bits of input in bits
        for (size_t i = 0; i != zN; ++i) {
            bits.push_back((input & (1 << i)) >> i);
        }
        return multiplyRaw(bits, getTInv());
    }

    /// This gives off-by-one representation with possible I term
    static ZPolynomial roundBigI(ZPolynomial input) {
        // Add small epsilon to ensure range in (-1, 0)
        // epsilon = 2^{-zN-1}
        auto eps = BigFixedPoint(BigInteger(1) << z_upper_roots_scale - zN - 1, z_upper_roots_scale, false);
        // do [\cdot ]_1
        std::vector<BigFixedPoint> output;
        for (auto i : input.getCoefficients()) {
            i = i - eps;
            // then do ceil
            if (i.getNeg()) {
                // make it positive
                auto ni      = -i;
                auto niFloor = (ni.getValue() >> ni.getLog2Scale());
                auto niFP    = BigFixedPoint(niFloor, 0, false);
                //std::cout << "ni: " << ni.toBinary() << " niFP: " << niFP.toBinary() << "\n";
                output.push_back(-(ni - niFP) + eps);
            }
            else {
                auto one    = BigFixedPoint(1, 0, false);
                auto iFloor = (i.getValue() >> i.getLog2Scale());
                auto iCeil  = BigFixedPoint(iFloor, 0, false) + one;
                output.push_back(i - iCeil + eps);
            }
        }
        return output;
    }

    // This gives redundant representation in [-1/2, 1/2)
    static ZPolynomial roundBigINatural(ZPolynomial input) {
        // Do []_1, which reduces to [-1/2, 1/2)
        std::vector<BigFixedPoint> output;
        for (auto i : input.getCoefficients()) {
            auto iRound = i.round();
            output.push_back(i - iRound);
        }
        return output;
    }

    static uint32_t decodeFromOffByOne(ZPolynomial input) {
        auto result = multiplyRaw(input, getT());

        // Now we have [m]_t + (X-2)I with I in {0, -1}
        // Need to remove the possible I
        // We need to add eps = 3 * 2^{-zN - 1}
        auto eps        = BigFixedPoint(BigInteger(3) << z_upper_roots_scale - zN - 1, z_upper_roots_scale, false);
        auto constCoeff = result.getCoefficients()[0] + eps;
        auto constCoeffInteger = constCoeff.getValue() >> constCoeff.getLog2Scale();
        // This is b1. OpenFHE count bits from 1...
        int32_t I = constCoeffInteger.GetBitAtIndex(2);
        auto BigI = BigFixedPoint(I, 0, false);
        result[1] += BigI;
        // Now we get pure [m]_t

        uint32_t ret = 0;
        for (size_t i = 0; i != zN; ++i) {
            // actually we should use rounding...
            // But for convenience let's just make it larger than 0
            // It becomes either 0 + small or 1 + small
            result[i] += eps;
            auto iInt = result[i].getValue() >> result[i].getLog2Scale();
            // OpenFHE count from 1???
            uint32_t iBit = iInt.GetBitAtIndex(1);
            //std::cout << "iInt " << result[i].toHexString() << " Bit " << i << "\n";
            ret += iBit << i;
        }
        return ret;
    }

    static uint32_t decode(ZPolynomial input) {
        auto poly = multiplyRaw(input, getT());

        // Do Euclidean division by X-2
        auto two                          = BigFixedPoint(2, 0, false).scaleTo(z_upper_roots_scale);
        std::vector<BigFixedPoint> result = poly.getCoefficients();
        for (size_t i = result.size() - 1; i >= 1; --i) {
            result[i - 1] += two * result[i];
            result.pop_back();
        }
        auto rounded = result[0].round();
        //std::cout << "Decoded X-2: " << result[0].toHexString() << "\n";
        return (rounded.getValue() >> rounded.getLog2Scale()).ConvertToInt();
    }

    static ZPolynomial add(ZPolynomial lhs, ZPolynomial rhs) {
        return roundBigI(addRaw(lhs, rhs));
    }

    static ZPolynomial multiply(ZPolynomial lhs, ZPolynomial rhs) {
        return roundBigI(multiplyRaw(multiplyRaw(lhs, rhs), getT()));
    }

    static ZPolynomial extractError(ZPolynomial input) {
        // Input: [m]_t / t + I + e
        // recoded: [m_t] / t
        auto decoded = decode(input);
        auto recoded = encode(decoded);
        // error: I + e
        ZPolynomial error;
        for (size_t i = 0; i != input.coefficients.size(); ++i) {
            error[i] = input[i] - recoded[i];
        }
        // This only works for off-by-one representation
        // Now there is a possible I term with I in {-1, 0}
        // We need to subtract eps = 2^{-zN - 1} to ensure < -1 or < 0
        //auto eps        = BigFixedPoint(BigInteger(1) << z_upper_roots_scale - zN - 1, z_upper_roots_scale, false);
        //auto constCoeff = error.getCoefficients()[0] - eps;
        //auto constCoeffInteger = constCoeff.getValue() >> constCoeff.getLog2Scale();
        //// This is b0. OpenFHE count bits from 1...
        //int32_t I = constCoeffInteger.GetBitAtIndex(1);
        //auto BigI = BigFixedPoint(I, 0, false);
        //error[0] += BigI;
        for (size_t i = 0; i != error.coefficients.size(); ++i) {
            auto errorRounded = error[i].round();
            error[i] -= errorRounded;
        }
        // Now we get pure e
        return error;
    }

    static ZPolynomial extractI(ZPolynomial input) {
        ZPolynomial Ipoly;
        // Input: [m]_t / t + I + e
        // recoded: [m_t] / t
        auto decoded = decode(input);
        auto recoded = encode(decoded);
        // error: I + e
        for (size_t i = 0; i != input.coefficients.size(); ++i) {
            auto diff      = input[i] - recoded[i];
            auto diffRound = diff.round();
            Ipoly[i]       = diffRound;
        }
        return Ipoly;
    }

    static ZPolynomial truncError(ZPolynomial input) {
        ZPolynomial output;
        auto error = extractError(input);
        for (size_t i = 0; i != error.coefficients.size(); ++i) {
            output[i] = input[i] - error[i];
        }
        return output;
    }

    static ZPolynomial extractStandardMessage(ZPolynomial input) {
        auto noError = truncError(input);
        auto decoded = decode(noError);
        auto recoded = encode(decoded);
        return recoded;
    }

    CSlots toCSlots() const;

    RPolynomial toRPolynomial() const;

private:
    std::vector<BigFixedPoint> coefficients;
};

struct RPolynomial {
public:
    RPolynomial() : coefficients(rN) {}
    RPolynomial(const std::vector<BigFixedPoint>& coeffs) : coefficients(coeffs) {}
    std::vector<BigFixedPoint> getCoefficients() const {
        return coefficients;
    }

    BigFixedPoint& operator[](size_t index) {
        return coefficients[index];
    }
    const BigFixedPoint& operator[](size_t index) const {
        return coefficients[index];
    }

    CSlots toCSlots() const;

    ZPolynomial toZPolynomial() const;

private:
    std::vector<BigFixedPoint> coefficients;
};

struct CSlots {
public:
    CSlots() : slots(rN / 2) {}
    CSlots(const std::vector<BigComplex>& slotVec) : slots(slotVec) {}
    std::vector<BigComplex> getSlots() const {
        return slots;
    }

    BigComplex& operator[](size_t index) {
        return slots[index];
    }
    const BigComplex& operator[](size_t index) const {
        return slots[index];
    }

    ZPolynomial toZPolynomial() const;
    RPolynomial toRPolynomial() const;

private:
    std::vector<BigComplex> slots;
};

//===
// Conversion between them
//===

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

#endif  // SRC_PKE_EXAMPLES_Z_ENCODE_H_