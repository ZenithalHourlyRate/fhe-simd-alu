#ifndef SRC_PKE_EXAMPLES_Z_ENCODE_UTILS_H_
#define SRC_PKE_EXAMPLES_Z_ENCODE_UTILS_H_

#include <cassert>
#include "math/z-constants.h"

namespace lbcrypto {

struct RPolynomial;
struct CSlots;

struct ZPolynomial {
public:
    ZPolynomial(uint32_t zN) : zN(zN), coefficients(zN) {}
    ZPolynomial(const std::vector<BigFixedPoint>& coeffs) : zN(coeffs.size()), coefficients(coeffs) {}
    std::vector<BigFixedPoint> getCoefficients() const {
        return coefficients;
    }

    uint32_t getZN() const {
        return zN;
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

    static ZPolynomial getT(uint32_t zN) {
        std::vector<BigFixedPoint> t(zN, BigFixedPoint::zero());
        auto two = BigFixedPoint::two();
        auto one = BigFixedPoint::one();
        t[0]     = -two;
        t[1]     = one;
        return ZPolynomial(t);
    }

    static ZPolynomial getTInv(uint32_t zN) {
        std::vector<BigFixedPoint> tInv;
        // denominator
        auto one = BigFixedPoint::one();
        auto d   = BigFixedPoint(BigInteger(1) << zN, 0, false).scaleTo(128);
        for (size_t i = 0; i != zN; ++i) {
            auto powerOf2 = BigFixedPoint(BigInteger(1) << (zN - 1 - i), 0, false).scaleTo(128);
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
        auto zN = a.getZN();
        std::vector<BigFixedPoint> result(zN, BigFixedPoint::zero());
        for (size_t i = 0; i != zN; ++i) {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    static ZPolynomial multiplyRaw(ZPolynomial a, ZPolynomial b) {
        auto zN = a.getZN();
        std::vector<BigFixedPoint> result(2 * zN - 1, BigFixedPoint::zero());
        for (size_t i = 0; i != zN; ++i) {
            for (size_t j = 0; j != zN; ++j) {
                result[i + j] += a[i] * b[j];
            }
        }
        auto two = BigFixedPoint::two();
        // now euclidean reduction mod X^zN - X + 2
        for (size_t i = result.size() - 1; i >= zN; --i) {
            result[i - zN + 1] += result[i];
            result[i - zN] += -two * result[i];
            result.pop_back();
        }
        assert(result.size() == zN);
        return result;
    }

    // Binary
    // TODO: support larger integer
    static ZPolynomial encodeBinary(uint32_t zN, uint64_t input) {
        auto one  = BigFixedPoint::one();
        auto zero = BigFixedPoint::zero();
        std::vector<BigFixedPoint> bits;
        // get bits of input in {0, 1}
        for (size_t i = 0; i != zN; ++i) {
            auto flag = (input & (1 << i)) >> i;
            if (flag) {
                bits.push_back(one);
            }
            else {
                bits.push_back(zero);
            }
        }
        return bits;
    }

    // Standard
    static ZPolynomial encode(uint32_t zN, uint64_t input) {
        return multiplyRaw(encodeBinary(zN, input), getTInv(zN));
    }

    // Round to [-1, 1)
    static ZPolynomial roundNOneToOne(ZPolynomial input) {
        std::vector<BigFixedPoint> output;
        auto two = BigFixedPoint::positive(2);
        for (auto i : input.getCoefficients()) {
            auto j      = i / two;
            auto jRound = j.round();
            auto jFrac  = j - jRound;
            auto iFrac  = jFrac * two;
            output.push_back(iFrac);
        }
        return ZPolynomial(output);
    }

    // Round to (-1, 0]
    static ZPolynomial roundNOneToZero(ZPolynomial input) {
        // Add small epsilon to ensure range in (-1, 0)
        // epsilon = 2^{-zN-1}
        auto zN  = input.getZN();
        auto eps = BigFixedPoint(BigInteger(1) << (128 - zN - 1), 128, false);
        // do [\cdot ]_1
        std::vector<BigFixedPoint> output;
        for (auto i : input.getCoefficients()) {
            i = i - eps;
            // then do ceil
            auto iCeil = i.ceil();
            output.push_back(i - iCeil + eps);
        }
        return output;
    }

    static ZPolynomial roundBigI(ZPolynomial input) {
        return roundNOneToZero(input);
    }

    static uint64_t decode(ZPolynomial input) {
        auto poly = multiplyRaw(input, getT(input.getZN()));

        // For each coefficient, do rounding
        std::vector<BigFixedPoint> result = poly.getCoefficients();
        for (size_t i = 0; i != result.size(); ++i) {
            auto rounded = result[i].round();
            result[i]    = rounded;
        }
        // Do Euclidean division by X-2
        auto two = BigFixedPoint::positive(2);
        for (size_t i = result.size() - 1; i >= 1; --i) {
            result[i - 1] += two * result[i];
            result.pop_back();
        }
        auto rounded = result[0].round();
        // Now do mod 2^32
        auto modValueBFP = BigFixedPoint(BigInteger(1) << 32, 0, false).scaleTo(128);
        rounded          = rounded - (rounded / modValueBFP).floor() * modValueBFP;
        //std::cout << "Decoded X-2: " << result[0].toHexString() << "\n";
        return (rounded.getValue() >> rounded.getLog2Scale()).ConvertToInt();
    }

    static ZPolynomial add(ZPolynomial lhs, ZPolynomial rhs) {
        return roundBigI(addRaw(lhs, rhs));
    }

    static ZPolynomial multiply(ZPolynomial lhs, ZPolynomial rhs) {
        return roundBigI(multiplyRaw(multiplyRaw(lhs, rhs), getT(lhs.getZN())));
    }

    static ZPolynomial extractError(ZPolynomial input) {
        // Input: [m]_t / t + I + e
        // recoded: [m_t] / t
        auto zN      = input.getZN();
        auto decoded = decode(input);
        auto recoded = encode(zN, decoded);
        // error: I + e
        ZPolynomial error(zN);
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
        ZPolynomial Ipoly(input.getZN());
        // Input: [m]_t / t + I + e
        // recoded: [m_t] / t
        auto decoded = decode(input);
        auto recoded = encode(input.getZN(), decoded);
        // error: I + e
        for (size_t i = 0; i != input.coefficients.size(); ++i) {
            auto diff      = input[i] - recoded[i];
            auto diffRound = diff.round();
            Ipoly[i]       = diffRound;
        }
        return Ipoly;
    }

    static ZPolynomial truncError(ZPolynomial input) {
        ZPolynomial output(input.getZN());
        auto error = extractError(input);
        for (size_t i = 0; i != error.coefficients.size(); ++i) {
            output[i] = input[i] - error[i];
        }
        return output;
    }

    static ZPolynomial extractStandardMessage(ZPolynomial input) {
        auto noError = truncError(input);
        auto decoded = decode(noError);
        auto recoded = encode(input.getZN(), decoded);
        return recoded;
    }

    CSlots toCSlots() const;

    RPolynomial toRPolynomial() const;

    RPolynomial interpretAsRPolynomial() const;

private:
    uint32_t zN;
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

    ZPolynomial interpretAsZPolynomial() const;

private:
    std::vector<BigFixedPoint> coefficients;
};

struct CSlots {
public:
    CSlots(uint32_t zN, uint32_t zSlots) : slots(rN / 2) {}
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
    uint32_t zN;
    uint32_t zSlots;
    std::vector<BigComplex> slots;
};

}  // namespace lbcrypto

#endif  // SRC_PKE_EXAMPLES_Z_ENCODE_UTILS_H_