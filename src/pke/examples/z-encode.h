#ifndef SRC_PKE_EXAMPLES_Z_ENCODE_H_
#define SRC_PKE_EXAMPLES_Z_ENCODE_H_

#include <cassert>
#include "openfhe.h"
#include "high-prec-complex.h"

std::vector<BigFixedPoint> encodeInZ(uint32_t input) {
    std::vector<double> bits;
    // get bits of input in bits
    for (size_t i = 0; i != zN; ++i) {
        bits.push_back((input & (1 << i)) >> i);
    }

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
    // multiply bits by tInv
    std::vector<BigFixedPoint> result(2 * zN - 1, BigFixedPoint(0));
    for (size_t i = 0; i != zN; ++i) {
        for (size_t j = 0; j != zN; ++j) {
            result[i + j] += bits[i] * tInv[j];
        }
    }
    auto two = BigFixedPoint(2, 0, false).scaleTo(z_upper_roots_scale);
    // now euclidean reduction mod X^zN - X + 2
    for (size_t i = result.size() - 1; i >= zN; --i) {
        result[i - zN + 1] += result[i];
        result[i - zN] += -two * result[i];
        result.pop_back();
    }
    return result;
}

std::vector<BigFixedPoint> roundInZ(std::vector<BigFixedPoint> input) {
    // Add small epsilon to ensure range in (-1, 0)
    // epsilon = 2^{-zN-1}
    auto eps = BigFixedPoint(BigInteger(1) << z_upper_roots_scale - zN - 1, z_upper_roots_scale, false);
    // do [\cdot ]_1
    std::vector<BigFixedPoint> output;
    for (auto i : input) {
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

uint32_t decodeFromZ(const std::vector<BigFixedPoint>& input) {
    std::vector<BigFixedPoint> t(zN, 0);
    auto two = BigFixedPoint(2, 0, false).scaleTo(z_upper_roots_scale);
    auto one = BigFixedPoint(1, 0, false).scaleTo(z_upper_roots_scale);
    t[0]     = -two;
    t[1]     = one;

    std::vector<BigFixedPoint> result(2 * zN - 1, 0);
    for (size_t i = 0; i != zN; ++i) {
        for (size_t j = 0; j != zN; ++j) {
            result[i + j] += input[i] * t[j];
        }
    }
    // now euclidean reduction mod X^zN - X + 2
    for (size_t i = result.size() - 1; i >= zN; --i) {
        result[i - zN + 1] += result[i];
        result[i - zN] += -two * result[i];
        result.pop_back();
    }
    // Now we have [m]_t + (X-2)I with I in {0, -1}
    // Need to remove the possible I
    // We need to add eps = 3 * 2^{-zN - 1}
    auto eps               = BigFixedPoint(BigInteger(3) << z_upper_roots_scale - zN - 1, z_upper_roots_scale, false);
    auto constCoeff        = result[0] + eps;
    auto constCoeffInteger = constCoeff.getValue() >> constCoeff.getLog2Scale();
    // This is b1. OpenFHE count bits from 1...
    int32_t I = constCoeffInteger.GetBitAtIndex(2);
    auto BigI = BigFixedPoint(I, 0, false);
    result[1] += BigI;
    // Now we get pure [m]_t

    int32_t ret = 0;
    for (size_t i = 0; i != zN; ++i) {
        // actually we should use rounding...
        // But for convenience let's just make it larger than 0
        // It becomes either 0 + small or 1 + small
        result[i] += eps;
        auto iInt = result[i].getValue() >> result[i].getLog2Scale();
        // OpenFHE count from 1???
        int64_t iBit = iInt.GetBitAtIndex(1);
        //std::cout << "iInt " << result[i].toHexString() << " Bit " << i << "\n";
        ret += iBit << i;
    }
    return ret;
}

std::vector<BigFixedPoint> multiplyInRE(const std::vector<BigFixedPoint>& lhs, const std::vector<BigFixedPoint>& rhs) {
    auto multResult = std::vector<BigFixedPoint>(lhs.size() + lhs.size() - 1, 0);
    for (size_t i = 0; i != lhs.size(); ++i) {
        for (size_t j = 0; j != rhs.size(); ++j) {
            multResult[i + j] += lhs[i] * rhs[j];
        }
    }
    auto one = BigFixedPoint(1, 0, false);
    auto two = BigFixedPoint(2, 0, false);
    // now euclidean reduction mod X^zN - X + 2
    for (size_t i = multResult.size() - 1; i >= zN; --i) {
        multResult[i - zN + 1] += multResult[i];
        multResult[i - zN] += -two * multResult[i];
        multResult.pop_back();
    }
    // then muliply by t
    std::vector<BigFixedPoint> t(zN, 0);
    t[0] = -two;
    t[1] = one;
    std::vector<BigFixedPoint> finalResult(zN + zN - 1, 0);
    for (size_t i = 0; i != zN; ++i) {
        for (size_t j = 0; j != zN; ++j) {
            finalResult[i + j] += multResult[i] * t[j];
        }
    }
    // now euclidean reduction mod X^zN - X + 2
    for (size_t i = finalResult.size() - 1; i >= zN; --i) {
        finalResult[i - zN + 1] += finalResult[i];
        finalResult[i - zN] += -two * finalResult[i];
        finalResult.pop_back();
    }
    return finalResult;
}

#endif  // SRC_PKE_EXAMPLES_Z_ENCODE_H_