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
    // do [\cdot ]_1
    std::vector<BigFixedPoint> output;
    for (auto i : input) {
        if (i.getNeg()) {
            // make it positive
            auto ni      = -i;
            auto niFloor = (ni.getValue() >> ni.getLog2Scale());
            auto niFP    = BigFixedPoint(niFloor, 0, false);
            //std::cout << "ni: " << ni.toBinary() << " niFP: " << niFP.toBinary() << "\n";
            output.push_back(-(ni - niFP));
        }
        else {
            auto iFloor = (i.getValue() >> i.getLog2Scale());
            auto iFP    = BigFixedPoint(iFloor, 0, false);
            output.push_back(i - iFP);
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
    int32_t ret = 0;
    for (size_t i = 0; i != zN; ++i) {
        auto iInt = result[i].getValue() >> result[i].getLog2Scale();
        // OpenFHE count from 1???
        int64_t iBit = iInt.GetBitAtIndex(1);
        bool iNeg    = result[i].getNeg();
        if (iNeg && iBit) {
            iBit = 0;
        }
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