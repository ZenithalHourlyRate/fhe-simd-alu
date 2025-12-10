#include "math/z-constants.h"

namespace lbcrypto {

BigCMatrix getZU() {
    // Vandermond matrix
    BigCMatrix zu(zN / 2, std::vector<BigComplex>(zN));
    for (size_t i = 0; i != zN / 2; ++i) {
        zu[i][0] = BigFixedPoint::one();
        for (size_t j = 1; j != zN; ++j) {
            zu[i][j] = zu[i][j - 1] * z_upper_roots[i];
        }
    }
    return zu;
}

BigCMatrix getZUInverse() {
    std::vector<std::vector<BigComplex>> zUInv(zN, std::vector<BigComplex>(zN / 2));

    // Build the inverse Vandermond matrix
    for (size_t j = 0; j != zN / 2; ++j) {
        auto xi = z_upper_roots[j];
        std::vector<BigComplex> xiPowers;
        auto one   = BigFixedPoint::one();
        auto zNbig = BigFixedPoint::positive(zN);
        xiPowers.push_back(one);
        for (size_t p = 1; p != zN; ++p) {
            xiPowers.push_back(xiPowers[p - 1] * xi);
        }
        // build power map
        for (size_t i = 0; i != zN; ++i) {
            // d = zN * xi^(zN-1) - 1
            auto d = (zNbig * xiPowers[zN - 1]) - one;
            if (i == 0) {
                zUInv[i][j] = (xiPowers[zN - 1] - one) / d;
            }
            else {
                zUInv[i][j] = xiPowers[zN - 1 - i] / d;
            }
        }
    }
    return zUInv;
}

BigCMatrix GetZU(uint32_t zN) {
    // Vandermond matrix
    BigCMatrix zu(zN / 2, std::vector<BigComplex>(zN));
    for (size_t i = 0; i != zN / 2; ++i) {
        zu[i][0] = BigFixedPoint::one();
        for (size_t j = 1; j != zN; ++j) {
            zu[i][j] = zu[i][j - 1] * Z_ROOTS_MAP.at(zN)[i];
        }
    }
    return zu;
}

BigCMatrix GetZUInverse(uint32_t zN) {
    std::vector<std::vector<BigComplex>> zUInv(zN, std::vector<BigComplex>(zN / 2));

    // Build the inverse Vandermond matrix
    for (size_t j = 0; j != zN / 2; ++j) {
        auto xi = Z_ROOTS_MAP.at(zN)[j];
        std::vector<BigComplex> xiPowers;
        auto one   = BigFixedPoint::one();
        auto zNbig = BigFixedPoint::positive(zN);
        xiPowers.push_back(one);
        for (size_t p = 1; p != zN; ++p) {
            xiPowers.push_back(xiPowers[p - 1] * xi);
        }
        // build power map
        for (size_t i = 0; i != zN; ++i) {
            // d = zN * xi^(zN-1) - 1
            auto d = (zNbig * xiPowers[zN - 1]) - one;
            if (i == 0) {
                zUInv[i][j] = (xiPowers[zN - 1] - one) / d;
            }
            else {
                zUInv[i][j] = xiPowers[zN - 1 - i] / d;
            }
        }
    }
    return zUInv;
}

std::vector<BigComplex> multU(const BigCMatrix& U, std::vector<BigFixedPoint> input) {
    assert(input.size() == U[0].size() && "Input size does not match expected size for U");

    // This is the vandermond matrix
    std::vector<BigComplex> result;

    size_t halfSize = input.size() / 2;
    for (size_t i = 0; i < halfSize; ++i) {
        BigComplex sum;
        for (size_t j = 0; j < input.size(); ++j) {
            sum = sum + U[i][j] * input[j];
        }
        result.push_back(sum);
    }
    return result;
}

std::vector<BigFixedPoint> multUInverse(const BigCMatrix& UInv, std::vector<BigComplex> input) {
    std::vector<BigFixedPoint> result;
    assert(input.size() == z_upper_roots.size() && "Input size does not match expected size for multiplyByZUInverse");

    for (size_t i = 0; i != zN; ++i) {
        BigComplex sum;
        for (size_t j = 0; j != zN / 2; ++j) {
            sum = sum + UInv[i][j] * input[j];
        }
        // z + conj(z) = 2*real(z)
        auto two = BigFixedPoint::positive(2);
        result.push_back(two * sum.getReal());
    }
    return result;
}

BigCMatrix getRU() {
    std::vector<std::vector<BigComplex>> U(rN / 2, std::vector<BigComplex>(rN));
    for (uint32_t i = 0; i < rN / 2; ++i) {
        auto zeta = r_roots[i];
        // build power map
        std::vector<BigComplex> zetaPows;
        auto one = BigFixedPoint::one();
        zetaPows.push_back(one);
        for (uint32_t p = 1; p != rN; ++p) {
            zetaPows.push_back(zetaPows[p - 1] * zeta);
        }

        for (uint32_t j = 0; j < rN; ++j) {
            U[i][j] = zetaPows[j];
        }
    }
    return U;
}

BigCMatrix getRUInverse() {
    auto U = getRU();
    BigCMatrix UInverse(rN, std::vector<BigComplex>(rN / 2));
    auto rNBig = BigFixedPoint::positive(rN);
    for (size_t i = 0; i != rN / 2; ++i) {
        for (size_t j = 0; j != rN; ++j) {
            UInverse[j][i] = U[i][j].conj() / rNBig;
        }
    }
    return UInverse;
}

}  // namespace lbcrypto