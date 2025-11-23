#include <cassert>
#include "openfhe.h"

using namespace lbcrypto;
using CiphertextT        = ConstCiphertext<DCRTPoly>;
using MutableCiphertextT = Ciphertext<DCRTPoly>;
using CCParamsT          = CCParams<CryptoContextBFVRNS>;
using CryptoContextT     = CryptoContext<DCRTPoly>;
using EvalKeyT           = EvalKey<DCRTPoly>;
using PlaintextT         = Plaintext;
using PrivateKeyT        = PrivateKey<DCRTPoly>;
using PublicKeyT         = PublicKey<DCRTPoly>;

// DecryptCore not accessible from CryptoContext
// so copy from @openfhe//src/pke/lib/schemerns/rns-pke.cpp
DCRTPoly DecryptCore(const std::vector<DCRTPoly>& cv, const PrivateKey<DCRTPoly> privateKey) {
    const DCRTPoly& s = privateKey->GetPrivateElement();

    size_t sizeQ  = s.GetParams()->GetParams().size();
    size_t sizeQl = cv[0].GetParams()->GetParams().size();

    size_t diffQl = sizeQ - sizeQl;

    auto scopy(s);
    scopy.DropLastElements(diffQl);

    DCRTPoly sPower(scopy);

    DCRTPoly b(cv[0]);
    b.SetFormat(Format::EVALUATION);

    DCRTPoly ci;
    for (size_t i = 1; i < cv.size(); i++) {
        ci = cv[i];
        ci.SetFormat(Format::EVALUATION);

        b += sPower * ci;
        sPower *= scopy;
    }
    return b;
}

std::vector<std::complex<double>> multiplyByUComplex(std::vector<double> input) {
    // now input.size() == slots * 2
    // as it is the special input polynomial
    std::vector<std::complex<double>> val(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        val[i].real(input[i]);
        val[i].imag(0);
    }
    size_t slots   = input.size() / 2;
    uint32_t m     = 4 * slots;
    uint32_t mmask = m - 1;  // assumes m is power of 2

    // computes indices for all primitive roots of unity
    std::vector<uint32_t> rotGroup(slots);
    uint32_t fivePows = 1;
    for (uint32_t i = 0; i < slots; ++i) {
        rotGroup[i] = fivePows;
        fivePows *= 5;
        fivePows &= mmask;
    }

    // computes all powers of a primitive root of unity exp(2 * M_PI/m)
    std::vector<std::complex<double>> ksiPows(m + 1);
    double ak = 2 * M_PI / m;
    for (uint32_t j = 0; j < m; ++j) {
        double angle = ak * j;
        ksiPows[j].real(std::cos(angle));
        ksiPows[j].imag(std::sin(angle));
    }
    ksiPows[m] = ksiPows[0];

    std::vector<std::vector<std::complex<double>>> U(slots, std::vector<std::complex<double>>(2 * slots));

    for (uint32_t i = 0; i < slots; ++i) {
        for (uint32_t j = 0; j < 2 * slots; ++j) {
            U[i][j] = ksiPows[(j * rotGroup[i]) & mmask];
        }
    }

    std::vector<std::complex<double>> result(slots);

    // do matrix multiplication of U * val
    for (uint32_t i = 0; i < slots; ++i) {
        result[i] = 0;
        for (uint32_t j = 0; j < 2 * slots; ++j) {
            result[i] += U[i][j] * val[j];
        }
    }
    return result;
}

std::vector<double> multiplyByU(std::vector<double> input) {
    auto result = multiplyByUComplex(input);
    std::vector<double> realResult;
    for (size_t i = 0; i < result.size(); ++i) {
        realResult.push_back(result[i].real());
    }
    return realResult;
}

// For X^8-X+2
static const inline std::vector<std::complex<double>> z_upper_roots_8 = {
    -1.0553707970856263 + 0.46044989545852666j, -0.3759904375006573 + 1.062443306305578j,
    0.4809269633899156 + 0.9626774673864164j, 0.9504342711963706 + 0.3496349308325642j};

const auto z_upper_roots = z_upper_roots_8;
const auto zN            = 8;

std::vector<std::complex<double>> multiplyByZUComplex(std::vector<double> input) {
    std::vector<std::complex<double>> val(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        val[i].real(input[i]);
        val[i].imag(0);
    }
    assert(input.size() == z_upper_roots.size() * 2 &&
           "Input size does not match expected size for multiplyByZUComplex");

    // This is the vandermond matrix
    std::vector<std::complex<double>> result;
    size_t halfSize = input.size() / 2;
    for (size_t i = 0; i < halfSize; ++i) {
        std::complex<double> sum = 0;
        for (size_t j = 0; j < input.size(); ++j) {
            sum += std::pow(z_upper_roots[i], static_cast<int>(j)) * val[j];
        }
        result.push_back(sum);
    }
    return result;
}

std::vector<double> multiplyByZUInverse(std::vector<std::complex<double>> input) {
    std::vector<double> result;
    assert(input.size() == z_upper_roots.size() && "Input size does not match expected size for multiplyByZUInverse");

    std::vector<std::vector<std::complex<double>>> zUInv(zN, std::vector<std::complex<double>>(zN / 2, 0));

    // Build the inverse Vandermond matrix
    for (size_t i = 0; i != zN; ++i) {
        for (size_t j = 0; j != zN / 2; ++j) {
            auto xi = z_upper_roots[j];
            // d = zN * xi^(zN-1) - 1
            auto d = (static_cast<double>(zN) * std::pow(xi, static_cast<int>(zN - 1))) - 1.0;
            if (i == 0) {
                zUInv[i][j] = (std::pow(xi, static_cast<int>(zN - 1)) - 1.0) / d;
            }
            else {
                zUInv[i][j] = std::pow(xi, static_cast<int>(zN - 1 - i)) / d;
            }
        }
    }

    for (size_t i = 0; i != zN; ++i) {
        std::complex<double> sum = 0;
        for (size_t j = 0; j != zN / 2; ++j) {
            sum += zUInv[i][j] * input[j];
        }
        // z + conj(z) = 2*real(z)
        result.push_back(2 * sum.real());
    }
    return result;
}

void test_zu() {
    std::vector<double> input = {1, 2, 3, 4, 5, 6, 7, 8};
    auto zuResult             = multiplyByZUComplex(input);
    auto zuInvResult          = multiplyByZUInverse(zuResult);

    std::cout << "Input: ";
    for (const auto& val : input) {
        std::cout << std::setprecision(20) << val << " ";
    }
    std::cout << "\n";

    std::cout << "After ZU and inverse ZU: ";
    for (const auto& val : zuInvResult) {
        std::cout << std::setprecision(20) << val << " ";
    }
    std::cout << "\n";
}

std::vector<double> encodeInRE(int32_t input) {
    std::vector<double> bits;
    // get bits of input in bits
    for (size_t i = 0; i != zN; ++i) {
        bits.push_back((input & (1 << i)) >> i);
    }
    std::vector<double> tInv;
    // denominator
    double d = std::pow(2.0, static_cast<int>(zN));
    for (size_t i = 0; i != zN; ++i) {
        if (i == 0) {
            tInv.push_back((1.0 - std::pow(2.0, static_cast<int>(zN - 1))) / d);
        }
        else {
            tInv.push_back(-std::pow(2.0, static_cast<int>(zN - 1 - i)) / d);
        }
    }
    // multiply bits by tInv
    std::vector<double> result(2 * zN - 1, 0.0);
    for (size_t i = 0; i != zN; ++i) {
        for (size_t j = 0; j != zN; ++j) {
            result[i + j] += bits[i] * tInv[j];
        }
    }
    // now euclidean reduction mod X^zN - X + 2
    for (size_t i = result.size() - 1; i >= zN; --i) {
        result[i - zN + 1] += result[i];
        result[i - zN] += -2.0 * result[i];
        result.pop_back();
    }
    return result;
}

int32_t decodeFromRE(const std::vector<double>& input) {
    std::vector<double> t(zN, 0);
    t[0] = -2.0;
    t[1] = 1.0;

    std::vector<double> result(2 * zN - 1, 0);
    for (size_t i = 0; i != zN; ++i) {
        for (size_t j = 0; j != zN; ++j) {
            result[i + j] += input[i] * t[j];
        }
    }
    // now euclidean reduction mod X^zN - X + 2
    for (size_t i = result.size() - 1; i >= zN; --i) {
        result[i - zN + 1] += result[i];
        result[i - zN] += -2.0 * result[i];
        result.pop_back();
    }
    int32_t ret = 0;
    for (size_t i = 0; i != zN; ++i) {
        ret += static_cast<int>(std::round(result[i])) << i;
    }
    return ret;
}

void test_encodeInRE() {
    auto input   = 255;
    auto encoded = encodeInRE(input);
    std::cout << "Input: " << input << "\nEncoded in RE: ";
    for (const auto& val : encoded) {
        std::cout << std::setprecision(20) << val << " ";
    }
    auto decoded = decodeFromRE(encoded);
    std::cout << "Decoded: " << decoded;
    std::cout << "\n";
}

void test2_encodeInRE() {
    auto input1   = 1;
    auto input2   = 255;
    auto encoded1 = encodeInRE(input1);
    auto encoded2 = encodeInRE(input2);
    std::vector<double> sum(encoded1.size(), 0.0);
    for (size_t i = 0; i != encoded1.size(); ++i) {
        sum[i] = encoded1[i] + encoded2[i];
    }
    std::cout << "Sum 255 + 1 in RE: ";
    for (const auto& val : sum) {
        std::cout << std::setprecision(20) << val << " ";
    }
}

void test3_encodeInRE() {
    auto input1     = 33;
    auto input2     = 57;
    auto encoded1   = encodeInRE(input1);
    auto encoded2   = encodeInRE(input2);
    auto multResult = std::vector<double>(encoded1.size() + encoded2.size() - 1, 0.0);
    for (size_t i = 0; i != encoded1.size(); ++i) {
        for (size_t j = 0; j != encoded2.size(); ++j) {
            multResult[i + j] += encoded1[i] * encoded2[j];
        }
    }
    // now euclidean reduction mod X^zN - X + 2
    for (size_t i = multResult.size() - 1; i >= zN; --i) {
        multResult[i - zN + 1] += multResult[i];
        multResult[i - zN] += -2.0 * multResult[i];
        multResult.pop_back();
    }
    // then muliply by t
    std::vector<double> t(zN, 0);
    t[0] = -2.0;
    t[1] = 1.0;
    std::vector<double> finalResult(zN + zN - 1, 0.0);
    for (size_t i = 0; i != zN; ++i) {
        for (size_t j = 0; j != zN; ++j) {
            finalResult[i + j] += multResult[i] * t[j];
        }
    }
    // now euclidean reduction mod X^zN - X + 2
    for (size_t i = finalResult.size() - 1; i >= zN; --i) {
        finalResult[i - zN + 1] += finalResult[i];
        finalResult[i - zN] += -2.0 * finalResult[i];
        finalResult.pop_back();
    }
    auto decoded = decodeFromRE(finalResult);
    std::cout << "Decoded 33 * 57 from RE: " << decoded << "\n";
}