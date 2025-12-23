//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*

Example for CKKS bootstrapping with full packing

*/

#include "math/chebyshev.h"
#include "openfhe.h"
#include "utils.h"
#include "scheme/ckksrns/z-fhe.h"
#include "math/dftransform-bigcomplex.h"

using namespace lbcrypto;

// related to __heir_debug2
CryptoContextT cc_global;
PublicKeyT pk_global;
PrivateKeyT sk_global;
size_t zN_global;
size_t zSlots_global;

std::vector<BigComplex> zC2SVals;

double __heir_debug2(CiphertextT ct, std::string msg) {
    auto b = DecryptCore(ct->GetElements(), sk_global);

    auto sfBigFP = ct->GetScalingFactorBFP();
    //std::cout << msg << "  Ciphertext Scaling Factor BFP: " << sfBigFP.toHexString() << std::endl;
    auto log2sf = std::log2(sfBigFP.convertToDouble());
    std::cout << msg << "  Scaling factor log2: " << std::setprecision(20) << log2sf << std::endl;
    auto q = ct->GetElements()[0].GetParams()->GetModulus();
    double log2q;
    if (q.GetMSB() > 128) {
        auto offset = q.GetMSB() - 128;
        q >>= offset;
        log2q = std::log2(q.ConvertToDouble()) + offset;
    }
    else {
        log2q = std::log2(q.ConvertToDouble());
    }
    std::cout << msg << "  q: " << log2q << std::endl;
    auto l = ct->GetElements()[0].GetParams()->GetParams().size();
    std::cout << msg << "  l: " << l - 1 << std::endl;

    ZEncodingParams params(ZMode, zN_global, zSlots_global);
    auto zEncode = std::make_shared<ZEncodingImpl>(b.GetParams(), b, sfBigFP, params);

    RPolynomial values = ZEncodingImpl::decodeR(zEncode);

    // This is for sparse LT
    ZEncodingParams paramsTwice(CMode, zN_global * zSlots_global * 2);
    auto zEncodeTwice       = std::make_shared<ZEncodingImpl>(b.GetParams(), b, sfBigFP, paramsTwice);
    RPolynomial valuesTwice = ZEncodingImpl::decodeR(zEncodeTwice);

    auto printZPoly = [&](const ZPolynomial zPoly) {
        auto decoded = ZPolynomial::decode(zPoly);
        std::cout << msg << "  zPoly Decoded: " << std::hex << decoded << std::dec << std::endl;
        auto I = ZPolynomial::extractI(zPoly);
        std::cout << msg << "  zPoly I: ";
        for (size_t i = 0; i != I.getCoefficients().size(); ++i) {
            std::cout << I[i].toHexString(16) << " ";
        }
        std::cout << std::endl;
        std::cout << msg << "  zPoly error log2Norm: " << ZPolynomial::extractError(zPoly).getLog2Norm() << std::endl;
    };

    [[maybe_unused]] auto extractErrorRoly = [&](const RPolynomial rPoly, size_t bits) {
        std::vector<BigFixedPoint> errors;
        auto scalarBFP = BigFixedPoint(BigInteger(1) << bits, 0, false).scaleTo(128);
        for (auto& coeff : rPoly.getCoefficients()) {
            auto integerPart = (coeff * scalarBFP).round() / scalarBFP;
            auto fracPart    = coeff - integerPart;
            errors.push_back(fracPart);
            //std::cout << msg << "  rPoly error: " << fracPart.toHexString(16) << std::endl;
        }
        // Use the ZPolynomial method to get log2Norm
        // Not thinking it is ZPolynomial here
        std::cout << msg << "  rPoly error log2Norm: " << ZPolynomial(errors).getLog2Norm() << std::endl;
    };

    auto printCSlots = [&](const CSlots cSlots) {
        uint64_t reconstructed = 0;
        double log2MaxError    = -1000.0;
        for (size_t i = 0; i != cSlots.getSlots().size(); ++i) {
            auto realPart = cSlots[i].getReal();
            auto fracPart = realPart - realPart.round();
            auto k        = static_cast<int>(std::round(realPart.convertToDouble()));
            reconstructed += (static_cast<uint64_t>(k) << (i));
            log2MaxError = std::max(log2MaxError, fracPart.log2Norm());
        }
        std::cout << msg << "  CSlots Reconstructed: " << std::hex << reconstructed << std::dec
                  << " with error: " << log2MaxError << std::endl;
    };

    enum class DecodeMode { RDecode, CSlotsDecode, CSlotsTwiceDecode, ZDecode };
    std::map<std::string, DecodeMode> decodeMap = {
        // RDecode
        {"ModRaise", DecodeMode::RDecode},
        {"C2R", DecodeMode::RDecode},
        // CSlotsDecode
        {"CSlotsDecode", DecodeMode::CSlotsDecode},
        // CSlotsTwiceDecode
        {"LUT", DecodeMode::CSlotsTwiceDecode},
        {"Normalize", DecodeMode::CSlotsTwiceDecode},
        {"Core", DecodeMode::CSlotsTwiceDecode},
        {"Z2C", DecodeMode::CSlotsTwiceDecode},
        {"Boolean", DecodeMode::CSlotsTwiceDecode},
        {"BooleanAgain", DecodeMode::CSlotsTwiceDecode},
        // ZDecode
        {"Input", DecodeMode::ZDecode},
        {"CMult", DecodeMode::ZDecode},
    };

    auto decodeModeIt = decodeMap.find(msg);
    if (decodeModeIt == decodeMap.end()) {
        return 0;
    }
    auto decodeMode = decodeModeIt->second;

    // Check the r encoding?
    if (decodeMode == DecodeMode::RDecode) {
        for (size_t i = 0; i != values.getCoefficients().size(); ++i) {
            std::cout << msg << "  values [" << i << "]: " << values[i].toHexString(ceil(log2sf / 4.0)) << std::endl;
        }
    }
    if (decodeMode == DecodeMode::CSlotsTwiceDecode) {
        auto cSlotsTwice = valuesTwice.toCSlots();
        for (size_t i = 0; i != cSlotsTwice.getSlots().size(); ++i) {
            std::cout << msg << "  complexValuesTwice [" << i << "]: " << cSlotsTwice[i].toHexString(ceil(log2sf / 4.0))
                      << std::endl;
        }
        printCSlots(cSlotsTwice);
    }
    if (decodeMode == DecodeMode::ZDecode) {
        ZPolynomial zValues = values.toCSlots().getZPolynomial(0);
        for (size_t i = 0; i != values.getCoefficients().size(); ++i) {
            std::cout << msg << "  zValues [" << i << "]: " << zValues[i].toHexString(ceil(log2sf / 4.0)) << std::endl;
        }
        printZPoly(zValues);
    }
    return 0;
}

std::vector<BigComplex> another_interpolate(size_t p, int order = 1) {
    auto f = [&](auto x) {
        return -(p - double(x)) / p;
    };

    auto omega = std::exp(2i * M_PI / double(p));

    std::vector<std::complex<double>> beta;

    for (size_t m = 0; m != p; ++m) {
        std::complex<double> ret = 0;
        for (size_t ell = 0; ell != p; ++ell) {
            ret += double(f(ell)) * std::pow(omega, -double(ell) * m);
        }
        ret /= double(p);
        beta.push_back(ret);
    }
    std::vector<std::complex<double>> alpha;

    if (order == 1) {
        for (size_t m = 0; m != p; ++m) {
            alpha.push_back((1.0 + double(m) / p) * beta[m]);
        }
        for (size_t m = p; m != 2 * p; ++m) {
            alpha.push_back((1.0 - double(m) / p) * beta[m - p]);
        }
    }
    if (order == 2) {
        for (size_t m = 0; m != p; ++m) {
            alpha.push_back((1.0 + double(m) * (double(m) + 3 * p) / (2.0 * p * p)) * beta[m]);
        }
        for (size_t m = p; m != 2 * p; ++m) {
            alpha.push_back((-double(m - p) * (double(m - p) + 2 * p) / (double(p) * p)) * beta[m - p]);
        }
        for (size_t m = 2 * p; m != 3 * p; ++m) {
            alpha.push_back((double(m - p - p) * (double(m - p - p) + p) / (2.0 * double(p) * p)) * beta[m - p - p]);
        }
    }

    std::vector<BigComplex> alphaBigComplex(alpha.size());
    for (size_t i = 0; i != alpha.size(); ++i) {
        alphaBigComplex[i] =
            BigComplex(BigFixedPoint::fromDouble(alpha[i].real()), BigFixedPoint::fromDouble(alpha[i].imag()));
    }
    return alphaBigComplex;
}

void MSBBootstrap(CiphertextT ct, LeveledZ z, AdvancedZ advZ, FHEZ fheZ, uint32_t oneHotBit) {
    // Normalize to [-1, 1] from [-16, 16]
    // This is required by Chebyshev
    // Multiply by 2 * (rN / N / 32) * Delta, so the result is m / 16 * Delta^2
    // NOTE: two here should be changed.
    //auto raisedSF = raised->GetScalingFactorBFP();
    //auto NBigFP   = BigFixedPoint::positive(N);
    //auto two      = BigFixedPoint::two();
    //auto rNBFP    = BigFixedPoint::positive(zN_global * zSlots_global);
    //auto BFP32    = BigFixedPoint::positive(32);
    //// Then C2S below will multiply by rN again because of the construction of U0HatT
    //BigFixedPoint normalizeFactorBFP = raisedSF * two * rNBFP / NBigFP / BFP32 / rNBFP;
    //BigInteger normalizeFactor       = (normalizeFactorBFP.round().getValue()) >> normalizeFactorBFP.getLog2Scale();
    //raised                           = z->EvalMultScalar(raised, normalizeFactor);
    //raised->SetScalingFactorBFP(raisedSF * raisedSF);
    //z->ModReduceInPlace(raised);
    //__heir_debug2(raised, "Normalize");

    // Note that there are other ways...some work first multiply by 1 / N
    // Then PartialSum
    // Then use CoeffsToSlots matrix to do the /16
    // //__heir_debug2(c2s, "C2S");

    //------------------------------------------------------------------------------
    // Running Approximate Mod Reduction
    //------------------------------------------------------------------------------

    //auto& coeff_exp = coeff_exp_16_big_complex_46;
    //auto res        = advZ->EvalChebyshevSeriesPS(c2s, coeff_exp);

    //// Double angle-iterations to get exp(2*Pi*i*x)
    //res = z->EvalMult(res, res);
    //z->ModReduceInPlace(res);
    //res = z->EvalMult(res, res);
    //z->ModReduceInPlace(res);

    //__heir_debug2(res, "Cheby1");

    //------------------------------------------------------------------------------
    // Running LUT
    //------------------------------------------------------------------------------

    //auto lutCoeffs = another_interpolate(1l << 8, 2);
    //auto powers    = advZ->EvalPowers(res, lutCoeffs);
    //auto lut       = advZ->EvalPolyWithPrecomp(powers, lutCoeffs);
    ////__heir_debug2(lut, "LUT");
}

void SimpleBootstrapExample();

int main(int argc, char* argv[]) {
    SimpleBootstrapExample();
}

void SimpleBootstrapExample() {
    CCParams<CryptoContextCKKSRNS> parameters;
    // A. Specify main parameters
    /*  A1) Secret key distribution
    * The secret key distribution for CKKS should either be SPARSE_TERNARY or UNIFORM_TERNARY.
    * The SPARSE_TERNARY distribution was used in the original CKKS paper,
    * but in this example, we use UNIFORM_TERNARY because this is included in the homomorphic
    * encryption standard.
    */
    SecretKeyDist secretKeyDist = lbcrypto::SPARSE_ENCAPSULATED;
    parameters.SetSecretKeyDist(secretKeyDist);

    /*  A2) Desired security level based on FHE standards.
    * In this example, we use the "NotSet" option, so the example can run more quickly with
    * a smaller ring dimension. Note that this should be used only in
    * non-production environments, or by experts who understand the security
    * implications of their choices. In production-like environments, we recommend using
    * HEStd_128_classic, HEStd_192_classic, or HEStd_256_classic for 128-bit, 192-bit,
    * or 256-bit security, respectively. If you choose one of these as your security level,
    * you do not need to set the ring dimension.
    */
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 10);
    //parameters.SetNumLargeDigits(6);

    ScalingTechnique rescaleTech = FLEXIBLEMANUAL;
    uint32_t dcrtBits            = 59;
    uint32_t firstMod            = 59;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    /*  A4) Multiplicative depth.
    * The goal of bootstrapping is to increase the number of available levels we have, or in other words,
    * to dynamically increase the multiplicative depth. However, the bootstrapping procedure itself
    * needs to consume a few levels to run. We compute the number of bootstrapping levels required
    * using GetBootstrapDepth, and add it to levelsAvailableAfterBootstrap to set our initial multiplicative
    * depth. We recommend using the input parameters below to get started.
    */
    std::vector<uint32_t> levelBudget = {2, 2};

    // Note that the actual number of levels avalailable after bootstrapping before next bootstrapping
    // will be levelsAvailableAfterBootstrap - 1 because an additional level
    // is used for scaling the ciphertext before next bootstrapping (in 64-bit CKKS bootstrapping)
    //uint32_t levelsAvailableAfterBootstrap = 10;
    //uint32_t depth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(20);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    //cc->Enable(ADVANCEDSHE);
    //cc->Enable(FHE);

    uint32_t ringDim = cc->GetRingDimension();
    // This is the maximum number of slots that can be used for full packing.
    //uint32_t numSlots = ringDim / 2;
    std::cout << "CKKS scheme ring dimension: " << ringDim << "\n\n";

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::vector<int> rotateIndices = {};
    for (int i = 1; i <= 31; ++i) {
        rotateIndices.push_back(i);
    }
    for (int i = 1; i <= 31; ++i) {
        rotateIndices.push_back(-i);
    }
    const uint32_t limit = ringDim / 32;
    for (uint32_t j = 1; j < limit; j <<= 1) {
        rotateIndices.push_back(j * 16);
    }
    cc->EvalRotateKeyGen(keyPair.secretKey, rotateIndices);
    // 2 * N - 1 for conjugdate
    cc->EvalAutomorphismKeyGen(keyPair.secretKey, {2 * ringDim - 1});
    EvalSparseEncapsulatedKeyGen(keyPair.secretKey);
    //cc->EvalBootstrapKeyGen(keyPair.secretKey, 16);

    std::cout << *(cc->GetCryptoParameters()) << std::endl;

    std::cout << *(std::static_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters())->GetParamsP())
              << " primes in the special prime modulus." << std::endl;

    LeveledZ z     = std::make_shared<LeveledZImpl>();
    AdvancedZ advZ = std::make_shared<AdvancedZImpl>(z);
    FHEZ fheZ      = std::make_shared<FHEZImpl>(z, advZ);

    uint32_t zN = 32;
    uint32_t N  = cc->GetCyclotomicOrder();
    //uint32_t zSlots = N / zN / 2;  // Maximal sparse packing
    uint32_t zSlots = 1;  // Maximal sparse packing
    zN_global       = zN;
    zSlots_global   = zSlots;

    auto cSlots  = zN * zSlots / 2;
    auto cSlots2 = cSlots * 2;
    // For regular encoding
    DiscreteFourierTransformBigComplex::Initialize(cSlots * 4, cSlots);
    // For encoding of bootstrapping related plaintext for sparse bootstrapping
    DiscreteFourierTransformBigComplex::Initialize(cSlots2 * 4, cSlots2);
    ZLinearTransform::Initialize(zN);

    fheZ->EvalBootstrapSetup(*cc, zN * zSlots / 2, levelBudget);

    cc_global = cc;
    pk_global = keyPair.publicKey;
    sk_global = keyPair.secretKey;

    auto zero      = EncryptZero(keyPair.publicKey);
    auto elemParam = zero->GetElements()[0].GetParams();
    auto sfq0      = BigFixedPoint::positive(elemParam->GetParams()[0]->GetModulus().ConvertToInt());

    {
        auto sfNow = sfq0;
        auto sfMax = sfNow;
        auto sfMin = sfNow;
        for (size_t i = elemParam->GetParams().size() - 1; i != size_t(-1); --i) {
            auto qi = BigFixedPoint::positive(elemParam->GetParams()[i]->GetModulus().ConvertToInt());
            sfNow   = sfNow * sfNow / qi;
            std::cout << "Level " << (elemParam->GetParams().size() - 1 - i)
                      << " Scaling Factor log2: " << std::log2(sfNow.convertToDouble()) << std::endl;
        }
    }

    //
    //__heir_debug2(zero, "Input");

    auto zPoly = ZPolynomial::encode(32, 0x1);
    // add some noise
    //for (size_t i = 0; i != zPoly.getCoefficients().size(); ++i) {
    //    zPoly[i] += BigFixedPoint::positive(i + 1) / BigFixedPoint::positive(1 << 25);
    //}
    auto singleCSlots = zPoly.toCSlots();
    // Multiply by N/(2 * n) for maximal sparse packing
    ZEncodingParams paramsMaximalSparse(ZMode, zN, zSlots);
    std::vector<BigComplex> scaledSlots;
    for (size_t i = 0; i != singleCSlots.getSlots().size(); ++i) {
        scaledSlots.push_back(singleCSlots[i] * BigFixedPoint::positive(N) / BigFixedPoint::positive(2 * zN));
    }
    std::vector<BigComplex> maximalSlots;
    for (size_t i = 0; i != zSlots; ++i) {
        for (size_t j = 0; j != singleCSlots.getSlots().size(); ++j) {
            maximalSlots.push_back(scaledSlots[j]);
        }
    }
    CSlots maximalCSlots(paramsMaximalSparse, maximalSlots);
    RPolynomial scaledRPoly = maximalCSlots.toRPolynomial();

    RPolynomial value1 = zPoly.toCSlots().toRPolynomial();
    Plaintext ptxt1    = ZEncodingImpl::encodeR(value1, elemParam, sfq0);

    RPolynomial value2 = ZPolynomial::encode(32, -1).toCSlots().toRPolynomial();
    Plaintext ptxt2    = ZEncodingImpl::encodeR(value2, elemParam, sfq0);

    /// TEST ENCODE
    auto encoded  = Encrypt(ptxt1, keyPair.publicKey);
    auto encoded2 = Encrypt(ptxt2, keyPair.publicKey);

    //__heir_debug2(encoded, "Encode");

    /// TEST ADD
    if (0) {
        auto ctAdd = z->EvalAdd(encoded, ptxt2);

        __heir_debug2(ctAdd, "Add");
    }

    /// TEST CT-PT-MULT
    RPolynomial t   = ZPolynomial::getT(32).toCSlots().toRPolynomial();
    Plaintext tPtxt = ZEncodingImpl::encodeR(t, elemParam, sfq0);

    if (0) {
        auto ctMul = z->EvalMult(encoded, tPtxt);
        __heir_debug2(ctMul, "Mult");
        z->ModReduceInPlace(ctMul);
        __heir_debug2(ctMul, "Mult");
    }

    /// TEST CT-CT-MULT
    Ciphertext<DCRTPoly> ct = encoded;
    if (0) {
        auto ctMul = zEvalMultFull(z, encoded, encoded2, tPtxt);
        z->ModReduceInPlace(ctMul, 2);
        __heir_debug2(ctMul, "CMult");
        ct = ctMul;
    }
    __heir_debug2(ct, "Input");

    if (0) {
        ct = fheZ->EvalArithToArith(ct);
        __heir_debug2(ct, "CMult");
    }

    Ciphertext<DCRTPoly> ct2;
    for (size_t i = 0; i != 1; ++i) {
        ct2 = fheZ->EvalArithToBoolean(ct);
    }
    __heir_debug2(ct2, "Boolean");
    std::cout << "Finished Warmup\n";
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i != 3; ++i) {
        ct2 = fheZ->EvalArithToBoolean(ct);
    }
    auto end      = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Average time for ArithToBoolean: " << double(duration) / 3.0 << " ms\n";

    Ciphertext<DCRTPoly> ct3;
    for (size_t i = 0; i != 1; ++i) {
        ct3 = fheZ->EvalArithToArith(ct);
    }
    __heir_debug2(ct3, "CMult");
    start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i != 3; ++i) {
        ct3 = fheZ->EvalArithToArith(ct);
    }
    end      = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Average time for ArithToArith: " << double(duration) / 3.0 << " ms\n";

    for (size_t i = 0; i != 1; ++i) {
        ct = fheZ->EvalBooleanToBoolean(ct2);
    }
    __heir_debug2(ct, "BooleanAgain");
    start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i != 3; ++i) {
        ct = fheZ->EvalBooleanToBoolean(ct2);
    }
    end      = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Average time for BooleanToBoolean: " << double(duration) / 3.0 << " ms\n";

    /// TEST Rotate
    if (0) {
        auto ctRot = cc->EvalRotate(encoded, 1);

        __heir_debug2(ctRot, "Rotate");
    }
}
