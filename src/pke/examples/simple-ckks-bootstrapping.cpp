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
#include "z-utils.h"

using namespace lbcrypto;

// related to __heir_debug2
CryptoContextT cc_global;
PublicKeyT pk_global;
PrivateKeyT sk_global;
size_t slots_global;

std::vector<BigComplex> zC2SVals;

double __heir_debug2(CiphertextT ct, std::string msg) {
    auto b = DecryptCore(ct->GetElements(), sk_global);

    auto sfBigFP = ct->GetScalingFactorBFP();
    std::cout << msg << "  Ciphertext Scaling Factor BFP: " << sfBigFP.toHexString() << std::endl;
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

    // valueSize = zN
    auto valueSize = slots_global * 2;

    auto zEncode = std::make_shared<ZEncodingImpl>(b.GetParams(), b, valueSize, sfBigFP);

    RPolynomial values = ZEncodingImpl::decodeR(zEncode);

    auto printZPoly = [&](const ZPolynomial zPoly) {
        auto decoded = ZPolynomial::decode(zPoly);
        std::cout << msg << "  zPoly Decoded: " << decoded << std::endl;
        auto I = ZPolynomial::extractI(zPoly);
        std::cout << msg << "  zPoly I: ";
        for (size_t i = 0; i != I.getCoefficients().size(); ++i) {
            std::cout << I[i].toHexString(16) << " ";
        }
        std::cout << std::endl;
        std::cout << msg << "  zPoly error log2Norm: " << ZPolynomial::extractError(zPoly).getLog2Norm() << std::endl;
    };

    auto extractErrorRoly = [&](const RPolynomial rPoly, size_t bits) {
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

    // Check the r encoding?
    if (msg == "S2RC" || msg == "MSB" || msg == "Low4" || msg == "ModRaise" || msg == "PSum" || msg == "Normalize" ||
        msg == "LUTR" || msg.find("MSB") != std::string::npos || msg == "Reconstructed") {
        for (size_t i = 0; i != values.getCoefficients().size(); ++i) {
            std::cout << msg << "  values [" << i << "]: " << values[i].toHexString(16) << std::endl;
        }
        if (msg == "S2RC" || msg == "MSB" || msg == "Reconstructed") {
            // Directly interpret as ZPolynomial
            ZPolynomial zPoly(values.getCoefficients());
            printZPoly(zPoly);
        }
        if (msg == "ModRaise") {
            std::cout << msg << "  values [" << 0 << "]: " << values[0].toString(24) << std::endl;
        }
        if (msg == "Low4" || msg == "ModRaise") {
            extractErrorRoly(values, 4);
        }
        if (msg == "Normalize") {
            extractErrorRoly(values, 8);
        }
    }
    auto cSlots = values.toCSlots();
    if (msg == "Upper" || msg == "Down" || msg == "Rotate" || msg == "MSBC2S" || msg == "LUT") {
        for (size_t i = 0; i != cSlots.getSlots().size(); ++i) {
            std::cout << msg << "  complexValues [" << i << "]: " << cSlots[i].toHexString(16) << std::endl;
        }
        if (msg == "Upper") {
            zC2SVals.clear();
            zC2SVals = cSlots.getSlots();
        }
        if (msg == "Down") {
            for (size_t i = 0; i != cSlots.getSlots().size(); ++i) {
                zC2SVals.push_back(cSlots[i]);
            }
            // Now construct ZPolynomial
            std::vector<BigFixedPoint> zCoeffs(zC2SVals.size());
            for (size_t i = 0; i != zC2SVals.size(); ++i) {
                zCoeffs[i] = zC2SVals[i].getReal();
            }
            printZPoly(zCoeffs);
        }
        if (msg == "LUT") {
            RPolynomial rPoly;
            for (size_t i = 0; i != cSlots.getSlots().size(); ++i) {
                rPoly[i] = cSlots[i].getReal();
            }
            extractErrorRoly(rPoly, 4);
        }
    }
    ZPolynomial zValues = values.toZPolynomial();
    if (msg == "Encode" || msg == "Input" || msg == "Add" || msg == "Mult" || msg == "CMult" || msg == "Z2S2Z") {
        for (size_t i = 0; i != values.getCoefficients().size(); ++i) {
            std::cout << msg << "  zValues [" << i << "]: " << zValues[i].toHexString(16) << std::endl;
        }
        printZPoly(zValues);
    }
#if 0
    if (msg == "Encode" && false) {
        std::vector<BigComplex> encodeZCoeffsToSlots(16);
        auto zUInverse = getZUInverse();

        unsigned halfSize = 16;
        std::vector<std::vector<BigComplex>> results1;
        std::vector<std::vector<BigComplex>> results2;

        for (size_t i = 0; i != halfSize; ++i) {
            auto diagonal = std::vector<BigComplex>(halfSize, BigComplex());
            for (size_t j = 0; j != halfSize; ++j) {
                diagonal[j] = zUInverse[j][(i + j) % halfSize];
            }
            results1.push_back(diagonal);
        }
        for (size_t i = 0; i != halfSize; ++i) {
            auto diagonal = std::vector<BigComplex>(halfSize, BigComplex());
            for (size_t j = 0; j != halfSize; ++j) {
                diagonal[j] = zUInverse[j + 16][(i + j) % halfSize];
            }
            results2.push_back(diagonal);
        }

        //auto ptxt1 = results1[0];
        //for (size_t i = 0; i != ptxt1.size(); ++i) {
        //    std::cout << msg << "  results1[0] [" << i << "]: " << ptxt1[i].toHexString(32) << std::endl;
        //}

        auto elementParams = cc_global->GetCryptoParameters()->GetElementParams();
        //auto scalingFactor = BigFixedPoint(BigInteger(1) << static_cast<uint32_t>(std::log2(sf)), 0, false);
        auto auxPtxts      = getZCoeffToSlotsAuxDCRTPoly(elementParams, scalingFactor);
        auto auxPtxts0     = auxPtxts[0];
        for (size_t i = 0; i != 16; ++i) {
            auto diagR     = auxPtxts0[i];
            auto diagRPoly = getFixedPointVecFromDCRTPoly(diagR, scalingFactor, 32);
            auto diagRInC  = RPolynomial(diagRPoly).toCSlots();
            results1[i]    = diagRInC.getSlots();
        }

        [[maybe_unused]] auto cVecMult = [](std::vector<BigComplex> a, std::vector<BigComplex> b) {
            std::vector<BigComplex> result(16);
            for (size_t i = 0; i != 16; ++i) {
                result[i] = a[i] * b[i];
            }
            return result;
        };
        [[maybe_unused]] auto cVecAdd = [](std::vector<BigComplex> a, std::vector<BigComplex> b) {
            std::vector<BigComplex> result(16);
            for (size_t i = 0; i != 16; ++i) {
                result[i] = a[i] + b[i];
            }
            return result;
        };
        [[maybe_unused]] auto cVecRotate = [](std::vector<BigComplex> a, size_t r) {
            std::vector<BigComplex> result(16);
            for (size_t i = 0; i != 16; ++i) {
                result[i] = a[(i + r) % 16];
            }
            return result;
        };

        //auto start = cVecMult(results1[0], cSlots.getSlots());
        //for (size_t i = 1; i != 16; ++i) {
        //    auto rotated = cVecRotate(cSlots.getSlots(), i);
        //    auto mult    = cVecMult(results1[i], rotated);
        //    start        = cVecAdd(start, mult);
        //}
        //for (size_t i = 0; i != 16; ++i) {
        //    encodeZCoeffsToSlots[i] = start[i];
        //}

        //for (size_t i = 0; i != zN; ++i) {
        //    for (size_t j = 0; j != zN / 2; ++j) {
        //        encodeZCoeffsToSlots[i] = encodeZCoeffsToSlots[i] + zUInverse[i][j] * cSlots[j];
        //    }
        //}
        //for (size_t i = 0; i != encodeZCoeffsToSlots.size(); ++i) {
        //    encodeZCoeffsToSlots[i] = encodeZCoeffsToSlots[i] + encodeZCoeffsToSlots[i].conj();
        //}
        //for (size_t i = 0; i != encodeZCoeffsToSlots.size(); ++i) {
        //    std::cout << msg << "  encodeZCoeffsToSlots [" << i << "]: " << encodeZCoeffsToSlots[i].toHexString(16)
        //              << std::endl;
        //}
    }
#endif
    return 0;
}

std::vector<BigComplex> another_interpolate(int order = 1) {
    size_t p = 16;

    auto f = [&](auto x) {
        return double(x) / p;
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

Ciphertext<DCRTPoly> MSBBootstrap(CiphertextT ct) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ct->GetCryptoParameters());

    auto paramsQ   = cryptoParams->GetElementParams()->GetParams();
    uint32_t sizeQ = paramsQ.size();
    std::vector<NativeInteger> moduli(sizeQ);
    std::vector<NativeInteger> roots(sizeQ);
    for (uint32_t i = 0; i < sizeQ; ++i) {
        moduli[i] = paramsQ[i]->GetModulus();
        roots[i]  = paramsQ[i]->GetRootOfUnity();
    }

    auto cc = ct->GetCryptoContext();
    auto M  = cc->GetCyclotomicOrder();
    auto N  = cc->GetRingDimension();

    auto elementParamsRaisedPtr = std::make_shared<ILDCRTParams<DCRTPoly::Integer>>(M, moduli, roots);

    //------------------------------------------------------------------------------
    // RAISING THE MODULUS
    //------------------------------------------------------------------------------

    auto raised = ct->Clone();
    auto algo   = cc->GetScheme();

    uint32_t L0 = cryptoParams->GetElementParams()->GetParams().size();

    if (cryptoParams->GetSecretKeyDist() == SPARSE_ENCAPSULATED) {
        auto evalKeyMap = cc->GetEvalAutomorphismKeyMap(raised->GetKeyTag());

        // transform from a denser secret to a sparser one
        raised = FHECKKSRNS::KeySwitchSparse(raised, evalKeyMap.at(2 * N - 4));

        // Only level 0 ciphertext used here. Other towers ignored to make CKKS bootstrapping faster.
        auto& ctxtDCRTs = raised->GetElements();

        for (auto& dcrt : ctxtDCRTs) {
            dcrt.SetFormat(COEFFICIENT);
            DCRTPoly tmp(dcrt.GetElementAtIndex(0), elementParamsRaisedPtr);
            tmp.SetFormat(EVALUATION);
            dcrt = std::move(tmp);
        }
        raised->SetLevel(L0 - ctxtDCRTs[0].GetNumOfElements());

        // go back to a denser secret
        algo->KeySwitchInPlace(raised, evalKeyMap.at(2 * N - 2));
    }
    else {
        // Only level 0 ciphertext used here. Other towers ignored to make CKKS bootstrapping faster.
        auto& ctxtDCRTs = raised->GetElements();

        for (auto& dcrt : ctxtDCRTs) {
            dcrt.SetFormat(COEFFICIENT);
            DCRTPoly tmp(dcrt.GetElementAtIndex(0), elementParamsRaisedPtr);
            tmp.SetFormat(EVALUATION);
            dcrt = std::move(tmp);
        }
        raised->SetLevel(L0 - ctxtDCRTs[0].GetNumOfElements());
    }

    raised->SetScalingFactorBFP(ct->GetScalingFactorBFP());
    //__heir_debug2(raised, "ModRaise");

    //------------------------------------------------------------------------------
    // SPARSELY PACKED CASE
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Running PartialSum
    //------------------------------------------------------------------------------

    // Here 32 is rPoly.size()
    const uint32_t limit = N / 32;
    for (uint32_t j = 1; j < limit; j <<= 1) {
        cc->EvalAddInPlace(raised, cc->EvalRotate(raised, j * 16));
    }
    // Now the message is multplied by N/32
    //__heir_debug2(raised, "PSum");

    // Normalize to [-1, 1] from [-16, 16]
    // Multiply by 2/N * Delta, so the result is m / 16 * Delta^2
    // NOTE: two here should be changed.
    auto raisedSF                    = raised->GetScalingFactorBFP();
    auto NBigFP                      = BigFixedPoint::positive(N);
    auto two                         = BigFixedPoint::two();
    BigFixedPoint normalizeFactorBFP = raisedSF * two / NBigFP;
    BigInteger normalizeFactor       = (normalizeFactorBFP.round().getValue()) >> normalizeFactorBFP.getLog2Scale();
    raised                           = gEvalMultScalar(raised, normalizeFactor);
    raised->SetScalingFactorBFP(raisedSF * raisedSF);
    gModReduceInPlace(raised);
    //__heir_debug2(raised, "Normalize");
    // Note that there are other ways...some work first multiply by 1 / N
    // Then PartialSum
    // Then use CoeffsToSlots matrix to do the /16

    //------------------------------------------------------------------------------
    // Running CoeffsToSlots
    //------------------------------------------------------------------------------

    auto raisedRC2S = RCoeffsToSlots(cc, raised);
    gModReduceInPlace(raisedRC2S[0]);
    gModReduceInPlace(raisedRC2S[1]);
    //__heir_debug2(raisedRC2S[0], "MSBC2S");
    //__heir_debug2(raisedRC2S[1], "MSBC2S");

    //------------------------------------------------------------------------------
    // Running Approximate Mod Reduction
    //------------------------------------------------------------------------------

    auto& coeff_exp = coeff_exp_16_big_complex_58;
    auto res        = EvalChebyshevSeriesPS(raisedRC2S[0], coeff_exp);
    auto resI       = EvalChebyshevSeriesPS(raisedRC2S[1], coeff_exp);

    // Double angle-iterations to get exp(2*Pi*i*x)
    res = gEvalMult(res, res);
    gModReduceInPlace(res);
    res = gEvalMult(res, res);
    gModReduceInPlace(res);

    resI = gEvalMult(resI, resI);
    gModReduceInPlace(resI);
    resI = gEvalMult(resI, resI);
    gModReduceInPlace(resI);

    //__heir_debug2(res, "Cheby1");

    //------------------------------------------------------------------------------
    // Running LUT
    //------------------------------------------------------------------------------

    auto lutCoeffs = another_interpolate(2);
    auto powers    = EvalPowers(res, lutCoeffs);
    auto lut       = EvalPolyWithPrecomp(powers, lutCoeffs);
    //__heir_debug2(lut, "LUT");

    auto powersI = EvalPowers(resI, lutCoeffs);
    auto lutI    = EvalPolyWithPrecomp(powersI, lutCoeffs);
    //__heir_debug2(lutI, "LUT");

    auto rCoeffLut = SlotsToRCoeffs(cc, lut, lutI);
    gModReduceInPlace(rCoeffLut);
    //__heir_debug2(rCoeffLut, "LUTR");
    return rCoeffLut;
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
    parameters.SetRingDim(1 << 12);
    //parameters.SetNumLargeDigits(6);

    /*  A3) Scaling parameters.
    * By default, we set the modulus sizes and rescaling technique to the following values
    * to obtain a good precision and performance tradeoff. We recommend keeping the parameters
    * below unless you are an FHE expert.
    */
#if NATIVEINT == 128
    ScalingTechnique rescaleTech = FIXEDMANUAL;
    uint32_t dcrtBits            = 78;
    uint32_t firstMod            = 89;
#else
    ScalingTechnique rescaleTech = FLEXIBLEMANUAL;
    uint32_t dcrtBits            = 59;
    uint32_t firstMod            = 60;
#endif

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
    std::vector<uint32_t> levelBudget = {1, 1};

    // Note that the actual number of levels avalailable after bootstrapping before next bootstrapping
    // will be levelsAvailableAfterBootstrap - 1 because an additional level
    // is used for scaling the ciphertext before next bootstrapping (in 64-bit CKKS bootstrapping)
    //uint32_t levelsAvailableAfterBootstrap = 10;
    //uint32_t depth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(25);

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

    //cc->EvalBootstrapSetup(levelBudget);

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    std::vector<int> rotateIndices = {};
    for (int i = 1; i <= 31; ++i) {
        rotateIndices.push_back(i);
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

    cc_global    = cc;
    pk_global    = keyPair.publicKey;
    sk_global    = keyPair.secretKey;
    slots_global = 32 / 2;

    auto sf = BigFixedPoint(BigInteger(1) << dcrtBits, 0, false).scaleTo(128);

    auto zero      = EncryptZero(keyPair.publicKey);
    auto elemParam = zero->GetElements()[0].GetParams();
    //
    //__heir_debug2(zero, "Input");

    RPolynomial value1 = ZPolynomial::encode(-1).toRPolynomial();
    Plaintext ptxt1    = ZEncodingImpl::encodeR(value1, elemParam, sf);

    RPolynomial value2 = ZPolynomial::encode(-3).toRPolynomial();
    Plaintext ptxt2    = ZEncodingImpl::encodeR(value2, elemParam, sf);

    /// TEST ENCODE
    auto encoded  = Encrypt(ptxt1, keyPair.publicKey);
    auto encoded2 = Encrypt(ptxt2, keyPair.publicKey);

    __heir_debug2(encoded, "Encode");

    /// TEST ADD
    if (0) {
        auto ctAdd = gEvalAdd(encoded, ptxt2);

        __heir_debug2(ctAdd, "Add");
    }

    /// TEST CT-PT-MULT
    RPolynomial t   = ZPolynomial::getT().toRPolynomial();
    Plaintext tPtxt = ZEncodingImpl::encodeR(t, elemParam, sf);

    if (0) {
        auto ctMul = gEvalMult(encoded, tPtxt);
        __heir_debug2(ctMul, "Mult");
        gModReduceInPlace(ctMul);
        __heir_debug2(ctMul, "Mult");
    }

    /// TEST CT-CT-MULT
    Ciphertext<DCRTPoly> ct;
    if (1) {
        auto ctMul = zEvalMultFull(encoded, encoded2, tPtxt);
        gModReduceInPlace(ctMul, 2);
        __heir_debug2(ctMul, "CMult");
        ct = ctMul;
    }

    /// TEST Rotate
    if (0) {
        auto ctRot = cc->EvalRotate(encoded, 1);

        __heir_debug2(ctRot, "Rotate");
    }

    /// TEST ZCoeffToSlots and SlotsToZCoeffs
    std::vector<Ciphertext<DCRTPoly>> zC2S;
    if (1) {
        zC2S = ZCoeffsToSlots(cc, ct);
        gModReduceInPlace(zC2S[0]);
        gModReduceInPlace(zC2S[1]);

        __heir_debug2(zC2S[0], "Upper");
        __heir_debug2(zC2S[1], "Down");

        // Now SF is sf
        if (0) {
            auto z2S2z = SlotsToZCoeffs(cc, zC2S[0], zC2S[1]);

            __heir_debug2(z2S2z, "Z2S2Z");
        }
    }

    // TEST SlotsToRCoeffs
    Ciphertext<DCRTPoly> s2rc;
    if (1) {
        auto z2S2r = SlotsToRCoeffs(cc, zC2S[0], zC2S[1]);
        gModReduceInPlace(z2S2r);
        s2rc = z2S2r;
        __heir_debug2(z2S2r, "S2RC");

        //auto z2S2r_1 = EvalMultScalar(z2S2r, 2);
        ////z2S2r_1->SetScalingFactor(sf * sf * 2);
        //__heir_debug2(z2S2r_1, "S2RC");

        //auto z2S2r_2 = EvalMultScalar(z2S2r, 4);
        ////z2S2r_2->SetScalingFactor(sf * sf * 4);
        //__heir_debug2(z2S2r_2, "S2RC");

        //auto z2S2r_3 = EvalMultScalar(z2S2r, 8);
        ////z2S2r_3->SetScalingFactor(sf * sf * 8);
        //__heir_debug2(z2S2r_3, "S2RC");
    }

    // TEST Scale to MSB
    Ciphertext<DCRTPoly> ctMSB;
    if (1) {
        auto q     = s2rc->GetElements()[0].GetModulus();
        auto sfNow = s2rc->GetScalingFactorBFP();
        auto qBFP  = BigFixedPoint(q, 0, false).scaleTo(128);
        auto two   = BigFixedPoint::two();
        // q / (2 * Delta)
        auto div       = (qBFP / sfNow / two).round();
        auto divScalar = div.getValue() >> div.getLog2Scale();
        auto ct2       = gEvalMultScalar(s2rc, divScalar);
        ct2->SetScalingFactorBFP(qBFP / two);
        // Reduce all the way to the bottom
        gModReduceInPlace(ct2, ct2->GetElements()[0].GetNumOfElements() - 1);
        //auto q0     = ct2->GetElements()[0].GetModulus();
        //auto q0BFP  = BigFixedPoint(q0, 0, false).scaleTo(128);
        //auto sfNow2 = q0BFP / two;
        //ct2->SetScalingFactorBFP(sfNow2);
        __heir_debug2(ct2, "MSB");
        ctMSB = ct2;
    }

    auto lowBitBTS = [](Ciphertext<DCRTPoly> input, unsigned bits) -> std::array<Ciphertext<DCRTPoly>, 2> {
        auto lowScalar = BigInteger(1) << bits;
        auto ctLow     = gEvalMultScalar(input, lowScalar);
        // Just a different interpretation...
        ctLow->SetScalingFactorBFP(input->GetScalingFactorBFP());
        auto ctLowBTS = MSBBootstrap(ctLow);

        // Adjust to original MSB representation
        auto q            = ctLowBTS->GetElements()[0].GetModulus();
        auto sfNow        = ctLowBTS->GetScalingFactorBFP();
        auto qBFP         = BigFixedPoint(q, 0, false).scaleTo(128);
        auto two          = BigFixedPoint::two();
        auto lowScalarBFP = BigFixedPoint(BigInteger(1) << bits, 0, false).scaleTo(128);
        // q / (2 * Delta * lowScalar)
        auto div       = (qBFP / sfNow / two / lowScalarBFP).round();
        auto divScalar = div.getValue() >> div.getLog2Scale();
        auto ct2       = gEvalMultScalar(ctLowBTS, divScalar);
        ct2->SetScalingFactorBFP(qBFP / two);
        // Reduce all the way to the bottom
        gModReduceInPlace(ct2, ct2->GetElements()[0].GetNumOfElements() - 1);
        return {ct2, ctLowBTS};
    };

    if (1) {
        std::vector<Ciphertext<DCRTPoly>> lowBTSs;
        std::vector<Ciphertext<DCRTPoly>> lowBTSHighs;
        for (size_t bits = 4; bits <= 32; bits += 4) {
            auto [lowBTS, lowBTSHigh] = lowBitBTS(ctMSB, 32 - bits);
            lowBTSs.push_back(lowBTS);
            lowBTSHighs.push_back(lowBTSHigh);
            gEvalSubInPlace(ctMSB, lowBTS);
            //__heir_debug2(ctMSB, "MSB" + std::to_string(bits));
        }
        auto ctNew = lowBTSs[0];
        for (size_t i = 1; i != lowBTSs.size(); ++i) {
            gEvalAddInPlace(ctNew, lowBTSs[i]);
        }
        __heir_debug2(ctNew, "Reconstructed");
    }

    // MSBBootstrap(cc, z2S2r);

    // SlotsToR-Coeffs

    // auto encoded2 = encodeREInCt(ciph, 2);
    // __heir_debug2(encoded2, "Encode");

    // // should use another empty ct
    // auto ZPtm = encodeZPtmInCt(ciph);

    // auto multResultHalf = cc->EvalMult(encoded, encoded2);
    // auto multResult     = cc->EvalMult(multResultHalf, ZPtm);

    // __heir_debug2(multResult, "Mult");
}
