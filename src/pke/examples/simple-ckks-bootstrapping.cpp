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

    auto sf      = ct->GetScalingFactor();
    auto log2sf  = std::log2(sf);
    auto sfBigFP = BigFixedPoint(BigInteger(1) << log2sf, 0, false).scaleTo(128);
    //std::cout << msg << "  Scaling factor: " << std::log2(sf) << std::endl;

    // valueSize = zN
    auto valueSize = slots_global * 2;

    RPolynomial values = getFixedPointVecFromDCRTPoly(b, sfBigFP, valueSize);

    // Check the slot encoding?
    if (msg == "S2RC" || msg == "MSB") {
        for (size_t i = 0; i != values.getCoefficients().size(); ++i) {
            std::cout << msg << "  values [" << i << "]: " << values[i].toHexString(16) << std::endl;
        }
    }
    auto cSlots = values.toCSlots();
    if (msg == "Upper" || msg == "Down" || msg == "Rotate") {
        for (size_t i = 0; i != cSlots.getSlots().size(); ++i) {
            std::cout << msg << "  complexValues [" << i << "]: " << cSlots[i].toHexString(32) << std::endl;
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
            ZPolynomial zPoly(zCoeffs);
            std::cout << msg << "  zValues error log2Norm: " << ZPolynomial::extractError(zPoly).getLog2Norm()
                      << std::endl;
        }
    }
    ZPolynomial zValues = values.toZPolynomial();
    if (msg == "Encode" || msg == "Input" || msg == "Add" || msg == "Mult" || msg == "CMult" || msg == "Z2S2Z") {
        for (size_t i = 0; i != values.getCoefficients().size(); ++i) {
            std::cout << msg << "  zValues [" << i << "]: " << zValues[i].toHexString(32) << std::endl;
        }
        auto decoded = ZPolynomial::decode(zValues);
        std::cout << msg << "  Decoded: " << decoded << std::endl;
        std::cout << msg << "  zValues error log2Norm: " << ZPolynomial::extractError(zValues).getLog2Norm()
                  << std::endl;
    }
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
        auto scalingFactor = BigFixedPoint(BigInteger(1) << static_cast<uint32_t>(std::log2(sf)), 0, false);
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

    //auto rounded   = roundInRe(zValues);
    //auto decodeInR = decodeFromRE(rounded);
    //std::cout << "Z: " << decodeInR << std::endl;

    //if (msg != "Input" && msg != "ModRaise") {
    //    for (size_t i = 0; i != 1; ++i) {
    //        std::cout << msg << "  complexValues [" << i << "]: " << complexValues[i] << std::endl;
    //    }
    //}
    return 0;
}

//void MSBBootstrap(CryptoContextT cc, CiphertextT ct) {
//    auto q      = ct->GetElements()[0].GetModulus();
//    auto sf     = ct->GetScalingFactor();
//    auto log2sf = std::log2(sf);
//    std::cout << "q: " << q.GetMSB() << " log2sf: " << log2sf << std::endl;
//    auto multBy = q >> (log2sf - 1);
//    auto ct2    = ct->Clone();
//    auto& cv    = ct2->GetElements();
//    cv[0] *= BigInteger(64);
//    cv[1] *= BigInteger(64);
//    __heir_debug2(ct2, "MSB");
//}

void SimpleBootstrapExample();

int main(int argc, char* argv[]) {
    //test4_encodeInRE();
    //test_UInverse();
    //test_zu();
    //test4_encodeInRE();
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
    SecretKeyDist secretKeyDist = lbcrypto::UNIFORM_TERNARY;
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
    parameters.SetNumLargeDigits(5);

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
    ScalingTechnique rescaleTech = FIXEDMANUAL;
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
    parameters.SetMultiplicativeDepth(4);

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
    cc->EvalRotateKeyGen(keyPair.secretKey, rotateIndices);
    cc->EvalAutomorphismKeyGen(keyPair.secretKey, {2 * ringDim - 1});
    //cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);

    cc_global    = cc;
    pk_global    = keyPair.publicKey;
    sk_global    = keyPair.secretKey;
    slots_global = 32 / 2;

    auto sf      = std::pow(2.0, dcrtBits);
    auto sfBigFP = BigFixedPoint(BigInteger(1) << dcrtBits, 0, false).scaleTo(128);

    auto zero      = EncryptZero(sf, keyPair.publicKey);
    auto elemParam = zero->GetElements()[0].GetParams();

    //__heir_debug2(zero, "Input");

    RPolynomial value1 = ZPolynomial::encode(255).toRPolynomial();
    DCRTPoly ptxt1     = getDCRTPolyFromFixedPointVec(value1.getCoefficients(), elemParam, sfBigFP);

    RPolynomial value2 = ZPolynomial::encode(2).toRPolynomial();
    DCRTPoly ptxt2     = getDCRTPolyFromFixedPointVec(value2.getCoefficients(), elemParam, sfBigFP);

    /// TEST ENCODE
    auto encoded  = EncryptDCRTPoly(ptxt1, sf, keyPair.publicKey);
    auto encoded2 = EncryptDCRTPoly(ptxt2, sf, keyPair.publicKey);

    __heir_debug2(encoded, "Encode");

    /// TEST ADD
    if (0) {
        auto ctAdd = EvalAddDCRTPoly(encoded, ptxt2);

        __heir_debug2(ctAdd, "Add");
    }

    /// TEST CT-PT-MULT
    RPolynomial t = ZPolynomial::getT().toRPolynomial();

    DCRTPoly tPtxt = getDCRTPolyFromFixedPointVec(t.getCoefficients(), elemParam, sfBigFP);

    if (0) {
        auto ctMul = EvalMultDCRTPoly(encoded, tPtxt);
        ctMul->SetScalingFactor(sf * sf);
        __heir_debug2(ctMul, "Mult");
        ModReduceCustomInPlace(ctMul);
        ctMul->SetScalingFactor(sf);
        __heir_debug2(ctMul, "Mult");
    }

    /// TEST CT-CT-MULT
    Ciphertext<DCRTPoly> ct;
    if (0) {
        auto ctMulRaw = cc->EvalMult(encoded, encoded2);
        ModReduceCustomInPlace(ctMulRaw);
        ctMulRaw->SetScalingFactor(sf);
        __heir_debug2(ctMulRaw, "CMult");

        auto ctMulRawT = EvalMultDCRTPoly(ctMulRaw, tPtxt);
        ctMulRawT->SetScalingFactor(sf * sf);
        __heir_debug2(ctMulRawT, "CMult");

        ModReduceCustomInPlace(ctMulRawT);
        ctMulRawT->SetScalingFactor(sf);
        ct = ctMulRawT;
    }

    /// TEST Rotate
    if (0) {
        auto ctRot = cc->EvalRotate(encoded, 1);

        __heir_debug2(ctRot, "Rotate");
    }

    /// TEST ZCoeffToSlots and SlotsToZCoeffs
    {
        auto zC2S = ZCoeffsToSlots(cc, encoded);
        zC2S[0]->SetScalingFactor(sf * sf);
        zC2S[1]->SetScalingFactor(sf * sf);

        __heir_debug2(zC2S[0], "Upper");
        __heir_debug2(zC2S[1], "Down");

        //ModReduceCustomInPlace(zC2S[0]);
        //ModReduceCustomInPlace(zC2S[1]);
        //// Now SF is sf

        //auto z2S2z = SlotsToZCoeffs(cc, zC2S[0], zC2S[1]);
        //z2S2z->SetScalingFactor(sf * sf * 32 * 32);

        //__heir_debug2(z2S2z, "Z2S2Z");
    }

    // auto zC2S = ZCoeffToSlots(cc, encoded, zero);

    // __heir_debug2(zC2S[0], "Upper");
    // __heir_debug2(zC2S[1], "Down");

    // auto z2S2r = SlotsToRCoeffs(cc, zC2S[0], zC2S[1], zero);

    // __heir_debug2(z2S2r, "S2RC");

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
