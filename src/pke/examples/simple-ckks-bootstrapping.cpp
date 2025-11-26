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

double __heir_debug2(CiphertextT ct, std::string msg) {
    auto b = DecryptCore(ct->GetElements(), sk_global);

    auto sf      = ct->GetScalingFactor();
    auto log2sf  = std::log2(sf);
    auto sfBigFP = BigFixedPoint(BigInteger(1) << log2sf, 0, false).scaleTo(z_upper_roots_scale);
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
    if (msg == "Upper" || msg == "Down" || msg == "Aux" || msg == "AuxI" || msg == "AuxJ") {
        for (size_t i = 0; i != cSlots.getSlots().size(); ++i) {
            std::cout << msg << "  complexValues [" << i << "]: " << cSlots[i].toHexString(16) << std::endl;
        }
    }
    ZPolynomial zValues = values.toZPolynomial();
    if (msg == "Encode" || msg == "Input" || msg == "Add" || msg == "Mult") {
        for (size_t i = 0; i != values.getCoefficients().size(); ++i) {
            std::cout << msg << "  zValues [" << i << "]: " << zValues[i].toHexString(16) << std::endl;
        }
        std::cout << msg << "  zValues error log2Norm: " << ZPolynomial::extractError(zValues).getLog2Norm()
                  << std::endl;
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

//std::array<std::vector<Ciphertext<DCRTPoly>>, 2> getAuxUInverseCt(CryptoContextT cc, CiphertextT zero) {
//    auto UT = getUT(zN * 2);
//    return getAuxInverseCt(cc, zero, UT);
//}
//
//std::array<std::vector<Ciphertext<DCRTPoly>>, 2> getAuxZUInverseCt(CryptoContextT cc, CiphertextT zero) {
//    auto zUInv = getZUInverse();
//    return getAuxInverseCt(cc, zero, zUInv);
//}

//std::array<std::vector<Ciphertext<DCRTPoly>>, 2> getAuxCt(CryptoContextT cc, CiphertextT zero, CMatrix T) {
//    // T is of shape (zN / 2) * zN
//    auto halfSize = T.size();
//
//    std::vector<Ciphertext<DCRTPoly>> results1;
//    std::vector<Ciphertext<DCRTPoly>> results2;
//
//    // up matrix
//    for (size_t i = 0; i != halfSize; ++i) {
//        auto newCt    = zero->Clone();
//        auto diagonal = std::vector<std::complex<double>>(halfSize, 0);
//        for (size_t j = 0; j != halfSize; ++j) {
//            diagonal[j] = T[j][(i + j) % halfSize];
//        }
//        auto diagonalInR = multiplyByUInverseComplex(diagonal);
//        auto finalPoly   = getPolyFromVec(diagonalInR, zero->GetElements()[0].GetParams(), zero->GetScalingFactor());
//        finalPoly.SetFormat(Format::EVALUATION);
//        auto& cv = newCt->GetElements();
//        cv[0] += finalPoly;
//        results1.push_back(newCt);
//    }
//    // down matrix
//    for (size_t i = 0; i != halfSize; ++i) {
//        auto newCt    = zero->Clone();
//        auto diagonal = std::vector<std::complex<double>>(halfSize, 0);
//        for (size_t j = 0; j != halfSize; ++j) {
//            diagonal[j] = T[j][(i + j) % halfSize + halfSize];
//        }
//        auto diagonalInR = multiplyByUInverseComplex(diagonal);
//        auto finalPoly   = getPolyFromVec(diagonalInR, zero->GetElements()[0].GetParams(), zero->GetScalingFactor());
//        finalPoly.SetFormat(Format::EVALUATION);
//        auto& cv = newCt->GetElements();
//        cv[0] += finalPoly;
//        results2.push_back(newCt);
//    }
//    return {results1, results2};
//}
//
// std::array<std::vector<Ciphertext<DCRTPoly>>, 2> getAuxUCt(CryptoContextT cc, CiphertextT zero) {
//     auto U = getU(zN * 2);
//     return getAuxCt(cc, zero, U);
// }
//
// std::array<std::vector<Ciphertext<DCRTPoly>>, 2> getAuxZUCt(CryptoContextT cc, CiphertextT zero) {
//     auto zU = getZU();
//     return getAuxCt(cc, zero, zU);
// }

std::vector<Ciphertext<DCRTPoly>> CoeffsToSlots(CryptoContextT cc, CiphertextT ct, CiphertextT zero,
                                                const std::array<std::vector<Ciphertext<DCRTPoly>>, 2>& auxCts) {
    // Halevi-Shoup
    // Z-CoeffToSlots
    auto startCt     = ct->Clone();
    auto resultUpper = zero->Clone();
    auto resultDown  = zero->Clone();
    for (size_t i = 0; i != auxCts[0].size(); ++i) {
        auto diagonalUpper = cc->EvalMult(startCt, auxCts[0][i]);
        resultUpper        = cc->EvalAdd(resultUpper, diagonalUpper);
        auto diagonalDown  = cc->EvalMult(startCt, auxCts[1][i]);
        resultDown         = cc->EvalAdd(resultDown, diagonalDown);
        // rotate one more
        startCt = cc->EvalRotate(startCt, 1);
    }
    cc->EvalAddInPlace(resultUpper, Conjugate(resultUpper, cc->GetEvalAutomorphismKeyMap(resultUpper->GetKeyTag())));
    cc->EvalAddInPlace(resultDown, Conjugate(resultDown, cc->GetEvalAutomorphismKeyMap(resultDown->GetKeyTag())));

    return {resultUpper, resultDown};
}

//std::vector<Ciphertext<DCRTPoly>> ZCoeffToSlots(CryptoContextT cc, CiphertextT ct, CiphertextT zero) {
//    return CoeffsToSlots(cc, ct, zero, getAuxZUInverseCt(cc, zero));
//}
//
//std::vector<Ciphertext<DCRTPoly>> RCoeffToSlots(CryptoContextT cc, CiphertextT ct, CiphertextT zero) {
//    return CoeffsToSlots(cc, ct, zero, getAuxUInverseCt(cc, zero));
//}

Ciphertext<DCRTPoly> SlotsToCoeffs(CryptoContextT cc, CiphertextT ctLeft, CiphertextT ctRight, CiphertextT zero,
                                   const std::array<std::vector<Ciphertext<DCRTPoly>>, 2>& auxCts) {
    // Halevi-Shoup
    // Z-CoeffToSlots
    auto startCtLeft  = ctLeft->Clone();
    auto startCtRight = ctRight->Clone();
    auto result       = zero->Clone();
    for (size_t i = 0; i != auxCts[0].size(); ++i) {
        auto diagonalLeft  = cc->EvalMult(startCtLeft, auxCts[0][i]);
        result             = cc->EvalAdd(result, diagonalLeft);
        auto diagonalRight = cc->EvalMult(startCtRight, auxCts[1][i]);
        result             = cc->EvalAdd(result, diagonalRight);
        // rotate one more
        startCtLeft  = cc->EvalRotate(startCtLeft, 1);
        startCtRight = cc->EvalRotate(startCtRight, 1);
    }
    return result;
}

//Ciphertext<DCRTPoly> SlotsToRCoeffs(CryptoContextT cc, CiphertextT ctLeft, CiphertextT ctRight, CiphertextT zero) {
//    return SlotsToCoeffs(cc, ctLeft, ctRight, zero, getAuxUCt(cc, zero));
//}
//
//Ciphertext<DCRTPoly> SlotsToZCoeffs(CryptoContextT cc, CiphertextT ctLeft, CiphertextT ctRight, CiphertextT zero) {
//    return SlotsToCoeffs(cc, ctLeft, ctRight, zero, getAuxZUCt(cc, zero));
//}

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
    // To reduce encryption noise
    //parameters.SetEncryptionTechnique(EncryptionTechnique::EXTENDED);

    /*  A3) Scaling parameters.
    * By default, we set the modulus sizes and rescaling technique to the following values
    * to obtain a good precision and performance tradeoff. We recommend keeping the parameters
    * below unless you are an FHE expert.
    */
#if NATIVEINT == 128
    ScalingTechnique rescaleTech = FIXEDAUTO;
    uint32_t dcrtBits            = 78;
    uint32_t firstMod            = 89;
#else
    ScalingTechnique rescaleTech = FIXEDAUTO;
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
    parameters.SetMultiplicativeDepth(3);

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
    cc->EvalRotateKeyGen(keyPair.secretKey, {1});
    cc->EvalAutomorphismKeyGen(keyPair.secretKey, {2 * ringDim - 1});
    //cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);

    cc_global    = cc;
    pk_global    = keyPair.publicKey;
    sk_global    = keyPair.secretKey;
    slots_global = 32 / 2;

    auto sf      = std::pow(2.0, firstMod);
    auto sfBigFP = BigFixedPoint(BigInteger(1) << firstMod, 0, false).scaleTo(z_upper_roots_scale);

    auto zero      = EncryptZero(sf, keyPair.publicKey);
    auto elemParam = zero->GetElements()[0].GetParams();

    __heir_debug2(zero, "Input");

    RPolynomial value1 = ZPolynomial::encode(255).toRPolynomial();
    DCRTPoly ptxt1     = getDCRTPolyFromFixedPointVec(value1.getCoefficients(), elemParam, sfBigFP);

    RPolynomial value2 = ZPolynomial::encode(1).toRPolynomial();
    DCRTPoly ptxt2     = getDCRTPolyFromFixedPointVec(value2.getCoefficients(), elemParam, sfBigFP);

    auto encoded = EncryptDCRTPoly(ptxt1, sf, keyPair.publicKey);

    __heir_debug2(encoded, "Encode");

    auto ctAdd = EvalAddDCRTPoly(encoded, ptxt2);

    __heir_debug2(ctAdd, "Add");

    RPolynomial t = ZPolynomial::getT().toRPolynomial();

    DCRTPoly tPtxt = getDCRTPolyFromFixedPointVec(t.getCoefficients(), elemParam, sfBigFP);

    //auto ctMul1 = EvalMultDCRTPoly(ctAdd, ptxt2);
    //cc->ModReduceInPlace(ctMul1);
    auto ctMul2 = EvalMultDCRTPoly(encoded, tPtxt);
    ctMul2->SetScalingFactor(sf * sf);
    __heir_debug2(ctMul2, "Mult");
    ModReduceCustomInPlace(ctMul2);
    ctMul2->SetScalingFactor(sf);
    __heir_debug2(ctMul2, "Mult");

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
