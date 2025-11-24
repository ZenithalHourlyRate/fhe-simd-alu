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
#include "re-utils.h"

using namespace lbcrypto;

// related to __heir_debug2
CryptoContextT cc_global;
PrivateKeyT sk_global;
size_t slots_global;

double __heir_debug2(CiphertextT ct, std::string msg) {
    auto b = DecryptCore(ct->GetElements(), sk_global);
    b.SetFormat(Format::COEFFICIENT);
    auto bigB       = b.CRTInterpolate();
    auto& bigCoeffs = bigB.GetValues();

    auto sf = ct->GetScalingFactor();
    //auto level = ct->GetLevel();

    //std::cout << msg << "  Scaling factor: " << std::log2(sf) << std::endl;

    const auto& q{b.GetParams()->GetModulus()};
    const auto& half{q >> 1};

    std::vector<double> values;
    for (size_t i = 0; i < b.GetRingDimension(); i++) {
        auto x   = bigCoeffs[i];
        bool neg = false;
        if (x > half) {
            x   = q - x;
            neg = true;
        }
        double value = x.ConvertToDouble() / sf;
        if (neg)
            value = -value;
        // generic for different packed cases
        if ((msg == "Input" || msg == "Encode" || msg == "Mult") &&
            (i % (b.GetRingDimension() / (slots_global * 2)) == 0)) {
            values.push_back(value);
        }
    }

    // Check the slot encoding?
    //for (size_t i = 0; i != values.size(); ++i) {
    //    std::cout << msg << "  values [" << i << "]: " << values[i] << std::endl;
    //}
    auto complexValues = multiplyByUComplex(values);
    auto zValues       = multiplyByZUInverse(complexValues);
    for (size_t i = 0; i != values.size(); ++i) {
        std::cout << msg << "  zValues [" << i << "]: " << zValues[i] << std::endl;
    }
    auto rounded   = roundInRe(zValues);
    auto decodeInR = decodeFromRE(rounded);
    std::cout << "Z: " << decodeInR << std::endl;

    //if (msg != "Input" && msg != "ModRaise") {
    //    for (size_t i = 0; i != 1; ++i) {
    //        std::cout << msg << "  complexValues [" << i << "]: " << complexValues[i] << std::endl;
    //    }
    //}
    return 0;
}

DCRTPoly getPolyFromVec(std::vector<double> input, const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                        double scalingFactor) {
    DCRTPoly newPoly(elementParams, Format::COEFFICIENT, true);
    auto bigPoly = newPoly.CRTInterpolate();
    auto ringDim = elementParams->GetRingDimension();
    for (size_t i = 0; i != input.size(); ++i) {
        double val = input[i] * scalingFactor;
        bool neg   = false;
        if (val < 0) {
            val = -val;
            neg = true;
        }
        int64_t intVal = static_cast<int64_t>(std::round(val));
        if (neg) {
            bigPoly[i * (ringDim / (input.size()))] = elementParams->GetModulus() - intVal;
        }
        else {
            bigPoly[i * (ringDim / (input.size()))] = intVal;
        }
    }
    DCRTPoly finalPoly(bigPoly, elementParams);
    finalPoly.SetFormat(Format::EVALUATION);
    return finalPoly;
}

CiphertextT encodeREVecInCt(CiphertextT ct, std::vector<double> input) {
    auto newCt         = ct->Clone();
    auto& cv0          = ct->GetElements()[0];
    double sf          = ct->GetScalingFactor();
    auto elementParams = cv0.GetParams();

    auto encodedC = multiplyByZUComplex(input);
    auto encodedR = multiplyByUInverseComplex(encodedC);

    auto finalPoly = getPolyFromVec(encodedR, elementParams, sf);
    finalPoly.SetFormat(Format::EVALUATION);
    auto& cv = newCt->GetElements();
    cv[0] += finalPoly;
    return newCt;
}

CiphertextT encodeREInCt(CiphertextT ct, int32_t input) {
    auto encodedZ = encodeInRE(input);
    return encodeREVecInCt(ct, encodedZ);
}

CiphertextT encodeZPtmInCt(CiphertextT ct) {
    // X-2 in calZ
    std::vector<double> encodedZPtm = {-2, 1, 0, 0, 0, 0, 0, 0};
    assert(encodedZPtm.size() == zN && "Encoded ZPtm size does not match zN");
    return encodeREVecInCt(ct, encodedZPtm);
}

CiphertextT multiplyByZptm(CryptoContextT cc, CiphertextT ct) {
    // X-2 in calZ
    std::vector<double> encodedZPtm = {-2, 1, 0, 0, 0, 0, 0, 0};
    auto encodedC                   = multiplyByZUComplex(encodedZPtm);
    auto encodedR                   = multiplyByUInverseComplex(encodedC);

    auto newCt         = ct->Clone();
    auto& cv0          = ct->GetElements()[0];
    double sf          = ct->GetScalingFactor();
    auto elementParams = cv0.GetParams();
    auto finalPoly     = getPolyFromVec(encodedR, elementParams, sf);
    finalPoly.SetFormat(Format::EVALUATION);
    auto& cv = newCt->GetElements();
    cv[0] *= finalPoly;
    cv[1] *= finalPoly;
    cc->ModReduceInPlace(newCt);
    return newCt;
}

void SimpleBootstrapExample();

int main(int argc, char* argv[]) {
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
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
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
    parameters.SetMultiplicativeDepth(5);

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);

    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    //cryptoContext->Enable(ADVANCEDSHE);
    //cryptoContext->Enable(FHE);

    uint32_t ringDim = cryptoContext->GetRingDimension();
    // This is the maximum number of slots that can be used for full packing.
    //uint32_t numSlots = ringDim / 2;
    std::cout << "CKKS scheme ring dimension: " << ringDim << "\n\n";

    //cryptoContext->EvalBootstrapSetup(levelBudget);

    auto keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeyGen(keyPair.secretKey);
    //cryptoContext->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);

    cc_global    = cryptoContext;
    sk_global    = keyPair.secretKey;
    slots_global = zN / 2;

    std::vector<double> x = {0, 0, 0, 0, 0, 0, 0, 0};
    size_t encodedLength  = x.size();

    // We start with a depleted ciphertext that has used up all of its levels.
    Plaintext ptxt = cryptoContext->MakeCKKSPackedPlaintext(x, 1);

    ptxt->SetLength(encodedLength);
    std::cout << "Input: " << ptxt << "\n";

    Ciphertext<DCRTPoly> ciph = cryptoContext->Encrypt(keyPair.publicKey, ptxt);

    __heir_debug2(ciph, "Input");

    auto encoded = encodeREInCt(ciph, 255);
    __heir_debug2(encoded, "Encode");
    auto encoded2 = encodeREInCt(ciph, 2);
    __heir_debug2(encoded2, "Encode");

    // should use another empty ct
    auto ZPtm = encodeZPtmInCt(ciph);

    auto multResultHalf = cryptoContext->EvalMult(encoded, encoded2);
    auto multResult     = cryptoContext->EvalMult(multResultHalf, ZPtm);

    __heir_debug2(multResult, "Mult");
}
