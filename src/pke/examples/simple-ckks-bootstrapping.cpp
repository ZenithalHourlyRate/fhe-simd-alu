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
size_t slots_global;
size_t zN_global;
size_t zSlots_global;

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

    ZEncodingParams params(ZMode, zN_global, zSlots_global);
    auto zEncode = std::make_shared<ZEncodingImpl>(b.GetParams(), b, sfBigFP, params);

    RPolynomial values = ZEncodingImpl::decodeR(zEncode);

    // This is for sparse LT
    ZEncodingParams paramsTwice(CMode, zN_global * zSlots_global * 2);
    auto zEncodeTwice       = std::make_shared<ZEncodingImpl>(b.GetParams(), b, sfBigFP, paramsTwice);
    RPolynomial valuesTwice = ZEncodingImpl::decodeR(zEncodeTwice);

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

    // Check the r encoding?
    if (msg == "S2RC" || msg == "MSB" || msg == "Low4" || msg == "ModRaise" || msg == "PSum" || msg == "Normalize" ||
        msg == "LUTR" || msg.find("MSB") != std::string::npos || msg == "Reconstructed" || msg == "Binary") {
        for (size_t i = 0; i != values.getCoefficients().size(); ++i) {
            std::cout << msg << "  values [" << i << "]: " << values[i].toHexString(16) << std::endl;
        }
    }
    if (msg == "Rotate") {
        auto cSlots = values.toCSlots();
        for (size_t i = 0; i != cSlots.getSlots().size(); ++i) {
            std::cout << msg << "  complexValues [" << i << "]: " << cSlots[i].toHexString(16) << std::endl;
        }
    }
    if (msg == "Upper" || msg == "C2S" || msg == "LUT") {
        auto cSlotsTwice = valuesTwice.toCSlots();
        for (size_t i = 0; i != cSlotsTwice.getSlots().size(); ++i) {
            std::cout << msg << "  complexValuesTwice [" << i << "]: " << cSlotsTwice[i].toHexString(16) << std::endl;
        }
    }
    if (msg == "Encode" || msg == "Input" || msg == "Add" || msg == "Mult" || msg == "CMult" || msg == "Z2S2Z") {
        ZPolynomial zValues = values.toCSlots().getZPolynomial(0);
        for (size_t i = 0; i != values.getCoefficients().size(); ++i) {
            std::cout << msg << "  zValues [" << i << "]: " << zValues[i].toHexString(16) << std::endl;
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

std::array<Ciphertext<DCRTPoly>, 2> MSBBootstrap(CiphertextT ct, LeveledZ z, AdvancedZ advZ, FHEZ fheZ,
                                                 uint32_t oneHotBit) {
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
    const uint32_t limit = N / (zN_global * zSlots_global);
    for (uint32_t j = 1; j < limit; j <<= 1) {
        cc->EvalAddInPlace(raised, cc->EvalRotate(raised, j * (zN_global * zSlots_global / 2)));
    }
    // Now the message is multplied by N/(rN)
    //__heir_debug2(raised, "PSum");

    // Normalize to [-1, 1] from [-16, 16]
    // This is required by Chebyshev
    // Multiply by 2 * (rN / N / 32) * Delta, so the result is m / 16 * Delta^2
    // NOTE: two here should be changed.
    auto raisedSF = raised->GetScalingFactorBFP();
    auto NBigFP   = BigFixedPoint::positive(N);
    auto two      = BigFixedPoint::two();
    auto rNBFP    = BigFixedPoint::positive(zN_global * zSlots_global);
    auto BFP32    = BigFixedPoint::positive(32);
    // Then C2S below will multiply by rN again because of the construction of U0HatT
    BigFixedPoint normalizeFactorBFP = raisedSF * two * rNBFP / NBigFP / BFP32 / rNBFP;
    BigInteger normalizeFactor       = (normalizeFactorBFP.round().getValue()) >> normalizeFactorBFP.getLog2Scale();
    raised                           = z->EvalMultScalar(raised, normalizeFactor);
    raised->SetScalingFactorBFP(raisedSF * raisedSF);
    z->ModReduceInPlace(raised);
    //__heir_debug2(raised, "Normalize");

    // Note that there are other ways...some work first multiply by 1 / N
    // Then PartialSum
    // Then use CoeffsToSlots matrix to do the /16

    //------------------------------------------------------------------------------
    // Running CoeffsToSlots
    //------------------------------------------------------------------------------

    auto cSlots = zN_global * zSlots_global / 2;

    auto& precomp = fheZ->GetBootPrecom(cSlots);
    // This is FFT
    //auto c2s           = fheZ->EvalCoeffsToSlots(precomp.m_U0hatTPreFFT, raised, cSlots);
    // This is LT
    auto c2s = fheZ->EvalLinearTransform(precomp.m_U0hatTPre, raised);
    z->EvalAddInPlace(c2s, Conjugate(c2s, cc->GetEvalAutomorphismKeyMap(c2s->GetKeyTag())));
    z->ModReduceInPlace(c2s);
    //__heir_debug2(c2s, "C2S");

    //------------------------------------------------------------------------------
    // Running Approximate Mod Reduction
    //------------------------------------------------------------------------------

    auto& coeff_exp = coeff_exp_16_big_complex_46;
    auto res        = advZ->EvalChebyshevSeriesPS(c2s, coeff_exp);

    //// Double angle-iterations to get exp(2*Pi*i*x)
    res = z->EvalMult(res, res);
    z->ModReduceInPlace(res);
    res = z->EvalMult(res, res);
    z->ModReduceInPlace(res);

    //__heir_debug2(res, "Cheby1");

    //------------------------------------------------------------------------------
    // Running LUT
    //------------------------------------------------------------------------------

    auto lutCoeffs = another_interpolate(1l << 8, 2);
    auto powers    = advZ->EvalPowers(res, lutCoeffs);
    auto lut       = advZ->EvalPolyWithPrecomp(powers, lutCoeffs);
    //__heir_debug2(lut, "LUT");

    //------------------------------------------------------------------------------
    // Masking and Rotating
    //------------------------------------------------------------------------------

    auto lutElementParams = lut->GetElements()[0].GetParams();
    // TODO: fix the encode One logic
    auto oneHotEncoding1 = ZEncodingImpl::encodeOneHotInC(zN_global, zSlots_global, oneHotBit, lutElementParams,
                                                          lut->GetScalingFactorBFP(), 256);
    auto oneHotEncoding2 = ZEncodingImpl::encodeOneHotInC(zN_global, zSlots_global, oneHotBit + 8, lutElementParams,
                                                          lut->GetScalingFactorBFP(), 256 * 256);

    auto masked1 = z->EvalMult(lut, oneHotEncoding1);
    masked1      = cc->EvalRotate(masked1, -8);
    z->ModReduceInPlace(masked1);

    auto masked2 = z->EvalMult(lut, oneHotEncoding2);
    masked2      = cc->EvalRotate(masked1, -16);
    z->ModReduceInPlace(masked2);

    //------------------------------------------------------------------------------
    // Running SlotsToCoeffs
    //------------------------------------------------------------------------------

    auto s2c1 = fheZ->EvalLinearTransform(precomp.m_U0Pre, masked1);
    // This is needed for sparsely packed case
    {
        z->EvalAddInPlace(s2c1, cc->EvalRotate(s2c1, zN_global * zSlots_global / 2));
    }
    z->ModReduceInPlace(s2c1);

    auto s2c2 = fheZ->EvalLinearTransform(precomp.m_U0Pre, masked2);
    // This is needed for sparsely packed case
    {
        z->EvalAddInPlace(s2c2, cc->EvalRotate(s2c2, zN_global * zSlots_global / 2));
    }
    z->ModReduceInPlace(s2c2);
    return {s2c1, s2c2};
}

Ciphertext<DCRTPoly> ToBottom(CiphertextT ct, LeveledZ z) {
    auto q     = ct->GetElements()[0].GetModulus();
    auto sfNow = ct->GetScalingFactorBFP();
    auto qBFP  = BigFixedPoint(q, 0, false).scaleTo(128);
    auto two   = BigFixedPoint::two();
    // q / (2 * Delta)
    auto div       = (qBFP / sfNow / two).round();
    auto divScalar = div.getValue() >> div.getLog2Scale();
    auto ct2       = z->EvalMultScalar(ct, divScalar);
    ct2->SetScalingFactorBFP(qBFP / two);
    // Reduce all the way to the bottom
    z->ModReduceInPlace(ct2, ct2->GetElements()[0].GetNumOfElements() - 1);
    // To make sure the scaling factor is exactly q0 / 2
    auto q0     = ct2->GetElements()[0].GetModulus();
    auto q0BFP  = BigFixedPoint(q0, 0, false).scaleTo(128);
    auto sfNow2 = q0BFP / two;
    ct2->SetScalingFactorBFP(sfNow2);
    return ct2;
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

    uint32_t zN     = 32;
    uint32_t zSlots = 1;
    zN_global       = zN;
    zSlots_global   = zSlots;

    auto cSlots  = zN * zSlots / 2;
    auto cSlots2 = cSlots * 2;
    // For regular encoding
    DiscreteFourierTransformBigComplex::Initialize(cSlots * 4, cSlots);
    // For encoding of bootstrapping related plaintext for sparse bootstrapping
    DiscreteFourierTransformBigComplex::Initialize(cSlots2 * 4, cSlots2);
    ZLinearTransform::Initialize(zN);

    fheZ->EvalBootstrapSetup(*cc, zN * zSlots / 2, {1, 1});
    auto precom = fheZ->GetBootPrecom(cSlots);

    cc_global    = cc;
    pk_global    = keyPair.publicKey;
    sk_global    = keyPair.secretKey;
    slots_global = 32 / 2;

    auto sf = BigFixedPoint(BigInteger(1) << dcrtBits, 0, false).scaleTo(128);

    auto zero      = EncryptZero(keyPair.publicKey);
    auto elemParam = zero->GetElements()[0].GetParams();
    //
    //__heir_debug2(zero, "Input");

    RPolynomial value1 = ZPolynomial::encode(32, -1).toCSlots().toRPolynomial();
    Plaintext ptxt1    = ZEncodingImpl::encodeR(value1, elemParam, sf);

    RPolynomial value2 = ZPolynomial::encode(32, 1).toCSlots().toRPolynomial();
    Plaintext ptxt2    = ZEncodingImpl::encodeR(value2, elemParam, sf);

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
    Plaintext tPtxt = ZEncodingImpl::encodeR(t, elemParam, sf);

    if (0) {
        auto ctMul = z->EvalMult(encoded, tPtxt);
        __heir_debug2(ctMul, "Mult");
        z->ModReduceInPlace(ctMul);
        __heir_debug2(ctMul, "Mult");
    }

    /// TEST CT-CT-MULT
    Ciphertext<DCRTPoly> ct;
    if (1) {
        auto ctMul = zEvalMultFull(z, encoded, encoded2, tPtxt);
        z->ModReduceInPlace(ctMul, 2);
        __heir_debug2(ctMul, "CMult");
        ct = ctMul;
    }

    ct = fheZ->EvalArithToArithHigh(ct);
    __heir_debug2(ct, "CMult");

    return;

    /// TEST Rotate
    if (0) {
        auto ctRot = cc->EvalRotate(encoded, 1);

        __heir_debug2(ctRot, "Rotate");
    }

    /// TEST ZCoeffToSlots
    Ciphertext<DCRTPoly> zC2S;
    if (1) {
        zC2S = fheZ->EvalZLinearTransform(precom.m_ZUInversePre, ct);
        z->EvalAddInPlace(zC2S, z->EvalConjugateInC(zC2S));
        z->ModReduceInPlace(zC2S);

        //__heir_debug2(zC2S, "Upper");
    }

    // TEST SlotsToRCoeffs
    Ciphertext<DCRTPoly> s2rc;
    if (1) {
        // This is FFT
        //auto z2S2r = fheZ->EvalSlotsToCoeffs(precom.m_U0PreFFT, zC2S, zN * zSlots / 2);
        // This is LT
        auto z2S2r = fheZ->EvalLinearTransform(precom.m_U0Pre, zC2S);
        // TODO: fix the scaling
        // This is needed for sparsely packed case
        {
            z->EvalAddInPlace(z2S2r, cc->EvalRotate(z2S2r, zN * zSlots / 2));
        }
        z->ModReduceInPlace(z2S2r);
        s2rc = z2S2r;
        __heir_debug2(z2S2r, "S2RC");
    }

    // TEST Scale to MSB
    Ciphertext<DCRTPoly> ctMSB = ToBottom(s2rc, z);
    __heir_debug2(ctMSB, "MSB");

    auto [ct8, ct8_2] = MSBBootstrap(ctMSB, z, advZ, fheZ, 8);
    __heir_debug2(ct8, "MSB");

    auto ctMSBSub8      = z->EvalSub(ctMSB, ToBottom(ct8, z));
    auto [ct16, ct16_2] = MSBBootstrap(ctMSBSub8, z, advZ, fheZ, 16);
    __heir_debug2(ct16, "MSB");

    auto ctMSBSub16     = z->EvalSub(z->EvalSub(ctMSB, ToBottom(ct8_2, z)), ToBottom(ct16, z));
    auto [ct24, ct24_2] = MSBBootstrap(ctMSBSub16, z, advZ, fheZ, 24);
    __heir_debug2(ct24, "MSB");

    auto ctMSBSub24     = z->EvalSub(z->EvalSub(ctMSB, ToBottom(ct16_2, z)), ToBottom(ct24, z));
    auto [ct32, ct32_2] = MSBBootstrap(ctMSBSub24, z, advZ, fheZ, 32);
    __heir_debug2(ct32, "MSB");
}
