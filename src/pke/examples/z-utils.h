#include <cassert>
#include "openfhe.h"
#include "math/z-encode.h"

using namespace lbcrypto;
using CiphertextT        = ConstCiphertext<DCRTPoly>;
using MutableCiphertextT = Ciphertext<DCRTPoly>;
using CCParamsT          = CCParams<CryptoContextBFVRNS>;
using CryptoContextT     = CryptoContext<DCRTPoly>;
using EvalKeyT           = EvalKey<DCRTPoly>;
using PlaintextT         = Plaintext;
using PrivateKeyT        = PrivateKey<DCRTPoly>;
using PublicKeyT         = PublicKey<DCRTPoly>;

//=============================================================================
// Encryption Utils
//=============================================================================

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

std::shared_ptr<std::vector<DCRTPoly>> EncryptZeroCore(const PublicKey<DCRTPoly> publicKey) {
    const auto cryptoParams =
        std::dynamic_pointer_cast<CryptoParametersRLWE<DCRTPoly>>(publicKey->GetCryptoParameters());

    const auto ns = cryptoParams->GetNoiseScale();

    auto elementParams = cryptoParams->GetElementParams();

    auto pk = publicKey->GetPublicElements();

    auto p0 = pk[0];
    auto p1 = pk[1];

    uint32_t sizeQ  = elementParams->GetParams().size();
    uint32_t sizePK = p0.GetParams()->GetParams().size();

    if (sizePK > sizeQ) {
        p0.DropLastElements(sizePK - sizeQ);
        p1.DropLastElements(sizePK - sizeQ);
    }

    DCRTPoly::TugType tug;
    DCRTPoly v = DCRTPoly(tug, elementParams, Format::EVALUATION);

    // noise generation with the discrete gaussian generator dgg
    auto& dgg = cryptoParams->GetDiscreteGaussianGenerator();
    DCRTPoly e0(dgg, elementParams, Format::EVALUATION);
    DCRTPoly e1(dgg, elementParams, Format::EVALUATION);

    DCRTPoly b(elementParams);
    DCRTPoly a(elementParams);

    b = p0 * v + ns * e0;
    a = p1 * v + ns * e1;

    return std::make_shared<std::vector<DCRTPoly>>(std::initializer_list<DCRTPoly>({std::move(b), std::move(a)}));
}

std::shared_ptr<std::vector<DCRTPoly>> EncryptZeroCore(const PrivateKey<DCRTPoly> privateKey) {
    const auto cryptoParams =
        std::dynamic_pointer_cast<CryptoParametersRLWE<DCRTPoly>>(privateKey->GetCryptoParameters());
    const auto elementParams = cryptoParams->GetElementParams();

    DCRTPoly::DugType dug;
    DCRTPoly a(dug, elementParams, Format::EVALUATION);

    DCRTPoly e(cryptoParams->GetDiscreteGaussianGenerator(), elementParams, Format::EVALUATION);
    NativeInteger ns = cryptoParams->GetNoiseScale();

    // {b = ns * e - a * s, a}
    DCRTPoly b(std::move((e *= ns) -= (a * privateKey->GetPrivateElement())));

    return std::make_shared<std::vector<DCRTPoly>>(std::initializer_list<DCRTPoly>({std::move(b), std::move(a)}));
}

Ciphertext<DCRTPoly> EncryptDCRTPoly(DCRTPoly ptxt, double scalingFactor, const PublicKey<DCRTPoly> publicKey) {
    auto ba = EncryptZeroCore(publicKey);
    (*ba)[0] += ptxt;

    auto ctxt = std::make_shared<CiphertextImpl<DCRTPoly>>(publicKey);
    ctxt->SetElements(std::move(*ba));
    ctxt->SetNoiseScaleDeg(1);
    ctxt->SetScalingFactor(scalingFactor);
    return ctxt;
}

Ciphertext<DCRTPoly> EncryptDCRTPoly(DCRTPoly ptxt, double scalingFactor, const PrivateKey<DCRTPoly> privateKey) {
    auto ba = EncryptZeroCore(privateKey);
    (*ba)[0] += ptxt;

    auto ctxt = std::make_shared<CiphertextImpl<DCRTPoly>>(privateKey);
    ctxt->SetElements(std::move(*ba));
    ctxt->SetNoiseScaleDeg(1);
    ctxt->SetScalingFactor(scalingFactor);
    return ctxt;
}

Ciphertext<DCRTPoly> EncryptZero(double scalingFactor, const PublicKey<DCRTPoly> publicKey) {
    return EncryptDCRTPoly(DCRTPoly(publicKey->GetCryptoParameters()->GetElementParams(), Format::EVALUATION, true),
                           scalingFactor, publicKey);
}

Ciphertext<DCRTPoly> EncryptZero(double scalingFactor, const PrivateKey<DCRTPoly> privateKey) {
    return EncryptDCRTPoly(DCRTPoly(privateKey->GetCryptoParameters()->GetElementParams(), Format::EVALUATION, true),
                           scalingFactor, privateKey);
}

Ciphertext<DCRTPoly> Conjugate(ConstCiphertext<DCRTPoly> ciphertext,
                               const std::map<uint32_t, EvalKey<DCRTPoly>>& evalKeyMap) {
    uint32_t N = ciphertext->GetElements()[0].GetRingDimension();
    std::vector<uint32_t> vec(N);
    PrecomputeAutoMap(N, 2 * N - 1, &vec);

    auto result = ciphertext->Clone();

    auto algo = ciphertext->GetCryptoContext()->GetScheme();
    algo->KeySwitchInPlace(result, evalKeyMap.at(2 * N - 1));

    auto& rcv = result->GetElements();
    rcv[0]    = rcv[0].AutomorphismTransform(2 * N - 1, vec);
    rcv[1]    = rcv[1].AutomorphismTransform(2 * N - 1, vec);
    return result;
}

//=============================================================================
// Encode Utils
//=============================================================================

DCRTPoly getDCRTPolyFromFixedPointVec(std::vector<BigFixedPoint> input,
                                      const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                      BigFixedPoint scalingFactor) {
    DCRTPoly newPoly(elementParams, Format::COEFFICIENT, true);
    auto bigPoly = newPoly.CRTInterpolate();
    auto ringDim = elementParams->GetRingDimension();
    auto q       = elementParams->GetModulus();
    for (size_t i = 0; i != input.size(); ++i) {
        auto val        = input[i] * scalingFactor;
        auto valInteger = (val.round()).getValue() >> val.getLog2Scale();
        if (val.getNeg()) {
            bigPoly[i * (ringDim / (input.size()))] = q - valInteger;
        }
        else {
            bigPoly[i * (ringDim / (input.size()))] = valInteger;
        }
    }
    DCRTPoly finalPoly(bigPoly, elementParams);
    finalPoly.SetFormat(Format::EVALUATION);
    return finalPoly;
}

std::vector<BigFixedPoint> getFixedPointVecFromDCRTPoly(DCRTPoly input, BigFixedPoint scalingFactor,
                                                        size_t outputSize) {
    input.SetFormat(Format::COEFFICIENT);
    auto bigPoly = input.CRTInterpolate();
    std::vector<BigFixedPoint> output;
    auto ringDim = input.GetParams()->GetRingDimension();
    auto q       = input.GetParams()->GetModulus();
    for (size_t i = 0; i != outputSize; ++i) {
        auto valInteger = bigPoly[i * (ringDim / outputSize)];
        bool neg        = false;
        if (valInteger > q / 2) {
            neg        = true;
            valInteger = q - valInteger;
        }
        auto valFixedPoint = BigFixedPoint(valInteger, 0, neg).scaleTo(z_upper_roots_scale);
        output.push_back(valFixedPoint / scalingFactor);
    }
    return output;
}

//=============================================================================
// Custom Evals
//=============================================================================

Ciphertext<DCRTPoly> EvalAddDCRTPoly(ConstCiphertext<DCRTPoly> ct, const DCRTPoly& ptxt) {
    assert(ct->GetElements().size() == ptxt.GetParams()->GetParams().size() &&
           "Ciphertext and Plaintext size mismatch in EvalAddDCRTPoly");
    assert(ptxt.GetFormat() == Format::EVALUATION && "Plaintext must be in EVALUATION format in EvalMultDCRTPoly");
    auto ctNew = ct->Clone();
    auto& b    = ctNew->GetElements()[0];
    b += ptxt;
    return ctNew;
}

Ciphertext<DCRTPoly> EvalMultDCRTPoly(ConstCiphertext<DCRTPoly> ct, const DCRTPoly& ptxt) {
    assert(ct->GetElements().size() == ptxt.GetParams()->GetParams().size() &&
           "Ciphertext and Plaintext size mismatch in EvalMultDCRTPoly");
    assert(ptxt.GetFormat() == Format::EVALUATION && "Plaintext must be in EVALUATION format in EvalMultDCRTPoly");
    auto ctNew = ct->Clone();
    for (auto& a : ctNew->GetElements()) {
        a *= ptxt;
    }
    return ctNew;
}

Ciphertext<DCRTPoly> EvalMultScalar(ConstCiphertext<DCRTPoly> ct, uint32_t scalar) {
    assert(ct->GetElements().size() == ptxt.GetParams()->GetParams().size() &&
           "Ciphertext and Plaintext size mismatch in EvalMultDCRTPoly");
    assert(ptxt.GetFormat() == Format::EVALUATION && "Plaintext must be in EVALUATION format in EvalMultDCRTPoly");
    auto ctNew = ct->Clone();
    for (auto& a : ctNew->GetElements()) {
        a *= BigInteger(scalar);
    }
    return ctNew;
}

void ModReduceCustomInPlace(Ciphertext<DCRTPoly>& ciphertext, size_t levels = 1) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());

    auto& cv = ciphertext->GetElements();

    size_t sizeQ  = cryptoParams->GetElementParams()->GetParams().size();
    size_t sizeQl = cv[0].GetNumOfElements();
    size_t diffQl = sizeQ - sizeQl;

    ciphertext->SetLevel(ciphertext->GetLevel() + levels);

    for (size_t i = 0; i < levels; ++i) {
        for (auto& dcrtpoly : cv)
            dcrtpoly.DropLastElementAndScale(cryptoParams->GetQlQlInvModqlDivqlModq(diffQl + i),
                                             cryptoParams->GetqlInvModq(diffQl + i));
        // We manually track scaling factor
        // double modReduceFactor = cryptoParams->GetModReduceFactor(sizeQl - 1 - i);
        // ciphertext->SetScalingFactor(ciphertext->GetScalingFactor() / modReduceFactor);
    }
}

//=============================================================================
// Linear Transform Auxiliaries
//=============================================================================

std::array<std::vector<DCRTPoly>, 2> getAuxLTDCRTPoly(BigCMatrix T, bool inverse,
                                                      const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                                      BigFixedPoint scalingFactor) {
    // Tinverse is of shape zN * (zN / 2)
    // T of is of (zN / 2) * zN
    auto halfSize = T.size();
    if (inverse) {
        halfSize = T[0].size();
    }
    auto rowOffset    = 0;
    auto columnOffset = 0;
    if (inverse) {
        rowOffset = halfSize;
    }
    else {
        columnOffset = halfSize;
    }

    std::vector<DCRTPoly> results1;
    std::vector<DCRTPoly> results2;

    for (size_t i = 0; i != halfSize; ++i) {
        auto diagonal = std::vector<BigComplex>(halfSize, BigComplex());
        for (size_t j = 0; j != halfSize; ++j) {
            diagonal[j] = T[j][(i + j) % halfSize];
        }
        auto diagonalInR = CSlots(diagonal).toRPolynomial();
        auto dcrtPoly    = getDCRTPolyFromFixedPointVec(diagonalInR.getCoefficients(), elementParams, scalingFactor);
        results1.push_back(dcrtPoly);
    }
    for (size_t i = 0; i != halfSize; ++i) {
        auto diagonal = std::vector<BigComplex>(halfSize, BigComplex());
        for (size_t j = 0; j != halfSize; ++j) {
            diagonal[j] = T[j + rowOffset][(i + j) % halfSize + columnOffset];
        }
        auto diagonalInR = CSlots(diagonal).toRPolynomial();
        auto dcrtPoly    = getDCRTPolyFromFixedPointVec(diagonalInR.getCoefficients(), elementParams, scalingFactor);
        results2.push_back(dcrtPoly);
    }
    return {results1, results2};
}

std::array<std::vector<DCRTPoly>, 2> getZCoeffToSlotsAuxDCRTPoly(
    const std::shared_ptr<typename DCRTPoly::Params>& elementParams, BigFixedPoint scalingFactor) {
    return getAuxLTDCRTPoly(getZUInverse(), true, elementParams, scalingFactor);
}

std::array<std::vector<DCRTPoly>, 2> getRCoeffToSlotsAuxDCRTPoly(
    const std::shared_ptr<typename DCRTPoly::Params>& elementParams, BigFixedPoint scalingFactor) {
    return getAuxLTDCRTPoly(getRUInverse(), true, elementParams, scalingFactor);
}

std::array<std::vector<DCRTPoly>, 2> getSlotsToZCoeffsAuxDCRTPoly(
    const std::shared_ptr<typename DCRTPoly::Params>& elementParams, BigFixedPoint scalingFactor) {
    return getAuxLTDCRTPoly(getZU(), false, elementParams, scalingFactor);
}

std::array<std::vector<DCRTPoly>, 2> getSlotsToRCoeffsAuxDCRTPoly(
    const std::shared_ptr<typename DCRTPoly::Params>& elementParams, BigFixedPoint scalingFactor) {
    return getAuxLTDCRTPoly(getRU(), false, elementParams, scalingFactor);
}

//=============================================================================
// Linear Transform Evals
//=============================================================================

// We now do not exploit sparse encoding

std::vector<Ciphertext<DCRTPoly>> CoeffsToSlots(CryptoContextT cc, CiphertextT ct,
                                                const std::array<std::vector<DCRTPoly>, 2>& auxPtxts) {
    // Halevi-Shoup
    // Peel the first loop
    auto resultUpper = EvalMultDCRTPoly(ct, auxPtxts[0][0]);
    auto resultDown  = EvalMultDCRTPoly(ct, auxPtxts[1][0]);
    auto startCt     = cc->EvalRotate(ct, 1);
    for (size_t i = 1; i != auxPtxts[0].size(); ++i) {
        auto diagonalUpper = EvalMultDCRTPoly(startCt, auxPtxts[0][i]);
        cc->EvalAddInPlace(resultUpper, diagonalUpper);
        auto diagonalDown = EvalMultDCRTPoly(startCt, auxPtxts[1][i]);
        cc->EvalAddInPlace(resultDown, diagonalDown);
        if (i + 1 < auxPtxts[0].size()) {
            //rotate one more
            startCt = cc->EvalRotate(ct, i + 1);
        }
    }
    cc->EvalAddInPlace(resultUpper, Conjugate(resultUpper, cc->GetEvalAutomorphismKeyMap(resultUpper->GetKeyTag())));
    cc->EvalAddInPlace(resultDown, Conjugate(resultDown, cc->GetEvalAutomorphismKeyMap(resultDown->GetKeyTag())));

    return {resultUpper, resultDown};
}

Ciphertext<DCRTPoly> SlotsToCoeffs(CryptoContextT cc, CiphertextT ctLeft, CiphertextT ctRight,
                                   const std::array<std::vector<DCRTPoly>, 2>& auxPtxts) {
    // Halevi-Shoup
    // Peel the first loop
    auto result = EvalMultDCRTPoly(ctLeft, auxPtxts[0][0]);
    cc->EvalAddInPlace(result, EvalMultDCRTPoly(ctRight, auxPtxts[1][0]));
    auto startCtLeft  = cc->EvalRotate(ctLeft, 1);
    auto startCtRight = cc->EvalRotate(ctRight, 1);
    for (size_t i = 1; i != auxPtxts[0].size(); ++i) {
        auto diagonalLeft = EvalMultDCRTPoly(startCtLeft, auxPtxts[0][i]);
        cc->EvalAddInPlace(result, diagonalLeft);
        auto diagonalRight = EvalMultDCRTPoly(startCtRight, auxPtxts[1][i]);
        cc->EvalAddInPlace(result, diagonalRight);
        // rotate one more
        if (i + 1 < auxPtxts[0].size()) {
            startCtLeft  = cc->EvalRotate(startCtLeft, 1);
            startCtRight = cc->EvalRotate(startCtRight, 1);
        }
    }
    return result;
}

// NOTE: should pre-compute ptxts outside
std::vector<Ciphertext<DCRTPoly>> ZCoeffsToSlots(CryptoContextT cc, CiphertextT ct) {
    auto elementParams = ct->GetElements()[0].GetParams();
    auto sfBigFP =
        BigFixedPoint(BigInteger(1) << static_cast<uint32_t>(std::log2(ct->GetScalingFactor())), 0, false).scaleTo(128);
    return CoeffsToSlots(cc, ct, getZCoeffToSlotsAuxDCRTPoly(elementParams, sfBigFP));
}

std::vector<Ciphertext<DCRTPoly>> RCoeffsToSlots(CryptoContextT cc, CiphertextT ct) {
    auto elementParams = ct->GetElements()[0].GetParams();
    auto sfBigFP =
        BigFixedPoint(BigInteger(1) << static_cast<uint32_t>(std::log2(ct->GetScalingFactor())), 0, false).scaleTo(128);
    return CoeffsToSlots(cc, ct, getRCoeffToSlotsAuxDCRTPoly(elementParams, sfBigFP));
}

Ciphertext<DCRTPoly> SlotsToZCoeffs(CryptoContextT cc, CiphertextT ctLeft, CiphertextT ctRight) {
    auto elementParams = ctLeft->GetElements()[0].GetParams();
    auto sfBigFP =
        BigFixedPoint(BigInteger(1) << static_cast<uint32_t>(std::log2(ctLeft->GetScalingFactor())), 0, false)
            .scaleTo(128);
    return SlotsToCoeffs(cc, ctLeft, ctRight, getSlotsToZCoeffsAuxDCRTPoly(elementParams, sfBigFP));
}

Ciphertext<DCRTPoly> SlotsToRCoeffs(CryptoContextT cc, CiphertextT ctLeft, CiphertextT ctRight) {
    auto elementParams = ctLeft->GetElements()[0].GetParams();
    auto sfBigFP =
        BigFixedPoint(BigInteger(1) << static_cast<uint32_t>(std::log2(ctLeft->GetScalingFactor())), 0, false)
            .scaleTo(128);
    return SlotsToCoeffs(cc, ctLeft, ctRight, getSlotsToRCoeffsAuxDCRTPoly(elementParams, sfBigFP));
}
