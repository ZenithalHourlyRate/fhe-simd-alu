#include "scheme/ckksrns/z-leveledshe.h"
#include "encoding/z-encode.h"

//=============================================================================
// Custom Evals
//=============================================================================

namespace lbcrypto {

//=============================================================================
// Generic methods
//=============================================================================

void gEvalAddInPlace(Ciphertext<DCRTPoly> ct, Plaintext ptxt) {
    if (ct->GetElements()[0].GetParams()->GetParams().size() !=
        ptxt->GetElement<DCRTPoly>().GetParams()->GetParams().size()) {
        OPENFHE_THROW("Ciphertext and Plaintext size mismatch in EvalAdd");
    }
    if (ptxt->GetElement<DCRTPoly>().GetFormat() != Format::EVALUATION) {
        OPENFHE_THROW("Plaintext must be in EVALUATION format in EvalMultDCRTPoly");
    }
    ZEncoding zEnc = std::dynamic_pointer_cast<ZEncodingImpl>(ptxt);
    if (!zEnc->GetScalingFactorBFP().almostEqual(ct->GetScalingFactorBFP())) {
        OPENFHE_THROW("Ciphertext and Plaintext scaling factor mismatch in EvalAdd");
    }
    auto zEncDCRTPoly = zEnc->GetElement<DCRTPoly>();

    auto& b = ct->GetElements()[0];
    b += zEncDCRTPoly;
}

Ciphertext<DCRTPoly> gEvalAdd(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt) {
    auto ctNew = ct->Clone();
    gEvalAddInPlace(ctNew, ptxt);
    return ctNew;
}

void gEvalAddInPlace(Ciphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2) {
    if (ct->GetElements()[0].GetParams()->GetParams().size() != ct2->GetElements()[0].GetParams()->GetParams().size()) {
        OPENFHE_THROW("Ciphertext size mismatch in EvalAddInplace");
    }
    if (!ct->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
        OPENFHE_THROW("Ciphertext scaling factor mismatch in EvalAddInplace");
    }
    auto& cv1  = ct->GetElements();
    auto& cv2  = ct2->GetElements();
    uint32_t n = cv1.size();
    for (uint32_t i = 0; i < n; ++i)
        cv1[i] += cv2[i];
}

Ciphertext<DCRTPoly> gEvalAdd(ConstCiphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2) {
    auto ctNew = ct->Clone();
    gEvalAddInPlace(ctNew, ct2);
    return ctNew;
}

void gEvalSubInPlace(Ciphertext<DCRTPoly> ct, Plaintext ptxt) {
    if (ct->GetElements()[0].GetParams()->GetParams().size() !=
        ptxt->GetElement<DCRTPoly>().GetParams()->GetParams().size()) {
        OPENFHE_THROW("Ciphertext and Plaintext size mismatch in EvalAdd");
    }
    if (ptxt->GetElement<DCRTPoly>().GetFormat() != Format::EVALUATION) {
        OPENFHE_THROW("Plaintext must be in EVALUATION format in EvalMultDCRTPoly");
    }
    ZEncoding zEnc = std::dynamic_pointer_cast<ZEncodingImpl>(ptxt);
    if (!zEnc->GetScalingFactorBFP().almostEqual(ct->GetScalingFactorBFP())) {
        OPENFHE_THROW("Ciphertext and Plaintext scaling factor mismatch in EvalAdd");
    }
    auto zEncDCRTPoly = zEnc->GetElement<DCRTPoly>();

    auto& b = ct->GetElements()[0];
    b -= zEncDCRTPoly;
}

Ciphertext<DCRTPoly> gEvalSub(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt) {
    auto ctNew = ct->Clone();
    gEvalAddInPlace(ctNew, ptxt);
    return ctNew;
}

void gEvalSubInPlace(Ciphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2) {
    if (ct->GetElements()[0].GetParams()->GetParams().size() != ct2->GetElements()[0].GetParams()->GetParams().size()) {
        OPENFHE_THROW("Ciphertext size mismatch in EvalAddInplace");
    }
    if (!ct->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
        OPENFHE_THROW("Ciphertext scaling factor mismatch in EvalAddInplace");
    }
    auto& cv1  = ct->GetElements();
    auto& cv2  = ct2->GetElements();
    uint32_t n = cv1.size();
    for (uint32_t i = 0; i < n; ++i)
        cv1[i] -= cv2[i];
}

Ciphertext<DCRTPoly> gEvalSub(ConstCiphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2) {
    auto ctNew = ct->Clone();
    gEvalSubInPlace(ctNew, ct2);
    return ctNew;
}

void gEvalMultInPlace(Ciphertext<DCRTPoly> ct, Plaintext ptxt) {
    ZEncoding zEnc    = std::dynamic_pointer_cast<ZEncodingImpl>(ptxt);
    auto zEncDCRTPoly = zEnc->GetElement<DCRTPoly>();
    if (ct->GetElements()[0].GetParams()->GetParams().size() !=
        ptxt->GetElement<DCRTPoly>().GetParams()->GetParams().size()) {
        OPENFHE_THROW("Ciphertext and Plaintext size mismatch in EvalMult");
    }
    if (zEncDCRTPoly.GetFormat() != Format::EVALUATION) {
        OPENFHE_THROW("Plaintext must be in EVALUATION format in EvalMult");
    }
    auto ctSFBFP   = ct->GetScalingFactorBFP();
    auto ptxtSFBFP = zEnc->GetScalingFactorBFP();

    for (auto& a : ct->GetElements()) {
        a *= zEncDCRTPoly;
    }
    ct->SetScalingFactorBFP(ctSFBFP * ptxtSFBFP);
}

Ciphertext<DCRTPoly> gEvalMult(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt) {
    auto ctNew = ct->Clone();
    gEvalMultInPlace(ctNew, ptxt);
    return ctNew;
}

Ciphertext<DCRTPoly> gEvalMult(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    auto sfBFP1 = ct1->GetScalingFactorBFP();
    auto sfBFP2 = ct2->GetScalingFactorBFP();
    auto cc     = ct1->GetCryptoContext();
    // We use cc here for automatic relinearization
    // But we do not rely on automatic rescaling
    auto ctNew = cc->EvalMult(ct1, ct2);
    ctNew->SetScalingFactorBFP(sfBFP1 * sfBFP2);
    return ctNew;
}

Ciphertext<DCRTPoly> gEvalMultWithAdjust(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    // GetLevel() at qL is 0...
    // GetLevel() at q0 is L...
    if (ct1->GetLevel() < ct2->GetLevel()) {
        auto ct1Adjusted = gAdjustCiphertext(ct1, ct2);
        return gEvalMult(ct1Adjusted, ct2);
    }
    else if (ct1->GetLevel() > ct2->GetLevel()) {
        auto ct2Adjusted = gAdjustCiphertext(ct2, ct1);
        return gEvalMult(ct1, ct2Adjusted);
    }
    else {
        if (!ct1->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
            OPENFHE_THROW("Scaling factors are not equal in gEvalMultWithAdjust");
        }
        return gEvalMult(ct1, ct2);
    }
}

Ciphertext<DCRTPoly> gEvalAddWithAdjust(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    // GetLevel() at qL is 0...
    // GetLevel() at q0 is L...
    if (ct1->GetLevel() < ct2->GetLevel()) {
        auto ct1Adjusted = gAdjustCiphertext(ct1, ct2);
        return gEvalAdd(ct1Adjusted, ct2);
    }
    else if (ct1->GetLevel() > ct2->GetLevel()) {
        auto ct2Adjusted = gAdjustCiphertext(ct2, ct1);
        return gEvalAdd(ct1, ct2Adjusted);
    }
    else {
        if (!ct1->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
            OPENFHE_THROW("Scaling factors are not equal in gEvalAddtWithAdjust");
        }
        return gEvalAdd(ct1, ct2);
    }
}

void gEvalAddWithAdjustInPlace(Ciphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    if (ct2->GetLevel() < ct1->GetLevel()) {
        auto ct2Adjusted = gAdjustCiphertext(ct2, ct1);
        gEvalAddInPlace(ct1, ct2Adjusted);
    }
    else if (ct1->GetLevel() < ct2->GetLevel()) {
        auto ct1Adjusted = gAdjustCiphertext(ct1, ct2);
        ct1->SetElements(ct1Adjusted->GetElements());
        ct1->SetLevel(ct1Adjusted->GetLevel());
        ct1->SetScalingFactorBFP(ct1Adjusted->GetScalingFactorBFP());
        gEvalAddInPlace(ct1, ct2);
    }
    else {
        if (!ct1->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
            OPENFHE_THROW("Scaling factors are not equal in gEvalAddWithAdjust");
        }
        gEvalAddInPlace(ct1, ct2);
    }
}

Ciphertext<DCRTPoly> gEvalSubWithAdjust(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    // GetLevel() at qL is 0...
    // GetLevel() at q0 is L...
    if (ct1->GetLevel() < ct2->GetLevel()) {
        auto ct1Adjusted = gAdjustCiphertext(ct1, ct2);
        return gEvalSub(ct1Adjusted, ct2);
    }
    else if (ct1->GetLevel() > ct2->GetLevel()) {
        auto ct2Adjusted = gAdjustCiphertext(ct2, ct1);
        return gEvalSub(ct1, ct2Adjusted);
    }
    else {
        if (!ct1->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
            OPENFHE_THROW("Scaling factors are not equal in gEvalSubWithAdjust");
        }
        return gEvalSub(ct1, ct2);
    }
}

void gEvalSubWithAdjustInPlace(Ciphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    if (ct2->GetLevel() < ct1->GetLevel()) {
        auto ct2Adjusted = gAdjustCiphertext(ct2, ct1);
        gEvalSubInPlace(ct1, ct2Adjusted);
    }
    else if (ct1->GetLevel() < ct2->GetLevel()) {
        auto ct1Adjusted = gAdjustCiphertext(ct1, ct2);
        ct1->SetElements(ct1Adjusted->GetElements());
        ct1->SetLevel(ct1Adjusted->GetLevel());
        ct1->SetScalingFactorBFP(ct1Adjusted->GetScalingFactorBFP());
        gEvalSubInPlace(ct1, ct2);
    }
    else {
        if (!ct1->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
            OPENFHE_THROW("Scaling factors are not equal in gEvalSubWithAdjustInPlace");
        }
        gEvalSubInPlace(ct1, ct2);
    }
}

void gEvalMultScalarInPlace(Ciphertext<DCRTPoly> ct, BigInteger scalar) {
    // NOTE: the limit here requires RNS with more than 40 bits moduli
    auto limit = BigInteger(1) << 40;
    if (scalar > limit) {
        auto elemParams = ct->GetElements()[0].GetParams();
        auto scalarBFP  = BigFixedPoint(scalar, 0, false).scaleTo(128);
        RPolynomial rPoly;
        rPoly[0]       = scalarBFP;
        Plaintext ptxt = ZEncodingImpl::encodeR(rPoly, elemParams, BigFixedPoint::one());
        auto ctOldBFP  = ct->GetScalingFactorBFP();
        gEvalMultInPlace(ct, ptxt);
        // gEvalMultInPlace updates scaling factor
        // We need to restore it
        ct->SetScalingFactorBFP(ctOldBFP);
    }
    else {
        for (auto& a : ct->GetElements()) {
            a *= scalar;
        }
    }
}

Ciphertext<DCRTPoly> gEvalMultScalar(ConstCiphertext<DCRTPoly> ct, BigInteger scalar) {
    auto ctNew = ct->Clone();
    gEvalMultScalarInPlace(ctNew, scalar);
    return ctNew;
}

void gLevelReduceInPlace(Ciphertext<DCRTPoly>& ciphertext, size_t levels) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());

    auto& cv = ciphertext->GetElements();

    ciphertext->SetLevel(ciphertext->GetLevel() + levels);
    for (size_t i = 0; i < levels; ++i) {
        for (auto& dcrtpoly : cv)
            dcrtpoly.DropLastElement();
    }
}

void gModReduceInPlace(Ciphertext<DCRTPoly>& ciphertext, size_t levels) {
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
        auto ql      = cryptoParams->GetElementParams()->GetParams()[sizeQl - 1 - i]->GetModulus();
        auto qlBigFP = BigFixedPoint(ql, 0, false).scaleTo(128);
        ciphertext->SetScalingFactorBFP(ciphertext->GetScalingFactorBFP() / qlBigFP);

        // Old code from OpenFHE
        // double modReduceFactor = cryptoParams->GetModReduceFactor(sizeQl - 1 - i);
        // ciphertext->SetScalingFactor(ciphertext->GetScalingFactor() / modReduceFactor);
    }
}

Ciphertext<DCRTPoly> gAdjustCiphertext(ConstCiphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ctTarget) {
    auto ctBFP       = ct->GetScalingFactorBFP();
    auto ctTargetBFP = ctTarget->GetScalingFactorBFP();

    auto ctBFPLog2       = std::log2(ctBFP.convertToDouble());
    auto ctTargetBFPLog2 = std::log2(ctBFP.convertToDouble());
    // The case of Noise Deg = 2 is not handled now.
    // Should track noise degree...
    if (ctBFPLog2 > 100 || ctTargetBFPLog2 > 100) {
        OPENFHE_THROW("Can not Adjust Ciphertext with large scaling factor");
    }

    auto sizeQl       = ct->GetElements()[0].GetNumOfElements();
    auto sizeQlTarget = ctTarget->GetElements()[0].GetNumOfElements();
    if (sizeQl == sizeQlTarget) {
        if (!ctBFP.almostEqual(ctTargetBFP)) {
            OPENFHE_THROW("Can not Adjust Ciphertext");
        }
        return ct->Clone();
    }
    if (sizeQl < sizeQlTarget) {
        OPENFHE_THROW("Can not Adjust Ciphertext to larger size");
    }
    auto ctNew = ct->Clone();

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ct->GetCryptoParameters());
    auto qlTargetPlusOne    = cryptoParams->GetElementParams()->GetParams()[sizeQlTarget]->GetModulus();
    auto qlTargetPlusOneBFP = BigFixedPoint(qlTargetPlusOne, 0, false).scaleTo(128);

    gLevelReduceInPlace(ctNew, sizeQl - sizeQlTarget - 1);

    auto adjustFactorBFP = ctTargetBFP * qlTargetPlusOneBFP / ctBFP;
    auto adjustFactor    = adjustFactorBFP.round().getValue() >> adjustFactorBFP.getLog2Scale();
    gEvalMultScalarInPlace(ctNew, adjustFactor);
    gModReduceInPlace(ctNew);
    ctNew->SetScalingFactorBFP(ctTargetBFP);
    return ctNew;
}

//=============================================================================
// Zencoding methods
//=============================================================================

Ciphertext<DCRTPoly> zEvalMultShort(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    return gEvalMult(ct1, ct2);
}
Ciphertext<DCRTPoly> zEvalMultFull(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2, Plaintext t) {
    return gEvalMult(gEvalMult(ct1, ct2), t);
}

Ciphertext<DCRTPoly> zEvalAdd(ConstCiphertext<DCRTPoly> ct, uint32_t ptxt) {
    auto elemParam     = ct->GetElements()[0].GetParams();
    auto sf            = ct->GetScalingFactorBFP();
    RPolynomial value1 = ZPolynomial::encode(ptxt).toRPolynomial();
    Plaintext ptxt1    = ZEncodingImpl::encodeR(value1, elemParam, sf);
    return gEvalAdd(ct, ptxt1);
}

Ciphertext<DCRTPoly> zEvalMult(ConstCiphertext<DCRTPoly> ct, uint32_t ptxt, BigFixedPoint scalingFactor) {
    auto elemParam     = ct->GetElements()[0].GetParams();
    auto sf            = ct->GetScalingFactorBFP();
    RPolynomial value1 = ZPolynomial::encode(ptxt).toRPolynomial();
    if (scalingFactor.equalZero()) {
        scalingFactor = sf;
    }
    Plaintext ptxt1 = ZEncodingImpl::encodeR(value1, elemParam, scalingFactor);
    return gEvalMult(ct, ptxt1);
}

void cEvalAddInPlace(Ciphertext<DCRTPoly> ct, BigComplex ptxt) {
    auto elemParam = ct->GetElements()[0].GetParams();
    auto sf        = ct->GetScalingFactorBFP();
    // TODO: remove the 16 requirement
    CSlots cslots(std::vector<BigComplex>(16, ptxt));
    RPolynomial value1 = cslots.toRPolynomial();
    Plaintext ptxt1    = ZEncodingImpl::encodeR(value1, elemParam, sf);
    gEvalAddInPlace(ct, ptxt1);
}

Ciphertext<DCRTPoly> cEvalAdd(ConstCiphertext<DCRTPoly> ct, BigComplex ptxt) {
    auto ctNew = ct->Clone();
    cEvalAddInPlace(ctNew, ptxt);
    return ctNew;
}

Ciphertext<DCRTPoly> cEvalMult(ConstCiphertext<DCRTPoly> ct, BigComplex ptxt, BigFixedPoint scalingFactor) {
    auto elemParam = ct->GetElements()[0].GetParams();
    if (scalingFactor.equalZero()) {
        scalingFactor = ct->GetScalingFactorBFP();
    }
    // TODO: remove the 16 requirement
    CSlots cslots(std::vector<BigComplex>(16, ptxt));
    //Is not this just a constant? We can optimize it.
    RPolynomial value1 = cslots.toRPolynomial();
    Plaintext ptxt1    = ZEncodingImpl::encodeR(value1, elemParam, scalingFactor);
    return gEvalMult(ct, ptxt1);
}

//=============================================================================
// Generic methods
//=============================================================================

void LeveledZImpl::EvalAddInPlace(Ciphertext<DCRTPoly> ct, Plaintext ptxt) {
    if (ct->GetElements()[0].GetParams()->GetParams().size() !=
        ptxt->GetElement<DCRTPoly>().GetParams()->GetParams().size()) {
        OPENFHE_THROW("Ciphertext and Plaintext size mismatch in EvalAdd");
    }
    if (ptxt->GetElement<DCRTPoly>().GetFormat() != Format::EVALUATION) {
        OPENFHE_THROW("Plaintext must be in EVALUATION format in EvalMultDCRTPoly");
    }
    ZEncoding zEnc = std::dynamic_pointer_cast<ZEncodingImpl>(ptxt);
    if (!zEnc->GetScalingFactorBFP().almostEqual(ct->GetScalingFactorBFP())) {
        OPENFHE_THROW("Ciphertext and Plaintext scaling factor mismatch in EvalAdd");
    }
    auto zEncDCRTPoly = zEnc->GetElement<DCRTPoly>();

    auto& b = ct->GetElements()[0];
    b += zEncDCRTPoly;
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalAdd(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt) {
    auto ctNew = ct->Clone();
    EvalAddInPlace(ctNew, ptxt);
    return ctNew;
}

void LeveledZImpl::EvalAddInPlace(Ciphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2) {
    if (ct->GetElements()[0].GetParams()->GetParams().size() != ct2->GetElements()[0].GetParams()->GetParams().size()) {
        OPENFHE_THROW("Ciphertext size mismatch in EvalAddInplace");
    }
    if (!ct->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
        OPENFHE_THROW("Ciphertext scaling factor mismatch in EvalAddInplace");
    }
    auto& cv1  = ct->GetElements();
    auto& cv2  = ct2->GetElements();
    uint32_t n = cv1.size();
    for (uint32_t i = 0; i < n; ++i)
        cv1[i] += cv2[i];
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalAdd(ConstCiphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2) {
    auto ctNew = ct->Clone();
    EvalAddInPlace(ctNew, ct2);
    return ctNew;
}

void LeveledZImpl::EvalSubInPlace(Ciphertext<DCRTPoly> ct, Plaintext ptxt) {
    if (ct->GetElements()[0].GetParams()->GetParams().size() !=
        ptxt->GetElement<DCRTPoly>().GetParams()->GetParams().size()) {
        OPENFHE_THROW("Ciphertext and Plaintext size mismatch in EvalAdd");
    }
    if (ptxt->GetElement<DCRTPoly>().GetFormat() != Format::EVALUATION) {
        OPENFHE_THROW("Plaintext must be in EVALUATION format in EvalMultDCRTPoly");
    }
    ZEncoding zEnc = std::dynamic_pointer_cast<ZEncodingImpl>(ptxt);
    if (!zEnc->GetScalingFactorBFP().almostEqual(ct->GetScalingFactorBFP())) {
        OPENFHE_THROW("Ciphertext and Plaintext scaling factor mismatch in EvalAdd");
    }
    auto zEncDCRTPoly = zEnc->GetElement<DCRTPoly>();

    auto& b = ct->GetElements()[0];
    b -= zEncDCRTPoly;
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalSub(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt) {
    auto ctNew = ct->Clone();
    EvalSubInPlace(ctNew, ptxt);
    return ctNew;
}

void LeveledZImpl::EvalSubInPlace(Ciphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2) {
    if (ct->GetElements()[0].GetParams()->GetParams().size() != ct2->GetElements()[0].GetParams()->GetParams().size()) {
        OPENFHE_THROW("Ciphertext size mismatch in EvalAddInplace");
    }
    if (!ct->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
        OPENFHE_THROW("Ciphertext scaling factor mismatch in EvalAddInplace");
    }
    auto& cv1  = ct->GetElements();
    auto& cv2  = ct2->GetElements();
    uint32_t n = cv1.size();
    for (uint32_t i = 0; i < n; ++i)
        cv1[i] -= cv2[i];
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalSub(ConstCiphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2) {
    auto ctNew = ct->Clone();
    EvalSubInPlace(ctNew, ct2);
    return ctNew;
}

void LeveledZImpl::EvalMultInPlace(Ciphertext<DCRTPoly> ct, Plaintext ptxt) {
    ZEncoding zEnc    = std::dynamic_pointer_cast<ZEncodingImpl>(ptxt);
    auto zEncDCRTPoly = zEnc->GetElement<DCRTPoly>();
    if (ct->GetElements()[0].GetParams()->GetParams().size() !=
        ptxt->GetElement<DCRTPoly>().GetParams()->GetParams().size()) {
        OPENFHE_THROW("Ciphertext and Plaintext size mismatch in EvalMult");
    }
    if (zEncDCRTPoly.GetFormat() != Format::EVALUATION) {
        OPENFHE_THROW("Plaintext must be in EVALUATION format in EvalMult");
    }
    auto ctSFBFP   = ct->GetScalingFactorBFP();
    auto ptxtSFBFP = zEnc->GetScalingFactorBFP();

    for (auto& a : ct->GetElements()) {
        a *= zEncDCRTPoly;
    }
    ct->SetScalingFactorBFP(ctSFBFP * ptxtSFBFP);
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalMult(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt) {
    auto ctNew = ct->Clone();
    gEvalMultInPlace(ctNew, ptxt);
    return ctNew;
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalMult(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    auto sfBFP1 = ct1->GetScalingFactorBFP();
    auto sfBFP2 = ct2->GetScalingFactorBFP();
    auto cc     = ct1->GetCryptoContext();
    // We use cc here for automatic relinearization
    // But we do not rely on automatic rescaling
    // TODO: remove use of cc
    auto ctNew = cc->EvalMult(ct1, ct2);
    ctNew->SetScalingFactorBFP(sfBFP1 * sfBFP2);
    return ctNew;
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalMultWithAdjust(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    // GetLevel() at qL is 0...
    // GetLevel() at q0 is L...
    if (ct1->GetLevel() < ct2->GetLevel()) {
        auto ct1Adjusted = gAdjustCiphertext(ct1, ct2);
        return gEvalMult(ct1Adjusted, ct2);
    }
    else if (ct1->GetLevel() > ct2->GetLevel()) {
        auto ct2Adjusted = gAdjustCiphertext(ct2, ct1);
        return gEvalMult(ct1, ct2Adjusted);
    }
    else {
        if (!ct1->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
            OPENFHE_THROW("Scaling factors are not equal in gEvalMultWithAdjust");
        }
        return gEvalMult(ct1, ct2);
    }
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalAddWithAdjust(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    // GetLevel() at qL is 0...
    // GetLevel() at q0 is L...
    if (ct1->GetLevel() < ct2->GetLevel()) {
        auto ct1Adjusted = gAdjustCiphertext(ct1, ct2);
        return gEvalAdd(ct1Adjusted, ct2);
    }
    else if (ct1->GetLevel() > ct2->GetLevel()) {
        auto ct2Adjusted = gAdjustCiphertext(ct2, ct1);
        return gEvalAdd(ct1, ct2Adjusted);
    }
    else {
        if (!ct1->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
            OPENFHE_THROW("Scaling factors are not equal in gEvalAddtWithAdjust");
        }
        return gEvalAdd(ct1, ct2);
    }
}

void LeveledZImpl::EvalAddWithAdjustInPlace(Ciphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    if (ct2->GetLevel() < ct1->GetLevel()) {
        auto ct2Adjusted = gAdjustCiphertext(ct2, ct1);
        EvalAddInPlace(ct1, ct2Adjusted);
    }
    else if (ct1->GetLevel() < ct2->GetLevel()) {
        auto ct1Adjusted = gAdjustCiphertext(ct1, ct2);
        ct1->SetElements(ct1Adjusted->GetElements());
        ct1->SetLevel(ct1Adjusted->GetLevel());
        ct1->SetScalingFactorBFP(ct1Adjusted->GetScalingFactorBFP());
        EvalAddInPlace(ct1, ct2);
    }
    else {
        if (!ct1->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
            OPENFHE_THROW("Scaling factors are not equal in gEvalAddWithAdjust");
        }
        EvalAddInPlace(ct1, ct2);
    }
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalSubWithAdjust(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    // GetLevel() at qL is 0...
    // GetLevel() at q0 is L...
    if (ct1->GetLevel() < ct2->GetLevel()) {
        auto ct1Adjusted = gAdjustCiphertext(ct1, ct2);
        return EvalSub(ct1Adjusted, ct2);
    }
    else if (ct1->GetLevel() > ct2->GetLevel()) {
        auto ct2Adjusted = gAdjustCiphertext(ct2, ct1);
        return EvalSub(ct1, ct2Adjusted);
    }
    else {
        if (!ct1->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
            OPENFHE_THROW("Scaling factors are not equal in gEvalSubWithAdjust");
        }
        return EvalSub(ct1, ct2);
    }
}

void LeveledZImpl::EvalSubWithAdjustInPlace(Ciphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    if (ct2->GetLevel() < ct1->GetLevel()) {
        auto ct2Adjusted = AdjustCiphertext(ct2, ct1);
        EvalSubInPlace(ct1, ct2Adjusted);
    }
    else if (ct1->GetLevel() < ct2->GetLevel()) {
        auto ct1Adjusted = AdjustCiphertext(ct1, ct2);
        ct1->SetElements(ct1Adjusted->GetElements());
        ct1->SetLevel(ct1Adjusted->GetLevel());
        ct1->SetScalingFactorBFP(ct1Adjusted->GetScalingFactorBFP());
        EvalSubInPlace(ct1, ct2);
    }
    else {
        if (!ct1->GetScalingFactorBFP().almostEqual(ct2->GetScalingFactorBFP())) {
            OPENFHE_THROW("Scaling factors are not equal in gEvalSubWithAdjustInPlace");
        }
        EvalSubInPlace(ct1, ct2);
    }
}

void LeveledZImpl::EvalMultScalarInPlace(Ciphertext<DCRTPoly> ct, BigInteger scalar) {
    // NOTE: the limit here requires RNS with more than 40 bits moduli
    auto limit = BigInteger(1) << 40;
    if (scalar > limit) {
        auto elemParams = ct->GetElements()[0].GetParams();
        auto scalarBFP  = BigFixedPoint(scalar, 0, false).scaleTo(128);
        RPolynomial rPoly;
        rPoly[0]       = scalarBFP;
        Plaintext ptxt = ZEncodingImpl::encodeR(rPoly, elemParams, BigFixedPoint::one());
        auto ctOldBFP  = ct->GetScalingFactorBFP();
        EvalMultInPlace(ct, ptxt);
        // gEvalMultInPlace updates scaling factor
        // We need to restore it
        ct->SetScalingFactorBFP(ctOldBFP);
    }
    else {
        for (auto& a : ct->GetElements()) {
            a *= scalar;
        }
    }
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalMultScalar(ConstCiphertext<DCRTPoly> ct, BigInteger scalar) {
    auto ctNew = ct->Clone();
    gEvalMultScalarInPlace(ctNew, scalar);
    return ctNew;
}

void LeveledZImpl::LevelReduceInPlace(Ciphertext<DCRTPoly>& ciphertext, size_t levels) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());

    auto& cv = ciphertext->GetElements();

    ciphertext->SetLevel(ciphertext->GetLevel() + levels);
    for (size_t i = 0; i < levels; ++i) {
        for (auto& dcrtpoly : cv)
            dcrtpoly.DropLastElement();
    }
}

void LeveledZImpl::ModReduceInPlace(Ciphertext<DCRTPoly>& ciphertext, size_t levels) {
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
        auto ql      = cryptoParams->GetElementParams()->GetParams()[sizeQl - 1 - i]->GetModulus();
        auto qlBigFP = BigFixedPoint(ql, 0, false).scaleTo(128);
        ciphertext->SetScalingFactorBFP(ciphertext->GetScalingFactorBFP() / qlBigFP);

        // Old code from OpenFHE
        // double modReduceFactor = cryptoParams->GetModReduceFactor(sizeQl - 1 - i);
        // ciphertext->SetScalingFactor(ciphertext->GetScalingFactor() / modReduceFactor);
    }
}

Ciphertext<DCRTPoly> LeveledZImpl::AdjustCiphertext(ConstCiphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ctTarget) {
    auto ctBFP       = ct->GetScalingFactorBFP();
    auto ctTargetBFP = ctTarget->GetScalingFactorBFP();

    auto ctBFPLog2       = std::log2(ctBFP.convertToDouble());
    auto ctTargetBFPLog2 = std::log2(ctBFP.convertToDouble());
    // The case of Noise Deg = 2 is not handled now.
    // Should track noise degree...
    if (ctBFPLog2 > 100 || ctTargetBFPLog2 > 100) {
        OPENFHE_THROW("Can not Adjust Ciphertext with large scaling factor");
    }

    auto sizeQl       = ct->GetElements()[0].GetNumOfElements();
    auto sizeQlTarget = ctTarget->GetElements()[0].GetNumOfElements();
    if (sizeQl == sizeQlTarget) {
        if (!ctBFP.almostEqual(ctTargetBFP)) {
            OPENFHE_THROW("Can not Adjust Ciphertext");
        }
        return ct->Clone();
    }
    if (sizeQl < sizeQlTarget) {
        OPENFHE_THROW("Can not Adjust Ciphertext to larger size");
    }
    auto ctNew = ct->Clone();

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ct->GetCryptoParameters());
    auto qlTargetPlusOne    = cryptoParams->GetElementParams()->GetParams()[sizeQlTarget]->GetModulus();
    auto qlTargetPlusOneBFP = BigFixedPoint(qlTargetPlusOne, 0, false).scaleTo(128);

    LevelReduceInPlace(ctNew, sizeQl - sizeQlTarget - 1);

    auto adjustFactorBFP = ctTargetBFP * qlTargetPlusOneBFP / ctBFP;
    auto adjustFactor    = adjustFactorBFP.round().getValue() >> adjustFactorBFP.getLog2Scale();
    EvalMultScalarInPlace(ctNew, adjustFactor);
    ModReduceInPlace(ctNew);
    ctNew->SetScalingFactorBFP(ctTargetBFP);
    return ctNew;
}

}  // namespace lbcrypto