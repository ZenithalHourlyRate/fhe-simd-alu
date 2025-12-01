#include "scheme/ckksrns/z-leveledshe.h"
#include "encoding/z-encode.h"

//=============================================================================
// Custom Evals
//=============================================================================

namespace lbcrypto {

Ciphertext<DCRTPoly> zEvalAdd(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt) {
    assert(ct->GetElements().size() == ptxt.GetParams()->GetParams().size() &&
           "Ciphertext and Plaintext size mismatch in EvalAdd");
    assert(ptxt.GetFormat() == Format::EVALUATION && "Plaintext must be in EVALUATION format in EvalMultDCRTPoly");
    ZEncoding zEnc = std::dynamic_pointer_cast<ZEncodingImpl>(ptxt);
    assert(zEnc->GetScalingFactorBFP() == ct->GetScalingFactorBFP() &&
           "Ciphertext and Plaintext scaling factor mismatch in EvalAdd");
    auto zEncDCRTPoly = zEnc->GetElement<DCRTPoly>();

    auto ctNew = ct->Clone();
    auto& b    = ctNew->GetElements()[0];
    b += zEncDCRTPoly;
    ctNew->SetScalingFactorBFP(ct->GetScalingFactorBFP());
    return ctNew;
}

Ciphertext<DCRTPoly> zEvalMult(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt) {
    ZEncoding zEnc    = std::dynamic_pointer_cast<ZEncodingImpl>(ptxt);
    auto zEncDCRTPoly = zEnc->GetElement<DCRTPoly>();
    assert(ct->GetElements().size() == zEncDCRTPoly.GetParams()->GetParams().size() &&
           "Ciphertext and Plaintext size mismatch in EvalMult");
    assert(zEncDCRTPoly.GetFormat() == Format::EVALUATION &&
           "Plaintext must be in EVALUATION format in EvalMultDCRTPoly");
    auto ctSFBFP   = ct->GetScalingFactorBFP();
    auto ptxtSFBFP = zEnc->GetScalingFactorBFP();

    auto ctNew = ct->Clone();
    for (auto& a : ctNew->GetElements()) {
        a *= zEncDCRTPoly;
    }
    ctNew->SetScalingFactorBFP(ctSFBFP * ptxtSFBFP);
    return ctNew;
}

Ciphertext<DCRTPoly> zEvalMult(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    auto sfBFP1 = ct1->GetScalingFactorBFP();
    auto sfBFP2 = ct2->GetScalingFactorBFP();
    auto cc     = ct1->GetCryptoContext();
    // We use cc here for automatic relinearization
    // But we do not rely on automatic rescaling
    auto ctNew = cc->EvalMult(ct1, ct2);
    ctNew->SetScalingFactorBFP(sfBFP1 * sfBFP2);
    return ctNew;
}

Ciphertext<DCRTPoly> zEvalMultScalar(ConstCiphertext<DCRTPoly> ct, BigInteger scalar) {
    assert(ct->GetElements().size() == ptxt.GetParams()->GetParams().size() &&
           "Ciphertext and Plaintext size mismatch in EvalMultDCRTPoly");
    assert(ptxt.GetFormat() == Format::EVALUATION && "Plaintext must be in EVALUATION format in EvalMultDCRTPoly");
    auto ctNew = ct->Clone();
    for (auto& a : ctNew->GetElements()) {
        a *= scalar;
    }
    return ctNew;
}

void zModReduceInPlace(Ciphertext<DCRTPoly>& ciphertext, size_t levels) {
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

}  // namespace lbcrypto