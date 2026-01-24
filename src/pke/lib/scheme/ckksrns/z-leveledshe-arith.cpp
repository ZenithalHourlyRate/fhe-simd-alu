#include "scheme/ckksrns/z-leveledshe.h"
#include "encoding/z-encoding.h"

namespace lbcrypto {

//
// Operations in Z
//

Ciphertext<DCRTPoly> LeveledZImpl::EvalAddInZ(ConstCiphertext<DCRTPoly> ct, BigInteger ptxt) {
    // Encode ptxt in Z encoding
    auto zN        = ct->GetZEncodingParams().getZN();
    auto elemParam = ct->GetElements()[0].GetParams();
    auto sf        = ct->GetScalingFactorBFP();
    auto plaintext = ZEncodingImpl::encodeArith(ptxt.ConvertToInt(), zN, elemParam, sf);
    return EvalAdd(ct, plaintext);
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalMultInZ(ConstCiphertext<DCRTPoly> ct, BigInteger ptxt) {
    // Encode ptxt in Z encoding
    auto zN        = ct->GetZEncodingParams().getZN();
    auto elemParam = ct->GetElements()[0].GetParams();
    auto sf        = ct->GetScalingFactorBFP();
    // encodeBinary here is crucial: making it become MultShort
    auto plaintext =
        ZEncodingImpl::encodeZ({ZPolynomial::encodeBinary(zN, ptxt.ConvertToInt())}, zN, 1, 0, elemParam, sf);
    // This is actually MultShort
    return EvalMult(ct, plaintext);
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalMultShortInZ(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    // Each input message is scaled by zSlots (maybe?)
    auto ct = EvalMultWithAdjust(ct1, ct2);
    // Now the message is scaled by zSlots^2 (maybe?)
    return ct;
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalMultTInZ(ConstCiphertext<DCRTPoly> ct) {
    // Multiply by tPtxt
    auto zN         = ct->GetZEncodingParams().getZN();
    auto elemParam  = ct->GetElements()[0].GetParams();
    auto sf         = ct->GetScalingFactorBFP();
    Plaintext tPtxt = GetTPlaintext(zN, sf, elemParam);

    return EvalMult(ct, tPtxt);
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalMultTInvInZ(ConstCiphertext<DCRTPoly> ct) {
    // Multiply by tPtxt
    auto zN            = ct->GetZEncodingParams().getZN();
    auto elemParam     = ct->GetElements()[0].GetParams();
    auto sf            = ct->GetScalingFactorBFP();
    Plaintext tInvPtxt = GetTInvPlaintext(zN, sf, elemParam);

    return EvalMult(ct, tInvPtxt);
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalMultFullInZ(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2) {
    auto ct = EvalMultWithAdjust(ct1, ct2);
    ModReduceInPlace(ct);

    // Adjust the ZDeg
    auto ct1Params = ct1->GetZEncodingParams();
    auto ct1ZDeg   = ct1->GetZEncodingParams().getZDeg();
    auto ct2ZDeg   = ct2->GetZEncodingParams().getZDeg();
    ct->SetZEncodingParams(ZEncodingParams(ZMode, ct1Params.getZN(), ct1Params.getZSlots(), ct1ZDeg + ct2ZDeg));

    // Multiply by tPtxt
    auto zN         = ct->GetZEncodingParams().getZN();
    auto elemParam  = ct->GetElements()[0].GetParams();
    auto sf         = ct->GetScalingFactorBFP();
    Plaintext tPtxt = GetTPlaintext(zN, sf, elemParam);
    return EvalMult(ct, tPtxt);
}

//
// Operations in C
//

void LeveledZImpl::EvalAddInPlaceInC(Ciphertext<DCRTPoly> ct, const BigComplex& ptxt) {
    auto elemParams = ct->GetElements()[0].GetParams();
    auto ptxtEnc    = GetBCInCPlaintext(ptxt, ct->GetScalingFactorBFP(), elemParams);
    EvalAddInPlace(ct, ptxtEnc);
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalAddInC(ConstCiphertext<DCRTPoly> ct, const BigComplex& ptxt) {
    auto elemParams = ct->GetElements()[0].GetParams();
    auto ptxtEnc    = GetBCInCPlaintext(ptxt, ct->GetScalingFactorBFP(), elemParams);
    return EvalAdd(ct, ptxtEnc);
}

void LeveledZImpl::EvalMultInPlaceInC(Ciphertext<DCRTPoly> ct, const BigComplex& ptxt, BigFixedPoint scalingFactor) {
    auto elemParams = ct->GetElements()[0].GetParams();
    if (scalingFactor.equalZero()) {
        scalingFactor = ct->GetScalingFactorBFP();
    }
    auto ptxtEnc = GetBCInCPlaintext(ptxt, scalingFactor, elemParams);
    EvalMultInPlace(ct, ptxtEnc);
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalMultInC(ConstCiphertext<DCRTPoly> ct, const BigComplex& ptxt,
                                               BigFixedPoint scalingFactor) {
    auto elemParams = ct->GetElements()[0].GetParams();
    if (scalingFactor.equalZero()) {
        scalingFactor = ct->GetScalingFactorBFP();
    }
    auto ptxtEnc = GetBCInCPlaintext(ptxt, scalingFactor, elemParams);
    return EvalMult(ct, ptxtEnc);
}

Ciphertext<DCRTPoly> LeveledZImpl::EvalConjugateInC(ConstCiphertext<DCRTPoly> ct) {
    uint32_t N = ct->GetElements()[0].GetRingDimension();
    std::vector<uint32_t> vec(N);
    PrecomputeAutoMap(N, 2 * N - 1, &vec);

    auto result = ct->Clone();

    auto cc         = ct->GetCryptoContext();
    auto algo       = ct->GetCryptoContext()->GetScheme();
    auto evalKeyMap = cc->GetEvalAutomorphismKeyMap(ct->GetKeyTag());
    algo->KeySwitchInPlace(result, evalKeyMap.at(2 * N - 1));

    auto& rcv = result->GetElements();
    rcv[0]    = rcv[0].AutomorphismTransform(2 * N - 1, vec);
    rcv[1]    = rcv[1].AutomorphismTransform(2 * N - 1, vec);
    return result;
}

}  // namespace lbcrypto