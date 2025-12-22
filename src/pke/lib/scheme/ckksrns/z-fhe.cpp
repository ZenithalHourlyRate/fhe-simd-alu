#include "scheme/ckksrns/z-fhe.h"
#include "scheme/ckksrns/ckksrns-fhe.h"

namespace lbcrypto {

Ciphertext<DCRTPoly> FHEZImpl::EvalTruncate(ConstCiphertext<DCRTPoly>& ct) const {
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

Ciphertext<DCRTPoly> FHEZImpl::EvalModRaise(ConstCiphertext<DCRTPoly>& ct) const {
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
        OPENFHE_THROW("Other Key unsupported");
    }

    raised->SetScalingFactorBFP(ct->GetScalingFactorBFP());
    return raised;
}

Ciphertext<DCRTPoly> FHEZImpl::EvalArithToArithHigh(ConstCiphertext<DCRTPoly>& ct, uint32_t cSlots) const {
    auto cc            = ct->GetCryptoContext();
    auto precomp       = GetBootPrecom(cSlots);
    bool isLTBootstrap = (precomp.m_paramsEnc.lvlb == 1) && (precomp.m_paramsDec.lvlb == 1);
    auto N             = cc->GetRingDimension();

    //------------------------------------------------------------------------------
    // Z-To-C
    //------------------------------------------------------------------------------

    auto z2c = EvalZLinearTransform(precomp.m_ZUInversePre, ct);
    z->EvalAddInPlace(z2c, z->EvalConjugateInC(z2c));
    z->ModReduceInPlace(z2c);

    //------------------------------------------------------------------------------
    // C-To-R
    //------------------------------------------------------------------------------

    auto c2r =
        isLTBootstrap ? EvalLinearTransform(precomp.m_U0Pre, z2c) : EvalCoeffsToSlots(precomp.m_U0PreFFT, z2c, cSlots);
    // TODO: fix the scaling
    // This is needed for sparsely packed case
    {
        z->EvalAddInPlace(c2r, cc->EvalRotate(c2r, cSlots));
    }
    z->ModReduceInPlace(c2r);

    //------------------------------------------------------------------------------
    // Truncate and ModRaise
    //------------------------------------------------------------------------------

    auto truncated = EvalTruncate(c2r);
    auto raised    = EvalModRaise(truncated);

    //------------------------------------------------------------------------------
    // SPARSELY PACKED CASE
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Running PartialSum
    //------------------------------------------------------------------------------

    const uint32_t limit = N / (cSlots * 2);
    for (uint32_t j = 1; j < limit; j <<= 1) {
        cc->EvalAddInPlace(raised, cc->EvalRotate(raised, j * (cSlots)));
    }
    // Now the message is multplied by N/(rN)

    //------------------------------------------------------------------------------
    // R-To-C
    //------------------------------------------------------------------------------

    // Then R2C here will multiply by rN because of the construction of U0HatT
    auto r2c = isLTBootstrap ? EvalLinearTransform(precomp.m_U0Pre, raised) :
                               EvalSlotsToCoeffs(precomp.m_U0PreFFT, raised, cSlots);
    z->EvalAddInPlace(r2c, z->EvalConjugateInC(r2c));
    z->ModReduceInPlace(r2c);
    // Now the message is N * m

    // Normalize to m
    {
        auto sf     = r2c->GetScalingFactorBFP();
        auto NBigFP = BigFixedPoint::positive(N);
        // Then C2S below will multiply by rN again because of the construction of U0HatT
        BigFixedPoint normalizeFactorBFP = sf / NBigFP;
        BigInteger normalizeFactor       = (normalizeFactorBFP.round().getValue()) >> normalizeFactorBFP.getLog2Scale();
        r2c                              = z->EvalMultScalar(r2c, normalizeFactor);
        r2c->SetScalingFactorBFP(sf * sf);
        z->ModReduceInPlace(r2c);
    }

    //------------------------------------------------------------------------------
    // C-To-Z
    //------------------------------------------------------------------------------

    auto c2z = EvalZLinearTransform(precomp.m_ZUPre, r2c);
    // TODO: fix the scaling
    // This is needed for sparsely packed case
    {
        z->EvalAddInPlace(c2z, cc->EvalRotate(c2z, cSlots));
    }
    z->ModReduceInPlace(c2z);
    return c2z;
}

}  // namespace lbcrypto