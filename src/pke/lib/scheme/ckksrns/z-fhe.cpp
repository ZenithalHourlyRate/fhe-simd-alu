#include "scheme/ckksrns/z-fhe.h"

namespace lbcrypto {

//------------------------------------------------------------------------------
// EVALUATION: CoeffsToSlots and SlotsToCoeffs
//------------------------------------------------------------------------------

Ciphertext<DCRTPoly> FHEZImpl::EvalLinearTransform(std::vector<ZBootstrapPlaintextCache>& A,
                                                   ConstCiphertext<DCRTPoly>& ct) const {
    // Computing the baby-step bStep and the giant-step gStep.
    const uint32_t slots = A.size();
    const auto& p        = GetBootPrecom(slots);
    const uint32_t bStep = (p.m_paramsEnc.g == 0) ? std::ceil(std::sqrt(slots)) : p.m_paramsEnc.g;
    const uint32_t gStep = std::ceil(static_cast<double>(slots) / bStep);

    auto cc     = ct->GetCryptoContext();
    auto digits = cc->EvalFastRotationPrecompute(ct);

    // hoisted automorphisms
    std::vector<Ciphertext<DCRTPoly>> fastRotation(bStep - 1);
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(bStep - 1))
    for (uint32_t j = 1; j < bStep; ++j)
        fastRotation[j - 1] = cc->EvalFastRotationExt(ct, j, digits, true);

    auto elementParams = fastRotation[0]->GetElements()[0].GetParams();
    auto sfBFP         = ct->GetScalingFactorBFP();

    const uint32_t M = cc->GetCyclotomicOrder();
    const uint32_t N = cc->GetRingDimension();
    std::vector<uint32_t> map(N);
    Ciphertext<DCRTPoly> result;
    DCRTPoly first;
    for (uint32_t j = 0; j < gStep; ++j) {
        auto inner = z->EvalMult(cc->KeySwitchExt(ct, true), A[bStep * j]->GetPlaintext(sfBFP, elementParams));
        for (uint32_t i = 1; i < bStep; ++i) {
            if (bStep * j + i < slots)
                z->EvalAddInPlace(
                    inner, z->EvalMult(fastRotation[i - 1], A[bStep * j + i]->GetPlaintext(sfBFP, elementParams)));
        }

        if (j == 0) {
            first         = cc->KeySwitchDownFirstElement(inner);
            auto elements = inner->GetElements();
            elements[0].SetValuesToZero();
            inner->SetElements(std::move(elements));
            result = std::move(inner);
        }
        else {
            inner = cc->KeySwitchDown(inner);
            // Find the automorphism index that corresponds to rotation index index.
            uint32_t autoIndex = FindAutomorphismIndex2nComplex(bStep * j, M);
            PrecomputeAutoMap(N, autoIndex, &map);
            first += inner->GetElements()[0].AutomorphismTransform(autoIndex, map);

            auto&& innerDigits = cc->EvalFastRotationPrecompute(inner);
            z->EvalAddInPlace(result, cc->EvalFastRotationExt(inner, bStep * j, innerDigits, false));
        }
    }
    result = cc->KeySwitchDown(result);
    result->GetElements()[0] += first;
    return result;
}

}  // namespace lbcrypto