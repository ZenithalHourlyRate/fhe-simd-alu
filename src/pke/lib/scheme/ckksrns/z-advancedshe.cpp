#include "cryptocontext.h"
#include "scheme/ckksrns/ckksrns-utils.h"
#include "scheme/ckksrns/z-leveledshe.h"
#include "scheme/ckksrns/z-advancedshe.h"

namespace lbcrypto {

std::shared_ptr<seriesPowers<DCRTPoly>> zInternalEvalChebyPolysPS(ConstCiphertext<DCRTPoly>& x, uint32_t degree) {
    auto degs  = ComputeDegreesPS(degree);
    uint32_t k = degs[0];
    uint32_t m = degs[1];

    auto cc = x->GetCryptoContext();
    std::vector<Ciphertext<DCRTPoly>> T(k);
    // no linear transformation is needed if a = -1, b = 1
    // T_1(y) = y
    T[0] = x->Clone();

    // Computes Chebyshev polynomials up to degree k
    // for y: T_1(y) = y, T_2(y), ... , T_k(y)
    // uses binary tree multiplication
    for (uint32_t i = 2; i <= k; ++i) {
        if (i & 0x1) {  // if i is odd
            // compute T_{2i+1}(y) = 2*T_i(y)*T_{i+1}(y) - y
            if (T[i / 2 - 1]->GetLevel() != T[i / 2]->GetLevel()) {
                // This is to ensure per-level scaling factor is the same
                auto TIAdjusted = gAdjustCiphertext(T[i / 2 - 1], T[i / 2]);
                T[i - 1]        = gEvalMult(TIAdjusted, T[i / 2]);
            }
            else {
                T[i - 1] = gEvalMult(T[i / 2 - 1], T[i / 2]);
            }
            gEvalAddInPlace(T[i - 1], T[i - 1]);
            gModReduceInPlace(T[i - 1]);
            // TODO: maybe we can hoist T0Adjust
            auto T0Adjusted = gAdjustCiphertext(T[0], T[i - 1]);
            gEvalSubInPlace(T[i - 1], T0Adjusted);
        }
        else {
            // compute T_{2i}(y) = 2*T_i(y)^2 - 1
            T[i - 1] = gEvalMult(T[i / 2 - 1], T[i / 2 - 1]);
            gEvalAddInPlace(T[i - 1], T[i - 1]);
            gModReduceInPlace(T[i - 1]);
            auto one = BigFixedPoint::one();
            cEvalAddInPlace(T[i - 1], -one);
        }
    }

    // Adjust them to have same depth and scaling factor
    for (uint32_t i = 1; i < k; ++i) {
        if (T[i - 1]->GetLevel() != T[k - 1]->GetLevel()) {
            T[i - 1] = gAdjustCiphertext(T[i - 1], T[k - 1]);
        }
        else {
            assert(T[i - 1]->GetScalingFactorBFP().almostEqual(T[k - 1]->GetScalingFactorBFP()) &&
                   "Scaling factors are not equal!");
        }
    }

    std::vector<Ciphertext<DCRTPoly>> T2(m);
    // T2[0] is used as a placeholder
    T2[0] = T.back();

    // computes T_{k(2*m - 1)}(y)
    auto T2km1 = T.back();

    for (uint32_t i = 1; i < m; ++i) {
        // Compute the Chebyshev polynomials T_k(y), T_{2k}(y), T_{4k}(y), ... , T_{2^{m-1}k}(y)
        T2[i] = gEvalMult(T2[i - 1], T2[i - 1]);
        gEvalAddInPlace(T2[i], T2[i]);
        gModReduceInPlace(T2[i]);
        auto one = BigFixedPoint::one();
        cEvalAddInPlace(T2[i], -one);

        // compute T_{k(2*m - 1)} = 2*T_{k(2^{m-1}-1)}(y)*T_{k*2^{m-1}}(y) - T_k(y)
        auto T2km1Adjusted = gAdjustCiphertext(T2km1, T2[i]);
        T2km1              = gEvalMult(T2km1Adjusted, T2[i]);
        gEvalAddInPlace(T2km1, T2km1);
        gModReduceInPlace(T2km1);
        auto T20Adjusted = gAdjustCiphertext(T2[0], T2km1);
        gEvalSubInPlace(T2km1, T20Adjusted);
    }

    return std::make_shared<seriesPowers<DCRTPoly>>(std::move(T), std::move(T2), std::move(T2km1), k, m);
}

}  // namespace lbcrypto