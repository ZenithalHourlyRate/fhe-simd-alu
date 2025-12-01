#include "cryptocontext.h"
#include "scheme/ckksrns/ckksrns-utils.h"
#include "scheme/ckksrns/z-advancedshe.h"

namespace lbcrypto {

std::shared_ptr<seriesPowers<DCRTPoly>> zInternalEvalChebyPolysPS(ConstCiphertext<DCRTPoly>& x, uint32_t degree,
                                                                  double a, double b) {
    auto degs  = ComputeDegreesPS(degree);
    uint32_t k = degs[0];
    uint32_t m = degs[1];

    // computes linear transformation y = -1 + 2 (x-a)/(b-a)
    // consumes one level when a <> -1 && b <> 1
    auto cc = x->GetCryptoContext();
    std::vector<Ciphertext<DCRTPoly>> T(k);
    if (!IsNotEqualNegOne(a) && !IsNotEqualOne(b)) {
        // no linear transformation is needed if a = -1, b = 1
        // T_1(y) = y
        T[0] = x->Clone();
    }
    else {
        // linear transformation is needed
        double alpha = 2 / (b - a);
        double beta  = a * alpha;

        T[0] = cc->EvalMult(x, alpha);
        cc->ModReduceInPlace(T[0]);
        cc->EvalAddInPlace(T[0], -1.0 - beta);
    }

    // Computes Chebyshev polynomials up to degree k
    // for y: T_1(y) = y, T_2(y), ... , T_k(y)
    // uses binary tree multiplication
    for (uint32_t i = 2; i <= k; ++i) {
        if (i & 0x1) {  // if i is odd
            // compute T_{2i+1}(y) = 2*T_i(y)*T_{i+1}(y) - y
            T[i - 1] = cc->EvalMult(T[i / 2 - 1], T[i / 2]);
            cc->EvalAddInPlaceNoCheck(T[i - 1], T[i - 1]);
            cc->ModReduceInPlace(T[i - 1]);
            cc->EvalSubInPlace(T[i - 1], T[0]);
        }
        else {
            // compute T_{2i}(y) = 2*T_i(y)^2 - 1
            T[i - 1] = cc->EvalSquare(T[i / 2 - 1]);
            cc->EvalAddInPlaceNoCheck(T[i - 1], T[i - 1]);
            cc->ModReduceInPlace(T[i - 1]);
            cc->EvalAddInPlace(T[i - 1], -1.0);
        }
    }

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(T[k - 1]->GetCryptoParameters());
    if (cryptoParams->GetScalingTechnique() == FIXEDMANUAL) {
        // brings all powers of x to the same level
        for (uint32_t i = 1; i < k; ++i)
            cc->LevelReduceInPlace(T[i - 1], nullptr, T[k - 1]->GetLevel() - T[i - 1]->GetLevel());
    }
    else {
        for (uint32_t i = 1; i < k; ++i)
            cc->GetScheme()->AdjustLevelsAndDepthInPlace(T[i - 1], T[k - 1]);
    }

    std::vector<Ciphertext<DCRTPoly>> T2(m);
    // T2[0] is used as a placeholder
    T2[0] = T.back();

    // computes T_{k(2*m - 1)}(y)
    auto T2km1 = T.back();

    for (uint32_t i = 1; i < m; ++i) {
        // Compute the Chebyshev polynomials T_k(y), T_{2k}(y), T_{4k}(y), ... , T_{2^{m-1}k}(y)
        T2[i] = cc->EvalSquare(T2[i - 1]);
        cc->EvalAddInPlaceNoCheck(T2[i], T2[i]);
        cc->ModReduceInPlace(T2[i]);
        cc->EvalAddInPlace(T2[i], -1.0);

        // compute T_{k(2*m - 1)} = 2*T_{k(2^{m-1}-1)}(y)*T_{k*2^{m-1}}(y) - T_k(y)
        T2km1 = cc->EvalMult(T2km1, T2[i]);
        cc->EvalAddInPlaceNoCheck(T2km1, T2km1);
        cc->ModReduceInPlace(T2km1);
        cc->EvalSubInPlace(T2km1, T2[0]);
    }

    return std::make_shared<seriesPowers<DCRTPoly>>(std::move(T), std::move(T2), std::move(T2km1), k, m);
}

}  // namespace lbcrypto