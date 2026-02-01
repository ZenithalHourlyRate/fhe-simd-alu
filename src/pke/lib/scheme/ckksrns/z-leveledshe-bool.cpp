#include "scheme/ckksrns/z-leveledshe.h"
#include "scheme/ckksrns/z-user.h"
#include "encoding/z-encoding.h"

namespace lbcrypto {

CiphertextGroup UserZImpl::EvalBooleanAND(CiphertextGroup ct1, CiphertextGroup ct2) {
    if (ct1.size() != ct2.size()) {
        OPENFHE_THROW("CiphertextGroup size mismatch");
    }
    std::vector<Ciphertext<DCRTPoly>> resVec;
    for (size_t i = 0; i != ct1.size(); ++i) {
        auto res = z->EvalMultWithAdjust(ct1[i], ct2[i]);
        z->ModReduceInPlace(res);
        resVec.push_back(res);
    }
    return resVec;
}

CiphertextGroup UserZImpl::EvalBooleanOR(CiphertextGroup ct1, CiphertextGroup ct2) {
    if (ct1.size() != ct2.size()) {
        OPENFHE_THROW("CiphertextGroup size mismatch");
    }
    std::vector<Ciphertext<DCRTPoly>> resVec;
    for (size_t i = 0; i != ct1.size(); ++i) {
        auto addRes = z->EvalAddWithAdjust(ct1[i], ct2[i]);
        auto res    = z->EvalMultWithAdjust(ct1[i], ct2[i]);
        z->ModReduceInPlace(res);
        auto orRes = z->EvalSubWithAdjust(addRes, res);
        resVec.push_back(orRes);
    }
    return resVec;
}

CiphertextGroup UserZImpl::EvalBooleanXOR(CiphertextGroup ct1, CiphertextGroup ct2) {
    if (ct1.size() != ct2.size()) {
        OPENFHE_THROW("CiphertextGroup size mismatch");
    }
    std::vector<Ciphertext<DCRTPoly>> resVec;
    for (size_t i = 0; i != ct1.size(); ++i) {
        auto addRes = z->EvalAddWithAdjust(ct1[i], ct2[i]);
        auto res    = z->EvalMultWithAdjust(ct1[i], ct2[i]);
        res         = z->EvalMultScalar(res, 2);
        z->ModReduceInPlace(res);
        auto orRes = z->EvalSubWithAdjust(addRes, res);
        resVec.push_back(orRes);
    }
    return resVec;
}

CiphertextGroup UserZImpl::EvalBooleanNOT(CiphertextGroup ct1) {
    std::vector<Ciphertext<DCRTPoly>> resVec;
    for (size_t i = 0; i != ct1.size(); ++i) {
        auto res = z->EvalNegate(ct1[i]);
        z->EvalAddInPlaceInC(res, BigFixedPoint::one());
        resVec.push_back(res);
    }
    return resVec;
}

}  // namespace lbcrypto