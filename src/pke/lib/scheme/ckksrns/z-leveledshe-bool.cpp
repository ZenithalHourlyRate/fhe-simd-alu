#include "scheme/ckksrns/z-leveledshe.h"
#include "encoding/z-encoding.h"

namespace lbcrypto {

CiphertextGroup LeveledZImpl::EvalBooleanAND(CiphertextGroup ct1, CiphertextGroup ct2) {
    if (ct1.size() != ct2.size()) {
        OPENFHE_THROW("CiphertextGroup size mismatch");
    }
    std::vector<Ciphertext<DCRTPoly>> resVec;
    for (size_t i = 0; i != ct1.size(); ++i) {
        auto res = EvalMultWithAdjust(ct1[i], ct2[i]);
        ModReduceInPlace(res);
        resVec.push_back(res);
    }
    return resVec;
}

CiphertextGroup LeveledZImpl::EvalBooleanOR(CiphertextGroup ct1, CiphertextGroup ct2) {
    if (ct1.size() != ct2.size()) {
        OPENFHE_THROW("CiphertextGroup size mismatch");
    }
    std::vector<Ciphertext<DCRTPoly>> resVec;
    for (size_t i = 0; i != ct1.size(); ++i) {
        auto addRes = EvalAddWithAdjust(ct1[i], ct2[i]);
        auto res    = EvalMultWithAdjust(ct1[i], ct2[i]);
        ModReduceInPlace(res);
        auto orRes = EvalSubWithAdjust(addRes, res);
        resVec.push_back(orRes);
    }
    return resVec;
}

CiphertextGroup LeveledZImpl::EvalBooleanXOR(CiphertextGroup ct1, CiphertextGroup ct2) {
    if (ct1.size() != ct2.size()) {
        OPENFHE_THROW("CiphertextGroup size mismatch");
    }
    std::vector<Ciphertext<DCRTPoly>> resVec;
    for (size_t i = 0; i != ct1.size(); ++i) {
        auto addRes = EvalAddWithAdjust(ct1[i], ct2[i]);
        auto res    = EvalMultWithAdjust(ct1[i], ct2[i]);
        res         = EvalMultScalar(res, 2);
        ModReduceInPlace(res);
        auto orRes = EvalSubWithAdjust(addRes, res);
        resVec.push_back(orRes);
    }
    return resVec;
}

CiphertextGroup LeveledZImpl::EvalBooleanNOT(CiphertextGroup ct1) {
    std::vector<Ciphertext<DCRTPoly>> resVec;
    for (size_t i = 0; i != ct1.size(); ++i) {
        auto res = EvalNegate(ct1[i]);
        EvalAddInPlaceInC(res, BigFixedPoint::one());
        resVec.push_back(res);
    }
    return resVec;
}

}  // namespace lbcrypto