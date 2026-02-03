#include "scheme/ckksrns/z-user-advanced.h"

namespace lbcrypto {

void UserZImpl::Setup(const PrivateKey<DCRTPoly> privateKey, uint32_t zN, uint32_t zSlot,
                      std::vector<uint32_t> leftShifts, std::vector<uint32_t> rightShifts, bool enableSignExtract) {
    auto cc   = privateKey->GetCryptoContext();
    auto algo = cc->GetScheme();

    std::vector<int32_t> indices;
    if (enableSignExtract) {
        // For full packing...
        indices.push_back(static_cast<int32_t>(zN / 2 - 1));
    }
    auto evalKeys = algo->EvalAtIndexKeyGen(privateKey, indices);
    cc->InsertEvalAutomorphismKey(evalKeys, privateKey->GetKeyTag());

    // TODO: add left/right shifts keys
}

CiphertextGroup UserZAdvancedImpl::EvalLessThan(Ciphertext<DCRTPoly> ct1, Ciphertext<DCRTPoly> ct2) {
    VERIFY_ARITH_ENCODING(ct1);
    VERIFY_ARITH_ENCODING(ct2);
    // Compute ct1 - ct2
    auto ctDiff     = userZ->EvalSubInZ(ct1, ct2);
    auto ctDiffBool = fheZ->EvalArithToBoolean(ctDiff);

    // Sign extract
    auto sign = userZ->EvalSignExtract(ctDiffBool);
    return sign;
}

//CiphertextGroup UserZAdvancedImpl::EvalEqualTo(Ciphertext<DCRTPoly> ct1, Ciphertext<DCRTPoly> ct2) {
//    // TODO: optimize it using reduce operation....
//    VERIFY_ARITH_ENCODING(ct1);
//    VERIFY_ARITH_ENCODING(ct2);
//    //// Compute ct1 - ct2
//    //auto ctDiff = userZ->EvalSubInZ(ct1, ct2);
//    //// Compute ct2 - ct1
//    //auto ctDiff2    = userZ->EvalSubInZ(ct2, ct1);
//
//    //auto ctDiffBool = fheZ->EvalArithToBoolean(CiphertextGroup({ctDiff, ctDiff2}));
//
//    //// Sign extract
//    //auto sign = userZ->EvalSignExtract(ctDiffBool.append(ctDiff2Bool));
//
//    //// NOT sign
//    //auto signNeg = userZ->EvalBooleanNOT(sign);
//    //return signNeg;
//}

}  // namespace lbcrypto