#ifndef SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_USER_H_
#define SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_USER_H_

#include "scheme/ckksrns/z-leveledshe.h"

namespace lbcrypto {

// User facing API
// LeveledSHE API are too internal...
class UserZImpl {
public:
    UserZImpl(LeveledZ z) : z(z) {}

    //
    // Operations in Z
    //
    Ciphertext<DCRTPoly> EvalAddInZ(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
    Ciphertext<DCRTPoly> EvalSubInZ(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
    Ciphertext<DCRTPoly> EvalNegateInZ(ConstCiphertext<DCRTPoly> ct1);

    // Here is short cut multiplication, ct1 * ct2 where one is in binary encoding
    Ciphertext<DCRTPoly> EvalMultShortInZ(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
    // Here is full multiplication, ct1 * ct2 * t
    Ciphertext<DCRTPoly> EvalMultFullInZ(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);

    // Helpers. Here ptxt will be ZEncoded
    // TODO: add vector<BigInteger> version
    Ciphertext<DCRTPoly> EvalAddInZ(ConstCiphertext<DCRTPoly> ct, BigInteger ptxt);
    Ciphertext<DCRTPoly> EvalMultInZ(ConstCiphertext<DCRTPoly> ct, BigInteger ptxt);

    // Conversion between [m]_t / t and [m]_t
    Ciphertext<DCRTPoly> EvalMultTInZ(ConstCiphertext<DCRTPoly> ct);
    Ciphertext<DCRTPoly> EvalMultTInvInZ(ConstCiphertext<DCRTPoly> ct);

    //
    // Boolean Mode Operations
    //

    // TODO: add ct-pt versions
    CiphertextGroup EvalBooleanAND(CiphertextGroup ct1, CiphertextGroup ct2);
    CiphertextGroup EvalBooleanOR(CiphertextGroup ct1, CiphertextGroup ct2);
    CiphertextGroup EvalBooleanXOR(CiphertextGroup ct1, CiphertextGroup ct2);
    CiphertextGroup EvalBooleanNOT(CiphertextGroup ct1);

    // Here left/right is defined in a big-endian manner
    CiphertextGroup EvalBooleanShiftLeft(CiphertextGroup ct1, uint64_t offset);
    CiphertextGroup EvalBooleanShiftRight(CiphertextGroup ct1, uint64_t offset);
    CiphertextGroup EvalBooleanRotateLeft(CiphertextGroup ct1, uint64_t offset);
    CiphertextGroup EvalBooleanRotateRight(CiphertextGroup ct1, uint64_t offset);

    Ciphertext<DCRTPoly> EvalSignExtract(CiphertextGroup ct1);

private:
    LeveledZ z;
};

using UserZ = std::shared_ptr<UserZImpl>;

}  // namespace lbcrypto

#endif  //SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_USER_H_