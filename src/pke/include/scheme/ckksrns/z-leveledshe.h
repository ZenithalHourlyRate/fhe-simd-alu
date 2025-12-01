#ifndef SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_LEVELEDSHE_H_
#define SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_LEVELEDSHE_H_

#include "math/hal/bigfixedpoint.h"
#include "cryptocontext.h"

namespace lbcrypto {

// These methods assume that the Plaintext is ZEncoded
Ciphertext<DCRTPoly> zEvalAdd(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt);
Ciphertext<DCRTPoly> zEvalMult(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt);
Ciphertext<DCRTPoly> zEvalMult(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
Ciphertext<DCRTPoly> zEvalMultScalar(ConstCiphertext<DCRTPoly> ct, BigInteger scalar);
void zModReduceInPlace(Ciphertext<DCRTPoly>& ciphertext, size_t levels = 1);

// Helpers. Here ptxt will be ZEncoded
Ciphertext<DCRTPoly> zEvalAdd(ConstCiphertext<DCRTPoly> ct, uint32_t ptxt);
Ciphertext<DCRTPoly> zEvalMult(ConstCiphertext<DCRTPoly> ct, uint32_t ptxt);

// Helpers. Here ptxt will be CEncoded
Ciphertext<DCRTPoly> cEvalAdd(ConstCiphertext<DCRTPoly> ct, uint32_t ptxt);
Ciphertext<DCRTPoly> cEvalMult(ConstCiphertext<DCRTPoly> ct, uint32_t ptxt);

// Helpers. Here ptxt will be REncoded
Ciphertext<DCRTPoly> rEvalAdd(ConstCiphertext<DCRTPoly> ct, uint32_t ptxt);
Ciphertext<DCRTPoly> rEvalMult(ConstCiphertext<DCRTPoly> ct, uint32_t ptxt);

// These methods assume that the Plaintext is REncoded

}  // namespace lbcrypto

#endif  //SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_LEVELEDSHE_H_