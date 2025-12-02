#ifndef SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_LEVELEDSHE_H_
#define SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_LEVELEDSHE_H_

#include "math/hal/bigfixedpoint.h"
#include "cryptocontext.h"

namespace lbcrypto {

// Generic methods
Ciphertext<DCRTPoly> gEvalAdd(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt);
Ciphertext<DCRTPoly> gEvalSub(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt);
Ciphertext<DCRTPoly> gEvalAdd(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
Ciphertext<DCRTPoly> gEvalSub(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
Ciphertext<DCRTPoly> gEvalMult(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt);
Ciphertext<DCRTPoly> gEvalMult(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
Ciphertext<DCRTPoly> gEvalMultScalar(ConstCiphertext<DCRTPoly> ct, BigInteger scalar);

void gEvalAddInPlace(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt);
void gEvalAddInPlace(Ciphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2);
void gEvalSubInPlace(Ciphertext<DCRTPoly> ct, Plaintext ptxt);
void gEvalSubInPlace(Ciphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2);
void gEvalMultScalarInPlace(Ciphertext<DCRTPoly> ct, BigInteger scalar);
void gLevelReduceInPlace(Ciphertext<DCRTPoly>& ciphertext, size_t levels = 1);
void gModReduceInPlace(Ciphertext<DCRTPoly>& ciphertext, size_t levels = 1);

Ciphertext<DCRTPoly> gAdjustCiphertext(ConstCiphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ctTarget);

// Automatic adjustment family
Ciphertext<DCRTPoly> gEvalMultWithAdjust(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
Ciphertext<DCRTPoly> gEvalAddWithAdjust(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
Ciphertext<DCRTPoly> gEvalSubWithAdjust(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);

void gEvalAddWithAdjustInPlace(Ciphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
void gEvalSubWithAdjustInPlace(Ciphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);

// Here is short cut multiplication, ct1 * ct2 where one is in binary encoding
Ciphertext<DCRTPoly> zEvalMultShort(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
// Here is full multiplication, ct1 * ct2 * t
Ciphertext<DCRTPoly> zEvalMultFull(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2, Plaintext t);

// Helpers. Here ptxt will be ZEncoded
// If scalingFactor is not given, we use ct
Ciphertext<DCRTPoly> zEvalAdd(ConstCiphertext<DCRTPoly> ct, uint32_t ptxt);
// If scalingFactor is not given, we use ct
Ciphertext<DCRTPoly> zEvalMult(ConstCiphertext<DCRTPoly> ct, uint32_t ptxt,
                               BigFixedPoint scalingFactor = BigFixedPoint::zero());

// Helpers. Here ptxt will be CEncoded
Ciphertext<DCRTPoly> cEvalAdd(ConstCiphertext<DCRTPoly> ct, BigComplex ptxt);
Ciphertext<DCRTPoly> cEvalMult(ConstCiphertext<DCRTPoly> ct, BigComplex ptxt,
                               BigFixedPoint scalingFactor = BigFixedPoint::zero());

void cEvalAddInPlace(Ciphertext<DCRTPoly> ct, BigComplex ptxt);

// Helpers. Here ptxt will be REncoded
// Ciphertext<DCRTPoly> rEvalAdd(ConstCiphertext<DCRTPoly> ct, RPolynomial ptxt);
// Ciphertext<DCRTPoly> rEvalMult(ConstCiphertext<DCRTPoly> ct, RPolynomial ptxt);

// These methods assume that the Plaintext is REncoded

}  // namespace lbcrypto

#endif  //SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_LEVELEDSHE_H_