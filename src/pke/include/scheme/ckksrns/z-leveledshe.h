#ifndef SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_LEVELEDSHE_H_
#define SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_LEVELEDSHE_H_

#include "math/hal/bigfixedpoint.h"
#include "cryptocontext.h"

namespace lbcrypto {

class LeveledZImpl;
using LeveledZ = std::shared_ptr<LeveledZImpl>;

// Here is short cut multiplication, ct1 * ct2 where one is in binary encoding
Ciphertext<DCRTPoly> zEvalMultShort(LeveledZ z, ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
// Here is full multiplication, ct1 * ct2 * t
Ciphertext<DCRTPoly> zEvalMultFull(LeveledZ z, ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2,
                                   Plaintext t);

// Helpers. Here ptxt will be ZEncoded
// If scalingFactor is not given, we use ct
Ciphertext<DCRTPoly> zEvalAdd(LeveledZ z, ConstCiphertext<DCRTPoly> ct, uint32_t ptxt);
// If scalingFactor is not given, we use ct
Ciphertext<DCRTPoly> zEvalMult(LeveledZ z, ConstCiphertext<DCRTPoly> ct, uint32_t ptxt,
                               BigFixedPoint scalingFactor = BigFixedPoint::zero());

class LeveledZImpl {
public:
    //
    // Generic methods involving Ciphertext and Plaintext
    //
    Ciphertext<DCRTPoly> EvalAdd(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt);
    Ciphertext<DCRTPoly> EvalSub(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt);
    Ciphertext<DCRTPoly> EvalAdd(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
    Ciphertext<DCRTPoly> EvalSub(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
    Ciphertext<DCRTPoly> EvalMult(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt);
    Ciphertext<DCRTPoly> EvalMult(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
    void EvalAddInPlace(Ciphertext<DCRTPoly> ct, Plaintext ptxt);
    void EvalAddInPlace(Ciphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2);
    void EvalSubInPlace(Ciphertext<DCRTPoly> ct, Plaintext ptxt);
    void EvalSubInPlace(Ciphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2);
    void EvalMultInPlace(Ciphertext<DCRTPoly> ct, Plaintext ptxt);

    // Special methods for MSB bootstrapping
    Ciphertext<DCRTPoly> EvalMultScalar(ConstCiphertext<DCRTPoly> ct, BigInteger scalar);
    void EvalMultScalarInPlace(Ciphertext<DCRTPoly> ct, BigInteger scalar);

    //
    // Level management
    //
private:
    // Should not allow user to call level reduce
    void LevelReduceInPlace(Ciphertext<DCRTPoly>& ciphertext, size_t levels = 1);

public:
    void ModReduceInPlace(Ciphertext<DCRTPoly>& ciphertext, size_t levels = 1);

    //
    // Cross level adjustment
    //
    Ciphertext<DCRTPoly> AdjustCiphertext(ConstCiphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ctTarget);

    // Automatic adjustment family
    Ciphertext<DCRTPoly> EvalMultWithAdjust(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
    Ciphertext<DCRTPoly> EvalAddWithAdjust(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
    Ciphertext<DCRTPoly> EvalSubWithAdjust(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);

    void EvalAddWithAdjustInPlace(Ciphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
    void EvalSubWithAdjustInPlace(Ciphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);

    //
    // Operations in Z
    //

    // Here is short cut multiplication, ct1 * ct2 where one is in binary encoding
    Ciphertext<DCRTPoly> EvalMultShortInZ(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
    // Here is full multiplication, ct1 * ct2 * t
    Ciphertext<DCRTPoly> EvalMultFullInZ(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);

    // Helpers. Here ptxt will be ZEncoded
    // If scalingFactor is not given, we use ct
    // TODO: support larger integer..
    Ciphertext<DCRTPoly> EvalAddInZ(ConstCiphertext<DCRTPoly> ct, uint32_t ptxt);
    // If scalingFactor is not given, we use ct
    Ciphertext<DCRTPoly> EvalMultInZ(ConstCiphertext<DCRTPoly> ct, uint32_t ptxt,
                                     BigFixedPoint scalingFactor = BigFixedPoint::zero());

    //
    // Operations in C
    //

    void EvalAddInPlaceInC(Ciphertext<DCRTPoly> ct, const BigComplex& ptxt);
    void EvalMultInPlaceInC(Ciphertext<DCRTPoly> ct, const BigComplex& ptxt,
                            BigFixedPoint scalingFactor = BigFixedPoint::zero());

    Ciphertext<DCRTPoly> EvalAddInC(ConstCiphertext<DCRTPoly> ct, const BigComplex& ptxt);

    Ciphertext<DCRTPoly> EvalMultInC(ConstCiphertext<DCRTPoly> ct, const BigComplex& ptxt,
                                     BigFixedPoint scalingFactor = BigFixedPoint::zero());

    Ciphertext<DCRTPoly> EvalConjugateInC(ConstCiphertext<DCRTPoly> ct);

private:
    // value, scalingFactor, modulus
    using BCInCPlaintextKey = std::tuple<BigComplex, BigFixedPoint, BigInteger>;
    struct BCInCPlaintextKeyCompare {
        bool operator()(const BCInCPlaintextKey& a, const BCInCPlaintextKey& b) const;
    };

    std::map<BCInCPlaintextKey, Plaintext, BCInCPlaintextKeyCompare> m_bcInCPlaintextCache;

    Plaintext GetBCInCPlaintext(const BigComplex& value, const BigFixedPoint& scalingFactor,
                                const std::shared_ptr<typename DCRTPoly::Params>& elementParams);
};

using LeveledZ = std::shared_ptr<LeveledZImpl>;

}  // namespace lbcrypto

#endif  //SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_LEVELEDSHE_H_