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
    Ciphertext<DCRTPoly> EvalNegate(ConstCiphertext<DCRTPoly> ct1);
    Ciphertext<DCRTPoly> EvalMult(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt);
    Ciphertext<DCRTPoly> EvalMult(ConstCiphertext<DCRTPoly> ct1, ConstCiphertext<DCRTPoly> ct2);
    Ciphertext<DCRTPoly> EvalSquare(ConstCiphertext<DCRTPoly> ct1);
    void EvalAddInPlace(Ciphertext<DCRTPoly> ct, Plaintext ptxt);
    void EvalAddInPlace(Ciphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2);
    void EvalSubInPlace(Ciphertext<DCRTPoly> ct, Plaintext ptxt);
    void EvalSubInPlace(Ciphertext<DCRTPoly> ct, ConstCiphertext<DCRTPoly> ct2);
    void EvalNegateInPlace(Ciphertext<DCRTPoly> ct);
    void EvalMultInPlace(Ciphertext<DCRTPoly> ct, Plaintext ptxt);
    // Buggy
    //void EvalSquareInPlace(Ciphertext<DCRTPoly> ct);

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
    // This wont change the zSlot scale
    Ciphertext<DCRTPoly> EvalAddInZ(ConstCiphertext<DCRTPoly> ct, BigInteger ptxt);
    // This wont change the zSlot scale
    Ciphertext<DCRTPoly> EvalMultInZ(ConstCiphertext<DCRTPoly> ct, BigInteger ptxt);

    // Conversion between [m]_t / t and [m]_t
    Ciphertext<DCRTPoly> EvalMultTInZ(ConstCiphertext<DCRTPoly> ct);
    Ciphertext<DCRTPoly> EvalMultTInvInZ(ConstCiphertext<DCRTPoly> ct);

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

    // zN, scaling factor, modulus
    using tPlaintextKey = std::tuple<BigInteger, BigFixedPoint, BigInteger>;
    struct tPlaintextKeyCompare {
        bool operator()(const tPlaintextKey& a, const tPlaintextKey& b) const;
    };

    std::map<tPlaintextKey, Plaintext, tPlaintextKeyCompare> m_tPlaintextCache;
    std::map<tPlaintextKey, Plaintext, tPlaintextKeyCompare> m_tInvPlaintextCache;

    Plaintext GetTPlaintext(uint32_t zN, const BigFixedPoint& scalingFactor,
                            const std::shared_ptr<typename DCRTPoly::Params>& elementParams);
    Plaintext GetTInvPlaintext(uint32_t zN, const BigFixedPoint& scalingFactor,
                               const std::shared_ptr<typename DCRTPoly::Params>& elementParams);
};

using LeveledZ = std::shared_ptr<LeveledZImpl>;

}  // namespace lbcrypto

#endif  //SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_LEVELEDSHE_H_