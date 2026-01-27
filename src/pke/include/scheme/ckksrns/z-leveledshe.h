#ifndef SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_LEVELEDSHE_H_
#define SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_LEVELEDSHE_H_

#include "math/hal/bigfixedpoint.h"
#include "cryptocontext.h"

namespace lbcrypto {

class LeveledZImpl;
using LeveledZ = std::shared_ptr<LeveledZImpl>;

class CiphertextGroup {
public:
    CiphertextGroup() = default;
    CiphertextGroup(Ciphertext<DCRTPoly> ct) : parts({ct}) {}
    CiphertextGroup(std::vector<Ciphertext<DCRTPoly>> cts) : parts(cts) {}

    std::vector<Ciphertext<DCRTPoly>> getParts() const {
        return parts;
    }

    Ciphertext<DCRTPoly>& operator[](size_t idx) {
        return parts[idx];
    }

    size_t size() const {
        return parts.size();
    }

    operator Ciphertext<DCRTPoly>() const {
        if (parts.size() != 1) {
            OPENFHE_THROW("Cannot convert CiphertextGroup with multiple parts to single Ciphertext");
        }
        return parts[0];
    }

    using MapFunc = std::function<Ciphertext<DCRTPoly>(ConstCiphertext<DCRTPoly>&)>;
    CiphertextGroup map(MapFunc func) const {
        std::vector<Ciphertext<DCRTPoly>> result;
        for (const auto& part : parts) {
            result.push_back(func(part));
        }
        return CiphertextGroup(result);
    }

    using MapWideFunc = std::function<CiphertextGroup(ConstCiphertext<DCRTPoly>&)>;
    CiphertextGroup mapWide(MapWideFunc func) const {
        std::vector<Ciphertext<DCRTPoly>> result;
        for (const auto& part : parts) {
            auto newParts = func(part);
            for (auto& newPart : newParts.getParts()) {
                result.push_back(newPart);
            }
        }
        return CiphertextGroup(result);
    }

private:
    std::vector<Ciphertext<DCRTPoly>> parts;
};

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
    Ciphertext<DCRTPoly> AdjustCiphertextToLevel(ConstCiphertext<DCRTPoly> ciphertext, size_t level);
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

    //
    // Boolean Mode Operations
    //

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