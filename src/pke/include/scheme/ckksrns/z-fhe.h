#ifndef SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_FHE_H_
#define SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_FHE_H_

#include "z-leveledshe.h"
#include "z-advancedshe.h"
#include "z-fhe-utils.h"

namespace lbcrypto {

class ZBootstrapPlaintextCacheImpl {
public:
    ZBootstrapPlaintextCacheImpl()  = default;
    ~ZBootstrapPlaintextCacheImpl() = default;

    ZBootstrapPlaintextCacheImpl(const BigCVector& val) : m_value(val) {}

    Plaintext GetPlaintext(const BigFixedPoint& scalingFactor,
                           const std::shared_ptr<typename DCRTPoly::Params>& elementParams);

    BigCVector& getValue() {
        return m_value;
    }

private:
    BigCVector m_value;

    // scalingFactor, modulus
    using Key = std::tuple<BigFixedPoint, BigInteger>;
    std::map<Key, Plaintext> m_cache;
};

using ZBootstrapPlaintextCache = std::shared_ptr<ZBootstrapPlaintextCacheImpl>;

class ZBootstrapPrecom {
public:
    ZBootstrapPrecom() = default;

    virtual ~ZBootstrapPrecom() = default;

    ZBootstrapPrecom(const ZBootstrapPrecom& rhs) = default;

    ZBootstrapPrecom(ZBootstrapPrecom&& rhs) noexcept = default;

    // level budget for homomorphic encoding, number of layers to collapse in one level,
    // number of layers remaining to be collapsed in one level to have exactly the number
    // of levels specified in the level budget, the number of rotations in one level,
    // the baby step and giant step in the baby-step giant-step strategy, the number of
    // rotations in the remaining level, the baby step and giant step in the baby-step
    // giant-step strategy for the remaining level
    struct ckks_boot_params m_paramsEnc;

    // level budget for homomorphic decoding, number of layers to collapse in one level,
    // number of layers remaining to be collapsed in one level to have exactly the number
    // of levels specified in the level budget, the number of rotations in one level,
    // the baby step and giant step in the baby-step giant-step strategy, the number of
    // rotations in the remaining level, the baby step and giant step in the baby-step
    // giant-step strategy for the remaining level
    struct ckks_boot_params m_paramsDec;

    // number of complex slots for which the bootstrapping is performed
    // In Z we have two kinds of slots: complex slots and Z slots
    uint32_t m_cSlots;

    // Linear map U0; used in decoding
    std::vector<ZBootstrapPlaintextCache> m_U0Pre;

    // Conj(U0^T); used in encoding
    std::vector<ZBootstrapPlaintextCache> m_U0hatTPre;

    // coefficients corresponding to U0; used in decoding
    std::vector<std::vector<ZBootstrapPlaintextCache>> m_U0PreFFT;

    // coefficients corresponding to conj(U0^T); used in encoding
    std::vector<std::vector<ZBootstrapPlaintextCache>> m_U0hatTPreFFT;

    //
    // coefficients for ZU and ZUInverse
    //
    // digonal, 1st digonal, ..., zN/2-1th diagonal, -1th diagonal, ..., -zN/2+1 th diagonal
    std::vector<ZBootstrapPlaintextCache> m_ZUPre;
    std::vector<ZBootstrapPlaintextCache> m_ZUInversePre;

    std::vector<ZBootstrapPlaintextCache> m_ZU0Pre;
    std::vector<ZBootstrapPlaintextCache> m_ZU1Pre;
    std::vector<ZBootstrapPlaintextCache> m_ZUInverse0Pre;
    std::vector<ZBootstrapPlaintextCache> m_ZUInverse1Pre;
};

class FHEZImpl {
public:
    FHEZImpl(LeveledZ z, AdvancedZ advZ) : z(z), advZ(advZ) {}

    void EvalBootstrapSetup(const CryptoContextImpl<DCRTPoly>& cc, uint32_t numCSlots,
                            std::vector<uint32_t> levelBudget, std::vector<uint32_t> dim1 = {0, 0});

    Ciphertext<DCRTPoly> EvalLinearTransform(std::vector<ZBootstrapPlaintextCache>& A,
                                             ConstCiphertext<DCRTPoly>& ct) const;

    Ciphertext<DCRTPoly> EvalZLinearTransform(std::vector<ZBootstrapPlaintextCache>& A,
                                              ConstCiphertext<DCRTPoly>& ct) const;

    Ciphertext<DCRTPoly> EvalCoeffsToSlots(const std::vector<std::vector<ZBootstrapPlaintextCache>>& A,
                                           ConstCiphertext<DCRTPoly>& ctxt, uint32_t cSlots) const;

    Ciphertext<DCRTPoly> EvalSlotsToCoeffs(const std::vector<std::vector<ZBootstrapPlaintextCache>>& A,
                                           ConstCiphertext<DCRTPoly>& ctxt, uint32_t cSlots) const;

    Ciphertext<DCRTPoly> EvalTruncate(ConstCiphertext<DCRTPoly>& ctxt) const;

    Ciphertext<DCRTPoly> EvalModRaise(ConstCiphertext<DCRTPoly>& ctxt) const;

    Ciphertext<DCRTPoly> EvalArithToArithHigh(ConstCiphertext<DCRTPoly>& ctxt, uint32_t cSlots) const;
    // Temporary:
public:
    ZBootstrapPrecom& GetBootPrecom(uint32_t slots) const {
        auto pair = m_bootPrecomMap.find(slots);
        if (pair != m_bootPrecomMap.end())
            return *(pair->second);
        OPENFHE_THROW("Precomputations for " + std::to_string(slots) + " slots not found.");
    }

private:
    //------------------------------------------------------------------------------
    // Precomputations for ZCoeffsToSlots and SlotsToZCoeffs
    //------------------------------------------------------------------------------
    std::vector<ZBootstrapPlaintextCache> EvalZLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                         const BigCMatrix& A, uint32_t zSlots) const;

    std::vector<ZBootstrapPlaintextCache> EvalZLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                         const BigCMatrix& A, const BigCMatrix& B,
                                                                         uint32_t zSlots) const;

    //------------------------------------------------------------------------------
    // Precomputations for CoeffsToSlots and SlotsToCoeffs
    //------------------------------------------------------------------------------

    std::vector<ZBootstrapPlaintextCache> EvalLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                        const BigCMatrix& A) const;

    std::vector<ZBootstrapPlaintextCache> EvalLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                        const BigCMatrix& A, const BigCMatrix& B,
                                                                        uint32_t orientation = 0) const;

    std::vector<std::vector<ZBootstrapPlaintextCache>> EvalCoeffsToSlotsPrecompute(
        const CryptoContextImpl<DCRTPoly>& cc, const BigCVector& pows, const std::vector<uint32_t>& rotGroup,
        bool flag_i) const;

    std::vector<std::vector<ZBootstrapPlaintextCache>> EvalSlotsToCoeffsPrecompute(
        const CryptoContextImpl<DCRTPoly>& cc, const BigCVector& pows, const std::vector<uint32_t>& rotGroup,
        bool flag_i) const;

private:
    // corresponds to probability of less than 2^{-128}
    static constexpr uint32_t K_SPARSE_ENCAPSULATED = 16;

private:
    LeveledZ z;
    AdvancedZ advZ;

    // key tuple is dim1, levelBudgetEnc, levelBudgetDec
    // key tuple is cSlots, zSlots,
    std::map<uint32_t, std::shared_ptr<ZBootstrapPrecom>> m_bootPrecomMap;
};

using FHEZ = std::shared_ptr<FHEZImpl>;

}  // namespace lbcrypto

#endif  //SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_FHE_H_