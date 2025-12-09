#include "encoding/z-encode.h"
#include "math/z-constants.h"
#include "scheme/ckksrns/z-fhe.h"
#include "math/dftransform-bigcomplex.h"

namespace lbcrypto {

//------------------------------------------------------------------------------
// Bootstrap Wrapper
//------------------------------------------------------------------------------

void FHEZImpl::EvalBootstrapSetup(const CryptoContextImpl<DCRTPoly>& cc, uint32_t numCSlots,
                                  std::vector<uint32_t> levelBudget, std::vector<uint32_t> dim1) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

    uint32_t M     = cc.GetCyclotomicOrder();
    uint32_t slots = (numCSlots == 0) ? M / 4 : numCSlots;

    m_bootPrecomMap[slots] = std::make_shared<ZBootstrapPrecom>();

    auto& precom     = m_bootPrecomMap[slots];
    precom->m_cSlots = slots;

    // even for the case of a single slot we need one level for rescaling
    uint32_t logSlots = (slots < 3) ? 1 : std::log2(slots);

    // Perform some checks on the level budget and compute parameters
    uint32_t newBudget0 = levelBudget[0];
    if (newBudget0 > logSlots) {
        std::cerr << "\nWarning, the level budget for encoding is too large. Setting it to " << logSlots << std::endl;
        newBudget0 = logSlots;
    }
    if (newBudget0 < 1) {
        std::cerr << "\nWarning, the level budget for encoding can not be zero. Setting it to 1" << std::endl;
        newBudget0 = 1;
    }
    uint32_t newBudget1 = levelBudget[1];
    if (newBudget1 > logSlots) {
        std::cerr << "\nWarning, the level budget for decoding is too large. Setting it to " << logSlots << std::endl;
        newBudget1 = logSlots;
    }
    if (newBudget1 < 1) {
        std::cerr << "\nWarning, the level budget for decoding can not be zero. Setting it to 1" << std::endl;
        newBudget1 = 1;
    }

    precom->m_paramsEnc = GetCollapsedFFTParams(slots, newBudget0, dim1[0]);
    precom->m_paramsDec = GetCollapsedFFTParams(slots, newBudget1, dim1[1]);

    uint32_t m     = 4 * slots;
    uint32_t mmask = m - 1;  // assumes m is power of 2
    bool isSparse  = (M != m);

    // computes indices for all primitive roots of unity
    std::vector<uint32_t> rotGroup(slots);
    uint32_t fivePows = 1;
    for (uint32_t i = 0; i < slots; ++i) {
        rotGroup[i] = fivePows;
        fivePows *= 5;
        fivePows &= mmask;
    }

    // computes all powers of a primitive root of unity exp(2 * M_PI/m)
    std::vector<BigComplex> ksiPows(m + 1);
    ksiPows[0] = BigComplex(BigFixedPoint::one(), BigFixedPoint::zero());
    ksiPows[1] = R_ROOT_MAP.at(m);
    for (uint32_t j = 2; j < m; ++j) {
        ksiPows[j] = ksiPows[j - 1] * ksiPows[1];
    }
    ksiPows[m] = ksiPows[0];

    //uint32_t approxModDepth = GetModDepthInternal(cryptoParams->GetSecretKeyDist());

    //uint32_t depthBT = approxModDepth + precom->m_paramsEnc.lvlb + precom->m_paramsDec.lvlb;

    // compute # of levels to remain when encoding the coefficients
    // for FLEXIBLEAUTOEXT we do not need extra modulus in auxiliary plaintexts
    //uint32_t L0 = cryptoParams->GetElementParams()->GetParams().size();

    //uint32_t lEnc = L0 - (precom->m_paramsEnc.lvlb + 1);
    //uint32_t lDec = L0 - depthBT;

    bool isLTBootstrap = (precom->m_paramsEnc.lvlb == 1) && (precom->m_paramsDec.lvlb == 1);

    auto I = BigComplex(BigFixedPoint::zero(), BigFixedPoint::one());

    if (isLTBootstrap) {
        if (isSparse) {
            BigCMatrix U0(slots, BigCVector(slots));
            BigCMatrix U0hatT(slots, BigCVector(slots));
            BigCMatrix U1(slots, BigCVector(slots));
            BigCMatrix U1hatT(slots, BigCVector(slots));
            for (uint32_t i = 0; i < slots; ++i) {
                for (uint32_t j = 0; j < slots; ++j) {
                    U0[i][j]     = ksiPows[(j * rotGroup[i]) & mmask];
                    U0hatT[j][i] = U0[i][j].conj();
                    U1[i][j]     = I * U0[i][j];
                    U1hatT[j][i] = U1[i][j].conj();
                }
            }
            precom->m_U0Pre     = EvalLinearTransformPrecompute(cc, U0, U1, 1);
            precom->m_U0hatTPre = EvalLinearTransformPrecompute(cc, U0hatT, U1hatT, 0);
        }
        else {
            BigCMatrix U0(slots, BigCVector(slots));
            BigCMatrix U0hatT(slots, BigCVector(slots));
            for (uint32_t i = 0; i < slots; ++i) {
                for (uint32_t j = 0; j < slots; ++j) {
                    U0[i][j]     = ksiPows[(j * rotGroup[i]) & mmask];
                    U0hatT[j][i] = U0[i][j].conj();
                }
            }
            precom->m_U0Pre     = EvalLinearTransformPrecompute(cc, U0);
            precom->m_U0hatTPre = EvalLinearTransformPrecompute(cc, U0hatT);
        }
    }
    else {
        precom->m_U0PreFFT     = EvalSlotsToCoeffsPrecompute(cc, ksiPows, rotGroup, false);
        precom->m_U0hatTPreFFT = EvalCoeffsToSlotsPrecompute(cc, ksiPows, rotGroup, false);
    }
}

//------------------------------------------------------------------------------
// Precomputations for CoeffsToSlots and SlotsToCoeffs
//------------------------------------------------------------------------------

std::shared_ptr<typename DCRTPoly::Params> FHEZImpl::GetParamsQP(
    uint32_t m, const std::vector<std::shared_ptr<ILNativeParams>>& paramsQ,
    const std::vector<std::shared_ptr<ILNativeParams>>& paramsP) const {
    uint32_t sizeQ = paramsQ.size();
    uint32_t sizeP = paramsP.size();
    std::vector<NativeInteger> moduli(sizeQ + sizeP);
    std::vector<NativeInteger> roots(sizeQ + sizeP);
    for (uint32_t i = 0; i < sizeQ; ++i) {
        moduli[i] = paramsQ[i]->GetModulus();
        roots[i]  = paramsQ[i]->GetRootOfUnity();
    }
    for (uint32_t i = 0; i < sizeP; ++i) {
        moduli[sizeQ + i] = paramsP[i]->GetModulus();
        roots[sizeQ + i]  = paramsP[i]->GetRootOfUnity();
    }
    return std::make_shared<ILDCRTParams<DCRTPoly::Integer>>(m, moduli, roots);
}

std::vector<ZBootstrapPlaintextCache> FHEZImpl::EvalLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                              const BigCMatrix& A) const {
    const int32_t slots = A.size();
    if (slots != static_cast<int32_t>(A[0].size()))
        OPENFHE_THROW("The matrix passed to EvalLTPrecompute is not square");

    auto g = GetBootPrecom(slots).m_paramsEnc.g;

    const int32_t step = (g == 0) ? std::ceil(std::sqrt(slots)) : g;

    std::vector<ZBootstrapPlaintextCache> result(slots);
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(slots))
#endif
    for (int32_t ji = 0; ji < slots; ++ji) {
        auto diag = ExtractShiftedDiagonal(A, ji);
        //for (auto& d : diag)
        //    d *= scale;
        result[ji] = Rotate(diag, -step * (ji / step));
    }
    return result;
}

std::vector<ZBootstrapPlaintextCache> FHEZImpl::EvalLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                              const BigCMatrix& A, const BigCMatrix& B,
                                                                              uint32_t orientation) const {
    const int32_t slots = static_cast<int32_t>(A.size());

    auto g = GetBootPrecom(slots).m_paramsEnc.g;

    const int32_t step = (g == 0) ? std::ceil(std::sqrt(slots)) : g;

    std::vector<ZBootstrapPlaintextCache> result(slots);

    if (orientation == 0) {
        // vertical concatenation - used during homomorphic encoding
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(slots))
#endif
        for (int32_t ji = 0; ji < slots; ++ji) {
            auto vecA = ExtractShiftedDiagonal(A, ji);
            auto vecB = ExtractShiftedDiagonal(B, ji);
            vecA.insert(vecA.end(), vecB.begin(), vecB.end());
            result[ji] = Rotate(vecA, -step * (ji / step));
        }
    }
    else {
        // horizontal concatenation - used during homomorphic decoding
        BigCMatrix newA(slots);

        //  A and B are concatenated horizontally
        for (int32_t i = 0; i < slots; ++i) {
            newA[i].reserve(A[i].size() + B[i].size());
            newA[i].insert(newA[i].end(), A[i].begin(), A[i].end());
            newA[i].insert(newA[i].end(), B[i].begin(), B[i].end());
        }

#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(slots))
#endif
        for (int32_t ji = 0; ji < slots; ++ji) {
            // shifted diagonal is computed for rectangular map newA of dimension
            // slots x 2*slots
            auto vec   = ExtractShiftedDiagonal(newA, ji);
            result[ji] = Rotate(vec, -step * (ji / step));
        }
    }

    return result;
}

std::vector<std::vector<ZBootstrapPlaintextCache>> FHEZImpl::EvalCoeffsToSlotsPrecompute(
    const CryptoContextImpl<DCRTPoly>& cc, const BigCVector& A, const std::vector<uint32_t>& rotGroup,
    bool flag_i) const {
    const uint32_t slots = rotGroup.size();

    const auto& p = GetBootPrecom(slots).m_paramsEnc;

    // result is the rotated plaintext version of the coefficients
    std::vector<std::vector<ZBootstrapPlaintextCache>> result(p.lvlb,
                                                              std::vector<ZBootstrapPlaintextCache>(p.numRotations));

    int32_t stop    = -1;
    int32_t flagRem = 0;
    if (p.remCollapse != 0) {
        stop    = 0;
        flagRem = 1;

        // remainder corresponds to index 0 in encoding and to last index in decoding
        result[0].resize(p.numRotationsRem);
    }

    auto M = cc.GetCyclotomicOrder();

    if (uint32_t M4 = M / 4; slots == M4) {
        //------------------------------------------------------------------------------
        // fully-packed mode
        //------------------------------------------------------------------------------

        auto coeff = CoeffEncodingCollapse(A, rotGroup, p.lvlb, flag_i);

        for (int32_t s = -1 + p.lvlb; s > stop; --s) {
            const int32_t rotScale = (1 << ((s - flagRem) * p.layersCollapse + p.remCollapse)) * p.g;
            const uint32_t limit   = p.b * p.g;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotations) {
                    // Maybe we need it one day???
                    if ((flagRem == 0) && (s == stop + 1)) {
                        // do the scaling only at the last set of coefficients
                        //for (auto& c : coeff[s][ij])
                        //    c *= scale;
                    }

                    auto rot = Rotate(coeff[s][ij], ReduceRotation(-rotScale * (ij / p.g), slots));

                    result[s][ij] = rot;
                }
            }
        }

        if (flagRem == 1) {
            const uint32_t limit = p.bRem * p.gRem;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotationsRem) {
                    //for (auto& c : coeff[stop][ij])
                    //    c *= scale;

                    auto rot = Rotate(coeff[stop][ij], ReduceRotation(-p.gRem * (ij / p.gRem), slots));

                    result[stop][ij] = rot;
                }
            }
        }
    }
    else {
        //------------------------------------------------------------------------------
        // sparsely-packed mode
        //------------------------------------------------------------------------------

        auto coeff  = CoeffEncodingCollapse(A, rotGroup, p.lvlb, false);
        auto coeffi = CoeffEncodingCollapse(A, rotGroup, p.lvlb, true);

        for (int32_t s = -1 + p.lvlb; s > stop; --s) {
            const int32_t rotScale = (1 << ((s - flagRem) * p.layersCollapse + p.remCollapse)) * p.g;
            const uint32_t limit   = p.b * p.g;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotations) {
                    // concatenate the coefficients horizontally on their third dimension, which corresponds to the # of slots
                    auto clearTmp   = coeff[s][ij];
                    auto& clearTmpi = coeffi[s][ij];
                    clearTmp.insert(clearTmp.end(), clearTmpi.begin(), clearTmpi.end());
                    // Maybe we need it one day
                    if ((flagRem == 0) && (s == stop + 1)) {
                        // do the scaling only at the last set of coefficients
                        //for (auto& c : clearTmp)
                        //    c *= scale;
                    }

                    auto rot = Rotate(clearTmp, ReduceRotation(-rotScale * (ij / p.g), M4));

                    result[s][ij] = rot;
                }
            }
        }

        if (flagRem == 1) {
            const uint32_t limit = p.bRem * p.gRem;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotationsRem) {
                    // concatenate the coefficients on their third dimension, which corresponds to the # of slots
                    auto clearTmp   = coeff[stop][ij];
                    auto& clearTmpi = coeffi[stop][ij];
                    clearTmp.insert(clearTmp.end(), clearTmpi.begin(), clearTmpi.end());
                    //for (auto& c : clearTmp)
                    //    c *= scale;

                    auto rot = Rotate(clearTmp, ReduceRotation(-p.gRem * (ij / p.gRem), M4));

                    result[stop][ij] = rot;
                }
            }
        }
    }
    return result;
}

std::vector<std::vector<ZBootstrapPlaintextCache>> FHEZImpl::EvalSlotsToCoeffsPrecompute(
    const CryptoContextImpl<DCRTPoly>& cc, const BigCVector& A, const std::vector<uint32_t>& rotGroup,
    bool flag_i) const {
    const uint32_t slots = rotGroup.size();

    const auto& p = GetBootPrecom(slots).m_paramsDec;

    const int32_t flagRem = (p.remCollapse == 0) ? 0 : 1;

    // result is the rotated plaintext version of coeff
    std::vector<std::vector<ZBootstrapPlaintextCache>> result(p.lvlb,
                                                              std::vector<ZBootstrapPlaintextCache>(p.numRotations));
    if (flagRem == 1) {
        // remainder corresponds to index 0 in encoding and to last index in decoding
        result[p.lvlb - 1].resize(p.numRotationsRem);
    }

    // make sure the plaintext is created only with the necessary amount of moduli

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());
    auto elementParams      = *(cryptoParams->GetElementParams());

    auto paramsQ          = elementParams.GetParams();
    auto paramsP          = cryptoParams->GetParamsP()->GetParams();
    auto elementParamsPtr = GetParamsQP(cc.GetCyclotomicOrder(), paramsQ, paramsP);

    if (uint32_t M4 = cc.GetCyclotomicOrder() / 4; M4 == slots) {
        // fully-packed
        auto coeff          = CoeffDecodingCollapse(A, rotGroup, p.lvlb, flag_i);
        const uint32_t smax = p.lvlb - flagRem;
        for (uint32_t s = 0; s < smax; ++s) {
            const int32_t rotScale = (1 << (s * p.layersCollapse)) * p.g;
            const uint32_t limit   = p.b * p.g;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotations) {
                    if ((flagRem == 0) && (s + 1 == smax)) {
                        // do the scaling only at the last set of coefficients
                        //for (auto& c : coeff[s][ij])
                        //    c *= scale;
                    }

                    auto rot = Rotate(coeff[s][ij], ReduceRotation(-rotScale * (ij / p.g), slots));

                    result[s][ij] = rot;
                }
            }
        }

        if (flagRem == 1) {
            const int32_t rotScale = (1 << (smax * p.layersCollapse)) * p.gRem;
            const uint32_t limit   = p.bRem * p.gRem;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotationsRem) {
                    //for (auto& c : coeff[smax][ij])
                    //    c *= scale;

                    auto rot = Rotate(coeff[smax][ij], ReduceRotation(-rotScale * (ij / p.g), slots));

                    result[smax][ij] = rot;
                }
            }
        }
    }
    else {
        //------------------------------------------------------------------------------
        // sparsely-packed mode
        //------------------------------------------------------------------------------

        auto coeff  = CoeffDecodingCollapse(A, rotGroup, p.lvlb, false);
        auto coeffi = CoeffDecodingCollapse(A, rotGroup, p.lvlb, true);

        const uint32_t smax = p.lvlb - flagRem;
        for (uint32_t s = 0; s < smax; ++s) {
            const int32_t rotScale = (1 << (s * p.layersCollapse)) * p.g;
            const uint32_t limit   = p.b * p.g;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotations) {
                    // concatenate the coefficients horizontally on their third dimension, which corresponds to the # of slots
                    auto clearTmp   = coeff[s][ij];
                    auto& clearTmpi = coeffi[s][ij];
                    clearTmp.insert(clearTmp.end(), clearTmpi.begin(), clearTmpi.end());
                    if ((flagRem == 0) && (s + 1 == smax)) {
                        // do the scaling only at the last set of coefficients
                        //for (auto& c : clearTmp)
                        //    c *= scale;
                    }

                    auto rot = Rotate(clearTmp, ReduceRotation(-rotScale * (ij / p.g), M4));

                    result[s][ij] = rot;
                }
            }
        }

        if (flagRem == 1) {
            const int32_t rotScale = (1 << (smax * p.layersCollapse)) * p.g;
            const uint32_t limit   = p.bRem * p.gRem;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotationsRem) {
                    // concatenate the coefficients on their third dimension, which corresponds to the # of slots
                    auto clearTmp   = coeff[smax][ij];
                    auto& clearTmpi = coeffi[smax][ij];
                    clearTmp.insert(clearTmp.end(), clearTmpi.begin(), clearTmpi.end());
                    //for (auto& c : clearTmp)
                    //    c *= scale;

                    auto rot = Rotate(clearTmp, ReduceRotation(-rotScale * (ij / p.g), M4));

                    result[smax][ij] = rot;
                }
            }
        }
    }
    return result;
}

Plaintext ZBootstrapPlaintextCache::GetPlaintext(const BigFixedPoint& scalingFactor,
                                                 const std::shared_ptr<typename DCRTPoly::Params>& elementParams) {
    auto q   = elementParams->GetModulus();
    auto key = std::make_tuple(scalingFactor, q);
    auto it  = m_cache.find(key);
    if (it != m_cache.end()) {
        return it->second;
    }

    // Construct plaintext
    auto slots = m_value.size();
    auto m     = m_value.size() * 4;
    auto n     = m_value.size() * 2;
    DiscreteFourierTransformBigComplex::Initialize(m, m / 4);

    BigCVector inverse = m_value;
    // Scale first
    for (auto& val : inverse) {
        val *= scalingFactor;
    }

    DiscreteFourierTransformBigComplex::FFTSpecialInv(inverse, m);

    std::vector<BigFixedPoint> rValues(2 * slots);
    for (size_t i = 0; i != inverse.size(); ++i) {
        rValues[i]         = inverse[i].getReal().round();
        rValues[i + slots] = inverse[i].getImag().round();
    }

    // The big one
    auto N = elementParams->GetRingDimension();
    BigVector V(N, q);
    for (size_t i = 0; i < n; ++i) {
        auto bfp     = rValues[i];
        auto integer = bfp.getValue() >> bfp.getLog2Scale();
        auto neg     = bfp.getNeg();
        if (neg) {
            V[i * N / n] = q.Sub(integer.Mod(q));
        }
        else {
            V[i * N / n] = integer.Mod(q);
        }
    }

    DCRTPoly::PolyLargeType polyLarge(std::make_shared<ILParamsImpl<DCRTPoly::Integer>>(2 * n, q, 1));
    polyLarge.SetValues(std::move(V), Format::COEFFICIENT);

    DCRTPoly poly(polyLarge, elementParams);
    poly.SetFormat(Format::EVALUATION);
    Plaintext ptxt = std::make_shared<ZEncodingImpl>(elementParams, poly, 0, scalingFactor);
    // Store in cache
    m_cache[key] = ptxt;
    return ptxt;
}

}  // namespace lbcrypto