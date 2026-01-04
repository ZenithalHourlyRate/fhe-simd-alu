#include "scheme/ckksrns/z-fhe.h"
#include "scheme/ckksrns/ckksrns-fhe.h"
#include "encoding/z-encoding.h"

double __heir_debug2(lbcrypto::ConstCiphertext<lbcrypto::DCRTPoly> ct, std::string msg) __attribute__((weak));

double __heir_debug2(lbcrypto::ConstCiphertext<lbcrypto::DCRTPoly> ct, std::string msg) {
    return 0.0;
}

namespace lbcrypto {

Ciphertext<DCRTPoly> FHEZImpl::EvalTruncate(ConstCiphertext<DCRTPoly>& ct) const {
    auto q     = ct->GetElements()[0].GetModulus();
    auto sfNow = ct->GetScalingFactorBFP();
    auto qBFP  = BigFixedPoint(q, 0, false).scaleTo(128);
    // q / Delta
    auto div       = (qBFP / sfNow).round();
    auto divScalar = div.getValue() >> div.getLog2Scale();
    auto ct2       = z->EvalMultScalar(ct, divScalar);
    ct2->SetScalingFactorBFP(qBFP);
    // Reduce all the way to the bottom
    z->ModReduceInPlace(ct2, ct2->GetElements()[0].GetNumOfElements() - 1);
    // To make sure the scaling factor is exactly q0 / 2
    auto q0     = ct2->GetElements()[0].GetModulus();
    auto q0BFP  = BigFixedPoint(q0, 0, false).scaleTo(128);
    auto sfNow2 = q0BFP;
    ct2->SetScalingFactorBFP(sfNow2);
    return ct2;
}

Ciphertext<DCRTPoly> FHEZImpl::EvalModRaise(ConstCiphertext<DCRTPoly>& ct) const {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ct->GetCryptoParameters());

    auto paramsQ   = cryptoParams->GetElementParams()->GetParams();
    uint32_t sizeQ = paramsQ.size();
    std::vector<NativeInteger> moduli(sizeQ);
    std::vector<NativeInteger> roots(sizeQ);
    for (uint32_t i = 0; i < sizeQ; ++i) {
        moduli[i] = paramsQ[i]->GetModulus();
        roots[i]  = paramsQ[i]->GetRootOfUnity();
    }

    auto cc = ct->GetCryptoContext();
    auto M  = cc->GetCyclotomicOrder();
    auto N  = cc->GetRingDimension();

    auto elementParamsRaisedPtr = std::make_shared<ILDCRTParams<DCRTPoly::Integer>>(M, moduli, roots);

    auto raised = ct->Clone();
    auto algo   = cc->GetScheme();

    uint32_t L0 = cryptoParams->GetElementParams()->GetParams().size();

    if (cryptoParams->GetSecretKeyDist() == SPARSE_ENCAPSULATED) {
        auto evalKeyMap = cc->GetEvalAutomorphismKeyMap(raised->GetKeyTag());

        // transform from a denser secret to a sparser one
        raised = FHECKKSRNS::KeySwitchSparse(raised, evalKeyMap.at(2 * N - 4));

        // Only level 0 ciphertext used here. Other towers ignored to make CKKS bootstrapping faster.
        auto& ctxtDCRTs = raised->GetElements();

        for (auto& dcrt : ctxtDCRTs) {
            dcrt.SetFormat(COEFFICIENT);
            DCRTPoly tmp(dcrt.GetElementAtIndex(0), elementParamsRaisedPtr);
            tmp.SetFormat(EVALUATION);
            dcrt = std::move(tmp);
        }
        raised->SetLevel(L0 - ctxtDCRTs[0].GetNumOfElements());

        // go back to a denser secret
        algo->KeySwitchInPlace(raised, evalKeyMap.at(2 * N - 2));
    }
    else {
        OPENFHE_THROW("Other Key unsupported");
    }

    raised->SetScalingFactorBFP(ct->GetScalingFactorBFP());
    return raised;
}

void FHEZImpl::EvalPartialSumInPlace(Ciphertext<DCRTPoly>& ct) const {
    auto cc     = ct->GetCryptoContext();
    auto cSlots = ct->GetZEncodingParams().getCSlots();
    auto N      = cc->GetRingDimension();

    const uint32_t limit = N / (cSlots * 2);
    for (uint32_t j = 1; j < limit; j <<= 1) {
        cc->EvalAddInPlace(ct, cc->EvalRotate(ct, j * (cSlots)));
    }
    // Now the message is multplied by N/(rN)
}

Ciphertext<DCRTPoly> FHEZImpl::EvalTruncateModRaisePartialSum(ConstCiphertext<DCRTPoly>& ct) const {
    auto ctNew = EvalModRaise(EvalTruncate(ct));
    EvalPartialSumInPlace(ctNew);
    return ctNew;
}

CiphertextGroup FHEZImpl::EvalZ2C(ConstCiphertext<DCRTPoly>& ct, Z2CScalingOption scalingOption, bool specialB0) const {
    auto cc       = ct->GetCryptoContext();
    auto cSlots   = ct->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    auto isSparse = precomp.m_isSparse;

    const CiphertextGroup::MapFunc postProcess = [&, this](ConstCiphertext<DCRTPoly>& z2c) -> Ciphertext<DCRTPoly> {
        // Take the Real
        auto ct = z->EvalAdd(z2c, z->EvalConjugateInC(z2c));
        z->ModReduceInPlace(ct);

        // For large zSlots, we do not scale down by zSlots during multiplying ZV
        auto zSlots = ct->GetZEncodingParams().getZSlots();
        if (precomp.m_zSlotsThresholdForScaling <= zSlots) {
            // Manually scale down by zSlots
            BigFixedPoint scaleDown = BigFixedPoint::one() / BigFixedPoint::positive(zSlots);
            if (scalingOption == SCALE_ZV_TWICE) {
                // Scale down by zSlots twice
                // To compensate for the zSlots scaling during encoding
                scaleDown = scaleDown / BigFixedPoint::positive(zSlots);
            }
            z->EvalMultInPlaceInC(ct, scaleDown);
            z->ModReduceInPlace(ct);
        }
        else {
            OPENFHE_THROW("Not implemented yet");
        }
        return ct;
    };

    // We force a sparse packing
    if (isSparse) {
        auto z2c = EvalZLinearTransform(specialB0 ? precomp.m_ZVSpecialB0Pre : precomp.m_ZVPre, ct);
        return postProcess(z2c);
    }
    else {
        auto z2c0 = EvalZLinearTransform(specialB0 ? precomp.m_ZV0SpecialB0Pre : precomp.m_ZV0Pre, ct);
        auto z2c1 = EvalZLinearTransform(precomp.m_ZV1Pre, ct);

        CiphertextGroup z2cGroup({z2c0, z2c1});
        return z2cGroup.map(postProcess);
    }
}

// The m0 + i * m1 preprocess is done elsewhere
Ciphertext<DCRTPoly> FHEZImpl::EvalC2R(ConstCiphertext<DCRTPoly>& ct) const {
    auto cc            = ct->GetCryptoContext();
    auto cSlots        = ct->GetZEncodingParams().getCSlots();
    auto precomp       = GetBootPrecom(cSlots);
    bool isSparse      = precomp.m_isSparse;
    bool isLTBootstrap = (precomp.m_paramsEnc.lvlb == 1) && (precomp.m_paramsDec.lvlb == 1);

    if (isSparse) {
        auto c2r = isLTBootstrap ? EvalLinearTransform(precomp.m_U0Pre, ct) : EvalSlotsToCoeffs(precomp.m_U0PreFFT, ct);
        // Trace
        z->EvalAddInPlace(c2r, cc->EvalRotate(c2r, cSlots));
        z->ModReduceInPlace(c2r);
        return c2r;
    }
    else {
        auto c2r =
            (isLTBootstrap) ? EvalLinearTransform(precomp.m_U0Pre, ct) : EvalSlotsToCoeffs(precomp.m_U0PreFFT, ct);
        z->ModReduceInPlace(c2r);
        return c2r;
    }
}

CiphertextGroup FHEZImpl::EvalR2C(ConstCiphertext<DCRTPoly>& ct, R2CScalingOption scalingOption) const {
    auto cc            = ct->GetCryptoContext();
    auto cSlots        = ct->GetZEncodingParams().getCSlots();
    auto precomp       = GetBootPrecom(cSlots);
    bool isLTBootstrap = (precomp.m_paramsEnc.lvlb == 1) && (precomp.m_paramsDec.lvlb == 1);
    auto isSparse      = precomp.m_isSparse;

    auto scaleDown = [&](Ciphertext<DCRTPoly> target, BigFixedPoint scale) {
        // Manually scale down by N
        auto sfNow     = target->GetScalingFactorBFP();
        auto scalarBFP = (sfNow * scale).round();
        auto scalar    = scalarBFP.getValue() >> scalarBFP.getLog2Scale();
        // Use EvalMultScalarInPlace to avoid precision loss
        // Should not use EvalMultInPlaceInC here as it will cause precision loss
        // TODO: should be OK to use EvalMultInPlaceInC?
        z->EvalMultScalarInPlace(target, scalar);
        target->SetScalingFactorBFP(sfNow * sfNow);
        z->ModReduceInPlace(target);
    };

    BigFixedPoint scaleN  = BigFixedPoint::one() / BigFixedPoint::positive(cc->GetRingDimension());
    BigFixedPoint scaleNK = BigFixedPoint::one() / (BigFixedPoint::positive(cc->GetRingDimension()) *
                                                    BigFixedPoint::positive(K_SPARSE_ENCAPSULATED));

    Ciphertext<DCRTPoly> r2c;

    // only one linear transform is needed as the other one can be derived
    if (scalingOption == SCALE_N || scalingOption == SCALE_NK) {
        r2c = (isLTBootstrap) ? EvalLinearTransform(precomp.m_U0hatTPre, ct) :
                                EvalCoeffsToSlots(precomp.m_U0hatTPreFFT, ct);
    }
    else if (scalingOption == SCALE_N_PRE) {
        r2c = (isLTBootstrap) ? EvalLinearTransform(precomp.m_U0hatTPreScaledN, ct) :
                                EvalCoeffsToSlots(precomp.m_U0hatTPreFFTScaledN, ct);
    }
    else {  // scalingOption == SCALE_NK_PRE
        r2c = (isLTBootstrap) ? EvalLinearTransform(precomp.m_U0hatTPreScaledNK, ct) :
                                EvalCoeffsToSlots(precomp.m_U0hatTPreFFTScaledNK, ct);
    }

    if (isSparse) {
        z->EvalAddInPlace(r2c, z->EvalConjugateInC(r2c));
        z->ModReduceInPlace(r2c);

        if (scalingOption == SCALE_N || scalingOption == SCALE_NK) {
            if (scalingOption == SCALE_N) {
                scaleDown(r2c, scaleN);
            }
            else {
                scaleDown(r2c, scaleNK);
            }
        }
        return r2c;
    }
    else {
        auto r2cConj = z->EvalConjugateInC(r2c);
        auto r2cI    = z->EvalSub(r2c, r2cConj);
        z->EvalAddInPlace(r2c, r2cConj);
        cc->GetScheme()->MultByMonomialInPlace(r2cI, 3 * cSlots);
        z->ModReduceInPlace(r2c);
        z->ModReduceInPlace(r2cI);

        if (scalingOption == SCALE_N || scalingOption == SCALE_NK) {
            if (scalingOption == SCALE_N) {
                scaleDown(r2c, scaleN);
                scaleDown(r2cI, scaleN);
            }
            else {
                scaleDown(r2c, scaleNK);
                scaleDown(r2cI, scaleNK);
            }
        }
        return std::vector{r2c, r2cI};
    }
}

Ciphertext<DCRTPoly> FHEZImpl::EvalC2Z(CiphertextGroup ct) const {
    auto cc       = ct[0]->GetCryptoContext();
    auto cSlots   = ct[0]->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    bool isSparse = precomp.m_isSparse;

    if (isSparse) {
        auto c2z = EvalZLinearTransform(precomp.m_ZUPre, ct[0]);
        z->EvalAddInPlace(c2z, cc->EvalRotate(c2z, cSlots));
        z->ModReduceInPlace(c2z);
        return c2z;
    }
    else {
        auto c2z0 = EvalZLinearTransform(precomp.m_ZU0Pre, ct[0]);
        auto c2z1 = EvalZLinearTransform(precomp.m_ZU1Pre, ct[1]);
        z->EvalAddInPlace(c2z0, c2z1);
        z->ModReduceInPlace(c2z0);
        return c2z0;
    }
}

Ciphertext<DCRTPoly> FHEZImpl::EvalZ2R(ConstCiphertext<DCRTPoly>& ct, Z2CScalingOption scalingOption) const {
    auto cc       = ct->GetCryptoContext();
    auto cSlots   = ct->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    auto isSparse = precomp.m_isSparse;

    if (isSparse) {
        return EvalC2R(EvalZ2C(ct, scalingOption)[0]);
    }
    else {
        auto z2cGroup = EvalZ2C(ct, scalingOption);
        auto z2c0     = z2cGroup[0];
        auto z2c1     = z2cGroup[1];
        cc->GetScheme()->MultByMonomialInPlace(z2c1, cSlots);
        cc->EvalAddInPlaceNoCheck(z2c0, z2c1);
        return EvalC2R(z2c0);
    }
}

Ciphertext<DCRTPoly> FHEZImpl::EvalR2Z(ConstCiphertext<DCRTPoly>& ct, R2CScalingOption scalingOption) const {
    return EvalC2Z(EvalR2C(ct, scalingOption));
}

Ciphertext<DCRTPoly> FHEZImpl::EvalArithToArithHigh(ConstCiphertext<DCRTPoly>& ct,
                                                    Z2CScalingOption scalingOption) const {
    auto cc      = ct->GetCryptoContext();
    auto cSlots  = ct->GetZEncodingParams().getCSlots();
    auto precomp = GetBootPrecom(cSlots);

    auto z2r = EvalZ2R(ct, scalingOption);

    //------------------------------------------------------------------------------
    // Truncate, ModRaise and PartialSum
    //------------------------------------------------------------------------------

    auto raised = EvalTruncateModRaisePartialSum(z2r);

    //------------------------------------------------------------------------------
    // R-To-Z
    //------------------------------------------------------------------------------

    // R2C then will multiply by rN then divide by N. Carried out in high precision
    auto r2z = EvalR2Z(raised, R2CScalingOption::SCALE_N);
    return r2z;
}

void FHEZImpl::ApplyDoubleAngleIterations(Ciphertext<DCRTPoly>& ct, uint32_t numIter) const {
    auto cc = ct->GetCryptoContext();
    for (int32_t i = 0; i != numIter; ++i) {
        ct = z->EvalMult(ct, ct);
        z->EvalAddInPlace(ct, z->EvalAddInC(ct, r_sparse_scalars[i]));
        z->ModReduceInPlace(ct);
    }
}

Ciphertext<DCRTPoly> FHEZImpl::EvalArithToArithNoise(ConstCiphertext<DCRTPoly>& ct,
                                                     Z2CScalingOption scalingOption) const {
    auto cc        = ct->GetCryptoContext();
    auto cSlots    = ct->GetZEncodingParams().getCSlots();
    auto precomp   = GetBootPrecom(cSlots);
    auto elemParam = ct->GetElements()[0].GetParams();
    auto sf        = ct->GetScalingFactorBFP();

    auto ctT = z->EvalMultTInZ(ct);
    z->ModReduceInPlace(ctT);

    auto z2r = EvalZ2R(ctT, scalingOption);
    //__heir_debug2(z2r, "Z2R");

    //------------------------------------------------------------------------------
    // Truncate, ModRaise and PartialSum
    //------------------------------------------------------------------------------

    auto raised = EvalTruncateModRaisePartialSum(z2r);
    // Now the message is multplied by N/(rN)
    //__heir_debug2(raised, "Raised");

    //------------------------------------------------------------------------------
    // R-To-C
    //------------------------------------------------------------------------------

    // R2C then will multiply by rN then divide by N * K; carried out in high precision
    auto r2c = EvalR2C(raised, R2CScalingOption::SCALE_NK);
    //__heir_debug2(r2c, "R2C");

    //------------------------------------------------------------------------------
    // Approximate Mod Reduction
    //------------------------------------------------------------------------------

    // Evaluate Chebyshev series for the sine wave
    auto& coeff_g0 = coeff_g0_big_complex_32;
    auto g0        = r2c.map([&](ConstCiphertext<DCRTPoly>& input) -> Ciphertext<DCRTPoly> {
        auto g0 = advZ->EvalChebyshevSeriesPS(input, coeff_g0);

        // Double-angle iterations
        uint32_t numIter = FHEZImpl::R_SPARSE;
        ApplyDoubleAngleIterations(g0, numIter);
        return g0;
    });

    //__heir_debug2(g0, "Sine");

    //------------------------------------------------------------------------------
    // C-To-Z
    //------------------------------------------------------------------------------

    auto r2z = EvalC2Z(g0);

    //------------------------------------------------------------------------------
    // Multiply by t^{-1} in Z
    //------------------------------------------------------------------------------

    r2z = z->EvalMultTInvInZ(r2z);
    z->ModReduceInPlace(r2z);

    return r2z;
}

Ciphertext<DCRTPoly> FHEZImpl::EvalArithToArith(ConstCiphertext<DCRTPoly>& ct, Z2CScalingOption scalingOption) const {
    auto high  = EvalArithToArithHigh(ct, scalingOption);
    auto noise = EvalArithToArithNoise(ct, scalingOption);
    return z->EvalSubWithAdjust(high, noise);
}

Ciphertext<DCRTPoly> FHEZImpl::EvalArithToBooleanSparse(ConstCiphertext<DCRTPoly>& ct, Z2CScalingOption scalingOption) {
    auto cc       = ct->GetCryptoContext();
    auto cSlots   = ct->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    bool isSparse = precomp.m_isSparse;
    if (!isSparse) {
        OPENFHE_THROW("ArithToBooleanSparse called on non-sparse ciphertext");
    }

    // Our core ct
    auto core = EvalZ2C(ct, scalingOption, /*specialB0=*/true)[0];
    // TODO: We need to keep core ct at bottom; rescale if necessary

    auto zN     = ct->GetZEncodingParams().getZN();
    auto zSlots = ct->GetZEncodingParams().getZSlots();
    uint32_t w  = precomp.m_w;
    // We ask zN to be multiple of w now...
    uint32_t numIter = static_cast<uint32_t>(std::ceil(static_cast<double>(zN) / (static_cast<double>(w))));

    std::vector<Ciphertext<DCRTPoly>> parts;

    // Iteratively process each low bits
    for (uint32_t iter = 0; iter != numIter; ++iter) {
        // Encode low-hot vector
        std::vector<BigComplex> oneHotVec(2 * cSlots, BigFixedPoint::zero());
        for (uint32_t j = 0; j != zSlots; ++j) {
            auto index = j * (zN / 2) + (iter * w);
            if (iter * w >= zN / 2) {
                index += (zSlots - 1) * (zN / 2);
            }
            for (uint32_t b = 0; b != w; ++b) {
                oneHotVec[index + b] = BigFixedPoint::one();
            }
        }

        ZEncodingParams oneHotZEncodeParams(CMode, zN * zSlots * 2);  // sparse packing
        RPolynomial oneHotPoly = CSlots(oneHotZEncodeParams, oneHotVec).toRPolynomial();

        // core may change over time
        auto elemParam       = core->GetElements()[0].GetParams();
        auto sf              = core->GetScalingFactorBFP();
        Plaintext oneHotPtxt = ZEncodingImpl::encodeR(oneHotPoly, elemParam, sf);

        auto coreMasked = z->EvalMult(core, oneHotPtxt);
        z->ModReduceInPlace(coreMasked);

        auto lut = internalBooleanToBooleanCustomLUTSparse(coreMasked, precomp.m_lutIDCoeffs);

        //------------------------------------------------------------------------------
        // Store the parts and remove it from core
        //------------------------------------------------------------------------------

        parts.push_back(lut);

        // remove part from core
        for (uint32_t nextIter = iter + 1; nextIter != numIter; ++nextIter) {
            int32_t diff = static_cast<int32_t>(iter) - static_cast<int32_t>(nextIter);
            // Note the rotation index is negative here
            int32_t rotationIndex = diff * w;
            if (rotationIndex <= precomp.m_cutoff) {
                // We do not remove them any more
                // Just treat the lower parts as noises
                break;
            }
            // Multiply by 1 / (2^{nextIter - iter} * p)
            auto scaled = z->EvalMultInC(lut, BigFixedPoint::pow2(diff * w));
            // If cross the half-way point, need to rotate more
            if (nextIter * w >= zN / 2 && iter * w < zN / 2) {
                rotationIndex -= static_cast<int32_t>((zSlots - 1) * zN / 2);
            }
            scaled = cc->EvalRotate(scaled, rotationIndex);
            z->ModReduceInPlace(scaled);
            z->EvalSubWithAdjustInPlace(core, scaled);
        }
        // Remove itself from core
        z->EvalSubWithAdjustInPlace(core, lut);
    }

    //------------------------------------------------------------------------------
    // Combine all parts
    //------------------------------------------------------------------------------

    for (size_t i = 1; i != parts.size(); ++i) {
        // They may have different scaling factors...
        z->EvalAddInPlace(parts[0], parts[i]);
    }

    return internalBooleanToBooleanCustomLUTSparse(parts[0], precomp.m_lutMSBCoeffs);
}

CiphertextGroup FHEZImpl::EvalArithToBooleanFull(ConstCiphertext<DCRTPoly>& ct, Z2CScalingOption scalingOption) {
    auto cc       = ct->GetCryptoContext();
    auto cSlots   = ct->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    bool isSparse = precomp.m_isSparse;
    if (isSparse) {
        OPENFHE_THROW("ArithToBooleanFull called on sparse ciphertext");
    }

    // Our core ct
    auto cores = EvalZ2C(ct, scalingOption, /*specialB0=*/true);
    auto core0 = cores[0];
    auto core1 = cores[1];
    // TODO: We need to keep core ct at bottom; rescale if necessary

    auto zN     = ct->GetZEncodingParams().getZN();
    auto zSlots = ct->GetZEncodingParams().getZSlots();
    uint32_t w  = precomp.m_w;
    // We ask zN to be multiple of w now...
    uint32_t numIter = static_cast<uint32_t>(std::ceil(static_cast<double>(zN) / (static_cast<double>(w))));

    std::vector<Ciphertext<DCRTPoly>> parts;

    // Iteratively process each low bits
    for (uint32_t iter = 0; iter != numIter; ++iter) {
        auto core       = core0;
        bool secondHalf = (iter * w >= zN / 2);
        if (secondHalf) {
            core = core1;
        }

        auto elemParam  = core->GetElements()[0].GetParams();
        auto sf         = core->GetScalingFactorBFP();
        auto maskPtxt   = getATBMaskFullPacking(iter, w, zN, zSlots, elemParam, sf);
        auto coreMasked = z->EvalMult(core, maskPtxt);
        z->ModReduceInPlace(coreMasked);

        auto lut = internalBooleanToBooleanCustomLUTSparse(coreMasked, precomp.m_lutIDCoeffs);

        //------------------------------------------------------------------------------
        // Store the parts and remove it from core
        //------------------------------------------------------------------------------

        parts.push_back(lut);

        // remove part from core
        for (uint32_t nextIter = iter + 1; nextIter != numIter; ++nextIter) {
            int32_t diff = static_cast<int32_t>(iter) - static_cast<int32_t>(nextIter);
            // Note the rotation index is negative here
            int32_t rotationIndex = diff * w;
            if (rotationIndex <= precomp.m_cutoff) {
                // We do not remove them any more
                // Just treat the lower parts as noises
                break;
            }
            // Multiply by 1 / (2^{nextIter - iter} * p)
            auto scaled = z->EvalMultInC(lut, BigFixedPoint::pow2(diff * w));

            bool targetSecondHalf = (nextIter * w >= zN / 2);

            // If cross the half-way point, switch to second ciphertext
            auto targetCore = core;
            if (targetSecondHalf && !secondHalf) {
                rotationIndex += zN / 2;
                targetCore = core1;
            }
            if (rotationIndex != 0) {
                scaled = cc->EvalRotate(scaled, rotationIndex);
            }
            z->ModReduceInPlace(scaled);
            z->EvalSubWithAdjustInPlace(targetCore, scaled);
        }
        // Remove itself from core
        z->EvalSubWithAdjustInPlace(core, lut);
    }

    //------------------------------------------------------------------------------
    // Combine all parts
    //------------------------------------------------------------------------------

    auto part0 = parts[0];
    auto part1 = parts[parts.size() / 2];
    for (size_t i = 1; i != parts.size() / 2; ++i) {
        // They may have different scaling factors...
        z->EvalAddInPlace(part0, parts[i]);
        z->EvalAddInPlace(part1, parts[parts.size() / 2 + i]);
    }
    CiphertextGroup partGroup({part0, part1});

    return internalBooleanToBooleanCustomLUTFull(partGroup, precomp.m_lutMSBCoeffs);
}

CiphertextGroup FHEZImpl::EvalArithToBooleanBatched(CiphertextGroup ctxts, Z2CScalingOption scalingOption) {
    auto cc       = ctxts[0]->GetCryptoContext();
    auto cSlots   = ctxts[0]->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    auto isSparse = precomp.m_isSparse;

    auto zN     = ctxts[0]->GetZEncodingParams().getZN();
    auto zSlots = ctxts[0]->GetZEncodingParams().getZSlots();
    uint32_t w  = precomp.m_w;
    // We ask zN to be multiple of w now...
    uint32_t numIter = static_cast<uint32_t>(std::ceil(static_cast<double>(zN) / (static_cast<double>(w))));

    if (isSparse) {
        OPENFHE_THROW("Batched EvalArithToBooleanBatched only supports full packing");
    }
    if (ctxts.getParts().size() != numIter) {
        OPENFHE_THROW("Batch size mismatch for batched EvalArithToBooleanBatched");
    }

    auto core =
        ctxts.mapWide([&](ConstCiphertext<DCRTPoly>& ct) { return EvalZ2C(ct, scalingOption, /*specialB0=*/true); });

    std::vector<Ciphertext<DCRTPoly>> coreFirstHalfVec;
    std::vector<Ciphertext<DCRTPoly>> coreSecondHalfVec;
    for (auto i = 0; i != ctxts.getParts().size(); ++i) {
        coreFirstHalfVec.push_back(core[2 * i]);
        coreSecondHalfVec.push_back(core[2 * i + 1]);
    }
    CiphertextGroup coreFirstHalf(coreFirstHalfVec);
    CiphertextGroup coreSecondHalf(coreSecondHalfVec);

    std::vector<Ciphertext<DCRTPoly>> luts;

    // Iteratively process each low bits
    for (uint32_t iter = 0; iter != numIter; ++iter) {
        bool secondHalf = (iter * w >= zN / 2);

        CiphertextGroup targetGroup = coreFirstHalf;
        if (secondHalf) {
            targetGroup = coreSecondHalf;
        }

        // core may change over time
        auto elemParam = targetGroup[0]->GetElements()[0].GetParams();
        auto sf        = targetGroup[0]->GetScalingFactorBFP();
        auto maskPtxt  = getATBMask(iter, w, zN, zSlots, elemParam, sf);

        auto masked = targetGroup.map([&](ConstCiphertext<DCRTPoly>& ct) {
            auto newCt = z->EvalMult(ct, maskPtxt);
            z->ModReduceInPlace(newCt);
            return newCt;
        });

        auto targetCombined0 = masked[0];
        auto targetCombined1 = masked[masked.size() / 2];
        for (size_t i = 1; i != masked.size() / 2; ++i) {
            auto rotated0 = cc->EvalRotate(masked[i], -static_cast<int32_t>(i * w));
            z->EvalAddInPlace(targetCombined0, rotated0);

            auto rotated1 = cc->EvalRotate(masked[masked.size() / 2 + i], -static_cast<int32_t>(i * w));
            z->EvalAddInPlace(targetCombined1, rotated1);
        }

        auto lut =
            internalBooleanToBooleanCustomLUTFull(std::vector{targetCombined0, targetCombined1}, precomp.m_lutIDCoeffs);
        auto lut0 = lut[0];
        auto lut1 = lut[1];

        //------------------------------------------------------------------------------
        // Store the parts and remove it from core
        //------------------------------------------------------------------------------

        luts.push_back(lut0);
        luts.push_back(lut1);

        elemParam = lut0->GetElements()[0].GetParams();
        sf        = lut0->GetScalingFactorBFP();

        // remove part from core
        for (uint32_t nextIter = iter + 1; nextIter != numIter; ++nextIter) {
            int32_t diff = static_cast<int32_t>(iter) - static_cast<int32_t>(nextIter);
            // Note the rotation index is negative here
            int32_t rotationIndex = diff * w;
            if (rotationIndex <= precomp.m_cutoff) {
                // We do not remove them any more
                // Just treat the lower parts as noises
                break;
            }
            if (nextIter * w >= zN / 2 && iter * w < zN / 2) {
                rotationIndex += zN / 2;
            }
            for (size_t j = 0; j != targetGroup.size() / 2; ++j) {
                // Multiply by 1 / (2^{nextIter - iter} * p)
                auto scaleDown      = BigFixedPoint::pow2(diff * w);
                auto recombMaskPtxt = getATBRecombMaskFullPacking(iter, w, zN, zSlots, elemParam, sf, scaleDown,
                                                                  -static_cast<int32_t>(j * w));

                {
                    auto scaled            = z->EvalMult(lut0, recombMaskPtxt);
                    auto targetRotateIndex = rotationIndex + static_cast<int32_t>(j * w);
                    if (targetRotateIndex != 0) {
                        scaled = cc->EvalRotate(scaled, targetRotateIndex);
                    }
                    z->ModReduceInPlace(scaled);

                    bool targetSecondHalf = (nextIter * w >= zN / 2);

                    auto coreCt = core[2 * j + targetSecondHalf];
                    z->EvalSubWithAdjustInPlace(coreCt, scaled);
                }
                {
                    auto scaled            = z->EvalMult(lut1, recombMaskPtxt);
                    auto targetRotateIndex = rotationIndex + static_cast<int32_t>(j * w);
                    if (targetRotateIndex != 0) {
                        scaled = cc->EvalRotate(scaled, targetRotateIndex);
                    }
                    z->ModReduceInPlace(scaled);

                    bool targetSecondHalf = (nextIter * w >= zN / 2);
                    auto coreCt           = core[2 * (targetGroup.size() / 2 + j) + targetSecondHalf];
                    z->EvalSubWithAdjustInPlace(coreCt, scaled);
                }
            }
        }
        // Remove itself from core
        for (size_t j = 0; j != targetGroup.size() / 2; ++j) {
            auto recombMaskPtxt = getATBRecombMaskFullPacking(iter, w, zN, zSlots, elemParam, sf, BigFixedPoint::one(),
                                                              -static_cast<int32_t>(j * w));
            {
                auto scaled = z->EvalMult(lut0, recombMaskPtxt);
                if (static_cast<int32_t>(j * w) != 0) {
                    scaled = cc->EvalRotate(scaled, static_cast<int32_t>(j * w));
                }
                z->ModReduceInPlace(scaled);

                bool targetSecondHalf = secondHalf;
                auto coreCt           = core[2 * j + targetSecondHalf];
                z->EvalSubWithAdjustInPlace(coreCt, scaled);
            }
            {
                auto scaled = z->EvalMult(lut1, recombMaskPtxt);
                if (static_cast<int32_t>(j * w) != 0) {
                    scaled = cc->EvalRotate(scaled, static_cast<int32_t>(j * w));
                }
                z->ModReduceInPlace(scaled);

                bool targetSecondHalf = secondHalf;
                auto coreCt           = core[2 * (targetGroup.size() / 2 + j) + targetSecondHalf];
                z->EvalSubWithAdjustInPlace(coreCt, scaled);
            }
        }
    }

    //------------------------------------------------------------------------------
    // Combine all LUTs
    //------------------------------------------------------------------------------

    // For each input ciphertext, extract MSBs from luts
    std::vector<Ciphertext<DCRTPoly>> MSBs;
    std::vector<Ciphertext<DCRTPoly>> MSBsSec;
    for (size_t j = 0; j != ctxts.size() / 2; ++j) {
        Ciphertext<DCRTPoly> MSBjFirst, MSBjSecond, MSBSecjFirst, MSBSecjSecond;
        for (size_t iter = 0; iter != numIter / 2; ++iter) {
            auto mask = getATBRecombMaskFullPacking(iter, w, zN, zSlots, luts[0]->GetElements()[0].GetParams(),
                                                    luts[0]->GetScalingFactorBFP(), BigFixedPoint::one(),
                                                    -static_cast<int32_t>(j * w));
            {
                auto maskedLut = z->EvalMult(luts[2 * iter], mask);
                z->ModReduceInPlace(maskedLut);
                if (iter == 0) {
                    MSBjFirst = maskedLut;
                }
                else {
                    z->EvalAddInPlace(MSBjFirst, maskedLut);
                }
            }
            {
                auto maskedLut = z->EvalMult(luts[2 * iter + 1], mask);
                z->ModReduceInPlace(maskedLut);
                if (iter == 0) {
                    MSBSecjFirst = maskedLut;
                }
                else {
                    z->EvalAddInPlace(MSBSecjFirst, maskedLut);
                }
            }
        }
        if (j != 0) {
            MSBjFirst    = cc->EvalRotate(MSBjFirst, static_cast<int32_t>(j * w));
            MSBSecjFirst = cc->EvalRotate(MSBSecjFirst, static_cast<int32_t>(j * w));
        }
        for (size_t iter = numIter / 2; iter != numIter; ++iter) {
            auto mask = getATBRecombMaskFullPacking(iter, w, zN, zSlots, luts[0]->GetElements()[0].GetParams(),
                                                    luts[0]->GetScalingFactorBFP(), BigFixedPoint::one(),
                                                    -static_cast<int32_t>(j * w));
            {
                auto maskedLut = z->EvalMult(luts[2 * iter], mask);
                z->ModReduceInPlace(maskedLut);
                if (iter == numIter / 2) {
                    MSBjSecond = maskedLut;
                }
                else {
                    z->EvalAddInPlace(MSBjSecond, maskedLut);
                }
            }
            {
                auto maskedLut = z->EvalMult(luts[2 * iter + 1], mask);
                z->ModReduceInPlace(maskedLut);
                if (iter == numIter / 2) {
                    MSBSecjSecond = maskedLut;
                }
                else {
                    z->EvalAddInPlace(MSBSecjSecond, maskedLut);
                }
            }
        }
        if (j != 0) {
            MSBjSecond    = cc->EvalRotate(MSBjSecond, static_cast<int32_t>(j * w));
            MSBSecjSecond = cc->EvalRotate(MSBSecjSecond, static_cast<int32_t>(j * w));
        }
        MSBs.push_back(MSBjFirst);
        MSBs.push_back(MSBjSecond);
        MSBsSec.push_back(MSBSecjFirst);
        MSBsSec.push_back(MSBSecjSecond);
    }

    MSBs.insert(MSBs.end(), MSBsSec.begin(), MSBsSec.end());

    std::vector<Ciphertext<DCRTPoly>> finalMSBs;
    for (size_t i = 0; i != MSBs.size() / 2; ++i) {
        auto boolGroup =
            internalBooleanToBooleanCustomLUTFull(std::vector{MSBs[2 * i], MSBs[2 * i + 1]}, precomp.m_lutMSBCoeffs);
        finalMSBs.push_back(boolGroup[0]);
        finalMSBs.push_back(boolGroup[1]);
    }
    return finalMSBs;
}

Ciphertext<DCRTPoly> FHEZImpl::internalBooleanToBooleanLTsSparse(ConstCiphertext<DCRTPoly>& ct) const {
    auto cc      = ct->GetCryptoContext();
    auto cSlots  = ct->GetZEncodingParams().getCSlots();
    auto precomp = GetBootPrecom(cSlots);

    auto c2r = EvalC2R(ct);

    //------------------------------------------------------------------------------
    // Truncate, ModRaise and PartialSum
    //------------------------------------------------------------------------------

    auto raised = EvalTruncateModRaisePartialSum(c2r);
    // Now the message is multplied by N/(rN)

    //------------------------------------------------------------------------------
    // R-To-C
    //------------------------------------------------------------------------------

    // R2C then will multiply by rN then divide by N * K because of scaling in pre-compute
    // can do so in low precision as we do not need high precision here
    // Since we call C2R without imaginary part, we only need the real part here
    auto r2c = EvalR2C(raised, SCALE_NK_PRE)[0];
    return r2c;
}

CiphertextGroup FHEZImpl::internalBooleanToBooleanLTsFull(CiphertextGroup ct) const {
    auto cc      = ct[0]->GetCryptoContext();
    auto cSlots  = ct[0]->GetZEncodingParams().getCSlots();
    auto precomp = GetBootPrecom(cSlots);

    auto ctComb = ct[0]->Clone();
    auto ct1I   = cc->GetScheme()->MultByMonomial(ct[1], cSlots);
    z->EvalAddInPlace(ctComb, ct1I);
    auto c2r = EvalC2R(ctComb);

    //------------------------------------------------------------------------------
    // Truncate, ModRaise and PartialSum
    //------------------------------------------------------------------------------

    auto raised = EvalTruncateModRaisePartialSum(c2r);

    //------------------------------------------------------------------------------
    // R-To-C
    //------------------------------------------------------------------------------

    // R2C then will multiply by rN then divide by N * K because of scaling in pre-compute
    // can do so in low precision as we do not need high precision here
    auto r2cGroup = EvalR2C(raised, SCALE_NK_PRE);
    return r2cGroup;
}

Ciphertext<DCRTPoly> FHEZImpl::internalBooleanToBooleanCustomLUTSparse(ConstCiphertext<DCRTPoly>& ct,
                                                                       const std::vector<BigComplex>& lutCoeffs) const {
    auto r2c = internalBooleanToBooleanLTsSparse(ct);

    //------------------------------------------------------------------------------
    // Exp
    //------------------------------------------------------------------------------

    auto& coeff_exp = coeff_exp_16_big_complex_46;
    auto res        = advZ->EvalChebyshevSeriesPS(r2c, coeff_exp);

    // Double angle-iterations to get exp(2*Pi*i*x)
    res = z->EvalSquare(res);
    z->ModReduceInPlace(res);
    res = z->EvalSquare(res);
    z->ModReduceInPlace(res);

    //------------------------------------------------------------------------------
    // Running LUT
    //------------------------------------------------------------------------------

    auto powers = advZ->EvalPowers(res, lutCoeffs);
    auto lut    = advZ->EvalPolyWithPrecomp(powers, lutCoeffs);
    // Take the real part
    z->EvalAddInPlace(lut, z->EvalConjugateInC(lut));
    return lut;
}

CiphertextGroup FHEZImpl::internalBooleanToBooleanCustomLUTFull(CiphertextGroup ctGroup,
                                                                const std::vector<BigComplex>& lutCoeffs) const {
    auto r2c = internalBooleanToBooleanLTsFull(ctGroup);

    return r2c.map([&](ConstCiphertext<DCRTPoly> ct) {
        //------------------------------------------------------------------------------
        // Exp
        //------------------------------------------------------------------------------

        auto& coeff_exp = coeff_exp_16_big_complex_46;
        auto res        = advZ->EvalChebyshevSeriesPS(ct, coeff_exp);

        // Double angle-iterations to get exp(2*Pi*i*x)
        res = z->EvalSquare(res);
        z->ModReduceInPlace(res);
        res = z->EvalSquare(res);
        z->ModReduceInPlace(res);

        //------------------------------------------------------------------------------
        // Running LUT
        //------------------------------------------------------------------------------

        auto powers = advZ->EvalPowers(res, lutCoeffs);
        auto lut    = advZ->EvalPolyWithPrecomp(powers, lutCoeffs);
        // Take the real part
        z->EvalAddInPlace(lut, z->EvalConjugateInC(lut));
        return lut;
    });
}

Ciphertext<DCRTPoly> FHEZImpl::EvalBooleanToBooleanSparse(ConstCiphertext<DCRTPoly>& ct) const {
    auto ctNew = ct->Clone();
    // Double scaling factor to divide by 2
    ctNew->SetScalingFactorBFP(ctNew->GetScalingFactorBFP() * BigFixedPoint::two());
    auto r2c = internalBooleanToBooleanLTsSparse(ctNew);

    //------------------------------------------------------------------------------
    // Cosine
    //------------------------------------------------------------------------------

    auto& coeff_cos = coeff_cos_16_big_complex_50;
    auto res        = advZ->EvalChebyshevSeriesPS(r2c, coeff_cos);

    // Double angle-iterations to get cos(pi*x)
    res = z->EvalSquare(res);
    z->EvalAddInPlace(res, res);
    z->EvalAddInPlaceInC(res, -BigFixedPoint::one());
    z->ModReduceInPlace(res);  // cos(pi x)
    res = z->EvalSquare(res);
    z->ModReduceInPlace(res);  // cos^2(pi x)

    //------------------------------------------------------------------------------
    // The LUT (1 - cos^2(pi x))
    // This is p == 2 and order == 1 of AKP25 for the identity map
    // f(0) = 0, f(1) = 1
    //------------------------------------------------------------------------------

    z->EvalNegateInPlace(res);                        // -cos^2(pi x)
    z->EvalAddInPlaceInC(res, BigFixedPoint::one());  // 1 - cos^2(pi x)

    return res;
}

CiphertextGroup FHEZImpl::EvalBooleanToBooleanFull(CiphertextGroup ct) const {
    auto ctNew0 = ct[0]->Clone();
    auto ctNew1 = ct[1]->Clone();
    // Double scaling factor to divide by 2
    ctNew0->SetScalingFactorBFP(ctNew0->GetScalingFactorBFP() * BigFixedPoint::two());
    ctNew1->SetScalingFactorBFP(ctNew1->GetScalingFactorBFP() * BigFixedPoint::two());
    auto r2c = internalBooleanToBooleanLTsFull(std::vector{ctNew0, ctNew1});

    return r2c.map([&](ConstCiphertext<DCRTPoly> ct) {
        //------------------------------------------------------------------------------
        // Cosine
        //------------------------------------------------------------------------------

        auto& coeff_cos = coeff_cos_16_big_complex_50;
        auto res        = advZ->EvalChebyshevSeriesPS(ct, coeff_cos);

        // Double angle-iterations to get cos(pi*x)
        res = z->EvalSquare(res);
        z->EvalAddInPlace(res, res);
        z->EvalAddInPlaceInC(res, -BigFixedPoint::one());
        z->ModReduceInPlace(res);  // cos(pi x)
        res = z->EvalSquare(res);
        z->ModReduceInPlace(res);  // cos^2(pi x)

        //------------------------------------------------------------------------------
        // The LUT (1 - cos^2(pi x))
        // This is p == 2 and order == 1 of AKP25 for the identity map
        // f(0) = 0, f(1) = 1
        //------------------------------------------------------------------------------

        z->EvalNegateInPlace(res);                        // -cos^2(pi x)
        z->EvalAddInPlaceInC(res, BigFixedPoint::one());  // 1 - cos^2(pi x)

        return res;
    });
}

}  // namespace lbcrypto