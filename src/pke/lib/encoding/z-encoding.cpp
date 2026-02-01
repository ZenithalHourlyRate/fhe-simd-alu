#include "encoding/z-encoding.h"

namespace lbcrypto {

//#define HIGH_PREC

#ifndef HIGH_PREC
ZEncoding ZEncodingImpl::encodeR(const RPolynomial& input,
                                 const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                 const BigFixedPoint& scalingFactor) {
    if (scalingFactor.log2Norm() > 60) {
        OPENFHE_THROW("Scaling factor too large in ZEncodingImpl::encodeR");
    }
    auto N           = elementParams->GetRingDimension();
    auto n           = input.getCoefficients().size();
    auto inputCoeffs = input.getCoefficients();
    std::vector<int64_t> roundedCoeffs;
    roundedCoeffs.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        roundedCoeffs.push_back((inputCoeffs[i] * scalingFactor).getRoundedInteger() *
                                (inputCoeffs[i].getNeg() ? -1 : 1));
    }
    DCRTPoly poly(elementParams, Format::COEFFICIENT, true);

    auto& mVectors = poly.GetAllElements();
    auto t         = poly.GetNumOfElements();
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(t))
    for (size_t i = 0; i < t; ++i) {
        auto& singlePoly = mVectors[i];
        auto qi          = mVectors[i].GetModulus();
        int64_t qiInt    = qi.template ConvertToInt<int64_t>();
        for (size_t j = 0; j < n; ++j) {
            auto signedRem = roundedCoeffs[j] % qiInt;
            if (signedRem < 0) {
                signedRem += qiInt;
            }
            singlePoly[j * N / n] = static_cast<uint64_t>(signedRem);
        }
    }
    poly.SetFormat(Format::EVALUATION);
    return std::make_shared<ZEncodingImpl>(elementParams, poly, scalingFactor, input.getZEncodingParams());
}

#else

ZEncoding ZEncodingImpl::encodeR(const RPolynomial& input,
                                 const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                 const BigFixedPoint& scalingFactor) {
    auto N              = elementParams->GetRingDimension();
    auto q              = elementParams->GetModulus();
    auto n              = input.getCoefficients().size();
    const auto& rCoeffs = input.getCoefficients();
    BigVector V(N, q);
    for (size_t i = 0; i < n; ++i) {
        auto bfp     = (rCoeffs[i] * scalingFactor);
        auto integer = bfp.getRoundedBigInteger();
        auto neg     = bfp.getNeg();
        if (neg) {
            V[i * N / n] = q.Sub(integer.Mod(q));
        }
        else {
            V[i * N / n] = integer.Mod(q);
        }
    }

    DCRTPoly::PolyLargeType polyLarge(std::make_shared<ILParamsImpl<DCRTPoly::Integer> >(2 * N, q, 1));
    polyLarge.SetValues(std::move(V), Format::COEFFICIENT);

    DCRTPoly poly(polyLarge, elementParams);
    poly.SetFormat(Format::EVALUATION);
    return std::make_shared<ZEncodingImpl>(elementParams, poly, scalingFactor, input.getZEncodingParams());
}
#endif

ZEncoding ZEncodingImpl::encodeC(const BigComplex& input,
                                 const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                 const BigFixedPoint& scalingFactor) {
    if (scalingFactor.log2Norm() > 60) {
        OPENFHE_THROW("Scaling factor too large in ZEncodingImpl::encodeC");
    }
    auto N              = elementParams->GetRingDimension();
    int64_t roundedReal = (input.getReal() * scalingFactor).getRoundedInteger() * (input.getReal().getNeg() ? -1 : 1);
    int64_t roundedImag = (input.getImag() * scalingFactor).getRoundedInteger() * (input.getImag().getNeg() ? -1 : 1);
    DCRTPoly poly(elementParams, Format::COEFFICIENT, true);

    auto& mVectors = poly.GetAllElements();
    auto t         = poly.GetNumOfElements();
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(t))
    for (size_t i = 0; i < t; ++i) {
        auto& singlePoly = mVectors[i];
        auto qi          = mVectors[i].GetModulus();
        int64_t qiInt    = qi.template ConvertToInt<int64_t>();

        // Real part in coeff 0
        {
            auto signedRem = roundedReal % qiInt;
            if (signedRem < 0) {
                signedRem += qiInt;
            }
            singlePoly[0] = static_cast<uint64_t>(signedRem);
        }
        // Imag part in coeff N/2
        {
            auto signedRem = roundedImag % qiInt;
            if (signedRem < 0) {
                signedRem += qiInt;
            }
            singlePoly[N / 2] = static_cast<uint64_t>(signedRem);
        }
    }
    poly.SetFormat(Format::EVALUATION);
    ZEncodingParams params(CMode, 2);
    return std::make_shared<ZEncodingImpl>(elementParams, poly, scalingFactor, params);

    /*
        auto n = elementParams->GetRingDimension();
        auto q = elementParams->GetModulus();

        // for a value in C, its representation in R is
        // 1. Real part in constant coeff
        // 2. Imaginary part in coeff of n/2
        auto scaledValueComplex = input * scalingFactor;
        BigVector V(n, q);
        // Real part in V[0]
        {
            auto realBFP     = scaledValueComplex.getReal();
            auto realInteger = realBFP.getRoundedBigInteger();
            auto realNeg     = realBFP.getNeg();
            if (realNeg) {
                V[0] = q.Sub(realInteger.Mod(q));
            }
            else {
                V[0] = realInteger.Mod(q);
            }
        }
        // Imag part in V[n//2]
        {
            auto imagBFP     = scaledValueComplex.getImag();
            auto imagInteger = imagBFP.getRoundedBigInteger();
            auto imagNeg     = imagBFP.getNeg();
            if (imagNeg) {
                V[n / 2] = q.Sub(imagInteger.Mod(q));
            }
            else {
                V[n / 2] = imagInteger.Mod(q);
            }
        }

        DCRTPoly::PolyLargeType polyLarge(std::make_shared<ILParamsImpl<DCRTPoly::Integer>>(2 * n, q, 1));
        polyLarge.SetValues(std::move(V), Format::COEFFICIENT);

        DCRTPoly poly(polyLarge, elementParams);
        poly.SetFormat(Format::EVALUATION);
        ZEncodingParams params(CMode, 2);
        return std::make_shared<ZEncodingImpl>(elementParams, poly, scalingFactor, params);
        */
}

}  // namespace lbcrypto