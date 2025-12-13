#ifndef SRC_PKE_EXAMPLES_Z_ENCODE_H_
#define SRC_PKE_EXAMPLES_Z_ENCODE_H_

#include <cassert>
#include "ciphertext-fwd.h"
#include "math/z-polynomial.h"
#include "plaintext.h"

namespace lbcrypto {

class ZEncodingImpl;

using ZEncoding = std::shared_ptr<ZEncodingImpl>;

// This is only interplay between Plaintext and RPolynomial
// RPolynomial to CSlots to ZPolynomial is handled elsewhere
class ZEncodingImpl : public PlaintextImpl {
private:
    BigFixedPoint m_scalingFactorBFP;
    ZEncodingParams m_zEncodingParams;

public:
    ZEncodingImpl(std::shared_ptr<DCRTPoly::Params> vp, const DCRTPoly& elements, BigFixedPoint scalingFactorBFP,
                  const ZEncodingParams& params)
        : PlaintextImpl(vp, nullptr, INVALID_ENCODING, INVALID_SCHEME) {
        encodedVectorDCRT  = elements;
        m_scalingFactorBFP = scalingFactorBFP;
        m_zEncodingParams  = params;
    }

    // Just to fullfill the PlaintextImpl interface
    virtual bool Encode() override {
        OPENFHE_THROW("not implemented");
    }
    virtual bool Decode() override {
        OPENFHE_THROW("not implemented");
    }
    virtual void PrintValue(std::ostream& out) const override {
        OPENFHE_THROW("Not implemented");
    }
    virtual bool CompareTo(const PlaintextImpl& other) const override {
        OPENFHE_THROW("Not implemented");
    }
    virtual size_t GetLength() const override {
        return m_zEncodingParams.getLength();
    }
    BigFixedPoint GetScalingFactorBFP() const {
        return m_scalingFactorBFP;
    }
    ZEncodingParams GetZEncodingParams() const {
        return m_zEncodingParams;
    }

    static ZEncoding encodeR(const RPolynomial& input, const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                             const BigFixedPoint& scalingFactor) {
        auto N              = elementParams->GetRingDimension();
        auto q              = elementParams->GetModulus();
        auto n              = input.getCoefficients().size();
        const auto& rCoeffs = input.getCoefficients();
        BigVector V(N, q);
        for (size_t i = 0; i < n; ++i) {
            auto bfp     = (rCoeffs[i] * scalingFactor).round();
            auto integer = bfp.getValue() >> bfp.getLog2Scale();
            auto neg     = bfp.getNeg();
            if (neg) {
                V[i * N / n] = q.Sub(integer.Mod(q));
            }
            else {
                V[i * N / n] = integer.Mod(q);
            }
        }

        DCRTPoly::PolyLargeType polyLarge(std::make_shared<ILParamsImpl<DCRTPoly::Integer>>(2 * N, q, 1));
        polyLarge.SetValues(std::move(V), Format::COEFFICIENT);

        DCRTPoly poly(polyLarge, elementParams);
        poly.SetFormat(Format::EVALUATION);
        return std::make_shared<ZEncodingImpl>(elementParams, poly, scalingFactor, input.getZEncodingParams());
    }

    static RPolynomial decodeR(ZEncoding input) {
        auto dcrtPoly      = input->GetElement<DCRTPoly>();
        auto outputSize    = input->GetLength();
        auto scalingFactor = input->GetScalingFactorBFP();
        dcrtPoly.SetFormat(Format::COEFFICIENT);
        auto bigPoly = dcrtPoly.CRTInterpolate();
        std::vector<BigFixedPoint> output;
        auto ringDim = dcrtPoly.GetParams()->GetRingDimension();
        auto q       = dcrtPoly.GetParams()->GetModulus();
        for (size_t i = 0; i != outputSize; ++i) {
            auto valInteger = bigPoly[i * (ringDim / outputSize)];
            bool neg        = false;
            if (valInteger > q / 2) {
                neg        = true;
                valInteger = q - valInteger;
            }
            auto valFixedPoint = BigFixedPoint(valInteger, 0, neg).scaleTo(128);
            output.push_back(valFixedPoint / scalingFactor);
        }
        return RPolynomial(input->GetZEncodingParams(), output);
    }

    static ZEncoding encodeC(const BigComplex& input, const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                             const BigFixedPoint& scalingFactor) {
        auto n = elementParams->GetRingDimension();
        auto q = elementParams->GetModulus();

        // for a value in C, its representation in R is
        // 1. Real part in constant coeff
        // 2. Imaginary part in coeff of n/2
        auto scaledValueComplex = input * scalingFactor;
        BigVector V(n, q);
        // Real part in V[0]
        {
            auto realBFP     = scaledValueComplex.getReal().round();
            auto realInteger = realBFP.getValue() >> realBFP.getLog2Scale();
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
            auto imagBFP     = scaledValueComplex.getImag().round();
            auto imagInteger = imagBFP.getValue() >> imagBFP.getLog2Scale();
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
    }

    static ZEncoding encodeC(const CSlots& cSlots, const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                             const BigFixedPoint& scalingFactor) {
        return encodeR(cSlots.toRPolynomial(), elementParams, scalingFactor);
    }

    static ZEncoding encodeOneHotInC(uint32_t zN, uint32_t zSlots, uint32_t bit,
                                     const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                     const BigFixedPoint& scalingFactor, uint32_t downScale) {
        // TODO: handle sparse packing and full packing
        ZEncodingParams params(CMode, zN * zSlots * 2);  // for sparse packing now...
        std::vector<BigComplex> vec(zN * zSlots, BigFixedPoint::zero());
        // TODO: fix it
        vec[bit - 1] = BigFixedPoint::one() / BigFixedPoint::positive(downScale);
        CSlots cSlots(params, vec);
        return encodeC(cSlots, elementParams, scalingFactor);
    }
};

}  // namespace lbcrypto

#endif  // SRC_PKE_EXAMPLES_Z_ENCODE_H_