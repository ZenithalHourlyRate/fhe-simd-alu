#ifndef SRC_PKE_EXAMPLES_Z_ENCODE_H_
#define SRC_PKE_EXAMPLES_Z_ENCODE_H_

#include <cassert>
#include "math/z-polynomial.h"
#include "plaintext.h"

namespace lbcrypto {

enum ZEncodingType {
    Z = 1,  // Representing ZPolynomial, at Z/arithmetic mode
    C       // Representing C values. Used at bootstrapping and C/binary mode
};

class ZEncodingImpl;

using ZEncoding = std::shared_ptr<ZEncodingImpl>;

// This is only interplay between Plaintext and RPolynomial
// RPolynomial to CSlots to ZPolynomial is handled elsewhere
class ZEncodingImpl : public PlaintextImpl {
private:
    BigFixedPoint m_scalingFactorBFP;
    uint32_t m_zN;
    uint32_t m_zSlots;

public:
    ZEncodingImpl(std::shared_ptr<DCRTPoly::Params> vp, const DCRTPoly& elements, uint32_t zN, uint32_t zSlots,
                  BigFixedPoint scalingFactorBFP)
        : PlaintextImpl(vp, nullptr, INVALID_ENCODING, INVALID_SCHEME) {
        encodedVectorDCRT  = elements;
        m_zN               = zN;
        m_zSlots           = zSlots;
        m_scalingFactorBFP = scalingFactorBFP;
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
        return m_zSlots * m_zN;
    }
    BigFixedPoint GetScalingFactorBFP() const {
        return m_scalingFactorBFP;
    }
    uint32_t GetZN() const {
        return m_zN;
    }
    uint32_t GetZSlots() const {
        return m_zSlots;
    }

    static ZEncoding encodeR(const RPolynomial& input, const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                             const BigFixedPoint& scalingFactor) {
        // The big one
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
        return std::make_shared<ZEncodingImpl>(elementParams, poly, input.getZN(), input.getZSlots(), scalingFactor);
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
        return {input->GetZN(), input->GetZSlots(), output};
    }
};

}  // namespace lbcrypto

#endif  // SRC_PKE_EXAMPLES_Z_ENCODE_H_