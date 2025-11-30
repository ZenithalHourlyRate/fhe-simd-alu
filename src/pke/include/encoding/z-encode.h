#ifndef SRC_PKE_EXAMPLES_Z_ENCODE_H_
#define SRC_PKE_EXAMPLES_Z_ENCODE_H_

#include <cassert>
#include "math/z-encode-utils.h"
#include "plaintext.h"

class ZEncodingImpl;

using ZEncoding = std::shared_ptr<ZEncodingImpl>;

// This is only interplay between Plaintext and RPolynomial
// RPolynomial to CSlots to ZPolynomial is handled elsewhere
class ZEncodingImpl : public PlaintextImpl {
private:
    BigFixedPoint m_scalingFactorBFP;
    size_t m_rPolySize;

public:
    ZEncodingImpl(std::shared_ptr<DCRTPoly::Params> vp, const DCRTPoly& elements, size_t rPolySize,
                  BigFixedPoint scalingFactorBFP)
        : PlaintextImpl(vp, nullptr, INVALID_ENCODING, INVALID_SCHEME) {
        encodedVectorDCRT  = elements;
        m_rPolySize        = rPolySize;
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
        return m_rPolySize;
    }
    BigFixedPoint GetScalingFactorBFP() const {
        return m_scalingFactorBFP;
    }

    static ZEncoding encodeR(RPolynomial input, const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                             BigFixedPoint scalingFactor) {
        DCRTPoly newPoly(elementParams, Format::COEFFICIENT, true);
        auto bigPoly = newPoly.CRTInterpolate();
        auto ringDim = elementParams->GetRingDimension();
        auto q       = elementParams->GetModulus();
        for (size_t i = 0; i != input.getCoefficients().size(); ++i) {
            auto val        = input[i] * scalingFactor;
            auto valInteger = (val.round()).getValue() >> val.getLog2Scale();
            if (val.getNeg()) {
                bigPoly[i * (ringDim / (input.getCoefficients().size()))] = q - valInteger;
            }
            else {
                bigPoly[i * (ringDim / (input.getCoefficients().size()))] = valInteger;
            }
        }
        DCRTPoly finalPoly(bigPoly, elementParams);
        finalPoly.SetFormat(Format::EVALUATION);
        return std::make_shared<ZEncodingImpl>(elementParams, finalPoly, input.getCoefficients().size(), scalingFactor);
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
            auto valFixedPoint = BigFixedPoint(valInteger, 0, neg).scaleTo(z_upper_roots_scale);
            output.push_back(valFixedPoint / scalingFactor);
        }
        return output;
    }
};

#endif  // SRC_PKE_EXAMPLES_Z_ENCODE_H_