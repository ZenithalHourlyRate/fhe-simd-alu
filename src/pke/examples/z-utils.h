#include <cassert>
#include "openfhe.h"
#include "encoding/z-encode.h"

using namespace lbcrypto;
using CiphertextT        = ConstCiphertext<DCRTPoly>;
using MutableCiphertextT = Ciphertext<DCRTPoly>;
using CCParamsT          = CCParams<CryptoContextBFVRNS>;
using CryptoContextT     = CryptoContext<DCRTPoly>;
using EvalKeyT           = EvalKey<DCRTPoly>;
using PlaintextT         = Plaintext;
using PrivateKeyT        = PrivateKey<DCRTPoly>;
using PublicKeyT         = PublicKey<DCRTPoly>;

//=============================================================================
// Encryption Utils
//=============================================================================

// DecryptCore not accessible from CryptoContext
// so copy from @openfhe//src/pke/lib/schemerns/rns-pke.cpp
DCRTPoly DecryptCore(const std::vector<DCRTPoly>& cv, const PrivateKey<DCRTPoly> privateKey) {
    const DCRTPoly& s = privateKey->GetPrivateElement();

    size_t sizeQ  = s.GetParams()->GetParams().size();
    size_t sizeQl = cv[0].GetParams()->GetParams().size();

    size_t diffQl = sizeQ - sizeQl;

    auto scopy(s);
    scopy.DropLastElements(diffQl);

    DCRTPoly sPower(scopy);

    DCRTPoly b(cv[0]);
    b.SetFormat(Format::EVALUATION);

    DCRTPoly ci;
    for (size_t i = 1; i < cv.size(); i++) {
        ci = cv[i];
        ci.SetFormat(Format::EVALUATION);

        b += sPower * ci;
        sPower *= scopy;
    }
    return b;
}

std::shared_ptr<std::vector<DCRTPoly>> EncryptZeroCore(const PublicKey<DCRTPoly> publicKey) {
    const auto cryptoParams =
        std::dynamic_pointer_cast<CryptoParametersRLWE<DCRTPoly>>(publicKey->GetCryptoParameters());

    const auto ns = cryptoParams->GetNoiseScale();

    auto elementParams = cryptoParams->GetElementParams();

    auto pk = publicKey->GetPublicElements();

    auto p0 = pk[0];
    auto p1 = pk[1];

    uint32_t sizeQ  = elementParams->GetParams().size();
    uint32_t sizePK = p0.GetParams()->GetParams().size();

    if (sizePK > sizeQ) {
        p0.DropLastElements(sizePK - sizeQ);
        p1.DropLastElements(sizePK - sizeQ);
    }

    DCRTPoly::TugType tug;
    DCRTPoly v = DCRTPoly(tug, elementParams, Format::EVALUATION);

    // noise generation with the discrete gaussian generator dgg
    auto& dgg = cryptoParams->GetDiscreteGaussianGenerator();
    DCRTPoly e0(dgg, elementParams, Format::EVALUATION);
    DCRTPoly e1(dgg, elementParams, Format::EVALUATION);

    DCRTPoly b(elementParams);
    DCRTPoly a(elementParams);

    b = p0 * v + ns * e0;
    a = p1 * v + ns * e1;

    return std::make_shared<std::vector<DCRTPoly>>(std::initializer_list<DCRTPoly>({std::move(b), std::move(a)}));
}

std::shared_ptr<std::vector<DCRTPoly>> EncryptZeroCore(const PrivateKey<DCRTPoly> privateKey) {
    const auto cryptoParams =
        std::dynamic_pointer_cast<CryptoParametersRLWE<DCRTPoly>>(privateKey->GetCryptoParameters());
    const auto elementParams = cryptoParams->GetElementParams();

    DCRTPoly::DugType dug;
    DCRTPoly a(dug, elementParams, Format::EVALUATION);

    DCRTPoly e(cryptoParams->GetDiscreteGaussianGenerator(), elementParams, Format::EVALUATION);
    NativeInteger ns = cryptoParams->GetNoiseScale();

    // {b = ns * e - a * s, a}
    DCRTPoly b(std::move((e *= ns) -= (a * privateKey->GetPrivateElement())));

    return std::make_shared<std::vector<DCRTPoly>>(std::initializer_list<DCRTPoly>({std::move(b), std::move(a)}));
}

Ciphertext<DCRTPoly> Encrypt(Plaintext ptxt, const PublicKey<DCRTPoly> publicKey) {
    ZEncoding zEnc    = std::dynamic_pointer_cast<ZEncodingImpl>(ptxt);
    auto zEncDCRTPoly = zEnc->GetElement<DCRTPoly>();

    auto ba = EncryptZeroCore(publicKey);
    (*ba)[0] += zEncDCRTPoly;

    auto ctxt = std::make_shared<CiphertextImpl<DCRTPoly>>(publicKey);
    ctxt->SetElements(std::move(*ba));
    ctxt->SetNoiseScaleDeg(1);
    ctxt->SetScalingFactorBFP(zEnc->GetScalingFactorBFP());
    return ctxt;
}

Ciphertext<DCRTPoly> EncryptZero(const PublicKey<DCRTPoly> publicKey) {
    auto ba   = EncryptZeroCore(publicKey);
    auto ctxt = std::make_shared<CiphertextImpl<DCRTPoly>>(publicKey);
    ctxt->SetElements(std::move(*ba));
    return ctxt;
}

Ciphertext<DCRTPoly> Conjugate(ConstCiphertext<DCRTPoly> ciphertext,
                               const std::map<uint32_t, EvalKey<DCRTPoly>>& evalKeyMap) {
    uint32_t N = ciphertext->GetElements()[0].GetRingDimension();
    std::vector<uint32_t> vec(N);
    PrecomputeAutoMap(N, 2 * N - 1, &vec);

    auto result = ciphertext->Clone();

    auto algo = ciphertext->GetCryptoContext()->GetScheme();
    algo->KeySwitchInPlace(result, evalKeyMap.at(2 * N - 1));

    auto& rcv = result->GetElements();
    rcv[0]    = rcv[0].AutomorphismTransform(2 * N - 1, vec);
    rcv[1]    = rcv[1].AutomorphismTransform(2 * N - 1, vec);
    return result;
}

//=============================================================================
// Custom Evals
//=============================================================================

Ciphertext<DCRTPoly> EvalAdd(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt) {
    assert(ct->GetElements().size() == ptxt.GetParams()->GetParams().size() &&
           "Ciphertext and Plaintext size mismatch in EvalAdd");
    assert(ptxt.GetFormat() == Format::EVALUATION && "Plaintext must be in EVALUATION format in EvalMultDCRTPoly");
    ZEncoding zEnc    = std::dynamic_pointer_cast<ZEncodingImpl>(ptxt);
    auto zEncDCRTPoly = zEnc->GetElement<DCRTPoly>();

    auto ctNew = ct->Clone();
    auto& b    = ctNew->GetElements()[0];
    b += zEncDCRTPoly;
    return ctNew;
}

Ciphertext<DCRTPoly> EvalMult(ConstCiphertext<DCRTPoly> ct, Plaintext ptxt) {
    ZEncoding zEnc    = std::dynamic_pointer_cast<ZEncodingImpl>(ptxt);
    auto zEncDCRTPoly = zEnc->GetElement<DCRTPoly>();
    assert(ct->GetElements().size() == zEncDCRTPoly.GetParams()->GetParams().size() &&
           "Ciphertext and Plaintext size mismatch in EvalMult");
    assert(zEncDCRTPoly.GetFormat() == Format::EVALUATION &&
           "Plaintext must be in EVALUATION format in EvalMultDCRTPoly");

    auto ctNew = ct->Clone();
    for (auto& a : ctNew->GetElements()) {
        a *= zEncDCRTPoly;
    }
    return ctNew;
}

Ciphertext<DCRTPoly> EvalMultScalar(ConstCiphertext<DCRTPoly> ct, BigInteger scalar) {
    assert(ct->GetElements().size() == ptxt.GetParams()->GetParams().size() &&
           "Ciphertext and Plaintext size mismatch in EvalMultDCRTPoly");
    assert(ptxt.GetFormat() == Format::EVALUATION && "Plaintext must be in EVALUATION format in EvalMultDCRTPoly");
    auto ctNew = ct->Clone();
    for (auto& a : ctNew->GetElements()) {
        a *= scalar;
    }
    return ctNew;
}

void ModReduceCustomInPlace(Ciphertext<DCRTPoly>& ciphertext, size_t levels = 1) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());

    auto& cv = ciphertext->GetElements();

    size_t sizeQ  = cryptoParams->GetElementParams()->GetParams().size();
    size_t sizeQl = cv[0].GetNumOfElements();
    size_t diffQl = sizeQ - sizeQl;

    ciphertext->SetLevel(ciphertext->GetLevel() + levels);

    for (size_t i = 0; i < levels; ++i) {
        for (auto& dcrtpoly : cv)
            dcrtpoly.DropLastElementAndScale(cryptoParams->GetQlQlInvModqlDivqlModq(diffQl + i),
                                             cryptoParams->GetqlInvModq(diffQl + i));
        // We manually track scaling factor
        auto ql      = cryptoParams->GetElementParams()->GetParams()[sizeQl - 1 - i]->GetModulus();
        auto qlBigFP = BigFixedPoint(ql, 0, false).scaleTo(128);
        ciphertext->SetScalingFactorBFP(ciphertext->GetScalingFactorBFP() / qlBigFP);

        // Old code from OpenFHE
        // double modReduceFactor = cryptoParams->GetModReduceFactor(sizeQl - 1 - i);
        // ciphertext->SetScalingFactor(ciphertext->GetScalingFactor() / modReduceFactor);
    }
}

//=============================================================================
// KeyGen Related
//=============================================================================

std::shared_ptr<std::map<uint32_t, EvalKey<DCRTPoly>>> EvalSparseEncapsulatedKeyGen(
    const PrivateKey<DCRTPoly> privateKey) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(privateKey->GetCryptoParameters());

    auto cc   = privateKey->GetCryptoContext();
    auto algo = cc->GetScheme();
    auto M    = cc->GetCyclotomicOrder();

    // computing all indices for baby-step giant-step procedure
    auto evalKeys = std::make_shared<std::map<uint32_t, EvalKey<DCRTPoly>>>();

    if (cryptoParams->GetSecretKeyDist() == SPARSE_ENCAPSULATED) {
        DCRTPoly::TugType tug;

        // sparse key used for the modraising step
        auto skNew = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
        skNew->SetPrivateElement(DCRTPoly(tug, cryptoParams->GetElementParams(), Format::EVALUATION, 32));

        // we reserve M-4 and M-2 for the sparse encapsulation switching keys
        // Even autorphism indices are not possible, so there will not be any conflict
        (*evalKeys)[M - 4] = lbcrypto::FHECKKSRNS::KeySwitchGenSparse(privateKey, skNew);
        (*evalKeys)[M - 2] = algo->KeySwitchGen(skNew, privateKey);
    }

    cc->InsertEvalAutomorphismKey(evalKeys, privateKey->GetKeyTag());
    return evalKeys;
}

//=============================================================================
// Linear Transform Auxiliaries
//=============================================================================

std::array<std::vector<Plaintext>, 2> getAuxLTPtxt(BigCMatrix T, bool inverse,
                                                   const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                                   BigFixedPoint scalingFactor) {
    // Tinverse is of shape zN * (zN / 2)
    // T of is of (zN / 2) * zN
    auto halfSize = T.size();
    if (inverse) {
        halfSize = T[0].size();
    }
    auto rowOffset    = 0;
    auto columnOffset = 0;
    if (inverse) {
        rowOffset = halfSize;
    }
    else {
        columnOffset = halfSize;
    }

    std::vector<Plaintext> results1;
    std::vector<Plaintext> results2;

    for (size_t i = 0; i != halfSize; ++i) {
        auto diagonal = std::vector<BigComplex>(halfSize, BigComplex());
        for (size_t j = 0; j != halfSize; ++j) {
            diagonal[j] = T[j][(i + j) % halfSize];
        }
        auto diagonalInR = CSlots(diagonal).toRPolynomial();
        Plaintext ptxt   = ZEncodingImpl::encodeR(diagonalInR.getCoefficients(), elementParams, scalingFactor);
        results1.push_back(ptxt);
    }
    for (size_t i = 0; i != halfSize; ++i) {
        auto diagonal = std::vector<BigComplex>(halfSize, BigComplex());
        for (size_t j = 0; j != halfSize; ++j) {
            diagonal[j] = T[j + rowOffset][(i + j) % halfSize + columnOffset];
        }
        auto diagonalInR = CSlots(diagonal).toRPolynomial();
        Plaintext ptxt   = ZEncodingImpl::encodeR(diagonalInR.getCoefficients(), elementParams, scalingFactor);
        results2.push_back(ptxt);
    }
    return {results1, results2};
}

std::array<std::vector<Plaintext>, 2> getZCoeffToSlotsAuxDCRTPoly(
    const std::shared_ptr<typename DCRTPoly::Params>& elementParams, BigFixedPoint scalingFactor) {
    return getAuxLTPtxt(getZUInverse(), true, elementParams, scalingFactor);
}

std::array<std::vector<Plaintext>, 2> getRCoeffToSlotsAuxDCRTPoly(
    const std::shared_ptr<typename DCRTPoly::Params>& elementParams, BigFixedPoint scalingFactor) {
    return getAuxLTPtxt(getRUInverse(), true, elementParams, scalingFactor);
}

std::array<std::vector<Plaintext>, 2> getSlotsToZCoeffsAuxDCRTPoly(
    const std::shared_ptr<typename DCRTPoly::Params>& elementParams, BigFixedPoint scalingFactor) {
    return getAuxLTPtxt(getZU(), false, elementParams, scalingFactor);
}

std::array<std::vector<Plaintext>, 2> getSlotsToRCoeffsAuxDCRTPoly(
    const std::shared_ptr<typename DCRTPoly::Params>& elementParams, BigFixedPoint scalingFactor) {
    return getAuxLTPtxt(getRU(), false, elementParams, scalingFactor);
}

//=============================================================================
// Linear Transform Evals
//=============================================================================

// We now do not exploit sparse encoding

std::vector<Ciphertext<DCRTPoly>> CoeffsToSlots(CryptoContextT cc, CiphertextT ct,
                                                const std::array<std::vector<Plaintext>, 2>& auxPtxts) {
    // Halevi-Shoup
    // Peel the first loop
    auto resultUpper = EvalMult(ct, auxPtxts[0][0]);
    auto resultDown  = EvalMult(ct, auxPtxts[1][0]);
    auto startCt     = cc->EvalRotate(ct, 1);
    for (size_t i = 1; i != auxPtxts[0].size(); ++i) {
        auto diagonalUpper = EvalMult(startCt, auxPtxts[0][i]);
        cc->EvalAddInPlace(resultUpper, diagonalUpper);
        auto diagonalDown = EvalMult(startCt, auxPtxts[1][i]);
        cc->EvalAddInPlace(resultDown, diagonalDown);
        if (i + 1 < auxPtxts[0].size()) {
            //rotate one more
            startCt = cc->EvalRotate(ct, i + 1);
        }
    }
    cc->EvalAddInPlace(resultUpper, Conjugate(resultUpper, cc->GetEvalAutomorphismKeyMap(resultUpper->GetKeyTag())));
    cc->EvalAddInPlace(resultDown, Conjugate(resultDown, cc->GetEvalAutomorphismKeyMap(resultDown->GetKeyTag())));

    return {resultUpper, resultDown};
}

Ciphertext<DCRTPoly> SlotsToCoeffs(CryptoContextT cc, CiphertextT ctLeft, CiphertextT ctRight,
                                   const std::array<std::vector<Plaintext>, 2>& auxPtxts) {
    // Halevi-Shoup
    // Peel the first loop
    auto result = EvalMult(ctLeft, auxPtxts[0][0]);
    cc->EvalAddInPlace(result, EvalMult(ctRight, auxPtxts[1][0]));
    auto startCtLeft  = cc->EvalRotate(ctLeft, 1);
    auto startCtRight = cc->EvalRotate(ctRight, 1);
    for (size_t i = 1; i != auxPtxts[0].size(); ++i) {
        auto diagonalLeft = EvalMult(startCtLeft, auxPtxts[0][i]);
        cc->EvalAddInPlace(result, diagonalLeft);
        auto diagonalRight = EvalMult(startCtRight, auxPtxts[1][i]);
        cc->EvalAddInPlace(result, diagonalRight);
        // rotate one more
        if (i + 1 < auxPtxts[0].size()) {
            startCtLeft  = cc->EvalRotate(startCtLeft, 1);
            startCtRight = cc->EvalRotate(startCtRight, 1);
        }
    }
    return result;
}

// NOTE: should pre-compute ptxts outside
std::vector<Ciphertext<DCRTPoly>> ZCoeffsToSlots(CryptoContextT cc, CiphertextT ct,
                                                 BigFixedPoint sf = BigFixedPoint::zero()) {
    auto elementParams = ct->GetElements()[0].GetParams();
    if (sf.equalZero()) {
        sf = ct->GetScalingFactorBFP();
    }
    return CoeffsToSlots(cc, ct, getZCoeffToSlotsAuxDCRTPoly(elementParams, sf));
}

std::vector<Ciphertext<DCRTPoly>> RCoeffsToSlots(CryptoContextT cc, CiphertextT ct,
                                                 BigFixedPoint sf = BigFixedPoint::zero()) {
    auto elementParams = ct->GetElements()[0].GetParams();
    if (sf.equalZero()) {
        sf = ct->GetScalingFactorBFP();
    }
    return CoeffsToSlots(cc, ct, getRCoeffToSlotsAuxDCRTPoly(elementParams, sf));
}

Ciphertext<DCRTPoly> SlotsToZCoeffs(CryptoContextT cc, CiphertextT ctLeft, CiphertextT ctRight,
                                    BigFixedPoint sf = BigFixedPoint::zero()) {
    auto elementParams = ctLeft->GetElements()[0].GetParams();
    if (sf.equalZero()) {
        sf = ctLeft->GetScalingFactorBFP();
    }
    return SlotsToCoeffs(cc, ctLeft, ctRight, getSlotsToZCoeffsAuxDCRTPoly(elementParams, sf));
}

Ciphertext<DCRTPoly> SlotsToRCoeffs(CryptoContextT cc, CiphertextT ctLeft, CiphertextT ctRight,
                                    BigFixedPoint sf = BigFixedPoint::zero()) {
    auto elementParams = ctLeft->GetElements()[0].GetParams();
    if (sf.equalZero()) {
        sf = ctLeft->GetScalingFactorBFP();
    }
    return SlotsToCoeffs(cc, ctLeft, ctRight, getSlotsToRCoeffsAuxDCRTPoly(elementParams, sf));
}

//=============================================================================
// Pre-defined Chebyshev Series Coefficients
//=============================================================================

// Coefficients for the function std::exp(1i * Pi/2.0 * x) in [-16, 16] of degree 46
// Need two double-angle iterations to get std::exp(1i * 2Pi * x)
static const inline std::vector<std::complex<double>> coeff_exp_16_double_46{
    0.22393566906777329084,
    4.72e-18 - 0.22176384914036376128i,
    0.24158307546266130639 - 7.09e-18i,
    -0.18331470851313935722i,
    0.28534623846463541552 + 5.91e-18i,
    -2.362e-18 - 0.092486179824488604084i,
    0.32214532018151836867 - 4.72e-18i,
    0.061326880477941263237i,
    0.28798365357787297780 - 4.72e-18i,
    1.417e-17 + 0.24466296846427090794i,
    0.112756709876059055264 - 7.087e-18i,
    0.33439190718203853914i,
    -0.17995397739265364678 + 2.36e-18i,
    9.45e-18 + 0.16254851699551037258i,
    -0.34811157721125479680 - 2.36e-18i,
    2.95e-18 - 0.22527723082929962395i,
    -0.079206690817228031509 - 3.543e-18i,
    9.45e-18 - 0.32612632178540534866i,
    0.36198254675123692214 - 4.13e-18i,
    -4.72e-18 + 0.19237548287066727482i,
    0.071116210979946081761 - 9.449e-18i,
    -7.09e-18 + 0.30556044798491310832i,
    -0.43951407397686892420 + 2.36e-18i,
    4.72e-18 - 0.46389876376571992367i,
    0.40955141151976887093 - 3.54e-18i,
    0.31828681535788988510i,
    -0.22366008829505132360 - 4.72e-18i,
    4.72e-18 - 0.14446909676096356123i,
    0.086745018497585937856 - 7.087e-18i,
    2.362e-18 + 0.048813481993878263254i,
    -0.025904132260782003483 - 2.362e-18i,
    2.3622e-18 - 0.0130280784432669702322i,
    0.0062348555293589751417 - 7.0865e-18i,
    -1.41731e-17 + 0.0028488507881147717878i,
    -0.00124638777412581667689 + 6.49599e-18i,
    -1.77163e-18 - 0.00052341839132928819674i,
    0.00021144315086689963591 - 2.36218e-18i,
    -4.724353e-18 + 0.000082321616249500458072i,
    -0.000030941853907365078754,
    -9.4487066e-18 - 0.0000112448146643207249387i,
    3.9566691190946415917e-6 + 4.7243533e-18i,
    -4.72435330e-18 + 1.34965353351963424859e-6i,
    -4.4681665430474632002e-7 - 7.08652994e-18i,
    3.54326497e-18 - 1.4370869470312489006e-7i,
    4.4978580236218894379e-8 + 2.362176648e-18i,
    -4.7243532963e-18 + 1.35960019696524527036e-8i,
    -4.3910916203385461663e-9 - 4.7243532963e-18i};
