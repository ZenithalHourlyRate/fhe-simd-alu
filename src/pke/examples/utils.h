#include <cassert>
#include "openfhe.h"
#include "encoding/z-encoding.h"
#include "scheme/ckksrns/z-leveledshe.h"
#include "scheme/ckksrns/z-advancedshe.h"

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
    ctxt->SetZEncodingParams(zEnc->GetZEncodingParams());
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