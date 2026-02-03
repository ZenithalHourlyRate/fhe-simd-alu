
#include "openfhe.h"
#include "utils.h"
#include "scheme/ckksrns/z-fhe.h"

using namespace lbcrypto;

double __heir_debug3(lbcrypto::Plaintext ptxt, std::string msg) {
    auto b       = ptxt->GetElement<DCRTPoly>();
    auto zEncode = std::dynamic_pointer_cast<ZEncodingImpl>(ptxt);

    auto sfBigFP = zEncode->GetScalingFactorBFP();
    //std::cout << msg << "  Ciphertext Scaling Factor BFP: " << sfBigFP.toHexString() << std::endl;
    auto log2sf = std::log2(sfBigFP.convertToDouble());
    std::cout << msg << "  Scaling factor log2: " << std::setprecision(20) << log2sf << std::endl;
    auto q = b.GetParams()->GetModulus();
    double log2q;
    if (q.GetMSB() > 128) {
        auto offset = q.GetMSB() - 128;
        q >>= offset;
        log2q = std::log2(q.ConvertToDouble()) + offset;
    }
    else {
        log2q = std::log2(q.ConvertToDouble());
    }
    std::cout << msg << "  q: " << log2q << std::endl;
    auto l = b.GetParams()->GetParams().size();
    std::cout << msg << "  l: " << l - 1 << std::endl;

    RPolynomial values   = ZEncodingImpl::decodeR(zEncode);
    auto cSlots          = values.toCSlots();
    auto maxSlotsToPrint = std::min(zSlots_global, size_t(2));
    for (size_t i = 0; i != maxSlotsToPrint; ++i) {
        auto value    = cSlots[i].getReal();
        auto p        = BigFixedPoint::positive(1 << zN_global);
        auto lutPart  = (value * p).round() / p;
        auto fracPart = value - lutPart;
        //double log2Error = std::log2(std::abs(fracPart.convertToDouble()));
        //std::cout << msg << "  cSlots Slot " << i << " Reconstructed: " << lutPart.toHexString() << " "
        //          << fracPart.toHexString(ceil(log2sf / 4.0)) << " error: " << std::setprecision(2) << log2Error
        //          << std::endl;
        std::cout << msg << "  cSlots Slot " << i << " " << cSlots[i].getReal().toHexString(ceil(log2sf / 4.0))
                  << std::endl;
    }
    if (zSlots_global > maxSlotsToPrint) {
        std::cout << msg << "  ... (total " << zSlots_global << " slots)" << std::endl;
    }
    return 0.0;
}
void SimpleBootstrapExample2(int zN) {
    CCParams<CryptoContextCKKSRNS> parameters;

    SecretKeyDist secretKeyDist = lbcrypto::SPARSE_ENCAPSULATED;
    parameters.SetSecretKeyDist(secretKeyDist);

    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 9);
    //parameters.SetNumLargeDigits(6);

    ScalingTechnique rescaleTech = FLEXIBLEMANUAL;
    uint32_t dcrtBits            = 30;
    // bit size for aux moduli in P
    // for HEXL acceleration. Extra 6 bit for SPARSE_ENCAPSULATED
    AUXMODSIZE_FLEXIBLEMANUAL = 45;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetNumLargeDigits(3);

    std::vector<uint32_t> levelBudget = {1, 1};

    uint32_t mulDepth = 18;
    parameters.SetMultiplicativeDepth(mulDepth);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    uint32_t ringDim = cc->GetRingDimension();
    std::cout << "CKKS scheme ring dimension: " << ringDim << "\n\n";

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    //std::cout << *cryptoParams << std::endl;
    //std::cout << *(cryptoParams->GetParamsP()) << " primes in the special prime modulus." << std::endl;
    double logQ = 0;
    double logP = 0;
    {
        auto moduliQ = cc->GetCryptoParameters()->GetElementParams()->GetModulus();
        auto moduliP = cryptoParams->GetParamsP()->GetModulus();
        logQ         = moduliQ.GetMSB();
        logP         = moduliP.GetMSB();
    }
    std::cout << "log2(Q) = " << logQ << " log2(P) = " << logP << " log2(QP) = " << logQ + logP << std::endl;

    LeveledZ z     = std::make_shared<LeveledZImpl>();
    AdvancedZ advZ = std::make_shared<AdvancedZImpl>(z);
    FHEZ fheZ      = std::make_shared<FHEZImpl>(z, advZ);

    //uint32_t zN     = 16;
    //uint32_t zSlots = cc->GetRingDimension() / zN;  // Full packing
    uint32_t zSlots = 1;  // 1 Slot
    //uint32_t zSlots = cc->GetRingDimension() / zN / 2;  // Maximal sparse packing
    //uint32_t zSlots = 32;
    zN_global     = zN;
    zSlots_global = zSlots;
    std::cout << "Bootstrapping parameters: zN = " << zN << ", zSlots = " << zSlots << std::endl;

    fheZ->EvalBootstrapSetup(*cc, zN, zSlots, levelBudget, {0, 0}, 4, -16, 1, 1);
    fheZ->EvalBootstrapKeyGen(keyPair.secretKey, zN, zSlots);

    cc_global = cc;
    pk_global = keyPair.publicKey;
    sk_global = keyPair.secretKey;

    auto elemParam = cc->GetCryptoParameters()->GetElementParams();
    auto sfq0      = cryptoParams->GetScalingFactorBFP(0);
    //auto lArith         = cryptoParams->GetMultiplicativeDepth() - 10;
    //auto lArith = 0;
    //auto elemParamArith = cryptoParams->GetParamsQl(lArith);
    //auto sfArith        = cryptoParams->GetScalingFactorBFP(lArith);
    auto elemParamArith = elemParam;
    auto sfArith        = sfq0;

    //auto lBool         = cryptoParams->GetMultiplicativeDepth() - 10;
    //auto elemParamBool = cryptoParams->GetParamsQl(lBool);
    //auto sfBool = cryptoParams->GetScalingFactorBFP(lBool);
    auto elemParamBool = elemParam;
    auto sfBool        = sfq0;

    std::vector<uint64_t> vec(zSlots, 0);
    for (size_t i = 0; i != zSlots; ++i) {
        vec[i] = i + 3;
    }
    Plaintext ptxt1 = ZEncodingImpl::encodeArith(vec, zN, zSlots, elemParamArith, sfArith);

    std::vector<uint64_t> vec2(zSlots, 0);
    for (size_t i = 0; i != zSlots; ++i) {
        //vec2[i] = -i - 3;
        //vec2[i] = i + 0xdeadbeaf;
        vec2[i] = i + 16;
    }
    Plaintext ptxt2 = ZEncodingImpl::encodeArith(vec2, zN, zSlots, elemParamArith, sfArith);

    /// TEST ENCODE
    auto encoded  = Encrypt(ptxt1, keyPair.publicKey);
    auto encoded2 = Encrypt(ptxt2, keyPair.publicKey);

    //__heir_debug2(encoded, "Input");

    /// TEST Boolean Encode
    auto ptxt3    = ZEncodingImpl::encodeBooleanSparse(vec, zN, zSlots, elemParamBool, sfBool);
    auto encoded3 = Encrypt(ptxt3, keyPair.publicKey);

    auto ptxt4    = ZEncodingImpl::encodeBooleanSparse(vec2, zN, zSlots, elemParamBool, sfBool);
    auto encoded4 = Encrypt(ptxt4, keyPair.publicKey);

#define DEBUG

#ifdef DEBUG
    __heir_debug2(encoded, "Input");
#endif

    // Mult
    if (1) {
        BENCHMARK(z->EvalMultFullInZ(encoded2, encoded), 3, "MultFull");
#ifdef DEBUG
        auto ct2 = z->EvalMultFullInZ(encoded, encoded);
        __heir_debug2(ct2, "CMult");
#endif
    }

    // Bool
    if (1) {
        BENCHMARK(z->EvalBooleanOR(encoded3, encoded4), 3, "BooleanOR");
#ifdef DEBUG
        auto ct2 = z->EvalBooleanOR(encoded3, encoded4);
        __heir_debug2(ct2[0], "BoolOR");
#endif
    }

    // BoolToArith
    if (0) {
        BENCHMARK(fheZ->EvalBooleanToArith(encoded3), 1, "BooleanToArith");
#ifdef DEBUG
        auto ct2 = z->EvalBooleanOR(encoded3, encoded4);
        auto ct3 = fheZ->EvalBooleanToArith(ct2);
        __heir_debug2(ct3, "B2A");
#endif
    }

    // BooleanToBoolean
    if (1) {
        BENCHMARK((fheZ->EvalBooleanToBooleanSparse(encoded3)), 1, "BooleanToBoolean");
    }

    // ArithToArithHigh
    if (0) {
        BENCHMARK(fheZ->EvalArithToArithHigh(encoded2), 1, "ArithToArithHigh");
    }

    // ArithToArithNoise
    if (0) {
        BENCHMARK(fheZ->EvalArithToArithNoise(encoded2), 1, "ArithToArithNoise");
    }

    // ArithToArith
    if (1) {
        BENCHMARK(fheZ->EvalArithToArith(encoded2), 1, "ArithToArith");
#ifdef DEBUG
        auto ct2 = fheZ->EvalArithToArithHigh(encoded2);
        __heir_debug2(ct2, "A2A");

        auto ct3 = z->EvalMultFullInZ(ct2, ct2);
        z->ModReduceInPlace(ct3);
        __heir_debug2(ct3, "CMult");

        auto ct4 = fheZ->EvalArithToArith(ct3);
        __heir_debug2(ct4, "CMult");

        auto ct5 = fheZ->EvalArithToArith(ct4);
        __heir_debug2(ct5, "CMult");

        //auto ct5 = z->EvalMultFullInZ(ct4, ct4);
        //z->ModReduceInPlace(ct5);
        //__heir_debug2(ct5, "CMult");

        //auto ct6 = fheZ->EvalArithToArith(ct5);
        //__heir_debug2(ct6, "A2A");

        //auto ct7 = fheZ->EvalArithToArith(ct6);
        //__heir_debug2(ct7, "A2A");
#endif
    }

    // ArithToBooleanBatched
    if (0) {
        std::vector<Ciphertext<DCRTPoly>> batchCts = {encoded, encoded2};
        auto batchSize                             = zN / 4;  // zN / w
        while (batchCts.size() < batchSize) {
            batchCts.push_back(encoded);
        }
#ifdef DEBUG
        auto ctGroupBool = fheZ->EvalArithToBooleanBatched(batchCts);
        __heir_debug2(ctGroupBool[0], "A2B0");
        __heir_debug2(ctGroupBool[1], "A2B1");
        __heir_debug2(ctGroupBool[2], "A2B0");
        __heir_debug2(ctGroupBool[3], "A2B1");
        __heir_debug2(ctGroupBool[14], "A2B0");
        __heir_debug2(ctGroupBool[15], "A2B1");
#endif
        BENCHMARK(fheZ->EvalArithToBooleanBatched(batchCts), 1, "ArithToBooleanBatched");
    }

    //ArithToBooleanFull
    if (1) {
#ifdef DEBUG
        auto ctBool = fheZ->EvalArithToBooleanSparse(encoded2);
        __heir_debug2(ctBool, "A2B");
#endif
        BENCHMARK(fheZ->EvalArithToBooleanSparse(encoded2), 1, "ArithToBoolean");
    }
}