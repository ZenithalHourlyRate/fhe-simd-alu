#include "openfhe.h"
#include "utils.h"
#include "scheme/ckksrns/z-fhe.h"

using namespace lbcrypto;

// related to __heir_debug2
CryptoContextT cc_global;
PublicKeyT pk_global;
PrivateKeyT sk_global;
size_t zN_global;
size_t zSlots_global;

std::vector<BigComplex> zC2SVals;

double __heir_debug2(CiphertextT ct, std::string msg) {
    auto b = DecryptCore(ct->GetElements(), sk_global);

    auto sfBigFP = ct->GetScalingFactorBFP();
    //std::cout << msg << "  Ciphertext Scaling Factor BFP: " << sfBigFP.toHexString() << std::endl;
    auto log2sf = std::log2(sfBigFP.convertToDouble());
    std::cout << msg << "  Scaling factor log2: " << std::setprecision(20) << log2sf << std::endl;
    auto q = ct->GetElements()[0].GetParams()->GetModulus();
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
    auto l = ct->GetElements()[0].GetParams()->GetParams().size();
    std::cout << msg << "  l: " << l - 1 << std::endl;

    ZEncodingParams params(ZMode, zN_global, zSlots_global);
    auto zEncode = std::make_shared<ZEncodingImpl>(b.GetParams(), b, sfBigFP, params);

    RPolynomial values = ZEncodingImpl::decodeR(zEncode);

    // This is for sparse LT
    ZEncodingParams paramsTwice(CMode, zN_global * zSlots_global * 2);
    auto zEncodeTwice       = std::make_shared<ZEncodingImpl>(b.GetParams(), b, sfBigFP, paramsTwice);
    RPolynomial valuesTwice = ZEncodingImpl::decodeR(zEncodeTwice);

    auto printZPoly = [&](const ZPolynomial zPoly) {
        auto decoded = ZPolynomial::decode(zPoly);
        std::cout << msg << "  zPoly Decoded: " << std::hex << decoded << std::dec
                  << "  zPoly error log2Norm: " << ZPolynomial::extractError(zPoly).getLog2Norm() << std::endl;
        //auto I = ZPolynomial::extractI(zPoly);
        //std::cout << msg << "  zPoly I: ";
        //for (size_t i = 0; i != I.getCoefficients().size(); ++i) {
        //    std::cout << I[i].toHexString(16) << " ";
        //}
        //std::cout << std::endl;
        //std::cout << msg << "  zPoly error log2Norm: " << ZPolynomial::extractError(zPoly).getLog2Norm() << std::endl;
    };

    enum class DecodeMode { RDecode, CSlotsDecode, CSlotsTwiceDecode, CSlotsTwiceBooleanMode, ZDecode };
    std::map<std::string, DecodeMode> decodeMap = {
        // RDecode
        {"ModRaise", DecodeMode::RDecode},
        {"C2R", DecodeMode::RDecode},
        {"Z2R", DecodeMode::RDecode},
        // CSlotsDecode
        {"CSlotsDecode", DecodeMode::CSlotsDecode},
        // CSlotsTwiceDecode
        {"LUT", DecodeMode::CSlotsTwiceDecode},
        {"Normalize", DecodeMode::CSlotsTwiceDecode},
        {"Core", DecodeMode::CSlotsTwiceDecode},
        {"Z2C", DecodeMode::CSlotsTwiceDecode},
        {"Z2C2", DecodeMode::CSlotsTwiceDecode},
        {"Boolean", DecodeMode::CSlotsTwiceBooleanMode},
        {"BooleanAgain", DecodeMode::CSlotsTwiceBooleanMode},
        // ZDecode
        {"Input", DecodeMode::ZDecode},
        {"CMult", DecodeMode::ZDecode},
        {"NoiseT", DecodeMode::ZDecode},
    };

    auto decodeModeIt = decodeMap.find(msg);
    if (decodeModeIt == decodeMap.end()) {
        return 0;
    }
    auto decodeMode = decodeModeIt->second;

    // Check the r encoding?
    if (decodeMode == DecodeMode::RDecode) {
        // Threshold print to avoid too much output
        for (size_t i = 0; i != std::min(values.getCoefficients().size(), 32ul); ++i) {
            std::cout << msg << "  values [" << i << "]: " << values[i].toHexString(ceil(log2sf / 4.0)) << std::endl;
        }
        if (values.getCoefficients().size() > 32) {
            std::cout << msg << "  ... (total " << values.getCoefficients().size() << " coefficients)" << std::endl;
        }
    }
    if (decodeMode == DecodeMode::CSlotsTwiceDecode || decodeMode == DecodeMode::CSlotsTwiceBooleanMode) {
        auto cSlotsTwice = valuesTwice.toCSlots();

        if (decodeMode == DecodeMode::CSlotsTwiceBooleanMode) {
            // Force re-interpretation at boolean mode
            CSlots cSlotsTwiceBooleanMode(params, cSlotsTwice.getSlots());
            auto maxSlotsToPrint = std::min(zSlots_global, size_t(8));
            for (size_t i = 0; i != maxSlotsToPrint; ++i) {
                auto [integerValue, log2Error] = cSlotsTwiceBooleanMode.getIntegerAndErrorAtBooleanMode(i);
                std::cout << msg << "  Boolean Mode Slot " << i << " Reconstructed: " << std::hex << integerValue
                          << std::dec << " with error log2: " << log2Error << std::endl;
            }
            if (zSlots_global > maxSlotsToPrint) {
                std::cout << msg << "  ... (total " << zSlots_global << " slots)" << std::endl;
            }
        }
        else {
            //auto exampleSlotIndex = zSlots_global - 1;
            auto exampleSlotIndex = 0;
            for (size_t i = 0; i != zN_global; ++i) {
                auto slotIndex = exampleSlotIndex;
                auto index     = slotIndex * (zN_global / 2) + i;
                if (i >= zN_global / 2) {
                    index += (zSlots_global - 1) * zN_global / 2;
                }
                auto p        = BigFixedPoint::positive(1 << zN_global);
                auto lutPart  = (cSlotsTwice[index].getReal() * p).round() / p;
                auto fracPart = cSlotsTwice[index].getReal() - lutPart;
                std::cout << msg << "  cSlotsTwice Slot Example " << exampleSlotIndex << " " << index << " : "
                          << lutPart.toHexString() << " " << fracPart.toHexString(ceil(log2sf / 4.0)) << std::endl;
            }
        }
    }
    if (decodeMode == DecodeMode::ZDecode) {
        auto cSlots = values.toCSlots();
        std::vector<ZPolynomial> zPolys(zSlots_global, ZPolynomial(zN_global));
        auto maxSlotsToPrint = std::min(zSlots_global, size_t(8));
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(maxSlotsToPrint))
        for (size_t i = 0; i != maxSlotsToPrint; ++i) {
            zPolys[i] = cSlots.getZPolynomial(i);
        }
        for (size_t i = 0; i != maxSlotsToPrint; ++i) {
            printZPoly(zPolys[i]);
        }
        if (zSlots_global > maxSlotsToPrint) {
            std::cout << msg << "  ... (total " << zSlots_global << " slots)" << std::endl;
        }

        //ZPolynomial zValues = values.toCSlots().getZPolynomial(0);
        //for (size_t i = 0; i != zValues.getCoefficients().size(); ++i) {
        //    std::cout << msg << "  zValues [" << i << "]: " << zValues[i].toHexString(ceil(log2sf / 4.0)) << std::endl;
        //}
    }
    return 0;
}

void SimpleBootstrapExample();

int main(int argc, char* argv[]) {
    SimpleBootstrapExample();
}

void SimpleBootstrapExample() {
    CCParams<CryptoContextCKKSRNS> parameters;

    SecretKeyDist secretKeyDist = lbcrypto::SPARSE_ENCAPSULATED;
    parameters.SetSecretKeyDist(secretKeyDist);

    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 12);
    //parameters.SetNumLargeDigits(6);

    ScalingTechnique rescaleTech = FLEXIBLEMANUAL;
    uint32_t dcrtBits            = 45;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);

    std::vector<uint32_t> levelBudget = {2, 2};

    parameters.SetMultiplicativeDepth(20);

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

    uint32_t zN     = 16;
    uint32_t zSlots = cc->GetRingDimension() / zN / 2;  // Maximal sparse packing
    //uint32_t zSlots = 32;
    zN_global     = zN;
    zSlots_global = zSlots;
    std::cout << "Bootstrapping parameters: zN = " << zN << ", zSlots = " << zSlots << std::endl;

    fheZ->EvalBootstrapSetup(*cc, zN, zSlots, levelBudget, {0, 0}, 4, -16);
    fheZ->EvalBootstrapKeyGen(keyPair.secretKey, zN, zSlots);

    cc_global = cc;
    pk_global = keyPair.publicKey;
    sk_global = keyPair.secretKey;

    auto elemParam = cc->GetCryptoParameters()->GetElementParams();
    auto sfq0      = cryptoParams->GetScalingFactorBFP(0);

    std::vector<uint64_t> vec(zSlots, 0);
    for (size_t i = 0; i != zSlots; ++i) {
        vec[i] = i + 1;
    }
    Plaintext ptxt1 = ZEncodingImpl::encodeArith(vec, zN, zSlots, elemParam, sfq0);

    /// TEST ENCODE
    auto encoded = Encrypt(ptxt1, keyPair.publicKey);

    __heir_debug2(encoded, "Input");

    Ciphertext<DCRTPoly> ct = encoded;

    //Ciphertext<DCRTPoly> ct2 = fheZ->EvalArithToBoolean(ct);
    //__heir_debug2(ct2, "Boolean");
    //BENCHMARK(fheZ->EvalArithToBoolean(ct), 3, "ArithToBoolean");

    Ciphertext<DCRTPoly> ct3 = fheZ->EvalArithToArith(ct);
    __heir_debug2(ct3, "CMult");
    //BENCHMARK(fheZ->EvalArithToArith(ct), 3, "ArithToArith");

    //BENCHMARK(fheZ->EvalBooleanToBoolean(ct2), 3, "BooleanToBoolean");
    //__heir_debug2(ct, "BooleanAgain");
}
