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
    std::cout << msg << "  zEncodingParams: " << ct->GetZEncodingParams().toString() << std::endl;

    // TODO: fix it
    auto zDeg = msg == "CMult" ? ct->GetZEncodingParams().getZDeg() : 1;
    ZEncodingParams params(ZMode, zN_global, zSlots_global, zDeg);
    auto zEncode = std::make_shared<ZEncodingImpl>(b.GetParams(), b, sfBigFP, params);

    RPolynomial values = ZEncodingImpl::decodeR(zEncode);

    // This is for sparse LT
    ZEncodingParams paramsTwice(CMode, zN_global * zSlots_global * 2);
    auto zEncodeTwice       = std::make_shared<ZEncodingImpl>(b.GetParams(), b, sfBigFP, paramsTwice);
    RPolynomial valuesTwice = ZEncodingImpl::decodeR(zEncodeTwice);

    auto printZPoly = [&](const ZPolynomial zPoly) {
        auto decoded = ZPolynomial::decode(zPoly);
        //for (size_t i = 0; i != zPoly.getCoefficients().size(); ++i) {
        //    auto coeff = zPoly[i];
        //    std::cout << msg << "  zPoly [" << i << "]: " << coeff.toHexString(ceil(log2sf / 4.0)) << std::endl;
        //}
        auto I         = ZPolynomial::extractI(zPoly);
        auto maxIValue = 1;  // for log2 calculation
        for (size_t i = 0; i != I.getCoefficients().size(); ++i) {
            auto Idouble = std::abs(I[i].round().convertToDouble());
            if (Idouble > maxIValue) {
                maxIValue = static_cast<int>(Idouble);
            }
        }
        std::cout << msg << "  zPoly Decoded: 0x" << std::hex << decoded << std::dec
                  << "  zPoly error log2Norm: " << std::setprecision(2) << std::fixed
                  << ZPolynomial::extractError(zPoly).getLog2Norm() << " logMaxI: " << std::log2(maxIValue)
                  << std::endl;

        //std::cout << msg << "  zPoly I: ";
        //std::cout << std::endl;
        //std::cout << msg << "  zPoly error log2Norm: " << ZPolynomial::extractError(zPoly).getLog2Norm() << std::endl;
    };

    enum class DecodeMode { RDecode, CSlotsDecode, CSlotsTwiceDecode, CSlotsTwiceBooleanMode, ZDecode };
    std::map<std::string, DecodeMode> decodeMap = {
        // RDecode
        {"ModRaise", DecodeMode::RDecode},
        {"C2R", DecodeMode::RDecode},
        {"Z2R", DecodeMode::RDecode},
        {"Raised", DecodeMode::RDecode},
        // CSlotsDecode
        {"CSlotsDecode", DecodeMode::CSlotsDecode},
        {"Core0", DecodeMode::CSlotsDecode},
        {"Core1", DecodeMode::CSlotsDecode},
        {"LUT0", DecodeMode::CSlotsDecode},
        {"LUT1", DecodeMode::CSlotsDecode},
        {"MSB0", DecodeMode::CSlotsDecode},
        {"MSB1", DecodeMode::CSlotsDecode},
        {"MSB2", DecodeMode::CSlotsDecode},
        {"MSB3", DecodeMode::CSlotsDecode},
        {"A2B0", DecodeMode::CSlotsDecode},
        {"A2B1", DecodeMode::CSlotsDecode},
        {"Comb0", DecodeMode::CSlotsDecode},
        {"Comb1", DecodeMode::CSlotsDecode},
        {"Scaled", DecodeMode::CSlotsDecode},
        {"BoolOR", DecodeMode::CSlotsDecode},
        // CSlotsTwiceDecode
        {"LUT", DecodeMode::CSlotsTwiceDecode},
        {"Normalize", DecodeMode::CSlotsTwiceDecode},
        {"Core", DecodeMode::CSlotsTwiceDecode},
        {"Z2C", DecodeMode::CSlotsTwiceDecode},
        {"R2C", DecodeMode::CSlotsTwiceDecode},
        {"Sine", DecodeMode::CSlotsTwiceDecode},
        {"Boolean", DecodeMode::CSlotsTwiceBooleanMode},
        {"BooleanAgain", DecodeMode::CSlotsTwiceBooleanMode},
        // ZDecode
        {"Input", DecodeMode::ZDecode},
        {"Sub", DecodeMode::ZDecode},
        {"CMult", DecodeMode::ZDecode},
        {"A2A", DecodeMode::ZDecode},
        {"B2A", DecodeMode::ZDecode},
    };

    auto decodeModeIt = decodeMap.find(msg);
    if (decodeModeIt == decodeMap.end()) {
        return 0;
    }
    auto decodeMode = decodeModeIt->second;

    // Check the r encoding?
    if (decodeMode == DecodeMode::RDecode) {
        // Threshold print to avoid too much output
        for (size_t i = 0; i != std::min(values.getCoefficients().size(), 8ul); ++i) {
            auto valueI      = values[i];
            auto p           = BigFixedPoint::positive(1 << 16);  // some random precision
            auto integerPart = (valueI * p).round() / p;
            auto fracPart    = valueI - integerPart;
            auto log2Error   = std::log2(std::abs(fracPart.convertToDouble()));
            if (msg == "Raised") {
                auto N  = cc_global->GetRingDimension();
                auto rN = zN_global * zSlots_global;
                log2Error -= std::log2(N / rN);
            }
            std::cout << msg << "  values [" << i << "]: " << values[i].toHexString(ceil(log2sf / 4.0))
                      << " error: " << log2Error << std::endl;
        }
        if (values.getCoefficients().size() > 4) {
            std::cout << msg << "  ... (total " << values.getCoefficients().size() << " coefficients)" << std::endl;
        }
    }
    if (decodeMode == DecodeMode::CSlotsDecode) {
        auto cSlots          = values.toCSlots();
        auto maxSlotsToPrint = std::min(zSlots_global, size_t(16));
        for (size_t i = 0; i != maxSlotsToPrint; ++i) {
            auto value       = cSlots[i].getReal();
            auto p           = BigFixedPoint::positive(1 << zN_global);
            auto lutPart     = (value * p).round() / p;
            auto fracPart    = value - lutPart;
            double log2Error = std::log2(std::abs(fracPart.convertToDouble()));
            std::cout << msg << "  cSlots Slot " << i << " Reconstructed: " << lutPart.toHexString() << " "
                      << fracPart.toHexString(ceil(log2sf / 4.0)) << " error: " << std::setprecision(2) << log2Error
                      << std::endl;
            //std::cout << msg << "  cSlots Slot " << i << " " << cSlots[i].getReal().toHexString(ceil(log2sf / 4.0))
            //          << std::endl;
        }
        if (zSlots_global > maxSlotsToPrint) {
            std::cout << msg << "  ... (total " << zSlots_global << " slots)" << std::endl;
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
                std::cout << msg << "  Boolean Mode Slot " << i << " Reconstructed: 0x" << std::hex << integerValue
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
                auto value = cSlotsTwice[index].getReal();
                if (msg == "R2C") {
                    // K_SPARSE_ENCAPSULATED scaling
                    value = value * BigFixedPoint::positive(16);
                }
                auto p           = BigFixedPoint::positive(1 << zN_global);
                auto lutPart     = (value * p).round() / p;
                auto fracPart    = value - lutPart;
                double log2Error = std::log2(std::abs(fracPart.convertToDouble()));
                std::cout << msg << "  cSlotsTwice Slot Example " << exampleSlotIndex << " " << index << " : "
                          << lutPart.toHexString() << " " << fracPart.toHexString(ceil(log2sf / 4.0))
                          << " error: " << log2Error << std::endl;
                //std::cout << msg << "  cSlotsTwice Slot Example " << exampleSlotIndex << " " << index << " : "
                //          << cSlotsTwice[index].getReal().toString() << std::endl;
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

void SimpleBootstrapExample(int zN);

int main(int argc, char* argv[]) {
    SimpleBootstrapExample(std::stoi(argv[1]));
}

void SimpleBootstrapExample(int zN) {
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
    AUXMODSIZE_FLEXIBLEMANUAL = 50;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetNumLargeDigits(3);

    std::vector<uint32_t> levelBudget = {2, 2};

    uint32_t mulDepth = 20;
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
    uint32_t zSlots = cc->GetRingDimension() / zN;  // Full packing
    //uint32_t zSlots = cc->GetRingDimension() / zN / 2;  // Maximal sparse packing
    //uint32_t zSlots = 32;
    zN_global     = zN;
    zSlots_global = zSlots;
    std::cout << "Bootstrapping parameters: zN = " << zN << ", zSlots = " << zSlots << std::endl;

    fheZ->EvalBootstrapSetup(*cc, zN, zSlots, levelBudget, {0, 0}, 4, -16, 1, 8, 1);
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
    auto ptxt3        = ZEncodingImpl::encodeBooleanFull(vec, zN, zSlots, elemParamBool, sfBool);
    Plaintext ptxt3_1 = ptxt3[0];
    Plaintext ptxt3_2 = ptxt3[1];
    auto encoded3_1   = Encrypt(ptxt3_1, keyPair.publicKey);
    auto encoded3_2   = Encrypt(ptxt3_2, keyPair.publicKey);
    CiphertextGroup encoded3({encoded3_1, encoded3_2});

    auto ptxt4        = ZEncodingImpl::encodeBooleanFull(vec2, zN, zSlots, elemParamBool, sfBool);
    Plaintext ptxt4_1 = ptxt4[0];
    Plaintext ptxt4_2 = ptxt4[1];
    auto encoded4_1   = Encrypt(ptxt4_1, keyPair.publicKey);
    auto encoded4_2   = Encrypt(ptxt4_2, keyPair.publicKey);
    CiphertextGroup encoded4({encoded4_1, encoded4_2});

    //#define DEBUG

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
    if (1) {
        BENCHMARK(fheZ->EvalBooleanToArith(encoded3), 1, "BooleanToArith");
#ifdef DEBUG
        auto ct2 = z->EvalBooleanOR(encoded3, encoded4);
        auto ct3 = fheZ->EvalBooleanToArith(ct2);
        __heir_debug2(ct3, "B2A");
#endif
    }

    // BooleanToBoolean
    if (1) {
        BENCHMARK((fheZ->EvalBooleanToBooleanFull(encoded3)), 1, "BooleanToBoolean");
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

        //auto ct3 = z->EvalMultFullInZ(ct2, ct2);
        //z->ModReduceInPlace(ct3);
        //__heir_debug2(ct3, "CMult");

        //auto ct4 = fheZ->EvalArithToArith(ct3);
        //__heir_debug2(ct4, "CMult");

        //auto ct5 = fheZ->EvalArithToArith(ct4);
        //__heir_debug2(ct5, "CMult");

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
    if (1) {
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
    if (0) {
#ifdef DEBUG
        auto ctGroupBool = fheZ->EvalArithToBooleanFull(encoded2);
        __heir_debug2(ctGroupBool[0], "A2B0");
        __heir_debug2(ctGroupBool[1], "A2B1");
#endif
        BENCHMARK(fheZ->EvalArithToBooleanFull(encoded2), 1, "ArithToBoolean");
    }
}
