#include "openfhe.h"
#include "scheme/ckksrns/z-user.h"
#include "utils.h"
#include "scheme/ckksrns/z-fhe.h"
#include "scheme/ckksrns/z-pke.h"

using namespace lbcrypto;

void NoiseTestExample(int zN, std::string directive);

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <zN> <test|verify|bench>" << std::endl;
    }
    auto zN = 64;
    if (argc >= 2) {
        zN = std::stoi(argv[1]);
    }
    auto directive = "test";
    ;
    if (argc >= 3) {
        directive = argv[2];
    }

    NoiseTestExample(zN, directive);
}

void NoiseTestExample(int zN, std::string directive) {
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecretKeyDist(lbcrypto::SPARSE_ENCAPSULATED);

    if (directive == "bench" || directive == "verify") {
        parameters.SetSecurityLevel(lbcrypto::HEStd_128_classic);
        parameters.SetRingDim(1 << 16);
    }
    else {  // test
        parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
        parameters.SetRingDim(1 << 10);
    }

    ScalingTechnique rescaleTech = FLEXIBLEMANUAL;
    uint32_t dcrtBits            = 43;
    // bit size for aux moduli in P
    // for HEXL acceleration. Extra 6 bit for SPARSE_ENCAPSULATED
    AUXMODSIZE_FLEXIBLEMANUAL = 50;
    std::cout << "Scaling Factor: " << dcrtBits << " bits\n";
    std::cout << "Auxiliary Prime Size: " << AUXMODSIZE_FLEXIBLEMANUAL << " bits\n";

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetNumLargeDigits(3);

    std::vector<uint32_t> levelBudget = {3, 2};
    std::cout << "Level Budget: ";
    for (auto lb : levelBudget) {
        std::cout << lb << " ";
    }
    std::cout << "\n";

    uint32_t mulDepth = 20;
    parameters.SetMultiplicativeDepth(mulDepth);
    std::cout << "Multiplicative Depth: " << mulDepth << "\n";

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    uint32_t ringDim = cc->GetRingDimension();
    std::cout << "CKKS scheme ring dimension: " << ringDim << "\n";

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
    UserZ u        = std::make_shared<UserZImpl>(z);
    AdvancedZ advZ = std::make_shared<AdvancedZImpl>(z);
    FHEZ fheZ      = std::make_shared<FHEZImpl>(z, advZ);
    PKEZ pkeZ      = std::make_shared<PKEZImpl>(keyPair.publicKey, keyPair.secretKey);
    pkeZ_global    = pkeZ;  // Temporary

    uint32_t zSlots = cc->GetRingDimension() / zN;  // Full packing
    std::cout << "Bootstrapping parameters: zN = " << zN << ", zSlots = " << zSlots << std::endl << std::endl;

    fheZ->EvalBootstrapSetup(*cc, zN, zSlots, levelBudget, {0, 0}, 4, -16, 1);
    fheZ->EvalBootstrapKeyGen(keyPair.secretKey, zN, zSlots);

    //------
    // User program
    //------

    auto elemParam = cc->GetCryptoParameters()->GetElementParams();
    auto sfq0      = cryptoParams->GetScalingFactorBFP(0);

    std::vector<BigInteger> vec(zSlots, 0);
    std::vector<BigInteger> vec2(zSlots, 0);
    for (size_t i = 0; i != zSlots; ++i) {
        vec[i]  = i + rand();  // between 0 and RAND_MAX = 2147483647
        vec2[i] = 2 * i + rand();
    }
    auto ptxt1 = ZEncodingImpl::encodeArith(vec, zN, zSlots, elemParam, sfq0);
    auto ptxt2 = ZEncodingImpl::encodeArith(vec2, zN, zSlots, elemParam, sfq0);

    /// TEST ENCODE
    auto ct  = pkeZ->Encrypt(ptxt1);
    auto ct2 = pkeZ->Encrypt(ptxt2);

    auto ctDec = pkeZ->Decrypt(ct);
    ctDec.print("Arith ct");
    auto ct2Dec = pkeZ->Decrypt(ct2);
    std::cout << std::endl;

    // Test ct-ct for freshly A2A ct

    auto ctA2A  = fheZ->EvalArithToArith(ct);
    auto ct2A2A = fheZ->EvalArithToArith(ct2);
    //
    auto ctA2ADec = pkeZ->Decrypt(ctA2A);
    ctA2ADec.print("Arith ct A2A");
    std::cout << std::endl;
    auto ct2A2ADec = pkeZ->Decrypt(ct2A2A);

    auto ctMult    = u->EvalMultFullInZ(ctA2A, ct2A2A);
    auto ctMultDec = pkeZ->Decrypt(ctMult);

    auto valueMult0         = ctMultDec[0];
    auto valueMultExpected0 = (ctA2ADec[0] * ct2A2ADec[0]) % (BigInteger(1) << zN);
    if (valueMult0 != valueMultExpected0) {
        std::cout << "Error in Mult!" << std::endl;
    }

    auto noiseGrowth    = ctMultDec.getLogMaxNoise() - std::max(ctA2ADec.getLogMaxNoise(), ct2A2ADec.getLogMaxNoise());
    auto overflowGrowth = ctMultDec.getLogMaxI() - std::max(ctA2ADec.getLogMaxI(), ct2A2ADec.getLogMaxI());
    std::cout << "Noise growth after Mult on A2A ct: " << noiseGrowth << " bits." << std::endl;
    std::cout << "Overflow growth after Mult on A2A ct: " << overflowGrowth << " bits." << std::endl;

    auto ctMultMult             = u->EvalMultFullInZ(ctMult, ctMult);
    auto ctMultMultDec          = pkeZ->Decrypt(ctMultMult);
    auto valueMultMult0         = ctMultMultDec[0];
    auto valueMultMultExpected0 = (valueMult0 * valueMult0) % (BigInteger(1) << zN);
    if (valueMultMult0 != valueMultMultExpected0) {
        std::cout << "Error in Mult on Mult!" << std::endl;
    }
    auto noiseGrowthMultMult    = ctMultMultDec.getLogMaxNoise() - ctMultDec.getLogMaxNoise();
    auto overflowGrowthMultMult = ctMultMultDec.getLogMaxI() - ctMultDec.getLogMaxI();
    std::cout << "Noise growth after Mult on Mult ct: " << noiseGrowthMultMult << " bits." << std::endl;
    std::cout << "Overflow growth after Mult on Mult ct: " << overflowGrowthMultMult << " bits." << std::endl;

    auto ctMultMultMult             = u->EvalMultFullInZ(ctMultMult, ctA2A);
    auto ctMultMultMultDec          = pkeZ->Decrypt(ctMultMultMult);
    auto valueMultMultMult0         = ctMultMultMultDec[0];
    auto valueMultMultMultExpected0 = (valueMultMult0 * vec[0]) % (BigInteger(1) << zN);
    if (valueMultMultMult0 != valueMultMultMultExpected0) {
        std::cout << "Error in Mult on Mult on Mult!" << std::endl;
    }
    auto noiseGrowthMultMultMult    = ctMultMultMultDec.getLogMaxNoise() - ctMultMultDec.getLogMaxNoise();
    auto overflowGrowthMultMultMult = ctMultMultMultDec.getLogMaxI() - ctMultMultDec.getLogMaxI();
    std::cout << "Noise growth after Mult on Mult on Mult ct: " << noiseGrowthMultMultMult << " bits." << std::endl;
    std::cout << "Overflow growth after Mult on Mult on Mult ct: " << overflowGrowthMultMultMult << " bits."
              << std::endl;

    // Test ct-pt mult, i.e. multshort

    auto ctMultShort = u->EvalMultPtInZ(ctA2A, BigInteger(vec2[0]));

    auto ctMultShortDec          = pkeZ->Decrypt(ctMultShort);
    auto valueMultShort0         = ctMultShortDec[0];
    auto valueMultShortExpected0 = (ctA2ADec[0] * BigInteger(vec2[0])) % (BigInteger(1) << zN);
    if (valueMultShort0 != valueMultShortExpected0) {
        std::cout << "Error in MultShort!" << std::endl;
    }

    auto noiseGrowthShort    = ctMultShortDec.getLogMaxNoise() - ctA2ADec.getLogMaxNoise();
    auto overflowGrowthShort = ctMultShortDec.getLogMaxI() - ctA2ADec.getLogMaxI();

    std::cout << "Noise growth after MultShort on A2A ct: " << noiseGrowthShort << " bits." << std::endl;
    std::cout << "Overflow growth after MultShort on A2A ct: " << overflowGrowthShort << " bits." << std::endl;

    auto ctMultShortMultShort             = u->EvalMultPtInZ(ctMultShort, BigInteger(vec2[0]));
    auto ctMultShortMultShortDec          = pkeZ->Decrypt(ctMultShortMultShort);
    auto valueMultShortMultShort0         = ctMultShortMultShortDec[0];
    auto valueMultShortMultShortExpected0 = (valueMultShort0 * BigInteger(vec2[0])) % (BigInteger(1) << zN);
    if (valueMultShortMultShort0 != valueMultShortMultShortExpected0) {
        std::cout << "Error in MultShort on MultShort!" << std::endl;
    }
    auto noiseGrowthShortMultShort    = ctMultShortMultShortDec.getLogMaxNoise() - ctMultShortDec.getLogMaxNoise();
    auto overflowGrowthShortMultShort = ctMultShortMultShortDec.getLogMaxI() - ctMultShortDec.getLogMaxI();
    std::cout << "Noise growth after MultShort on MultShort ct: " << noiseGrowthShortMultShort << " bits." << std::endl;
    std::cout << "Overflow growth after MultShort on MultShort ct: " << overflowGrowthShortMultShort << " bits."
              << std::endl;

    auto ctMultShortMultShortMultShort     = u->EvalMultPtInZ(ctMultShortMultShort, BigInteger(vec2[0]));
    auto ctMultShortMultShortMultShortDec  = pkeZ->Decrypt(ctMultShortMultShortMultShort);
    auto valueMultShortMultShortMultShort0 = ctMultShortMultShortMultShortDec[0];
    auto valueMultShortMultShortMultShortExpected0 =
        (valueMultShortMultShort0 * BigInteger(vec2[0])) % (BigInteger(1) << zN);
    if (valueMultShortMultShortMultShort0 != valueMultShortMultShortMultShortExpected0) {
        std::cout << "Error in MultShort on MultShort on MultShort!" << std::endl;
    }
    auto noiseGrowthShortMultShortMultShort =
        ctMultShortMultShortMultShortDec.getLogMaxNoise() - ctMultShortMultShortDec.getLogMaxNoise();
    auto overflowGrowthShortMultShortMultShort =
        ctMultShortMultShortMultShortDec.getLogMaxI() - ctMultShortMultShortDec.getLogMaxI();
    std::cout << "Noise growth after MultShort on MultShort on MultShort ct: " << noiseGrowthShortMultShortMultShort
              << " bits." << std::endl;
    std::cout << "Overflow growth after MultShort on MultShort on MultShort ct: "
              << overflowGrowthShortMultShortMultShort << " bits." << std::endl;
}