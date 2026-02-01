#include "openfhe.h"
#include "utils.h"
#include "scheme/ckksrns/z-fhe.h"
#include "scheme/ckksrns/z-pke.h"

using namespace lbcrypto;

void SimpleBootstrapExample(int zN);

int main(int argc, char* argv[]) {
    SimpleBootstrapExample(std::stoi(argv[1]));
}

void SimpleBootstrapExample(int zN) {
    CCParams<CryptoContextCKKSRNS> parameters;

    SecretKeyDist secretKeyDist = lbcrypto::SPARSE_ENCAPSULATED;
    parameters.SetSecretKeyDist(secretKeyDist);

    parameters.SetSecurityLevel(HEStd_NotSet);
    // 1 << 12 is buggy? WHY?
    parameters.SetRingDim(1 << 10);
    //parameters.SetNumLargeDigits(6);

    ScalingTechnique rescaleTech = FLEXIBLEMANUAL;
    uint32_t dcrtBits            = 43;
    // bit size for aux moduli in P
    // for HEXL acceleration. Extra 6 bit for SPARSE_ENCAPSULATED
    AUXMODSIZE_FLEXIBLEMANUAL = 50;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetNumLargeDigits(3);

    std::vector<uint32_t> levelBudget = {3, 2};

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
    PKEZ pkeZ      = std::make_shared<PKEZImpl>(keyPair.publicKey, keyPair.secretKey);
    pkeZ_global    = pkeZ;  // Temporary

    uint32_t zSlots = cc->GetRingDimension() / zN;  // Full packing
    std::cout << "Bootstrapping parameters: zN = " << zN << ", zSlots = " << zSlots << std::endl;

    fheZ->EvalBootstrapSetup(*cc, zN, zSlots, levelBudget, {0, 0}, 4, -16, 1);
    fheZ->EvalBootstrapKeyGen(keyPair.secretKey, zN, zSlots);

    auto elemParam = cc->GetCryptoParameters()->GetElementParams();
    auto sfq0      = cryptoParams->GetScalingFactorBFP(0);

    std::vector<uint64_t> vec(zSlots, 0);
    for (size_t i = 0; i != zSlots; ++i) {
        //vec[i] = i + 3;
        vec[i] = i + 0xdeadbeaf;
    }
    auto ptxt1 = ZEncodingImpl::encodeArith(vec, zN, zSlots, elemParam, sfq0);

    std::vector<uint64_t> vec2(zSlots, 0);
    for (size_t i = 0; i != zSlots; ++i) {
        //vec2[i] = -i - 3;
        vec2[i] = i + 0xf0f0f0ff;
        //vec2[i] = i + 16;
    }
    auto ptxt2 = ZEncodingImpl::encodeArith(vec2, zN, zSlots, elemParam, sfq0);

    /// TEST ENCODE
    auto ct  = pkeZ->Encrypt(ptxt1);
    auto ct2 = pkeZ->Encrypt(ptxt2);

    auto ctDec = pkeZ->Decrypt(ct);
    ctDec.print("ct");
    auto ct2Dec = pkeZ->Decrypt(ct2);
    ctDec.print("ct2");

    /// TEST Boolean Encode
    auto ptxt3          = ZEncodingImpl::encodeBooleanFull(vec, zN, zSlots, elemParam, sfq0);
    CiphertextGroup ct3 = pkeZ->Encrypt(ptxt3);

    auto ptxt4          = ZEncodingImpl::encodeBooleanFull(vec2, zN, zSlots, elemParam, sfq0);
    CiphertextGroup ct4 = pkeZ->Encrypt(ptxt4);

    auto ct3Dec = pkeZ->Decrypt(ct3);
    ct3Dec.print("ct3");
    auto ct4Dec = pkeZ->Decrypt(ct4);
    ct3Dec.print("ct4");

    //#define DEBUG

#ifdef DEBUG
    __heir_debug2(encoded, "Input");
#endif

    // Mult
    if (1) {
        //BENCHMARK(z->EvalMultFullInZ(encoded2, encoded), 3, "MultFull");
        auto ctRes    = z->EvalMultFullInZ(ct, ct);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        ctResDec.print("Mult Result");
#ifdef DEBUG
        auto ct2 = z->EvalMultFullInZ(encoded, encoded);
        __heir_debug2(ct2, "CMult");
#endif
    }

    // Bool
    if (1) {
        auto ctRes    = z->EvalBooleanOR(ct3, ct4);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        ctResDec.print("BoolOR Result");
        //BENCHMARK(z->EvalBooleanOR(encoded3, encoded4), 3, "BooleanOR");
#ifdef DEBUG
        auto ct2 = z->EvalBooleanOR(encoded3, encoded4);
        __heir_debug2(ct2[0], "BoolOR");
#endif
    }

    // BoolToArith
    if (0) {
        BENCHMARK(fheZ->EvalBooleanToArith(ct3), 1, "BooleanToArith");
#ifdef DEBUG
        auto ct2 = z->EvalBooleanOR(encoded3, encoded4);
        auto ct3 = fheZ->EvalBooleanToArith(ct2);
        __heir_debug2(ct3, "B2A");
#endif
    }

    // BooleanToBoolean
    if (0) {
        BENCHMARK((fheZ->EvalBooleanToBooleanFull(ct3)), 1, "BooleanToBoolean");
    }

    // ArithToArithHigh
    if (0) {
        BENCHMARK(fheZ->EvalArithToArithHigh(ct2), 1, "ArithToArithHigh");
    }

    // ArithToArithNoise
    if (1) {
        //BENCHMARK(fheZ->EvalArithToArithNoise(ct2), 1, "ArithToArithNoise");
        auto ctRes  = fheZ->EvalArithToArithHigh(ct2);
        auto a2aDec = pkeZ->Decrypt(ctRes);
        a2aDec.print("A2Ae");
        a2aDec.printNoiseComparison(ct2Dec, "A2Ae vs ct2");

        auto ctRes2 = z->EvalMultFullInZ(ctRes, ctRes);
        z->ModReduceInPlace(ctRes2);
        ctRes2 = z->EvalMultFullInZ(ctRes2, ctRes);
        z->ModReduceInPlace(ctRes2);
        auto a2aMulDec = pkeZ->Decrypt(ctRes2);
        a2aMulDec.print("A2Ae Mult");
        a2aMulDec.printNoiseComparison(ct2Dec, "A2Ae Mult vs ct2");

        auto ctRes3     = fheZ->EvalArithToArith(ctRes2);
        auto a2aMul2Dec = pkeZ->Decrypt(ctRes3);
        a2aMul2Dec.print("A2Ae Mult 2");
        a2aMul2Dec.printNoiseComparison(ct2Dec, "A2Ae Mult 2 vs ct2");

        //auto ctRes4     = fheZ->EvalArithToArith(ctRes3);
        //auto a2aMul3Dec = pkeZ->Decrypt(ctRes4);
        //a2aMul3Dec.print("A2Ae Mult 3");
        //a2aMul3Dec.printNoiseComparison(ct2Dec, "A2Ae Mult 3 vs ct2");
    }

    // ArithToArith
    if (1) {
        //BENCHMARK(fheZ->EvalArithToArith(encoded2), 1, "ArithToArith");
    }

    // ArithToBooleanBatched
    if (0) {
        std::vector<Ciphertext<DCRTPoly>> batchCts = {ct, ct2};
        auto batchSize                             = zN / 4;  // zN / w
        while (batchCts.size() < batchSize) {
            batchCts.push_back(ct);
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
        BENCHMARK(fheZ->EvalArithToBooleanFull(ct2), 1, "ArithToBoolean");
    }
}