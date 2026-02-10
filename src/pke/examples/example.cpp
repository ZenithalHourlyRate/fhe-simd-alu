#include "openfhe.h"
#include "scheme/ckksrns/z-user.h"
#include "scheme/ckksrns/z-user-advanced.h"
#include "scheme/ckksrns/z-fhe.h"
#include "scheme/ckksrns/z-pke.h"

using namespace lbcrypto;

void SimpleExample();

int main(int argc, char* argv[]) {
    SimpleExample();
}

void SimpleExample() {
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecretKeyDist(lbcrypto::SPARSE_ENCAPSULATED);

    // parameters.SetSecurityLevel(lbcrypto::HEStd_128_classic);
    // parameters.SetRingDim(1 << 16);
    // Toy parameter for easy testing; set to above parameter for real security
    parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
    parameters.SetRingDim(1 << 10);

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
    auto moduliQ            = cc->GetCryptoParameters()->GetElementParams()->GetModulus();
    auto moduliP            = cryptoParams->GetParamsP()->GetModulus();
    auto logQ               = moduliQ.GetMSB();
    auto logP               = moduliP.GetMSB();
    std::cout << "log2(Q) = " << logQ << " log2(P) = " << logP << " log2(QP) = " << logQ + logP << std::endl;

    LeveledZ z         = std::make_shared<LeveledZImpl>();
    UserZ u            = std::make_shared<UserZImpl>(z);
    AdvancedZ advZ     = std::make_shared<AdvancedZImpl>(z);
    FHEZ fheZ          = std::make_shared<FHEZImpl>(z, advZ);
    PKEZ pkeZ          = std::make_shared<PKEZImpl>(keyPair.publicKey, keyPair.secretKey);
    UserZAdvanced uAdv = std::make_shared<UserZAdvancedImpl>(z, fheZ, u);

    uint32_t zN     = 64;
    uint32_t zSlots = cc->GetRingDimension() / zN;  // Full packing
    std::cout << "Bootstrapping parameters: zN = " << zN << ", zSlots = " << zSlots << std::endl << std::endl;

    // dim1 for BSGS related to CKKS bootstrapping, set to 0 for default
    // w = 4, the A2B LUT parameter
    // a2b-cutoff = -24
    // lutOrder = 1, the A2B LUT parameter
    fheZ->EvalBootstrapSetup(*cc, zN, zSlots, levelBudget, {0, 0}, 4, -24, 1);
    fheZ->EvalBootstrapKeyGen(keyPair.secretKey, zN, zSlots);

    // Enable sign extraction
    // zSlotRotates = {1, -1} for one left and one right rotate in the slot dimension, which is useful for many applications
    // leftShifts = {8} for shift left by 8 bits, rightShifts = {8} for shift right by 8 bits, leftRotates = {8} for rotate left by 8 bits, rightRotates = {8} for rotate right by 8 bits
    u->Setup(keyPair.secretKey, zN, zSlots, true, {1, -1}, {8}, {8}, {0}, {8});

    auto elemParam = cc->GetCryptoParameters()->GetElementParams();
    auto sfq0      = cryptoParams->GetScalingFactorBFP(0);

    std::vector<BigInteger> vec(zSlots, 0);
    for (size_t i = 0; i != zSlots; ++i) {
        vec[i] = i + 0xdeadbeaf;
    }
    auto ptxt1 = ZEncodingImpl::encodeArith(vec, zN, zSlots, elemParam, sfq0);

    std::vector<BigInteger> vec2(zSlots, 0);
    for (size_t i = 0; i != zSlots; ++i) {
        vec2[i] = i + 0xf0f0f0ff;
    }
    auto ptxt2 = ZEncodingImpl::encodeArith(vec2, zN, zSlots, elemParam, sfq0);
    std::vector<BigInteger> vec3(zSlots, 0);
    for (size_t i = 0; i != zSlots; ++i) {
        vec3[i] = i + 3;
    }

    /// TEST ENCODE
    auto ct  = pkeZ->Encrypt(ptxt1);
    auto ct2 = pkeZ->Encrypt(ptxt2);

    auto ctDec = pkeZ->Decrypt(ct);
    ctDec.print("Arith ct");
    auto ct2Dec = pkeZ->Decrypt(ct2);
    ct2Dec.print("Arith ct2");

    /// TEST Boolean Encode
    auto ptxt3          = ZEncodingImpl::encodeBooleanFull(vec, zN, zSlots, elemParam, sfq0);
    CiphertextGroup ct3 = pkeZ->Encrypt(ptxt3);
    auto ct3Dec         = pkeZ->Decrypt(ct3);
    ct3Dec.print("Bool ct3");

    auto ptxt4          = ZEncodingImpl::encodeBooleanFull(vec2, zN, zSlots, elemParam, sfq0);
    CiphertextGroup ct4 = pkeZ->Encrypt(ptxt4);
    auto ct4Dec         = pkeZ->Decrypt(ct4);
    ct4Dec.print("Bool ct4");

    std::cout << std::endl;

    // Add
    if (1) {
        auto ctRes    = u->EvalAddPtVecInZ(ct, vec2);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        if (ctResDec[0] != (ctDec[0] + ct2Dec[0]) % (BigInteger(1) << zN)) {
            std::cout << "Error in Add!" << std::endl;
        }
        ctResDec.print("Add");
        ctResDec.printNoiseComparison(ctDec, "Add");
        std::cout << std::endl;
    }

    // Sub
    if (1) {
        auto ctRes    = u->EvalSubInZ(ct2, ct);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        if (ctResDec[0] != (ct2Dec[0] - ctDec[0] + (BigInteger(1) << zN)) % (BigInteger(1) << zN)) {
            std::cout << "Error in Sub!" << std::endl;
        }
        ctResDec.print("Sub");
        ctResDec.printNoiseComparison(ct2Dec, "Sub");
        std::cout << std::endl;
    }

    // MultFull
    if (1) {
        auto ctRes    = u->EvalMultFullInZ(ct, ct);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        if (ctResDec[0] != (ctDec[0] * ctDec[0]) % (BigInteger(1) << zN)) {
            std::cout << "Error in MultFull!" << std::endl;
        }
        ctResDec.print("MultFull");
        ctResDec.printNoiseComparison(ctDec, "MultFull");
        std::cout << std::endl;
    }

    // MultShort ct-pt
    if (1) {
        auto ctRes    = u->EvalMultPtVecInZ(ct, vec3);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        if (ctResDec[0] != (ctDec[0] * BigInteger(vec3[0])) % (BigInteger(1) << zN)) {
            std::cout << "Error in MultShort!" << std::endl;
        }
        ctResDec.print("MultShort");
        ctResDec.printNoiseComparison(ctDec, "MultShort");
        std::cout << std::endl;
    }

    // BoolAND
    if (1) {
        auto ctRes    = u->EvalBooleanAND(ct3, ct4);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        if (ctResDec[0].ConvertToInt() != (ct3Dec[0].ConvertToInt() & ct4Dec[0].ConvertToInt())) {
            std::cout << "Error in BooleanAND!" << std::endl;
        }
        ctResDec.print("BooleanAND");
        ctResDec.printNoiseComparison(ctDec, "BooleanAND");
        std::cout << std::endl;
    }

    // BoolAND with ptxt
    if (1) {
        auto ctRes    = u->EvalBooleanANDPt(ct3, 0xff);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        if (ctResDec[0].ConvertToInt() != (ct3Dec[0].ConvertToInt() & 0xff)) {
            std::cout << "Error in BooleanAND pt!" << std::endl;
        }
        ctResDec.print("BooleanAND pt");
        ctResDec.printNoiseComparison(ctDec, "BooleanAND pt");
        std::cout << std::endl;
    }

    // zSlot Rotate Right
    if (1) {
        auto ctRes    = uAdv->EvalRotateInZ(ct, 1);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        if (ctResDec[0] != ct3Dec[1]) {
            std::cout << "Error in zSlot Rotate Right!" << std::endl;
        }
        ctResDec.print("zSlot Rotate Right");
        ctResDec.printNoiseComparison(ct3Dec, "zSlot Rotate Right");
        std::cout << std::endl;
    }

    // zSlot Rotate Left
    if (1) {
        auto ctRes    = uAdv->EvalRotateInZ(ct, -1);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        if (ctResDec[0] != ct3Dec[zSlots - 1]) {
            std::cout << "Error in zSlot Rotate Left!" << std::endl;
        }
        ctResDec.print("zSlot Rotate Left");
        ctResDec.printNoiseComparison(ct3Dec, "zSlot Rotate Left");
        std::cout << std::endl;
    }

    // ShiftLeft
    if (1) {
        uint64_t shift      = 8;
        auto ctRes          = u->EvalBooleanShiftLeft(ct3, shift);
        auto ctResDec       = pkeZ->Decrypt(ctRes);
        BigInteger expected = (vec[0] << shift) % (BigInteger(1) << zN);
        if (ctResDec[0] != expected) {
            std::cout << "Error in ShiftLeft!" << std::endl;
        }
        ctResDec.print("ShiftLeft");
        ctResDec.printNoiseComparison(ct3Dec, "ShiftLeft");
        std::cout << std::endl;
    }

    // ShiftRight
    if (1) {
        uint64_t shift      = 8;
        auto ctRes          = u->EvalBooleanShiftRight(ct3, shift);
        auto ctResDec       = pkeZ->Decrypt(ctRes);
        BigInteger expected = (vec[0] >> shift) % (BigInteger(1) << zN);
        if (ctResDec[0] != expected) {
            std::cout << "Error in ShiftRight!" << std::endl;
        }
        ctResDec.print("ShiftRight");
        ctResDec.printNoiseComparison(ct3Dec, "ShiftRight");
        std::cout << std::endl;
    }

    // RotateRight
    if (1) {
        uint64_t shift = 8;
        auto ctRes     = u->EvalBooleanRotateRight(ct3, shift);
        auto ctResDec  = pkeZ->Decrypt(ctRes);
        BigInteger expected;
        expected = BigInteger((vec[0] >> shift) | (vec[0] << (zN - shift))) % (BigInteger(1) << zN);
        if (ctResDec[0] != expected) {
            std::cout << "Error in RotateRight!" << std::endl;
        }
        ctResDec.print("RotateRight");
        ctResDec.printNoiseComparison(ct3Dec, "RotateRight");
        std::cout << std::endl;
    }

    // Compare
    if (1) {
        std::vector<BigInteger> compareVec(zSlots, 0);
        for (size_t i = 0; i != zSlots; ++i) {
            if (i % 2 == 0)
                compareVec[i] = vec[i] + BigInteger(10);
            else
                compareVec[i] = vec[i] - BigInteger(10);
        }
        auto ptxtCompare = ZEncodingImpl::encodeArith(compareVec, zN, zSlots, elemParam, sfq0);
        auto ctCompare   = pkeZ->Encrypt(ptxtCompare);
        auto ctRes       = uAdv->EvalLessThan(ct, ctCompare);
        auto ctResDec    = pkeZ->Decrypt(ctRes);
        for (size_t i = 0; i != zSlots; ++i) {
            BigInteger expected = (vec[i] < compareVec[i]) ? BigInteger(1) : BigInteger(0);
            if (ctResDec[i] != expected) {
                std::cout << "Error in Compare at slot " << i << "!" << std::endl;
            }
        }
        ctResDec.print("Compare");
        ctResDec.printNoiseComparison(ctDec, "Compare");
        std::cout << std::endl;
    }

    // BoolToArith
    if (1) {
        auto ctRes    = fheZ->EvalBooleanToArith(ct3);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        if (!ctResDec.valuesEqual(ct3Dec)) {
            std::cout << "Error in B2A!" << std::endl;
        }
        ctResDec.printNoiseComparison(ctDec, "B2A");
        std::cout << std::endl;
    }

    // BooleanToBoolean
    if (1) {
        auto ctRes    = fheZ->EvalBooleanToBooleanFull(ct3);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        if (!ctResDec.valuesEqual(ct3Dec)) {
            std::cout << "Error in B2B!" << std::endl;
        }
        ctResDec.printNoiseComparison(ct3Dec, "B2B");
        std::cout << std::endl;
    }

    // ArithToArith
    if (1) {
        auto ctRes    = fheZ->EvalArithToArith(ct2);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        if (!ctResDec.valuesEqual(ct2Dec)) {
            std::cout << "Error in A2A!" << std::endl;
        }
        ctResDec.printNoiseComparison(ct2Dec, "A2A");
        std::cout << std::endl;
    }

    // ArithToBooleanBatched
    if (1) {
        std::vector<Ciphertext<DCRTPoly>> batchCts = {ct, ct2};
        auto batchSize                             = zN / 4;  // zN / w
        while (batchCts.size() < batchSize) {
            batchCts.push_back(ct);
        }

        auto ctGroupBool = fheZ->EvalArithToBooleanBatched(batchCts);
        auto ctResDec0   = pkeZ->Decrypt(CiphertextGroup({ctGroupBool[0], ctGroupBool[1]}));
        auto ctResDec1   = pkeZ->Decrypt(CiphertextGroup({ctGroupBool[2], ctGroupBool[3]}));
        if (!ctResDec0.valuesEqual(ctDec)) {
            std::cout << "Error in B-A2B!" << std::endl;
        }
        if (!ctResDec1.valuesEqual(ct2Dec)) {
            std::cout << "Error in B-A2B!" << std::endl;
        }
        ctResDec0.printNoiseComparison(ctDec, "B-A2B-0");
        ctResDec1.printNoiseComparison(ct2Dec, "B-A2B-1");
        ctResDec0.print("B-A2B-0");
        std::cout << std::endl;
    }

    //ArithToBooleanFull
    if (1) {
        auto ctRes    = fheZ->EvalArithToBooleanFull(ct2);
        auto ctResDec = pkeZ->Decrypt(ctRes);
        if (!ctResDec.valuesEqual(ct2Dec)) {
            std::cout << "Error in A2B!" << std::endl;
        }
        ctResDec.printNoiseComparison(ct2Dec, "A2B");
        ctResDec.print("A2B");
        std::cout << std::endl;
    }
    return;
}