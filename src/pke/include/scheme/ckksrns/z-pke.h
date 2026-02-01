#ifndef SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_DEBUG_H_
#define SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_DEBUG_H_

#include "encoding/z-encoding.h"
#include "scheme/ckksrns/z-fhe.h"

namespace lbcrypto {

struct ZDecryptResult {
    ZDecryptResult() : logMaxNoise(0), logMaxI(0) {}
    ZDecryptResult(std::vector<BigInteger> vals, double logNoise, double logI)
        : values(vals), logMaxNoise(logNoise), logMaxI(logI) {}

    size_t size() const {
        return values.size();
    }

    BigInteger operator[](size_t index) const {
        return values[index];
    }

    double getLogMaxNoise() const {
        return logMaxNoise;
    }

    double getLogMaxI() const {
        return logMaxI;
    }

    void print(std::string msg, size_t maxSlotsToPrint = 4) const;

    bool valuesEqual(const ZDecryptResult& other) const;

    double noiseComparison(const ZDecryptResult& other) const {
        return logMaxNoise - other.logMaxNoise;
    }

    void printNoiseComparison(const ZDecryptResult& other, std::string msg) const {
        double diff = noiseComparison(other);
        std::cout << msg << " Current Noise " << logMaxNoise << ", Difference: " << diff << std::endl;
    }

private:
    std::vector<BigInteger> values;
    double logMaxNoise;
    double logMaxI;
};

class PKEZImpl {
public:
    PKEZImpl(PublicKey<DCRTPoly> pk, PrivateKey<DCRTPoly> sk) : pk(pk), sk(sk) {};

    void debug(Ciphertext<DCRTPoly> ct, std::string msg);

    ZDecryptResult Decrypt(CiphertextGroup ct);

    Ciphertext<DCRTPoly> Encrypt(ZEncoding ptxt);
    CiphertextGroup Encrypt(std::vector<ZEncoding> ptxts);

private:
    PublicKey<DCRTPoly> pk;
    PrivateKey<DCRTPoly> sk;
};

using PKEZ = std::shared_ptr<PKEZImpl>;

extern PKEZ pkeZ_global;

}  // namespace lbcrypto

#endif  // SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_DEBUG_H_