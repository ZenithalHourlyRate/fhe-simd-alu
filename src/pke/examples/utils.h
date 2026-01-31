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

#define BENCHMARK(x, times, str)                                                           \
    {                                                                                      \
        auto startW = std::chrono::high_resolution_clock::now();                           \
        (x);                                                                               \
        auto endW                           = std::chrono::high_resolution_clock::now();   \
        std::chrono::duration<double> diffW = endW - startW;                               \
        std::cout << "Finished Warmup for " str " " << diffW.count() << " s" << std::endl; \
        auto start = std::chrono::high_resolution_clock::now();                            \
        for (size_t i = 0; i != (times); ++i)                                              \
            (x);                                                                           \
        auto end                           = std::chrono::high_resolution_clock::now();    \
        std::chrono::duration<double> diff = end - start;                                  \
        std::cout << "Time for " str " : " << diff.count() / (times) << " s" << std::endl; \
    }                                                                                      \
    while (0)

static inline std::string bigIntegerToHexString(BigInteger value) {
    std::stringstream ss;
    NTL::ZZ v = value;

    if (v == 0) {
        return "0x0";
    }

    long nbytes = NumBytes(v);
    std::vector<unsigned char> buf(nbytes);
    BytesFromZZ(buf.data(), v, nbytes);
    std::reverse(buf.begin(), buf.end());

    std::ostringstream oss;
    oss << "0x";

    for (unsigned char b : buf) {
        oss << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(b);
    }
    return oss.str();
};