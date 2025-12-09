#ifndef SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_ADVANCEDSHE_H_
#define SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_ADVANCEDSHE_H_

#include "math/hal/bigfixedpoint.h"
#include "cryptocontext.h"
#include "z-leveledshe.h"

namespace lbcrypto {

//=============================================================================
// Pre-defined Chebyshev Series Coefficients
//=============================================================================

// Coefficients for the function std::exp(1i * Pi/2.0 * x) in [-16, 16] of degree 58
// Need two double-angle iterations to get std::exp(1i * 2Pi * x)
// This gives -50 bit precision
static const inline std::vector<BigComplex> coeff_exp_16_big_complex_58 = {
    BigComplex(BigFixedPoint(BigInteger("76201359508406144086108162329195173168"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, true),
               BigFixedPoint(BigInteger("75462327482981010929732591801834267864"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("82206460726474050179385148728921771719"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("62378762904272904571622519573888623801"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("97098293416732629097600305219456855499"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("31471416178152430358852601914596611953"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("109620372043870639737320939173783651183"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, true),
               BigFixedPoint(BigInteger("20868456044911447829692299230268924039"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("97995759274018037306194755132115244659"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, true),
               BigFixedPoint(BigInteger("83254494006925109073099868507868876972"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("38369120122842855712446721811003270438"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("113787669655110863451684674786681491532"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("61235165364009198111122209480690318755"), 128, true),
               BigFixedPoint(BigInteger("0"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("55312394102720757898644017695158058423"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("118456231446026759188540160451912549635"), 128, true),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("76657869319988643285788486120579931992"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("26952640227261194763905894118976297191"), 128, true),
               BigFixedPoint(BigInteger("0"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, true),
               BigFixedPoint(BigInteger("110975036692357318664215943016268684310"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("123176277792580374930872323081376113165"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, true),
               BigFixedPoint(BigInteger("65461984648789268549752845531057520576"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("24199592598704844843639017755231204452"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("103976832477728479972166724489186183378"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("149558889387913469944174054915466899868"), 128, true),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("157856569345896333137463064213636837696"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("139363123687758077275148430791977316063"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, true),
               BigFixedPoint(BigInteger("108307390889710558441771144735781715323"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("76107584230786263566850937818652147411"), 128, true),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("49160286192750883756182833588143160122"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("29517800212959227901643325698970465645"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, true),
               BigFixedPoint(BigInteger("16610367190529583400159946573291359467"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("8814719438731980311580025064597388180"), 128, true),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("4433225369106595419301399779194470128"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("2121611396940476463520449109282130128"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("969413689184254142896588577455705847"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("424123781880824174246087041777004704"), 128, true),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("178110049091499884110719962349032111"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("71950375846091926345332846912452789"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, true),
               BigFixedPoint(BigInteger("28012594426554458521141725683489357"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("10528967282478530758575404449167301"), 128, true),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("3826412157972211940085931935062259"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("1346384697861443672052779089819636"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, true),
               BigFixedPoint(BigInteger("459263441386523328060768631772001"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("152043263138788980611771564757969"), 128, true),
               BigFixedPoint(BigInteger("0"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("48903733353037260271522149472554"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("15297061185064203581813928727122"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, true),
               BigFixedPoint(BigInteger("4657530499369695341684273052327"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("1381491350394687807178958024933"), 128, true),
               BigFixedPoint(BigInteger("0"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, true),
               BigFixedPoint(BigInteger("399506577472584567012933042613"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("112719644763759278673254371797"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("31050750072564917010226630651"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("8356424074880620222773761537"), 128, true),
               BigFixedPoint(BigInteger("0"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("2198404879529214920726341624"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("565694514958724793776479624"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("142455153436637585591268534"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("35125212779865136607144117"), 128, true),
               BigFixedPoint(BigInteger("0"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("8484235440507313142482802"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("2008720393085612905928457"), 128, false),
               BigFixedPoint(BigInteger("0"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("0"), 128, false),
               BigFixedPoint(BigInteger("465076165811325908974862"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("111346816677428262205798"), 128, true),
               BigFixedPoint(BigInteger("0"), 128, true)),
};
class AdvancedZImpl {
public:
    AdvancedZImpl(LeveledZ z) : z(z) {}

    Ciphertext<DCRTPoly> EvalChebyshevSeriesPS(ConstCiphertext<DCRTPoly>& x, const std::vector<BigComplex>& coeffs);

    std::shared_ptr<seriesPowers<DCRTPoly>> EvalPowers(ConstCiphertext<DCRTPoly>& ciphertext,
                                                       const std::vector<BigComplex>& coefficients);

    Ciphertext<DCRTPoly> EvalPolyWithPrecomp(std::shared_ptr<seriesPowers<DCRTPoly>> ctxtPowers,
                                             const std::vector<BigComplex>& coeffs);

private:
    // WSum related
    Ciphertext<DCRTPoly> EvalPartialLinearWSum(const std::vector<Ciphertext<DCRTPoly>>& ciphertexts,
                                               const std::vector<BigComplex>& constants, uint32_t limit = 0);

    // ChebyshevPS related
    std::shared_ptr<seriesPowers<DCRTPoly>> internalEvalChebyPolysPS(ConstCiphertext<DCRTPoly>& x, uint32_t degree);

    Ciphertext<DCRTPoly> InnerEvalChebyshevPS(ConstCiphertext<DCRTPoly>& x, const std::vector<BigComplex>& coefficients,
                                              uint32_t k, uint32_t m, const std::vector<Ciphertext<DCRTPoly>>& T,
                                              const std::vector<Ciphertext<DCRTPoly>>& T2);

    Ciphertext<DCRTPoly> internalEvalChebyshevSeriesPSWithPrecomp(
        const std::shared_ptr<seriesPowers<DCRTPoly>>& ctxtPolys, const std::vector<BigComplex>& coefficients);

    // PS related
    std::shared_ptr<seriesPowers<DCRTPoly>> internalEvalPowersPS(ConstCiphertext<DCRTPoly>& x, uint32_t degree);

    Ciphertext<DCRTPoly> InnerEvalPolyPS(ConstCiphertext<DCRTPoly>& x, const std::vector<BigComplex>& coefficients,
                                         uint32_t k, uint32_t m, const std::vector<Ciphertext<DCRTPoly>>& powers,
                                         const std::vector<Ciphertext<DCRTPoly>>& powers2);

    Ciphertext<DCRTPoly> internalEvalPolyPSWithPrecomp(const std::shared_ptr<seriesPowers<DCRTPoly>>& ctxtPowers,
                                                       const std::vector<BigComplex>& coefficients);

private:
    LeveledZ z;
};

using AdvancedZ = std::shared_ptr<AdvancedZImpl>;

}  // namespace lbcrypto

#endif  //SRC_PKE_INCLUDE_SCHEME_CKKSRNS_Z_ADVANCEDSHE_H_