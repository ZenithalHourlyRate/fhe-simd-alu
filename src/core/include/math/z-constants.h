#ifndef SRC_CORE_INCLUDE_MATH_Z_CONSTANTS_H_
#define SRC_CORE_INCLUDE_MATH_Z_CONSTANTS_H_

#include "math/hal/bigfixedpoint.h"

// For X^32 - X + 2
const std::vector<BigComplex> z_upper_roots_32 = {
    BigComplex(BigFixedPoint(BigInteger("350551088990192644005587638168552186653"), 128, true),
               BigFixedPoint(BigInteger("34900550222895179714397787198203513173"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("336695400736306456376811924693546769860"), 128, true),
               BigFixedPoint(BigInteger("103299331694394466680808556291558998110"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("309551080397728457701091292115470551977"), 128, true),
               BigFixedPoint(BigInteger("167549268673586454371916803328721791984"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("270228898653084472294462899074132324165"), 128, true),
               BigFixedPoint(BigInteger("225075857771326758427915586988199792655"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("220337572168790987178561844883347289138"), 128, true),
               BigFixedPoint(BigInteger("273583686243471418078046036292938539652"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("161917343781254382163963077301550563567"), 128, true),
               BigFixedPoint(BigInteger("311152443074692736824352186014335408156"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("97355433476106157288077219655957918362"), 128, true),
               BigFixedPoint(BigInteger("336318161343136424694695901841775800828"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("29286437520563535328487200882786010496"), 128, true),
               BigFixedPoint(BigInteger("348136658982339264967867138036106806354"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("39518997925195344591398501225959432883"), 128, false),
               BigFixedPoint(BigInteger("346226947218601574485487145087988588362"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("106274150951931075782750588501284147648"), 128, false),
               BigFixedPoint(BigInteger("330793168550893306693876614221731825723"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("168302822651450193478236284374557061220"), 128, false),
               BigFixedPoint(BigInteger("302624017298379579243622911980893710731"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("223172712788269601050257487780680679410"), 128, false),
               BigFixedPoint(BigInteger("263067441411603613475403140531675225854"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("268838277499134743505006047266580888328"), 128, false),
               BigFixedPoint(BigInteger("213972753153500930100328268852044496823"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("303792157693636518457724129356882969974"), 128, false),
               BigFixedPoint(BigInteger("157576849671356512677734982596091307782"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("327191119747388911900252417088460946142"), 128, false),
               BigFixedPoint(BigInteger("96289389574348118233366899415301954126"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("338833016467020703571417641180937488614"), 128, false),
               BigFixedPoint(BigInteger("32362949412358133459608716350368962602"), 128, false)),
};

const auto z_upper_roots_32_scale = 128;

const auto z_upper_roots       = z_upper_roots_32;
const auto z_upper_roots_scale = z_upper_roots_32_scale;
const size_t zN                = 32;

BigCMatrix getZU();

BigCMatrix getZUInverse();

// Now find the primitive roots for x^N+1

// this is zeta_0 to zeta_15 for N=32
// other parts are just conjugates
// we use the zeta notataion as in original CKKS paper
const std::vector<BigComplex> r_roots_32 = {
    BigComplex(BigFixedPoint(BigInteger("338643814315582355912937410983851649519"), 128, false),
               BigFixedPoint(BigInteger("33353504510164655995148692934415085169"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("300102255270364911803307903244229695504"), 128, false),
               BigFixedPoint(BigInteger("160407997365957196655387116530659854223"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("263041826724899848421390153303749767787"), 128, true),
               BigFixedPoint(BigInteger("215872848293952797667059038254236574542"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("325629922445073537381874007153931917509"), 128, false),
               BigFixedPoint(BigInteger("98778757057029163154891549982764432298"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("33353504510164655995148692934415085169"), 128, false),
               BigFixedPoint(BigInteger("338643814315582355912937410983851649519"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("160407997365957196655387116530659854223"), 128, false),
               BigFixedPoint(BigInteger("300102255270364911803307903244229695504"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("215872848293952797667059038254236574542"), 128, false),
               BigFixedPoint(BigInteger("263041826724899848421390153303749767787"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("98778757057029163154891549982764432298"), 128, true),
               BigFixedPoint(BigInteger("325629922445073537381874007153931917509"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("338643814315582355912937410983851649519"), 128, true),
               BigFixedPoint(BigInteger("33353504510164655995148692934415085169"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("300102255270364911803307903244229695504"), 128, true),
               BigFixedPoint(BigInteger("160407997365957196655387116530659854223"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("263041826724899848421390153303749767787"), 128, false),
               BigFixedPoint(BigInteger("215872848293952797667059038254236574542"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("325629922445073537381874007153931917509"), 128, true),
               BigFixedPoint(BigInteger("98778757057029163154891549982764432298"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("33353504510164655995148692934415085169"), 128, true),
               BigFixedPoint(BigInteger("338643814315582355912937410983851649519"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("160407997365957196655387116530659854223"), 128, true),
               BigFixedPoint(BigInteger("300102255270364911803307903244229695504"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("215872848293952797667059038254236574542"), 128, true),
               BigFixedPoint(BigInteger("263041826724899848421390153303749767787"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("98778757057029163154891549982764432298"), 128, false),
               BigFixedPoint(BigInteger("325629922445073537381874007153931917509"), 128, false)),
};

const auto r_roots_32_scale = 128;

const auto r_roots       = r_roots_32;
const auto r_roots_scale = r_roots_32_scale;
const auto rN            = 32;

BigCMatrix getRU();

BigCMatrix getRUInverse();

// Multiplies input vector by U matrix
std::vector<BigComplex> multU(const BigCMatrix& U, std::vector<BigFixedPoint> input);

std::vector<BigFixedPoint> multUInverse(const BigCMatrix& UInv, std::vector<BigComplex> input);

#endif