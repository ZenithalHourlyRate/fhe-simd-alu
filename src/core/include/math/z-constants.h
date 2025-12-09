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

const auto r_root_32 = BigComplex(BigFixedPoint(BigInteger("338643814315582355912937410983851649519"), 128, false),
                                  BigFixedPoint(BigInteger("33353504510164655995148692934415085169"), 128, false));
const auto r_root_64 = BigComplex(BigFixedPoint(BigInteger("339872481907374595982655373237523399998"), 128, false),
                                  BigFixedPoint(BigInteger("16696864359439569162511425913059001623"), 128, false));

const auto r_roots_32_scale = 128;

const auto r_roots       = r_roots_32;
const auto r_roots_scale = r_roots_32_scale;
const auto rN            = 32;

BigCMatrix getRU();

BigCMatrix getRUInverse();

// Multiplies input vector by U matrix
std::vector<BigComplex> multU(const BigCMatrix& U, std::vector<BigFixedPoint> input);

std::vector<BigFixedPoint> multUInverse(const BigCMatrix& UInv, std::vector<BigComplex> input);

// Primitive roots of unity for various m
// Here we define exp(2 pi i / m) for m = 8, 16, ..., 262144
// They are roots for X^{m/2} + 1, and offers m/4 C slots
// The minimal m is 8 because the minimal C slot we have is 2 (sparse packing of 1 complex number)
// If we account for Z slot, minimal m is much larger depending on the choice
const auto R_ROOT_M8      = BigComplex(BigFixedPoint(BigInteger("240615969168004511545033772477625056927"), 128, false),
                                       BigFixedPoint(BigInteger("240615969168004511545033772477625056927"), 128, false));
const auto R_ROOT_M16     = BigComplex(BigFixedPoint(BigInteger("314379914072750776175968446001973125505"), 128, false),
                                       BigFixedPoint(BigInteger("130220424146621615521356287377630227996"), 128, false));
const auto R_ROOT_M32     = BigComplex(BigFixedPoint(BigInteger("333743936656827580393945988193699131951"), 128, false),
                                       BigFixedPoint(BigInteger("66385796539016198537004233668372411924"), 128, false));
const auto R_ROOT_M64     = BigComplex(BigFixedPoint(BigInteger("338643814315582355912937410983851649519"), 128, false),
                                       BigFixedPoint(BigInteger("33353504510164655995148692934415085169"), 128, false));
const auto R_ROOT_M128    = BigComplex(BigFixedPoint(BigInteger("339872481907374595982655373237523399998"), 128, false),
                                       BigFixedPoint(BigInteger("16696864359439569162511425913059001623"), 128, false));
const auto R_ROOT_M256    = BigComplex(BigFixedPoint(BigInteger("340179880234010501256410555240815138610"), 128, false),
                                       BigFixedPoint(BigInteger("8350947328924239869306493724758181497"), 128, false));
const auto R_ROOT_M512    = BigComplex(BigFixedPoint(BigInteger("340256744284537774291616451567166656226"), 128, false),
                                       BigFixedPoint(BigInteger("4175788093625692072829030064408394408"), 128, false));
const auto R_ROOT_M1024   = BigComplex(BigFixedPoint(BigInteger("340275961201545362514682404085958700136"), 128, false),
                                       BigFixedPoint(BigInteger("2087933351568136955647909863034040375"), 128, false));
const auto R_ROOT_M2048   = BigComplex(BigFixedPoint(BigInteger("340280765487321862450779202347741264351"), 128, false),
                                       BigFixedPoint(BigInteger("1043971588913162969962874901364098823"), 128, false));
const auto R_ROOT_M4096   = BigComplex(BigFixedPoint(BigInteger("340281966562298792572160534645360975660"), 128, false),
                                       BigFixedPoint(BigInteger("521986408598802148050493966777353034"), 128, false));
const auto R_ROOT_M8192   = BigComplex(BigFixedPoint(BigInteger("340282266831263825696362301831656808750"), 128, false),
                                       BigFixedPoint(BigInteger("260993281067212527299757014848791288"), 128, false));
const auto R_ROOT_M16384  = BigComplex(BigFixedPoint(BigInteger("340282341898518884018790830072534181986"), 128, false),
                                       BigFixedPoint(BigInteger("130496650129583753759140694874978295"), 128, false));
const auto R_ROOT_M32768  = BigComplex(BigFixedPoint(BigInteger("340282360665333511102050687223426738058"), 128, false),
                                       BigFixedPoint(BigInteger("65248326264289096219790833612199600"), 128, false));
const auto R_ROOT_M65536  = BigComplex(BigFixedPoint(BigInteger("340282365357037221779282487371359009126"), 128, false),
                                       BigFixedPoint(BigInteger("32624163282081701561065576602700777"), 128, false));
const auto R_ROOT_M131072 = BigComplex(BigFixedPoint(BigInteger("340282366529963152817741505908074724901"), 128, false),
                                       BigFixedPoint(BigInteger("16312081659782994994230389606789418"), 128, false));
const auto R_ROOT_M262144 = BigComplex(BigFixedPoint(BigInteger("340282366823194635787928202577525532155"), 128, false),
                                       BigFixedPoint(BigInteger("8156040832234265524836811571533929"), 128, false));

const std::map<uint32_t, BigComplex> R_ROOT_MAP = {
    {8, R_ROOT_M8},         {16, R_ROOT_M16},       {32, R_ROOT_M32},         {64, R_ROOT_M64},
    {128, R_ROOT_M128},     {256, R_ROOT_M256},     {512, R_ROOT_M512},       {1024, R_ROOT_M1024},
    {2048, R_ROOT_M2048},   {4096, R_ROOT_M4096},   {8192, R_ROOT_M8192},     {16384, R_ROOT_M16384},
    {32768, R_ROOT_M32768}, {65536, R_ROOT_M65536}, {131072, R_ROOT_M131072}, {262144, R_ROOT_M262144},
};

#endif