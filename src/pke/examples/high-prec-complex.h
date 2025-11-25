
#include <cassert>
#include "openfhe.h"

using namespace lbcrypto;

struct BigFixedPoint {
public:
    BigFixedPoint() : value(0), log2Scale(0), neg(false) {}
    BigFixedPoint(const BigInteger& val, int log2S, bool neg) : value(val), log2Scale(log2S), neg(neg) {}
    BigFixedPoint(int a) {
        log2Scale = 0;
        if (a < 0) {
            value = -a;
            neg   = true;
        }
        else {
            value = a;
            neg   = false;
        }
    }
    explicit BigFixedPoint(double a) {
        neg = false;
        if (a < 0) {
            neg = true;
            a   = -a;
        }
        // find log2 scale to represent fractional part
        log2Scale = 0;
        double intpart;
        double fracpart = std::modf(a, &intpart);
        while (fracpart != 0.0 && log2Scale < 128) {
            a *= 2.0;
            fracpart = std::modf(a, &intpart);
            log2Scale++;
        }
        value = BigInteger(static_cast<int64_t>(std::round(a)));
    }

    BigInteger getValue() const {
        return value;
    }
    int getLog2Scale() const {
        return log2Scale;
    }
    bool getNeg() const {
        return neg;
    }

    BigFixedPoint scaleUp(int log2S) const {
        return BigFixedPoint(getValue() << log2S, getLog2Scale() + log2S, getNeg());
    }
    BigFixedPoint scaleDown(int log2S) const {
        return BigFixedPoint(getValue() >> log2S, getLog2Scale() - log2S, getNeg());
    }
    BigFixedPoint scaleTo(int log2S) const {
        if (log2S > log2Scale) {
            return scaleUp(log2S - log2Scale);
        }
        if (log2S < log2Scale) {
            return scaleDown(log2Scale - log2S);
        }
        return *this;
    }

    double convertToDouble() const {
        double val = value.ConvertToDouble();
        val /= std::pow(2.0, log2Scale);
        if (neg) {
            val = -val;
        }
        return val;
    }

    std::string toString(int prec = 32) const {
        std::string result;

        if (neg && value != 0) {
            result += "-";
        }

        // Calculate the integer and fractional parts
        BigInteger scaledValue    = value;
        BigInteger integerPart    = scaledValue >> log2Scale;                  // Divide by 2^log2Scale
        BigInteger fractionalPart = scaledValue - (integerPart << log2Scale);  // Remainder

        // Convert integer part to string
        result += integerPart.ToString();

        // Convert fractional part to decimal
        if (fractionalPart != 0 && log2Scale > 0) {
            result += ".";

            // Convert fractional part (which is in base-2) to base-10
            BigInteger frac            = fractionalPart;
            const int maxDecimalDigits = prec;  // Maximum decimal digits to display

            for (int i = 0; i < maxDecimalDigits && frac != 0; i++) {
                frac *= 10;                            // Multiply by 10 to get next decimal digit
                BigInteger digit = frac >> log2Scale;  // Extract the digit
                result += digit.ToString();
                frac = frac - (digit << log2Scale);  // Remainder for next iteration

                // Early exit if remainder becomes zero
                if (frac == 0) {
                    break;
                }
            }
        }
        else if (log2Scale > 0) {
            // Add .0 if there's fractional precision but value is integer
            result += ".0";
        }

        return result;
    }

    std::string toBinary(int prec = 64) const {
        std::string result;

        if (neg && value != 0) {
            result += "-";
        }

        if (value == 0) {
            return "0";
        }

        // Convert integer part to binary manually
        BigInteger integerPart = value >> log2Scale;
        if (integerPart != 0) {
            // Manual binary conversion for integer part
            BigInteger temp = integerPart;
            std::string binaryStr;
            while (temp != 0) {
                // OpenFHE count from 1???
                if (temp.GetBitAtIndex(1)) {
                    binaryStr = "1" + binaryStr;
                }
                else {
                    binaryStr = "0" + binaryStr;
                }
                temp >>= 1;
            }
            result += binaryStr;
        }
        else {
            result += "0";
        }

        // Convert fractional part
        if (log2Scale > 0) {
            BigInteger fractionalPart = value - (integerPart << log2Scale);

            if (fractionalPart != 0) {
                result += ".";

                BigInteger frac = fractionalPart;
                for (int i = 0; i < log2Scale && i < prec && frac != 0; i++) {
                    frac <<= 1;  // Shift left to check next bit
                    if (frac >= (BigInteger(1) << log2Scale)) {
                        result += "1";
                        frac = frac - (BigInteger(1) << log2Scale);
                    }
                    else {
                        result += "0";
                    }

                    if (frac == 0) {
                        break;
                    }
                }
            }
            else {
                result += ".0";
            }
        }
        return result;
    }

private:
    BigInteger value;  // unsigned big
    int log2Scale;     // scale = 2^log2Scale
    bool neg;
    // we actually encode neg * value * 2^{-log2Scale}
};

BigFixedPoint operator+(const BigFixedPoint& a, const BigFixedPoint& b) {
    if (a.getLog2Scale() == b.getLog2Scale()) {
        auto aValue = a.getValue();
        auto bValue = b.getValue();
        auto aNeg   = a.getNeg();
        auto bNeg   = b.getNeg();
        if (aNeg == bNeg) {
            return BigFixedPoint(aValue + bValue, a.getLog2Scale(), aNeg);
        }
        else {
            if (aValue > bValue) {
                return BigFixedPoint(aValue - bValue, a.getLog2Scale(), aNeg);
            }
            else if (bValue > aValue) {
                return BigFixedPoint(bValue - aValue, a.getLog2Scale(), bNeg);
            }
            else {
                // equal values, return zero
                return BigFixedPoint(BigInteger(0), a.getLog2Scale(), false);
            }
        }
    }
    else if (a.getLog2Scale() > b.getLog2Scale()) {
        int diff           = a.getLog2Scale() - b.getLog2Scale();
        BigInteger scaledB = b.getValue() << diff;
        BigFixedPoint bigScaledB(scaledB, a.getLog2Scale(), b.getNeg());
        return a + bigScaledB;
    }
    else {
        int diff           = b.getLog2Scale() - a.getLog2Scale();
        BigInteger scaledA = a.getValue() << diff;
        BigFixedPoint bigScaledA(scaledA, b.getLog2Scale(), a.getNeg());
        return bigScaledA + b;
    }
}

BigFixedPoint operator-(const BigFixedPoint& a, const BigFixedPoint& b) {
    BigFixedPoint negB = b;
    negB               = BigFixedPoint(negB.getValue(), negB.getLog2Scale(), !negB.getNeg());
    return a + negB;
}

// make it scale preserving
BigFixedPoint operator*(const BigFixedPoint& a, const BigFixedPoint& b) {
    BigInteger prod = a.getValue() * b.getValue();
    int log2S       = a.getLog2Scale() + b.getLog2Scale();
    bool neg        = (a.getNeg() != b.getNeg());
    return BigFixedPoint(prod, log2S, neg).scaleTo(std::max(a.getLog2Scale(), b.getLog2Scale()));
}

BigFixedPoint operator/(const BigFixedPoint& a, const BigFixedPoint& b) {
    // Scale up the dividend to preserve precision
    BigInteger scaledA = a.getValue() << b.getLog2Scale();
    BigInteger quot    = scaledA / b.getValue();
    int log2S          = a.getLog2Scale();
    bool neg           = (a.getNeg() != b.getNeg());
    return BigFixedPoint(quot, log2S, neg);
}

struct BigComplex {
public:
    BigComplex() : real(), imag() {}
    BigComplex(const BigFixedPoint& re, const BigFixedPoint& im) : real(re), imag(im) {}
    BigComplex(const BigFixedPoint& re) : real(re), imag(BigFixedPoint(0, 0, false).scaleTo(re.getLog2Scale())) {}
    BigFixedPoint getReal() const {
        return real;
    }
    BigFixedPoint getImag() const {
        return imag;
    }

    BigComplex scaleTo(int log2S) const {
        return BigComplex(real.scaleTo(log2S), imag.scaleTo(log2S));
    }
    BigComplex conj() const {
        return BigComplex(real, BigFixedPoint(imag.getValue(), imag.getLog2Scale(), !imag.getNeg()));
    }

    std::complex<double> convertToComplex() const {
        return {real.convertToDouble(), imag.convertToDouble()};
    }

    std::string toString(int prec = 32) const {
        return "(" + real.toString(prec) + ", " + imag.toString(prec) + ")";
    }

    std::string toBinary(int prec = 64) const {
        return "(" + real.toBinary(prec) + ", " + imag.toBinary(prec) + ")";
    }

private:
    BigFixedPoint real;
    BigFixedPoint imag;
};
BigComplex operator+(const BigComplex& a, const BigComplex& b) {
    return BigComplex(a.getReal() + b.getReal(), a.getImag() + b.getImag());
}
BigComplex operator-(const BigComplex& a, const BigComplex& b) {
    return BigComplex(a.getReal() - b.getReal(), a.getImag() - b.getImag());
}
BigComplex operator*(const BigComplex& a, const BigComplex& b) {
    BigFixedPoint rePart = (a.getReal() * b.getReal()) - (a.getImag() * b.getImag());
    BigFixedPoint imPart = (a.getReal() * b.getImag()) + (a.getImag() * b.getReal());
    return BigComplex(rePart, imPart);
}
BigComplex operator/(const BigComplex& a, const BigComplex& b) {
    BigFixedPoint denom  = (b.getReal() * b.getReal()) + (b.getImag() * b.getImag());
    BigFixedPoint rePart = ((a.getReal() * b.getReal()) + (a.getImag() * b.getImag())) / denom;
    BigFixedPoint imPart = ((a.getImag() * b.getReal()) - (a.getReal() * b.getImag())) / denom;
    return BigComplex(rePart, imPart);
}

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

using BigCMatrix = std::vector<std::vector<BigComplex>>;

BigCMatrix getZU() {
    // Vandermond matrix
    BigCMatrix zu(zN / 2, std::vector<BigComplex>(zN));
    for (size_t i = 0; i != zN / 2; ++i) {
        zu[i][0] = BigFixedPoint(1, 0, false);
        for (size_t j = 1; j != zN; ++j) {
            zu[i][j] = zu[i][j - 1] * z_upper_roots[i];
        }
    }
    return zu;
}

BigCMatrix getZUInverse() {
    std::vector<std::vector<BigComplex>> zUInv(zN, std::vector<BigComplex>(zN / 2));

    // Build the inverse Vandermond matrix
    for (size_t j = 0; j != zN / 2; ++j) {
        auto xi = z_upper_roots[j];
        std::vector<BigComplex> xiPowers;
        auto one   = BigFixedPoint(1, 0, false);
        auto zNbig = BigFixedPoint(zN, 0, false);
        xiPowers.push_back(one);
        for (size_t p = 1; p != zN; ++p) {
            xiPowers.push_back(xiPowers[p - 1] * xi);
        }
        // build power map
        for (size_t i = 0; i != zN; ++i) {
            // d = zN * xi^(zN-1) - 1
            auto d = (zNbig * xiPowers[zN - 1]) - one;
            if (i == 0) {
                zUInv[i][j] = (xiPowers[zN - 1] - one) / d;
            }
            else {
                zUInv[i][j] = xiPowers[zN - 1 - i] / d;
            }
        }
    }
    return zUInv;
}

std::vector<BigComplex> multU(const BigCMatrix& U, std::vector<BigFixedPoint> input) {
    assert(input.size() == U[0].size() && "Input size does not match expected size for U");

    // This is the vandermond matrix
    std::vector<BigComplex> result;

    size_t halfSize = input.size() / 2;
    for (size_t i = 0; i < halfSize; ++i) {
        BigComplex sum;
        for (size_t j = 0; j < input.size(); ++j) {
            sum = sum + U[i][j] * input[j];
        }
        result.push_back(sum);
    }
    return result;
}

std::vector<BigFixedPoint> multUInverse(const BigCMatrix& UInv, std::vector<BigComplex> input) {
    std::vector<BigFixedPoint> result;
    assert(input.size() == z_upper_roots.size() && "Input size does not match expected size for multiplyByZUInverse");

    for (size_t i = 0; i != zN; ++i) {
        BigComplex sum;
        for (size_t j = 0; j != zN / 2; ++j) {
            sum = sum + UInv[i][j] * input[j];
        }
        // z + conj(z) = 2*real(z)
        auto two = BigFixedPoint(2, 0, false);
        result.push_back(two * sum.getReal());
    }
    return result;
}

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

BigCMatrix getRU() {
    std::vector<std::vector<BigComplex>> U(rN / 2, std::vector<BigComplex>(rN));
    for (uint32_t i = 0; i < rN / 2; ++i) {
        auto zeta = r_roots[i];
        // build power map
        std::vector<BigComplex> zetaPows;
        auto one = BigFixedPoint(1, 0, false);
        zetaPows.push_back(one);
        for (uint32_t p = 1; p != rN; ++p) {
            zetaPows.push_back(zetaPows[p - 1] * zeta);
        }

        for (uint32_t j = 0; j < rN; ++j) {
            U[i][j] = zetaPows[j];
        }
    }
    return U;
}

BigCMatrix getRUInverse() {
    auto U = getRU();
    BigCMatrix UInverse(rN, std::vector<BigComplex>(rN / 2));
    for (size_t i = 0; i != rN / 2; ++i) {
        auto rNBig = BigFixedPoint(rN, 0, false);
        for (size_t j = 0; j != rN; ++j) {
            UInverse[j][i] = U[i][j].conj() / rNBig;
        }
    }
    return UInverse;
}