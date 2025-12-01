from mpmath import * # pyright: ignore[reportMissingImports]

"""
std::vector<double> EvalChebyshevCoefficients(std::function<double(double)> func, double a, double b, uint32_t degree) {
    if (degree == 0)
        OPENFHE_THROW("The degree of approximation can not be zero");
    // the number of coefficients to be generated should be degree+1 as zero is also included
    size_t coeffTotal{degree + 1};
    double bMinusA = 0.5 * (b - a);
    double bPlusA  = 0.5 * (b + a);
    double PiByDeg = M_PI / static_cast<double>(coeffTotal);
    std::vector<double> functionPoints(coeffTotal);
    for (size_t i = 0; i < coeffTotal; ++i)
        functionPoints[i] = func(std::cos(PiByDeg * (i + 0.5)) * bMinusA + bPlusA);

    double multFactor = 2.0 / static_cast<double>(coeffTotal);
    std::vector<double> coefficients(coeffTotal);
    for (size_t i = 0; i < coeffTotal; ++i) {
        for (size_t j = 0; j < coeffTotal; ++j)
            coefficients[i] += functionPoints[j] * std::cos(PiByDeg * i * (j + 0.5));
        coefficients[i] *= multFactor;
    }
    return coefficients;
}
"""

mp.dps = 100  # Set decimal places for high precision

def EvalChebyshevCoefficients(f, a, b, degree):
    """
    Compute the Chebyshev coefficients for a given function f on the interval [a, b]
    using mpmath for high-precision arithmetic.
    """
    # translate the above C++ code to Python using mpmath
    if degree == 0:
        raise ValueError("The degree of approximation can not be zero")
    coeffTotal = degree + 1
    bMinusA = mp.mpf(0.5) * (mp.mpf(b) - mp.mpf(a))
    bPlusA = mp.mpf(0.5) * (mp.mpf(b) + mp.mpf(a))
    PiByDeg = mp.pi / mp.mpf(coeffTotal)
    functionPoints = [mp.mpf(0)] * coeffTotal
    for i in range(coeffTotal):
        functionPoints[i] = f(mp.cos(PiByDeg * (mp.mpf(i) + mp.mpf(0.5))) * bMinusA + bPlusA)
    multFactor = mp.mpf(2) / mp.mpf(coeffTotal)
    coefficients = [mp.mpf(0)] * coeffTotal
    for i in range(coeffTotal):
        for j in range(coeffTotal):
            coefficients[i] += functionPoints[j] * mp.cos(PiByDeg * mp.mpf(i) * (mp.mpf(j) + mp.mpf(0.5)))
        coefficients[i] *= multFactor
    return coefficients


"""
// A cleartext version of CryptoContext<...>::EvalChebyshevFunction(...)
std::vector<double> EvalChebyshevFunctionPtxt(std::function<double(double)> func, const std::vector<double>& ptxt,
                                              double a, double b, size_t degree) {
    auto coeffs = EvalChebyshevCoefficients(func, a, b, degree);

    // The standard practice is to halve the 1st coefficient.
    // See, for example, the Chebyshev Series section at
    // https://www.cfm.brown.edu/people/dobrush/am34/Mathematica/ch5/chebyshev.html
    // and derivation of Eq. (6) in https://arxiv.org/pdf/1810.04282.
    // The halving requirement follows from the discrete orthogonality relation for Chebyshev polynomials,
    // i.e., Eq. (4) in https://arxiv.org/pdf/1810.04282.
    coeffs[0] /= 2.0;

    // Special case for trivial case of a degee-0 approximation
    if (degree == 0)
        return std::vector<double>(ptxt.size(), coeffs[0]);

    // If [a,b] is different than [-1,1] then need to scale the input
    double scaleFactor = 2.0 / (b - a);
    double offset      = (b + a) * scaleFactor / -2.0;

    std::vector<double> result(ptxt.size());
    for (size_t i = 0; i < ptxt.size(); i++) {
        double x  = ptxt[i] * scaleFactor + offset;
        double x2 = 2 * x;

        double t_prev = 1.0;  // T0(x) = 1
        double t_j    = x;    // T1(x) = x
        double y      = coeffs[0] + coeffs[1] * x;
        // Use the recursive formula T_{i+1}(X) = 2x T_i(x) - T_{i-1}(x)
        for (size_t j = 2; j < coeffs.size(); j++) {
            // Compute T_j(x) and add it to the approximation
            double t_next = x2 * t_j - t_prev;
            t_prev        = t_j;
            t_j           = t_next;
            y += coeffs[j] * t_next;
        }
        result[i] = y;
    }
    return result;
}
"""

def EvalChebyshevFunctionPtxt(f, ptxt, a, b, degree):
    """
    Evaluate the Chebyshev function approximation on plaintext data.
    """
    coeffs = EvalChebyshevCoefficients(f, a, b, degree)
    coeffs[0] /= mp.mpf(2)

    if degree == 0:
        return [coeffs[0]] * len(ptxt)

    scaleFactor = mp.mpf(2) / (mp.mpf(b) - mp.mpf(a))
    offset = (mp.mpf(b) + mp.mpf(a)) * scaleFactor / mp.mpf(-2)

    result = [mp.mpf(0)] * len(ptxt)
    for i in range(len(ptxt)):
        x = mp.mpf(ptxt[i]) * scaleFactor + offset
        x2 = mp.mpf(2) * x

        t_prev = mp.mpf(1)  # T0(x) = 1
        t_j = x             # T1(x) = x
        y = coeffs[0] + coeffs[1] * x
        for j in range(2, len(coeffs)):
            t_next = x2 * t_j - t_prev
            t_prev = t_j
            t_j = t_next
            y += coeffs[j] * t_next
        result[i] = y
    return result

def g0(x):
    """
    auto K = 16.0;
    auto R = 3.0;
    auto f = [&](double x) {
        return (1.0 / std::pow(2 * M_PI, std::pow(2.0, -R))) * std::cos(2 * M_PI * (K * x - 0.25) / std::pow(2.0, R));
    };
    """
    K = mp.mpf(16)
    R = mp.mpf(3)
    return (mp.mpf(1) / mp.power(mp.mpf(2) * mp.pi, mp.power(mp.mpf(2), -R))) * mp.cos(
        mp.mpf(2) * mp.pi * (K * x - mp.mpf(0.25)) / mp.power(mp.mpf(2), R)
    )

def g0_cheby_coeffs():
    print("Computing Chebyshev coefficients for g0...")
    coeffs = EvalChebyshevCoefficients(g0, mp.mpf(-1), mp.mpf(1), 48)
    print("Chebyshev Coefficients for g0:")
    for i in range(len(coeffs)):
        print(f"  coeffs[{i}] = {mp.nstr(coeffs[i], 50)}")

def g0_cheby(x):
    return EvalChebyshevFunctionPtxt(g0, [x], mp.mpf(-1), mp.mpf(1), 48)[0]

# estimate the errors of the approximation on a dense grid
def estimateError(f, g, a, b, num_points=100):
    max_error = mp.mpf(0)
    for i in range(num_points):
        x = a + (b - a) * mp.mpf(i) / mp.mpf(num_points - 1)
        error = mp.fabs(f(x) - g(x))
        if error > max_error:
            max_error = error
    return max_error

def g0_error():
    print("Estimating error for g0 approximation...")
    error = estimateError(g0, g0_cheby, mp.mpf(-1), mp.mpf(1))
    print(f"Maximum error bits in g0 approximation: {mp.log(error) / mp.log(2)}")

#if __name__ == "__main__":
#    g0_cheby_coeffs()
#    g0_error()

def expHalf(x):
    """
    std::exp(1j * Pi/2.0 * x) in [-16, 16] of degree 46
    """
    return mp.exp(mp.mpc(0, mp.pi / mp.mpf(2) * x))

EXP_HALF_DEGREE=58
def expHalf_cheby_coeffs():
    print("Computing Chebyshev coefficients for expHalf...")
    coeffs = EvalChebyshevCoefficients(expHalf, mp.mpf(-16), mp.mpf(16), EXP_HALF_DEGREE)
    print("Chebyshev Coefficients for expHalf:")
    for i in range(len(coeffs)):
        print(f"  coeffs[{i}] = {mp.nstr(coeffs[i], 50)}")

def expHalf_cheby(x):
    return EvalChebyshevFunctionPtxt(expHalf, [x], mp.mpf(-16), mp.mpf(16), EXP_HALF_DEGREE)[0]

def expHalf_error():
    print("Estimating error for expHalf approximation...")
    error = estimateError(expHalf, expHalf_cheby, mp.mpf(-16), mp.mpf(16))
    print(f"Maximum error bits in expHalf approximation: {mp.log(error) / mp.log(2)}")

def expHalf_plot():
    import matplotlib.pyplot as plt
    import numpy as np

    xs = np.linspace(-16, 16, 100)
    exp_half_vals = [expHalf(mp.mpf(x)) for x in xs]
    exp_half_cheby_vals = [expHalf_cheby(mp.mpf(x)) for x in xs]

    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    plt.plot(xs, [mp.re(z) for z in exp_half_vals], label='Re(expHalf)', color='blue')
    plt.plot(xs, [mp.re(z) for z in exp_half_cheby_vals], label='Re(expHalf_cheby)', linestyle='dashed', color='orange')
    plt.title('Real Part of expHalf vs Chebyshev Approximation')
    plt.xlabel('x')
    plt.ylabel('Real Part')
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(xs, [mp.im(z) for z in exp_half_vals], label='Im(expHalf)', color='blue')
    plt.plot(xs, [mp.im(z) for z in exp_half_cheby_vals], label='Im(expHalf_cheby)', linestyle='dashed', color='orange')
    plt.title('Imaginary Part of expHalf vs Chebyshev Approximation')
    plt.xlabel('x')
    plt.ylabel('Imaginary Part')
    plt.legend()

    plt.tight_layout()
    #plt.show()
    plt.savefig("expHalf_cheby_comparison.png")

if __name__ == "__main__":
    #expHalf_cheby_coeffs()
    expHalf_error()
    #expHalf_plot()

def printC(z):
    a = mp.re(z)
    b = mp.im(z)

    aNeg = (a < 0)
    bNeg = (b < 0)

    # c = abs(round(a * 2^128)), d similarly
    scale_int = 1 << 128
    c = abs(int(mp.nint(a * scale_int)))
    d = abs(int(mp.nint(b * scale_int)))

    aNeg = "true" if aNeg else "false"
    bNeg = "true" if bNeg else "false"

    # Print in the requested format
    # (Python bools print as True/False, which is usually fine for code-gen style)
    print(f"BigComplex(BigFixedPoint(BigInteger(\"{c}\"), 128, {aNeg}), BigFixedPoint(BigInteger(\"{d}\"), 128, {bNeg})),")

def printExpHalf():
    coeffs = EvalChebyshevCoefficients(expHalf, mp.mpf(-16), mp.mpf(16), EXP_HALF_DEGREE)
    print(f"coeff_exp_16_big_complex_{EXP_HALF_DEGREE} = ", "{")
    for coeff in coeffs:
        printC(coeff)
    print("};")
#
#if __name__ == "__main__":
#    printExpHalf()

def expFull(x):
    """
    std::exp(1j * 2 * Pi * x) in [-16, 16]
    """
    return mp.exp(mp.mpc(0, mp.pi * mp.mpf(2) * x))

EXP_FULL_DEGREE=58
def expFull_cheby_coeffs():
    print("Computing Chebyshev coefficients for expFull...")
    coeffs = EvalChebyshevCoefficients(expFull, mp.mpf(-16), mp.mpf(16), EXP_FULL_DEGREE)
    print("Chebyshev Coefficients for expFull:")
    for i in range(len(coeffs)):
        print(f"  coeffs[{i}] = {mp.nstr(coeffs[i], 50)}")
def expFull_cheby(x):
    return EvalChebyshevFunctionPtxt(expFull, [x], mp.mpf(-16), mp.mpf(16), EXP_FULL_DEGREE)[0]
def expFull_error():
    print("Estimating error for expFull approximation...")
    error = estimateError(expFull, expFull_cheby, mp.mpf(-16), mp.mpf(16))
    print(f"Maximum error bits in expFull approximation: {mp.log(error) / mp.log(2)}")

def expFull_plot():
    import matplotlib.pyplot as plt
    import numpy as np

    xs = np.linspace(-16, 16, 100)
    exp_full_vals = [expFull(mp.mpf(x)) for x in xs]
    exp_full_cheby_vals = [expFull_cheby(mp.mpf(x)) for x in xs]

    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    plt.plot(xs, [mp.re(z) for z in exp_full_vals], label='Re(expFull)', color='blue')
    plt.plot(xs, [mp.re(z) for z in exp_full_cheby_vals], label='Re(expFull_cheby)', linestyle='dashed', color='orange')
    plt.title('Real Part of expFull vs Chebyshev Approximation')
    plt.xlabel('x')
    plt.ylabel('Real Part')
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(xs, [mp.im(z) for z in exp_full_vals], label='Im(expFull)', color='blue')
    plt.plot(xs, [mp.im(z) for z in exp_full_cheby_vals], label='Im(expFull_cheby)', linestyle='dashed', color='orange')
    plt.title('Imaginary Part of expFull vs Chebyshev Approximation')
    plt.xlabel('x')
    plt.ylabel('Imaginary Part')
    plt.legend()

    plt.tight_layout()
    #plt.show()
    plt.savefig("expFull_cheby_comparison.png")

if __name__ == "__main__":
    #expFull_cheby_coeffs()
    #expFull_error()
    #expFull_plot()
    pass