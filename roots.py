#!/usr/bin/env python3
"""
Compute the roots of X^32 - X + 2 in the upper half-plane using mpmath
with high precision, and print them in the form:

    BigComplex(
        BigFixedPoint(str(c), 64, aNeg),
        BigFixedPoint(str(d), 64, bNeg)
    )

where for a complex root z = a + i b:

    c    = abs(round(a * 2^64))
    d    = abs(round(b * 2^64))
    aNeg = (a < 0)
    bNeg = (b < 0)
"""

from mpmath import *

# Set high precision: 80 decimal digits (~265 bits), > 64-bit fractional precision
mp.dps = 360

def main(n):
    # Polynomial: X^n - X + 2
    # Coefficients for mp.polyroots: highest degree first
    # X^n + 0*X^{n-1} + ... + 0*X^2 - X + 2
    coeffs = [mp.mpf(1)] + [mp.mpf(0)] * (n-2) + [mp.mpf(-1), mp.mpf(2)]

    # Find all roots with high precision
    roots = mp.polyroots(coeffs, maxsteps=400, error=False)

    # Keep only roots with positive imaginary part (upper half-plane).
    # Use a tiny epsilon to avoid numerical noise.
    eps = mp.mpf("1e-40")
    upper_roots = [z for z in roots if mp.im(z) > eps]

    # Sort by real part (then imaginary) for a stable, readable order
    upper_roots.sort(key=lambda z: (mp.re(z), mp.im(z)))

    # 2^64 as an exact integer
    scale_int = 1 << 128

    for z in upper_roots:
        a = mp.re(z)
        b = mp.im(z)

        aNeg = (a < 0)
        bNeg = (b < 0)

        # c = abs(round(a * 2^128)), d similarly
        c = abs(int(mp.nint(a * scale_int)))
        d = abs(int(mp.nint(b * scale_int)))

        aNeg = "true" if aNeg else "false"
        bNeg = "true" if bNeg else "false"

        # Print in the requested format
        # (Python bools print as True/False, which is usually fine for code-gen style)
        print(
            f"BigComplex("
            f"BigFixedPoint(BigInteger(\"{c}\"), 128, {aNeg}), "
            f"BigFixedPoint(BigInteger(\"{d}\"), 128, {bNeg})),"
        )

def main2(n):
    # Polynomial: X^n + 1
    # Coefficients for mp.polyroots: highest degree first
    # X^n + 0*X^{n-1} + ... + 0*X^2 - X + 2
    coeffs = [mp.mpf(1)] + [mp.mpf(0)] * (n-1) + [mp.mpf(1)]
    m = 2 * n
    # 2n-th primitive roots
    omega = exp(2j * pi / m)

    all_roots = [omega]

    # 2^128 as an exact integer
    scale_int = 1 << 128

    for z in all_roots:
        a = mp.re(z)
        b = mp.im(z)

        aNeg = (a < 0)
        bNeg = (b < 0)

        # c = abs(round(a * 2^128)), d similarly
        c = abs(int(mp.nint(a * scale_int)))
        d = abs(int(mp.nint(b * scale_int)))

        aNeg = "true" if aNeg else "false"
        bNeg = "true" if bNeg else "false"

        # Print in the requested format
        # (Python bools print as True/False, which is usually fine for code-gen style)
        print(
            f"const auto R_ROOT_M{m} = "
            f"BigComplex("
            f"BigFixedPoint(BigInteger(\"{c}\"), 128, {aNeg}), "
            f"BigFixedPoint(BigInteger(\"{d}\"), 128, {bNeg}));"
        )

if __name__ == "__main__":
    for deg in [8, 16, 32, 64, 128, 256]:
        print(f"const BigCVector Z_ROOTS_N{deg} = ", "{")
        main(deg)
        print("};")
    #for deg in range(1, 18):
    #    main2(2 ** deg)