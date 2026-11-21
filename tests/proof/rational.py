#!/usr/bin/env python3
"""Proves what source/rational.cpp's bounds on a strip's integrals rest on, that a computer can.

    python3 -I tests/proof/rational.py

- The pairing identity. For Bernstein polynomials of degree n and i < j,

      B_i B_j' - B_i' B_j = (j - i) C(n, i) C(n, j) / C(2 n - 2, i + j - 1) B_(i+j-1)^(2 n - 2),

  which Kernel::Derived works out a piece's W X' - W' X by and Fanning its X x X': a polynomial
  identity of degree 2 n - 2, so that it holds once it holds at 2 n - 1 points. Both sides are
  times C(2 n - 2, i + j - 1) to clear the division, and taken at the integers 2 to 2 n, where
  every term is an integer: every pair of every degree to MaxDegree is checked exactly.

- The Bernstein ellipse's bound on a polynomial. On the ellipse of rho about a piece of [0, 1]
  with middle c and half axis R, |B_i(z)| summed is (|z| + |1 - z|)^n by the binomial theorem, and
  |z| + |1 - z| is at most |c| + |1 - c| + 2 R, which is max(1, |2 c - 1|) + 2 R: each of |z| and
  |1 - z| is within R of |c| and |1 - c| by the triangle inequality, every point of the ellipse
  being within its half axis of its middle. Checked here on a grid of pieces, rhos and points of
  their ellipses, in exact rational arithmetic, against |z| and |1 - z| bounded from above and the
  bound from below; a grid is no proof, the two inequalities above are.

- How far a root lies from an ellipse inside its own. A root r lies on the ellipse about the piece
  whose foci sum to 2 A_r, A_r the half sum of its distances from the piece's ends; a point z of a
  confocal ellipse of half axis R < A_r has |r - f_k| <= |r - z| + |z - f_k| for each focus, and
  summed, 2 A_r <= 2 |r - z| + 2 R: |r - z| >= A_r - R. No more is needed; it is stated here as
  the bound on the weight from below rests on it.

- The constants: Close to the MostCuts stays a normal double, and the stack of pieces holds every
  cut and every halving: the cuts at most the roots of the speed's numerator, 2 MaxDegree - 2, and
  of the weight, MaxDegree, and a piece more per halving, at most Halved.

The exit status is 1 when a check fails.
"""

import pathlib
import re
import sys
from fractions import Fraction
from math import comb

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def constant(text, name):
    """A constant's value as the source writes it, a hexadecimal float or an integer."""
    found = re.search(r"constexpr \w+(?:::\w+)?\s+" + name + r"\s*=\s*([^;]+);", text)
    value = found[1].strip()
    if value.startswith("0x"):
        return float.fromhex(value)
    return value


def bernstein(n, i, t):
    """B_i^n(t), exactly, t an integer."""
    return comb(n, i) * t**i * (1 - t) ** (n - i)


def slope(n, i, t):
    """B_i^n'(t), exactly, t an integer."""
    value = 0
    if i > 0:
        value += comb(n, i) * i * t ** (i - 1) * (1 - t) ** (n - i)
    if i < n:
        value -= comb(n, i) * (n - i) * t**i * (1 - t) ** (n - i - 1)
    return value


def paired(most):
    """Whether the pairing identity holds for every degree from 1 to the most and every pair."""
    for n in range(1, most + 1):
        twice = 2 * n - 2
        points = range(2, 2 * n + 1)
        for i in range(n):
            for j in range(i + 1, n + 1):
                below = comb(twice, i + j - 1)
                scale = (j - i) * comb(n, i) * comb(n, j)
                for t in points:
                    lhs = below * (bernstein(n, i, t) * slope(n, j, t)
                                   - slope(n, i, t) * bernstein(n, j, t))
                    rhs = scale * bernstein(twice, i + j - 1, t)
                    if lhs != rhs:
                        return False
    return True


def bounded():
    """Whether |z| + |1 - z| <= max(1, |2 c - 1|) + 2 R on a grid of pieces' ellipses."""
    for start, end in ((Fraction(0), Fraction(1)), (Fraction(1, 3), Fraction(1, 2)),
                       (Fraction(-2), Fraction(-1)), (Fraction(3, 4), Fraction(5, 2))):
        half = (end - start) / 2
        middle = start + half
        for rho in (Fraction(5, 4), Fraction(3), Fraction(64)):
            reach = half * (rho + 1 / rho) / 2
            minor = half * (rho - 1 / rho) / 2
            bound = max(Fraction(1), abs(2 * middle - 1)) + 2 * reach
            for k in range(16):
                # A point of the ellipse at rational cosine and sine (a Pythagorean pair).
                cosine = Fraction(k * k - 1, k * k + 1) if k else Fraction(-1)
                sine = Fraction(2 * k, k * k + 1) if k else Fraction(0)
                x = middle + reach * cosine
                y = minor * sine
                near = (x * x + y * y)
                far = ((1 - x) ** 2 + y * y)
                # |z| + |1 - z| <= bound iff, bound being positive, the squares obey it.
                if not (near <= bound**2 and far <= bound**2):
                    return False
                rest = bound**2 - near - far
                if rest < 0 or 4 * near * far > rest**2:
                    return False
    return True


def main():
    failed = False

    def report(name, held, told=""):
        nonlocal failed
        failed = failed or not held
        print(f"{name:60} {'held' if held else 'BROKEN':7} {told}", flush=True)

    spline = (ROOT / "include" / "geoxy" / "spline.hpp").read_text()
    rational = (ROOT / "source" / "rational.cpp").read_text()
    hodograph = (ROOT / "source" / "hodograph.hpp").read_text()
    terms = (ROOT / "source" / "kernel" / "piece.hpp").read_text()
    most = int(constant(spline, "MaxDegree"))
    report(f"pairing identity, every degree to {most}", paired(most))
    report("|z| + |1 - z| on a piece's ellipse within max(1, |2 c - 1|) + 2 R", bounded())
    close = constant(rational, "Close")
    deepest = int(constant(rational, "MostCuts"))
    report(f"Close^{deepest} a normal double", close**deepest >= 2.0**-1022, f"{close**deepest!r}")
    halved = int(constant(hodograph, "Halved"))
    kernel = 3 * most
    report("the pieces' stack holds every cut and halving",
           halved + kernel + 4 >= (2 * most - 2) + most + 1 + halved,
           f"{halved + kernel + 4} for {(2 * most - 2) + most + 1 + halved}")
    if "Terms = 3 * static_cast<std::size_t>(MaxDegree)" not in terms:
        report("Kernel::Terms is 3 MaxDegree", False)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
