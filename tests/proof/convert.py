#!/usr/bin/env python3
"""Proves the degree and radian conversions correctly rounded for every double.

    python3 -I tests/proof/convert.py

Angle::Converted (source/angle.cpp) multiplies a double by C = pi / 180 or 180 / pi, held in three
doubles, and decides the rounding by an exact sign: its answer is within 2^-100 ulp of the exact
product. That is the nearest double unless the exact product comes within 2^-100 ulp of a halfway
point between doubles. This shows none does, nor within 2^-61 ulp, in exact whole-number and
rational arithmetic:

- x = X 2^e with X a whole number below 2^53, and x C / ulp is X c or X c / 2, c the significand of
  C in [1, 2); on the subnormal grid it is X c 2^t. A halfway point is X a an odd whole number,
  a = 2c, c or c 2^(t+1).
- By Lagrange's theorem on best approximations, |X a - n| >= |q a - p| for every whole n and every
  X below q', where p / q is the last convergent of a's continued fraction with q <= 2^53 and q'
  the next denominator. So no X below 2^53 brings X a nearer a whole number, odd or not.
- pi comes from Machin's formula, to 2^-1190; the convergents up to 2^53 and the bounds need far
  less.

It also checks source/degrees.hpp's three doubles for each C are the nearest splits of C. The exit
status is 1 when a check fails.
"""

import math
import pathlib
import re
import sys
from fractions import Fraction

HERE = pathlib.Path(__file__).resolve().parent
BITS = 1200
DECIDED = Fraction(1, 2**100)
MARGIN = Fraction(1, 2**61)


def pi_scaled(bits):
    """floor(pi 2^bits) or one below, by Machin's formula in whole numbers."""
    one = 1 << (bits + 16)

    def arctan_inverse(n):
        total, term, k, sign = 0, one // n, 1, 1
        while term:
            total += sign * (term // k)
            term //= n * n
            k += 2
            sign = -sign
        return total
    return (16 * arctan_inverse(5) - 4 * arctan_inverse(239)) >> 16


def nearest_double(value):
    """The double nearest a nonzero rational in the normal range, as an exact rational."""
    exponent = math.floor(math.log2(abs(value)))
    while Fraction(2) ** exponent > abs(value):
        exponent -= 1
    while Fraction(2) ** (exponent + 1) <= abs(value):
        exponent += 1
    ulp = Fraction(2) ** (exponent - 52)
    scaled = value / ulp
    whole = math.floor(scaled)
    if scaled - whole > Fraction(1, 2) or (scaled - whole == Fraction(1, 2) and whole % 2):
        whole += 1
    return whole * ulp


def bound(alpha, limit):
    """|q alpha - p| for the last convergent p / q of alpha with q <= limit."""
    h0, h1, k0, k1 = 0, 1, 1, 0
    x = alpha
    while True:
        a = math.floor(x)
        h0, h1 = h1, a * h1 + h0
        k0, k1 = k1, a * k1 + k0
        if k1 > limit:
            return abs(k0 * alpha - h0)
        x = 1 / (x - a)


def significand(value):
    """value's significand, in [1, 2)."""
    while value >= 2:
        value /= 2
    while value < 1:
        value *= 2
    return value


def main():
    failed = False

    def report(name, held, told=""):
        nonlocal failed
        failed = failed or not held
        print(f"{name:52} {'held' if held else 'BROKEN':7} {told}", flush=True)

    pi = Fraction(pi_scaled(BITS), 1 << BITS)
    table = (HERE.parent.parent / "source" / "degrees.hpp").read_text()
    for name, constant in (("ToRadians", pi / 180), ("ToDegrees", 180 / pi)):
        pair = re.search(name + r" = \{\.High = (\S+),\s*\.Low\s*= (\S+)\}", table)
        third = re.search(name + r"Third = (\S+);", table)
        parts = [float.fromhex(pair[1]), float.fromhex(pair[2].rstrip(",")), float.fromhex(third[1])]
        left = constant
        nearest = True
        for part in parts:
            nearest = nearest and Fraction(part) == nearest_double(left)
            left -= Fraction(part)
        report(f"{name}: the three doubles are C's nearest splits", nearest)
        report(f"{name}: they hold C to 2^-150 of it", abs(left) < constant / 2**150,
               f"2^{math.log2(abs(left) / constant):.1f}")

        c = significand(constant)
        alphas = [2 * c, c] + [c * Fraction(2) ** t for t in range(-60, 1)]
        nearest_approach = min(bound(alpha, 1 << 53) for alpha in alphas) / 2
        report(f"{name}: no product within 2^-61 ulp of halfway", nearest_approach > MARGIN,
               f"2^{math.log2(nearest_approach):.2f} ulp")
        report(f"{name}: the 2^-100 ulp decision is on the right side",
               nearest_approach > DECIDED)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
