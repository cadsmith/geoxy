#!/usr/bin/env python3
"""Proves source/legendre.hpp's Gauss-Legendre rules: each node and weight the nearest double.

    python3 -I tests/proof/legendre.py

The rule of n = 2 h points over -1 to 1 has the roots of the Legendre polynomial P_n for nodes,
h of them positive, and weights w = 2 (1 - x^2) / (n P_(n-1)(x))^2 at each root x, P_n being 0
there. In exact rational arithmetic, by P_0 = 1, P_1 = x and j P_j = (2 j - 1) x P_(j-1) -
(j - 1) P_(j-2):

- A node is the nearest double to a root when P_n changes sign across the doubles' halfway
  points on either side of it. The h nodes of a rule are distinct, so their h halfway intervals
  are apart, and P_n has h positive roots: each interval holds exactly one, and each node is the
  nearest double to its own root.
- The root is then bisected to within 2^-160. On that interval |P_(n-1)'| is at most (n - 1)^2,
  by Markov's inequality, |P_(n-1)| being at most 1 on -1 to 1, so P_(n-1) is known there to
  within (n - 1)^2 times the interval's width; 1 - x^2 is known from the interval's ends. The
  weight lies between the bounds that gives, and both must lie inside the tabulated weight's
  halfway interval.

The exit status is 1 when a check fails.
"""

import math
import pathlib
import re
import sys
from fractions import Fraction

HERE = pathlib.Path(__file__).resolve().parent
PAIRS = 20
WIDTH = Fraction(1, 2**160)


def legendre(degree, x):
    """P_(degree - 1)(x) and P_degree(x), exactly."""
    previous, current = Fraction(1), x
    for j in range(2, degree + 1):
        previous, current = current, ((2 * j - 1) * x * current - (j - 1) * previous) / j
    return previous, current


def halfway(value):
    """The halfway points between a positive normal double and its neighbours."""
    exponent = math.frexp(value)[1] - 1
    above = Fraction(2) ** (exponent - 52)
    below = above / 2 if value == 2.0**exponent else above
    return Fraction(value) - below / 2, Fraction(value) + above / 2


def table(text, name):
    """An array of doubles read from the header, std::numbers::inv_sqrt3 the double nearest
    1 / sqrt(3), as the standard library rounds it."""
    body = re.search(name + r" = \{(.*?)\};", text, re.S)[1]
    body = body.replace("std::numbers::inv_sqrt3", "0.5773502691896257")
    return [float(word) for word in body.replace(",", " ").split()]


def weighed(degree, low, high):
    """Bounds on the weight at the root that lies between two rationals."""
    reach = (degree - 1) ** 2 * (high - low)
    ends = [legendre(degree, x)[0] for x in (low, high)]
    small = min(abs(end) for end in ends) - reach
    large = max(abs(end) for end in ends) + reach
    square = degree * degree
    return (2 * (1 - high * high) / (square * large * large),
            2 * (1 - low * low) / (square * small * small))


def main():
    failed = False

    def report(name, held, told=""):
        nonlocal failed
        failed = failed or not held
        print(f"{name:52} {'held' if held else 'BROKEN':7} {told}", flush=True)

    text = (HERE.parent.parent / "source" / "legendre.hpp").read_text()
    nodes = table(text, "Nodes")
    weights = table(text, "Weights")
    report("the table holds every rule of 2 to 40 points",
           len(nodes) == len(weights) == PAIRS * (PAIRS + 1) // 2)
    for pairs in range(1, PAIRS + 1):
        degree = 2 * pairs
        first = pairs * (pairs - 1) // 2
        rule = nodes[first:first + pairs]
        nearest = all(a < b for a, b in zip([0.0] + rule, rule + [1.0]))
        weighted = True
        for node, weight in zip(rule, weights[first:first + pairs]):
            low, high = halfway(node)
            sign = legendre(degree, low)[1] > 0
            nearest = nearest and (legendre(degree, high)[1] > 0) != sign
            nearest = nearest and legendre(degree, low)[1] != 0 != legendre(degree, high)[1]
            while high - low > WIDTH:
                middle = (low + high) / 2
                if (legendre(degree, middle)[1] > 0) == sign:
                    low = middle
                else:
                    high = middle
            least, most = weighed(degree, low, high)
            below, above = halfway(weight)
            weighted = weighted and below < least <= most < above
        report(f"{degree} points: each node the nearest double", nearest)
        report(f"{degree} points: each weight the nearest double", weighted)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
