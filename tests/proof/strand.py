#!/usr/bin/env python3
"""Machine-checks the bounds source/swift/strand.cpp's kernel of a line and a spline of any degree
adds to the line kernel's: the slope of the line's equation worked out in doubles from its rounded
Bernstein coefficients, and the running bound Isolated keeps on them, for every degree it takes.

    python3 -I tests/proof/strand.py [z3]

A rounding is lemma.py's: |r - z| <= u |z| for a sum, a difference, or a product clear of
underflow.

The slope: each coefficient c_k of degree n, within E of the exact C_k and no larger than S, its
difference rnd(c_k+1 - c_k) times n rounded (slope step): within 2 n E + 4 n u S (1 + u) of n
(C_k+1 - C_k), and no larger than 2 n S (1 + u)^2. De Casteljau's n - 1 rounds on them, each
twin.py's casteljau round, add three steps of their largest a round: f' at a place within 2 n E +
(4 n + 6 n (n - 1)) u S (1 + u)^2 + 6 n (n - 1) u E of the exact one, which the kernel's (2 n E +
(6 n^2 + 4 n) u S) Grow covers for every degree to sixteen, the u E terms within Grow's share of 2 n
E (slack, checked in exact rationals): fifteen for a span of strand.cpp's, sixteen for circuit.hpp's
equation on a span of degree eight.

The running bound: swift.py's for a cubic, here for every m = n k to sixteen times Deepest
(running).

Left by hand, as for the line kernel: that a spline's point times its weight, less the line's
point times the weight, is the spline function whose coefficient at a pole is its weight times the
line's equation there; that a span's Bernstein coefficients are that function's blossoms at the
span's ends, as section.py bounds them; and that a rational span's derivative is no larger than 2 n
w_max over w_min times the diagonal of its poles' hull.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from expand import Fraction, U  # pylint: disable=wrong-import-position
# pylint: disable-next=wrong-import-position
from lemma import Lemma, run  # pylint: disable=wrong-import-position

GROW = 1 + Fraction(1, 2 ** 40)
H = Fraction(1, 2 ** 1075)
DEEPEST = 14
STRANDS = 16


def step(both=True):
    """One slope: c and d within E of the exact C and D, each no larger than S; g = rnd(d - c), s =
    rnd(n g) for a whole n of at least 1. Then s is within 2 n E + 4 n u S (1 + u) of n (D - C)."""
    lemma = Lemma("slope step", "c", "d", "C", "D", "S", "E", "n", "g", "s")
    lemma.fact("(>= E 0.0)")
    lemma.fact("(>= n 1.0)")
    lemma.fact("(<= (absr c) S)")
    lemma.fact("(<= (absr d) S)")
    lemma.fact("(<= (absr (- c C)) E)")
    lemma.fact("(<= (absr (- d D)) E)")
    lemma.fact("(<= (absr (- g (- d c))) (* u (absr (- d c))))")
    if both:
        lemma.fact("(<= (absr (- s (* n g))) (* u (absr (* n g))))")
    else:
        lemma.fact("(= s (* n g))")
    lemma.goal("(<= (absr (- s (* n (- D C)))) (+ (* 2.0 n E) (* 4.0 n u S (+ 1.0 u))))")
    return lemma


def slack():
    """The kernel's slope bound covers the slope steps and n - 1 casteljau rounds, each adding three
    steps of the slopes' largest, 2 n (S + E) (1 + u)^2, for every degree from two to fifteen, the
    bound's own roundings, three of them, taken down."""
    held = True
    for order in range(2, STRANDS + 1):
        n = Fraction(order)
        largest = 2 * n * (1 + U) ** 2
        # The error per unit of S and of E, as the rounds carry them.
        on_size = 4 * n * U * (1 + U) + 3 * (n - 1) * U * largest
        on_error = 2 * n + 3 * (n - 1) * U * 2 * n * (1 + U) ** 2
        given = GROW * (1 - U) ** 3
        held = held and on_size <= (6 * n * n + 4 * n) * U * given
        held = held and on_error <= 2 * n * given
    return held, f"every degree to {STRANDS}"


def running():
    """swift.py's running bound for every m = n k to STRANDS times Deepest rounds."""
    held = True
    for rounds in range(1, STRANDS * DEEPEST + 1):
        given = GROW * (1 - U) ** 3
        held = held and (1 + rounds * U) * (1 + U) ** rounds <= given
        held = held and rounds * U * (1 + U) ** rounds <= 2 * rounds * U * given
        held = held and rounds * H <= (GROW * (1 - U) ** 3 - (1 + rounds * U) * (1 + U) ** rounds) \
            * Fraction(1, 2 ** 1000)
    return held, f"every m to {STRANDS * DEEPEST} rounds"


MUTANTS = (("slope step without the scaling's rounding", lambda: _unscaled()),)


def _unscaled():
    """slope step, its goal's scaling term halved though the product rounds."""
    lemma = step()
    lemma.goals = ["(<= (absr (- s (* n (- D C)))) (+ (* 2.0 n E) (* 2.0 n u S (+ 1.0 u))))"]
    return lemma


LEMMAS = (step,)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, checks=(slack, running), mutants=MUTANTS))
