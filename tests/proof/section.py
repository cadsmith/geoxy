#!/usr/bin/env python3
"""Machine-checks with Z3 the bound source/swift/section.cpp's spans of a polynomial spline rest
on, and nearest.cpp's branch and bound with them: a round of de Boor's algorithm in doubles, a mix
of two cells by a fraction of knot differences worked out in doubles, adds at most ten steps u of
a double at the largest size S of the cells and the exact ones to what the cells were off by,
which it carries as it was. Over n rounds, a span's Bézier control point, or a point At works out
within the span, is within 10 n u S of the true one, S no more than the largest coordinate M of
the span's poles grown by that: section.cpp's slack of 11 n u M.

    python3 -I tests/proof/section.py [z3]

A rounding is lemma.py's: |r - z| <= u |z| for a sum, a difference, or a product clear of
underflow; and for a quotient, which IEEE 754 rounds to nearest as it does a product, clear of
underflow as Sectioned checks a fraction is. Each lemma states the operations of section.cpp as
written, one fact a rounding; they are not bound to program.py's graphs, so translation.py does
not check them against the compiler's IR. Each mutant shrinks a bound, and Z3 must refute it.

Left to the comments of section.cpp, worked by hand: that every exact fraction of a round lies
from 0 to 1, the value a blossom is taken at lying between the two knots the round divides by, so
that every exact cell is a convex combination of the span's poles; that the fraction worked out
lies from 0 to 1 too, rounding being monotone; and that n rounds, each adding at most ten steps
of its cells' largest size, leave a control point within 10 n u S, every cell no larger than M
grown by that.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from lemma import Lemma, rational, run  # pylint: disable=wrong-import-position
# pylint: disable-next=wrong-import-position
from twin import casteljau  # pylint: disable=wrong-import-position

GROW = 1.0 + 2.0 ** -40
FRACTION = 3.0 + 2.0 ** -20


def fraction(share=FRACTION):
    """A round's fraction: ahead = rnd(x - s), width = rnd(e - s), q = rnd(ahead / width), for
    knots s <= x <= e, s < e, the exact fraction f = (x - s) / (e - s), from 0 to 1. Each
    rounding is a relative error of at most u: ahead = (x - s)(1 + a), width = (e - s)(1 + b) and
    q = (ahead / width)(1 + c), so q (1 + b) = f g for g = (1 + a)(1 + c), within 2 u + u^2 of 1;
    and q is within (3 + 2^-20) u f of f."""
    lemma = Lemma("blossom fraction", "f", "q", "a", "b", "c", "g")
    lemma.fact("(<= 0.0 f)")
    lemma.fact("(<= f 1.0)")
    for name in ("a", "b", "c"):
        lemma.fact(f"(<= (absr {name}) u)")
    lemma.fact("(= g (* (+ 1.0 a) (+ 1.0 c)))")
    lemma.fact("(<= (absr (- g 1.0)) (+ (* 2.0 u) (* u u)))")
    lemma.fact("(= (* q (+ 1.0 b)) (* f g))")
    lemma.goal(f"(<= (absr (- q f)) (* {rational(share)} u f))")
    return lemma


def product():
    """(1 + a)(1 + c) is within 2 u + u^2 of 1 for |a|, |c| <= u: fraction's g."""
    lemma = Lemma("blossom product", "a", "c", "g")
    lemma.fact("(<= (absr a) u)")
    lemma.fact("(<= (absr c) u)")
    lemma.fact("(= g (* (+ 1.0 a) (+ 1.0 c)))")
    lemma.goal("(<= (absr (- g 1.0)) (+ (* 2.0 u) (* u u)))")
    return lemma


def mix(share=10.0):
    """A round's mix: rest = rnd(1 - q), r = rnd(rnd(c rest) + rnd(d q)), for the fraction q
    worked out, from 0 to 1 as rounding is monotone and x - s is no more than e - s: twin.py's
    casteljau round at q, which puts r within E + 3 u S Grow of C (1 - q) + D q for cells c and d
    within E of the exact C and D, all four no larger than S. Moving q to the exact fraction f
    moves that by (D - C)(q - f), at most 2 S (3 + 2^-20) u: r is within E + 10 u S Grow of the
    exact mix C (1 - f) + D f, the cells' error carried as it was."""
    lemma = Lemma("blossom mix", "C", "D", "f", "q", "S", "E", "r")
    lemma.fact("(<= 0.0 f)")
    lemma.fact("(<= f 1.0)")
    lemma.fact("(>= E 0.0)")
    lemma.fact(f"(<= (absr (- q f)) (* {rational(FRACTION)} u f))")
    lemma.fact("(<= (absr C) S)")
    lemma.fact("(<= (absr D) S)")
    lemma.fact(f"(<= (absr (- r (+ (* C (- 1.0 q)) (* D q)))) (+ E (* 3.0 u S {rational(GROW)})))")
    lemma.goal(f"(<= (absr (- r (+ (* C (- 1.0 f)) (* D f)))) "
               f"(+ E (* {rational(share)} u S {rational(GROW)})))")
    return lemma


MUTANTS = (("fraction within 2.9 u", lambda: fraction(share=2.9)),
           ("mix within 8.9 u", lambda: mix(share=8.9)))

LEMMAS = (product, fraction, casteljau, mix)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, mutants=MUTANTS))
