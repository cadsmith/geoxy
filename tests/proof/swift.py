#!/usr/bin/env python3
"""Machine-checks the running bound source/swift/bernstein.cpp's Isolated keeps on the coefficients
of a line's equation on a curve as it halves: every coefficient within (first + n k 2 u largest)
Grow of the exact one after k halvings of a degree n curve.

    python3 -I tests/proof/swift.py [z3]

A rounding is lemma.py's: |r - z| <= u |z| for a sum, which underflow leaves exact; halving it
by 0.5 is exact but where it falls among the doubles below the normal ones, off by h = 2^-1075.

One halving is n rounds of de Casteljau's at 1/2, each coefficient rnd(a + b) 0.5 of two of the
round before: a convex combination of the exact ones, so no larger than largest + first, and
within E + u (largest + first + E) + h of it for inputs within E (halving round). Over m = n k
rounds, E_m <= (first + m u (largest + first))(1 + u)^m + m h, under the bound Isolated works out,
its three roundings taken down, for every m to n Deepest = 3 * 14, term by term; the m h is under
2^-40 of first, which carries sure.hpp's Tiny, 2^-1000 (running).

Left by hand: that Isolated's tests against the reach T beyond that bound tell f's sign on a
stretch, its hull holding f; Newton's step and the Taylor point are twin.py's taylor lemma's.

The lemma states the code's operations, one fact a rounding; it is not bound to program.py's
graphs. The mutant drops the round's own rounding, and Z3 must refute it. The exit status is 1
when a lemma is not proved or a mutant is.
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
HIGHEST = 3


def halving(own=True):
    """One round: a and b within E of the exact A and B, each at most S; s = rnd(a + b), w = s 0.5,
    exact or within h of it. Then w is within E + u (S + E) + h of (A + B) / 2."""
    lemma = Lemma("halving round", "a", "b", "A", "B", "S", "E", "s", "w")
    lemma.fact("(>= E 0.0)")
    lemma.fact("(<= (absr A) S)")
    lemma.fact("(<= (absr B) S)")
    lemma.fact("(<= (absr (- a A)) E)")
    lemma.fact("(<= (absr (- b B)) E)")
    if own:
        lemma.fact("(<= (absr (- s (+ a b))) (* u (absr (+ a b))))")
    else:
        lemma.fact("(= s (+ a b))")
    lemma.fact("(<= (absr (- w (* 0.5 s))) h)")
    lemma.goal("(<= (absr (- w (* 0.5 (+ A B)))) (+ E (* u (+ S E)) h))")
    return lemma


def running():
    """The bound Isolated works out, rnd(rnd(first + rnd(rnd(rnd(n k) 2 u) largest)) Grow), at
    least (first + 2 m u largest) Grow (1 - u)^3, m = n k exact, covers E_m term by term: first's
    share (1 + m u)(1 + u)^m and largest's m u (1 + u)^m; and m h is under 2^-40 2^-1000."""
    held = True
    for rounds in range(1, HIGHEST * DEEPEST + 1):
        given = GROW * (1 - U) ** 3
        held = held and (1 + rounds * U) * (1 + U) ** rounds <= given
        held = held and rounds * U * (1 + U) ** rounds <= 2 * rounds * U * given
        held = held and rounds * H <= (GROW * (1 - U) ** 3 - (1 + rounds * U) * (1 + U) ** rounds) \
            * Fraction(1, 2 ** 1000)
    return held, f"every m to {HIGHEST * DEEPEST} rounds"


MUTANTS = (("halving round without its rounding's u", lambda: _without()),)


def _without():
    """halving round, its goal's u (S + E) dropped though the sum rounds."""
    lemma = halving()
    lemma.goals = ["(<= (absr (- w (* 0.5 (+ A B)))) (+ E h))"]
    return lemma


LEMMAS = (halving,)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, checks=(running,), mutants=MUTANTS))
