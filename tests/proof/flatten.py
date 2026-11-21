#!/usr/bin/env python3
"""Machine-checks with Z3 the slack and skew source/flatten.cpp sizes an arc's chords by (Leeways,
Wheeled, Chords), and works out in exact rationals how far its runs of turned points stray.

    python3 -I tests/proof/flatten.py [z3]

A rounding is lemma.py's: |r - z| <= u |z| and |r - z| <= u |r| for a sum, a difference, a product
clear of underflow or a quotient. Each lemma states the operations of the code as written, one fact
a rounding; they are not bound to program.py's graphs, so translation.py does not check them
against the compiler's IR. Each mutant shrinks a bound, and Z3 must refute it; each check's bound
is shown to be within four times what it bounds. The exit status is 1 when a lemma is not proved,
a check fails or a mutant is proved.

Wheeled gives an arc's point k, k = m 32 + j and j < 32, as c + rnd(v r) for v = Rolled(w_m, p_j):
w_m = Direction(start + sweep (k / n)), as At takes it, and p_j the step t = Direction(sweep / n)
turned j times by Rolled from (1, 0). Direction rounds each component of its unit vector to the
nearest double, so w_m and t are within u of theirs, in length too. By hand: Rolled(a, b) is the
complex product a b, a rotation of a by b's angle times |b|, plus each component's rounding, which
rolled bounds by (2u + u^2)(|a_x b_x| + |a_y b_y|); squared and summed over both components that
is at most 2 (2u + u^2)^2 |a|^2 |b|^2, as 4 |a_x a_y b_x b_y| <= |a|^2 |b|^2. So with a within E_a
of a unit vector A and b within E_b of B, |Rolled(a, b) - A B| <= E_a |b| + E_b + sqrt(2) (2u +
u^2) |a| |b|: powers works out E_j for p_j, from E_1 = u, and the point's direction within Q of the
exact one, and runs that Q |r| and the product's rounding, u |v| |r|, stay within 2^-45 |r|.

The exact place p_j stands for is w_m's angle plus j steps of t's angle, sweep / n rounded once.
The point is c + v r rounded at each coordinate, |c_i| + |r| (1 + 2^-40) bounding it: half the
step between the doubles there, at most, which Halfsteps takes; a sum that underflows is exact.

The angles: anchor bounds how far w_m's angle, start + rnd(sweep rnd(k / n)) rounded, lies from
start + sweep k / n: u (|start| + 3 |sweep|)(1 + 2u). A run's own steps are sweep / n rounded,
within u of an even step each, and a run has 31 at most, so its last point is within u |sweep| more
of its place; and the end, start + sweep rounded, is within u |start + sweep| of the arc's: skew
shows two anchors' errors and a run's drift within Leeways' skew, 2^-51 (|start| + 2 |sweep|), and
the end's gap, swept at the radius, is within Leeways' |r| 2^-52 (|start| + |sweep|). Degrees and
radians alike: Direction turns the unit exactly, and InRadians rounds twice at most, which the
factors of two leave room for.

Left by hand: a chord between points each within e of the ends of a chord of the exact circle
strays no more than e farther from the circle than that one, the point at each share of the one
within e of the other's; the exact chord of a step s strays 2 r sin^2(s / 4). Chords takes 1 -
2^-20 of the tolerance less the slack, for the roundings of the square root and the arc sine, and
the sweep 2^-40 longer. A bulge's points rest on the library's arc tangent and sine, as Along's
do: its slack, 2^-40 of its reach, is measured (tests/geoxy/flatten.cpp), not proved.
"""

import fractions
import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from lemma import Lemma, run  # pylint: disable=wrong-import-position

Fraction = fractions.Fraction
U = Fraction(1, 2 ** 53)
ROOT_TWO = Fraction(14143, 10000)
RELAY = 32


def rounds(lemma, result, exact):
    """result is exact rounded to nearest: a sum, a difference, a product clear of underflow or a
    quotient."""
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {exact})))")
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {result})))")


def rolled(share="(+ (* 2.0 u) (* u u))"):
    """One component of Rolled, rnd(rnd(a c) - rnd(b d)), within (2u + u^2)(|a c| + |b d|) of a c
    - b d: each product off by u of itself, the difference by u of the rounded products'."""
    lemma = Lemma("rolled", "a", "b", "c", "d", "p", "q", "x")
    rounds(lemma, "p", "(* a c)")
    rounds(lemma, "q", "(* b d)")
    rounds(lemma, "x", "(- p q)")
    lemma.goal(f"(<= (absr (- x (- (* a c) (* b d))))"
               f" (* {share} (+ (absr (* a c)) (absr (* b d)))))")
    return lemma


def anchor(times="3.0"):
    """w_m's angle: f = rnd(k / n), the share q in [0, 1]; g = rnd(sweep f); a = rnd(start + g),
    within u (|start| + 3 |sweep|)(1 + 2u) of start + sweep q."""
    lemma = Lemma("anchor", "q", "f", "g", "a", "start", "sweep")
    lemma.fact("(<= 0.0 q)")
    lemma.fact("(<= q 1.0)")
    rounds(lemma, "f", "q")
    rounds(lemma, "g", "(* sweep f)")
    rounds(lemma, "a", "(+ start g)")
    lemma.goal(f"(<= (absr (- a (+ start (* sweep q))))"
               f" (* u (+ (absr start) (* {times} (absr sweep))) (+ 1.0 (* 2.0 u))))")
    return lemma


def skew(scale="(/ 1.0 2251799813685248.0)"):
    """Two anchors' errors, A each, and a run's drift, u |sweep|, sum to no more than Leeways' skew,
    2^-51 (|start| + 2 |sweep|), taken as worked out: its product and sum each rounded."""
    lemma = Lemma("skew", "start", "sweep", "A", "s", "k")
    lemma.fact("(>= start 0.0)")
    lemma.fact("(>= sweep 0.0)")
    lemma.fact("(= A (* u (+ start (* 3.0 sweep)) (+ 1.0 (* 2.0 u))))")
    rounds(lemma, "s", "(+ start (* 2.0 sweep))")
    rounds(lemma, "k", f"(* {scale} s)")
    lemma.goal("(<= (+ (* 2.0 A) (* u sweep)) k)")
    return lemma


def powers():
    """E_j for p_j, from E_1 = u, each turn by t adding E_j u + u + sqrt(2) (2u + u^2)(1 +
    E_j)(1 + u); then the point's direction, Rolled(w_m, p_j), within Q_j = u (1 + E_j) + E_j +
    sqrt(2) (2u + u^2)(1 + u)(1 + E_j) of the exact, and the radius's share of the slack, Q_j + u
    (1 + Q_j), within 2^-45 for every j a run takes, and past 2^-47."""
    share = 2 * U + U * U
    error = U
    worst = Fraction(0)
    for _ in range(1, RELAY):
        point = U * (1 + error) + error + ROOT_TWO * share * (1 + U) * (1 + error)
        worst = max(worst, point + U * (1 + point))
        error = error * (1 + U) + U + ROOT_TWO * share * (1 + error) * (1 + U)
    held = Fraction(1, 2 ** 47) < worst <= Fraction(1, 2 ** 45)
    return held, f"radius's share {float(worst / U):.2f} u, within 2^-45 = 256 u"


MUTANTS = (("rolled at 1.5 u", lambda: rolled(share="(* 1.5 u)")),
           ("anchor at 2 |sweep|", lambda: anchor(times="2.0")),
           ("skew at 2^-52", lambda: skew(scale="(/ 1.0 4503599627370496.0)")))

LEMMAS = (rolled, anchor, skew)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, checks=(powers,), mutants=MUTANTS))
