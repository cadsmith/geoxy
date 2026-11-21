#!/usr/bin/env python3
"""Machine-checks with Z3 the bounds source/swift/twin.cpp's kernel for two Bézier curves rests on:
a round of de Casteljau's algorithm in doubles, the band Krawczyk's operator leaves about a root,
the point worked out by Taylor from a grid place, and the least sine between two disks of
derivatives.

    python3 -I tests/proof/twin.py [z3]

A rounding is lemma.py's: |r - z| <= u |z| for a sum, a difference, or a product clear of
underflow. Each lemma states the operations of twin.cpp as written, one fact a rounding; they are
not bound to program.py's graphs, so translation.py does not check them against the compiler's
IR. Each mutant drops a term the code keeps, and Z3 must refute it: so the proofs need every part
of the bounds they are given. The exit status is 1 when a lemma is not proved or a mutant is.

Left to the comments of twin.cpp, worked by hand: that Krawczyk's operator inside its box proves
one root there and only one (Krawczyk 1969, Moore 1977), with the gap's and the derivatives'
ranges each covering what rounding leaves; the halving's and the fat line's slack; the disk a
derivative keeps to on a stretch, its Bernstein coefficients on every part of it within the same
disk; and the reasoning that the thorough way's tolerance decides nothing about a crossing the
kernel is sure of.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from lemma import Lemma, rational, run  # pylint: disable=wrong-import-position

GROW = 1.0 + 2.0 ** -40


def rounds(lemma, result, exact):
    """result is exact rounded to nearest: a sum, a difference, or a product clear of
    underflow."""
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {exact})))")


def casteljau(share=3.0, propagated=True):
    """Valued, one round: rest = rnd(1 - t), p = rnd(c rest), q = rnd(d t), r = rnd(p + q), for
    coefficients c and d no larger than S and t from 0 to 1, each within E of the exact
    coefficients C and D. Then r is within E + 3 u S Grow of C (1 - t) + D t: the roundings add
    u S ((1 - t) + (1 - t) + t + 1) and their products' u^2 terms, which Grow covers; the exact
    weights make the round a convex combination, which carries the coefficients' error no
    further."""
    lemma = Lemma("casteljau round", "c", "d", "C", "D", "t", "S", "E", "rest", "p", "q", "r")
    lemma.fact("(<= 0.0 t)")
    lemma.fact("(<= t 1.0)")
    lemma.fact("(<= (absr c) S)")
    lemma.fact("(<= (absr d) S)")
    lemma.fact("(<= (absr (- c C)) E)")
    lemma.fact("(<= (absr (- d D)) E)")
    rounds(lemma, "rest", "(- 1.0 t)")
    rounds(lemma, "p", "(* c rest)")
    rounds(lemma, "q", "(* d t)")
    rounds(lemma, "r", "(+ p q)")
    carried = "E" if propagated else "0.0"
    lemma.goal(f"(<= (absr (- r (+ (* C (- 1.0 t)) (* D t)))) "
               f"(+ {carried} (* {share} u S {rational(GROW)})))")
    return lemma


def band(spread=True):
    """Krawczyk's band: for a place x of the box where |F(x)| <= l, x - x* = C (x - x*) + Y F(x)
    with |C| <= Ĉ entrywise, so d = |x - x*| has d <= Ĉ d + g, g = l times each row of Y's length.
    Where Ĉ r < r for the box's radii r > 0, as the operator inside the box makes it, the
    determinant of I - Ĉ is above zero and d <= (I - Ĉ)^-1 g, entry by entry."""
    lemma = Lemma("band", "c11", "c12", "c21", "c22", "r1", "r2", "g1", "g2", "d1", "d2")
    for name in ("c11", "c12", "c21", "c22", "g1", "g2", "d1", "d2"):
        lemma.fact(f"(>= {name} 0.0)")
    lemma.fact("(> r1 0.0)")
    lemma.fact("(> r2 0.0)")
    lemma.fact("(< (+ (* c11 r1) (* c12 r2)) r1)")
    lemma.fact("(< (+ (* c21 r1) (* c22 r2)) r2)")
    lemma.fact("(<= d1 (+ (* c11 d1) (* c12 d2) g1))")
    lemma.fact("(<= d2 (+ (* c21 d1) (* c22 d2) g2))")
    rest = "(- (* (- 1.0 c11) (- 1.0 c22)) (* c12 c21))"
    across = "(* c12 g2)" if spread else "0.0"
    lemma.goal(f"(> {rest} 0.0)")
    lemma.goal(f"(<= (* d1 {rest}) (+ (* (- 1.0 c22) g1) {across}))")
    lemma.goal(f"(<= (* d2 {rest}) (+ (* c21 g1) (* (- 1.0 c11) g2)))")
    return lemma


def taylor(slipped=True):
    """The point by Taylor from the grid place s0: A(s*) = A(s0) + A'(s0) D + R, |R| <= B D^2 / 2
    for D = s* - s0 and B the most |A''| may be, a coordinate; D = m + e, m the step and |e| no
    more than its reach; A'(s0) worked out as a, within the slip of the true one, which is no more
    than the speed. Then A(s*) is within speed reach + slip |m| + B (|m| + reach)^2 / 2 of
    A(s0) + a m."""
    lemma = Lemma("taylor", "at", "a", "m", "e", "R", "B", "speed", "slip", "reach", "start")
    lemma.fact("(>= B 0.0)")
    lemma.fact("(>= slip 0.0)")
    lemma.fact("(<= (absr at) speed)")
    lemma.fact("(<= (absr (- at a)) slip)")
    lemma.fact("(<= (absr e) reach)")
    lemma.fact("(<= (absr R) (* B (+ m e) (+ m e) 0.5))")
    slip = "(* slip (absr m))" if slipped else "0.0"
    truth = "(+ start (* at (+ m e)) R)"
    lemma.goal(f"(<= (absr (- {truth} (+ start (* a m)))) "
               f"(+ (* speed reach) {slip} (* B (+ (absr m) reach) (+ (absr m) reach) 0.5)))")
    return lemma


def lagrange():
    """|x cross y| <= |x| |y|: Lagrange's identity, |x|^2 |y|^2 = (x . y)^2 + (x cross y)^2."""
    lemma = Lemma("lagrange", "x1", "x2", "y1", "y2", "nx", "ny")
    lemma.fact("(>= nx 0.0)")
    lemma.fact("(>= ny 0.0)")
    lemma.fact("(= (* nx nx) (+ (* x1 x1) (* x2 x2)))")
    lemma.fact("(= (* ny ny) (+ (* y1 y1) (* y2 y2)))")
    lemma.goal("(<= (absr (- (* x1 y2) (* x2 y1))) (* nx ny))")
    return lemma


def sine(both=True):
    """LeastSine: for a within ra of a0 and b within rb of b0, |a x b| >= |a0 x b0| - (|a0| rb
    + |b0| ra + ra rb), each cross product of a vector and an offset no more than their lengths'
    product by lagrange; and |a| |b| <= (|a0| + ra) (|b0| + rb). So the sine between any such a and
    b is at least the first over the second, where the first is above zero, every pair turning the
    same way."""
    lemma = Lemma("sine", "w0", "wb", "wa", "wab", "na", "nb", "ra", "rb", "size")
    lemma.fact("(>= ra 0.0)")
    lemma.fact("(>= rb 0.0)")
    lemma.fact("(>= na 0.0)")
    lemma.fact("(>= nb 0.0)")
    lemma.fact("(<= (absr wb) (* na rb))")
    lemma.fact("(<= (absr wa) (* nb ra))")
    lemma.fact("(<= (absr wab) (* ra rb))")
    lemma.fact("(<= size (* (+ na ra) (+ nb rb)))")
    lemma.fact("(>= size 0.0)")
    least = "(- (absr w0) (+ (* na rb) (* nb ra) (* ra rb)))" if both else \
        "(- (absr w0) (+ (* na rb) (* nb ra)))"
    lemma.goal(f"(>= (absr (+ w0 wb wa wab)) {least})")
    lemma.goal(f"(<= size (* (+ na ra) (+ nb rb)))")
    return lemma


MUTANTS = (("casteljau round at 2.5 u", lambda: casteljau(share=2.5)),
           ("casteljau round without the carried error", lambda: casteljau(propagated=False)),
           ("band without the other row's share", lambda: band(spread=False)),
           ("taylor without the slip", lambda: taylor(slipped=False)),
           ("sine without ra rb", lambda: sine(both=False)))

LEMMAS = (casteljau, band, taylor, lagrange, sine)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, mutants=MUTANTS))
