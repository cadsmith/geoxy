#!/usr/bin/env python3
"""Machine-checks with Z3 the bounds source/swift/circuit.hpp's kernel for an arc or an ellipse met
with a Bézier curve or another conic rests on: the reach the conic's equation on the other curve
must clear for that curve to lie farther than the reach from the conic, the error of a coefficient
of that equation worked out in doubles, the sine at a crossing, the margin a direction on the unit
circle keeps from a sweep's end, and the quarter of the unit circle as a rational quadratic piece;
and that a coefficient's 16 u holds for a span of every degree to Broadest, eight.

    python3 -I tests/proof/ring.py [z3]

A rounding is lemma.py's: |r - z| <= u |z| for a sum, a difference, a product or a quotient clear
of underflow. Each lemma states the operations of circuit.hpp as written, one fact a rounding; they
are not bound to program.py's graphs, so translation.py does not check them against the
compiler's IR. Each mutant drops a term the code keeps, and Z3 must refute it: so the proofs need
every part of the bounds they are given. The exit status is 1 when a lemma is not proved or a
mutant is.

Left to the comments of circuit.hpp, worked by hand: that a point at distance s from an ellipse is
within s / kappa of it in the unit circle's plane, kappa = 1 / sigma_max of its map, and a point of
the ellipse within s of a crossing is within s / sigma_min of it on the unit circle; the mean value
theorem's bound on the one double-double Newton step, as the kernel of a line and a Bézier curve
has it; that the isolation's stretches hold every place where the other curve comes within the
reach; and the reasoning that the thorough way's tolerance decides nothing about a crossing the
kernel is sure of.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from expand import Fraction, U  # pylint: disable=wrong-import-position
# pylint: disable-next=wrong-import-position
from lemma import Lemma, rational, run  # pylint: disable=wrong-import-position

GROW = 1.0 + 2.0 ** -40
BROADEST = 8


def rounds(lemma, result, exact):
    """result is exact rounded to nearest: a sum, a difference, a product or a quotient clear of
    underflow."""
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {exact})))")


def reach(heavy=True):
    """The reach: a point p of the other curve at its place t has H = W^2 (|v|^2 - rho^2), v = A
    (p - o) and W its weight, no more than W_max; its distance from the carrier is at least kappa
    | |v| - rho |, and |v| no more than U. Where |H| kappa > reach W_max^2 (U + rho), that distance
    is above the reach."""
    lemma = Lemma("reach", "v", "rho", "big", "kappa", "w", "most", "level", "far", "gap")
    lemma.fact("(>= v 0.0)")
    lemma.fact("(> rho 0.0)")
    lemma.fact("(<= v big)")
    lemma.fact("(> kappa 0.0)")
    lemma.fact("(> w 0.0)")
    lemma.fact("(<= w most)")
    lemma.fact("(>= far 0.0)")
    lemma.fact("(= level (* w w (- (* v v) (* rho rho))))")
    lemma.fact("(>= gap (* kappa (absr (- v rho))))")
    weight = "(* most most)" if heavy else "1.0"
    lemma.fact(f"(> (* (absr level) kappa) (* far {weight} (+ big rho)))")
    lemma.goal("(> gap far)")
    return lemma


def product(both=True):
    """A product of two terms each within its error of the exact one: |a b - A B| <= |a| e_B +
    |b| e_A + e_A e_B, the error a coefficient's term carries before its own roundings."""
    lemma = Lemma("product", "a", "b", "A", "B", "ea", "eb")
    lemma.fact("(>= ea 0.0)")
    lemma.fact("(>= eb 0.0)")
    lemma.fact("(<= (absr (- a A)) ea)")
    lemma.fact("(<= (absr (- b B)) eb)")
    cross = "(* ea eb)" if both else "0.0"
    lemma.goal(f"(<= (absr (- (* a b) (* A B))) (+ (* (absr a) eb) (* (absr b) ea) {cross}))")
    return lemma


def term(share=3.0):
    """A coefficient's term, worked out as the code does: p = rnd(x y), q = rnd(z w), d = rnd(p +
    q), t = rnd(m d) for a multiplier m >= 0. It is within 3 u Grow of its size m (|x y| + |z w|)
    of the exact term: its two products, its sum and the multiplier's product."""
    lemma = Lemma("term", "x", "y", "z", "w", "m", "p", "q", "d", "t")
    lemma.fact("(>= m 0.0)")
    rounds(lemma, "p", "(* x y)")
    rounds(lemma, "q", "(* z w)")
    rounds(lemma, "d", "(+ p q)")
    rounds(lemma, "t", "(* m d)")
    lemma.goal(f"(<= (absr (- t (* m (+ (* x y) (* z w))))) "
               f"(* {share} u {rational(GROW)} m (+ (absr (* x y)) (absr (* z w)))))")
    return lemma


def summed(share=3.0):
    """A coefficient's running sum of four terms, as a cubic piece's middle one has: each t_i
    within e_i of its exact T_i, no larger than M_i; s1 = rnd(t0 + t1), s2 = rnd(s1 + t2), s3 =
    rnd(s2 + t3). It is within the e_i and 3 u Grow of the M_i and e_i of the exact sum, so that
    with term's 3 u each the sum is within 6 u Grow of the terms' sizes, and the 16 u circuit.hpp
    charges covers it and finished's 4 u."""
    lemma = Lemma("summed", "t0", "t1", "t2", "t3", "T0", "T1", "T2", "T3",
                  "e0", "e1", "e2", "e3", "M0", "M1", "M2", "M3", "s1", "s2", "s3")
    for index in range(4):
        lemma.fact(f"(>= e{index} 0.0)")
        lemma.fact(f"(<= (absr (- t{index} T{index})) e{index})")
        lemma.fact(f"(<= (absr T{index}) M{index})")
    rounds(lemma, "s1", "(+ t0 t1)")
    rounds(lemma, "s2", "(+ s1 t2)")
    rounds(lemma, "s3", "(+ s2 t3)")
    lemma.goal(f"(<= (absr (- s3 (+ T0 T1 T2 T3))) "
               f"(+ e0 e1 e2 e3 (* {share} u {rational(GROW)} (+ M0 M1 M2 M3 e0 e1 e2 e3))))")
    return lemma


def running(grown=True):
    """One more term of a coefficient's running sum, for a span of any degree: the sum so far s
    within D of its exact S, which is no larger than P; the next term t within e of its exact T, no
    larger than m; r = rnd(s + t). Then r is within D + e + u (P + D + m + e) of S + T: the errors
    carried, and the sum's rounding of a size no more than the exact parts' and their errors."""
    lemma = Lemma("running", "s", "S", "t", "T", "D", "e", "P", "m", "r")
    for name in ("D", "e", "m"):
        lemma.fact(f"(>= {name} 0.0)")
    lemma.fact("(<= (absr (- s S)) D)")
    lemma.fact("(<= (absr S) P)")
    lemma.fact("(<= (absr (- t T)) e)")
    lemma.fact("(<= (absr T) m)")
    rounds(lemma, "r", "(+ s t)")
    rounding = "(* u (+ P D m e))" if grown else "0.0"
    lemma.goal(f"(<= (absr (- r (+ S T))) (+ D e {rounding}))")
    return lemma


def wide():
    """A coefficient of the equation on a span of any degree n to Broadest sums n + 1 terms, each
    within 3 u Grow of its size m_i (term), their sizes M in all. By running, n times, the sum is
    within D of the exact one: the terms' errors, 3 u Grow M in all, and each rounding's u times
    the sizes so far and the new term's, no more than M, and the errors, no more than D and 3 u
    Grow M; so D (1 - n u) is no more than 3 u Grow M + n u M (1 + 3 u Grow). Finished adds 4 u
    of the sum's size, no more than M + D, and of the weights' term. All of it within the 16 u
    circuit.hpp charges, for every n to Broadest; and the weights' sum, each product of a binomial
    weight and two weights rounded twice and n + 1 of them summed, within (1 + u)^(n + 2) - 1 of
    its size, which with finished's 4 u is within it too. Checked in exact rationals."""
    grow = Fraction(GROW)
    held = True
    for order in range(1, BROADEST + 1):
        term = 3 * U * grow
        error = (term + order * U * (1 + term)) / (1 - order * U)
        held = held and error + 4 * U * (1 + error) <= 16 * U
        weights = (1 + U) ** (order + 2) - 1
        held = held and weights + 4 * U * (1 + weights) <= 16 * U
    return held, f"every degree to {BROADEST}"


def finished(share=4.0):
    """A coefficient finished: v = rnd(rnd(s - rnd(l w)) / c), the sum s within e_s of its exact
    S, rho^2's double l within e_l of L, the weights' product w >= 0 exact and the binomial c >= 1.
    It is within (e_s + e_l w + 4 u (|s| + |l| w)) / c of (S - L w) / c: the product's, the
    difference's and the quotient's roundings."""
    lemma = Lemma("finished", "s", "S", "l", "L", "w", "c", "es", "el", "lw", "gap", "v")
    for name in ("es", "el", "w"):
        lemma.fact(f"(>= {name} 0.0)")
    lemma.fact("(>= c 1.0)")
    lemma.fact("(<= (absr (- s S)) es)")
    lemma.fact("(<= (absr (- l L)) el)")
    rounds(lemma, "lw", "(* l w)")
    rounds(lemma, "gap", "(- s lw)")
    rounds(lemma, "v", "(/ gap c)")
    lemma.goal(f"(<= (* c (absr (- v (/ (- S (* L w)) c)))) "
               f"(+ es (* el w) (* {share} u (+ (absr s) (* (absr l) w)))))")
    return lemma


def sine(heavy=True):
    """The sine at a crossing: H' = W^2 (grad Q . p') where Q is nothing, |grad Q| = 2 |A^T v| no
    more than 2 |A| rho, |p'| no more than the speed, W no more than W_max; where |H'| >= g, the
    sine |grad Q . p'| / (|grad Q| |p'|) is at least g / (W_max^2 2 |A| rho speed)."""
    lemma = Lemma("sine", "dot", "w", "most", "slope", "g", "normal", "norm", "rho", "speed", "way")
    for name in ("normal", "norm", "rho", "speed", "way"):
        lemma.fact(f"(> {name} 0.0)")
    lemma.fact("(> w 0.0)")
    lemma.fact("(<= w most)")
    lemma.fact("(= slope (* w w dot))")
    lemma.fact("(>= (absr slope) g)")
    lemma.fact("(>= g 0.0)")
    lemma.fact("(<= normal (* 2.0 norm rho))")
    lemma.fact("(<= way speed)")
    weight = "(* most most)" if heavy else "1.0"
    lemma.goal(f"(>= (* (absr dot) {weight} 2.0 norm rho speed) (* g normal way))")
    return lemma


def margin(sized=True):
    """A direction's side of a sweep's end s, a unit vector: where |s x e| > eta and f lies
    within eta of e, s x f has e's sign, as s x (f - e) is no more than |s| |f - e| by Lagrange's
    identity, twin.py's lagrange."""
    lemma = Lemma("margin", "s1", "s2", "e1", "e2", "f1", "f2", "eta")
    lemma.fact("(>= eta 0.0)")
    lemma.fact("(<= (absr (- (* s1 (- f2 e2)) (* s2 (- f1 e1)))) eta)")
    bound = "eta" if sized else "(* 0.5 eta)"
    lemma.fact(f"(> (absr (- (* s1 e2) (* s2 e1))) {bound})")
    lemma.goal("(=> (> (- (* s1 e2) (* s2 e1)) 0.0) (> (- (* s1 f2) (* s2 f1)) 0.0))")
    lemma.goal("(=> (< (- (* s1 e2) (* s2 e1)) 0.0) (< (- (* s1 f2) (* s2 f1)) 0.0))")
    return lemma


def quarter(weighted=True):
    """The quarter of the unit circle from angle 0 to pi / 2: control points (1, 0), (1, 1) and
    (0, 1) of weights 1, 1 and 2 give the numerator (1 - t^2, 2 t) and the weight 1 + t^2, whose
    point lies on the unit circle."""
    lemma = Lemma("quarter", "t", "nx", "ny", "w")
    last = "2.0" if weighted else "1.0"
    lemma.fact("(= nx (+ (* (- 1.0 t) (- 1.0 t)) (* 2.0 t (- 1.0 t))))")
    lemma.fact(f"(= ny (+ (* 2.0 t (- 1.0 t)) (* t t {last})))")
    lemma.fact(f"(= w (+ (* (- 1.0 t) (- 1.0 t)) (* 2.0 t (- 1.0 t)) (* t t {last})))")
    lemma.goal("(= nx (- 1.0 (* t t)))")
    lemma.goal("(= ny (* 2.0 t))")
    lemma.goal("(= w (+ 1.0 (* t t)))")
    lemma.goal("(= (+ (* nx nx) (* ny ny)) (* w w))")
    return lemma


MUTANTS = (("reach without the weight", lambda: reach(heavy=False)),
           ("product without ea eb", lambda: product(both=False)),
           ("term at 2 u", lambda: term(share=2.0)),
           ("summed at 2 u", lambda: summed(share=2.0)),
           ("finished at 1 u", lambda: finished(share=1.0)),
           ("running without the sum's rounding", lambda: running(grown=False)),
           ("sine without the weight", lambda: sine(heavy=False)),
           ("margin at half eta", lambda: margin(sized=False)),
           ("quarter of weights 1", lambda: quarter(weighted=False)))

LEMMAS = (reach, product, term, summed, running, finished, sine, margin, quarter)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, checks=(wide,), mutants=MUTANTS))
