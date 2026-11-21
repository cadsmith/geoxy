#!/usr/bin/env python3
"""Machine-checks both predicates' first stages with Z3, through underflow, overflow and NaN.

    python3 -I tests/proof/first.py [z3]

The first stages are inline in source/exact/orient.hpp and source/exact/circle.hpp: a determinant
rounded, a bound, one comparison. Their lemmas prove that a comparison that passes gives the sign
of the exact determinant, for every pair of finite doubles: the orientation test's bound,
OrientBoundA (|left + right| + Normal), and the circle test's, CircleBoundA permanent +
CircleBoundUnderflow (lifts + 1), each with its products allowed to underflow (lemma.py). Then,
in IEEE 754 arithmetic itself, that an infinity or a NaN anywhere makes the comparison fail.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
import lemma  # pylint: disable=wrong-import-position
import program  # pylint: disable=wrong-import-position
from lemma import (CIRCLE_BOUND_A, CIRCLE_BOUND_UNDERFLOW, NORMAL, ORIENT_BOUND_A, Floating,
                   Lemma, rational, run)

ALPHA = "(+ (* 4.0 u) (* 10.0 u u))"


def tails():
    """A product of two exact differences against the product of their roundings:
    |X Y - x y| <= (2u + u^2) |x y|, where |X - x| <= u |x| and |Y - y| <= u |y|."""
    lemma = Lemma("tails", "x", "y", "X", "Y")
    lemma.fact("(<= (absr (- X x)) (* u (absr x)))")
    lemma.fact("(<= (absr (- Y y)) (* u (absr y)))")
    lemma.goal("(<= (absr (- (* X Y) (* x y))) (* (+ (* 2.0 u) (* u u)) (absr (* x y))))")
    return lemma


def orient_filter():
    """The orientation test's first stage: left = fl((ax - cx)(by - cy)), right likewise, the
    determinant fl(left - right) and the bound fl(OrientBoundA fl(|fl(left + right)| + Normal)).
    P1 and P2 are the exact products of the exact differences, tails' facts binding them."""
    lemma = Lemma("orient filter", "P1", "d1", "left", "P2", "d2", "right", "det", "sum", "w",
                  "bound")
    for exact, product in (("P1", "d1"), ("P2", "d2")):
        lemma.fact(f"(<= (absr (- {exact} {product})) (* (+ (* 2.0 u) (* u u)) (absr {product})))")
    lemma.multiplied("left", "d1")
    lemma.multiplied("right", "d2")
    lemma.rounded("det", "(- left right)")
    lemma.rounded("sum", "(+ left right)")
    lemma.rounded("w", f"(+ (absr sum) {rational(NORMAL)})")
    lemma.multiplied("bound", f"(* {rational(ORIENT_BOUND_A)} w)")
    lemma.exceeds("det", "bound")
    lemma.same_sign("(- P1 P2)", "det")
    v = program.orient_first()
    lemma.code({"d1": v["left"], "left": v["left"], "d2": v["right"], "right": v["right"],
                "det": v["det"], "sum": v["sum"], "w": v["w"], "bound": v["bound"]})
    return lemma


def amgm():
    """|b c| <= (b^2 + c^2) / 2: a minor's product is bounded by the two corners' lifts."""
    lemma = Lemma("amgm", "b", "c")
    lemma.goal("(<= (absr (* b c)) (/ (+ (* b b) (* c c)) 2.0))")
    return lemma


def circle_lift():
    """A lift of the first stage: fl(fl(x x) + fl(y y)), against the exact differences' X^2 + Y^2,
    tails' facts binding the squares."""
    lemma = Lemma("circle lift", "x2", "y2", "X2", "Y2", "sxx", "syy", "lift")
    lemma.fact("(>= x2 0.0)")
    lemma.fact("(>= y2 0.0)")
    lemma.fact("(<= (absr (- X2 x2)) (* (+ (* 2.0 u) (* u u)) x2))")
    lemma.fact("(<= (absr (- Y2 y2)) (* (+ (* 2.0 u) (* u u)) y2))")
    lemma.multiplied("sxx", "x2")
    lemma.multiplied("syy", "y2")
    lemma.rounded("lift", "(+ sxx syy)")
    lemma.goal(f"(<= (absr (- (+ X2 Y2) lift)) (+ (* {ALPHA} lift) (* 3.0 h)))")
    lemma.goal("(<= (+ x2 y2) (+ (* (+ 1.0 (* 3.0 u)) lift) (* 3.0 h)))")
    lemma.goal("(>= lift 0.0)")
    v = program.circle_first()
    lemma.code({"x2": v["axx"], "y2": v["ayy"], "sxx": v["axx"], "syy": v["ayy"],
                "lift": v["alift"]})
    return lemma


def circle_minor():
    """A minor of the first stage, fl(fl(bx cy) - fl(cx by)), and its permanent part,
    fl(|fl(bx cy)| + |fl(cx by)|), against the exact differences' minor."""
    lemma = Lemma("circle minor", "bc", "cb", "BC", "CB", "p1", "p2", "m", "pm")
    lemma.fact("(<= (absr (- BC bc)) (* (+ (* 2.0 u) (* u u)) (absr bc)))")
    lemma.fact("(<= (absr (- CB cb)) (* (+ (* 2.0 u) (* u u)) (absr cb)))")
    lemma.multiplied("p1", "bc")
    lemma.multiplied("p2", "cb")
    lemma.rounded("m", "(- p1 p2)")
    lemma.rounded("pm", "(+ (absr p1) (absr p2))")
    exact = "(- BC CB)"
    lemma.goal(f"(<= (absr (- {exact} m)) (+ (* {ALPHA} pm) (* 3.0 h)))")
    lemma.goal(f"(<= (absr {exact}) (+ (* (+ 1.0 (* 6.0 u) (* 14.0 u u)) pm) (* 3.0 h)))")
    lemma.goal("(<= (absr m) (* (+ 1.0 (* 2.0 u) (* 2.0 u u)) pm))")
    lemma.goal("(<= pm (+ (* (+ 1.0 (* 3.0 u)) (+ (absr bc) (absr cb))) (* 3.0 h)))")
    v = program.circle_first()
    lemma.code({"bc": v["ap1"], "cb": v["ap2"], "p1": v["ap1"], "p2": v["ap2"], "m": v["am"],
                "pm": v["apm"]})
    return lemma


def affine_product():
    """|a| <= alpha A + beta and |b| <= gamma B + delta give
    |a b| <= alpha gamma A B + alpha delta A + beta gamma B + beta delta."""
    lemma = Lemma("affine product", "a", "b", "A", "B", "alpha", "beta", "gamma", "delta")
    for value in ("A", "B", "alpha", "beta", "gamma", "delta"):
        lemma.fact(f"(>= {value} 0.0)")
    lemma.fact("(<= (absr a) (+ (* alpha A) beta))")
    lemma.fact("(<= (absr b) (+ (* gamma B) delta))")
    lemma.goal("(<= (absr (* a b))"
               " (+ (* alpha gamma A B) (* alpha delta A) (* beta gamma B) (* beta delta)))")
    return lemma


def perturbation():
    """L M - l m = (L - l) M + l (M - m): a term's error split between its lift and its minor."""
    lemma = Lemma("perturbation", "L", "M", "l", "m")
    lemma.goal("(= (- (* L M) (* l m)) (+ (* (- L l) M) (* l (- M m))))")
    return lemma


def circle_corner():
    """A corner's term of the first stage, t = fl(lift m), and of its permanent,
    q = fl(pm lift), against the exact term T. LP = lift pm and LM = lift m are the exact
    products; d1 = (L - lift) M and d2 = lift (M - m) the perturbation's parts, bounded by the
    affine product lemma from circle lift's and circle minor's conclusions."""
    lemma = Lemma("circle corner", "T", "LM", "LP", "d1", "d2", "t", "q", "lift", "pm")
    alpha = ALPHA
    big = "(+ 1.0 (* 6.0 u) (* 14.0 u u))"
    lemma.fact("(>= lift 0.0)")
    lemma.fact("(>= pm 0.0)")
    lemma.fact("(>= LP 0.0)")
    lemma.fact("(= (- T LM) (+ d1 d2))")
    lemma.fact(f"(<= (absr d1) (+ (* {alpha} {big} LP) (* {alpha} 3.0 h lift)"
               f" (* 3.0 h {big} pm) (* 9.0 h h)))")
    lemma.fact(f"(<= (absr d2) (+ (* {alpha} LP) (* 3.0 h lift)))")
    lemma.fact("(<= (absr LM) (* (+ 1.0 (* 2.0 u) (* 2.0 u u)) LP))")
    lemma.multiplied("t", "LM")
    lemma.multiplied("q", "LP")
    lemma.goal("(<= (absr (- T t))"
               " (+ (* (+ (* 9.0 u) (* 50.0 u u)) LP) (* 4.0 h lift) (* 4.0 h pm) (* 2.0 h)))")
    lemma.goal("(<= (absr t) (+ (* (+ 1.0 (* 3.0 u) (* 5.0 u u)) LP) h))")
    lemma.goal("(<= LP (+ (* (+ 1.0 u (* 2.0 u u)) q) h))")
    lemma.goal("(>= q 0.0)")
    v = program.circle_first()
    lemma.code({"LM": v["at"], "t": v["at"], "LP": v["aq"], "q": v["aq"]})
    return lemma


def circle_permanent_minor():
    """A corner's permanent part against the other two lifts:
    pm <= (1 + 7u) (lift_b + lift_c) / 2 + 7h, by amgm and circle lift."""
    lemma = Lemma("permanent minor", "p1", "p2", "bc", "cb", "pm", "xb", "yb", "xc", "yc", "lb",
                  "lc")
    for value in ("xb", "yb", "xc", "yc", "lb", "lc"):
        lemma.fact(f"(>= {value} 0.0)")
    lemma.fact("(<= (absr bc) (/ (+ xb yc) 2.0))")
    lemma.fact("(<= (absr cb) (/ (+ xc yb) 2.0))")
    lemma.fact("(<= (+ xb yb) (+ (* (+ 1.0 (* 3.0 u)) lb) (* 3.0 h)))")
    lemma.fact("(<= (+ xc yc) (+ (* (+ 1.0 (* 3.0 u)) lc) (* 3.0 h)))")
    lemma.fact("(<= pm (+ (* (+ 1.0 (* 3.0 u)) (+ (absr bc) (absr cb))) (* 3.0 h)))")
    lemma.goal("(<= pm (+ (* (/ (+ 1.0 (* 7.0 u)) 2.0) (+ lb lc)) (* 7.0 h)))")
    return lemma


def circle_filter():
    """The circle test's first stage: det = fl(fl(t_a + t_b) + t_c), permanent
    fl(fl(q_a + q_b) + q_c), lifts fl(fl(l_a + l_b) + l_c), and the bound
    fl(fl(CircleBoundA permanent) + fl(CircleBoundUnderflow fl(lifts + 1))), with each corner's
    facts from circle corner and permanent minor."""
    names = []
    for corner in "abc":
        names += [f"T{corner}", f"t{corner}", f"LP{corner}", f"q{corner}", f"l{corner}",
                  f"pm{corner}"]
    lemma = Lemma("circle filter", *names, "s", "det", "ps", "perm", "s1", "lifts", "w", "k1",
                  "k2", "bound")
    for corner, (one, two) in zip("abc", ("bc", "ca", "ab")):
        lemma.fact(f"(>= l{corner} 0.0)")
        lemma.fact(f"(>= pm{corner} 0.0)")
        lemma.fact(f"(>= LP{corner} 0.0)")
        lemma.fact(f"(>= q{corner} 0.0)")
        lemma.fact(f"(<= (absr (- T{corner} t{corner})) (+ (* (+ (* 9.0 u) (* 50.0 u u))"
                   f" LP{corner}) (* 4.0 h l{corner}) (* 4.0 h pm{corner}) (* 2.0 h)))")
        lemma.fact(f"(<= (absr t{corner}) (+ (* (+ 1.0 (* 3.0 u) (* 5.0 u u)) LP{corner}) h))")
        lemma.fact(f"(<= LP{corner} (+ (* (+ 1.0 u (* 2.0 u u)) q{corner}) h))")
        lemma.fact(f"(<= pm{corner} (+ (* (/ (+ 1.0 (* 7.0 u)) 2.0) (+ l{one} l{two}))"
                   f" (* 7.0 h)))")
    lemma.rounded("s", "(+ ta tb)")
    lemma.rounded("det", "(+ s tc)")
    lemma.rounded("ps", "(+ qa qb)")
    lemma.rounded("perm", "(+ ps qc)")
    lemma.rounded("s1", "(+ la lb)")
    lemma.rounded("lifts", "(+ s1 lc)")
    lemma.rounded("w", "(+ lifts 1.0)")
    lemma.multiplied("k1", f"(* {rational(CIRCLE_BOUND_A)} perm)")
    lemma.multiplied("k2", f"(* {rational(CIRCLE_BOUND_UNDERFLOW)} w)")
    lemma.rounded("bound", "(+ k1 k2)")
    lemma.exceeds("det", "bound")
    lemma.same_sign("(+ Ta Tb Tc)", "det")
    v = program.circle_first()
    binding = {name: v[name] for name in ("s", "det", "ps", "perm", "s1", "lifts", "w", "k1",
                                          "k2", "bound")}
    for corner in "abc":
        binding.update({f"t{corner}": v[f"{corner}t"], f"q{corner}": v[f"{corner}q"],
                        f"l{corner}": v[f"{corner}lift"]})
    lemma.code(binding)
    return lemma


def nonfinite(name):
    return f"(or (fp.isInfinite {name}) (fp.isNaN {name}))"


def propagation():
    """IEEE 754: an infinity or a NaN never yields a finite sum, difference, product or
    magnitude, so one anywhere in the bound's operands leaves the bound infinite or NaN."""
    lemma = Floating("propagation", "a", "b")
    lemma.fact(nonfinite("a"))
    for operation in ("(fp.add RNE a b)", "(fp.add RNE b a)", "(fp.sub RNE a b)",
                      "(fp.sub RNE b a)", "(fp.mul RNE a b)", "(fp.mul RNE b a)", "(fp.abs a)"):
        lemma.goal(nonfinite(operation))
    return lemma


def comparison():
    """IEEE 754: |x| > bound fails for a bound infinite and not below zero, or NaN, whatever x is:
    so a non-finite bound never decides."""
    lemma = Floating("comparison", "x", "bound")
    lemma.fact(nonfinite("bound"))
    lemma.fact("(not (fp.isNegative bound))")
    lemma.goal("(not (fp.gt (fp.abs x) bound))")
    return lemma


def rounding():
    """IEEE 754: the rounded difference of two finite doubles is at most their rounded sum of
    magnitudes, |fl(p1 - p2)| <= fl(|p1| + |p2|), overflow included: the minor never exceeds its
    permanent part. Monotone rounding gives the rest alike, |fl(lift m)| <= fl(lift pm) and the
    determinant's partial sums against the permanent's, so where the permanent is finite the
    determinant is, and the determinant needs no overflow check of its own."""
    lemma = Floating("minor and permanent", "p1", "p2")
    lemma.fact(f"(not {nonfinite('p1')})")
    lemma.fact(f"(not {nonfinite('p2')})")
    lemma.goal("(fp.leq (fp.abs (fp.sub RNE p1 p2)) (fp.add RNE (fp.abs p1) (fp.abs p2)))")
    return lemma


def mutant(lemma_maker, **constants):
    """A lemma rebuilt with some of the code's constants replaced."""

    def make():
        saved = {name: globals()[name] for name in constants}
        globals().update(constants)
        try:
            return lemma_maker()
        finally:
            globals().update(saved)

    return make


MUTANTS = (("orient without u_N", mutant(orient_filter, NORMAL=0.0)),
           ("orient bound 2.9u", mutant(orient_filter, ORIENT_BOUND_A=2.9 * lemma.EPSILON)),
           ("circle without the underflow term",
            mutant(circle_filter, CIRCLE_BOUND_UNDERFLOW=0.0)),
           ("circle underflow term 2^-1072", mutant(circle_filter,
                                                     CIRCLE_BOUND_UNDERFLOW=2.0 ** -1072)),
           ("circle bound 9u", mutant(circle_filter, CIRCLE_BOUND_A=9.0 * lemma.EPSILON)))

LEMMAS = (tails, orient_filter, amgm, circle_lift, circle_minor, affine_product, perturbation,
          circle_corner, circle_permanent_minor, circle_filter, propagation, comparison,
          rounding)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, mutants=MUTANTS))
