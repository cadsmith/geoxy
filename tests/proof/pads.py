#!/usr/bin/env python3
"""Machine-checks with Z3 the pads source/extents.cpp grows a Bézier curve's box by, the grid
padding source/grid.cpp grows every box of a grid by, and what of source/hodograph.cpp's
Gauss-Legendre bound can be proved.

    python3 -I tests/proof/pads.py [z3]

A rounding is lemma.py's: |r - z| <= u |z| and |r - z| <= u |r| for a sum, a difference, or a
product clear of underflow; a product that underflows is off by h = 2^-1075 at the most. Each
lemma states the operations of the code as written, one fact a rounding; they are not bound to
program.py's graphs, so translation.py does not check them against the compiler's IR. Each mutant
shrinks a bound or drops a term the code keeps, and Z3 must refute it: so the proofs need every
part of the bounds they are given. The exit status is 1 when a lemma is not proved or a mutant
is.

- Wandering, extents.cpp: one level of de Casteljau's algorithm by Lerp, either branch, carries
  its inputs' error E to E + 3 u (S + E) Grow for coordinates within S, and h Grow more where
  its product underflows; d levels from exact control points then stay within the pad (4 d + 1)
  u S Grow and d 2^-1074, the pad's own rounding taken, for the degrees At works out, 2 and 3.
  The d 2^-1074 was not in the pad: below 2^-1025 the pad rounded to nothing, and At's point
  passed the box by a step of a double, as at a quadratic of a few steps of 2^-1074.
- Straying, extents.cpp: de Boor's fractions round, a mix's two factors within 6.0001 u of the
  exact ones and summing within u of 1, and over d rounds the coefficients over the poles are off
  by 6.0002 d u in all; the point by (w_max / w_min) times that times the poles' spread, twice
  under 12.5 d of it for At's point and Reckoned's. That the point is off by sum w_i D_i (x_i - x)
  / sum w_i c_i is algebra, and the mixes' own roundings (4 d + 2) steps of the largest pole by
  Farouki and Rajan's argument: both by hand. The pad left the fractions out, and At's point came
  345 steps of a double past the box of a spline whose weights differ 7e5 times.
- Gridwork::Grown, grid.cpp: Overshot and Undershot, in IEEE 754 itself, step past every double they
  are given, so a box spread holds the box grown exactly; Above is monotone; and two boxes grown
  that do not overlap lie farther apart than Reached's clearance, as Intersect works the gap out, for a
  curve of the drawing at its tolerance, at a finer one, and probed at a coarser one.
- GaussLegendre, hodograph.cpp: the Bernstein ellipse's points lie within its half axis of its
  middle, so the speed on it is at most the most GaussLegendre works out, and the loop stops where
  the bound with every rounding taken is within the share, while the bound from GaussLegendre's
  most and rho is within 5/4 of the exact one, a hypothesis. Trefethen's theorem 19.3 itself is
  cited, not proved.

Left to the comments, worked by hand or measured: Turnpoint's argument that At at a computed root
is within second order of the true extreme, the roots found as the stable quadratic formula
finds them; and Bowing's pad, by measure, as At's bulge rests on the library's sine and
arctangent; and that Roots finds a cubic's turns at every scale, its coefficients scaled by a
power of two past 2^500 or under 2^-500: below 2^-540 they underflowed, and a box missed its
curve's bulge.
"""

import fractions
import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from expand import Fraction, Poly, U, rounding  # pylint: disable=wrong-import-position
# pylint: disable-next=wrong-import-position
from lemma import Floating, Lemma, rational, run  # pylint: disable=wrong-import-position

GROW = 1.0 + 2.0 ** -40
NODE = Fraction(60001, 10000) * U
EPSILON = 2.0 ** -53


def rounds(lemma, result, exact):
    """result is exact rounded to nearest: a sum, a difference, or a product clear of
    underflow."""
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {exact})))")
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {result})))")


def multiplies(lemma, result, exact):
    """result is the product exact rounded to nearest, gradual underflow allowed: relative clear
    of underflow, within h below it, and zero for zero."""
    relative = (f"(and (>= (absr {exact}) uN)"
                f" (<= (absr (- {result} {exact})) (* u (absr {exact})))"
                f" (<= (absr (- {result} {exact})) (* u (absr {result}))))")
    absolute = (f"(and (< (absr {exact}) uN) (<= (absr (- {result} {exact})) h)"
                f" (<= (absr {result}) uN))")
    lemma.fact(f"(or {relative} {absolute})")
    lemma.fact(f"(=> (= {exact} 0.0) (= {result} 0.0))")


def lerp_shape(half=True):
    """Lerp, one level of de Casteljau, from x0 toward x1 by a fraction f up to a half: c =
    rnd(x1 - x0), then p = rnd(c f) and r = rnd(x0 + p). Past a half Lerp takes c' = rnd(x0 -
    x1) and r = rnd(x1 - rnd(c' (1 - f))), 1 - f exact by Sterbenz: rounding to nearest being
    symmetric, that is the same steps from x1 toward x0 by 1 - f, again up to a half. With x0
    and x1 within T, the product's operand c f is within (1 + u) T, the rounding of the
    difference moves it by u T at the most, and the exact place x0 + f (x1 - x0) is within T."""
    lemma = Lemma("lerp shape", "x0", "x1", "f", "T", "c")
    lemma.fact("(<= 0.0 f)")
    lemma.fact("(<= f 0.5)" if half else "(<= f 1.0)")
    lemma.fact("(<= (absr x0) T)")
    lemma.fact("(<= (absr x1) T)")
    rounds(lemma, "c", "(- x1 x0)")
    lemma.goal("(<= (absr (* c f)) (* (+ 1.0 u) T))")
    lemma.goal("(<= (absr (* f (- c (- x1 x0)))) (* u T))")
    lemma.goal("(<= (absr (+ x0 (* f (- x1 x0)))) T)")
    return lemma


def lerp_round(share=3.0, underflow=True):
    """Lerp's two last roundings, the product's operand P = c f and its exact part m = x0 + f
    (x1 - x0) as lerp shape bounds them, P = (m - x0) + X for X = f (c - (x1 - x0)): p = rnd(P)
    and r = rnd(x0 + p) are within 3 u T Grow of m, and h Grow more where the product
    underflows: u T for the sum, u (1 + u) T for the product, u T for the difference, and their
    products' u^2 terms, which Grow covers."""
    lemma = Lemma("lerp round", "x0", "m", "P", "X", "p", "r", "T")
    lemma.fact("(>= T 0.0)")
    lemma.fact("(<= (absr P) (* (+ 1.0 u) T))")
    lemma.fact("(<= (absr X) (* u T))")
    lemma.fact("(<= (absr m) T)")
    lemma.fact("(= P (+ (- m x0) X))")
    if underflow:
        multiplies(lemma, "p", "P")
    else:
        rounds(lemma, "p", "P")
    rounds(lemma, "r", "(+ x0 p)")
    slack = f" (* h {rational(GROW)})" if underflow else ""
    lemma.goal(f"(<= (absr (- r m)) (+ (* {share} u T {rational(GROW)}){slack}))")
    return lemma


def lerp_carry(carried=True):
    """Lerp's exact weights carry its inputs' error as a convex combination does: x0 and x1
    within E of A0 and A1, x0 + f (x1 - x0) is within E of A0 + f (A1 - A0) for f from 0 to 1.
    With |A0| and |A1| at most S, |x0| and |x1| are at most T = S + E: so one level is within E +
    3 u (S + E) Grow, and h Grow more where its product underflows."""
    lemma = Lemma("lerp carry", "x0", "x1", "A0", "A1", "f", "E")
    lemma.fact("(<= 0.0 f)")
    lemma.fact("(<= f 1.0)")
    lemma.fact("(<= (absr (- x0 A0)) E)")
    lemma.fact("(<= (absr (- x1 A1)) E)")
    bound = "E" if carried else "(* 0.5 E)"
    lemma.goal(f"(<= (absr (- (+ x0 (* f (- x1 x0))) (+ A0 (* f (- A1 A0))))) {bound})")
    return lemma


def wandering(degree=3, steps=None):
    """Wandering: d levels of lerp from exact control points, E_0 = 0 and E_(k+1) = E_k + 3 u (S
    + E_k) Grow + h Grow, stay within the pad Wandering works out, rnd(rnd(S rnd((4 d + 1) u
    Grow)) + d 2^-1074), the constant exact and each rounding down by u at the most. S is the
    largest control coordinate along the axis, a bound on every exact value de Casteljau's levels
    reach, each a convex combination of the controls."""
    steps = 4 * degree + 1 if steps is None else steps
    grow = rational(GROW)
    names = [f"E{level}" for level in range(degree + 1)]
    lemma = Lemma(f"wandering degree {degree}", "S", "pad", *names)
    lemma.fact("(>= S 0.0)")
    lemma.fact("(= E0 0.0)")
    for level in range(degree):
        here, after = names[level], names[level + 1]
        lemma.fact(f"(<= {after} (+ {here} (* 3.0 u (+ S {here}) {grow}) (* h {grow})))")
    constant = fractions.Fraction(steps) * fractions.Fraction(EPSILON) * fractions.Fraction(GROW)
    if fractions.Fraction(float(constant)) != constant:
        raise ValueError("the pad's constant is not a double")
    lemma.fact(f"(>= pad (* (+ (* S {rational(float(constant))} (- 1.0 u)) (* {degree}.0 2.0 h))"
               f" (- 1.0 u)))")
    lemma.goal(f"(<= {names[-1]} pad)")
    return lemma


def above_bits(name):
    """Above(magnitude) as real.hpp works it: the double whose bits follow the magnitude's,
    less the magnitude."""
    return (f"(fp.sub RNE ((_ to_fp 11 53) (bvadd (fp.to_ieee_bv {name}) #x0000000000000001))"
            f" {name})")


def beyond():
    """IEEE 754: Overshot(s) = s + Above(|s|) is above s for every finite double s, or not
    finite. So it is at least the double after s, which every number s is the rounding of lies
    below: Spread's box holds the box grown exactly, Undershot likewise from beneath."""
    lemma = Floating("beyond and below", "s")
    lemma.fact("(not (or (fp.isInfinite s) (fp.isNaN s)))")
    step = above_bits("(fp.abs s)")
    up = f"(fp.add RNE s {step})"
    down = f"(fp.sub RNE s {step})"
    lemma.goal(f"(or (fp.isInfinite {up}) (fp.isNaN {up}) (fp.gt {up} s))")
    lemma.goal(f"(or (fp.isInfinite {down}) (fp.isNaN {down}) (fp.lt {down} s))")
    return lemma


def monotone():
    """IEEE 754: Above is monotone on finite magnitudes, so 16 Above of the larger of two sizes
    is the larger of their 16 Above: Reached's grain is at most the larger of the two curves'
    grains, each read from the larger of its box's reach and its numbers' size."""
    lemma = Floating("above monotone", "m", "n")
    for name in ("m", "n"):
        lemma.fact(f"(not (or (fp.isInfinite {name}) (fp.isNaN {name})))")
        lemma.fact(f"(not (fp.isNegative {name}))")
    lemma.fact("(fp.leq m n)")
    lemma.fact(f"(not (fp.isInfinite {above_bits('n')}))")
    lemma.goal(f"(fp.leq {above_bits('m')} {above_bits('n')})")
    return lemma


def apart(ratio=2.0):
    """Grown and Probed against Reached: one box's right edge x1, the other's left edge x2. The
    first curve is grown by m1 >= 32 (L1 + D)(1 - u)^2, L1 the larger of the tolerance t it is
    met within and its grain, D = t - T >= 0 the excess over the drawing's tolerance T, Probed's
    sum rounded; the second by 32 L2T, its L at T, and its L at t is at most L2T + D. A pair of
    the drawing met at T or a finer tolerance is the case D = 0, L1 its L at T. Spread's edges
    lie past the exact ones, x1 + m1 < X1 and X2 < x2 - 32 L2T, and the boxes grown do not
    overlap, X1 < X2. Reached's clearance c, 4 or 16 of its length, is at most 16 of the larger
    L at t. Then the gap Intersect works out, g = rnd(x2 - x1), is more than c: the pair is far, and
    meets in nothing. Twice the clearance is covered without Probed's 32 D: D is at most t, so
    at most L1, and half of the first margin takes it."""
    lemma = Lemma("apart", "x1", "x2", "X1", "X2", "L1", "L2", "L2T", "D", "t", "m1", "c", "g")
    lemma.fact("(> L1 0.0)")
    lemma.fact("(> L2T 0.0)")
    lemma.fact("(>= D 0.0)")
    lemma.fact("(>= L1 t)")
    lemma.fact("(>= L2 t)")
    lemma.fact("(<= D t)")
    lemma.fact("(<= L2 (+ L2T D))")
    lemma.fact("(>= m1 (* 32.0 (+ L1 D) (- 1.0 u) (- 1.0 u)))")
    lemma.fact("(> X1 (+ x1 m1))")
    lemma.fact("(< X2 (- x2 (* 32.0 L2T)))")
    lemma.fact("(< X1 X2)")
    lemma.fact("(> c 0.0)")
    lemma.fact(f"(or (<= (* {ratio} c) (* 32.0 L1)) (<= (* {ratio} c) (* 32.0 L2)))")
    rounds(lemma, "g", "(- x2 x1)")
    lemma.goal("(> g c)")
    return lemma


def ellipse(reach=True):
    """The Bernstein ellipse of a piece of half width w about its middle m, z = m + (w / 2)(rho
    e^(i a) + e^(-i a) / rho): |z - m|^2 = (w / 2)^2 (rho^2 + rho^-2 + 2 cos 2a), at most (w (rho +
    1 / rho) / 2)^2, GaussLegendre's reach squared; k stands for cos 2a."""
    lemma = Lemma("bernstein ellipse", "w", "rho", "k", "d")
    lemma.fact("(> w 0.0)")
    lemma.fact("(> rho 1.0)")
    lemma.fact("(<= (absr k) 1.0)")
    lemma.fact("(= (* d rho rho)"
               " (* (/ 1.0 4.0) w w (+ (* rho rho rho rho) 1.0 (* 2.0 k rho rho))))")
    most = "(* (/ 1.0 4.0) w w (+ (* rho rho rho rho) 1.0 (* 2.0 rho rho)))" if reach else \
        "(* (/ 1.0 4.0) w w (* rho rho rho rho))"
    lemma.goal(f"(<= (* d rho rho) {most})")
    return lemma


def speed(geometric=True):
    """The speed's continuation over the ellipse, Lead times the square root of |z - r| |z -
    conj r| for each root r: each factor is at most |m - r| + a, a the half axis, as m is real so
    |m - r| = |m - conj r|, and |z - m| <= a; so the speed is at most GaussLegendre's most."""
    lemma = Lemma("speed on the ellipse", "p", "q", "s", "a", "D")
    lemma.fact("(>= p 0.0)")
    lemma.fact("(>= q 0.0)")
    lemma.fact("(>= s 0.0)")
    lemma.fact("(<= p (+ a D))")
    lemma.fact("(<= q (+ a D))")
    lemma.fact("(= (* s s) (* p q))" if geometric else "(= s (+ p q))")
    lemma.goal("(<= s (+ a D))")
    return lemma


SHARE_GROWTH = fractions.Fraction(1 + 2 ** -44)


def share(slack=0.75):
    """GaussLegendre's loop: it stops at the first pairs whose bound, as worked out, is within
    three quarters of the share, or at the last rule. Each of the bound's roundings, at most 160 in
    all over 20 pairs, is within (1 + u), (1 + u)^160 under 1 + 2^-44; with GaussLegendre's most
    and rho each such that the exact bound is at most 5/4 of the bound worked out from them
    exactly, a hypothesis, the exact bound is within the share."""
    lemma = Lemma("ruled share", "computed", "exact", "share")
    lemma.fact("(> share 0.0)")
    lemma.fact(f"(<= computed (* {slack} share))")
    lemma.fact(f"(<= exact (* computed (/ 5.0 4.0) {rational(float(SHARE_GROWTH))}))")
    lemma.goal("(<= exact share)")
    return lemma


def growth_check():
    """(1 + u)^160 is under 1 + 2^-44, exactly."""
    held = (1 + fractions.Fraction(EPSILON)) ** 160 < SHARE_GROWTH
    return held, "(1 + u)^160 < 1 + 2^-44"


def unpadded():
    """lerp round, underflow allowed but its h dropped from the bound."""
    lemma = lerp_round()
    lemma.goals = [goal.replace(f" (* h {rational(GROW)})))", "))") for goal in lemma.goals]
    return lemma


def quadratic():
    """wandering for degree 2."""
    return wandering(degree=2)


def cubic():
    """wandering for degree 3."""
    return wandering(degree=3)


def fraction_node():
    """Straying, a mix of de Boor's algorithm: the fraction a = (t - s) / (e - s) for s <= t <= e
    rounds in three, a' = a (1 + d1)(1 + d3) / (1 + d2), and the rest r' = rnd(1 - a') in one more.
    Then |a' - a| <= c a for c the expansion's bound of (1 + d1)(1 + d3) / (1 + d2) - 1, and |a' -
    a| + |r' - (1 - a)| <= 2 c a + u (1 - a + c a), at most 6.0001 u over 0 <= a <= 1, the larger
    of its ends; and r' + a' = 1 + d4 (1 - a') is within u of 1."""
    ratio = rounding("d1") * rounding("d3") * (1 + Poly.var("d2")) - 1
    bounds = {"d1": U, "d3": U, "d2": U / (1 - U)}
    spread = ratio.bound(bounds)
    node = max(U, 2 * spread + U * spread)
    return node <= NODE, f"c {float(spread / U):.6f} u, node {float(node / U):.6f} u"


def fraction_path():
    """Straying, the fractions over d rounds: each node's coefficients c' over the poles, all at
    least zero, are r' times the left node's and a' times the right's; the exact ones sum to 1 at
    every node. So the sum D of |c' - c| grows by at most (r' + a') D + 6.0001 u a round, D_d <=
    6.0001 u ((1 + u)^d - 1) / u, and the c' sum to at least (1 - u)^d. The point x' = sum w_i c'_i
    x_i / sum w_i c'_i is off from the exact x by sum w_i (c'_i - c_i)(x_i - x) / sum w_i c'_i, at
    most (w_max / w_min) D_d spread / (1 - u)^d; At's point at a parameter, and Reckoned's at a
    turn, each so: twice that, the ratio and the spread each rounded, is under 12.5 d (w_max /
    w_min) spread u for every degree to 31."""
    worst = Fraction(0)
    for degree in range(1, 32):
        grown = NODE * ((1 + U) ** degree - 1) / U
        moved = 2 * grown / (1 - U) ** degree * (1 + U) ** 3
        worst = max(worst, moved / (degree * U))
    return worst <= Fraction(25, 2), f"{float(worst):.4f} d u (w_max / w_min) spread"


MUTANTS = (("lerp shape past a half", lambda: lerp_shape(half=False)),
           ("lerp round at 2.5 u", lambda: lerp_round(share=2.5)),
           ("lerp round without underflow's h", unpadded),
           ("lerp carry at half the error", lambda: lerp_carry(carried=False)),
           ("wandering cubic at 8 steps", lambda: wandering(degree=3, steps=8)),
           ("apart at a clearance of 32 L", lambda: apart(ratio=1.0)),
           ("ellipse within rho^2 alone", lambda: ellipse(reach=False)),
           ("speed by the sum of factors", lambda: speed(geometric=False)),
           ("ruled share at 0.8", lambda: share(slack=0.8)))

LEMMAS = (lerp_shape, lerp_round, lerp_carry, quadratic, cubic, beyond, monotone, apart, ellipse,
          speed, share)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, checks=(growth_check, fraction_node, fraction_path),
                 mutants=MUTANTS))
