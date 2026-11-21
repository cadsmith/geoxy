#!/usr/bin/env python3
"""Machine-checks with Z3 the bounds source/abut.cpp's two passes rest on, Squarely for a level
line and an upright one, Rimmed for a level or upright line and an arc: that every test they make
on a number Intersect works out by its own roundings is told the way Intersect tells it.

    python3 -I tests/proof/abut.py [z3]

A rounding is lemma.py's: |r - z| <= u |z| and |r - z| <= u |r| for a sum, a difference, or a
product clear of underflow; a product, quotient or square that may underflow is off by h =
2^-1075 more. A square root is the real root rounded once. u is 2^-53, S the size each pass works
out, L the length. Squarely's lines are longer than 4 L, past 2^-108 as L is at least Fine, 2^-110,
and S is below 2^47 L, 2^97 as L is at most Coarse: every product of a gap and a way is a normal
double. Rimmed's radius is past 1024 L and 2^-32 S, so S is past 2^-100 and every product of a
line's way with itself a normal double too; S is max(the line's reach, |p| + r, |q| + r) rounded,
so |p| + r and |q| + r are at most S (1 + 2 u).

Squarely: its gap test past 4 L (1 + 2^-40) passes Intersect's Beside on the same rounded gaps (beside,
level length), and Away finds such an end farther than 4 L (away); Endly's point, the first line's
start moved along its way by a quotient of rounded products, is within 11.01 u S of the carriers'
crossing (quotient terms, quotient rounding, quotient chain, moved), so the end on the other's
carrier lies within L / 2 of it and Away finds it within as much of the other line (near end,
along); each place lies within 0 and 1 by 2^-46, and Intersect's rounded place inside [0, 1] (placed);
and the lines are never shallow at a tolerance's angle within Steep / 2 (steep lines).

Rimmed: Roundly's foot is within 11.01 u S of (p, c) (the same chain), the carrier's distance
from the center within 17.1 u S of |g| (length sum, length main, scaled square, scaled sum, scaled
root, length tiny: Length within 2.01 u of the exact length and 2^-540; hypotenuse, distance); the
test r^2 - g^2 > r^2 / 32 made here (rest test) passes Roundly's own: no touch, a sine past Steep,
the carrier cutting the circle (cut deep); the radicand is within 37.3 u S r of r^2 - g^2 (squares,
radicand), its root within 106.6 u S (root, under the 108.3 claimed), each crossing and its place
along the line, taken along the way, within 128 u S of the true one (level length, half, scaled,
crossing, crossing point; 124 u S at most). Spotted's root is within a fifth of offRoot of the true
one, and its crossing within 0.6 of the bound it gives Certain, as Certain asks, that bound past
2^-99 of the value, the double within 1.01 u S of the true crossing certified or not (gap high,
spotted rest, spotted bounds, spotted scale, spotted root, spotted share, spotted reach, spotted
size, spotted part, spotted value, spotted sum, spotted near). Then every number Roundly compares
is told within the slack 2^-45 S here: an end's distance from a crossing against L / 2 (told near,
told far, triangle, within), a line's end's from the circle against 4 L (off told), an arc's end's
from the line against 4 L (along, away far), a crossing's place along the line (inside), and its
turn from an arc's end (sine product, cut, toward, products, swept).

Assumed, from source/lean.hpp's comments, which tests/proof proves elsewhere: lean Times within
2^-102 |a b| and Plus within 2^-104 (|a| + |b|) of their high parts, so r^2 - g^2 within 2^-101
(r_H^2 + g_H^2); lean Root within 2^-103 s of the root of its input, its low part at most 1.52 u s;
TwoSum and TwoProduct exact. std::sin at most 1 in magnitude, and below 0.0626 at 1/16. Left to
abut.cpp's comments, by hand: which branch Intersect takes given these tests (Aligned finding no stretch
for two square lines, the order of hits, an end's kind), and that a coordinate's mirror and the
swap of x and y round alike, which takes each lemma's case to the others.

Each lemma states the operations as written, one fact a rounding; they are not bound to
program.py's graphs. Each mutant shrinks a bound or drops a term the code keeps, and Z3 must refute
it. The exit status is 1 when a lemma is not proved or a mutant is.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from lemma import Lemma, rational, run  # pylint: disable=wrong-import-position

MARGIN = 1.0 + 2.0 ** -40
TIGHT = 1.0 - 2.0 ** -48
LOOSE = 1.0 + 2.0 ** -48
GRAIN = 2.0 ** -98
SLACK = 2.0 ** -45
FOOT = "11.01"
OFF = "17.1"
ROOT = "106.6"
CROSSING = "128.0"
CHAIN = "4.0002"
ACROSS = "21.1"


def power(exponent):
    """2^exponent as an exact SMT-LIB rational."""
    return rational(2.0 ** exponent)


def rounds(lemma, result, exact):
    """result is exact rounded to nearest: a sum, a difference, or a product clear of
    underflow."""
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {exact})))")
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {result})))")


def multiplies(lemma, result, exact):
    """result is the product exact rounded to nearest, underflow allowed."""
    lemma.fact(f"(<= (absr (- {result} {exact})) (+ (* u (absr {exact})) h))")


def divides(lemma, result, top, bottom):
    """result is top / bottom rounded to nearest, underflow allowed, bottom not zero."""
    lemma.fact(f"(not (= {bottom} 0.0))")
    lemma.fact(f"(<= (absr (- (* {result} {bottom}) {top})) "
               f"(+ (* u (absr {top})) (* h (absr {bottom}))))")


def rooted(lemma, result, real, radicand):
    """result is the square root of radicand rounded to nearest, real the root itself."""
    lemma.fact(f"(>= {real} 0.0)")
    lemma.fact(f"(= (* {real} {real}) {radicand})")
    lemma.fact(f"(<= (absr (- {result} {real})) (* u {real}))")


def measured(lemma, result, exact):
    """result is Length of a vector whose exact length is exact, as length main and length scaled
    prove."""
    lemma.fact(f"(<= (absr (- {result} {exact})) (+ (* 2.01 u {exact}) {power(-540)}))")


def level_length():
    """The length of a way (d, 0), or the root of a line's square: sqrt(rnd(d d)) rounded, d^2 a
    normal double, within 1.51 u |d| of |d|. (It is |d| exactly, but this is all that is used.)"""
    lemma = Lemma("level length", "d", "square", "w", "length")
    lemma.fact(f"(>= (absr d) {power(-108)})")
    rounds(lemma, "square", "(* d d)")
    rooted(lemma, "length", "w", "square")
    lemma.goal("(<= (absr (- length (absr d))) (* 1.51 u (absr d)))")
    return lemma


def beside(margin=MARGIN):
    """Squarely's gap test, |g| past clear = rnd(4 L Margin), and Intersect's Beside on the same
    rounded gap g and the other's way d: P = rnd(g d) against least = rnd(4 L |d|'), |d|' the
    way's length as Length gives it. P is past least, of g's sign times d's, for both ends: so
    where Squarely tells a line beside the other's carrier, Beside does."""
    lemma = Lemma("beside", "L", "g", "d", "clear", "P", "length", "least")
    lemma.fact("(> L 0.0)")
    lemma.fact("(> (absr d) 0.0)")
    rounds(lemma, "clear", f"(* 4.0 L {rational(margin)})")
    lemma.fact("(> (absr g) clear)")
    rounds(lemma, "P", "(* g d)")
    lemma.fact("(<= (absr (- length (absr d))) (* 1.51 u (absr d)))")
    rounds(lemma, "least", "(* 4.0 L length)")
    lemma.goal("(> (absr P) least)")
    return lemma


def away(margin=MARGIN):
    """Intersect's Away of an end from the other's run, its gap from the carrier g as Squarely rounds
    it: the distance to its foot is Length of (a, g), a the gap along, whose square is at least
    g^2 however the foot rounds. Past clear, Away is past 4 L: the end is far."""
    lemma = Lemma("away", "L", "g", "clear", "sa", "sb", "s", "w", "length")
    lemma.fact("(> L 0.0)")
    rounds(lemma, "clear", f"(* 4.0 L {rational(margin)})")
    lemma.fact("(> (absr g) clear)")
    lemma.fact("(>= sa 0.0)")
    rounds(lemma, "sb", "(* g g)")
    rounds(lemma, "s", "(+ sa sb)")
    rooted(lemma, "length", "w", "s")
    lemma.goal("(> length (* 4.0 L))")
    return lemma


def quotient_terms():
    """A place by a quotient of rounded products, n = rnd(Ar w) and D = rnd(d w), d w normal: Q
    = d n / D is within 2.0001 u |Ar| of Ar, besides an underflow of n over |w| >= 2^-110."""
    lemma = Lemma("quotient terms", "Q", "D", "d", "n", "Ar", "w")
    lemma.fact("(= (* Q D) (* d n))")
    lemma.fact("(not (= D 0.0))")
    multiplies(lemma, "n", "(* Ar w)")
    rounds(lemma, "D", "(* d w)")
    lemma.fact(f"(>= (absr w) {power(-110)})")
    lemma.goal(f"(<= (absr (- Q Ar)) (+ (* 2.0001 u (absr Ar)) (* {power(111)} h)))")
    return lemma


def quotient_rounding():
    """The quotient t = rnd(n / D) taken along d: k = d t is within u |Q| of Q = d n / D, besides
    an underflow, |d| <= 2^99."""
    lemma = Lemma("quotient rounding", "Q", "D", "d", "n", "t", "k")
    lemma.fact("(= (* Q D) (* d n))")
    divides(lemma, "t", "n", "D")
    lemma.fact("(= k (* d t))")
    lemma.fact(f"(<= (absr d) {power(99)})")
    lemma.goal(f"(<= (absr (- k Q)) (+ (* u (absr Q)) (* {power(99)} h)))")
    return lemma


def quotient_chain():
    """The place along d, k = d t, from the rounded gap Ar = rnd(A): within 4.0002 u |A| of the
    exact gap A, by quotient terms and quotient rounding."""
    lemma = Lemma("quotient chain", "A", "Ar", "Q", "k")
    rounds(lemma, "Ar", "A")
    lemma.fact(f"(<= (absr (- Q Ar)) (+ (* 2.0001 u (absr Ar)) (* {power(111)} h)))")
    lemma.fact(f"(<= (absr (- k Q)) (+ (* u (absr Q)) (* {power(99)} h)))")
    lemma.goal(f"(<= (absr (- k A)) (+ (* {CHAIN} u (absr A)) (* {power(112)} h)))")
    return lemma


def moved(bound=FOOT):
    """A start a0 moved by m = rnd(k) and rounded, out = rnd(a0 + m): within 11.01 u S of a0 + A,
    the point sought, a0 within S and a0 + A within S (1 + 2 u). Endly's point, and Roundly's foot
    of the center on the carrier, are each such a coordinate; the other is the carrier's own,
    exactly."""
    lemma = Lemma("moved", "S", "a0", "A", "k", "m", "out")
    lemma.fact(f"(>= S {power(-110)})")
    lemma.fact("(<= (absr a0) S)")
    lemma.fact("(<= (absr (+ a0 A)) (* S (+ 1.0 (* 2.0 u))))")
    lemma.fact(f"(<= (absr (- k A)) (+ (* {CHAIN} u (absr A)) (* {power(112)} h)))")
    multiplies(lemma, "m", "k")
    rounds(lemma, "out", "(+ a0 m)")
    lemma.goal(f"(<= (absr (- out (+ a0 A))) (* {bound} u S))")
    return lemma


def near_end():
    """Endly's distance from the end on the other's carrier, the true crossing (x, y), to its
    point: Length of (rnd(v), 0), v within 11.01 u S, at most L / 2 as L is past 2^-47 S. Away of
    that end from the other line is the Length of such a vector too, by along: not far."""
    lemma = Lemma("near end", "S", "L", "v", "rv", "length")
    lemma.fact(f"(>= S {power(-110)})")
    lemma.fact(f"(> L (* {power(-47)} S))")
    lemma.fact(f"(<= (absr v) (* {FOOT} u S))")
    rounds(lemma, "rv", "v")
    measured(lemma, "length", "(absr rv)")
    lemma.goal("(<= length (/ L 2.0))")
    return lemma


def along():
    """Intersect's foot on a run of way dx > 0 (a mirror takes the other sign), its place clamped to
    0 and 1: fR = rnd(a0 + rnd(dx tau)), tau the quotient's place k / dx or the end it passes. For
    a point a0 + A along the carrier, the true foot on the line from a0 to a0 + B is F = a0 +
    clamp(A, 0, B). fR is within 11.01 u S of the point where it lies on the line, and within 21.1
    u S of F wherever it lies, the point within S (1 + 5 u)."""
    lemma = Lemma("along", "S", "a0", "B", "dx", "A", "k", "mR", "fR")
    lemma.fact(f"(>= S {power(-110)})")
    lemma.fact("(> B 0.0)")
    lemma.fact("(> dx 0.0)")
    rounds(lemma, "dx", "B")
    lemma.fact("(<= (absr a0) S)")
    lemma.fact("(<= (absr (+ a0 B)) S)")
    lemma.fact("(<= (absr (+ a0 A)) (* S (+ 1.0 (* 5.0 u))))")
    lemma.fact(f"(<= (absr (- k A)) (+ (* {CHAIN} u (absr A)) (* {power(112)} h)))")
    lemma.fact("(or (and (>= k 0.0) (<= k dx) (<= (absr (- mR k)) (+ (* u (absr k)) h)))"
               " (and (< k 0.0) (= mR 0.0)) (and (> k dx) (= mR dx)))")
    rounds(lemma, "fR", "(+ a0 mR)")
    foot = "(+ a0 (ite (< A 0.0) 0.0 (ite (> A B) B A)))"
    lemma.goal(f"(=> (and (>= A 0.0) (<= A B)) (<= (absr (- fR (+ a0 A))) (* {FOOT} u S)))")
    lemma.goal(f"(<= (absr (- fR {foot})) (* {ACROSS} u S))")
    return lemma


def placed(margin=MARGIN):
    """Squarely's crossing on a line from a0 to a0 + B, B > 0 and its way dx = rnd(B) (a mirror
    takes the other sign): both ends' gaps g0 = rnd(-A) and g1 = rnd(B - A) past clear and of
    opposite signs, the place A / B lies within 0 and 1 by 2^-46, and Intersect's rounded place, k /
    dx, within [0, 1]."""
    lemma = Lemma("placed", "S", "L", "a0", "B", "dx", "A", "k", "clear", "g0", "g1")
    lemma.fact(f"(>= S {power(-110)})")
    lemma.fact(f"(> L (* {power(-47)} S))")
    lemma.fact("(> B 0.0)")
    lemma.fact("(> dx 0.0)")
    rounds(lemma, "dx", "B")
    lemma.fact("(<= (absr a0) S)")
    lemma.fact("(<= (absr (+ a0 B)) S)")
    lemma.fact("(<= (absr (+ a0 A)) S)")
    rounds(lemma, "clear", f"(* 4.0 L {rational(margin)})")
    rounds(lemma, "g0", "(- 0.0 A)")
    rounds(lemma, "g1", "(- B A)")
    lemma.fact("(< g0 (- 0.0 clear))")
    lemma.fact("(> g1 clear)")
    lemma.fact(f"(<= (absr (- k A)) (+ (* {CHAIN} u (absr A)) (* {power(112)} h)))")
    lemma.goal(f"(>= A (* {power(-46)} B))")
    lemma.goal(f"(>= (- B A) (* {power(-46)} B))")
    lemma.goal("(> k 0.0)")
    lemma.goal("(< k dx)")
    return lemma


def steep_lines():
    """Intersect's shallow test on a level line and an upright one, at a tolerance's angle within
    Steep / 2: |rnd(d1 d2)| against rnd(rnd(2 sin(angle) l1) l2), lengths l1 and l2 as Length
    gives them, the sine under 0.0626. The lines are never shallow."""
    lemma = Lemma("steep lines", "d1", "d2", "l1", "l2", "cross", "sine", "t2", "t3")
    lemma.fact(f"(>= (absr d1) {power(-108)})")
    lemma.fact(f"(>= (absr d2) {power(-108)})")
    lemma.fact("(<= (absr (- l1 (absr d1))) (* 1.51 u (absr d1)))")
    lemma.fact("(<= (absr (- l2 (absr d2))) (* 1.51 u (absr d2)))")
    rounds(lemma, "cross", "(* d1 d2)")
    lemma.fact("(>= sine 0.0)")
    lemma.fact("(<= sine 0.0626)")
    rounds(lemma, "t2", "(* 2.0 sine l1)")
    rounds(lemma, "t3", "(* t2 l2)")
    lemma.goal("(> (absr cross) t3)")
    return lemma


def length_sum():
    """Length of (a, b) where its largest part is from 2^-500 to 2^500: rnd(rnd(a^2) + rnd(b^2)),
    the smaller square allowed to underflow, is within 2.0001 u of a^2 + b^2."""
    lemma = Lemma("length sum", "a2", "b2", "sa", "sb", "s")
    lemma.fact("(>= a2 0.0)")
    lemma.fact("(>= b2 0.0)")
    lemma.fact(f"(>= (+ a2 b2) {power(-1000)})")
    multiplies(lemma, "sa", "a2")
    multiplies(lemma, "sb", "b2")
    rounds(lemma, "s", "(+ sa sb)")
    lemma.goal("(<= (absr (- s (+ a2 b2))) (* 2.0001 u (+ a2 b2)))")
    return lemma


def length_main():
    """And its root rounded: within 2.01 u of the exact length T."""
    lemma = Lemma("length main", "T2", "s", "w", "length", "T")
    lemma.fact("(>= T2 0.0)")
    lemma.fact("(<= (absr (- s T2)) (* 2.0001 u T2))")
    rooted(lemma, "length", "w", "s")
    lemma.fact("(>= T 0.0)")
    lemma.fact("(= (* T T) T2)")
    lemma.goal("(<= (absr (- length T)) (* 2.01 u T))")
    return lemma


def scaled_square():
    """Length of (a, b) where its largest part m = |a| is under 2^-500: a / m is 1 exactly, and
    beta = b / m rounds to s, underflow allowed; s^2 is within (2 u + u^2) beta^2 + 3 h of
    beta^2."""
    lemma = Lemma("scaled square", "beta", "s")
    lemma.fact("(<= (absr beta) 1.0)")
    multiplies(lemma, "s", "beta")
    lemma.goal("(<= (absr (- (* s s) (* beta beta))) "
               "(+ (* (+ (* 2.0 u) (* u u)) (* beta beta)) (* 3.0 h)))")
    return lemma


def scaled_sum():
    """Then the sum of squares, dot = rnd(1 + rnd(s^2)), B2 = beta^2 and S2 = s^2: within 2.51 u
    T2 + 6 h of T2 = 1 + B2."""
    lemma = Lemma("scaled sum", "B2", "S2", "q", "dot")
    lemma.fact("(>= B2 0.0)")
    lemma.fact("(<= B2 1.0)")
    lemma.fact("(>= S2 0.0)")
    lemma.fact("(<= (absr (- S2 B2)) (+ (* (+ (* 2.0 u) (* u u)) B2) (* 3.0 h)))")
    lemma.fact("(>= q 0.0)")
    multiplies(lemma, "q", "S2")
    rounds(lemma, "dot", "(+ 1.0 q)")
    lemma.goal("(<= (absr (- dot (+ 1.0 B2))) (+ (* 2.51 u (+ 1.0 B2)) (* 6.0 h)))")
    return lemma


def scaled_root():
    """Then its root rounded, r1: within 2.27 u tau + 7 h of tau, the root of T2, at least 1."""
    lemma = Lemma("scaled root", "T2", "dot", "tau", "w", "r1")
    lemma.fact("(>= T2 1.0)")
    lemma.fact("(<= (absr (- dot T2)) (+ (* 2.51 u T2) (* 6.0 h)))")
    lemma.fact("(>= tau 0.0)")
    lemma.fact("(= (* tau tau) T2)")
    rooted(lemma, "r1", "w", "dot")
    lemma.goal("(<= (absr (- r1 tau)) (+ (* 2.27 u tau) (* 7.0 h)))")
    return lemma


def length_tiny():
    """And length = rnd(m r1), underflow allowed: within 2^-540 of the exact length m tau, m being
    under 2^-500."""
    lemma = Lemma("length tiny", "m", "tau", "r1", "length")
    lemma.fact("(> m 0.0)")
    lemma.fact(f"(< m {power(-500)})")
    lemma.fact("(>= tau 1.0)")
    lemma.fact("(<= tau 1.5)")
    lemma.fact("(<= (absr (- r1 tau)) (+ (* 2.27 u tau) (* 7.0 h)))")
    multiplies(lemma, "length", "(* m r1)")
    lemma.goal(f"(<= (absr (- length (* m tau))) {power(-540)})")
    return lemma


def hypotenuse():
    """The exact length T of (a, b) lies between |b| and |a| + |b|."""
    lemma = Lemma("hypotenuse", "a", "b", "T")
    lemma.fact("(>= T 0.0)")
    lemma.fact("(= (* T T) (+ (* a a) (* b b)))")
    lemma.goal("(<= (absr b) T)")
    lemma.goal("(<= T (+ (absr a) (absr b)))")
    return lemma


def distance(bound=OFF):
    """Roundly's off: Length of (rnd(p - fx), rnd(q - c)), the foot fx within 11.01 u S of p, c
    the carrier's coordinate: within 17.1 u S of |g|, g = c - q exact."""
    lemma = Lemma("distance", "S", "p", "fx", "c", "q", "a", "b", "T", "off")
    lemma.fact(f"(>= S {power(-100)})")
    lemma.fact(f"(<= (absr (- fx p)) (* {FOOT} u S))")
    lemma.fact("(<= (absr c) S)")
    lemma.fact("(<= (absr q) (* S (+ 1.0 (* 2.0 u))))")
    rounds(lemma, "a", "(- p fx)")
    rounds(lemma, "b", "(- q c)")
    lemma.fact("(<= (absr b) T)")
    lemma.fact("(<= T (+ (absr a) (absr b)))")
    measured(lemma, "off", "T")
    lemma.goal(f"(<= (absr (- off (absr (- q c)))) (* {bound} u S))")
    return lemma


def rest_test(grain=GRAIN):
    """Spotted's test, rnd(rnd(rest_H Tight) - offRest) > rnd(r_H^2 Loose) / 32, offRest =
    rnd(Grain rnd(r_H^2 + rnd(g_H^2))): r^2 - g^2, which rest holds within 2^-101 (r_H^2 +
    g_H^2) by lean Times and Plus, is past r^2 / 32. r^2 = r_H^2 + r_L exactly, TwoProduct's."""
    lemma = Lemma("rest test", "sqH", "sqL", "R2", "G2", "GH2", "restH", "restL", "gg", "s1",
                  "offRest", "t1", "t2", "t3", "t4")
    lemma.fact(f"(>= sqH {power(-200)})")
    lemma.fact("(= R2 (+ sqH sqL))")
    lemma.fact("(<= (absr sqL) (* u sqH))")
    lemma.fact("(>= G2 0.0)")
    lemma.fact("(>= GH2 0.0)")
    lemma.fact(f"(<= (absr (- (+ restH restL) (- R2 G2))) (* {power(-101)} (+ sqH GH2)))")
    lemma.fact("(<= (absr restL) (* u (absr restH)))")
    multiplies(lemma, "gg", "GH2")
    rounds(lemma, "s1", "(+ sqH gg)")
    rounds(lemma, "offRest", f"(* {rational(grain)} s1)")
    rounds(lemma, "t1", f"(* restH {rational(TIGHT)})")
    rounds(lemma, "t2", "(- t1 offRest)")
    rounds(lemma, "t3", f"(* sqH {rational(LOOSE)})")
    lemma.fact("(= (* t4 32.0) t3)")
    lemma.fact("(> t2 t4)")
    lemma.goal("(> (- R2 G2) (/ R2 32.0))")
    return lemma


def cut_deep():
    """Roundly's tests on the carrier, r^2 - G^2 past r^2 / 32, off within 17.1 u S of G = |g|, r
    past 1024 L and 2^-32 S: r past 16 clear; |rnd(off - r)| past clear = 4 L, so no touch;
    rnd(1 - rnd(rnd(off / r)^2)) past Steep^2 = 1/64; and off under r, so two crossings."""
    lemma = Lemma("cut deep", "r", "G", "off", "L", "S", "gap", "ratio", "square", "one")
    lemma.fact("(> r 0.0)")
    lemma.fact("(>= G 0.0)")
    lemma.fact("(>= off 0.0)")
    lemma.fact("(> L 0.0)")
    lemma.fact("(> (* 32.0 (- (* r r) (* G G))) (* r r))")
    lemma.fact(f"(<= (absr (- off G)) (* {OFF} u S))")
    lemma.fact("(> r (* 1024.0 L))")
    lemma.fact(f"(>= L (* {power(-42)} S))")
    rounds(lemma, "gap", "(- off r)")
    divides(lemma, "ratio", "off", "r")
    lemma.fact("(>= square 0.0)")
    multiplies(lemma, "square", "(* ratio ratio)")
    rounds(lemma, "one", "(- 1.0 square)")
    lemma.goal("(> r (* 64.0 L))")
    lemma.goal("(> (absr gap) (* 4.0 L))")
    lemma.goal("(> one (/ 1.0 64.0))")
    lemma.goal("(< off r)")
    return lemma


def squares():
    """The squares Roundly's radicand is made of, r^2 - G^2 past r^2 / 32, off within 17.1 u S of
    G = |g|, r past 2^-32 S: off^2 within 33.69 u S r of G^2, at most r^2; and r^2 at most (1 +
    2 u) S r."""
    lemma = Lemma("squares", "r", "G", "off", "S")
    lemma.fact(f"(>= r {power(-100)})")
    lemma.fact("(>= G 0.0)")
    lemma.fact("(>= off 0.0)")
    lemma.fact("(> (* 32.0 (- (* r r) (* G G))) (* r r))")
    lemma.fact(f"(<= (absr (- off G)) (* {OFF} u S))")
    lemma.fact("(<= r (* S (+ 1.0 (* 2.0 u))))")
    lemma.fact(f"(> r (* {power(-32)} S))")
    lemma.goal("(<= (absr (- (* off off) (* G G))) (* 33.69 u S r))")
    lemma.goal("(<= (* off off) (* r r))")
    lemma.goal("(<= (* r r) (* (+ 1.0 (* 2.0 u)) S r))")
    return lemma


def radicand(bound="37.3"):
    """Roundly's radicand rnd(rnd(r^2) - rnd(off^2)), R2 = r^2, O2 = off^2, G2 = G^2 and X = S r
    as squares binds them: within 37.3 u S r of r^2 - G^2."""
    lemma = Lemma("radicand", "R2", "O2", "G2", "X", "rr", "oo", "rest")
    lemma.fact(f"(>= R2 {power(-200)})")
    lemma.fact("(>= O2 0.0)")
    lemma.fact("(<= (absr (- O2 G2)) (* 33.69 u X))")
    lemma.fact("(<= O2 R2)")
    lemma.fact("(<= R2 (* (+ 1.0 (* 2.0 u)) X))")
    rounds(lemma, "rr", "R2")
    multiplies(lemma, "oo", "O2")
    rounds(lemma, "rest", "(- rr oo)")
    lemma.goal(f"(<= (absr (- rest (- R2 G2))) (* {bound} u X))")
    return lemma


def root(bound=ROOT):
    """Roundly's root of the radicand rounded: within 106.6 u S of the true root rho, rho^2 = r^2 -
    g^2 past r^2 / 32, the radicand within 37.3 u S r of it, r past 2^-32 S."""
    lemma = Lemma("root", "r", "S", "rest", "rho", "w", "root")
    lemma.fact("(> r 0.0)")
    lemma.fact("(<= r (* S (+ 1.0 (* 2.0 u))))")
    lemma.fact(f"(> r (* {power(-32)} S))")
    lemma.fact("(>= rho 0.0)")
    lemma.fact("(> (* 32.0 rho rho) (* r r))")
    lemma.fact("(<= (* rho rho) (* r r))")
    lemma.fact("(<= (absr (- rest (* rho rho))) (* 37.3 u S r))")
    rooted(lemma, "root", "w", "rest")
    lemma.goal(f"(<= (absr (- root rho)) (* {bound} u S))")
    return lemma


def half():
    """Roundly's half, rnd(root / length), length the root of the line's square as level length
    gives it: taken along the way, kh = |dx| half, it is within 2.52 u of the root."""
    lemma = Lemma("half", "dx", "length", "half", "root", "kh")
    lemma.fact(f"(>= (absr dx) {power(-108)})")
    lemma.fact(f"(<= (absr dx) {power(93)})")
    lemma.fact("(<= (absr (- length (absr dx))) (* 1.51 u (absr dx)))")
    lemma.fact("(>= root 0.0)")
    divides(lemma, "half", "root", "length")
    lemma.fact("(= kh (* (absr dx) half))")
    lemma.goal(f"(<= (absr (- kh root)) (+ (* 2.52 u root) (* {power(94)} h)))")
    return lemma


def scaled():
    """A rounding taken along dx stays a rounding: |s - z| <= u |z| gives |dx s - dx z| <= u |dx
    z|."""
    lemma = Lemma("scaled", "s", "z", "dx", "sigma", "zeta")
    lemma.fact("(<= (absr (- s z)) (* u (absr z)))")
    lemma.fact("(= sigma (* dx s))")
    lemma.fact("(= zeta (* dx z))")
    lemma.goal("(<= (absr (- sigma zeta)) (* u (absr zeta)))")
    return lemma


def crossing(bound="121.0"):
    """Roundly's crossing a0 + rnd(dx rnd(middle -+ half)), its place first: taken along dx,
    middle is km within 4.0002 u |A| of A, the gap from a0 to p, by the quotient chain, and half
    is K within 2.52 u of the root, itself within 106.6 u S of the true one P, signed; their sum
    rounds (scaled) to sigma, the place taken along dx: within 121 u S of A + P."""
    lemma = Lemma("crossing", "S", "a0", "A", "P", "Pc", "km", "K", "sigma")
    lemma.fact(f"(>= S {power(-100)})")
    lemma.fact("(<= (absr a0) S)")
    lemma.fact("(<= (absr (+ a0 A)) (* S (+ 1.0 (* 2.0 u))))")
    lemma.fact("(<= (absr P) (* S (+ 1.0 (* 2.0 u))))")
    lemma.fact(f"(<= (absr (- km A)) (+ (* {CHAIN} u (absr A)) (* {power(112)} h)))")
    lemma.fact(f"(<= (absr (- Pc P)) (* {ROOT} u S))")
    lemma.fact(f"(<= (absr (- K Pc)) (+ (* 2.52 u (absr Pc)) (* {power(94)} h)))")
    lemma.fact("(<= (absr (- sigma (+ km K))) (* u (absr (+ km K))))")
    lemma.goal(f"(<= (absr (- sigma (+ A P))) (* {bound} u S))")
    return lemma


def crossing_point(bound=CROSSING):
    """Then the crossing a0 + rnd(sigma), rounded: within 128 u S of the true one a0 + A + P,
    which lies within S (1 + 2 u); and so is its place along the line, taken along dx."""
    lemma = Lemma("crossing point", "S", "a0", "A", "P", "sigma", "m", "out")
    lemma.fact(f"(>= S {power(-100)})")
    lemma.fact("(<= (absr a0) S)")
    lemma.fact("(<= (absr (+ a0 A P)) (* S (+ 1.0 (* 2.0 u))))")
    lemma.fact("(<= (absr (- sigma (+ A P))) (* 121.0 u S))")
    multiplies(lemma, "m", "sigma")
    rounds(lemma, "out", "(+ a0 m)")
    lemma.goal(f"(<= (absr (- out (+ a0 A P))) (* {bound} u S))")
    return lemma


def spotted_rest():
    """Spotted's radicand rest, against r^2 - g^2 = rad: their roots rhoR and rho differ by at most
    |rest - rad| / rhoR."""
    lemma = Lemma("spotted rest", "rest", "rad", "rhoR", "rho", "E")
    lemma.fact("(> rhoR 0.0)")
    lemma.fact("(>= rho 0.0)")
    lemma.fact("(= (* rhoR rhoR) rest)")
    lemma.fact("(= (* rho rho) rad)")
    lemma.fact("(<= (absr (- rest rad)) E)")
    lemma.goal("(<= (* (absr (- rhoR rho)) rhoR) E)")
    return lemma


def gap_high():
    """Spotted's gap g = g_H + g_L exact, |g_L| <= u |g_H|, g^2 under 31/32 of r^2 <= r_H^2 (1 +
    u) by rest test: g_H^2 is at most r_H^2."""
    lemma = Lemma("gap high", "gH", "gL", "sqH")
    lemma.fact("(> sqH 0.0)")
    lemma.fact("(<= (absr gL) (* u (absr gH)))")
    lemma.fact("(< (* 32.0 (* (+ gH gL) (+ gH gL))) (* 31.0 sqH (+ 1.0 u)))")
    lemma.goal("(<= (* gH gH) sqH)")
    return lemma


def spotted_bounds():
    """Spotted's offRest = rnd(Grain rnd(r_H^2 + rnd(g_H^2))) is at least 7.9 times E = 2^-101
    (r_H^2 + g_H^2), how far lean double-double leaves its radicand rest from rad = r^2 - g^2,
    and at most 2^-97 (1 + 4 u) r_H^2; rest is past r_H^2 / 33 and at most r_H^2 (1 + 2 u); and
    lean Root's root rhoR, within 2^-103 of its high part s and its low part within 1.52 u s, is
    within s 2^-52 of s."""
    lemma = Lemma("spotted bounds", "sqH", "GH2", "gg", "s1", "offRest", "rad", "rest", "rootH",
                  "rootL", "rhoR")
    lemma.fact(f"(>= sqH {power(-200)})")
    lemma.fact("(>= GH2 0.0)")
    lemma.fact("(<= GH2 sqH)")
    multiplies(lemma, "gg", "GH2")
    rounds(lemma, "s1", "(+ sqH gg)")
    rounds(lemma, "offRest", f"(* {rational(GRAIN)} s1)")
    lemma.fact("(> (* 32.0 rad) (* sqH (- 1.0 u)))")
    lemma.fact("(<= rad (* sqH (+ 1.0 u)))")
    lemma.fact(f"(<= (absr (- rest rad)) (* {power(-101)} (+ sqH GH2)))")
    lemma.fact("(> rootH 0.0)")
    lemma.fact(f"(<= (absr (- (+ rootH rootL) rhoR)) (* {power(-103)} rootH))")
    lemma.fact("(<= (absr rootL) (* 1.52 u rootH))")
    lemma.goal(f"(<= (* 7.9 {power(-101)} (+ sqH GH2)) offRest)")
    lemma.goal(f"(<= offRest (* {power(-97)} (+ 1.0 (* 4.0 u)) sqH))")
    lemma.goal("(> (* 33.0 rest) sqH)")
    lemma.goal("(<= rest (* sqH (+ 1.0 (* 2.0 u))))")
    lemma.goal(f"(<= (absr (- rhoR rootH)) (* rootH {power(-52)}))")
    return lemma


def spotted_scale():
    """The root's error e, e rhoR at most E, E at most offRest / 7.9, by spotted rest and spotted
    bounds: e is at most 0.13 offRest inv, inv = rnd(1 / s)."""
    lemma = Lemma("spotted scale", "e", "E", "rhoR", "rootH", "inv", "offRest", "oi")
    lemma.fact("(>= e 0.0)")
    lemma.fact("(>= offRest 0.0)")
    lemma.fact("(> rootH 0.0)")
    lemma.fact("(<= (* e rhoR) E)")
    lemma.fact("(<= (* 7.9 E) offRest)")
    lemma.fact(f"(<= (absr (- rhoR rootH)) (* rootH {power(-52)}))")
    lemma.fact("(<= (absr (- (* inv rootH) 1.0)) u)")
    lemma.fact("(= oi (* offRest inv))")
    lemma.goal("(<= e (* 0.13 oi))")
    return lemma


def spotted_root():
    """Spotted's root against the true one: within 2^-103 s of rhoR, by lean Root, and rhoR
    within e of the true root, at most 0.13 offRest inv by spotted scale: within a fifth of
    offRoot = rnd(rnd(rnd(offRest inv) Loose) + rnd(Grain s))."""
    lemma = Lemma("spotted root", "rootH", "rootL", "rhoR", "rho", "oi", "q1", "q2", "q3",
                  "offRoot")
    lemma.fact("(> rootH 0.0)")
    lemma.fact("(>= oi 0.0)")
    lemma.fact(f"(<= (absr (- (+ rootH rootL) rhoR)) (* {power(-103)} rootH))")
    lemma.fact("(<= (absr (- rhoR rho)) (* 0.13 oi))")
    rounds(lemma, "q1", "oi")
    rounds(lemma, "q2", f"(* q1 {rational(LOOSE)})")
    rounds(lemma, "q3", f"(* {rational(GRAIN)} rootH)")
    rounds(lemma, "offRoot", "(+ q2 q3)")
    lemma.goal("(<= (absr (- (+ rootH rootL) rho)) (* 0.2 offRoot))")
    return lemma


def spotted_share():
    """offRest inv, inv = rnd(1 / s), times rhoR, within s 2^-52 of s: at most 2^-97 1.0001 r_H^2,
    offRest being at most 2^-97 (1 + 4 u) r_H^2."""
    lemma = Lemma("spotted share", "sqH", "offRest", "inv", "rootH", "rhoR", "oi")
    lemma.fact("(> sqH 0.0)")
    lemma.fact("(>= offRest 0.0)")
    lemma.fact(f"(<= offRest (* {power(-97)} (+ 1.0 (* 4.0 u)) sqH))")
    lemma.fact("(> rootH 0.0)")
    lemma.fact("(<= (absr (- (* inv rootH) 1.0)) u)")
    lemma.fact(f"(<= (absr (- rhoR rootH)) (* rootH {power(-52)}))")
    lemma.fact("(= oi (* offRest inv))")
    lemma.goal(f"(<= (* oi rhoR) (* {power(-97)} 1.0001 sqH))")
    return lemma


def spotted_reach():
    """offRest inv is under 2^-94 S, and rhoR at most S (1 + 5 u): rhoR^2 = rest is past r_H^2 /
    33 and at most r_H^2 (1 + 2 u), r_H^2 at most S^2 (1 + 6 u)."""
    lemma = Lemma("spotted reach", "S", "sqH", "rhoR", "oi")
    lemma.fact("(> sqH 0.0)")
    lemma.fact("(>= S 0.0)")
    lemma.fact("(<= sqH (* S S (+ 1.0 (* 6.0 u))))")
    lemma.fact("(> rhoR 0.0)")
    lemma.fact("(> (* 33.0 rhoR rhoR) sqH)")
    lemma.fact("(<= (* rhoR rhoR) (* sqH (+ 1.0 (* 2.0 u))))")
    lemma.fact("(>= oi 0.0)")
    lemma.fact(f"(<= (* oi rhoR) (* {power(-97)} 1.0001 sqH))")
    lemma.goal(f"(<= oi (* {power(-94)} S))")
    lemma.goal("(<= rhoR (* S (+ 1.0 (* 5.0 u))))")
    return lemma


def spotted_size():
    """offRoot = rnd(rnd(rnd(offRest inv) Loose) + rnd(Grain s)) is under 2^-93 S."""
    lemma = Lemma("spotted size", "S", "rhoR", "rootH", "oi", "q1", "q2", "q3", "offRoot")
    lemma.fact("(>= S 0.0)")
    lemma.fact("(> rootH 0.0)")
    lemma.fact("(>= oi 0.0)")
    lemma.fact(f"(<= oi (* {power(-94)} S))")
    lemma.fact("(<= rhoR (* S (+ 1.0 (* 5.0 u))))")
    lemma.fact(f"(<= (absr (- rhoR rootH)) (* rootH {power(-52)}))")
    rounds(lemma, "q1", "oi")
    rounds(lemma, "q2", f"(* q1 {rational(LOOSE)})")
    rounds(lemma, "q3", f"(* {rational(GRAIN)} rootH)")
    rounds(lemma, "offRoot", "(+ q2 q3)")
    lemma.goal(f"(<= offRoot (* {power(-93)} S))")
    return lemma


def spotted_part():
    """Spotted's crossing, the center's coordinate c plus the signed root Rh + Rl, |Rl| <= 1.52 u
    |Rh|: TwoSum's sumH + sumL exact and part = rnd(sumL + Rl), each within W = |c| + |Rh|."""
    lemma = Lemma("spotted part", "c", "Rh", "Rl", "sumH", "sumL", "part")
    lemma.fact("(<= (absr Rl) (* 1.52 u (absr Rh)))")
    lemma.fact("(= (+ sumH sumL) (+ c Rh))")
    rounds(lemma, "sumH", "(+ c Rh)")
    rounds(lemma, "part", "(+ sumL Rl)")
    size = "(+ (absr c) (absr Rh))"
    lemma.goal(f"(<= (absr part) (* 2.53 u {size}))")
    lemma.goal(f"(<= (absr (- part (+ sumL Rl))) (* 2.6 u u {size}))")
    lemma.goal(f"(<= (absr sumH) (* (+ 1.0 u) {size}))")
    return lemma


def spotted_value():
    """Then value = rnd(sumH + part) and the rest Certain is given, rnd(rnd(sumH - value) +
    part): their sum within 5 u^2 W of sumH + part."""
    lemma = Lemma("spotted value", "W", "sumH", "part", "value", "diff", "rest")
    lemma.fact("(>= W 0.0)")
    lemma.fact("(<= (absr part) (* 2.53 u W))")
    lemma.fact("(<= (absr sumH) (* (+ 1.0 u) W))")
    rounds(lemma, "value", "(+ sumH part)")
    rounds(lemma, "diff", "(- sumH value)")
    rounds(lemma, "rest", "(+ diff part)")
    lemma.goal("(<= (absr (- (+ value rest) (+ sumH part))) (* 5.0 u u W))")
    return lemma


def spotted_sum(doubled=True):
    """Then the crossing is within 0.6 of bound = rnd(2 offRoot + rnd(Grain rnd(|c| + |Rh|))) of
    the true one X = c + rho, as Certain asks, the root within a fifth of offRoot by spotted
    root; and the bound is past 2^-99 |value|."""
    lemma = Lemma("spotted sum", "c", "Rh", "Rl", "rho", "offRoot", "sumH", "sumL", "part",
                  "value", "rest", "size", "g", "bound")
    lemma.fact("(>= offRoot 0.0)")
    lemma.fact("(<= (absr (- (+ Rh Rl) rho)) (* 0.2 offRoot))")
    lemma.fact("(= (+ sumH sumL) (+ c Rh))")
    size = "(+ (absr c) (absr Rh))"
    lemma.fact(f"(<= (absr (- part (+ sumL Rl))) (* 2.6 u u {size}))")
    lemma.fact(f"(<= (absr part) (* 2.53 u {size}))")
    lemma.fact(f"(<= (absr sumH) (* (+ 1.0 u) {size}))")
    lemma.fact(f"(<= (absr (- (+ value rest) (+ sumH part))) (* 5.0 u u {size}))")
    rounds(lemma, "value", "(+ sumH part)")
    rounds(lemma, "size", size)
    rounds(lemma, "g", f"(* {rational(GRAIN)} size)")
    rounds(lemma, "bound", "(+ (* 2.0 offRoot) g)" if doubled else "g")
    lemma.goal("(<= (absr (- (+ value rest) (+ c rho))) (* 0.6 bound))")
    lemma.goal(f"(>= bound (* {power(-99)} (absr value)))")
    return lemma


def spotted_near():
    """And value is within 1.01 u S of X, certified or not, X within S (1 + 2 u) and offRoot
    under 2^-93 S."""
    lemma = Lemma("spotted near", "S", "c", "rho", "Rh", "Rl", "offRoot", "y", "e", "value")
    lemma.fact(f"(>= S {power(-100)})")
    lemma.fact("(<= (absr (+ c rho)) (* S (+ 1.0 (* 2.0 u))))")
    lemma.fact("(<= (absr c) (* S (+ 1.0 (* 2.0 u))))")
    lemma.fact("(<= (absr Rh) (* S (+ 1.0 (* 2.0 u))))")
    lemma.fact("(>= offRoot 0.0)")
    lemma.fact(f"(<= offRoot (* {power(-93)} S))")
    lemma.fact("(<= (absr (- (+ Rh Rl) rho)) (* 0.2 offRoot))")
    lemma.fact("(<= (absr (- y (+ c Rh Rl))) (* 2.6 u u (+ (absr c) (absr Rh))))")
    rounds(lemma, "value", "y")
    lemma.goal("(<= (absr (- value (+ c rho))) (* 1.01 u S))")
    return lemma


def told_near():
    """Told's near test on a distance squared, rnd(rnd(a^2) + rnd(b^2)) <= rnd(rnd(near^2) Tight),
    underflow allowed, near past 2^-300: the exact distance d of (a, b) is at most near."""
    lemma = Lemma("told near", "a2", "b2", "sa", "sb", "s", "near", "nn", "limit", "d")
    lemma.fact("(>= a2 0.0)")
    lemma.fact("(>= b2 0.0)")
    multiplies(lemma, "sa", "a2")
    multiplies(lemma, "sb", "b2")
    rounds(lemma, "s", "(+ sa sb)")
    lemma.fact(f"(>= near {power(-300)})")
    rounds(lemma, "nn", "(* near near)")
    rounds(lemma, "limit", f"(* nn {rational(TIGHT)})")
    lemma.fact("(<= s limit)")
    lemma.fact("(>= d 0.0)")
    lemma.fact("(= (* d d) (+ a2 b2))")
    lemma.goal("(<= d near)")
    return lemma


def told_far():
    """Told's far test, rnd(rnd(a^2) + rnd(b^2)) > rnd(rnd(far^2) Loose): the exact distance d is
    past far."""
    lemma = Lemma("told far", "a2", "b2", "sa", "sb", "s", "far", "ff", "limit", "d")
    lemma.fact("(>= a2 0.0)")
    lemma.fact("(>= b2 0.0)")
    multiplies(lemma, "sa", "a2")
    multiplies(lemma, "sb", "b2")
    rounds(lemma, "s", "(+ sa sb)")
    lemma.fact(f"(>= far {power(-300)})")
    rounds(lemma, "ff", "(* far far)")
    rounds(lemma, "limit", f"(* ff {rational(LOOSE)})")
    lemma.fact("(> s limit)")
    lemma.fact("(>= d 0.0)")
    lemma.fact("(= (* d d) (+ a2 b2))")
    lemma.goal("(> d far)")
    return lemma


def triangle():
    """Two vectors (a, y) and (b, y) differ in length by |a - b| at most."""
    lemma = Lemma("triangle", "a", "b", "y", "d1", "d2")
    lemma.fact("(>= d1 0.0)")
    lemma.fact("(>= d2 0.0)")
    lemma.fact("(= (* d1 d1) (+ (* a a) (* y y)))")
    lemma.fact("(= (* d2 d2) (+ (* b b) (* y y)))")
    lemma.goal("(<= (absr (- d2 d1)) (absr (- a b)))")
    return lemma


def within(slack=SLACK):
    """An end's distance from a crossing against half the tolerance: here from the spot sx,
    within 1.01 u S of the true crossing X, its square told near (L / 2 - slack) or far (L / 2 +
    slack); Roundly's Distance from its crossing xR, within 128 u S of X, against L / 2. The two
    share the carrier's coordinate, and the end's gap across it; told near, Roundly's is at most
    L / 2, told far past it. L is at most S (1 + 2 u) / 1024."""
    lemma = Lemma("within", "S", "L", "near", "far", "sx", "xR", "X", "ex", "ax", "bx", "dA",
                  "dR", "length")
    lemma.fact(f"(>= S {power(-100)})")
    lemma.fact("(> L 0.0)")
    lemma.fact("(<= (* 1024.0 L) (* S (+ 1.0 (* 2.0 u))))")
    rounds(lemma, "near", f"(- (/ L 2.0) (* {rational(slack)} S))")
    rounds(lemma, "far", f"(+ (/ L 2.0) (* {rational(slack)} S))")
    lemma.fact("(<= (absr (- sx X)) (* 1.01 u S))")
    lemma.fact(f"(<= (absr (- xR X)) (* {CROSSING} u S))")
    rounds(lemma, "ax", "(- sx ex)")
    rounds(lemma, "bx", "(- xR ex)")
    lemma.fact("(>= dA (absr ax))")
    lemma.fact("(>= dR (absr bx))")
    lemma.fact("(<= (absr (- dR dA)) (absr (- ax bx)))")
    measured(lemma, "length", "dR")
    lemma.goal("(=> (<= dA near) (<= length (/ L 2.0)))")
    lemma.goal("(=> (> dA far) (> length (/ L 2.0)))")
    return lemma


def off_told():
    """A line's end against the circle: its exact distance d from the center, of the components
    Roundly and Told both round, told within r - 4 L - slack (deep) or past r + 4 L + slack
    (wide); Roundly's Off, |rnd(Length - r)|, is then past 4 L, the end far."""
    lemma = Lemma("off told", "S", "L", "r", "deep", "nd", "wide", "nw", "d", "length", "off")
    lemma.fact(f"(>= S {power(-100)})")
    lemma.fact("(> L 0.0)")
    lemma.fact("(> r (* 1024.0 L))")
    lemma.fact("(<= r (* S (+ 1.0 (* 2.0 u))))")
    rounds(lemma, "deep", "(- r (* 4.0 L))")
    rounds(lemma, "nd", f"(- deep (* {rational(SLACK)} S))")
    rounds(lemma, "wide", "(+ r (* 4.0 L))")
    rounds(lemma, "nw", f"(+ wide (* {rational(SLACK)} S))")
    lemma.fact("(>= d 0.0)")
    measured(lemma, "length", "d")
    rounds(lemma, "off", "(- length r)")
    lemma.goal("(=> (<= d nd) (> (absr off) (* 4.0 L)))")
    lemma.goal("(=> (> d nw) (> (absr off) (* 4.0 L)))")
    return lemma


def away_far(slack=SLACK):
    """An arc's end against the line: its exact distance dA from the line's true foot fA, told
    past 4 L + slack; Intersect's Away from its foot fR, within 21.1 u S of fA by along, the two
    sharing the gap across the carrier: past 4 L, the end far."""
    lemma = Lemma("away far", "S", "L", "far", "fA", "fR", "cx", "ax", "bx", "dA", "dR",
                  "length")
    lemma.fact(f"(>= S {power(-100)})")
    lemma.fact("(> L 0.0)")
    lemma.fact("(<= (* 1024.0 L) (* S (+ 1.0 (* 2.0 u))))")
    rounds(lemma, "far", f"(+ (* 4.0 L) (* {rational(slack)} S))")
    lemma.fact(f"(<= (absr (- fR fA)) (* {ACROSS} u S))")
    rounds(lemma, "ax", "(- fA cx)")
    rounds(lemma, "bx", "(- fR cx)")
    lemma.fact("(>= dA (absr ax))")
    lemma.fact("(>= dR (absr bx))")
    lemma.fact("(<= (absr (- dR dA)) (absr (- ax bx)))")
    measured(lemma, "length", "dR")
    lemma.goal("(=> (> dA far) (> length (* 4.0 L)))")
    return lemma


def inside(slack=SLACK):
    """A crossing's place on the line from a0 to a0 + B, dx = rnd(B) > 0 (a mirror takes the other
    sign): Inside on the spot sx, within 1.01 u S of the true crossing X = a0 + A, against Roundly's
    place, which taken along dx is within 128 u S of A. Told within by the slack, the place is
    within (0, 1); told off, outside it."""
    lemma = Lemma("inside", "S", "a0", "B", "dx", "A", "sx", "sigma", "i0", "i1", "o0", "o1")
    lemma.fact(f"(>= S {power(-100)})")
    lemma.fact("(> B 0.0)")
    lemma.fact("(> dx 0.0)")
    rounds(lemma, "dx", "B")
    lemma.fact("(<= (absr a0) S)")
    lemma.fact("(<= (absr (+ a0 B)) S)")
    lemma.fact("(<= (absr (- sx (+ a0 A))) (* 1.01 u S))")
    lemma.fact(f"(<= (absr (- sigma A)) (* {CROSSING} u S))")
    rounds(lemma, "i0", "(- sx a0)")
    rounds(lemma, "i1", "(- (+ a0 B) sx)")
    rounds(lemma, "o0", f"(- a0 (* {rational(slack)} S))")
    rounds(lemma, "o1", f"(+ (+ a0 B) (* {rational(slack)} S))")
    spare = f"(* {rational(slack)} S)"
    lemma.goal(f"(=> (and (> i0 {spare}) (> i1 {spare})) (and (> sigma 0.0) (< sigma dx)))")
    lemma.goal("(=> (or (< sx o0) (> sx o1)) (or (< sigma 0.0) (> sigma dx)))")
    return lemma


def sine_product():
    """|sine r| <= r for |sine| <= 1."""
    lemma = Lemma("sine product", "sine", "r", "Y")
    lemma.fact("(> r 0.0)")
    lemma.fact("(<= (absr sine) 1.0)")
    lemma.fact("(= Y (* sine r))")
    lemma.goal("(<= (absr Y) r)")
    return lemma


def cut():
    """An arc's end from its center, rnd(rnd(q + rnd(Y)) - q), Y = sine r at most r by sine
    product, |q| + r <= S (1 + 2 u): at most r + 4 u S + h."""
    lemma = Lemma("cut", "S", "q", "r", "Y", "y1", "cy", "start")
    lemma.fact(f"(>= r {power(-100)})")
    lemma.fact("(<= (+ (absr q) r) (* S (+ 1.0 (* 2.0 u))))")
    lemma.fact("(<= (absr Y) r)")
    multiplies(lemma, "y1", "Y")
    rounds(lemma, "cy", "(+ q y1)")
    rounds(lemma, "start", "(- cy q)")
    lemma.goal("(<= (absr start) (+ r (* 4.0 u S) h))")
    return lemma


def toward():
    """The two towards' parts along the carrier, tA = rnd(sx - p) and tR = rnd(xR - p), the true
    crossing X within r of the center's p, the spot within 1.01 u S of it and Roundly's crossing
    within 128 u S, r past 2^-32 S: they differ by 131.1 u S at most, each at most r (1 + 2^-13)."""
    lemma = Lemma("toward", "S", "r", "p", "X", "sx", "xR", "tA", "tR")
    lemma.fact(f"(>= r {power(-100)})")
    lemma.fact(f"(> r (* {power(-32)} S))")
    lemma.fact("(<= r (* S (+ 1.0 (* 2.0 u))))")
    lemma.fact("(<= (absr (- X p)) r)")
    lemma.fact("(<= (absr (- sx X)) (* 1.01 u S))")
    lemma.fact(f"(<= (absr (- xR X)) (* {CROSSING} u S))")
    rounds(lemma, "tA", "(- sx p)")
    rounds(lemma, "tR", "(- xR p)")
    lemma.goal("(<= (absr (- tA tR)) (* 131.1 u S))")
    lemma.goal(f"(<= (absr tA) (* r (+ 1.0 {power(-13)})))")
    lemma.goal(f"(<= (absr tR) (* r (+ 1.0 {power(-13)})))")
    return lemma


def products():
    """The terms that differ, start tA and start tR, the arc's end from its center at most r + 4
    u S + h by cut: within 131.2 u S r of each other, each at most 1.0002 r^2; and r^2 at most (1
    + 2 u) S r."""
    lemma = Lemma("products", "S", "r", "start", "tA", "tR", "P1", "P2")
    lemma.fact(f"(>= r {power(-100)})")
    lemma.fact(f"(> r (* {power(-32)} S))")
    lemma.fact("(<= r (* S (+ 1.0 (* 2.0 u))))")
    lemma.fact("(<= (absr start) (+ r (* 4.0 u S) h))")
    lemma.fact("(<= (absr (- tA tR)) (* 131.1 u S))")
    lemma.fact(f"(<= (absr tA) (* r (+ 1.0 {power(-13)})))")
    lemma.fact(f"(<= (absr tR) (* r (+ 1.0 {power(-13)})))")
    lemma.fact("(= P1 (* start tA))")
    lemma.fact("(= P2 (* start tR))")
    lemma.goal("(<= (absr (- P1 P2)) (* 131.2 u S r))")
    lemma.goal("(<= (absr P1) (* 1.0002 r r))")
    lemma.goal("(<= (absr P2) (* 1.0002 r r))")
    lemma.goal("(<= (* r r) (* (+ 1.0 (* 2.0 u)) S r))")
    return lemma


def swept(slack=SLACK):
    """A crossing's turn from an arc's end: Cross(start, toward), toward from the center to the
    spot here and to Roundly's crossing there. The term across the carrier, Z = rnd(start_x
    toward_y), is the same double in both; the other, rnd(start_y toward_x), differs by the
    crossings' gap, P1 and P2 exact as products binds them, R2 = r^2 and X = S r. Told past the
    edge, rnd(slack S r), here, Roundly's turn has the same sign: so Swept tells within and off as
    Intersect's does, at a sweep up to a half turn or past it."""
    lemma = Lemma("swept", "X", "R2", "P1", "P2", "Z", "p1", "p2", "CA", "CR", "edge")
    lemma.fact(f"(>= R2 {power(-200)})")
    lemma.fact("(<= R2 (* (+ 1.0 (* 2.0 u)) X))")
    lemma.fact("(<= (absr (- P1 P2)) (* 131.2 u X))")
    lemma.fact("(<= (absr P1) (* 1.0002 R2))")
    lemma.fact("(<= (absr P2) (* 1.0002 R2))")
    multiplies(lemma, "p1", "P1")
    multiplies(lemma, "p2", "P2")
    rounds(lemma, "CA", "(- Z p1)")
    rounds(lemma, "CR", "(- Z p2)")
    rounds(lemma, "edge", f"(* {rational(slack)} X)")
    lemma.goal("(=> (> CA edge) (> CR 0.0))")
    lemma.goal("(=> (< CA (- 0.0 edge)) (< CR 0.0))")
    return lemma


MUTANTS = (("beside without the margin", lambda: beside(margin=1.0)),
           ("away without the margin", lambda: away(margin=1.0)),
           ("moved at 10.9 u S", lambda: moved(bound="10.9")),
           ("distance at 16.9 u S", lambda: distance(bound="16.9")),
           ("rest test at Grain 2^-104", lambda: rest_test(grain=2.0 ** -104)),
           ("radicand at 35 u S r", lambda: radicand(bound="35.0")),
           ("root at 105 u S", lambda: root(bound="105.0")),
           ("crossing at 118 u S", lambda: crossing(bound="118.0")),
           ("crossing point at 123 u S", lambda: crossing_point(bound="123.0")),
           ("spotted sum without 2 offRoot", lambda: spotted_sum(doubled=False)),
           ("within at slack 2^-47 S", lambda: within(slack=2.0 ** -47)),
           ("away far at slack 2^-50 S", lambda: away_far(slack=2.0 ** -50)),
           ("inside at slack 2^-47 S", lambda: inside(slack=2.0 ** -47)),
           ("swept at slack 2^-47 S", lambda: swept(slack=2.0 ** -47)))

LEMMAS = (level_length, beside, away, quotient_terms, quotient_rounding, quotient_chain, moved,
          near_end, along, placed, steep_lines, length_sum, length_main, scaled_square,
          scaled_sum, scaled_root, length_tiny, hypotenuse, distance, rest_test, cut_deep,
          squares, radicand, root, half, scaled, crossing, crossing_point, gap_high,
          spotted_rest, spotted_bounds, spotted_scale, spotted_root, spotted_share,
          spotted_reach, spotted_size, spotted_part, spotted_value, spotted_sum, spotted_near,
          told_near, told_far, triangle, within, off_told, away_far, inside, sine_product, cut,
          toward, products, swept)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, mutants=MUTANTS))
