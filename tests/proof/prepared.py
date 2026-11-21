#!/usr/bin/env python3
"""Machine-checks the bounds source/prepared.cpp's line pass certifies with: that Judged tells two
lines cross, or plainly do not, only where Straight tells them so, and that Crossed's bound, given
Certain, is under 0.6 of what Certain is given, so the point it certifies is the double nearest the
true crossing; and that Ends, the chord form's last sum, keeps within its bound.

    python3 -I tests/proof/prepared.py [z3]

A rounding is lemma.py's: |r - z| <= u |r| and |r - z| <= u |z| for a sum, a difference, a
product or a quotient clear of underflow, a quotient's written as |r d - n| <= u |n|. Rounding is
monotone (rounding.v's monotone), and 1 is a double: a quotient below 1 rounds to at most 1. The
lean toolkit's bounds are lean.py's, taken as facts here by name.

Judged against Straight (the same doubles: the ways, lengths, cross products and the gap's):
- plain: lap under rnd(rnd(L / 2 Tight) size) makes Straight's rnd(lap / size) under L / 2, and
  lap over rnd(rnd(32 L Loose) size) makes it over 32 L;
- placed: |n| reach over rnd(rnd(wide d) Loose), all rounded, makes Straight's place rnd(n / d)
  times reach over wide; |n - d| reach over that and 2^-50 |n| reach makes its rnd(place - 1)
  times reach over wide too;
- hit: then the place is above 0 exactly where n is, and below 1 exactly where n is below d.

Crossed: the lean cross products n' and d' within 2^-100 of their sizes Pn and Pd (lean.py's
cross); the quotient q' = n' / d' carries their errors to |t - q'| |d| <= En + |q'| Ed (carry);
the reciprocal of d's high part, I, is within 2^-50 of 1 / |d| where |d| >= Pd / 8 (inverse), as
the sine above Steep makes it; Over's q within 2^-101 |q| of q' (lean.py's over). Then the bound
Crossed works out, every rounding taken down, covers |t - q|: term by term, its Grain being four
times what each term needs, and the h's of underflow far below what is left, for the t within
[2^-45, 1] and the |d| above 2^-209 that Judged leaves (crossed_terms). Moved along the step, the
point is within 0.6 of the bound Certain is given, which is at least 2^-99 of the point
(crossed_certain).

Chorded's bound on the discriminant, OffD, covers its error, its parts' carried by sure.py's
carried and its products' and sum's own by lean.py's (chorded disc).

Left by hand, in prepared.cpp's comments: Lagrange's |wx ay| + |wy ax| <= |w| |a| with the plain
cross product's three roundings, which make |d| >= Pd / 8 (1 - 2^-48) of a sine above Steep (twin.py
proves Lagrange's identity); that a place clear of 0 and 1 by Wide tolerances is at least 2^-45 and
the ways at least 2^-103 long, from Reached's length; and the chord form's composite bounds, Cut's
and Circled's, whose steps lean.py proves.

Each lemma states the operations of prepared.cpp as written, one fact a rounding; they are not bound
to program.py's graphs. Each mutant drops a factor the code keeps, and Z3 must refute it. The exit
status is 1 when a lemma is not proved or a mutant is.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from expand import Fraction, U, verdict  # pylint: disable=wrong-import-position
# pylint: disable-next=wrong-import-position
from lemma import Lemma, rational, run  # pylint: disable=wrong-import-position

LOOSE = 1.0 + 2.0 ** -48
TIGHT = 1.0 - 2.0 ** -48
GRAIN = 2.0 ** -98
H = Fraction(1, 2 ** 1075)


def rounds(lemma, result, exact):
    """result is exact rounded to nearest, clear of underflow."""
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {exact})))")
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {result})))")


def divides(lemma, result, top, bottom):
    """result is top / bottom rounded to nearest, bottom above zero: |result bottom - top| <= u
    |top|, the quotient's rounding times bottom."""
    lemma.fact(f"(<= (absr (- (* {result} {bottom}) {top})) (* u (absr {top})))")


def plain_under(tight=TIGHT):
    """Judged's lap under half a tolerance: lap < rnd(rnd(L / 2 Tight) size), L / 2 exact, makes
    Straight's rnd(lap / size) under L / 2."""
    lemma = Lemma("plain under", "lap", "size", "half", "under", "product", "shadow")
    lemma.fact("(>= lap 0.0)")
    lemma.fact("(> size 0.0)")
    lemma.fact("(> half 0.0)")
    rounds(lemma, "under", f"(* half {rational(tight)})")
    rounds(lemma, "product", "(* under size)")
    lemma.fact("(< lap product)")
    divides(lemma, "shadow", "lap", "size")
    lemma.goal("(< shadow half)")
    return lemma


def plain_over(loose=LOOSE):
    """Judged's lap over 32 tolerances: lap > rnd(rnd(32 L Loose) size) makes Straight's rnd(lap
    / size) over 32 L."""
    lemma = Lemma("plain over", "lap", "size", "length", "over", "product", "shadow")
    lemma.fact("(> size 0.0)")
    lemma.fact("(> length 0.0)")
    rounds(lemma, "over", f"(* 32.0 length {rational(loose)})")
    rounds(lemma, "product", "(* over size)")
    lemma.fact("(> lap product)")
    divides(lemma, "shadow", "lap", "size")
    lemma.goal("(> shadow (* 32.0 length))")
    return lemma


def placed_start(loose=LOOSE):
    """Judged's place clear of the start: far = rnd(|n| reach) over edge = rnd(rnd(wide d) Loose),
    d = |turn| above zero, makes Straight's place p = rnd(n / d) give rnd(|p| reach) over wide."""
    lemma = Lemma("placed start", "n", "d", "reach", "wide", "far", "base", "edge", "place",
                  "product")
    lemma.fact("(> d 0.0)")
    lemma.fact("(> reach 0.0)")
    lemma.fact("(> wide 0.0)")
    rounds(lemma, "far", "(* (absr n) reach)")
    rounds(lemma, "base", "(* wide d)")
    rounds(lemma, "edge", f"(* base {rational(loose)})")
    lemma.fact("(> far edge)")
    divides(lemma, "place", "n", "d")
    rounds(lemma, "product", "(* (absr place) reach)")
    lemma.goal("(> product wide)")
    return lemma


def placed_end(share=2.0 ** -50):
    """Judged's place clear of the end: rnd(rnd(|n - d|) reach) over rnd(edge + rnd(far 2^-50))
    makes Straight's rnd(|rnd(p - 1)| reach) over wide, p = rnd(n / d): the 2^-50 |n| reach takes
    the rounding of n / d before 1 is taken from it. Products are named as Judged names them: far
    = rnd(|n| reach), edge = rnd(rnd(wide d) Loose)."""
    lemma = Lemma("placed end", "n", "d", "reach", "wide", "far", "base", "edge", "apart",
                  "gap", "slack", "limit", "place", "less", "product")
    lemma.fact("(> d 0.0)")
    lemma.fact("(> reach 0.0)")
    lemma.fact("(> wide 0.0)")
    rounds(lemma, "far", "(* (absr n) reach)")
    rounds(lemma, "base", "(* wide d)")
    rounds(lemma, "edge", f"(* base {rational(LOOSE)})")
    rounds(lemma, "apart", "(- n d)")
    rounds(lemma, "gap", "(* (absr apart) reach)")
    lemma.fact(f"(= slack (* far {rational(share)}))")
    rounds(lemma, "limit", "(+ edge slack)")
    lemma.fact("(> gap limit)")
    divides(lemma, "place", "n", "d")
    rounds(lemma, "less", "(- place 1.0)")
    rounds(lemma, "product", "(* (absr less) reach)")
    lemma.goal("(> product wide)")
    return lemma


def hit():
    """Judged's hit, n above 0 and below d, is Straight's, the place p = rnd(n / d) above 0 and
    below 1, where both tests above clear it of 0 and 1: p is not 1, and rounding is monotone."""
    lemma = Lemma("hit", "n", "d", "place")
    lemma.fact("(> d 0.0)")
    divides(lemma, "place", "n", "d")
    lemma.fact("(=> (< n d) (<= place 1.0))")
    lemma.fact("(=> (>= n d) (>= place 1.0))")
    lemma.fact("(not (= place 1.0))")
    lemma.fact("(not (= n 0.0))")
    lemma.goal("(= (> n 0.0) (> place 0.0))")
    lemma.goal("(= (and (> n 0.0) (< n d)) (and (> place 0.0) (< place 1.0)))")
    return lemma


def carry():
    """The quotient's carry: t d = n and q e = m for values n, d within En and Ed of m, e make |t -
    q| |d| <= En + |q| Ed."""
    lemma = Lemma("carry", "t", "d", "n", "q", "e", "m", "En", "Ed")
    lemma.fact("(= (* t d) n)")
    lemma.fact("(= (* q e) m)")
    lemma.fact("(not (= d 0.0))")
    lemma.fact("(<= (absr (- n m)) En)")
    lemma.fact("(<= (absr (- d e)) Ed)")
    lemma.goal("(<= (* (absr (- t q)) (absr d)) (+ En (* (absr q) Ed)))")
    return lemma


def inverse(steep=True):
    """The reciprocal kept: d' within 2^-100 Pd and 8 h of d, its high part within u of it, I =
    rnd(1 / high) with |I high - 1| <= u, and |d| >= Pd / 8 (1 - 2^-48), at least 2^-209: then I
    |d| is within 2^-50 of 1, d taken above zero as the sign changes nothing."""
    lemma = Lemma("inverse", "I", "d", "dH", "dp", "Pd")
    lemma.fact("(> d 0.0)")
    lemma.fact("(> dH 0.0)")
    lemma.fact("(>= Pd 0.0)")
    lemma.fact("(<= (absr (- dp dH)) (* u dH))")
    lemma.fact(f"(<= (absr (- dp d)) (+ (* {rational(2.0 ** -100)} Pd) (* 8.0 h)))")
    if steep:
        lemma.fact(f"(>= (* 8.0 d) (* Pd {rational(TIGHT)}))")
    lemma.fact(f"(>= d {rational(2.0 ** -209)})")
    lemma.fact("(<= (absr (- (* I dH) 1.0)) u)")
    lemma.goal(f"(>= (* I d) (- 1.0 {rational(2.0 ** -50)}))")
    lemma.goal(f"(<= (* I d) (+ 1.0 {rational(2.0 ** -50)}))")
    return lemma


def crossed_terms():
    """Crossed's bound covers |t - q|, term by term. By carry and inverse, |t - q'| <= (En + |q'|
    Ed) I / (1 - 2^-50), En = 2^-100 Pn + 8 h, Ed = 2^-100 Pd + 8 h; by over, |q - q'| <= 2^-101 S +
    3 h (1 + I) and |q'| <= c S + 6 h (1 + I), c = 1 + 4.1 u + 2^-101, S = |q's high part|. The
    bound, its five roundings each taken down, is at least (1 - u)^5 (Grain Loose (Pn I + S Pd I)
    + Grain S). Each of Pn I, S Pd I and S needs no more than a quarter of what the bound gives it;
    every h term is at most 2^-850, with I <= 2^210, Pd I <= 8.1, S <= 1.01, under the three
    quarters left of Grain Pn I, Pn I at least 2^-46 for t at least 2^-45."""
    lean, low = Fraction(1, 2 ** 100), Fraction(1, 2 ** 101)
    grain, loose = Fraction(GRAIN), Fraction(LOOSE)
    given = (1 - U) ** 5
    divided = 1 / (1 - Fraction(1, 2 ** 50))
    reach = 1 + Fraction(41, 10) * U + low
    held = lean * divided <= grain * loose * given / 4
    held = held and lean * reach * divided <= grain * loose * given / 4
    held = held and low <= grain * given / 4
    inverse_most, product_most, size_most = Fraction(2 ** 210), Fraction(81, 10), Fraction(101, 100)
    unders = (8 * H * inverse_most * divided
              + 8 * H * (reach * size_most + 6 * H * (1 + inverse_most)) * inverse_most * divided
              + 6 * H * (1 + inverse_most) * lean * product_most * divided
              + 3 * H * (1 + inverse_most))
    held = held and unders <= Fraction(1, 2 ** 850)
    left = grain * loose * given * 3 / 4 * Fraction(1, 2 ** 46)
    held = held and unders <= left
    return verdict("crossed terms", held, f"h terms {float(unders):.3g}, left {float(left):.3g}")


def crossed_certain(share=0.6):
    """Crossed's certificate: the coordinate start + s t, s = sH + sL the exact step, |s| <= (1 +
    u) |sH|; the place within the bound b of q, so s q within (1 + u) |sH| b of s t, written z =
    |sH| b; Moved within 50 u^2 Size and 3 h of start + s q; the bound Certain is given B =
    rnd(rnd(2 z) + rnd(Grain Size)), 2 z exact. Then the true coordinate is within 0.6 B of value +
    rest, and B is at least 2^-99 |value| for Size >= |value| (1 - 2^-60) and |value| at least
    Tiny."""
    lemma = Lemma("crossed certain", "z", "size", "value", "error", "first", "second", "bound")
    lemma.fact("(>= z 0.0)")
    lemma.fact(f"(>= size (* (absr value) (- 1.0 {rational(2.0 ** -60)})))")
    lemma.fact(f"(>= (absr value) (/ 1.0 {2 ** 960}.0))")
    lemma.fact("(<= error (+ (* (+ 1.0 u) z) (* 50.0 u u size) (* 3.0 h)))")
    lemma.fact("(= first (* 2.0 z))")
    rounds(lemma, "second", f"(* {rational(GRAIN)} size)")
    rounds(lemma, "bound", "(+ first second)")
    lemma.goal(f"(<= error (* {share} bound))")
    lemma.goal(f"(>= bound (* {rational(2.0 ** -99)} (absr value)))")
    return lemma


def ends(share=0.25):
    """Ends: M's coordinate as Moved gives it, value mV and rest mR within 32.5 u^2 mS and 3 h of
    its parts' exact M, mS = |start| + |pH| at least |mV| (1 - 2^-60); H's in double-double within
    17.3 u^2 |X| and 3 h, X the exact high product, |X| <= (1 + u) |hH| + 3 h (lean.py's times by
    a quotient: chord.hpp's comment said 13.1); TwoSum of mV and hH exact, s + t; the low part l =
    rnd(t + rnd(mR + hL)); value = rnd(s + l), rest = rnd(rnd(s - value) + l). Then value + rest
    is within a quarter of Grain (mS + |hH|) and 9 h of M + H, the rest of the bound Certain is
    given left to the places' share, Off, at most half of it."""
    lemma = Lemma("ends", "M", "mV", "mR", "mS", "Hx", "hH", "hL", "X", "s", "t", "c", "l",
                  "value", "w", "rest")
    lemma.fact(f"(>= mS (* (absr mV) (- 1.0 {rational(2.0 ** -60)})))")
    lemma.fact("(<= (absr (- (+ mV mR) M)) (+ (* 32.5 u u mS) (* 3.0 h)))")
    lemma.fact("(<= (absr mR) (+ (* 2.0 u (absr mV)) (* 4.0 h)))")
    lemma.fact("(<= (absr (- (+ hH hL) Hx)) (+ (* 17.3 u u (absr X)) (* 3.0 h)))")
    lemma.fact("(<= (absr X) (+ (* (+ 1.0 u) (absr hH)) (* 3.0 h)))")
    lemma.fact("(<= (absr hL) (+ (* u (absr hH)) (* 3.0 h)))")
    lemma.fact("(= (+ s t) (+ mV hH))")
    rounds(lemma, "s", "(+ mV hH)")
    rounds(lemma, "c", "(+ mR hL)")
    rounds(lemma, "l", "(+ t c)")
    rounds(lemma, "value", "(+ s l)")
    rounds(lemma, "w", "(- s value)")
    rounds(lemma, "rest", "(+ w l)")
    lemma.goal(f"(<= (absr (- (+ value rest) (+ M Hx))) "
               f"(+ (* {share} {rational(GRAIN)} (+ mS (absr hH))) (* 9.0 h)))")
    return lemma


def chorded_disc(kept=True):
    """Chorded's OffD: b, a and c in lean double-double within eB, eA and eC of the exact ones,
    eB at most offB / 4 (Dot's 2^-100 of b's size, offB its Grain), eA at most offA / 4 (1 + 3 u),
    eC at most offC / 2 (Dot and Plus within 2^-100 and 2^-104 of their sizes); |b|, |a| and |c|
    within (1 + u) of their high parts B, A and C. By sure.py's carried, b^2 is within 2 |b| eB +
    eB^2 of the exact square and a c within |a| eC + |c| eA + eA eC of the exact product; Times
    leaves lead and rest within 2^-102 of B^2 and A C, at most lead's and rest's high parts LH
    and RH (1 + 2 u), and Plus their difference within 2^-104 (LH + RH). The products of two
    numbers are named: T1 = B offB, T2 = offB^2, T3 = A offC, T4 = C offA, T5 = offA offC. Then
    OffD, its sum rounded term by term in the code's order and taken by Loose, covers the
    discriminant's error, and 3 h for each Times, LH + RH past 2^-800 as the lengths and radii past
    2^-103 that Reached leaves make it."""
    lemma = Lemma("chorded disc", "T1", "T2", "T3", "T4", "T5", "LH", "RH", "squared", "product",
                  "error", "p1", "p2", "p3", "p4", "p5", "s2", "s3", "s4", "s5", "grown", "tail",
                  "off")
    for name in ("T1", "T2", "T3", "T4", "T5", "LH", "RH", "squared", "product"):
        lemma.fact(f"(>= {name} 0.0)")
    lemma.fact("(<= squared (* LH (+ 1.0 (* 2.0 u))))")
    lemma.fact("(<= product (* RH (+ 1.0 (* 2.0 u))))")
    carried = ("(+ (* 2.0 (+ 1.0 u) T1 0.25) (* T2 0.0625) (* (+ 1.0 u) T3 0.5)"
               " (* (+ 1.0 u) T4 0.25 (+ 1.0 (* 3.0 u))) (* T5 0.125 (+ 1.0 (* 3.0 u))))")
    own = (f"(+ (* {rational(2.0 ** -102)} (+ squared product)) "
           f"(* {rational(2.0 ** -104)} (+ LH RH)) (* 6.0 h))")
    lemma.fact(f"(<= error (+ {carried} {own}))")
    rounds(lemma, "p1", "(* 2.0 T1)")
    rounds(lemma, "p2", "T2")
    rounds(lemma, "p3", "T3")
    if kept:
        rounds(lemma, "p4", "T4")
    else:
        lemma.fact("(= p4 0.0)")
    rounds(lemma, "p5", "T5")
    rounds(lemma, "s2", "(+ p1 p2)")
    rounds(lemma, "s3", "(+ s2 p3)")
    rounds(lemma, "s4", "(+ s3 p4)")
    rounds(lemma, "s5", "(+ s4 p5)")
    rounds(lemma, "grown", f"(* s5 {rational(LOOSE)})")
    rounds(lemma, "tail", f"(* {rational(GRAIN)} (+ LH RH))")
    rounds(lemma, "off", "(+ grown tail)")
    lemma.fact(f"(>= (+ LH RH) {rational(2.0 ** -800)})")
    lemma.goal("(<= error off)")
    return lemma


MUTANTS = (("plain under without Tight", lambda: plain_under(tight=1.0)),
           ("plain over without Loose", lambda: plain_over(loose=1.0)),
           ("placed start without Loose", lambda: placed_start(loose=1.0)),
           ("placed end without 2^-50 |n|", lambda: placed_end(share=0.0)),
           ("inverse without the sine", lambda: inverse(steep=False)),
           ("crossed certain at 0.5", lambda: crossed_certain(share=0.5)),
           ("ends at 2^-8 Grain", lambda: ends(share=2.0 ** -8)),
           ("chorded disc without |c| offA", lambda: chorded_disc(kept=False)))

LEMMAS = (plain_under, plain_over, placed_start, placed_end, hit, carry, inverse, crossed_certain,
          ends, chorded_disc)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, checks=(crossed_terms(),), mutants=MUTANTS))
