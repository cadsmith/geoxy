#!/usr/bin/env python3
"""Machine-checks the bounds source/expansion.cpp's Rounded settles a quotient by: a number read in
two doubles by Rounding within 2^-104 (1 + 2^-50) of it however its parts cancel, Estimated's
double-double quotient within 2^-96 of its terms' size and its value's, or within 2^-88 of its
value where the top and the root's share cancel, each under half that, and Settled's certificate: a
double it settles is the one nearest the quotient, never a tie.

    python3 -I tests/proof/expansion.py [z3]

Rounding<true> reads the number h, held exactly in limbs, as H, the double nearest it, |h - H| <= u
|H| clear of the subnormals, and R, the 64 bits of the rest r = h - H below H's last bit, less a
whole step where H was rounded up, converted in two roundings: the bits below those 64, under 2^-116
|H|, left out, so |r - R| <= 2 u |R| + 2^-116 |H|. That the two hold H + R within 2^-104 (1 +
2^-50) of h is a lemma for Z3; without the nearest double, or within 2^-105, Z3 refutes it.

Estimated's two branches are bounded by expansion (expand.py): each double-double the exact value
times 1 + e, the e's in boxes, the computed quotient over the exact one less 1 a polynomial in them.
The factors are proved elsewhere: Rounding's here; Twofold::Add within 2^-100 of its terms' high
parts and Twofold::Divide and Root within 16 u^2 of their results (sure.py's sum, divide, root);
Twofold::Multiply within 2^-102 of its factors' high parts (lean.py's times); each with a few
h = 2^-1075 where a part underflows, which Estimated keeps far below by settling nothing whose terms
or value lie under 2^-900. A square root of 1 + b, |b| <= D, is 1 + p, |p| <= D / 2 (1 + D). In the
plain branch the top and the share are of one sign, or one is zero, so |t| + |s| is |t + s|; in the
other they are of opposite signs, so the conjugate's |t - s| is |t| + |s|, and the numerator t^2 -
s^2 = top^2 - bottom^2 root is exact. The bound's own roundings, at most five, are taken down.

Settled's certificate is a lemma for Z3: the estimate within half its bound, rnd(gap - 2 |L|) > 2 B
puts the quotient strictly within half the smaller step of the high part, so that the high part is
the double nearest it, as rounding to nearest takes it, and it is no halfway point.

The statements follow expansion.cpp's operations as written; they are not bound to program.py's
graphs, so translation.py does not check them against the compiler's IR. Each mutant drops a fact
the proof needs or shrinks a bound, and Z3, or a corner of the boxes worked in exact rationals, must
refute it. The exit status is 1 when a bound is not proved or a mutant is.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from expand import Fraction, Poly, U, peak, verdict
# pylint: disable-next=wrong-import-position
from lemma import Lemma, rational, run  # pylint: disable=wrong-import-position

DOUBLED = Fraction(1, 2 ** 104) * (1 + Fraction(1, 2 ** 50))
SUM = Fraction(1, 2 ** 100) / (1 - U)
TIMES = Fraction(1, 2 ** 102) / (1 - U) ** 2 + Fraction(1, 2 ** 170)
QUOTIENT = 16 * U * U / (1 - 16 * U * U) + Fraction(1, 2 ** 170)
ROOT = 16 * U * U / (1 - 16 * U * U)
HALF = DOUBLED / 2 * (1 + DOUBLED)


def doubled(nearest=True, grain=DOUBLED):
    """Doubled: the number, top the double nearest it, within u |top| of it, the rest = number -
    top read from its 64 leading bits within 2 u |read| and 2^-116 |top|; the two doubles' sum top
    + read within grain of the number."""
    lemma = Lemma("doubled", "number", "top", "rest", "read")
    lemma.fact("(= rest (- number top))")
    if nearest:
        lemma.fact("(<= (absr rest) (* u (absr top)))")
    tail = rational(float(Fraction(1, 2 ** 116)))
    lemma.fact(f"(<= (absr (- rest read)) (+ (* 2.0 u (absr read)) (* {tail} (absr top))))")
    lemma.goal(f"(<= (absr (- number (+ top read))) (* {rational(float(grain))} (absr number)))")
    return lemma


def settled(share=0.5):
    """Settled: v the quotient within share B of H + L, the estimate's two doubles, gap the
    smaller step either side of H, a positive double, t = rnd(gap - 2 |L|) > 2 B. Then v lies
    strictly within half gap of H."""
    lemma = Lemma("settled", "v", "H", "L", "B", "gap", "t")
    lemma.fact("(> gap 0.0)")
    lemma.fact("(>= B 0.0)")
    lemma.fact(f"(<= (absr (- v (+ H L))) (* {rational(share)} B))")
    lemma.fact("(<= (absr (- t (- gap (* 2.0 (absr L))))) (* u (absr t)))")
    lemma.fact("(> t (* 2.0 B))")
    lemma.goal("(< (absr (- v H)) (* 0.5 gap))")
    return lemma


def share_model():
    """The root's share over its exact value: bottom, the radicand's root and the product, each
    read with its error, the sign's product exact."""
    var = Poly.var
    return (1 + var("b1")) * (1 + var("p")) * (1 + var("r1")) * (1 + var("m1"))


def share_bounds():
    """The boxes of the share's errors."""
    return {"b1": DOUBLED, "p": HALF, "r1": ROOT, "m1": TIMES}


def plain_model():
    """The plain branch, Divide(Add(t, s), b), over the exact quotient: tau and 1 - tau the top's
    and the share's parts of their sum, of one sign; d3 for 1 / (1 + the divisor's error)."""
    var = Poly.var
    tau = var("tau")
    summed = tau * (1 + var("a1")) + (1 - tau) * share_model()
    return summed * (1 + var("a4")) * (1 + var("a5")) * (1 + var("d3")) - 1


def plain_bounds():
    """The boxes of the plain branch's errors."""
    bounds = share_bounds()
    bounds.update({"tau": Fraction(1), "a1": DOUBLED, "a4": SUM, "a5": QUOTIENT,
                   "d3": DOUBLED / (1 - DOUBLED)})
    return bounds


def plain(share=Fraction(1, 2)):
    """The plain branch within share of its computed bound, 2^-96 (scale + |value|): the error e of
    the quotient over its value, against that bound taken down over its roundings and its
    terms' errors, in units of the value."""
    error = plain_model().bound(plain_bounds())
    terms = max(DOUBLED, (share_model() - 1).bound(share_bounds()))
    scale = (1 - terms) * (1 - U) ** 3 / ((1 + U) * (1 + DOUBLED))
    value = (1 - error) / (1 + U)
    least = Fraction(1, 2 ** 96) * (1 - U) * (scale + value)
    held = error <= share * least
    return verdict("estimate, plain", held,
                   f"{float(error / U ** 2):.1f} u^2 of the value, the bound "
                   f"{float(least / error):.1f} times that")


def plain_mutant():
    """The plain branch's error past 2^-100 of the value at a corner: the bound is not sixteen
    times looser than needed."""
    bounds = plain_bounds()
    refuted = peak(plain_model(), bounds, {"tau": Fraction(1, 2)}) > Fraction(1, 2 ** 100)
    return verdict("mutant: estimate at 2^-100", refuted, "refuted" if refuted else "PROVED")


def cancel_model():
    """The cancelling branch's divisor, Multiply(Difference(t, s), b), over its exact value: tau and
    1 - tau the top's and the share's parts of |t| + |s|, of opposite signs."""
    var = Poly.var
    tau = var("tau")
    summed = tau * (1 + var("a1")) + (1 - tau) * share_model()
    return summed * (1 + var("a4")) * (1 + var("m3")) * (1 + var("a3")) - 1


def cancel(share=Fraction(1, 2)):
    """The cancelling branch within share of its bound, 2^-88 |value|: Doubled(gap) over the
    divisor, the divisor's error delta taken through 1 / (1 + delta), and Divide's."""
    bounds = share_bounds()
    bounds.update({"tau": Fraction(1), "a1": DOUBLED, "a4": SUM, "m3": TIMES, "a3": DOUBLED})
    delta = cancel_model().bound(bounds)
    inverse = delta / (1 - delta)
    var = Poly.var
    error = ((1 + var("g1")) * (1 + var("a5")) * (1 + var("d")) - 1).bound(
        {"g1": DOUBLED, "a5": QUOTIENT, "d": inverse})
    least = Fraction(1, 2 ** 88) * (1 - error) / (1 + U)
    held = error <= share * least
    return verdict("estimate, cancelling", held,
                   f"{float(error / U ** 2):.1f} u^2 of the value, the bound "
                   f"{float(least / error):.0f} times that")


MUTANTS = (("doubled without the nearest double", lambda: doubled(nearest=False)),
           ("doubled at 2^-105", lambda: doubled(grain=Fraction(1, 2 ** 105))),
           ("settled without its rounding's share", lambda: settled(share=1.0)))

LEMMAS = (doubled, settled)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, checks=(plain(), plain_mutant(), cancel()), mutants=MUTANTS))
