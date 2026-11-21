#!/usr/bin/env python3
"""Machine-checks with Z3 the bounds source/sure.hpp's numbers rest on: the certificate, the sum
of two double-doubles, the product's last step, and how a product's bound carries its operands'.

    python3 -I tests/proof/sure.py [z3]

A rounding is lemma.py's: |r - z| <= u |r| and |r - z| <= u |z| for a sum or difference, which
underflow leaves exact, and for a product clear of underflow; a product that underflows is off by
h = 2^-1075 at the most. TwoSum's double and roundoff sum to the exact value (rounding.v's
two_sum). A rounding is monotone (rounding.v's monotone), which turns the certificate's interval
into the double nearest every value in it, ties to even as the hardware rounds.

The quotient and the square root, Twofold::Divide and Twofold::Root, are bounded by expansion
(expand.py, as lean.py does): within 15 u^2 and 6.63 u^2 of their results, under the sixteenth of
QuotientGrain their bounds take. That a quotient's and a root's bound carries its operands' is
prepared.py's carry and abut.py's root step; the bounds' own few roundings, by Grow, are by hand.

Each lemma states the operations of sure.hpp and twofold.hpp as written, one fact a rounding;
they are not yet bound to program.py's graphs, so translation.py does not check them against the
compiler's IR. Each mutant shrinks a bound or drops a term the code keeps, and Z3 must refute it:
so the proofs need every part of the bounds they are given. The exit status is 1 when a lemma is
not proved or a mutant is.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from expand import Fraction, Poly, U, errors, peak, rounding, shares, verdict
# pylint: disable-next=wrong-import-position
from lemma import Lemma, rational, run  # pylint: disable=wrong-import-position

GROW = 1.0 + 2.0 ** -40
REST = 2.0 ** -50
SUM_GRAIN = 2.0 ** -100


def rounds(lemma, result, exact):
    """result is exact rounded to nearest: a sum, a difference, or a product clear of
    underflow."""
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {result})))")
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {exact})))")


def multiplies(lemma, result, exact):
    """result is the product exact rounded to nearest, underflow allowed, its sign kept."""
    relative = (f"(and (<= (absr (- {result} {exact})) (* u (absr {result})))"
                f" (<= (absr (- {result} {exact})) (* u (absr {exact}))))")
    absolute = f"(and (< (absr {exact}) uN) (<= (absr (- {result} {exact})) h))"
    lemma.fact(f"(or {relative} {absolute})")
    lemma.fact(f"(=> (>= {exact} 0.0) (>= {result} 0.0))")


def certificate(grow=GROW, rest=REST):
    """Rounded: a number H + L, |L| <= u |H| as a double-double is, its exact value within B of
    it. r = rnd(H + L) is H or a neighbour, within a factor of two of H, so H - r is exact
    (Sterbenz) and e = rnd((H - r) + L); b = rnd(rnd(B Grow) + rnd(|e| 2^-50)), p = rnd(e + b)
    and m = rnd(e - b). Then every value the bound allows lies between r + m and r + p: where both
    round to r, as the certificate asks, so does every value between, rounding being monotone."""
    lemma = Lemma("certificate", "H", "L", "B", "r", "e", "g", "k", "b", "p", "m", "v")
    lemma.fact("(<= (absr L) (* u (absr H)))")
    lemma.fact("(>= B 0.0)")
    rounds(lemma, "r", "(+ H L)")
    rounds(lemma, "e", "(+ (- H r) L)")
    # A double rounds to itself, and a sum that underflows is exact (rounding.v's sum_relative):
    # H + 0, e + 0, e - 0, 0 + b and 0 - b are what they add, and e, p and m are exact below the
    # normal doubles.
    lemma.fact("(=> (= L 0.0) (= r H))")
    lemma.fact("(=> (= b 0.0) (and (= p e) (= m e)))")
    lemma.fact("(=> (= e 0.0) (and (= p b) (= m (- b))))")
    lemma.fact("(=> (< (absr (+ (- H r) L)) uN) (= e (+ (- H r) L)))")
    lemma.fact("(=> (< (absr (+ e b)) uN) (= p (+ e b)))")
    lemma.fact("(=> (< (absr (- e b)) uN) (= m (- e b)))")
    # Rounding is monotone (rounding.v's monotone) and B and g are doubles: B Grow >= B and g + k >=
    # g round to no less.
    lemma.fact("(>= g B)")
    lemma.fact("(>= b g)")
    multiplies(lemma, "g", f"(* B {rational(grow)})")
    multiplies(lemma, "k", f"(* (absr e) {rational(rest)})")
    rounds(lemma, "b", "(+ g k)")
    rounds(lemma, "p", "(+ e b)")
    rounds(lemma, "m", "(- e b)")
    lemma.fact("(<= (absr (- v (+ H L))) B)")
    lemma.goal("(<= (+ r m) v)")
    lemma.goal("(<= v (+ r p))")
    return lemma


def total(grain=SUM_GRAIN):
    """operator+: Twofold::Add of two double-doubles, each |L| <= u |H|. TwoSum's s + t = xH + yH
    exactly; a = rnd(t + xL), c = rnd(a + yL); FastTwoSum: z = rnd(s + c), d = rnd(z - s), l =
    rnd(c - d). Then z + l is within SumGrain (|xH| + |yH|) of the exact sum, the share of the
    bound the sum's own rounding takes."""
    lemma = Lemma("sum", "xH", "xL", "yH", "yL", "s", "t", "a", "c", "z", "d", "l")
    lemma.fact("(<= (absr xL) (* u (absr xH)))")
    lemma.fact("(<= (absr yL) (* u (absr yH)))")
    rounds(lemma, "s", "(+ xH yH)")
    lemma.fact("(= (+ s t) (+ xH yH))")
    rounds(lemma, "a", "(+ t xL)")
    rounds(lemma, "c", "(+ a yL)")
    rounds(lemma, "z", "(+ s c)")
    rounds(lemma, "d", "(- z s)")
    rounds(lemma, "l", "(- c d)")
    lemma.goal(f"(<= (absr (- (+ z l) (+ xH xL yH yL))) "
               f"(* {rational(grain)} (+ (absr xH) (absr yH))))")
    return lemma


def product_tail(grain=2.0 ** -103):
    """operator*, its last step: Twofold::Multiply ends in FastTwoSum of p, TwoProduct's rounded
    xH yH, and c = rnd(rnd(q + xH yL) + xL yH), its terms' sum: z = rnd(p + c), d = rnd(z - p), l =
    rnd(c - d), with |c| <= 4 u |p|. Then z + l is within 2^-103 |p| of p + c. That p + c is
    itself within 2^-102 |xH| |yH| of the exact product, xL yL and three roundings of terms under
    u |xH yH| dropped, is lean.py's times, by expansion: Z3 does not settle it within ten
    minutes."""
    lemma = Lemma("product tail", "p", "c", "z", "d", "l")
    lemma.fact("(<= (absr c) (* 4.0 u (absr p)))")
    rounds(lemma, "z", "(+ p c)")
    rounds(lemma, "d", "(- z p)")
    rounds(lemma, "l", "(- c d)")
    lemma.goal(f"(<= (absr (- (+ z l) (+ p c))) (* {rational(grain)} (absr p)))")
    return lemma


def carried():
    """operator*'s bound for its operands' own: values X and Y within ex and ey of x and y have
    |X Y - x y| <= |x| ey + |y| ex + ex ey."""
    lemma = Lemma("carried", "x", "y", "X", "Y", "ex", "ey")
    lemma.fact("(<= (absr (- X x)) ex)")
    lemma.fact("(<= (absr (- Y y)) ey)")
    lemma.goal("(<= (absr (- (* X Y) (* x y))) (+ (* (absr x) ey) (* (absr y) ex) (* ex ey)))")
    return lemma


def divide_model():
    """Twofold::Divide of x by y, divided by Q0 = xH / yH, the units of xH within: g = rnd(xH /
    yH) = Q0 (1 + e1); Multiply of -g and y: TwoProduct p + q = -g yH exactly, k = rnd(-g yL) and
    h z1 where it underflows, the zero terms adding nothing, c = rnd(q + k), FastTwoSum exact;
    Add of x and that: TwoSum s + t = xH + mH, c2 = rnd(rnd(t + xL) + mL), rest = rnd(s + c2); out
    = g + rnd(rest / yH), FastTwoSum exact. xL = xH a and yL = yH b, within u. Returns g / Q0, the
    correction / Q0 and the error times (1 + b)."""
    var = Poly.var
    guess = rounding("e1")
    lead = -guess * rounding("e2")
    low = -guess - lead
    taken = -guess * var("b") * rounding("f1") + var("z1")
    tail = (low + taken) * rounding("f2")
    middle = (lead + tail) * rounding("e3")
    middleLow = lead + tail - middle
    head = (1 + middle) * rounding("e4")
    headLow = 1 + middle - head
    carried = ((headLow + var("a")) * rounding("f3") + middleLow) * rounding("f4")
    rest = (head + carried) * rounding("e5")
    step = rest * rounding("e6")
    return guess, step, (guess + step) * (1 + var("b")) - (1 + var("a"))


def divide(grain=16 * U * U):
    """Divide within 16 u^2 of its result and 2 h, the QuotientGrain's sixteenth: the error over
    1 - u for the factor 1 + b, against |g| (1 - u) at least, the correction at most 4 u |g| so
    that FastTwoSum is exact, and its own FastTwoSum's: c within 3 u of p."""
    guess, step, error = divide_model()
    bounds = errors("e1", "e2", "e3", "e4", "e5", "e6", "f1", "f2", "f3", "f4", "a", "b")
    bounds["z1"] = Fraction(1)
    share, free = shares(error, ("z1",), bounds)
    most = free / (1 - U) / (1 - U)
    corrected, _ = shares(step, ("z1",), bounds)
    held = most <= grain and step.parts(("z1",))[0].bound(bounds) <= 4 * U
    held = held and share["z1"] <= 2
    del guess, corrected
    return verdict("divide", held, f"{float(most / U ** 2):.2f} u^2 |result|")


def divide_mutant():
    """Divide's error past 2 u^2 at a corner: the bound is not eight times looser."""
    _, _, error = divide_model()
    bounds = errors("e1", "e2", "e3", "e4", "e5", "e6", "f1", "f2", "f3", "f4", "a", "b")
    refuted = peak(error, bounds, {"z1": 0}) > 2 * U * U * (1 + U)
    return verdict("mutant: divide at 2 u^2", refuted, "refuted" if refuted else "PROVED")


def root_model():
    """Twofold::Root of v, divided by r = sqrt(vH): g = r (1 + e1); Multiply of -g and g:
    TwoProduct p + q = -g^2 exactly, the zero terms adding nothing, FastTwoSum exact; Add of v and
    that: TwoSum s + t = vH + p, c2 = rnd(rnd(t + vL) + q), rest = rnd(s + c2); out = g + rnd(rest /
    rnd(2 g)), 2 g exact, FastTwoSum exact, 1 / (1 + e1) = 1 + d. vL = vH t, within u; the exact
    root r (1 + t / 2 + p), |p| <= t^2 / 8 (1 + 2^-50). Returns the error / r."""
    var = Poly.var
    guess = rounding("e1")
    lead = -(guess * guess) * rounding("e2")
    low = -(guess * guess) - lead
    head = (1 + lead) * rounding("e3")
    headLow = 1 + lead - head
    carried = ((headLow + var("t")) * rounding("f1") + low) * rounding("f2")
    rest = (head + carried) * rounding("e4")
    step = rest * Fraction(1, 2) * (1 + var("d")) * rounding("e5")
    return guess + step - (1 + var("t") * Fraction(1, 2) + var("p"))


def root(grain=16 * U * U):
    """Root within 16 u^2 of its result, the QuotientGrain's sixteenth, against g (1 - u)."""
    error = root_model()
    bounds = errors("e1", "e2", "e3", "e4", "e5", "f1", "f2", "t")
    bounds["d"] = U / (1 - U)
    bounds["p"] = U * U / 8 * (1 + Fraction(1, 2 ** 50))
    most = error.bound(bounds) / (1 - U)
    return verdict("root", most <= grain, f"{float(most / U ** 2):.2f} u^2 |result|")


def root_mutant():
    """Root's error past u^2 / 2 at a corner."""
    error = root_model()
    bounds = errors("e1", "e2", "e3", "e4", "e5", "f1", "f2", "t")
    fixed = {"d": 1 / (1 + U) - 1, "p": -U * U / 8, "e1": U}
    refuted = peak(error, bounds, fixed) > U * U / 2
    return verdict("mutant: root at u^2 / 2", refuted, "refuted" if refuted else "PROVED")


MUTANTS = (("certificate without the growth", lambda: certificate(grow=1.0)),
           ("certificate without the rest's share", lambda: certificate(rest=0.0)),
           ("sum grain 2^-106", lambda: total(grain=2.0 ** -106)),
           ("product tail grain 2^-110", lambda: product_tail(grain=2.0 ** -110)))

LEMMAS = (certificate, total, product_tail, carried)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, checks=(divide(), divide_mutant(), root(), root_mutant()),
                 mutants=MUTANTS))
