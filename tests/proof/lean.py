#!/usr/bin/env python3
"""Machine-checks the bounds source/lean.hpp's toolkit gives its callers: the lean cross and dot
products within 2^-100 of their terms' size (16 u^2 found, 32.2 claimed), the sum within 2^-104
and the product within 2^-102 of their high parts, a start moved along a step within 50 u^2 of its
terms' size (32.5 found), the quotient by a fused multiply-add's remainder within 2^-101 of its
high part (22 u^2 found, 22.2 claimed) and its low part within 4.1 u, the square root within
2^-103 of its high part (7.63 u^2 found) and its low part within 1.52 u, and the certificate:
a double Certain certifies is the one nearest the number, never a tie.

    python3 -I tests/proof/lean.py [z3]

The sum and the certificate are lemmas for Z3, a rounding being lemma.py's: |r - z| <= u |r| and
|r - z| <= u |z| for a sum or a difference, which underflow leaves exact; TwoSum's double and
roundoff sum to the exact value (rounding.v's two_sum); and a double rounds to value only where it
lies within the half gaps either side of value, as rounding to nearest does (Flocq's round_N_pt).

The others are proved by expansion (expand.py): each rounding of an exact z a factor 1 + e, |e| <=
u, the computed value less the exact one a polynomial in the e's and the inputs' low parts, each
in a box, bounded term by term. TwoProduct is exact for a product at least 2^-969 (rounding.v's
two_product_fma) and FastTwoSum where its smaller term is the smaller (fast_two_sum), checked over
the boxes. A product of a low part may underflow, off by h = 2^-1075 at the most, each such h
bounded at its own scale: the bounds carry a few h, which the callers keep far below their terms
(prepared.py). The square root's exact value is sqrt(v) (1 + t / 2 + p) for v + v t, |t| <= u, |p|
<= t^2 / 8 (1 + 2^-50), by Taylor's theorem with Lagrange's remainder.

Left by hand: a high product under 2^-969, where TwoProduct's roundoff may itself round, makes
every part of Cross, Dot, Times and Moved a number under 2^-960, each rounding within 2^-1012;
the callers' terms are above 2^-400, so that is far below their bounds.

The statements follow lean.hpp's operations as written; they are not bound to program.py's graphs,
so translation.py does not check them against the compiler's IR. Each mutant shrinks a bound, and
Z3, or a corner of the boxes worked in exact rationals, must refute it: so the bounds proved are
not loose past what the code may reach. The exit status is 1 when a bound is not proved or a
mutant is.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from expand import Fraction, Poly, U, errors, peak, rounding, shares, verdict
# pylint: disable-next=wrong-import-position
from lemma import Lemma, rational, run  # pylint: disable=wrong-import-position

H = Fraction(1, 2 ** 1075)


def rounds(lemma, result, exact):
    """result is a sum or difference exact rounded to nearest; underflow leaves it exact."""
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {result})))")
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {exact})))")


def two_sum(lemma, big, small, exact):
    """TwoSum: big is exact rounded, big + small is exact (rounding.v's two_sum)."""
    rounds(lemma, big, exact)
    lemma.fact(f"(= (+ {big} {small}) {exact})")


def products_model(sign):
    """Cross (sign -1) and Dot (sign 1) of exact displacements a and b: X = a0H b1H and Y = a1H b0H
    for the cross product, a0H b0H and a1H b1H for the dot, the scales; each low part a share of
    its high one, within u. TwoProduct splits X into pX = X (1 + e1) and lX = -X e1, and h z1 where
    it underflows, Y likewise; TwoSum splits pX +- pY into hH and hL; the four terms of a low part
    with a high one, m1 to m4, round to q1 to q4, with h z3 to z6 where they underflow; the low
    parts' products n1 and n2 are dropped; the roundoffs and the q's sum in the code's order, d0,
    s1, r2, s2, low, and w = rnd(hL + low); TwoSum of hH and w is exact. Returns the error."""
    var = Poly.var
    big, other = var("X"), var("Y")
    lead = big * rounding("e1") + var("z1")
    trail = other * rounding("e2") + var("z2")
    leadLow = big - lead
    trailLow = other - trail
    head = (lead + sign * trail) * rounding("e3")
    headLow = (lead + sign * trail) - head
    terms = [big * var("b1"), big * var("a0"), other * var("b0"), other * var("a1")]
    rounded = [term * rounding(f"f{index}") + var(f"z{index + 3}")
               for index, term in enumerate(terms)]
    first = (leadLow + sign * trailLow) * rounding("g1")
    second = (rounded[0] + rounded[1]) * rounding("g2")
    third = (first + second) * rounding("g3")
    fourth = (rounded[2] + rounded[3]) * rounding("g4")
    low = (third + sign * fourth) * rounding("g5")
    tail = (headLow + low) * rounding("g6")
    exact = (big * (1 + var("a0")) * (1 + var("b1"))
             + sign * other * (1 + var("a1")) * (1 + var("b0")))
    return head + tail - exact


PRODUCTS_UNDER = ("z1", "z2", "z3", "z4", "z5", "z6")


def products(sign=-1, grain=Fraction(1, 2 ** 100)):
    """The lean cross or dot product within 2^-100 (|pX| + |pY|) and 8 h of the exact one: the
    shares of X and Y, each over 1 - u for |X| <= |pX| / (1 - u), at most 2^-100, and the
    underflows' at most 8 in all."""
    name = "cross" if sign < 0 else "dot"
    error = products_model(sign)
    bounds = errors("e1", "e2", "e3", "f0", "f1", "f2", "f3", "g1", "g2", "g3", "g4", "g5", "g6",
                    "a0", "a1", "b0", "b1")
    bounds.update({under: Fraction(1) for under in PRODUCTS_UNDER})
    share, rest = shares(error, ("X", "Y") + PRODUCTS_UNDER, bounds)
    most = max(share["X"], share["Y"]) / (1 - U)
    held = rest == 0 and most <= grain and sum(share[z] for z in PRODUCTS_UNDER) <= 8
    return verdict(name, held, f"{float(most / U ** 2):.2f} u^2 P")


def products_mutant(sign=-1):
    """The cross or dot product's error past 2^-104 P at a corner of the boxes, X and Y of one
    size, cancelling for the cross product."""
    error = products_model(sign)
    bounds = errors("e1", "e2", "e3", "f0", "f1", "f2", "f3", "g1", "g2", "g3", "g4", "g5", "g6",
                    "a0", "a1", "b0", "b1")
    fixed = {"X": 1, "Y": 1 if sign < 0 else -1}
    fixed.update({under: 0 for under in PRODUCTS_UNDER})
    refuted = peak(error, bounds, fixed) > Fraction(1, 2 ** 104) * 2 / (1 - U)
    name = "cross" if sign < 0 else "dot"
    return verdict(f"mutant: {name} at 2^-104", refuted, "refuted" if refuted else "PROVED")


def plus(grain=2.0 ** -104):
    """Plus: TwoSum of the high parts, the low parts summed and added to its roundoff, TwoSum
    again; within 2^-104 (|aH| + |bH|) of the exact sum, its low part within u of its high."""
    lemma = Lemma("plus", "aH", "aL", "bH", "bL", "hH", "hL", "c", "w", "oH", "oL")
    lemma.fact("(<= (absr aL) (* u (absr aH)))")
    lemma.fact("(<= (absr bL) (* u (absr bH)))")
    two_sum(lemma, "hH", "hL", "(+ aH bH)")
    rounds(lemma, "c", "(+ aL bL)")
    rounds(lemma, "w", "(+ hL c)")
    two_sum(lemma, "oH", "oL", "(+ hH w)")
    lemma.goal(f"(<= (absr (- (+ oH oL) (+ aH aL bH bL))) "
               f"(* {rational(grain)} (+ (absr aH) (absr bH))))")
    lemma.goal("(<= (absr oL) (* u (absr oH)))")
    return lemma


def times_model():
    """Times, Twofold::Multiply, X = aH bH the scale, each low part a share of its high one within
    scale: TwoProduct splits X into p = X (1 + e1) and q = -X e1, exactly, X being at least 2^-969;
    the terms m1 = aH bL and m2 = aL bH round to k1 and k2, with h z1 and h z2 where they
    underflow; c = rnd(rnd(q + k1) + k2); the low parts' product dropped; FastTwoSum of p and c.
    Returns p, c and the error of p + c."""
    var = Poly.var
    big = var("X")
    lead = big * rounding("e1")
    low = big - lead
    first = big * var("b") * rounding("f1") + var("z1")
    second = big * var("a") * rounding("f2") + var("z2")
    tail = ((low + first) * rounding("g1") + second) * rounding("g2")
    return lead, tail, lead + tail - big * (1 + var("a")) * (1 + var("b"))


def times(grain=Fraction(1, 2 ** 102), wide=U, name="times"):
    """Times within 2^-102 |X| and 3 h of the exact product, for factors whose low parts are
    within u, and within 17.3 u^2 |X| where one is Over's, its low part within 4.1 u, as Ends and
    Cut take it: FastTwoSum is exact, c being at most 3.01 u |X| and 2 h, under |p|, for |X| at
    least 2^-969; so the low part is within u of the high."""
    lead, tail, error = times_model()
    bounds = errors("e1", "f1", "f2", "g1", "g2", "a", "b")
    bounds["a"] = wide
    bounds.update({"z1": Fraction(1), "z2": Fraction(1)})
    share, rest = shares(error, ("X", "z1", "z2"), bounds)
    reach, _ = shares(tail, ("X", "z1", "z2"), bounds)
    least = Fraction(1, 2 ** 969)
    fast = reach["X"] + (reach["z1"] + reach["z2"]) * H / least <= 1 - U
    held = rest == 0 and fast and share["X"] <= grain and share["z1"] + share["z2"] <= 3
    del lead
    return verdict(name, held, f"{float(share['X'] / U ** 2):.2f} u^2 |X|, "
                                  f"c {float(reach['X'] / U):.2f} u |X|")


def times_mutant():
    """Times's error past 2^-105 |X| at a corner of the boxes."""
    _, _, error = times_model()
    bounds = errors("e1", "f1", "f2", "g1", "g2", "a", "b")
    refuted = peak(error, bounds, {"X": 1, "z1": 0, "z2": 0}) > Fraction(1, 2 ** 105)
    return verdict("mutant: times at 2^-105", refuted, "refuted" if refuted else "PROVED")


def moved_model():
    """Moved, its scales the start and X = sH aH, X at least 2^-969 or zero: TwoProduct splits X
    into pH = X (1 + e1) and pL = -X e1; TwoSum splits start + pH into sH and sL; the step's low
    part a share s of its high, within u, the amount's a, within 4.1 u, as Over leaves it; the
    terms sH aL and sL aH round, with h z1 and h z2 where they underflow; the low parts sum in the
    code's order; value = rnd(sH + low), rest = rnd(rnd(sH - value) + low). Returns the error of
    value + rest."""
    var = Poly.var
    start, big = var("S"), var("X")
    lead = big * rounding("e1")
    leadLow = big - lead
    head = (start + lead) * rounding("e2")
    headLow = (start + lead) - head
    first = big * var("a") * rounding("f1") + var("z1")
    second = big * var("s") * rounding("f2") + var("z2")
    low = ((headLow + leadLow) * rounding("g1") + (first + second) * rounding("g2")) \
        * rounding("g3")
    value = (head + low) * rounding("g4")
    rest = ((head - value) * rounding("g5") + low) * rounding("g6")
    return value + rest - (start + big * (1 + var("a")) * (1 + var("s")))


def moved(grain=50 * U * U):
    """Moved within 50 u^2 (|start| + |pH|) and 3 h of the exact point, |X| <= |pH| / (1 - u)."""
    error = moved_model()
    bounds = errors("e1", "e2", "f1", "f2", "g1", "g2", "g3", "g4", "g5", "g6", "s")
    bounds["a"] = Fraction(41, 10) * U
    bounds.update({"z1": Fraction(1), "z2": Fraction(1)})
    share, rest = shares(error, ("S", "X", "z1", "z2"), bounds)
    most = max(share["S"], share["X"] / (1 - U))
    held = rest == 0 and most <= grain and share["z1"] + share["z2"] <= 3
    return verdict("moved", held, f"{float(most / U ** 2):.2f} u^2 (|start| + |pH|)")


def moved_mutant():
    """Moved's error past 10 u^2 (|start| + |pH|) at a corner of the boxes."""
    error = moved_model()
    bounds = errors("e1", "e2", "f1", "f2", "g1", "g2", "g3", "g4", "g5", "g6", "s")
    bounds["a"] = Fraction(41, 10) * U
    fixed = {"S": Fraction(1, 2), "X": 1, "z1": 0, "z2": 0}
    refuted = peak(error, bounds, fixed) > 10 * U * U * Fraction(3, 2) / (1 - U)
    return verdict("mutant: moved at 10 u^2", refuted, "refuted" if refuted else "PROVED")


def certain(share=0.6, least=2.0 ** -99):
    """Certain: value a double at least Tiny, so normal, its neighbours up and down 2 up and 2
    down away, each half gap between u |value| / 2 and u |value|; the number N within 0.6 bound
    of value + rest, the bound at least 2^-99 |value|. A double rounds to value only where it lies
    within the half gaps either side of it, as rounding to nearest does: so value + rnd(rest +
    bound) and value + rnd(rest - bound) rounding to value put both within them. Then N lies
    strictly within them: value is the double nearest N, and N is no tie."""
    lemma = Lemma("certain", "value", "rest", "bound", "N", "t1", "t2", "up", "down")
    lemma.fact(f"(>= (absr value) (/ 1.0 {2 ** 960}.0))")
    for gap in ("up", "down"):
        lemma.fact(f"(<= (* 0.5 u (absr value)) {gap})")
        lemma.fact(f"(<= {gap} (* u (absr value)))")
    lemma.fact(f"(>= bound (* {rational(least)} (absr value)))")
    lemma.fact(f"(<= (absr (- N (+ value rest))) (* {share} bound))")
    rounds(lemma, "t1", "(+ rest bound)")
    rounds(lemma, "t2", "(- rest bound)")
    for test in ("t1", "t2"):
        lemma.fact(f"(<= (- down) {test})")
        lemma.fact(f"(<= {test} up)")
    lemma.goal("(< N (+ value up))")
    lemma.goal("(> N (- value down))")
    return lemma


def over_model():
    """Over, divided by the high parts' quotient Q0 = tH / bH: the inputs t = tH (1 + a) and b =
    bH (1 + b), their low parts within u; inverse = (1 + e1) / bH; high = tH inverse (1 + e2) and
    h z2; rest = rnd(tH - high bH), the fused multiply-add, (1 + e3) and h z3 at tH's scale; s1 =
    rnd(rest + tL); k = rnd(high bL) and h z5 at tH's scale; s2 = rnd(s1 - k); low = rnd(s2
    inverse) and h z7. Returns high / Q0, low / Q0 and the exact quotient's error times (1 + b),
    each as a polynomial."""
    var = Poly.var
    top, bottom = var("a"), var("b")
    high = (rounding("e1") * rounding("e2")) + var("z2")
    rest = -(high - 1) * rounding("e3") + var("z3")
    first = (rest + top) * rounding("e4")
    taken = high * bottom * rounding("e5") + var("z5")
    second = (first - taken) * rounding("e6")
    low = second * rounding("e1") * rounding("e7") + var("z7")
    error = (high + low) * (1 + bottom) - (1 + top)
    return high, low, error


OVER_UNDER = ("z2", "z3", "z5", "z7")


def over(grain=Fraction(1, 2 ** 101), most=Fraction(41, 10) * U):
    """Over's bound: |high + low - t / b| <= 2^-101 |high| + 3 h (1 + |inverse|) and |low| <= 4.1
    u |high| + 3 h (1 + |inverse|): the free part bounded by expansion, over 1 - u for the
    factor 1 + b and over (1 - u)^2 for |Q0| <= (|high| + h) / (1 - u)^2; each underflow's share
    at its own scale, z2 and z7 at the quotient's, z3 and z5 at tH's, which is the quotient's
    times 1 / |bH|, within |inverse| (1 + u)."""
    _, low, error = over_model()
    bounds = errors("e1", "e2", "e3", "e4", "e5", "e6", "e7", "a", "b")
    bounds.update({name: Fraction(1) for name in OVER_UNDER})
    least = (1 - U) ** 2
    shrink = 1 / (1 - U)
    told = []
    held = True
    for value, limit, scale in ((error, grain, shrink), (low, most, 1)):
        share, free = shares(value, OVER_UNDER, bounds)
        held = held and free * scale <= limit * least
        held = held and (share["z2"] + share["z7"]) * scale + free <= 3
        held = held and (share["z3"] + share["z5"]) * scale * (1 + U) <= 3
        told.append(float(free * scale / least))
    return verdict("over", held, f"error {told[0] / float(U) ** 2:.2f} u^2 |high|, "
                                 f"low {told[1] / float(U):.2f} u |high|")


def over_mutant():
    """Over's error is past 2^-105 |Q0| at a corner of the boxes: the bound is no looser than
    four times what the code reaches."""
    _, _, error = over_model()
    values = {"e1": U, "e2": U, "e3": -U, "e4": U, "e5": -U, "e6": U, "e7": U, "a": -U, "b": U}
    values.update({name: 0 for name in OVER_UNDER})
    refuted = abs(error.value(values)) / (1 + U) > Fraction(1, 2 ** 105)
    return verdict("mutant: over at 2^-105", refuted, "refuted" if refuted else "PROVED")


def root_model():
    """Root, divided by r = sqrt(vH): vL = vH t, high = r (1 + e1), the root correctly rounded;
    inverse = (1 + e2) / high, 1 / (1 + e1) = 1 + d; rest = rnd(vH - high^2), the fused
    multiply-add, (1 + e3) and h z3 at vH's scale; sum = rnd(rest + vL); low = rnd(sum (0.5
    inverse)) and h z5 at r's scale, 0.5 inverse exact; FastTwoSum of high and low exact. The
    exact root is r (1 + t / 2 + p). Returns high / r, low / r and the error / r."""
    var = Poly.var
    high = rounding("e1")
    rest = -(high * high - 1) * rounding("e3") + var("z3")
    total = (rest + var("t")) * rounding("e4")
    low = total * Fraction(1, 2) * rounding("e2") * (1 + var("d")) * rounding("e5") + var("z5")
    exact = 1 + var("t") * Fraction(1, 2) + var("p")
    return high, low, (high + low) - exact


def root(grain=Fraction(1, 2 ** 103), most=Fraction(152, 100) * U):
    """Root's bound: within 2^-103 |high| and 2 h of sqrt(vH + vL), its low part at most 1.52 u
    |high| and h before FastTwoSum leaves it within u of the high: z3 at vH's scale, which is r
    times r, at most 2^-450 for the vH at least 2^-900 the callers keep, z5 at r's."""
    high, low, error = root_model()
    bounds = errors("e1", "e2", "e3", "e4", "e5", "t")
    bounds["d"] = U / (1 - U)
    bounds["p"] = U * U / 8 * (1 + Fraction(1, 2 ** 50))
    bounds["z3"] = Fraction(1, 2 ** 450)
    bounds["z5"] = Fraction(0)
    del high
    least = 1 - U
    free = error.bound(bounds)
    held = free <= grain * least
    _, shares = error.parts(("z5",))
    held = held and shares["z5"].bound(bounds) <= 1
    lowFree = low.bound(bounds)
    held = held and lowFree <= most * least
    told = f"error {float(free / U ** 2):.3f} u^2 r, low {float(lowFree / U):.3f} u r"
    return verdict("root", held, told)


def root_mutant():
    """Root's error is past 2^-106 r at a corner: the bound is no looser than eight times what
    the code may reach."""
    _, _, error = root_model()
    values = {"e1": U, "e2": U, "e3": U, "e4": U, "e5": U, "t": U, "z3": 0, "z5": 0}
    values["d"] = 1 / (1 + U) - 1
    values["p"] = -U * U / 8
    refuted = abs(error.value(values)) > Fraction(1, 2 ** 106)
    return verdict("mutant: root at 2^-106", refuted, "refuted" if refuted else "PROVED")


MUTANTS = (("plus at 2^-106", lambda: plus(2.0 ** -106)),
           ("certain at 1.01 of the bound", lambda: certain(share=1.01)),
           ("certain with the bound at 2^-107", lambda: certain(least=2.0 ** -107)))

LEMMAS = (plus, certain)

CHECKS = (products(-1), products_mutant(-1), products(1), products_mutant(1), times(),
          times(Fraction(1731, 100) * U * U, Fraction(41, 10) * U, "times by a quotient"),
          times_mutant(), moved(), moved_mutant())

if __name__ == "__main__":
    sys.exit(run(LEMMAS, checks=CHECKS + (over(), over_mutant(), root(), root_mutant()),
                 mutants=MUTANTS))
