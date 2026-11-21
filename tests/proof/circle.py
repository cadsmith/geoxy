#!/usr/bin/env python3
"""Machine-checks the circle test's compensated second stage and its decisions with Z3.

    python3 -I tests/proof/circle.py [z3]

The proof on Compensated in source/exact/circle.cpp is split into lemmas, each a few operations of
the code as written. A rounded result r of an exact value z is bound to it by the two axioms of
round to nearest without underflow or overflow, |r - z| <= u |r| and |r - z| <= u |z|, u = 2^-53;
an error-free transformation adds that the double and its roundoff sum to z exactly. Each lemma's
conclusions are the next one's facts, and Z3 must find every conclusion's negation unsatisfiable,
and the facts alone satisfiable, so that no lemma holds for want of a case.
The decisions read the code's constants as the doubles the compiler makes of them.

That no product underflows or overflows inside the precondition, where the axioms hold, is a count
of binary places, checked first. Z3 names the first argument, z3 on the path by default; the exit
status is 1 when a lemma is not proved.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
import program  # pylint: disable=wrong-import-position
from lemma import (CIRCLE_BOUND_B, CIRCLE_BOUND_C, CIRCLE_BOUND_COMPENSATED, RESULT_BOUND, Lemma,
                   growth, power, rational, run)

THIRD_BOUND = CIRCLE_BOUND_C + CIRCLE_BOUND_COMPENSATED


def lift():
    """LiftNear: x^2 + y^2 to within 4u^2(1 + 6u) of it, its low part at most 2u(1 + 5u) of it."""
    lemma = Lemma("lift", "X2", "Y2", "xB", "xS", "yB", "yS", "hB", "hS", "s", "low")
    lemma.fact("(>= X2 0.0)")
    lemma.fact("(>= Y2 0.0)")
    lemma.split("xB", "xS", "X2")
    lemma.split("yB", "yS", "Y2")
    lemma.split("hB", "hS", "(+ xB yB)")
    lemma.rounded("s", "(+ xS yS)")
    lemma.rounded("low", "(+ s hS)")
    whole = "(+ X2 Y2)"
    lemma.goal(f"(<= (absr (- {whole} (+ hB low))) (* 4.0 u u (+ 1.0 (* 6.0 u)) {whole}))")
    lemma.goal(f"(<= (absr low) (* 2.0 u (+ 1.0 (* 5.0 u)) {whole}))")
    lemma.goal(f"(and (<= 0.0 hB) (<= hB (* (+ 1.0 u) (+ 1.0 u) {whole})))")
    v = program.compensated()
    lemma.code({"X2": v["axB"], "xB": v["axB"], "xS": v["axS"], "Y2": v["ayB"], "yB": v["ayB"],
                "yS": v["ayS"], "hB": v["ahB"], "hS": v["ahS"], "s": v["als"], "low": v["aLo"]})
    return lemma


def minor():
    """MinorNear: bx cy - by cx to within 4u^2(1 + 6u) of N = |bx cy| + |by cx|."""
    lemma = Lemma("minor", "T1", "T2", "f1", "e1", "f2", "e2", "mB", "mS", "s", "low")
    lemma.split("f1", "e1", "T1")
    lemma.split("f2", "e2", "T2")
    lemma.split("mB", "mS", "(- f1 f2)")
    lemma.rounded("s", "(- e1 e2)")
    lemma.rounded("low", "(+ s mS)")
    whole = "(+ (absr T1) (absr T2))"
    lemma.goal(f"(<= (absr (- (- T1 T2) (+ mB low))) (* 4.0 u u (+ 1.0 (* 6.0 u)) {whole}))")
    lemma.goal(f"(<= (absr low) (* 2.0 u (+ 1.0 (* 5.0 u)) {whole}))")
    lemma.goal(f"(<= (absr mB) (* (+ 1.0 u) (+ 1.0 u) {whole}))")
    lemma.goal(f"(<= (absr (- T1 T2)) {whole})")
    v = program.compensated()
    lemma.code({"T1": v["af1"], "f1": v["af1"], "e1": v["ae1"], "T2": v["af2"], "f2": v["af2"],
                "e2": v["ae2"], "mB": v["amB"], "mS": v["amS"], "s": v["ams"], "low": v["aMlo"]})
    return lemma


def product():
    """|a| <= alpha A and |b| <= beta B give |a b| <= alpha beta A B: the bound of every product of
    two bounded values below, proved once for any alpha and beta."""
    lemma = Lemma("product", "a", "b", "A", "B", "alpha", "beta")
    for value in ("A", "B", "alpha", "beta"):
        lemma.fact(f"(>= {value} 0.0)")
    lemma.fact("(<= (absr a) (* alpha A))")
    lemma.fact("(<= (absr b) (* beta B))")
    lemma.goal("(<= (absr (* a b)) (* alpha beta A B))")
    return lemma


def term_algebra():
    """A lift times a minor, split by the parts the code keeps: with eL = L - H - Lo and
    eM = M - Mh - Mlo, L M - (H Mh + H Mlo + Lo Mh) = Lo Mlo + eL (M - eM) + L eM."""
    lemma = Lemma("term algebra", "L", "H", "Lo", "M", "Mh", "Mlo")
    eL = "(- L H Lo)"
    eM = "(- M Mh Mlo)"
    lemma.goal(f"(= (- (* L M) (+ (* H Mh) (* H Mlo) (* Lo Mh)))"
               f" (+ (* Lo Mlo) (* {eL} (- M {eM})) (* L {eM})))")
    return lemma


def term():
    """TermNear: a lift times a minor to within 25u^2(1 + 10u) L N. The products are values of
    their own, bounded by the product lemma from the lift's and the minor's bounds: |H| <= (1 + u)^2
    L, |Lo| <= 2u(1 + 5u) L, |eL| <= 4u^2(1 + 6u) L, and alike for the minor against N, with
    |M - eM| <= (1 + 4u^2(1 + 6u)) N."""
    lemma = Lemma("term", "LN", "LM", "HMh", "HMlo", "LoMh", "LoMlo", "eLM", "LeM", "p", "pe",
                  "c1", "c2", "c3", "low")
    high = "(* (+ 1.0 u) (+ 1.0 u))"
    part = "(* 2.0 u (+ 1.0 (* 5.0 u)))"
    error = "(* 4.0 u u (+ 1.0 (* 6.0 u)))"
    lemma.fact("(>= LN 0.0)")
    lemma.fact(f"(<= (absr HMh) (* {high} {high} LN))")
    lemma.fact(f"(<= (absr HMlo) (* {high} {part} LN))")
    lemma.fact(f"(<= (absr LoMh) (* {part} {high} LN))")
    lemma.fact(f"(<= (absr LoMlo) (* {part} {part} LN))")
    lemma.fact(f"(<= (absr eLM) (* {error} (+ 1.0 {error}) LN))")
    lemma.fact(f"(<= (absr LeM) (* 1.0 {error} LN))")
    lemma.fact("(= (- LM (+ HMh HMlo LoMh)) (+ LoMlo eLM LeM))")
    lemma.split("p", "pe", "HMh")
    lemma.rounded("c1", "HMlo")
    lemma.rounded("c2", "LoMh")
    lemma.rounded("c3", "(+ c1 c2)")
    lemma.rounded("low", "(+ pe c3)")
    lemma.goal("(<= (absr (- LM (+ p low))) (* 25.0 u u (+ 1.0 (* 10.0 u)) LN))")
    lemma.goal("(<= (absr low) (* 5.0 u (+ 1.0 (* 11.0 u)) LN))")
    lemma.goal("(<= (absr p) (* (+ 1.0 (* 6.0 u)) LN))")
    v = program.compensated()
    lemma.code({"HMh": v["ap"], "p": v["ap"], "pe": v["ape"], "HMlo": v["ac1"], "c1": v["ac1"],
                "LoMh": v["ac2"], "c2": v["ac2"], "c3": v["ac3"], "low": v["alo"]})
    return lemma


def total():
    """Compensated's sum: E within u |E| + 46u^2(1 + 15u) P' of D, P' = sum L N."""
    names = []
    for index in (1, 2, 3):
        names += [f"T{index}", f"Q{index}", f"p{index}", f"lo{index}"]
    lemma = Lemma("sum", *names, "o", "os", "t", "ts", "a1", "a2", "a3", "low", "E")
    for index in (1, 2, 3):
        lemma.fact(f"(>= Q{index} 0.0)")
        lemma.fact(f"(<= (absr (- T{index} (+ p{index} lo{index})))"
                   f" (* 25.0 u u (+ 1.0 (* 10.0 u)) Q{index}))")
        lemma.fact(f"(<= (absr lo{index}) (* 5.0 u (+ 1.0 (* 11.0 u)) Q{index}))")
        lemma.fact(f"(<= (absr p{index}) (* (+ 1.0 (* 6.0 u)) Q{index}))")
    lemma.split("o", "os", "(+ p1 p2)")
    lemma.split("t", "ts", "(+ o p3)")
    lemma.rounded("a1", "(+ lo1 lo2)")
    lemma.rounded("a2", "(+ a1 lo3)")
    lemma.rounded("a3", "(+ os ts)")
    lemma.rounded("low", "(+ a2 a3)")
    lemma.rounded("E", "(+ t low)")
    lemma.goal("(<= (absr (- (+ T1 T2 T3) E))"
               " (+ (* u (absr E)) (* 46.0 u u (+ 1.0 (* 15.0 u)) (+ Q1 Q2 Q3))))")
    v = program.compensated()
    binding = {name: v[name] for name in ("o", "os", "t", "ts", "a1", "a2", "a3", "low", "E")}
    binding.update({"p1": v["ap"], "p2": v["bp"], "p3": v["cp"], "lo1": v["alo"],
                    "lo2": v["blo"], "lo3": v["clo"]})
    lemma.code(binding)
    return lemma


def corner_permanent():
    """A corner's term of the first stage's permanent: L N at most (1 + u)^5 of it."""
    lemma = Lemma("corner permanent", "X2", "Y2", "fx", "fy", "A", "T1", "T2", "F1", "F2", "S",
                  "k")
    lemma.fact("(>= X2 0.0)")
    lemma.fact("(>= Y2 0.0)")
    lemma.rounded("fx", "X2")
    lemma.rounded("fy", "Y2")
    lemma.rounded("A", "(+ fx fy)")
    lemma.rounded("F1", "T1")
    lemma.rounded("F2", "T2")
    lemma.rounded("S", "(+ (absr F1) (absr F2))")
    lemma.rounded("k", "(* S A)")
    lemma.goal(f"(<= (* (+ X2 Y2) (+ (absr T1) (absr T2))) (* {power('(+ 1.0 u)', 5)} k))")
    lemma.goal("(>= k 0.0)")
    v = program.circle_first()
    lemma.code({"X2": v["axx"], "fx": v["axx"], "Y2": v["ayy"], "fy": v["ayy"], "A": v["alift"],
                "T1": v["ap1"], "F1": v["ap1"], "T2": v["ap2"], "F2": v["ap2"], "S": v["apm"],
                "k": v["aq"]})
    return lemma


def permanent():
    """The permanent the first stage rounds: P' at most (1 + u)^7 of it."""
    lemma = Lemma("permanent", "Q1", "Q2", "Q3", "k1", "k2", "k3", "s", "P")
    for index in (1, 2, 3):
        lemma.fact(f"(>= Q{index} 0.0)")
        lemma.fact(f"(>= k{index} 0.0)")
        lemma.fact(f"(<= Q{index} (* {power('(+ 1.0 u)', 5)} k{index}))")
    lemma.rounded("s", "(+ k1 k2)")
    lemma.rounded("P", "(+ s k3)")
    lemma.goal(f"(<= (+ Q1 Q2 Q3) (* {power('(+ 1.0 u)', 7)} P))")
    lemma.goal(f"(<= (* 46.0 u u (+ 1.0 (* 15.0 u)) (+ Q1 Q2 Q3)) (* 47.0 u u P))")
    v = program.circle_first()
    lemma.code({"k1": v["aq"], "k2": v["bq"], "k3": v["cq"], "s": v["ps"], "P": v["perm"]})
    return lemma


def exact_decision():
    """Exact differences: |E| >= CircleBoundCompensated P decides D's sign."""
    lemma = Lemma("exact decision", "D", "E", "Pp", "P", "bound")
    lemma.fact("(>= Pp 0.0)")
    lemma.fact(f"(<= Pp (* {power('(+ 1.0 u)', 7)} P))")
    lemma.fact("(<= (absr (- D E)) (+ (* u (absr E)) (* 46.0 u u (+ 1.0 (* 15.0 u)) Pp)))")
    lemma.rounded("bound", f"(* {rational(CIRCLE_BOUND_COMPENSATED)} P)")
    lemma.decided("E", "bound")
    lemma.same_sign("D", "E")
    v = program.compensated()
    lemma.code({"P": program.given("permanent"), "bound": v["exact bound"], "E": v["E"]})
    return lemma


def lift_tails():
    """The exact differences' lift: within (2u + u^2) L of the rounded one's, split by order."""
    lemma = Lemma("lift tails", "x", "y", "tx", "ty")
    lemma.fact("(<= (absr tx) (* u (absr x)))")
    lemma.fact("(<= (absr ty) (* u (absr y)))")
    whole = "(+ (* x x) (* y y))"
    first = "(* 2.0 (+ (* x tx) (* y ty)))"
    second = "(+ (* tx tx) (* ty ty))"
    exact = "(+ (* (+ x tx) (+ x tx)) (* (+ y ty) (+ y ty)))"
    lemma.goal(f"(= (- {exact} {whole}) (+ {first} {second}))")
    lemma.goal(f"(<= (absr {first}) (* 2.0 u {whole}))")
    lemma.goal(f"(<= (absr {second}) (* u u {whole}))")
    lemma.goal(f"(<= (+ (absr (* x tx)) (absr (* y ty))) (* u {whole}))")
    return lemma


def minor_tails():
    """The exact differences' minor, split by order: (b + tb) x (c + tc) - b x c is
    (bx tcy + cy tbx) - (by tcx + cx tby), of first order, plus tbx tcy - tby tcx."""
    lemma = Lemma("minor tails", "bx", "by", "cx", "cy", "tbx", "tby", "tcx", "tcy")
    rounded = "(- (* bx cy) (* by cx))"
    exact = "(- (* (+ bx tbx) (+ cy tcy)) (* (+ by tby) (+ cx tcx)))"
    first = "(- (+ (* bx tcy) (* cy tbx)) (+ (* by tcx) (* cx tby)))"
    second = "(- (* tbx tcy) (* tby tcx))"
    lemma.goal(f"(= (- {exact} {rounded}) (+ {first} {second}))")
    return lemma


def minor_orders():
    """The two orders' bounds, 2u N and u^2 N. Each product is a value of its own, bounded by the
    product lemma from |tail| <= u |difference|: |bx tcy| <= u |bx cy|, |tbx tcy| <= u^2 |bx cy|,
    and alike."""
    lemma = Lemma("minor orders", "T1", "T2", "bxtcy", "cytbx", "bytcx", "cxtby", "tbxtcy",
                  "tbytcx")
    lemma.fact("(<= (absr bxtcy) (* u (absr T1)))")
    lemma.fact("(<= (absr cytbx) (* u (absr T1)))")
    lemma.fact("(<= (absr bytcx) (* u (absr T2)))")
    lemma.fact("(<= (absr cxtby) (* u (absr T2)))")
    lemma.fact("(<= (absr tbxtcy) (* u u (absr T1)))")
    lemma.fact("(<= (absr tbytcx) (* u u (absr T2)))")
    whole = "(+ (absr T1) (absr T2))"
    lemma.goal(f"(<= (+ (absr bxtcy) (absr cytbx) (absr bytcx) (absr cxtby)) (* 2.0 u {whole}))")
    lemma.goal(f"(<= (absr (- tbxtcy tbytcx)) (* u u {whole}))")
    return lemma


def orders():
    """A corner's exact term: its first order part, and the rest at most (6u^2 + 4u^3 + u^4) L N."""
    lemma = Lemma("orders", "L", "N", "M", "dL1", "dL2", "dM1", "dM2")
    lemma.fact("(>= L 0.0)")
    lemma.fact("(>= N 0.0)")
    lemma.fact("(<= (absr M) N)")
    lemma.fact("(<= (absr dL1) (* 2.0 u L))")
    lemma.fact("(<= (absr dL2) (* u u L))")
    lemma.fact("(<= (absr dM1) (* 2.0 u N))")
    lemma.fact("(<= (absr dM2) (* u u N))")
    exact = "(* (+ L dL1 dL2) (+ M dM1 dM2))"
    first = "(+ (* L dM1) (* dL1 M))"
    lemma.goal(f"(<= (absr (- {exact} (* L M) {first}))"
               f" (* (+ (* 6.0 u u) (* 4.0 u u u) (* u u u u)) L N))")
    lemma.goal(f"(<= (absr (- {exact} (* L M))) (* {growth(4)} L N))")
    lemma.goal(f"(<= (absr {first}) (* 4.0 u L N))")
    return lemma


def rounded_decision():
    """Rounded differences: |E| >= CircleBoundB P decides the determinant's sign."""
    lemma = Lemma("second decision", "det", "D", "E", "Pp", "P", "bound")
    lemma.fact("(>= Pp 0.0)")
    lemma.fact(f"(<= Pp (* {power('(+ 1.0 u)', 7)} P))")
    lemma.fact(f"(<= (absr (- det D)) (* {growth(4)} Pp))")
    lemma.fact("(<= (absr (- D E)) (+ (* u (absr E)) (* 46.0 u u (+ 1.0 (* 15.0 u)) Pp)))")
    lemma.rounded("bound", f"(* {rational(CIRCLE_BOUND_B)} P)")
    lemma.decided("E", "bound")
    lemma.same_sign("det", "E")
    v = program.compensated()
    lemma.code({"P": program.given("permanent"), "bound": v["second bound"], "E": v["E"]})
    return lemma


def cross_correction():
    """Correction's cross part, (bx tcy + cy tbx) - (by tcx + cx tby), three roundings deep."""
    lemma = Lemma("cross correction", "P1", "P2", "P3", "P4", "r1", "r2", "r3", "r4", "s1", "s2",
                  "q")
    for index in (1, 2, 3, 4):
        lemma.rounded(f"r{index}", f"P{index}")
    lemma.rounded("s1", "(+ r1 r2)")
    lemma.rounded("s2", "(+ r3 r4)")
    lemma.rounded("q", "(- s1 s2)")
    lemma.goal(f"(<= (absr (- q (- (+ P1 P2) (+ P3 P4))))"
               f" (* {growth(3)} (+ (absr P1) (absr P2) (absr P3) (absr P4))))")
    v = program.compensated()
    products = (program.mul(v["bx"], v["cytail"]), program.mul(v["cy"], v["bxtail"]),
                program.mul(v["by"], v["cxtail"]), program.mul(v["cx"], v["bytail"]))
    binding = {f"{kind}{index}": product for index, product in enumerate(products, 1)
               for kind in "Pr"}
    binding.update({"s1": program.add(products[0], products[1]),
                    "s2": program.add(products[2], products[3]), "q": v["aq"]})
    lemma.code(binding)
    return lemma


def two_rounded():
    """A sum of two rounded products, two roundings deep: the third stage's rounded lift and its
    2 (x tx + y ty) before the doubling."""
    lemma = Lemma("two rounded", "P1", "P2", "r1", "r2", "sum")
    lemma.rounded("r1", "P1")
    lemma.rounded("r2", "P2")
    lemma.rounded("sum", "(+ r1 r2)")
    lemma.goal(f"(<= (absr (- sum (+ P1 P2))) (* {growth(2)} (+ (absr P1) (absr P2))))")
    v = program.compensated()
    places = (((("ax", "ax"), ("ay", "ay")), "aLh"),
              ((("ax", "axtail"), ("ay", "aytail")), "aw"))
    for (x, y), name in places:
        first = program.mul(v[x[0]], v[x[1]])
        second = program.mul(v[y[0]], v[y[1]])
        lemma.code({"P1": first, "r1": first, "P2": second, "r2": second, "sum": v[name]})
    return lemma


def two_rounded_difference():
    """A difference of two rounded products, two roundings deep: the third stage's rounded
    minor."""
    lemma = Lemma("two rounded minor", "P1", "P2", "r1", "r2", "difference")
    lemma.rounded("r1", "P1")
    lemma.rounded("r2", "P2")
    lemma.rounded("difference", "(- r1 r2)")
    lemma.goal(f"(<= (absr (- difference (- P1 P2))) (* {growth(2)} (+ (absr P1) (absr P2))))")
    v = program.compensated()
    first = program.mul(v["bx"], v["cy"])
    second = program.mul(v["by"], v["cx"])
    lemma.code({"P1": first, "r1": first, "P2": second, "r2": second, "difference": v["am"]})
    return lemma


def first_part():
    """Correction's first part, the rounded lift times the cross part."""
    lemma = Lemma("first part", "L", "N", "Lh", "Q", "q", "S", "f")
    lemma.fact("(>= L 0.0)")
    lemma.fact("(>= N 0.0)")
    lemma.fact("(<= S (* 2.0 u N))")
    lemma.fact("(<= (absr Q) S)")
    lemma.fact(f"(<= (absr (- Lh L)) (* {growth(2)} L))")
    lemma.fact(f"(<= (absr (- q Q)) (* {growth(3)} S))")
    lemma.rounded("f", "(* Lh q)")
    lemma.goal(f"(<= (absr (- f (* L Q))) (* {growth(6)} 2.0 u L N))")
    lemma.goal("(<= (absr (* L Q)) (* 2.0 u L N))")
    v = program.compensated()
    lemma.code({"Lh": v["aLh"], "q": v["aq"], "f": v["af"]})
    return lemma


def second_part():
    """Correction's second part, 2 (x tx + y ty) times the rounded minor."""
    lemma = Lemma("second part", "L", "N", "W", "w", "V", "M", "m", "g")
    lemma.fact("(>= L 0.0)")
    lemma.fact("(>= N 0.0)")
    lemma.fact("(<= V (* u L))")
    lemma.fact("(<= (absr W) V)")
    lemma.fact("(<= (absr M) N)")
    lemma.fact(f"(<= (absr (- w W)) (* {growth(2)} V))")
    lemma.fact(f"(<= (absr (- m M)) (* {growth(2)} N))")
    lemma.rounded("g", "(* (* 2.0 w) m)")
    lemma.goal(f"(<= (absr (- g (* 2.0 W M))) (* {growth(5)} 2.0 u L N))")
    lemma.goal("(<= (absr (* 2.0 W M)) (* 2.0 u L N))")
    v = program.compensated()
    lemma.code({"w": v["aw"], "m": v["am"], "g": v["ag"]})
    return lemma


def corner_correction():
    """A corner's correction, the two parts summed: (1 + u)^7 - 1 of 4u L N off."""
    lemma = Lemma("corner correction", "Q", "F", "G", "f", "g", "c")
    lemma.fact("(>= Q 0.0)")
    lemma.fact(f"(<= (absr (- f F)) (* {growth(6)} 2.0 u Q))")
    lemma.fact(f"(<= (absr (- g G)) (* {growth(5)} 2.0 u Q))")
    lemma.fact("(<= (absr F) (* 2.0 u Q))")
    lemma.fact("(<= (absr G) (* 2.0 u Q))")
    lemma.rounded("c", "(+ f g)")
    lemma.goal(f"(<= (absr (- c (+ F G))) (* {growth(7)} 4.0 u Q))")
    v = program.compensated()
    lemma.code({"f": v["af"], "g": v["ag"], "c": v["acorr"]})
    return lemma


def corrections():
    """The three corrections summed: (1 + u)^9 - 1 of 4u P' off, the 36u^2 of the proof."""
    lemma = Lemma("corrections", "Q1", "Q2", "Q3", "F1", "F2", "F3", "c1", "c2", "c3", "s", "C")
    for index in (1, 2, 3):
        lemma.fact(f"(>= Q{index} 0.0)")
        lemma.fact(f"(<= (absr (- c{index} F{index})) (* {growth(7)} 4.0 u Q{index}))")
        lemma.fact(f"(<= (absr F{index}) (* 4.0 u Q{index}))")
    lemma.rounded("s", "(+ c1 c2)")
    lemma.rounded("C", "(+ s c3)")
    lemma.goal(f"(<= (absr (- C (+ F1 F2 F3))) (* {growth(9)} 4.0 u (+ Q1 Q2 Q3)))")
    lemma.goal(f"(<= (* {growth(9)} 4.0 u) (* 36.0 u u (+ 1.0 (* 10.0 u))))")
    v = program.compensated()
    lemma.code({"c1": v["acorr"], "c2": v["bcorr"], "c3": v["ccorr"], "s": v["cs"],
                "C": v["Corr"]})
    return lemma


def third_decision():
    """The third stage: |C| >= (CircleBoundC + CircleBoundCompensated) P + ResultBound |E|."""
    lemma = Lemma("third decision", "det", "D", "E", "F", "R", "Corr", "Pp", "P", "C", "k", "r",
                  "bound")
    lemma.fact("(>= Pp 0.0)")
    lemma.fact(f"(<= Pp (* {power('(+ 1.0 u)', 7)} P))")
    lemma.fact("(= det (+ D F R))")
    lemma.fact("(<= (absr R) (* (+ (* 6.0 u u) (* 4.0 u u u) (* u u u u)) Pp))")
    lemma.fact("(<= (absr (- D E)) (+ (* u (absr E)) (* 46.0 u u (+ 1.0 (* 15.0 u)) Pp)))")
    lemma.fact(f"(<= (absr (- Corr F)) (* {growth(9)} 4.0 u Pp))")
    lemma.rounded("C", "(+ E Corr)")
    lemma.rounded("k", f"(* {rational(THIRD_BOUND)} P)")
    lemma.rounded("r", f"(* {rational(RESULT_BOUND)} (absr E))")
    lemma.rounded("bound", "(+ k r)")
    lemma.decided("C", "bound")
    lemma.same_sign("det", "C")
    v = program.compensated()
    lemma.code({"E": v["E"], "Corr": v["Corr"], "C": v["C"], "P": program.given("permanent"),
                "k": v["third k"], "r": v["third r"], "bound": v["third bound"]})
    return lemma


LEMMAS = (lift, minor, product, term_algebra, term, total, corner_permanent, permanent,
          exact_decision, lift_tails, minor_tails, minor_orders, orders, rounded_decision,
          cross_correction, two_rounded, two_rounded_difference, first_part, second_part,
          corner_correction, corrections, third_decision)


def ranges():
    """Inside the precondition no product the stages form underflows or overflows.

    Every coordinate is zero or of a magnitude from 2^-150 up to 2^150, so a whole number of
    2^-202, its lowest place, and below 2^150. So is every difference and every roundoff of one,
    a difference of such numbers, below 2^151. A product of k of them is a whole number of
    2^(-202 k), its roundoff too, and every sum of such numbers, rounded or not: the stages form
    products of at most four, so every value that is not zero is at least 2^-808, and TwoProduct's
    roundoff, a whole number of a place at or above 2^-1074, is exact. The largest are sums of at
    most 48 products of four, below 2^(4 151 + 6).
    """
    lowest = -150 - 52
    highest = 151
    smallest = 4 * lowest
    largest = (4 * highest) + 6
    return smallest >= -1022 and largest < 1024, f"places from 2^{smallest} to 2^{largest}"


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


MUTANTS = (("exact threshold 40u^2",
            mutant(exact_decision, CIRCLE_BOUND_COMPENSATED=40.0 * 2.0 ** -106)),
           ("second bound (4 + 80u)u",
            mutant(rounded_decision, CIRCLE_BOUND_B=(4.0 + (80.0 * 2.0 ** -53)) * 2.0 ** -53)),
           ("third bound 85u^2", mutant(third_decision, THIRD_BOUND=85.0 * 2.0 ** -106)))

if __name__ == "__main__":
    sys.exit(run(LEMMAS, (ranges,), MUTANTS))
