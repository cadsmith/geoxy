#!/usr/bin/env python3
"""Machine-checks with Z3 that source/crossing.cpp's Bowline and source/crossings/bows.cpp's Hopped
choose, of two bulges' circles' two crossings, the one the exact way chooses: the nearer the hit,
by the distances of the two points each rounded once, the first where they are as near; and how
far Hopped's one Newton step may leave its crossing.

    python3 -I tests/proof/bows.py [z3]

Crossbows (source/crossings/bows.cpp) works each crossing's coordinates out of sure.hpp's numbers
only: each a value H + L, |L| <= u |H|, within its bound B of the exact coordinate X, as sure.py
proves of every operation it takes. The rounded coordinate R, the double nearest X, is within u |X|
of it. So R is within (1 + u) B + u (2 + u) |H| of H (coordinate lemma), under twice the 4 (B +
2^-50 |H|) Bowline adds to its slack for it.

A point's distance from the hit moves by no more than the point does, and the point by no more
than the sum of its coordinates' moves: the rounded point's distance t is within e, the sum of its
coordinates' gaps, of the high parts' distance g. Distance rounds g within 4 u g: the difference,
the two squares, their sum and the root, each once. Bowline takes a crossing only where the
distances it works out, d, part by more than its slack, at least twice both points' gaps and 2^-30
of the larger d; every number in it at least 0, its own roundings, a sum of eleven terms, within 11
u of it, which the factor two leaves room for. Then the rounded points' distances part the same
way, strictly: the exact way, which takes the first where they are as near, takes the same one
(choice lemma).

Hopped (source/crossings/bows.cpp) takes one Newton step from a guess on the two bulges'
quadratics, G_k(Y + D) = G_k(Y) + J_k D + 4 b_k |D|^2 exactly, each worked out of sure numbers: the
step d = -J^-1 G(Y) within reach r of its every coordinate, J the exact gradients. The crossings'
steps are the fixed points of T(D) = d - J^-1 (4 b1, 4 b2) |D|^2; with N bounding J^-1's rows and
p = N max |b_k| (pull), |T(D) - d| is at most 4 p |D|_2^2 <= 8 p m^2 for m = |D|_inf, and T moves
two steps apart by at most 16 p R times their distance within R of the guess. Where 64 p r <= 1, T
maps the ball of 2 r into itself and contracts it, so by Banach's theorem it holds one crossing,
within 32 p r^2 of d; and T contracts every ball of radius under 1 / (16 p), which so holds no
other (quadratic lemma, the theorem by hand). The point is nearer that crossing than the other,
both rounded, where twice its gap from the guess, 3 r and the rounding of every coordinate, at most
2^-50 of their sizes, fall short of 1 / (16 p) (reach lemma).

A rounding is lemma.py's. The lemmas are not bound to program.py's graphs. Each mutant drops a part
of the slack, and Z3 must refute it. The exit status is 1 when a lemma is not proved or a mutant
is.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from lemma import Lemma, run  # pylint: disable=wrong-import-position

REST = "(/ 1.0 1125899906842624.0)"  # 2^-50, as Bowline weighs a coordinate's high part.
WIDE = "(/ 1.0 1073741824.0)"  # 2^-30, as Bowline weighs the larger distance.


def coordinate():
    """A coordinate H + L within B of the exact X, |L| <= u |H|; R the double nearest X. Then
    |R - H| <= (1 + u) B + u (2 + u) |H|, under 2 (B + 2^-50 |H|)."""
    lemma = Lemma("coordinate", "H", "L", "B", "X", "R")
    lemma.fact("(<= (absr L) (* u (absr H)))")
    lemma.fact("(>= B 0.0)")
    lemma.fact("(<= (absr (- X (+ H L))) B)")
    lemma.fact("(<= (absr (- R X)) (* u (absr X)))")
    lemma.goal(f"(<= (absr (- R H)) (* 2.0 (+ B (* {REST} (absr H)))))")
    return lemma


def choice(wide=WIDE, gaps=True):
    """Two points, the rounded ones' distances t and s from the hit, within e and f of the high
    parts' distances g and k; the distances worked out d and c, within 4 u of g and k; the slack
    W, the code's sum rounded down by at most 11 u, at least 2 (e + f) and 2^-30 max(d, c) once
    rounded. Where d - c > W, t > s; where c - d > W, s > t."""
    lemma = Lemma("choice", "t", "s", "g", "k", "e", "f", "d", "c", "W", "M")
    for name in ("g", "k", "e", "f"):
        lemma.fact(f"(>= {name} 0.0)")
    lemma.fact("(<= (absr (- t g)) e)")
    lemma.fact("(<= (absr (- s k)) f)")
    lemma.fact("(<= (absr (- d g)) (* 4.0 u g))")
    lemma.fact("(<= (absr (- c k)) (* 4.0 u k))")
    lemma.fact("(= M (ite (>= d c) d c))")
    slack = f"(* {wide} M)"
    if gaps:
        slack = f"(+ (* 2.0 (+ e f)) {slack})"
    lemma.fact(f"(>= W (* (- 1.0 (* 11.0 u)) {slack}))")
    lemma.goal("(=> (> (- d c) W) (> t s))")
    lemma.goal("(=> (> (- c d) W) (> s t))")
    return lemma


def quadratic(bound="(* 32.0 p r r)"):
    """With 64 p r <= 1, d within r of nothing, |T(D) - d| <= 8 p m^2: a step within 2 r maps
    within 2 r, steps within 2 r move at most half as far as they lie apart, and the fixed point is
    within 32 p r^2 of d. Written in the infinity norm: m the step's largest coordinate, t the
    image's, the second goal the Lipschitz factor 16 p R at R = 2 r."""
    lemma = Lemma("quadratic", "p", "r", "m", "t", "e")
    lemma.fact("(>= p 0.0)")
    lemma.fact("(>= r 0.0)")
    lemma.fact("(<= (* 64.0 p r) 1.0)")
    lemma.fact("(>= m 0.0)")
    lemma.fact("(<= m (* 2.0 r))")
    lemma.fact("(<= e (* 8.0 p m m))")
    lemma.fact("(<= t (+ r e))")
    lemma.goal("(<= t (* 2.0 r))")
    lemma.goal("(<= (* 16.0 p 2.0 r) 0.5)")
    lemma.goal(f"(<= e {bound})")
    return lemma


def reach():
    """The point's distance a from the crossing found, at most its distance g from the guess and
    the crossing's sqrt(2) 2 r from that, the other's distance o at least F = 1 / (16 p) less g;
    each rounded point's distance moved by at most w. Where 2 g + 3 r + w < F, as Hopped asks with
    w its 2^-50 of every coordinate's size, the rounded crossing found is strictly nearer."""
    lemma = Lemma("reach", "g", "r", "w", "F", "a", "o", "A", "O")
    for name in ("g", "r", "w", "F"):
        lemma.fact(f"(>= {name} 0.0)")
    lemma.fact("(<= a (+ g (* 2.83 r)))")
    lemma.fact("(>= o (- F g))")
    lemma.fact("(<= (absr (- A a)) (* 0.5 w))")
    lemma.fact("(<= (absr (- O o)) (* 0.5 w))")
    lemma.fact("(< (+ (* 2.0 g) (* 3.0 r) w) F)")
    lemma.goal("(< A O)")
    return lemma


MUTANTS = (
    ("quadratic at 16 p r^2", lambda: quadratic("(* 16.0 p r r)")),
    ("choice without the larger distance's share", lambda: choice(wide="0.0")),
    ("choice without the points' gaps", lambda: choice(gaps=False)),
)

LEMMAS = (coordinate, choice, quadratic, reach)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, mutants=MUTANTS))
