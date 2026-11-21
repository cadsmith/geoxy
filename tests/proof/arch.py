#!/usr/bin/env python3
"""Machine-checks with Z3 the bounds source/swift/vault.cpp's kernel for a quadratic Bézier curve met
with a spline rests on: the curve's parabola as the implicit equation L1^2 - 4 L0 L2 of the
barycentric coordinates of its control triangle, its place affine in the point, its gradient
along the curve, how far apart two of its points lie by their places, the reach the equation on a
span must clear, a linear form's error worked out in doubles, and the sine at a crossing.

    python3 -I tests/proof/arch.py [z3]

A rounding is lemma.py's: |r - z| <= u |z| for a sum, a difference or a product clear of
underflow. Each lemma states the operations of vault.cpp as written, one fact a rounding; they are
not bound to program.py's graphs, so translation.py does not check them against the compiler's
IR. Each mutant drops a term the code keeps or changes a sign, and Z3 must refute it: so the
proofs need every part of the bounds they are given. The exit status is 1 when a lemma is not
proved or a mutant is.

The coefficients of the equation on a span are sums of C(n, i) C(n, j) (l1_i l1_j - 4 l0_i l2_j)
over C(2 n, k), each term two products and a difference and the binomial weight's product, as a
ring's |G_i . G_j| terms are: ring.py's term, running, finished (its rho^2 nothing, so its
product and difference exact) and wide hold them within the 16 u vault.cpp charges for every
degree to eight.

Left to the comments of vault.cpp, worked by hand: that |grad Q| over a span is no more than at the
hull's corners, each |L_k| being convex; the mean value theorem's bound on the one double-double
Newton step, as the ring kernel's; that the isolation's stretches hold every place where the span
comes within the reach; and the reasoning that the thorough way's tolerance decides nothing about
a crossing the kernel is sure of.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from lemma import Lemma, rational, run  # pylint: disable=wrong-import-position

GROW = 1.0 + 2.0 ** -40


def rounds(lemma, result, exact):
    """result is exact rounded to nearest: a sum, a difference or a product clear of underflow."""
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {exact})))")


def cross(one, other):
    """The cross product of two vectors given by their coordinates' names."""
    return f"(- (* {one[0]} {other[1]}) (* {one[1]} {other[0]}))"


def curve(lemma, place="s"):
    """The quadratic Bézier curve's point at a place less its middle control point, g = (1 - s)^2
    (-e2) + s^2 e0, for e2 = P1 - P0 and e0 = P2 - P1, the coordinates gx and gy."""
    lemma.fact(f"(= gx (+ (* (- 1.0 {place}) (- 1.0 {place}) (- ax)) (* {place} {place} bx)))")
    lemma.fact(f"(= gy (+ (* (- 1.0 {place}) (- 1.0 {place}) (- ay)) (* {place} {place} by)))")


def implicit(middle=True):
    """The barycentric coordinates about P1 at the curve's point: L0 = e0 x g = (1 - s)^2 D, L2 =
    e2 x g = s^2 D and L1 = D - e1 x g = 2 s (1 - s) D for e1 = e2 + e0 and D = e2 x e0; so the
    curve's points lie on L1^2 = 4 L0 L2."""
    lemma = Lemma("implicit", "s", "ax", "ay", "bx", "by", "gx", "gy", "d")
    curve(lemma)
    lemma.fact(f"(= d {cross(('ax', 'ay'), ('bx', 'by'))})")
    first = cross(("bx", "by"), ("gx", "gy"))
    last = cross(("ax", "ay"), ("gx", "gy"))
    line = cross(("(+ ax bx)", "(+ ay by)"), ("gx", "gy"))
    center = f"(- d {line})" if middle else f"(+ d {line})"
    lemma.goal(f"(= {first} (* (- 1.0 s) (- 1.0 s) d))")
    lemma.goal(f"(= {last} (* s s d))")
    lemma.goal(f"(= {center} (* 2.0 s (- 1.0 s) d))")
    lemma.goal(f"(= (* {center} {center}) (* 4.0 {first} {last}))")
    return lemma


def place(halved=True):
    """The curve's place at its point, affine in it: 2 D (s - 1/2) = (e2 - e0) x g."""
    lemma = Lemma("place", "s", "ax", "ay", "bx", "by", "gx", "gy", "d")
    curve(lemma)
    lemma.fact(f"(= d {cross(('ax', 'ay'), ('bx', 'by'))})")
    shift = "0.5" if halved else "0.0"
    lemma.goal(f"(= (* 2.0 d (- s {shift})) {cross(('(- ax bx)', '(- ay by)'), ('gx', 'gy'))})")
    return lemma


def gradient(halved=True):
    """Q's gradient at the curve's point: 2 L1 grad L1 - 4 (L0 grad L2 + L2 grad L0), grad L0 =
    (-e0_y, e0_x), grad L2 = (-e2_y, e2_x) and grad L1 = (e1_y, -e1_x), is 2 D times B'(s) turned
    a quarter clockwise, B'(s) = 2 ((1 - s) e2 + s e0): its length 2 |D| |B'|, and its dot with
    any vector v -2 D B' x v."""
    lemma = Lemma("gradient", "s", "ax", "ay", "bx", "by", "gx", "gy", "d", "vx", "vy")
    curve(lemma)
    lemma.fact(f"(= d {cross(('ax', 'ay'), ('bx', 'by'))})")
    first = cross(("bx", "by"), ("gx", "gy"))
    last = cross(("ax", "ay"), ("gx", "gy"))
    center = f"(- d {cross(('(+ ax bx)', '(+ ay by)'), ('gx', 'gy'))})"
    across = f"(- (* 2.0 {center} (+ ay by)) (* 4.0 (+ (* {first} (- ay)) (* {last} (- by)))))"
    upward = f"(- (* 2.0 {center} (- (+ ax bx))) (* 4.0 (+ (* {first} ax) (* {last} bx))))"
    scale = "2.0" if halved else "1.0"
    wayx = f"(* {scale} (+ (* (- 1.0 s) ax) (* s bx)))"
    wayy = f"(* {scale} (+ (* (- 1.0 s) ay) (* s by)))"
    lemma.goal(f"(= {across} (* 2.0 d {wayy}))")
    lemma.goal(f"(= {upward} (* -2.0 d {wayx}))")
    lemma.goal(f"(= (+ (* {across} vx) (* {upward} vy)) (* -2.0 d {cross((wayx, wayy), ('vx', 'vy'))}))")
    return lemma


def chord(halved=True):
    """Two points of a quadratic Bézier curve, a coordinate at a time: B(t) - B(s) = (t - s)
    B'((s + t) / 2), so that their distance is |t - s| times the speed halfway."""
    lemma = Lemma("chord", "s", "t", "p0", "p1", "p2")
    point = "(+ (* (- 1.0 {0}) (- 1.0 {0}) p0) (* 2.0 {0} (- 1.0 {0}) p1) (* {0} {0} p2))"
    middle = "(/ (+ s t) 2.0)" if halved else "s"
    speed = f"(* 2.0 (+ (* (- 1.0 {middle}) (- p1 p0)) (* {middle} (- p2 p1))))"
    lemma.goal(f"(= (- {point.format('t')} {point.format('s')}) (* (- t s) {speed}))")
    return lemma


def least(crossed=True):
    """The least speed: B'(m) / 2 = (1 - m) e2 + m e0 and e0 - e2 have the cross product D,
    whatever m, so by Lagrange's identity |B'(m) / 2| |e0 - e2| >= |D|."""
    lemma = Lemma("least", "m", "ax", "ay", "bx", "by", "d")
    lemma.fact(f"(= d {cross(('ax', 'ay'), ('bx', 'by'))})")
    wayx = "(+ (* (- 1.0 m) ax) (* m bx))"
    wayy = "(+ (* (- 1.0 m) ay) (* m by))"
    lemma.goal(f"(= {cross((wayx, wayy), ('(- bx ax)', '(- by ay)'))} d)")
    size = "(* d d)" if crossed else "(* 2.0 d d)"
    lemma.goal(f"(>= (* (+ (* {wayx} {wayx}) (* {wayy} {wayy})) "
               f"(+ (* (- bx ax) (- bx ax)) (* (- by ay) (- by ay)))) {size})")
    return lemma


def scaled(sized=True):
    """A product bounded by its factors' bounds: |a b| <= |a| B where |b| <= B, and |b c| <= B C
    where also |c| <= C."""
    lemma = Lemma("scaled", "a", "b", "c", "B", "C")
    lemma.fact("(<= (absr b) B)")
    lemma.fact("(<= (absr c) C)")
    lemma.goal("(<= (absr (* a b)) (* (absr a) B))")
    bound = "(* B C)" if sized else "(* 0.5 B C)"
    lemma.goal(f"(<= (absr (* b c)) {bound})")
    return lemma


def expand(paired=True):
    """The parabola's equation at a point q = p + Delta, each L_k(q) = l_k + d_k: Q(q) = 0 is
    Q(p) + 2 l1 d1 + d1^2 - 4 (l0 d2 + l2 d0 + d0 d2) = 0, the products named."""
    lemma = Lemma("expand", "l0", "l1", "l2", "d0", "d1", "d2", "level", "p1", "p02", "p20",
                  "s1", "s02")
    lemma.fact("(= (* (+ l1 d1) (+ l1 d1)) (* 4.0 (+ l0 d0) (+ l2 d2)))")
    lemma.fact("(= level (- (* l1 l1) (* 4.0 l0 l2)))")
    lemma.fact("(= p1 (* l1 d1))")
    lemma.fact("(= p02 (* l0 d2))")
    lemma.fact("(= p20 (* l2 d0))")
    lemma.fact("(= s1 (* d1 d1))")
    lemma.fact("(= s02 (* d0 d2))")
    twice = "2.0" if paired else "1.0"
    lemma.goal(f"(= (+ level (* {twice} p1) s1) (* 4.0 (+ p02 p20 s02)))")
    return lemma


def reach(curved=True):
    """A point p at distance no more than r from a point q of the parabola: each L_k is affine, so
    L_k(q) = l_k + d_k, l_k = L_k(p) and |d_k| <= |grad L_k| r = N_k r; by expand and scaled's
    bounds on the products, b1 = |l1| N1 r, b02 = |l0| N2 r, b20 = |l2| N0 r, c1 = N1^2 r^2 and
    c02 = N0 N2 r^2, |Q(p)| is no more than (2 |l1| N1 + 4 |l0| N2 + 4 |l2| N0) r + (N1^2 + 4 N0
    N2) r^2. Where |Q(p)| is above that, p lies farther than r from the parabola."""
    lemma = Lemma("reach", "level", "p1", "p02", "p20", "s1", "s02", "b1", "b02", "b20", "c1",
                  "c02")
    lemma.fact("(= (+ level (* 2.0 p1) s1) (* 4.0 (+ p02 p20 s02)))")
    for name, bound in (("p1", "b1"), ("p02", "b02"), ("p20", "b20"), ("s1", "c1"),
                        ("s02", "c02")):
        lemma.fact(f"(<= (absr {name}) {bound})")
    square = "(+ c1 (* 4.0 c02))" if curved else "0.0"
    lemma.goal(f"(<= (absr level) (+ (* 2.0 b1) (* 4.0 b02) (* 4.0 b20) {square}))")
    return lemma


def rounding(weighed=True):
    """A linear form at a moved control point, as Springing works it out from doubles: c = rnd(rnd(ex
    gy) - rnd(ey gx)) and l = rnd(rnd(a w) - c), the exact products P = ex gy, R = ey gx and A = a
    w. Then c is within 4 u Grow (|P| + |R|) of P - R, and l within 4 u Grow (|P| + |R| + |A|) of
    A - (P - R)."""
    lemma = Lemma("rounding", "P", "R", "A", "p", "q", "c", "m", "l")
    rounds(lemma, "p", "P")
    rounds(lemma, "q", "R")
    rounds(lemma, "c", "(- p q)")
    rounds(lemma, "m", "A")
    rounds(lemma, "l", "(- m c)")
    grow = rational(GROW)
    lemma.goal(f"(<= (absr (- c (- P R))) (* 4.0 u {grow} (+ (absr P) (absr R))))")
    weight = "(absr A)" if weighed else "0.0"
    lemma.goal(f"(<= (absr (- l (- A (- P R)))) (* 4.0 u {grow} (+ (absr P) (absr R) {weight})))")
    return lemma


def apart():
    """The inputs' errors split: e x g - E x G = e x (g - G) + (e - E) x G, and a w - D W = (a - D) w
    + D (w - W)."""
    lemma = Lemma("apart", "ex", "ey", "Ex", "Ey", "gx", "gy", "Gx", "Gy", "a", "A", "w", "W")
    lemma.goal("(= (- (- (* ex gy) (* ey gx)) (- (* Ex Gy) (* Ey Gx))) "
               "(- (+ (* ex (- gy Gy)) (* (- ex Ex) Gy)) (+ (* ey (- gx Gx)) (* (- ey Ey) Gx))))")
    lemma.goal("(= (- (* a w) (* A W)) (+ (* (- a A) w) (* A (- w W))))")
    return lemma


def carried(sloped=True):
    """What the inputs' errors carry into e x G and D W, by apart's split: the rounded side e within
    its slip h of the exact E, the moved point g within (ox, oy) of the true G, the area a within
    its drift f of D and the weight w within the heft k of W; each product bounded by scaled, with
    |G| no more than |g| + o and |D| no more than |a| + f."""
    lemma = Lemma("carried", "t1", "t2", "t3", "t4", "t5", "t6", "b1", "b2", "b3", "b4", "b5",
                  "b6", "split", "weighed")
    lemma.fact("(= split (- (+ t1 t2) (+ t3 t4)))")
    lemma.fact("(= weighed (+ t5 t6))")
    for index in range(1, 7):
        lemma.fact(f"(<= (absr t{index}) b{index})")
    across = "(+ b1 b2 b3 b4)" if sloped else "(+ b1 b2 b3)"
    lemma.goal(f"(<= (absr split) {across})")
    lemma.goal("(<= (absr weighed) (+ b5 b6))")
    return lemma


def sine(heavy=True):
    """The sine at a crossing: H' = W^2 (grad Q . p') where Q is nothing, grad Q . p' = -2 D B' x
    p' (gradient), |B'| no more than Fast, |p'| no more than the speed, W no more than W_max; where
    |H'| >= g, the sine |B' x p'| / (|B'| |p'|) is at least g / (W_max^2 2 |D| Fast speed)."""
    lemma = Lemma("sine", "crossed", "w", "most", "slope", "g", "d", "way", "fast", "pace",
                  "speed")
    for name in ("way", "fast", "pace", "speed"):
        lemma.fact(f"(> {name} 0.0)")
    lemma.fact("(> w 0.0)")
    lemma.fact("(<= w most)")
    lemma.fact("(not (= d 0.0))")
    lemma.fact("(= slope (* w w -2.0 d crossed))")
    lemma.fact("(>= (absr slope) g)")
    lemma.fact("(>= g 0.0)")
    lemma.fact("(<= way fast)")
    lemma.fact("(<= pace speed)")
    weight = "(* most most)" if heavy else "1.0"
    lemma.goal(f"(>= (* (absr crossed) {weight} 2.0 (absr d) fast speed) (* g way pace))")
    return lemma


MUTANTS = (("implicit with L1 = D + e1 x g", lambda: implicit(middle=False)),
           ("place without the half", lambda: place(halved=False)),
           ("gradient of B' / 2", lambda: gradient(halved=False)),
           ("chord at the first place", lambda: chord(halved=False)),
           ("least at twice D^2", lambda: least(crossed=False)),
           ("scaled at half", lambda: scaled(sized=False)),
           ("expand with l1 d1 once", lambda: expand(paired=False)),
           ("reach without the second order", lambda: reach(curved=False)),
           ("rounding without the weighted area", lambda: rounding(weighed=False)),
           ("carried without a product", lambda: carried(sloped=False)),
           ("sine without the weight", lambda: sine(heavy=False)))

LEMMAS = (implicit, place, gradient, chord, least, scaled, expand, reach, rounding, apart,
          carried, sine)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, mutants=MUTANTS))
