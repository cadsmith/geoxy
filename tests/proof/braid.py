#!/usr/bin/env python3
"""Machine-checks the bounds source/swift/braid.cpp's kernel for a cubic Bézier curve met with a
spline rests on: the curve's implicit cubic as the determinant of its Bézout matrix about the first
control point, the matrix's kernel and the place of a point of the cubic as a quotient of two linear
forms, the determinant's minors and partial derivatives, how far the cubic's equation may move
between a point and one within r of it, a coefficient's term worked out in doubles and its 16 u
for every degree to five, and the sine at a crossing.

    python3 -I tests/proof/braid.py [z3]

The identities are checked in exact rationals, as polynomials in the control points' coordinates
less the first's and the curve's place, every coefficient of their difference nothing. The bounds
are Z3's: a rounding is lemma.py's, |r - z| <= u |z| for a sum, a difference or a product clear of
underflow. Each lemma states the operations of braid.cpp as written, one fact a rounding; they are
not bound to program.py's graphs, so translation.py does not check them against the compiler's
IR. Each mutant drops a term the code keeps, and Z3 must refute it: so the proofs need every part
of the bounds they are given. The exit status is 1 when a check fails, a lemma is not proved or a
mutant is.

The entries on a span are worked out as the parabola's lines are, arch.py's rounding, apart and
carried; a minor's coefficients as a ring's, ring.py's term, running and wide, its finished's
quotient left out; a product of two of them with ring.py's product. Left to the comments of
braid.cpp, worked by hand: that the equation on a span is W^3 f(N / W), the determinant of the
entries at the homogeneous point, each homogeneous of degree one; that an entry at a point of the
span is a convex mix of its coefficients over the weights; that the quotient's bound holds every
place of the cubic within the near reach of a crossing; the mean value theorem's bound on the one
double-double Newton step, as the ring kernel's; that the isolation's stretches hold every place
where the span comes within the reach; and the reasoning that the thorough way's tolerance decides
nothing about a crossing the kernel is sure of.
"""

import pathlib
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

# pylint: disable-next=wrong-import-position
from expand import Fraction, U  # pylint: disable=wrong-import-position
# pylint: disable-next=wrong-import-position
from lemma import Lemma, rational, run  # pylint: disable=wrong-import-position

GROW = 1.0 + 2.0 ** -40
TIGHTEST = 5
NAMES = ("ax", "ay", "bx", "by", "cx", "cy", "s")


def rounds(lemma, result, exact):
    """result is exact rounded to nearest: a sum, a difference or a product clear of underflow."""
    lemma.fact(f"(<= (absr (- {result} {exact})) (* u (absr {exact})))")


class Poly:
    """A polynomial in NAMES with exact rational coefficients, by its monomials' exponents."""

    def __init__(self, terms=None):
        self.terms = {key: value for key, value in (terms or {}).items() if value != 0}

    @staticmethod
    def constant(value):
        """A number as a polynomial."""
        return Poly({(0,) * len(NAMES): Fraction(value)})

    @staticmethod
    def variable(name):
        """One of NAMES as a polynomial."""
        powers = [0] * len(NAMES)
        powers[NAMES.index(name)] = 1
        return Poly({tuple(powers): Fraction(1)})

    def __add__(self, other):
        other = other if isinstance(other, Poly) else Poly.constant(other)
        terms = dict(self.terms)
        for key, value in other.terms.items():
            terms[key] = terms.get(key, 0) + value
        return Poly(terms)

    __radd__ = __add__

    def __neg__(self):
        return Poly({key: -value for key, value in self.terms.items()})

    def __sub__(self, other):
        return self + (-(other if isinstance(other, Poly) else Poly.constant(other)))

    def __rsub__(self, other):
        return Poly.constant(other) - self

    def __mul__(self, other):
        other = other if isinstance(other, Poly) else Poly.constant(other)
        terms = {}
        for one, first in self.terms.items():
            for two, second in other.terms.items():
                key = tuple(left + right for left, right in zip(one, two))
                terms[key] = terms.get(key, 0) + (first * second)
        return Poly(terms)

    __rmul__ = __mul__

    def nothing(self):
        """Whether every coefficient is nothing."""
        return not self.terms


def cross(one, other):
    """The cross product of two vectors of polynomials."""
    return (one[0] * other[1]) - (one[1] * other[0])


def curve():
    """The cubic's control points less the first, a1, a2 and a3, and its point at s less P0."""
    first = (Poly.variable("ax"), Poly.variable("ay"))
    second = (Poly.variable("bx"), Poly.variable("by"))
    third = (Poly.variable("cx"), Poly.variable("cy"))
    place = Poly.variable("s")
    rest = 1 - place
    weights = (3 * place * rest * rest, 3 * place * place * rest, place * place * place)
    point = tuple(sum((weight * pole[axis] for weight, pole in zip(weights, (first, second, third))),
                      Poly()) for axis in range(2))
    return (first, second, third), point


def entries(poles, gap):
    """The Bézout matrix's entries a, b, c, d, e and g at a gap from P0, each K - S x g, as braid.cpp's
    Braiding has them."""
    first, second, third = poles

    def scaled(factor, side):
        return (factor * side[0], factor * side[1])

    def apart(head, tail):
        return (head[0] - tail[0], head[1] - tail[1])

    nine = scaled(9, apart(second, first))
    sides = (scaled(3, apart(third, second)), scaled(3, apart(third, first)), third,
             (third[0] + nine[0], third[1] + nine[1]), scaled(3, second), scaled(3, first))
    areas = (3 * cross(third, second), 3 * cross(third, first), Poly(), 9 * cross(second, first),
             Poly(), Poly())
    return [area - cross(side, gap) for area, side in zip(areas, sides)]


def determinant(cells):
    """a (d g - e^2) - b (b g - c e) + c (b e - c d), the symmetric matrix's determinant."""
    a, b, c, d, e, g = cells
    return (a * ((d * g) - (e * e))) - (b * ((b * g) - (c * e))) + (c * ((b * e) - (c * d)))


def implicit():
    """The cubic's points lie on f = 0, at every place: f(B(s) - P0) is nothing as a polynomial,
    and (s^2, s (1 - s), (1 - s)^2) is the matrix's kernel there."""
    poles, point = curve()
    cells = entries(poles, point)
    a, b, c, d, e, g = cells
    place = Poly.variable("s")
    kernel = (place * place, place * (1 - place), (1 - place) * (1 - place))
    rows = ((a, b, c), (b, d, e), (c, e, g))
    held = determinant(cells).nothing() and all(
        sum((cell * part for cell, part in zip(row, kernel)), Poly()).nothing() for row in rows)
    return held, "f(B(s)) and M(B(s)) v(s) nothing"


def inverse():
    """The place of a point of the cubic: N = 3 A21 a + A13 b + A32 c and D = 3 A21 b + A13 d + A32
    e give s (D - N) = D at B(s), lambda = (3 A21, A13, A32) killing the column of c, e and g."""
    poles, point = curve()
    first, second, third = poles
    a, b, c, d, e, g = entries(poles, point)
    lambdas = (3 * cross(second, first), -cross(third, first), cross(third, second))
    top = (lambdas[0] * a) + (lambdas[1] * b) + (lambdas[2] * c)
    bottom = (lambdas[0] * b) + (lambdas[1] * d) + (lambdas[2] * e)
    killed = (lambdas[0] * c) + (lambdas[1] * e) + (lambdas[2] * g)
    place = Poly.variable("s")
    held = killed.nothing() and ((place * (bottom - top)) - bottom).nothing()
    return held, "s (D - N) = D at B(s), the third column killed"


def partials():
    """The determinant as its monomials a d g - a e^2 - b^2 g + 2 b c e - c^2 d, and its partial
    derivatives as Pitched bounds them, each the difference of f at an entry moved by one less f,
    the parts of the second order and on taken off."""
    names = ("ax", "ay", "bx", "by", "cx", "cy")
    cells = [Poly.variable(name) for name in names]
    a, b, c, d, e, g = cells
    monomials = (a * d * g) - (a * e * e) - (b * b * g) + (2 * b * c * e) - (c * c * d)
    held = (determinant(cells) - monomials).nothing()
    derivatives = ((d * g) - (e * e), 2 * ((c * e) - (b * g)), 2 * ((b * e) - (c * d)),
                   (a * g) - (c * c), 2 * ((b * c) - (a * e)), (a * d) - (b * b))
    step = Poly.variable("s")
    for index, derivative in enumerate(derivatives):
        moved = list(cells)
        moved[index] = moved[index] + step
        change = determinant(moved) - determinant(cells)
        linear = Poly({key: value for key, value in change.terms.items()
                       if key[NAMES.index("s")] == 1})
        held = held and (linear - (derivative * step)).nothing()
    return held, "the monomials and the six partial derivatives"


def expanded(paired=True):
    """A product of three entries at a point q = p + Delta, each moved: (x + p)(y + q)(z + w) - x y
    z is the sum of the seven products of moves and values, the first order's three, the
    second's three and the third's one."""
    lemma = Lemma("expanded", "x", "y", "z", "p", "q", "w")
    second = "(+ (* p q z) (* p y w) (* x q w))" if paired else "(+ (* p q z) (* p y w))"
    lemma.goal("(= (- (* (+ x p) (+ y q) (+ z w)) (* x y z)) "
               f"(+ (* p y z) (* x q z) (* x y w) {second} (* p q w)))")
    return lemma


def scaled(sized=True):
    """A product bounded by its factors' bounds: |a b| <= A B where |a| <= A and |b| <= B, and so a
    product of three by applying it twice."""
    lemma = Lemma("scaled", "a", "b", "A", "B")
    lemma.fact("(<= (absr a) A)")
    lemma.fact("(<= (absr b) B)")
    bound = "(* A B)" if sized else "(* 0.5 A B)"
    lemma.goal(f"(<= (absr (* a b)) {bound})")
    return lemma


def monomial(cubed=True):
    """A product of three entries at a point q = p + Delta, each within its side's length times r
    of its value at p, |x| <= X and |p| <= Nx r and so on: by expanded and scaled, each of the seven
    products t_i within its factors' bounds b_i, the change no more than their sum, which gathered
    puts as U1 r + U2 r^2 + U3 r^3. f's five monomials so bounded, |f(p)| = |f(p) - f(q)| is no
    more than their sum, a point of the cubic q within r of p."""
    names = [f"t{index}" for index in range(7)]
    bounds = [f"b{index}" for index in range(7)]
    lemma = Lemma("monomial", "change", *names, *bounds)
    lemma.fact(f"(= change (+ {' '.join(names)}))")
    for name, bound in zip(names, bounds):
        lemma.fact(f"(<= (absr {name}) {bound})")
    kept = bounds if cubed else bounds[:-1]
    lemma.goal(f"(<= (absr change) (+ {' '.join(kept)}))")
    return lemma


def gathered(cubed=True):
    """The seven products' bounds gathered by the order of r: Nx r Y Z + X Ny r Z + X Y Nz r + Nx r
    Ny r Z + Nx r Y Nz r + X Ny r Nz r + Nx r Ny r Nz r = (Nx Y Z + X Ny Z + X Y Nz) r + (Nx Ny Z +
    Nx Y Nz + X Ny Nz) r^2 + Nx Ny Nz r^3, as Teased sums them."""
    lemma = Lemma("gathered", "X", "Y", "Z", "nx", "ny", "nz", "r")
    bounds = ("(* nx r Y Z)", "(* X ny r Z)", "(* X Y nz r)", "(* nx r ny r Z)",
              "(* nx r Y nz r)", "(* X ny r nz r)", "(* nx r ny r nz r)")
    first = "(+ (* nx Y Z) (* X ny Z) (* X Y nz))"
    second = "(+ (* nx ny Z) (* nx Y nz) (* X ny nz))"
    third = "(* nx ny nz r r r)" if cubed else "0.0"
    lemma.goal(f"(= (+ {' '.join(bounds)}) (+ (* {first} r) (* {second} r r) {third}))")
    return lemma


def summed(share=3.0):
    """A coefficient's term on a span before its binomial weight, as Springing works it out: three
    products, a difference and a sum, s = rnd(rnd(p0 - p1) + p2), p_i = rnd(P_i) the rounded
    products of the exact P_i = x_i y_i. It is within 3 u Grow of the sizes |P0| + |P1| + |P2| of
    the exact P0 - P1 + P2."""
    lemma = Lemma("summed", "P0", "P1", "P2", "p0", "p1", "p2", "d", "s")
    for index in range(3):
        rounds(lemma, f"p{index}", f"P{index}")
    rounds(lemma, "d", "(- p0 p1)")
    rounds(lemma, "s", "(+ d p2)")
    size = "(+ (absr P0) (absr P1) (absr P2))"
    lemma.goal(f"(<= (absr (- s (+ (- P0 P1) P2))) (* {share} u {rational(GROW)} {size}))")
    return lemma


def term(share=4.0):
    """The term weighed by its binomial weight m >= 0, t = rnd(m s), s within summed's 3 u Grow of
    the sizes M of its exact S: within 4 u Grow of m M of m S."""
    lemma = Lemma("term", "s", "S", "M", "m", "t")
    lemma.fact("(>= m 0.0)")
    lemma.fact("(>= M 0.0)")
    lemma.fact("(<= (absr S) M)")
    lemma.fact(f"(<= (absr (- s S)) (* 3.0 u {rational(GROW)} M))")
    rounds(lemma, "t", "(* m s)")
    lemma.goal(f"(<= (absr (- t (* m S))) (* {share} u {rational(GROW)} m M))")
    return lemma


def wide():
    """A coefficient of the equation on a span of degree n to Tightest sums no more than n + 1
    terms, each within 4 u Grow of its size m_i (term), their sizes M in all. By ring.py's running,
    n times, the sum is within D of the exact one, D (1 - n u) no more than 4 u Grow M + n u M (1 +
    4 u Grow); the quotient by the binomial C(3 n, k) adds u of the sum's size, no more than M + D.
    All of it within the 16 u Springing charges, for every n to Tightest. Checked in exact
    rationals."""
    grow = Fraction(GROW)
    held = True
    for order in range(1, TIGHTEST + 1):
        share = 4 * U * grow
        error = (share + order * U * (1 + share)) / (1 - order * U)
        held = held and error + U * (1 + error) <= 16 * U
    return held, f"every degree to {TIGHTEST}"


def sine(heavy=True):
    """The sine at a crossing: H' = W^3 (grad f . p') where f is nothing, |grad f| no more than G,
    |p'| no more than the speed, W no more than W_max; where |H'| >= g, the sine |grad f . p'| /
    (|grad f| |p'|) is at least g / (W_max^3 G speed)."""
    lemma = Lemma("sine", "dot", "w", "most", "slope", "g", "normal", "graded", "speed", "way")
    for name in ("normal", "graded", "speed", "way"):
        lemma.fact(f"(> {name} 0.0)")
    lemma.fact("(> w 0.0)")
    lemma.fact("(<= w most)")
    lemma.fact("(= slope (* w w w dot))")
    lemma.fact("(>= (absr slope) g)")
    lemma.fact("(>= g 0.0)")
    lemma.fact("(<= normal graded)")
    lemma.fact("(<= way speed)")
    weight = "(* most most most)" if heavy else "(* most most)"
    lemma.goal(f"(>= (* (absr dot) {weight} graded speed) (* g normal way))")
    return lemma


MUTANTS = (("expanded without a second order product", lambda: expanded(paired=False)),
           ("scaled at half", lambda: scaled(sized=False)),
           ("monomial without the third order's product", lambda: monomial(cubed=False)),
           ("gathered without the third order", lambda: gathered(cubed=False)),
           ("summed at 2 u", lambda: summed(share=2.0)),
           ("term at 3 u", lambda: term(share=3.0)),
           ("sine at W_max^2", lambda: sine(heavy=False)))

LEMMAS = (expanded, scaled, monomial, gathered, summed, term, sine)

if __name__ == "__main__":
    sys.exit(run(LEMMAS, checks=(implicit, inverse, partials, wide), mutants=MUTANTS))
