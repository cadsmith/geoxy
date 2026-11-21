#!/usr/bin/env python3
"""Works out, at 60 significant digits, the lengths and areas of the splines tests/geoxy/splines.hpp
holds, and writes its rows.

    python3 tools/reference/measures.py [first] > rows.txt

The splines are drawn from a fixed seed: circles and ellipses written as NURBS of integer weights,
and arcs of them with their weights a power of two apart; weights 1 to 1e16 apart, and 1e100 and
1e-100 beside 1; random rational splines of degrees 1 to 7 and polynomial ones of degrees 4 to 12,
clamped, uneven and unclamped; the image of a cusped cubic under a map of the projective plane; and
rational curves that run back along a line.

Each span is turned into its rational Bezier form in exact rational arithmetic, the blossoms of the
doubles given worked out as fractions. On a span, with X, Y and W its homogeneous coordinates in its
own parameter t from 0 to 1, the speed is |W X' - W' X, W Y' - W' Y| / W^2 and the area's rate
(X - x W)(Y - y W)' - (Y - y W)(X - x W)' over 2 W^2, x and y the curve's start, every numerator a
polynomial with exact coefficients. They are integrated by mpmath at 60 digits and as many more as
the weights lie powers of ten apart, so that a rush as near 1 as a heavy weight puts it is told
apart from 1; cut at the real part of every root of the speed's numerator and of W near the span,
and at points graded geometrically toward each; by tanh-sinh and by Gauss-Legendre, which must
agree to 38 digits. The area's scale, the integral of its rate's length, cut where the rate is 0
too, is wanted as a double alone, by tanh-sinh.

Each row is the spline, its length and its area to 30 digits, and the scale as a double, from the
row first on, 0 by default. Run with a Python that has mpmath.
"""

import random
import sys
from decimal import Decimal, localcontext
from fractions import Fraction

import mpmath
from mpmath import mp, mpf

mp.dps = 60


def spans(degree, knots, poles, weights):
    """Each span that is not empty, as its rational Bezier control points, homogeneous, exact."""
    knots = [Fraction(knot) for knot in knots]
    weights = [Fraction(weight) for weight in weights] if weights else [Fraction(1)] * len(poles)
    lifted = [(Fraction(x) * w, Fraction(y) * w, w) for (x, y), w in zip(poles, weights)]
    out = []
    for span in range(degree, len(poles)):
        if not knots[span] < knots[span + 1]:
            continue
        controls = []
        for ends in range(degree + 1):
            values = [knots[span]] * (degree - ends) + [knots[span + 1]] * ends
            cells = [lifted[span - degree + index] for index in range(degree + 1)]
            for level in range(1, degree + 1):
                for index in range(degree, level - 1, -1):
                    start = knots[span - degree + index]
                    end = knots[span + 1 + index - level]
                    share = (values[level - 1] - start) / (end - start)
                    cells[index] = tuple((1 - share) * cells[index - 1][c] + share * cells[index][c]
                                         for c in range(3))
            controls.append(cells[degree])
        out.append(controls)
    return out


def binomial(n, k):
    """n choose k."""
    out = 1
    for index in range(k):
        out = out * (n - index) // (index + 1)
    return out


def powers(bernstein):
    """A polynomial's Bernstein coefficients over 0 to 1 as its coefficients in powers of t."""
    degree = len(bernstein) - 1
    out = [Fraction(0)] * (degree + 1)
    for index, value in enumerate(bernstein):
        for step in range(degree - index + 1):
            out[index + step] += (value * binomial(degree, index) * binomial(degree - index, step)
                                  * (-1) ** step)
    return out


def product(lhs, rhs):
    """Two polynomials' product."""
    out = [Fraction(0)] * (len(lhs) + len(rhs) - 1)
    for i, a in enumerate(lhs):
        for j, b in enumerate(rhs):
            out[i + j] += a * b
    return out


def difference(lhs, rhs):
    """Two polynomials' difference."""
    size = max(len(lhs), len(rhs))
    lhs = lhs + [Fraction(0)] * (size - len(lhs))
    rhs = rhs + [Fraction(0)] * (size - len(rhs))
    return [a - b for a, b in zip(lhs, rhs)]


def slope(values):
    """A polynomial's derivative."""
    return [values[k] * k for k in range(1, len(values))] or [Fraction(0)]


def trimmed(values):
    """A polynomial without the zeros at its top."""
    values = list(values)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return values


def lifted(values):
    """A polynomial's exact coefficients as mpf."""
    return [mpf(value.numerator) / value.denominator for value in values]


def horner(values, along):
    """A polynomial at a point."""
    out = mpf(0)
    for value in reversed(values):
        out = out * along + value
    return out


def roots(values):
    """A polynomial's roots, its coefficients in powers, the constant first."""
    values = list(values)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    if len(values) <= 1:
        return []
    try:
        return mpmath.polyroots(list(reversed(values)), maxsteps=400, extraprec=400)
    except mpmath.libmp.NoConvergence:
        return mpmath.polyroots(list(reversed(values)), maxsteps=4000, extraprec=1500)


def cuts(singular):
    """Where to cut 0 to 1: each root's real part inside it, and points graded toward each root
    near it by 4 at a time, from as near as the working precision tells apart from 0 and 1."""
    points = {mpf(0), mpf(1)}
    for root in singular:
        real = mpmath.re(root)
        near = min(max(real, mpf(0)), mpf(1))
        reach = abs(root - near)
        if 0 < real < 1:
            points.add(real)
        if reach > mpf(1) / 4:
            continue
        step = max(reach, mpf(10) ** (10 - mp.dps))
        while step < 1:
            for point in (near - step, near + step):
                if 0 < point < 1:
                    points.add(point)
            step *= 4
    return sorted(points)


def integrated(function, points):
    """An integral by tanh-sinh and by Gauss-Legendre, which must agree."""
    first = mp.quad(function, points, method="tanh-sinh")
    second = mp.quad(function, points, method="gauss-legendre")
    if abs(first - second) > mpf(10) ** -38 * (abs(first) + mpf(10) ** -300):
        raise ValueError(f"the rules disagree: {first} {second}")
    return first


def measured(degree, knots, poles, weights):
    """A spline's length, area and the area's scale, worked out with as many more digits as its
    weights lie powers of ten apart, so that a place as near 1 as a heavy weight puts the curve's
    rush is told apart from 1."""
    spread = int(mpmath.log10(mpf(max(weights))) - mpmath.log10(mpf(min(weights)))) if weights else 0
    with mp.workdps(60 + spread):
        return worked(degree, knots, poles, weights)


def worked(degree, knots, poles, weights):
    """A spline's length, area and the area's scale, at the working precision."""
    pieces = spans(degree, knots, poles, weights)
    x0, y0, w0 = pieces[0][0]
    origin = (x0 / w0, y0 / w0)
    length = area = scale = mpf(0)
    for controls in pieces:
        x = powers([c[0] for c in controls])
        y = powers([c[1] for c in controls])
        w = powers([c[2] for c in controls])
        across = trimmed(difference(product(w, slope(x)), product(slope(w), x)))
        rising = trimmed(difference(product(w, slope(y)), product(slope(w), y)))
        moved = difference(x, [c * origin[0] for c in w])
        raised = difference(y, [c * origin[1] for c in w])
        turning = trimmed(difference(product(moved, slope(raised)), product(raised, slope(moved))))
        across, rising, turning, weight = (lifted(p) for p in (across, rising, turning, trimmed(w)))
        size = max(len(across), len(rising))
        speeds = roots([mpmath.mpc(across[k] if k < len(across) else 0,
                                   rising[k] if k < len(rising) else 0) for k in range(size)])
        poles_ = roots(weight)
        length += integrated(lambda t: mpmath.hypot(horner(across, t), horner(rising, t))
                             / horner(weight, t) ** 2, cuts(list(speeds) + list(poles_)))
        rate = lambda t: horner(turning, t) / horner(weight, t) ** 2  # noqa: E731
        area += integrated(rate, cuts(poles_)) / 2
        zeros = [r for r in roots(turning) if abs(mpmath.im(r)) < mpf(10) ** -30]
        scale += mp.quad(lambda t: abs(rate(t)), cuts(list(poles_) + zeros)) / 2
    return length, area, scale


def drawn():
    """The splines, with their categories."""
    random.seed(20261114)
    out = []

    def add(degree, knots, poles, weights):
        out.append((degree, [float(k) for k in knots], [(float(x), float(y)) for x, y in poles],
                    [float(w) for w in weights]))

    def clamped(degree, count, uneven):
        inner = count - degree - 1
        cuts_ = (sorted(random.uniform(0, 1) for _ in range(inner)) if uneven
                 else [(i + 1) / (inner + 1) for i in range(inner)])
        return [0.0] * (degree + 1) + cuts_ + [1.0] * (degree + 1)

    circle = [(9, -12), (25, 0), (9, 12), (0, 18.75), (-9, 12), (-25, 0), (-9, -12), (0, -18.75),
              (9, -12)]
    rim = [25, 3, 1, 4, 25, 3, 1, 4, 25]
    ring = [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]
    add(2, ring, circle, rim)
    for shift, size in ((1e6, 1.0), (0.0, 1e-6), (-3e8, 2.0**-30), (0.0, 1e100), (5.0, 2.0**-500)):
        add(2, ring, [(x * size + shift, y * size - shift) for x, y in circle], rim)
    add(2, ring, [(2 * x, y) for x, y in circle], rim)
    add(2, ring, [(x + 0.5 * y, y) for x, y in circle], rim)
    arc = [(9, -12), (25, 0), (9, 12)]
    for power in range(0, 18, 3):
        rho = 2.0**power
        add(2, [0, 0, 0, 1, 1, 1], arc, [25, 3 * rho, rho * rho])
        add(2, [0, 0, 0, 1, 1, 1], arc, [25 * rho * rho, 3 * rho, 1])
    add(2, [0, 0, 0, 1, 1, 1], arc, [25, 3 * 2.0**100, 2.0**200])
    for power in range(0, 17, 3):
        w = 10.0**power
        add(2, [0, 0, 0, 1, 1, 1], [(0, 0), (1, 1), (2, 0)], [1, w, 1])
        add(2, [0, 0, 0, 1, 1, 1], [(0, 0), (1, 1), (2, 0)], [w, 1, 1])
        add(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0, 0), (1, 2), (3, 2), (4, 0)], [1, w, w, 1])
        add(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0, 0), (1, 2), (3, 2), (4, 0)], [1, w, 1 / w, 1])
        add(3, [0, 0, 0, 0, 0.3, 0.6, 1, 1, 1, 1],
            [(0, 0), (1, 2), (2, -1), (4, 3), (5, 0), (6, 2)], [1, w, 1, 1 / w, w, 1])
    for w in (1e100, 1e-100):
        add(2, [0, 0, 0, 1, 1, 1], [(0, 0), (1, 1), (2, 0)], [1, w, 1])
        add(3, [0, 0, 0, 0, 1, 1, 1, 1], [(0, 0), (1, 2), (3, 2), (4, 0)], [1, w, 1 / w, 1])
        add(3, [0, 0, 0, 0, 0.5, 1, 1, 1, 1], [(0, 0), (1, 2), (3, 2), (4, 0), (5, 1)],
            [w, 1, 1, 1, 1 / w])
    for degree in range(1, 8):
        for trial in range(4):
            count = degree + 1 + random.randint(0, 4)
            knots = ([float(k) for k in range(count + degree + 1)] if trial == 3
                     else clamped(degree, count, trial == 1))
            poles = [(random.uniform(-5, 5), random.uniform(-5, 5)) for _ in range(count)]
            weights = ([random.uniform(0.1, 10) for _ in range(count)] if trial < 2
                       else [10.0 ** random.uniform(-8, 8) for _ in range(count)])
            add(degree, knots, poles, weights)
    for degree in (4, 5, 6, 7, 9, 12):
        for trial in range(3):
            count = degree + 1 + random.randint(0, 4)
            knots = ([float(k) for k in range(count + degree + 1)] if trial == 2
                     else clamped(degree, count, trial == 1))
            add(degree, knots, [(random.uniform(-5, 5), random.uniform(-5, 5))
                                for _ in range(count)], [])
    cusp = [(0, 0), (1, 1), (0, 1), (1, 0)]
    cusps = 0
    while cusps < 6:
        m = [[random.randint(-3, 3) for _ in range(3)] for _ in range(2)]
        m.append([random.randint(-1, 1), random.randint(-1, 1), random.randint(3, 6)])
        h = [(m[0][0] * x + m[0][1] * y + m[0][2], m[1][0] * x + m[1][1] * y + m[1][2],
              m[2][0] * x + m[2][1] * y + m[2][2]) for x, y in cusp]
        if len({(a / c, b / c) for a, b, c in h}) == 4:
            add(3, [0, 0, 0, 0, 1, 1, 1, 1], [(a / c, b / c) for a, b, c in h], [c for _, _, c in h])
            cusps += 1
    for _ in range(4):
        xs = [random.uniform(-3, 3) for _ in range(4)]
        add(3, [0, 0, 0, 0, 1, 1, 1, 1], [(x, 0.5 * x) for x in xs],
            [10.0 ** random.uniform(-3, 3) for _ in range(4)])
        xs = [random.uniform(-3, 3) for _ in range(3)]
        add(2, [0, 0, 0, 1, 1, 1], [(x, -2 * x + 1) for x in xs],
            [10.0 ** random.uniform(-2, 2) for _ in range(3)])
    return out


def written(value):
    """A value to 30 significant digits."""
    number = Decimal(mpmath.nstr(value, 45, min_fixed=-mp.inf, max_fixed=mp.inf))
    if number == 0:
        return "0"
    with localcontext() as context:
        context.prec = 80
        if Decimal("1e-3") <= abs(number) < Decimal("1e6"):
            return format(number.quantize(Decimal(1).scaleb(number.adjusted() - 29)), "f")
        mantissa, exponent = format(number, ".29e").split("e")
        return mantissa + "e" + str(int(exponent))


def double(value):
    """A double as C reads it back exactly."""
    text = repr(float(value))
    return text[:-2] if text.endswith(".0") else text


def main():
    first = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    for degree, knots, poles, weights in drawn()[first:]:
        length, area, scale = measured(degree, knots, poles, weights)
        print("{.Degree  = %d,\n .Knots   = {%s},\n .Points  = {%s},\n .Weights = {%s},\n"
              " .Length  = %s,\n .Area    = %s,\n .Scale   = %s},"
              % (degree, ", ".join(double(k) for k in knots),
                 ", ".join(double(c) for p in poles for c in p),
                 ", ".join(double(w) for w in weights), written(length), written(area),
                 double(scale)), flush=True)


if __name__ == "__main__":
    main()
