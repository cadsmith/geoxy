#!/usr/bin/env python3
"""Works out, at 400 significant digits, where the lines, circular arcs and polylines of a case file
truly meet.

    python3 tools/reference/truth.py CASES TRUTH

Each case's curves are read as the doubles every library is given: a line by its ends, a ray by its
start and way, an endless line by a point and its way, an arc by its center, radius, start and sweep
(radians, or degrees with a d after them), a polyline by its vertices and bulges. A polyline is its
segments: a line where the bulge is under 1e-6 in magnitude or the ends are one point, as geoxy's
Straight has it, else the arc of the bulge's exact circle, whose center and radius are worked out
from the vertices and the bulge, the arc the circle's points on the bulge's side of the chord. Their
carriers' crossings are worked out in closed form, a line and a line by their determinant, a line
and a circle by its quadratic, two circles by their radical line, at 400 significant digits in
mpmath, so that a coordinate truly zero, of which a square root leaves a trace, comes out under the
least double and rounds to zero; each kept where it lies on both curves, within a segment's ends or
an arc's sweep, to 10^-90 of their size; the double root of a quadratic whose discriminant is zero
to the last digit a touch. A point two segments of a polyline both pass through, as a vertex, is
written once. Each point is written as the double nearest it, rounded once, and with 40 of its
digits.

Cases with any other curve, and pairs with two pieces whose carriers are one line or one circle, are
left out: a shared stretch is not a point. Where curves only come within a tolerance of meeting,
there is no point to work out, and nothing is written.

mpmath is needed: python3 -m venv venv && venv/bin/pip install mpmath.
"""

import sys
from fractions import Fraction

import mpmath
from mpmath import mp, mpf

mp.dps = 120
WIDE = 400  # The digits a case file is worked out at; 120 for a caller that imports what is here.
SLACK = mpf(10) ** -90
STRAIGHT = mpf(1e-6)  # The bulge under which a polyline's segment is a line, as geoxy's Straight.
STRAIGHTS = ("line", "ray", "xline")


def parsed(text):
    """A curve from its text: (kind, numbers), angles in radians; None for any other kind."""
    words = text.split()
    kind, values = words[0], words[1:]
    if kind in STRAIGHTS:
        return kind, [mpf(float(value)) for value in values[:4]]
    if kind == "polyline":
        count = int(values[1])
        numbers = [mpf(float(value)) for value in values[2:2 + 3 * count]]
        return kind, (values[0] == "closed", [numbers[3 * i:3 * i + 3] for i in range(count)])
    if kind == "arc":
        numbers = [mpf(float(value)) for value in values[:3]]
        for value in values[3:5]:
            if value.endswith("d"):
                numbers.append(mpf(float(value[:-1])) * mp.pi / 180)
            else:
                numbers.append(mpf(float(value)))
        return kind, numbers
    return None


def run(curve):
    """A straight curve's point, way and stretch of its parameter."""
    kind, values = curve
    if kind == "line":
        x0, y0, x1, y1 = values
        return (x0, y0), (x1 - x0, y1 - y0), (mpf(0), mpf(1))
    low = mpf(0) if kind == "ray" else -mp.inf
    return (values[0], values[1]), (values[2], values[3]), (low, mp.inf)


def pieces(curve):
    """A curve's pieces: a polyline's segments, each a line or a bulge, any other curve itself. A
    bulge is its circle's center and radius, then its vertices and its bulge, worked out from them:
    the center off the chord's middle by (1 - b^2) / (4 b) of the chord turned a quarter left."""
    kind, values = curve
    if kind != "polyline":
        return [curve]
    closed, vertices = values
    count = len(vertices) if closed and len(vertices) > 1 else len(vertices) - 1
    found = []
    for index in range(max(count, 0)):
        (x0, y0, bulge), (x1, y1, _) = vertices[index], vertices[(index + 1) % len(vertices)]
        if not abs(bulge) >= STRAIGHT or (x0, y0) == (x1, y1):
            found.append(("line", [x0, y0, x1, y1]))
            continue
        ends = [rational(value) for value in (x0, y0, x1, y1, bulge)]
        offset = (1 - ends[4] * ends[4]) / (4 * ends[4])
        cx = (ends[0] + ends[2]) / 2 - (ends[3] - ends[1]) * offset
        cy = (ends[1] + ends[3]) / 2 + (ends[2] - ends[0]) * offset
        square = (ends[0] - cx) ** 2 + (ends[1] - cy) ** 2
        found.append(("bulge", [real(cx), real(cy), mpmath.sqrt(real(square)), x0, y0, x1, y1,
                                bulge, (cx, cy, square)]))
    return found


def rational(value):
    """A double, held as an mpf, as the exact fraction it is."""
    return Fraction(float(value))


def real(value):
    """A fraction at the working precision."""
    return mpf(value.numerator) / value.denominator


def circled(curve):
    """An arc's or a bulge's circle exactly: its center and squared radius, as fractions."""
    if curve[0] == "bulge":
        return curve[1][8]
    cx, cy, radius = (rational(value) for value in curve[1][:3])
    return cx, cy, radius * radius


def on(curve, point, along):
    """Whether a point of a curve's carrier lies on the curve: within a run's stretch, an arc's
    sweep or on a bulge's side of its chord, closed, to SLACK of the curve's size."""
    kind, values = curve
    if kind == "bulge":
        _, _, radius, x0, y0, x1, y1, bulge, _ = values
        side = (x1 - x0) * (point[1] - y0) - (y1 - y0) * (point[0] - x0)
        size = 1 + abs(x0) + abs(y0) + abs(x1) + abs(y1) + radius
        return side * mpmath.sign(bulge) <= SLACK * size * size
    if kind != "arc":
        low, high = run(curve)[2]
        (px, py), (wx, wy), _ = run(curve)
        size = 1 + abs(px) + abs(py) + abs(wx) + abs(wy)
        slack = SLACK * size / mpmath.sqrt(wx * wx + wy * wy)
        return low - slack <= along <= high + slack
    cx, cy, _, start, sweep = values
    if abs(sweep) >= 2 * mp.pi:
        return True
    angle = mpmath.atan2(point[1] - cy, point[0] - cx)
    turned = (angle - start) % (2 * mp.pi) if sweep > 0 else (start - angle) % (2 * mp.pi)
    slack = SLACK * (1 + abs(start) + abs(sweep))
    return turned <= abs(sweep) + slack or turned >= 2 * mp.pi - slack


def crossings(first, second):
    """The points the two carriers share: (point, place on a straight first, touch); None where
    they are one carrier."""
    straight = [curve[0] in STRAIGHTS for curve in (first, second)]
    if all(straight):
        (px, py), (ux, uy), _ = run(first)
        (qx, qy), (vx, vy), _ = run(second)
        turn = ux * vy - uy * vx
        if turn == 0:
            return None
        along = ((qx - px) * vy - (qy - py) * vx) / turn
        beyond = ((qx - px) * uy - (qy - py) * ux) / turn
        return [((px + ux * along, py + uy * along), (along, beyond), False)]
    if any(straight):
        line, arc = (first, second) if straight[0] else (second, first)
        (px, py), (ux, uy), _ = run(line)
        px, py, ux, uy = (rational(value) for value in (px, py, ux, uy))
        cx, cy, square = circled(arc)
        gx, gy = px - cx, py - cy
        a = ux * ux + uy * uy
        b = ux * gx + uy * gy
        c = gx * gx + gy * gy - square
        d = b * b - a * c
        if d < 0 or a == 0:
            return []
        roots = [mpf(0)] if d == 0 else [mpmath.sqrt(real(d)), -mpmath.sqrt(real(d))]
        found = []
        for root in roots:
            along = (-real(b) + root) / real(a)
            point = (real(px) + real(ux) * along, real(py) + real(uy) * along)
            found.append((point, (along, along), d == 0))
        return found
    (ax, ay, inner), (bx, by, outer) = circled(first), circled(second)
    gx, gy = bx - ax, by - ay
    apart = gx * gx + gy * gy
    if apart == 0:
        return None
    chord = apart + inner - outer
    big = 4 * apart * inner - chord * chord
    if big < 0:
        return []
    roots = [mpf(0)] if big == 0 else [mpmath.sqrt(real(big)), -mpmath.sqrt(real(big))]
    found = []
    for root in roots:
        x = real(ax) + (real(chord * gx) - real(gy) * root) / (2 * real(apart))
        y = real(ay) + (real(chord * gy) + real(gx) * root) / (2 * real(apart))
        found.append(((x, y), (mpf(0), mpf(0)), big == 0))
    return found


def piecewise(first, second):
    """Where two pieces truly meet: each point and whether a touch; None where they share a
    carrier."""
    found = crossings(first, second)
    if found is None:
        return None
    points = []
    for point, (along, beyond), touch in found:
        other = beyond if first[0] in STRAIGHTS else along
        if on(first, point, along) and on(second, point, other):
            points.append((point, touch))
    return points


def met(first, second):
    """Where two curves truly meet: each point and whether a touch, a point two pieces both pass
    through once; None where two pieces share a carrier."""
    points = []
    for one in pieces(first):
        for other in pieces(second):
            found = piecewise(one, other)
            if found is None:
                return None
            for point, touch in found:
                size = 1 + abs(point[0]) + abs(point[1])
                if all(mpmath.hypot(point[0] - x, point[1] - y) > SLACK * size
                       for (x, y), _ in points):
                    points.append((point, touch))
    return points


def nearest(value):
    """The double nearest a number: Python's conversion of its digits rounds once, correctly."""
    return float(mpmath.nstr(value, 100, strip_zeros=False, min_fixed=-mp.inf, max_fixed=mp.inf))


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    mp.dps = WIDE
    out = []
    name, curves = None, {}
    for line in open(sys.argv[1], encoding="utf-8"):
        word, _, rest = line.strip().partition(" ")
        if word == "case":
            name, curves = rest, {}
        elif word in ("first", "second"):
            curves[word] = parsed(rest)
        elif word == "end" and name is not None:
            first, second = curves.get("first"), curves.get("second")
            if first and second:
                points = met(first, second)
                if points is not None:
                    out.append(f"case {name}")
                    for (x, y), touch in points:
                        tail = " touch" if touch else ""
                        out.append(
                            f"truth {nearest(x)!r} {nearest(y)!r}{tail} "
                            f"{mpmath.nstr(x, 40)} {mpmath.nstr(y, 40)}"
                        )
                    out.append("end")
            name = None
    with open(sys.argv[2], "w", encoding="utf-8") as truth:
        truth.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
