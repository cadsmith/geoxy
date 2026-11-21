#!/usr/bin/env python3
"""Exact answers to what a point and a box ask of splines whose weights lie far apart.

    python3 tools/reference/stiff.py [seed [rounds]] > stiff.txt

For each of six shapes of spline and each weight from 2^-1020 to 2^1020, random points and boxes
are asked of the spline: its nearest place (N), the feet of the perpendiculars (F), the places its
tangent passes through the point (T), and the stretches within the box (C). A line is a question
and its answer: the kind, the shape, the power of two, the spline (degree, knots, poles with their
weights), the question, how many steps of a double the answer moves when the question moves one
(1e300 where a place comes or goes), and the answer, every number as C writes a double in sixteens.
A place is its point, its parameter and the point's speed by the parameter; a stretch its two
parameters and the speeds there.

Everything is rational arithmetic on the doubles given, with nothing rounded before the answer is
written: roots are isolated by Sturm sequences and refined by bisection to a relative 2^-200, from
the start of a piece or from its end, whichever is nearer, so a root at 1 - 2^-1000 is told as
exactly as one at 2^-1000. A question that is not generic, a root at an end or a repeated one, is
left out. tests/geoxy/stiff.hpp holds a few of these lines a shape and weight, those the question's
last digit moves less than two steps, written as the shape's name, the power, the question and
the answer.

It needs mpmath, for the square roots of the distances and speeds it writes.
"""

import random
import sys
from fractions import Fraction as F
from math import gcd

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 23)


def trim(p):
    while len(p) > 1 and p[-1] == 0:
        p = p[:-1]
    return p


def add(a, b):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])


def scale(a, k):
    return trim([k * v for v in a])


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def der(a):
    return trim([i * a[i] for i in range(1, len(a))] or [F(0)])


def val(p, x):
    out = F(0)
    for c in reversed(p):
        out = out * x + c
    return out


def whole(p):
    """The polynomial times a positive number that makes its coefficients whole and coprime."""
    lcm = 1
    for c in p:
        lcm = lcm * c.denominator // gcd(lcm, c.denominator)
    ints = [int(c * lcm) for c in p]
    g = 0
    for v in ints:
        g = gcd(g, abs(v))
    return [F(v // g) for v in ints] if g else [F(0)]


def rem(a, b):
    a = list(a)
    while len(a) >= len(b) and trim(a) != [0]:
        k = a[-1] / b[-1]
        shift = len(a) - len(b)
        for i, c in enumerate(b):
            a[shift + i] -= k * c
        a = a[:-1]
        if not a:
            return [F(0)]
    return trim(a)


def sturm(p):
    """The Sturm sequence of a polynomial, or None where it has a repeated root."""
    seq = [whole(p), whole(der(p))]
    while len(seq[-1]) > 1:
        r = rem(seq[-2], seq[-1])
        if r == [0]:
            return None
        seq.append(whole(scale(r, -1)))
    return seq if seq[-1] != [0] else None


def changes(seq, x):
    signs = [v for v in (val(p, x) for p in seq) if v != 0]
    return sum(1 for a, b in zip(signs, signs[1:]) if (a > 0) != (b > 0))


class Odd(Exception):
    """A case that is not generic: a root at an end or at a place tried, or a repeated one."""


def isolate(seq, lo, hi):
    """Stretches (lo, hi] of one root each, lo > 0, within (lo, hi]."""
    if val(seq[0], lo) == 0 or val(seq[0], hi) == 0:
        raise Odd
    n = changes(seq, lo) - changes(seq, hi)
    if n == 0:
        return []
    if n == 1 and lo > 0:
        return [(lo, hi)]
    if lo == 0:
        mid = hi / 2 ** 64
        if mid < F(1, 2 ** 8000):
            raise Odd
    elif hi / lo > 4:
        steps = (hi / lo).numerator.bit_length() - (hi / lo).denominator.bit_length()
        mid = lo * 2 ** max(1, steps // 2)
    else:
        mid = (lo + hi) / 2
    return isolate(seq, lo, mid) + isolate(seq, mid, hi)


def refine(p, lo, hi):
    low = val(p, lo) > 0
    while (hi - lo) * 2 ** 200 > lo:
        mid = (lo + hi) / 2
        v = val(p, mid)
        if v == 0:
            raise Odd
        if (v > 0) == low:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def flip(p):
    """p(1 - s) as a polynomial in s."""
    out = [F(0)]
    power = [F(1)]
    for c in p:
        out = add(out, scale(power, c))
        power = mul(power, [F(1), F(-1)])
    return out


def roots(p):
    """The roots of a polynomial in the open stretch 0 to 1, in order, each as a rational within a
    relative 2^-200 of it, counted from the nearer end."""
    p = trim(p)
    if p == [0]:
        raise Odd
    if val(p, F(0)) == 0 or val(p, F(1)) == 0 or val(p, F(1, 2)) == 0:
        raise Odd
    seq = sturm(p)
    if seq is None:
        raise Odd
    found = [refine(seq[0], lo, hi) for lo, hi in isolate(seq, F(0), F(1, 2))]
    q = flip(p)
    back = sturm(q)
    late = [1 - refine(back[0], lo, hi) for lo, hi in isolate(back, F(0), F(1, 2))]
    return found + sorted(late)


def blossom(spline, span, args):
    degree, knots, poles = spline
    d = [poles[span - degree + j] for j in range(degree + 1)]
    for r in range(1, degree + 1):
        u = args[r - 1]
        for j in range(degree, r - 1, -1):
            a = knots[span - degree + j]
            b = knots[span + j - r + 1]
            alpha = (u - a) / (b - a)
            d[j] = tuple((1 - alpha) * x + alpha * y for x, y in zip(d[j - 1], d[j]))
    return d[degree]


def pieces(spline):
    """Each span that is not empty as (start, end, X, Y, W), the three as polynomials of the span's
    own parameter from 0 to 1."""
    degree, knots, poles = spline
    out = []
    for span in range(degree, len(poles)):
        a, b = knots[span], knots[span + 1]
        if a == b:
            continue
        bezier = [blossom(spline, span, [a] * (degree - j) + [b] * j) for j in range(degree + 1)]
        coords = []
        for axis in range(3):
            poly = [F(0)]
            for i in range(degree + 1):
                # C(n, i) t^i (1 - t)^(n - i)
                term = [F(1)]
                for _ in range(i):
                    term = mul(term, [F(0), F(1)])
                for _ in range(degree - i):
                    term = mul(term, [F(1), F(-1)])
                comb = 1
                for j in range(i):
                    comb = comb * (degree - j) // (j + 1)
                poly = add(poly, scale(term, comb * bezier[i][axis]))
            coords.append(poly)
        out.append((a, b, coords[0], coords[1], coords[2]))
    return out


def point(piece, u):
    _, _, X, Y, W = piece
    w = val(W, u)
    return val(X, u) / w, val(Y, u) / w


def speed(piece, u):
    """How fast the spline's point moves by its parameter at a place of a piece, as a float."""
    a, b, X, Y, W = piece
    w = val(W, u)
    dx = (w * val(der(X), u) - val(der(W), u) * val(X, u)) / (w * w)
    dy = (w * val(der(Y), u) - val(der(W), u) * val(Y, u)) / (w * w)
    square = (dx * dx + dy * dy) / ((b - a) * (b - a))
    from mpmath import mp, mpf, sqrt
    mp.dps = 40
    return float(sqrt(mpf(square.numerator) / mpf(square.denominator)))


def along(piece, u):
    a, b = piece[0], piece[1]
    return a + (b - a) * u


def hexed(v):
    return float.hex(float(v))


def seen(piece, p):
    _, _, X, Y, W = piece
    return sub(X, scale(W, p[0])), sub(Y, scale(W, p[1])), W


def feet(spline, p):
    out = []
    for piece in pieces(spline):
        U, V, W = seen(piece, p)
        n = sub(mul(W, add(mul(U, der(U)), mul(V, der(V)))), mul(der(W), add(mul(U, U), mul(V, V))))
        out += [(piece, u) for u in roots(n)]
    return out


def touches(spline, p):
    out = []
    for piece in pieces(spline):
        U, V, _ = seen(piece, p)
        out += [(piece, u) for u in roots(sub(mul(U, der(V)), mul(V, der(U))))]
    return out


def told(places):
    return " ".join(f"{hexed(point(piece, u)[0])} {hexed(point(piece, u)[1])} {hexed(along(piece, u))} "
                    f"{float.hex(speed(piece, u))}" for piece, u in places)


def squared(piece, u, p):
    x, y = point(piece, u)
    return (x - p[0]) ** 2 + (y - p[1]) ** 2


STEP = F(1, 2 ** 52)


def swayed(found, spline, p):
    """How many steps a place moves, at most, when the point asked from moves one step: how finely
    the question can be answered at all. 1e300 where a place comes or goes."""
    base = [point(*place) for place in found(spline, p)]
    worst = 0.0
    for dx, dy in ((1, 0), (0, 1)):
        try:
            moved = [point(*place) for place in found(spline, (p[0] + dx * STEP, p[1] + dy * STEP))]
        except Odd:
            return 1e300
        if len(moved) != len(base):
            return 1e300
        for a, b in zip(base, moved):
            away = max(abs(a[0] - b[0]), abs(a[1] - b[1])) / STEP
            worst = max(worst, float(away) if away < 10 ** 290 else 1e300)
    return worst


def closest(spline, p):
    return nearest(spline, p, True)


def nearest(spline, p, place=False):
    places = feet(spline, p)
    for piece in pieces(spline):
        places += [(piece, F(0)), (piece, F(1))]
    places.sort(key=lambda place: squared(place[0], place[1], p))
    best, next_ = squared(*places[0], p), None
    for place in places[1:]:
        if along(*place) != along(*places[0]):
            next_ = squared(*place, p)
            break
    sure = 0 if next_ is not None and (next_ - best) < next_ / 10 ** 6 else 1
    if place:
        return places[:1]
    from mpmath import mp, mpf, sqrt
    mp.dps = 60
    distance = sqrt(mpf(best.numerator) / mpf(best.denominator))
    return f"{float.hex(float(distance))} {sure} {told(places[:1])}"


def clipped(spline, box):
    return clip(spline, box, True)


def clip(spline, box, ends=False):
    (x0, y0), (x1, y1) = box
    cuts = []
    for piece in pieces(spline):
        _, _, X, Y, W = piece
        for poly, level in ((X, x0), (X, x1), (Y, y0), (Y, y1)):
            cuts += [(piece, u) for u in roots(sub(poly, scale(W, level)))]
    ps = pieces(spline)
    marks = [(ps[0], F(0))] + sorted(cuts, key=lambda c: along(*c)) + [(ps[-1], F(1))]
    for a, b in zip(marks, marks[1:]):
        if along(*a) == along(*b):
            raise Odd
    ranges = []
    for a, b in zip(marks, marks[1:]):
        middle = (along(*a) + along(*b)) / 2
        piece = next(q for q in ps if q[0] <= middle <= q[1])
        x, y = point(piece, (middle - piece[0]) / (piece[1] - piece[0]))
        if x in (x0, x1) or y in (y0, y1):
            raise Odd
        if x0 < x < x1 and y0 < y < y1:
            if ranges and along(*ranges[-1][1]) == along(*a):
                ranges[-1] = (ranges[-1][0], b)
            else:
                ranges.append((a, b))
    if ends:
        return [end for stretch in ranges for end in stretch]
    out = [str(len(ranges))]
    for a, b in ranges:
        out.append(f"{hexed(along(*a))} {hexed(along(*b))} {float.hex(speed(*a))} {float.hex(speed(*b))}")
    return " ".join(out)


def family(name, k):
    m = 2.0 ** k
    if name == "quadratic":
        return 2, [0, 0, 0, 1, 1, 1], [(0, 0, 1), (1, 1, m), (2, 0, 1)]
    if name == "cubic":
        return 3, [0, 0, 0, 0, 1, 1, 1, 1], [(0, 0, 1), (0.5, 1.5, m), (1.5, 1.5, m), (2, 0, 1)]
    if name == "spans":
        return 3, [0, 0, 0, 0, 1, 2, 2, 2, 2], [(0, 0, 1), (0.5, 1.5, m), (1.5, 1.5, 1), (2.5, 0.25, 1 / m), (3, 1, 1)]
    if name == "ends":
        return 3, [0, 0, 0, 0, 1, 1, 1, 1], [(0, 0, m), (0.5, 1.5, 1), (1.5, 1.5, 1), (2, 0, 1 / m)]
    if name == "odd":
        return 2, [0, 0, 0, 1, 1, 1], [(0.1, 0.2, 1), (1.3, 1.1, 1.37 * m), (2.2, -0.3, 0.83)]
    return 2, [0, 0, 0, 1, 2, 3, 3, 3], [(0, 0, 1), (1, 1, m), (2, 0, 1), (3, 1, m), (4, 0, 1)]


def exact(shape):
    degree, knots, poles = shape
    return degree, [F(k) for k in knots], [(F(x) * F(w), F(y) * F(w), F(w)) for x, y, w in poles]


def header(op, name, k, shape):
    degree, knots, poles = shape
    return (f"{op} {name} {k} {degree} {len(knots)} {' '.join(float.hex(float(v)) for v in knots)} "
            f"{len(poles)} {' '.join(float.hex(float(v)) for pole in poles for v in pole)}")


SPREADS = [0, 1, 2, 3, 4, 8, 16, 30, 45, 52, 60, 100, 200, 400, 700, 1000, 1020]
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 6
for name in ("quadratic", "cubic", "spans", "waves", "ends", "odd"):
    for size in SPREADS:
        for k in ((size, -size) if size else (0,)):
            shape = family(name, k)
            spline = exact(shape)
            reach = 4.0 if name == "waves" else 3.0 if name == "spans" else 2.0
            for _ in range(ROUNDS):
                p = (random.uniform(-1, reach + 1), random.uniform(-1, 2.5))
                pf = (F(p[0]), F(p[1]))
                at = f"{float.hex(p[0])} {float.hex(p[1])}"
                try:
                    print(header("N", name, k, shape), at, float.hex(swayed(closest, spline, pf)), nearest(spline, pf))
                except Odd:
                    pass
                try:
                    places = feet(spline, pf)
                    print(header("F", name, k, shape), at, float.hex(swayed(feet, spline, pf)), len(places), told(places))
                except Odd:
                    pass
                try:
                    places = touches(spline, pf)
                    print(header("T", name, k, shape), at, float.hex(swayed(touches, spline, pf)), len(places), told(places))
                except Odd:
                    pass
                x0, y0 = random.uniform(-0.5, reach - 0.8), random.uniform(-0.5, 1.0)
                x1, y1 = x0 + random.uniform(0.3, reach), y0 + random.uniform(0.2, 1.5)
                try:
                    box = ((F(x0), F(y0)), (F(x1), F(y1)))
                    base = [point(*end) for end in clipped(spline, box)]
                    sway = 0.0
                    try:
                        moved = [point(*end) for end in clipped(spline, tuple((c[0] + STEP, c[1] + STEP) for c in box))]
                        if len(moved) != len(base):
                            sway = 1e300
                        for a, b in zip(base, moved):
                            sway = max(sway, float(max(abs(a[0] - b[0]), abs(a[1] - b[1])) / STEP))
                    except Odd:
                        sway = 1e300
                    print(header("C", name, k, shape), " ".join(float.hex(v) for v in (x0, y0, x1, y1)),
                          float.hex(sway), clip(spline, box))
                except Odd:
                    pass
            sys.stdout.flush()
