#!/usr/bin/env python3
"""Forges geoxy's own suite of where two curves meet, every answer known without geoxy.

    python3 -I tools/suites/forge.py > tests/suites/geoxy/cases/intersect.txt
    python3 -I tools/suites/forge.py --check tests/suites/geoxy/cases/intersect.txt

With --check, nothing is written: the exit status is 1 if the file is not what the forge makes,
so that the suite in the tree never drifts from the forge and its oracle.

Each case is built from exact rationals: small integers, Pythagorean circles, bulges whose arcs
have rational centers, Bézier pieces cut by exact de Casteljau steps. The oracle beside this
file works out where the curves meet exactly. Cases built to sit a fraction of the tolerance
from a tangency or an end take their answer from the rule the tolerance states: within it, one
touch or one end; beyond it, nothing or two crossings. Each base case is written nine times: as
built, turned a quarter, mirrored, each at the origin, at 5e6 and at 1e8, with tolerances that
are powers of two so that every shifted coordinate is a double exactly.

The hard family pushes where kernels break: crossings at angles near none, secants a few
tolerances inside a tangency, contact of order three and four, crossings clustered closer than
the tolerance, cusps, a loop's double point, seams, a radius of 2^20, a sweep of 2^-11, an
ellipse 2^19 times as long as wide, stretches shared over a sliver.

Floats touch a case in two places only, and both are guarded. Every input must be a double
exactly, or the forge stops: geoxy would read another curve than the oracle solved. An arc is
written by its angles, so geoxy's ends may sit a few steps of a double from the oracle's; a case
with an answer near such an end is refused unless that drift is under 1/64 of the tolerance.
The answers are rounded once, to the nearest double, when written. Each hit's slack is the
tolerance's own: for a crossing, how far along it both curves stay within the tolerance, the
tolerance over the sine of half the angle they meet at; for a touch, its band. Hits are merged
only where they are one place on both curves, and a hit is a crossing only where the curve
passes from one side of the other to its other side, as geoxy's Intersect states it.
"""

import contextlib
import io
import itertools
import math
import pathlib
import sys
from dataclasses import dataclass, field
from fractions import Fraction as F

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import oracle  # noqa: E402  pylint: disable=wrong-import-position

SIZE = 20
SCALES = ((F(0), F(0), F(1, 2 ** 30)), (F(5_000_000), F(-3_000_000), F(1, 2 ** 24)),
          (F(100_000_000), F(100_000_000), F(1, 2 ** 20)))


@dataclass
class Frame:
    """Where a base case is put: a quarter turn or a mirror, then a shift; and its tolerance."""
    dx: F
    dy: F
    turn: int
    mirror: bool
    tolerance: F

    def vector(self, v):
        """A vector turned and mirrored as the frame says, without the shift."""
        x, y = F(v[0]), F(v[1])
        for _ in range(self.turn):
            x, y = -y, x
        return (x, -y) if self.mirror else (x, y)

    def point(self, p):
        """A point put in the frame: turned or mirrored, then shifted."""
        x, y = self.vector(p)
        return (x + self.dx, y + self.dy)

    @property
    def flips(self):
        """Whether the frame mirrors, which reverses the sense of an arc."""
        return self.mirror


@dataclass
class Shape:
    """A curve as geoxy's text states it, its oracle pieces and where it ends."""
    text: str
    pieces: list
    ends: list = field(default_factory=list)
    drift: float = 0.0  # how far geoxy's ends may be from the oracle's, read from float angles


def r(value):
    """An input number as the shortest text that reads back as the same double: it must be a
    double exactly, or geoxy would read another curve than the one the oracle solves."""
    if F(float(value)) != F(value):
        raise ValueError(f"{value} is not a double exactly")
    return repr(float(value))


def pt(p):
    """An input point as two numbers in text."""
    return f"{r(p[0])} {r(p[1])}"


def answer(p):
    """An answer's point as two numbers in text, each the double nearest it."""
    return f"{float(p[0])!r} {float(p[1])!r}"


def segment(frame, a, b):
    """The line segment from a to b, put in the frame."""
    a, b = frame.point(a), frame.point(b)
    way = (b[0] - a[0], b[1] - a[1])
    return Shape(f"line {pt(a)} {pt(b)}", [oracle.line(a, way, F(0), F(1))], [a, b])


def ray(frame, a, w):
    """The ray from a along w, put in the frame."""
    a, w = frame.point(a), frame.vector(w)
    return Shape(f"ray {pt(a)} {pt(w)}", [oracle.line(a, w, F(0), None)], [a])


def xline(frame, a, w):
    """The endless line through a along w, put in the frame."""
    a, w = frame.point(a), frame.vector(w)
    return Shape(f"xline {pt(a)} {pt(w)}", [oracle.line(a, w, None, None)], [])


def angles(start, end, counterclockwise, full):
    """The start angle and the signed sweep of an arc, a whole turn if it is full."""
    first = math.atan2(float(start[1]), float(start[0]))
    if full:
        return first, (2 * math.pi if counterclockwise else -2 * math.pi)
    last = math.atan2(float(end[1]), float(end[0]))
    whole = 2 * math.pi
    sweep = (last - first) % whole if counterclockwise else -((first - last) % whole)
    return first, sweep


def arc_piece(center, radius, start, end, counterclockwise, full):
    """The oracle's piece for a circle's arc."""
    return oracle.conic(center, (radius, F(0)), (F(0), radius), start, end, counterclockwise, full)


def arc(frame, center, radius, start, end, counterclockwise=True, full=False):
    """A circle's arc from the unit vector start to end, both rational."""
    center = frame.point(center)
    start, end = frame.vector(start), frame.vector(end)
    counterclockwise = counterclockwise != frame.flips
    first, sweep = angles(start, end, counterclockwise, full)
    text = f"arc {pt(center)} {r(radius)} {first!r} {sweep!r}"
    piece = arc_piece(center, F(radius), start, end, counterclockwise, full)
    ends = [] if full else [(center[0] + radius * start[0], center[1] + radius * start[1]),
                            (center[0] + radius * end[0], center[1] + radius * end[1])]
    return Shape(text, [piece], ends, drifted(radius, first, sweep))


def drifted(reach, first, sweep):
    """How far an arc's ends, read from float angles, may be from its exact ends: four steps of
    a double in each angle, its start and its end, at the arc's reach from its center."""
    return float(reach) * 4 * (math.ulp(abs(first)) + math.ulp(abs(first) + abs(sweep)))


def ellipse(frame, origin, xaxis, yaxis, start, end, counterclockwise=True, full=False):
    """An ellipse's arc: origin + xaxis cos + yaxis sin, from the unit vector start to end."""
    origin = frame.point(origin)
    xaxis, yaxis = frame.vector(xaxis), frame.vector(yaxis)
    first, sweep = angles(start, end, counterclockwise, full)
    text = (f"ellipse {pt(xaxis)} {pt(yaxis)} {pt(origin)} {first!r} {sweep!r}")
    piece = oracle.conic(origin, xaxis, yaxis, start, end, counterclockwise, full)
    at = lambda q: (origin[0] + xaxis[0] * q[0] + yaxis[0] * q[1],  # noqa: E731
                    origin[1] + xaxis[1] * q[0] + yaxis[1] * q[1])
    reach = max(math.hypot(float(xaxis[0]), float(xaxis[1])),
                math.hypot(float(yaxis[0]), float(yaxis[1])))
    return Shape(text, [piece], [] if full else [at(start), at(end)],
                 drifted(reach, first, sweep))


def exact_root(value):
    """The square root of a rational that is a perfect square, else a ValueError."""
    value = F(value)
    num, den = math.isqrt(value.numerator), math.isqrt(value.denominator)
    if num * num != value.numerator or den * den != value.denominator:
        raise ValueError(f"no rational root of {value}")
    return F(num, den)


def polyline(frame, vertices, closed=False):
    """A polyline of (point, bulge) vertices, each bulged chord of rational length."""
    placed = [(frame.point(p), F(float(b)) * (-1 if frame.flips else 1)) for p, b in vertices]
    words = [f"{pt(p)} {r(b)}" for p, b in placed]
    if len(placed) == 1:
        lone = placed[0][0]
        point = oracle.Piece(x=[lone[0]], y=[lone[1]], w=[F(1)], lo=F(0), hi=F(0))
        return Shape(f"polyline open 1 {words[0]}", [point], [lone])
    text = f"polyline {'closed' if closed else 'open'} {len(placed)} {' '.join(words)}"
    pieces = []
    count = len(placed) if closed else len(placed) - 1
    for index in range(count):
        (a, bulge), (b, _) = placed[index], placed[(index + 1) % len(placed)]
        d = (b[0] - a[0], b[1] - a[1])
        if bulge == 0:
            pieces.append(oracle.line(a, d, F(0), F(1)))
            continue
        middle = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
        factor = (1 - bulge * bulge) / (4 * bulge)
        center = (middle[0] - d[1] * factor, middle[1] + d[0] * factor)
        radius = exact_root(d[0] ** 2 + d[1] ** 2) * (1 + bulge * bulge) / (4 * abs(bulge))
        start = ((a[0] - center[0]) / radius, (a[1] - center[1]) / radius)
        end = ((b[0] - center[0]) / radius, (b[1] - center[1]) / radius)
        pieces.append(arc_piece(center, radius, start, end, bulge > 0, False))
    ends = [] if closed else [placed[0][0], placed[-1][0]]
    return Shape(text, pieces, ends)


def bezier(frame, points):
    """A quadratic or cubic Bézier curve through these control points, put in the frame."""
    points = [frame.point(p) for p in points]
    word = {3: "quadratic", 4: "cubic"}[len(points)]
    return Shape(f"{word} {' '.join(pt(p) for p in points)}", [oracle.bezier(points)],
                 [points[0], points[-1]])


def conic(frame, points, weight):
    """A rational quadratic Bézier curve, a conic's arc, through three control points, the
    middle one weighted, put in the frame: in geoxy's words a spline of degree two, rational."""
    points = [frame.point(p) for p in points]
    text = (f"spline 2 6 0 0 0 1 1 1 3 {' '.join(pt(p) for p in points)} "
            f"3 1.0 {r(weight)} 1.0")
    return Shape(text, [oracle.bezier(points, [F(1), F(weight), F(1)])], [points[0], points[-1]])


def spline(frame, pieces):
    """A cubic spline of Bézier pieces that share their ends, every inner knot threefold."""
    placed = [[frame.point(p) for p in piece] for piece in pieces]
    poles = [placed[0][0]]
    for piece in placed:
        poles.extend(piece[1:])
    knots = [0, 0, 0, 0]
    for index in range(1, len(placed)):
        knots.extend([index] * 3)
    knots.extend([len(placed)] * 4)
    text = (f"spline 3 {len(knots)} {' '.join(str(k) for k in knots)} {len(poles)} "
            f"{' '.join(pt(p) for p in poles)} 0")
    return Shape(text, [oracle.bezier(piece) for piece in placed], [placed[0][0], placed[-1][-1]])


def cut(points, along):
    """A Bézier curve cut at a parameter by de Casteljau: its two halves, exactly."""
    left, right, work = [points[0]], [points[-1]], list(points)
    while len(work) > 1:
        work = [((1 - along) * a[0] + along * b[0], (1 - along) * a[1] + along * b[1])
                for a, b in zip(work, work[1:])]
        left.append(work[0])
        right.append(work[-1])
    return left, list(reversed(right))


# What the oracle finds, merged as geoxy's Intersect states it.

def expected(first, second, tolerance):
    """The hits and shared stretches of two shapes, merged as geoxy's Intersect states it: each hit
    with its kind and the sine of the angle the curves meet at, None where it is not known.

    A hit within the tolerance of a curve's end is an end. A case where that, or two distinct
    hits merging, turns on less than a factor of four in the tolerance is refused: its answer
    would be a matter of rounding, not of the rule. So is one where a curve's end comes within
    four tolerances of the other curve without meeting it: its answer is the tolerance's rule
    for ends, which the oracle, knowing only where curves meet, does not judge."""
    hits, shared = [], []
    for i, a in enumerate(first.pieces):
        for j, b in enumerate(second.pieces):
            found = oracle.meet_pieces(a, b)
            for (p, tangent), sine, places in zip(found.hits, found.sines, found.places):
                keys = (None, None) if places is None else (
                    spot(first, i, p, places[0][0]), spot(second, j, p, places[1][0]))
                side = None if places is None else turned(places[0][1], places[1][1])
                hits.append((p, tangent, sine, keys, side))
            shared.extend(found.shared)
    merged = []
    for p, tangent, sine, keys, side in hits:
        same = next((m for m in merged if oracle.near(m[0], p) and alike(m[3], keys)), None)
        if same is None:
            merged.append([p, tangent, sine, keys, {side}])
        else:
            same[1] = same[1] or tangent
            same[2] = sine if same[2] is None else (same[2] if sine is None else min(same[2], sine))
            same[4].add(side)
    merged = [[p, tangent or 0 in sides or {1, -1} <= sides, sine, keys, sides]
              for p, tangent, sine, keys, sides in merged]
    joined = join(shared)
    ends = first.ends + second.ends
    out = []
    for p, tangent, sine, _, _ in merged:
        if (any(e is not None and oracle.near(p, e) for pair in joined for e in pair)
                or inside_shared(p, joined, first)):
            continue
        reach = min((squared(p, e) for e in ends), default=None)
        if reach is not None and tolerance ** 2 / 16 < reach < 16 * tolerance ** 2:
            raise ValueError(f"a hit sits near the tolerance from an end: {answer(p)}")
        ending = reach is not None and reach <= tolerance ** 2 / 16
        kind = "End" if ending else ("Touch" if tangent else "Cross")
        out.append((p, kind, 0.0 if tangent else sine))
    for shape, other in ((first, second), (second, first)):
        for e in shape.ends:
            if any(squared(p, e) <= tolerance ** 2 / 16 for p, _, _ in out):
                continue
            near = min(oracle.reach2(piece, e) for piece in other.pieces)
            if oracle.ON ** 2 < near < 16 * tolerance ** 2:
                raise ValueError(f"an end comes within a few tolerances of the other curve "
                                 f"without meeting it: {answer(e)}")
    for one, two in itertools.combinations(out, 2):
        if not oracle.near(one[0], two[0]) and squared(one[0], two[0]) < 16 * tolerance ** 2:
            raise ValueError(f"two hits near the tolerance apart: {answer(one[0])}")
    return out, joined


def turned(way, other):
    """The side a curve running one way passes to across another running the other: 1, -1, 0
    where they run along each other, None where a way is not known."""
    if way is None or other is None:
        return None
    turn = oracle.cross(way, other)
    return (turn > 0) - (turn < 0)


def ends_of(piece):
    """Where a piece starts and where it ends, None where it runs on or is whole."""
    if piece.arc is not None:
        both = oracle.ends(piece)
        return (both[0], both[1]) if both else (None, None)
    return tuple(None if t is None else oracle.point_at(piece, t) for t in (piece.lo, piece.hi))


def spot(shape, index, p, local):
    """A place on a shape, the same for the two pieces that meet at a joint: the joint's
    number there, else the piece's number and the place on it."""
    start, end = ends_of(shape.pieces[index])
    count = len(shape.pieces)
    closed = not shape.ends and count > 1
    if start is not None and oracle.near(start, p):
        return ("joint", index)
    if end is not None and oracle.near(end, p):
        return ("joint", (index + 1) % count if closed else index + 1)
    return (index, local)


def alike(keys, others):
    """Whether two hits' places are the same on both shapes; an unknown place matches any."""
    def same(one, two):
        """Whether two places on one shape are the same."""
        if one is None or two is None:
            return True
        if one[0] != two[0]:
            return False
        a, b = one[1], two[1]
        if isinstance(a, tuple):
            return all(abs(x - y) <= oracle.ON for x, y in zip(a, b))
        return abs(a - b) <= oracle.ON
    return same(keys[0], others[0]) and same(keys[1], others[1])


def squared(p, q):
    """The square of the distance between two points, exactly."""
    return (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2


def join(shared):
    """Shared stretches joined where one ends where the next starts."""
    joined = [list(pair) for pair in shared]
    changed = True
    while changed:
        changed = False
        for i, one in enumerate(joined):
            for j, two in enumerate(joined):
                if (i != j and one[1] is not None and two[0] is not None
                        and oracle.near(one[1], two[0])):
                    joined[i] = [one[0], two[1]]
                    del joined[j]
                    changed = True
                    break
            if changed:
                break
    return joined


def inside_shared(p, joined, first):
    """Whether a point lies within a shared stretch on the first shape's straight pieces."""
    for a, b in joined:
        for piece in first.pieces:
            if piece.straight is not None and oracle.on_straight(piece, p):
                start, way, _, _ = piece.straight
                along = [None if q is None else oracle.dot(oracle.minus(q, start), way)
                         for q in (a, b)]
                u = oracle.dot(oracle.minus(p, start), way)
                low = -oracle.INF if along[0] is None else along[0]
                high = oracle.INF if along[1] is None else along[1]
                low, high = min(low, high), max(low, high)
                if low - oracle.ON <= u <= high + oracle.ON and on_line(piece, p):
                    return True
    return False


def on_line(piece, p):
    """Whether a point lies on a straight piece's carrier."""
    start, way, _, _ = piece.straight
    return abs(oracle.cross(oracle.minus(p, start), way)) <= oracle.ON


# The families of base cases, each given a frame.

def lines(f):
    """Base cases of lines, rays and endless lines."""
    yield "line crosses line", segment(f, (0, 0), (4, 4)), segment(f, (0, 4), (4, 0))
    yield "lines apart", segment(f, (0, 0), (4, 0)), segment(f, (0, 1), (4, 2))
    yield "lines would cross past an end", segment(f, (0, 0), (1, 1)), segment(f, (3, 0), (4, -1))
    yield "parallel lines", segment(f, (0, 0), (4, 0)), segment(f, (0, 1), (4, 1))
    yield "collinear lines share a stretch", segment(f, (0, 0), (4, 0)), segment(f, (2, 0), (6, 0))
    yield "collinear lines, one within", segment(f, (0, 0), (8, 0)), segment(f, (2, 0), (5, 0))
    yield "collinear lines, against", segment(f, (0, 0), (4, 0)), segment(f, (6, 0), (2, 0))
    yield "collinear lines end to end", segment(f, (0, 0), (4, 0)), segment(f, (4, 0), (7, 0))
    yield "collinear lines apart", segment(f, (0, 0), (2, 0)), segment(f, (3, 0), (5, 0))
    yield "line ends on a line", segment(f, (0, 0), (4, 0)), segment(f, (1, 0), (1, 3))
    yield "lines meet at their ends", segment(f, (0, 0), (2, 2)), segment(f, (2, 2), (4, 0))
    yield "ray crosses line", ray(f, (0, 0), (1, 2)), segment(f, (-3, 4), (5, 4))
    yield "ray points away", ray(f, (0, 0), (-1, -2)), segment(f, (-3, 4), (5, 4))
    yield "xline crosses line", xline(f, (0, 0), (1, 2)), segment(f, (-3, -4), (5, -4))
    yield "xline along a line", xline(f, (0, 0), (2, 0)), segment(f, (1, 0), (3, 0))
    yield "ray along an xline", ray(f, (1, 1), (1, 1)), xline(f, (0, 0), (3, 3))


def circles(f):
    """Base cases of lines and circles with circular arcs."""
    circle = arc(f, (0, 0), 5, (1, 0), (1, 0), full=True)
    yield "line crosses circle", segment(f, (-8, 3), (8, 3)), circle
    yield "line touches circle", segment(f, (-8, 5), (8, 5)), circle
    yield "line misses circle", segment(f, (-8, 6), (8, 6)), circle
    yield "line ends on circle", segment(f, (0, 0), (3, 4)), circle
    yield "line from circle to circle", segment(f, (-3, 4), (4, 3)), circle
    upper = arc(f, (0, 0), 5, (1, 0), (-1, 0))
    yield "line crosses half circle", segment(f, (-8, 3), (8, 3)), upper
    yield "line misses half circle", segment(f, (-8, -3), (8, -3)), upper
    yield "line through arc's ends", segment(f, (-8, 0), (8, 0)), upper
    yield "xline crosses arc", xline(f, (0, 4), (1, 0)), \
        arc(f, (0, 0), 5, (F(4, 5), F(-3, 5)), (F(-4, 5), F(-3, 5)))
    yield "ray from inside circle", ray(f, (0, 0), (3, 4)), circle
    yield "circles cross", circle, arc(f, (6, 0), 5, (1, 0), (1, 0), full=True)
    yield "circles touch outside", arc(f, (0, 0), 3, (1, 0), (1, 0), full=True), \
        arc(f, (5, 0), 2, (1, 0), (1, 0), full=True)
    yield "circles touch inside", circle, arc(f, (3, 0), 2, (1, 0), (1, 0), full=True)
    yield "circles apart", circle, arc(f, (12, 0), 5, (1, 0), (1, 0), full=True)
    yield "concentric circles", circle, arc(f, (0, 0), 3, (1, 0), (1, 0), full=True)
    yield "circle and its arc", circle, upper
    yield "arcs of one circle overlap", upper, arc(f, (0, 0), 5, (0, 1), (0, -1))
    yield "arcs of one circle, against", upper, \
        arc(f, (0, 0), 5, (0, -1), (0, 1), counterclockwise=False)
    yield "arcs of one circle meet at an end", upper, arc(f, (0, 0), 5, (-1, 0), (0, -1))
    yield "arcs of one circle apart", arc(f, (0, 0), 5, (F(3, 5), F(4, 5)), (F(-3, 5), F(4, 5))), \
        arc(f, (0, 0), 5, (F(-3, 5), F(-4, 5)), (F(3, 5), F(-4, 5)))
    yield "arcs of one circle overlap twice", \
        arc(f, (0, 0), 5, (F(4, 5), F(3, 5)), (F(4, 5), F(-3, 5))), \
        arc(f, (0, 0), 5, (F(-4, 5), F(-3, 5)), (F(-4, 5), F(3, 5)))


def ellipses(f):
    """Base cases of ellipses with lines, circles and each other."""
    wide = ellipse(f, (0, 0), (5, 0), (0, 3), (1, 0), (1, 0), full=True)
    tall = ellipse(f, (0, 0), (3, 0), (0, 5), (1, 0), (1, 0), full=True)
    yield "ellipses cross four times", wide, tall
    yield "line touches ellipse", segment(f, (-8, 3), (8, 3)), wide
    yield "line crosses ellipse", segment(f, (-8, 1), (8, 1)), wide
    yield "circle crosses ellipse", wide, arc(f, (0, 0), 4, (1, 0), (1, 0), full=True)
    yield "circle touches ellipse", wide, arc(f, (0, 0), 3, (1, 0), (1, 0), full=True)
    sheared = ellipse(f, (1, 1), (4, 1), (1, 3), (1, 0), (1, 0), full=True)
    yield "sheared ellipse and line", sheared, segment(f, (-8, 2), (9, 2))
    yield "sheared ellipse and circle", sheared, arc(f, (1, 1), 3, (1, 0), (1, 0), full=True)
    half = ellipse(f, (0, 0), (5, 0), (0, 3), (1, 0), (-1, 0))
    yield "half ellipse and line", half, segment(f, (-8, -1), (8, 1))
    yield "ellipse and its half", wide, half


def polylines(f):
    """Base cases of polylines with and without bulges."""
    square = polyline(f, [((0, 0), 0), ((4, 0), 0), ((4, 4), 0), ((0, 4), 0)], closed=True)
    yield "square and line", square, segment(f, (-1, 1), (5, 3))
    yield "line along a square's side", square, segment(f, (1, 0), (6, 0))
    yield "line through a square's corner", square, segment(f, (-1, -1), (5, 5))
    bulged = polyline(f, [((0, 0), 0), ((4, 0), 1), ((4, 4), 0), ((0, 4), F(1, 2))], closed=True)
    yield "bulged square and line", bulged, segment(f, (-1, 2), (8, 2))
    yield "bulge on its circle", bulged, arc(f, (4, 2), 2, (1, 0), (1, 0), full=True)
    yield "bulged square and circle", bulged, arc(f, (2, 2), 3, (1, 0), (1, 0), full=True)
    yield "open polyline and line", polyline(f, [((0, 0), 0), ((2, 2), F(1, 3)), ((5, 2), 0)]), \
        segment(f, (0, 2), (6, 2))
    yield "lone vertex on a line", polyline(f, [((2, 0), 0)]), segment(f, (0, 0), (4, 0))


def beziers(f):
    """Base cases of Bézier curves and splines."""
    hump = [(0, 0), (1, 3), (3, 2), (4, 0)]
    wave = [(0, 0), (1, 3), (3, -2), (4, 1)]
    arch = [(0, 0), (1, 2), (2, 0)]
    yield "cubic crosses line", bezier(f, hump), segment(f, (-1, 1), (5, 1))
    yield "quadratic touches line", bezier(f, arch), segment(f, (-1, 1), (3, 1))
    yield "quadratic ends on line", bezier(f, arch), segment(f, (-1, 0), (3, 0))
    yield "quadratics touch", bezier(f, arch), bezier(f, [(0, 2), (1, 0), (2, 2)])
    yield "cubics cross", bezier(f, hump), bezier(f, wave)
    yield "cubic crosses circle", bezier(f, wave), arc(f, (2, 0), 2, (1, 0), (1, 0), full=True)
    yield "cubic crosses ellipse", bezier(f, hump), \
        ellipse(f, (2, 1), (3, 0), (0, 1), (1, 0), (1, 0), full=True)
    first, second = cut(hump, F(1, 2))
    yield "cubic and its first half", bezier(f, hump), bezier(f, first)
    yield "cubic and its second half", bezier(f, hump), bezier(f, second)
    left, _ = cut(hump, F(3, 4))
    _, right = cut(hump, F(1, 4))
    yield "two parts of a cubic overlap", bezier(f, right), bezier(f, left)
    yield "cubic lying on a line", bezier(f, [(1, 0), (2, 0), (3, 0), (5, 0)]), \
        segment(f, (0, 0), (4, 0))
    yield "spline crosses line at a knot", spline(f, [first, second]), segment(f, (2, -1), (2, 4))
    yield "spline and the cubic it is", spline(f, [first, second]), bezier(f, hump)
    yield "spline crosses circle", spline(f, [first, second]), \
        arc(f, (2, 0), 2, (1, 0), (1, 0), full=True)


def near(f):
    """Cases a fraction of the tolerance from a tangency or an end: the tolerance's rule."""
    t = f.tolerance
    circle = arc(f, (0, 0), 5, (1, 0), (1, 0), full=True)
    for shift, kind in ((t / 2, "Touch"), (-t / 2, "Touch"), (2 * t, None)):
        yield (f"line {float(shift / t)} tolerances from touching a circle",
               segment(f, (-8, 5 + shift), (8, 5 + shift)), circle,
               [(f.point((0, 5 + shift / 2)), kind)] if kind else [])
    for shift, kind in ((t / 2, "End"), (2 * t, None)):
        yield (f"line stopping {float(shift / t)} tolerances short of a circle",
               segment(f, (0, 0), (0, 5 - shift)), circle,
               [(f.point((0, 5 - shift / 2)), kind)] if kind else [])
    for shift, kind in ((t / 2, "Touch"), (2 * t, None)):
        yield (f"circles {float(shift / t)} tolerances from touching",
               arc(f, (0, 0), 3, (1, 0), (1, 0), full=True),
               arc(f, (5 + shift, 0), 2, (1, 0), (1, 0), full=True),
               [(f.point((3 + shift / 2, 0)), kind)] if kind else [])
    for shift, kind in ((t / 2, "Touch"), (2 * t, None)):
        yield (f"quadratic {float(shift / t)} tolerances from touching a line",
               bezier(f, [(0, 0), (1, 2), (2, 0)]), segment(f, (-1, 1 + shift), (3, 1 + shift)),
               [(f.point((1, 1 + shift / 2)), kind)] if kind else [])


def graph(frame, px, py, span):
    """A cubic Bézier curve through the point (x(s), y(s)) for s across a span, x and y
    polynomials of degree at most three from the constant up: its poles by Bernstein's basis,
    exactly."""
    lo, hi = F(span[0]), F(span[1])
    axes = []
    for poly in (px, py):
        power = [F(0)] * 4
        for k, c in enumerate(poly):
            for i, v in enumerate(oracle.power([lo, hi - lo], k)):
                power[i] += F(c) * v
        axes.append([sum(F(math.comb(i, j), math.comb(3, j)) * power[j] for j in range(i + 1))
                     for i in range(4)])
    return bezier(frame, list(zip(axes[0], axes[1])))


def unit(s):
    """The rational unit vector at twice the angle whose tangent is s."""
    s = F(s)
    return ((1 - s * s) / (1 + s * s), 2 * s / (1 + s * s))


def contact(tolerance, gap, order):
    """How far along a contact of an order the curves stay within the tolerance, twice over:
    where they part as gap times the distance along to the order."""
    return 2 * (float(tolerance) / gap) ** (1 / order)


def hard(f):
    """Cases that push a kernel: crossings at angles near none, secants a few tolerances from a
    tangency, contact of high order, crossings clustered closer than the tolerance, cusps,
    loops, seams, curves too small or too large, and stretches shared over a sliver. Each is
    sized from the frame's tolerance to stay clear of it by a factor of four: within it is one
    answer, beyond it another, never rounding's choice."""
    t = f.tolerance
    bits = -int(math.log2(t))
    circle = arc(f, (0, 0), 5, (1, 0), (1, 0), full=True)
    for k in (8, bits - 2):
        e = F(1, 2 ** k)
        yield (f"lines crossing at 2^-{k} radians", segment(f, (0, 0), (16, 0)),
               segment(f, (0, -8 * e), (16, 8 * e)), {})
    arch = bezier(f, [(0, 0), (1, 2), (2, 0)])
    wide = ellipse(f, (0, 0), (5, 0), (0, 3), (1, 0), (1, 0), full=True)
    for d, name in ((4 * t, "4 tolerances"), (1024 * t, "1024 tolerances")):
        yield (f"line crossing a circle {name} inside its tangent",
               segment(f, (-8, 5 - d), (8, 5 - d)), circle, {})
        yield (f"line crossing a quadratic {name} below its top",
               arch, segment(f, (-1, 1 - d), (3, 1 - d)), {})
        yield (f"line crossing an ellipse {name} inside its tangent",
               wide, segment(f, (-8, 3 - d), (8, 3 - d)), {})
        yield (f"circles crossing {name} past touching inside", circle,
               arc(f, (2 + d, 0), 3, (1, 0), (1, 0), full=True), {})
        yield (f"circles crossing {name} past touching outside",
               arc(f, (0, 0), 3, (1, 0), (1, 0), full=True),
               arc(f, (5 - d, 0), 2, (1, 0), (1, 0), full=True), {})
    flat = segment(f, (-2, 0), (2, 0))
    yield ("line along a cubic's inflection, contact of order three",
           graph(f, [0, 1], [0, 0, 0, 1], (-F(3, 2), F(3, 2))), flat,
           {"band": contact(t, 1, 3)})
    yield ("cubic crossing a line three times 2^-4 apart",
           graph(f, [0, 1], [0, -F(1, 256), 0, 1], (-F(3, 2), F(3, 2))), flat, {})
    # Crossings at 0 and +-sqrt(e), where the curves part by 3/8 e^1.5 halfway, far within the
    # tolerance: one touch. e as large as a double holds beside the frame's shift.
    e = F(1, 2 ** (bits + 4))
    yield ("cubic crossing a line three times closer than the tolerance parts them, one touch",
           graph(f, [0, 1], [0, -e, 0, 1], (-F(3, 2), F(3, 2))), flat,
           {"hits": [(f.point((0, 0)), "Touch", contact(t, 1, 3) + 2 * math.sqrt(e))]})
    yield ("cubics tangent where they cross",
           graph(f, [0, 1], [0, -1, 0, 1], (-F(3, 2), F(3, 2))),
           graph(f, [0, -1, 0, 1], [0, 1], (-F(3, 2), F(3, 2))), {"band": contact(t, 2, 3)})
    yield ("ellipse and circle osculating, contact of order four, at the ellipse's seam",
           ellipse(f, (0, 0), (4, 0), (0, 2), (1, 0), (1, 0), full=True),
           arc(f, (3, 0), 1, (1, 0), (1, 0), full=True), {"band": contact(t, 3 / 32, 4)})
    sliver = F(1, 1024)
    yield ("circles touching inside, radii 2^-10 apart, at their seams", circle,
           arc(f, (sliver, 0), 5 - sliver, (1, 0), (1, 0), full=True),
           {"band": contact(t, float(sliver) / 50, 2)})
    cusp = graph(f, [0, 0, 1], [0, 0, 0, 1], (-F(3, 2), F(3, 2)))
    yield "line crossing a cusped cubic", cusp, segment(f, (F(1, 4), -1), (F(1, 4), 1)), {}
    yield ("line along a cusp's tangent", cusp, segment(f, (-1, 0), (3, 0)),
           {"band": 4 * float(t) ** (2 / 3)})
    yield ("line across a cusp", cusp, segment(f, (0, -1), (0, 1)), {"band": 4 * float(t)})
    loop = graph(f, [0, 0, 1], [0, -1, 0, 1], (-F(3, 2), F(3, 2)))
    yield "line through a loop's double point", loop, segment(f, (1, -2), (1, 2)), {}
    yield "line through a loop's double point and its nose", loop, segment(f, (-1, 0), (3, 0)), {}
    yield "line through a circle's seam", circle, segment(f, (-1, -6), (6, 1)), {}
    yield ("line through an ellipse's seam", wide,
           segment(f, (F(25, 4), -F(3, 4)), (-F(5, 4), F(15, 4))), {})
    hump = [(0, 0), (1, 3), (3, 2), (4, 0)]
    first, second = cut(hump, F(1, 2))
    yield ("line tangent at a spline's knot", spline(f, [first, second]),
           segment(f, (-1, F(19, 8)), (5, F(11, 8))), {})
    square = polyline(f, [((0, 0), 0), ((4, 0), 0), ((4, 4), 0), ((0, 4), 0)], closed=True)
    yield ("line touching a square's corner from outside", square, segment(f, (2, 6), (6, 2)), {})
    across = segment(f, (0, 0), (4, 0))
    yield ("segment shorter than the tolerance across a line", across,
           segment(f, (2, -t / 8), (2, t / 8)), {})
    yield "segment of no length on a line", across, segment(f, (2, 0), (2, 0)), {}
    big = 2 ** 20
    yield ("arc of radius 2^20 crossed near its middle",
           arc(f, (0, -big), big, unit(-F(1, 2 ** 18))[::-1], unit(F(1, 2 ** 18))[::-1],
               counterclockwise=False),
           segment(f, (F(1, 2), -1), (F(1, 2), 1)), {})
    yield ("arc of sweep 2^-11 crossed", arc(f, (0, 0), 5, (1, 0), unit(F(1, 4096))),
           segment(f, (4, F(4, 4096)), (6, F(6, 4096))), {})
    thin = ellipse(f, (0, 0), (8, 0), (0, F(1, 2 ** 16)), (1, 0), (1, 0), full=True)
    yield "ellipse 2^19 times as long as wide, crossed", thin, segment(f, (1, -1), (1, 1)), {}
    yield ("ellipse 2^19 times as long as wide, along its axis", thin,
           segment(f, (-10, 0), (10, 0)), {})
    yield ("lines sharing 16 tolerances", across, segment(f, (4 - 16 * t, 0), (8, 0)), {})
    yield ("arcs of one circle sharing a sliver", arc(f, (0, 0), 5, unit(F(1, 4096)), (-1, 0)),
           arc(f, (0, 0), 5, (0, -1), unit(F(2, 4096))), {})


FAMILIES = (lines, circles, ellipses, polylines, beziers)


def frames():
    """Each frame a base case is written in, with a name for where it sits."""
    for dx, dy, tolerance in SCALES:
        for turn, mirror, name in ((0, False, "as built"), (1, False, "turned"),
                                   (0, True, "mirrored")):
            yield f"{name} at {float(dx):g}", Frame(dx, dy, turn, mirror, tolerance)


def slacked(hits, tolerance, band=None):
    """Each hit with how far off geoxy's may be: a crossing or an end anywhere both curves are
    within the tolerance, the tolerance over the sine of half the angle they meet at; a touch
    within its band, given or of a contact whose curvatures differ by up to 2/SIZE."""
    out = []
    for p, kind, sine in hits:
        if kind == "Touch":
            slack = band if band is not None else math.sqrt(8 * float(tolerance) * SIZE)
        elif sine:
            slack = float(tolerance) / math.sin(math.asin(sine) / 2)
        else:
            slack = float(tolerance)
        out.append((p, kind, slack))
    return out


def emit(name, tolerance, shapes, hits, shared):
    """Prints one case in the suite's text, unless an arc's ends, read from float angles, may
    sit too far from the oracle's near an answer: then the case is refused."""
    for shape in shapes:
        near_end = any(squared(p, e) <= (shape.drift + 4 * float(tolerance)) ** 2
                       for p, _, _ in hits for e in shape.ends)
        if near_end and shape.drift > float(tolerance) / 64:
            raise ValueError(f"{name}: an arc's float angles move its ends {shape.drift:.3g}")
    print(f"case geoxy {name}")
    print(f"tolerance {r(tolerance)}")
    print(f"first {shapes[0].text}")
    print(f"second {shapes[1].text}")
    for p, kind, slack in hits:
        print(f"hit {answer(p)} {kind} {slack!r}")
    for p, q in shared:
        print(f"shared {'inf inf' if p is None else answer(p)} "
              f"{'inf inf' if q is None else answer(q)}")
    print("end")


def main():
    """Prints the whole suite."""
    print("# geoxy's own suite of where two curves meet, forged by tools/suites/forge.py. MIT.")
    for where, frame in frames():
        t = frame.tolerance
        for family in FAMILIES:
            for title, first, second in family(frame):
                hits, shared = expected(first, second, t)
                emit(f"{family.__name__}: {title}, {where}", t, (first, second),
                     slacked(hits, t), shared)
        for title, first, second, hits in near(frame):
            emit(f"near: {title}, {where}", t, (first, second),
                 slacked([(p, kind, None) for p, kind in hits], t), [])
        for title, first, second, options in hard(frame):
            hits, shared = expected(first, second, t)
            if "hits" in options:
                hits = [(p, kind, None) for p, kind, _ in options["hits"]]
            out = slacked(hits, t, options.get("band"))
            if "hits" in options:
                out = [(p, kind, given or slack) for (p, kind, slack), (_, _, given)
                       in zip(out, options["hits"])]
            emit(f"hard: {title}, {where}", t, (first, second), out, shared)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--check":
        made = io.StringIO()
        with contextlib.redirect_stdout(made):
            main()
        with open(sys.argv[2], encoding="utf-8") as kept:
            fresh = kept.read() == made.getvalue()
        print(f"{sys.argv[2]}: {'as forged' if fresh else 'not what the forge makes'}")
        sys.exit(0 if fresh else 1)
    main()
