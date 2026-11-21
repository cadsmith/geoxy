#!/usr/bin/env python3
"""Forges stress cases of lines and circular arcs, each with where they truly meet.

    python3 tools/reference/stress.py > tests/suites/stress/cases/stress.txt
    python3 tools/reference/stress.py --check tests/suites/stress/cases/stress.txt

Each case is two curves in the text geoxy's `FromText` reads, drawn to be hard for double
arithmetic: lines near parallel, lines grazing circles, circles nearly touching, radii past a
million, shapes a million or a billion units from the origin. The numbers written are the doubles
every library is given; where those curves truly meet is then worked out from them at 120
significant digits, as tools/reference/truth.py does, and each point written as the double nearest
it: the one answer that is right to the last bit.

Each point is a crossing, and is written with how far it moves when the curves move by a few steps
of a double, as they do spelled another way: much, where they cross at a slight angle.

A case is kept only where nothing is left to a tolerance: every crossing of the carriers lies well
within both curves or well off one of them, no end of either lies near the other, carriers that
cross twice part by ten thousand tolerances between their crossings, and nothing touches. So any
library that answers with a tolerance, whatever it is, owes the same points. The cases are seeded:
every run forges the same file.

mpmath is needed: python3 -m venv venv && venv/bin/pip install mpmath.
"""

import importlib.util
import math
import pathlib
import random
import sys

import mpmath
from mpmath import mp, mpf

sys.dont_write_bytecode = True
HERE = pathlib.Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("truth", HERE / "truth.py")
truth = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(truth)

TOLERANCE = 1e-12  # Of the case's size: what every case is met within.
CLEAR = 1e4  # How many tolerances apart things must lie to owe nothing to a tolerance.
MARGIN = mpf(10) ** -3  # How far inside, or outside, a crossing must lie: of a line, in radians.
EACH = 60  # Cases of each family.


def line(start, end):
    """A line's text."""
    return f"line {start[0]!r} {start[1]!r} {end[0]!r} {end[1]!r}"


def arc(center, radius, start, sweep, degrees):
    """An arc's text, its angles in degrees or in radians."""
    if degrees:
        start, sweep = f"{math.degrees(start)!r}d", f"{math.degrees(sweep)!r}d"
    else:
        start, sweep = repr(start), repr(sweep)
    return f"arc {center[0]!r} {center[1]!r} {radius!r} {start} {sweep}"


def about(rng, angle, least, most):
    """A start and a sweep whose arc holds an angle well within it, turning either way."""
    span = rng.uniform(least, most)
    before = rng.uniform(0.2, 0.8) * span
    if rng.random() < 0.5:
        return angle - before, span
    return angle + before, -span


def whole(rng):
    """A whole turn from somewhere."""
    return rng.uniform(-math.pi, math.pi), 2 * math.pi


def moved(point, offset):
    """A point carried off by an offset, rounded as a double is."""
    return point[0] + offset[0], point[1] + offset[1]


def offsets(rng, reach):
    """Where a case lies: a few units from the origin, or a thousand, a million, a billion off."""
    far = rng.choice(reach) or rng.uniform(1, 9)
    turn = rng.uniform(0, 2 * math.pi)
    return far * math.cos(turn), far * math.sin(turn)


def crossing_lines(rng, angle, reach):
    """Two lines crossing at an angle, the crossing well within both."""
    scale = 10.0 ** rng.uniform(-2, 3)
    offset = offsets(rng, reach)
    turn = rng.uniform(0, 2 * math.pi)
    lines = []
    for way in (turn, turn + angle):
        back, on = scale * rng.uniform(0.3, 2), scale * rng.uniform(0.3, 2)
        lines.append(line(moved((-back * math.cos(way), -back * math.sin(way)), offset),
                          moved((on * math.cos(way), on * math.sin(way)), offset)))
    return lines


def lines_plain(rng):
    """Lines crossing at a plain angle, near the origin and far from it."""
    return crossing_lines(rng, rng.uniform(0.3, math.pi - 0.3), (0, 0, 1e3, 1e6, 1e9))


def lines_narrow(rng):
    """Lines crossing at a slight angle, down to a ten millionth of a radian."""
    return crossing_lines(rng, 10.0 ** rng.uniform(-7, -2), (0, 0, 1e3))


def through_circle(rng, depth, radius, reach):
    """A line cutting a circle to a depth, of its radius, and an arc of it about one crossing."""
    offset = offsets(rng, reach)
    turn = rng.uniform(0, 2 * math.pi)
    off = radius * (1 - depth)
    half = radius * math.sqrt(max(0.0, 1 - (1 - depth) ** 2))
    across = (math.cos(turn), math.sin(turn))
    along = (-across[1], across[0])
    foot = (off * across[0], off * across[1])
    # The line runs through both crossings, or stops between them, well clear of the circle.
    back = half * rng.uniform(1.5, 4)
    on = half * rng.uniform(1.5, 4) if rng.random() < 0.6 else half * rng.uniform(-0.6, 0.6)
    start = (foot[0] - back * along[0], foot[1] - back * along[1])
    end = (foot[0] + on * along[0], foot[1] + on * along[1])
    side = rng.choice((-1, 1))
    angle = math.atan2(foot[1] + side * half * along[1], foot[0] + side * half * along[0])
    least = 0.4 * math.sqrt(depth)
    first, sweep = whole(rng) if rng.random() < 0.3 else about(rng, angle, least, 5)
    return [line(moved(start, offset), moved(end, offset)),
            arc(moved((0.0, 0.0), offset), radius, first, sweep, rng.random() < 0.5)]


def circle_plain(rng):
    """A line cutting a circle deep, near the origin and far from it."""
    reach = (0, 0, 1e3, 1e6, 1e9)
    return through_circle(rng, rng.uniform(0.1, 1.6), 10.0 ** rng.uniform(-2, 3), reach)


def circle_grazed(rng):
    """A line grazing a circle, cutting in by a hundred millionth of its radius at the least."""
    return through_circle(rng, 10.0 ** rng.uniform(-8, -3), 10.0 ** rng.uniform(0, 3), (0,))


def circle_huge(rng):
    """A line cutting a circle of a radius up to ten billion, near where the line lies."""
    return through_circle(rng, 10.0 ** rng.uniform(-8, -1), 10.0 ** rng.uniform(6, 10), (0,))


def two_circles(rng, radius, other, apart, reach):
    """Two circles a distance apart, and arcs of them about their crossings."""
    offset = offsets(rng, reach)
    turn = rng.uniform(0, 2 * math.pi)
    toward = (math.cos(turn), math.sin(turn))
    centre = (apart * toward[0], apart * toward[1])
    along = (apart * apart + radius * radius - other * other) / (2 * apart)
    height = math.sqrt(max(0.0, radius * radius - along * along))
    side = rng.choice((-1, 1))
    point = (along * toward[0] - side * height * toward[1],
             along * toward[1] + side * height * toward[0])
    arcs = []
    for middle, size in (((0.0, 0.0), radius), (centre, other)):
        angle = math.atan2(point[1] - middle[1], point[0] - middle[0])
        least = 0.4 * height / size
        first, sweep = whole(rng) if rng.random() < 0.3 else about(rng, angle, least, 5)
        arcs.append(arc(moved(middle, offset), size, first, sweep, rng.random() < 0.5))
    return arcs


def circles_plain(rng):
    """Two circles crossing deep, near the origin and far from it."""
    radius, other = 10.0 ** rng.uniform(-1, 3), 10.0 ** rng.uniform(-1, 3)
    low, high = abs(radius - other), radius + other
    apart = low + (high - low) * rng.uniform(0.15, 0.85)
    return two_circles(rng, radius, other, apart, (0, 0, 1e3, 1e6, 1e9))


def circles_grazed(rng):
    """Two circles nearly touching, from outside or from inside, by a hundred millionth of the
    smaller at the least."""
    radius, other = 10.0 ** rng.uniform(0, 2), 10.0 ** rng.uniform(0, 2)
    depth = 10.0 ** rng.uniform(-8, -3) * min(radius, other)
    inside = rng.random() < 0.5 and abs(radius - other) > 100 * depth
    apart = abs(radius - other) + depth if inside else radius + other - depth
    return two_circles(rng, radius, other, apart, (0,))


FAMILIES = (
    ("lines", "lines crossing at a plain angle", lines_plain),
    ("narrow lines", "lines crossing at a slight angle", lines_narrow),
    ("line and circle", "a line cutting a circle deep", circle_plain),
    ("grazed circle", "a line grazing a circle", circle_grazed),
    ("huge circle", "a line cutting a circle of a huge radius", circle_huge),
    ("circles", "two circles crossing deep", circles_plain),
    ("grazing circles", "two circles nearly touching", circles_grazed),
)


def placed(curve, point, along):
    """Whether a point of a curve's carrier lies well within the curve, 1, well off it, -1, or
    near an end of it, 0."""
    kind, values = curve
    if kind == "line":
        if MARGIN <= along <= 1 - MARGIN:
            return 1
        return -1 if along <= -MARGIN or along >= 1 + MARGIN else 0
    cx, cy, _, start, sweep = values
    span = abs(sweep)
    if span >= 2 * mp.pi:
        return 1
    angle = mpmath.atan2(point[1] - cy, point[0] - cx)
    turned = (angle - start) % (2 * mp.pi) if sweep > 0 else (start - angle) % (2 * mp.pi)
    if MARGIN <= turned <= span - MARGIN:
        return 1
    return -1 if span + MARGIN <= turned <= 2 * mp.pi - MARGIN else 0


def ends(curve):
    """A curve's ends: a line's two, an arc's two, none for a whole circle."""
    kind, values = curve
    if kind == "line":
        return [(values[0], values[1]), (values[2], values[3])]
    cx, cy, r, start, sweep = values
    if abs(sweep) >= 2 * mp.pi:
        return []
    turns = (start, start + sweep)
    return [(cx + r * mpmath.cos(angle), cy + r * mpmath.sin(angle)) for angle in turns]


def off(curve, point):
    """How far a point is from a curve's carrier."""
    kind, values = curve
    if kind == "line":
        x0, y0, x1, y1 = values
        wx, wy = x1 - x0, y1 - y0
        return abs((point[0] - x0) * wy - (point[1] - y0) * wx) / mpmath.sqrt(wx * wx + wy * wy)
    cx, cy, r = values[:3]
    return abs(mpmath.sqrt((point[0] - cx) ** 2 + (point[1] - cy) ** 2) - r)


def heading(curve, point):
    """The way a curve runs at a point of it: a line's own, a circle's across its radius."""
    kind, values = curve
    if kind == "line":
        return values[2] - values[0], values[3] - values[1]
    return -(point[1] - values[1]), point[0] - values[0]


def slack(first, second, point):
    """How far a crossing moves when its curves move by a few steps of a double, as when they are
    spelled another way: 64 steps at its size, over the sine of the angle the curves cross at."""
    (ux, uy), (vx, vy) = heading(first, point), heading(second, point)
    sine = abs(ux * vy - uy * vx) / (mpmath.hypot(ux, uy) * mpmath.hypot(vx, vy))
    size = max(1.0, abs(float(point[0])), abs(float(point[1])))
    return float(64 * sys.float_info.epsilon * size / sine)


def parting(first, second):
    """How far two carriers that cross twice part between their crossings: a line's depth into a
    circle; of two circles, how far one reaches into the other, or out of it."""
    if first[0] != "arc" or second[0] != "arc":
        straight, round_ = (first, second) if first[0] != "arc" else (second, first)
        cx, cy, r = round_[1][:3]
        x0, y0, x1, y1 = straight[1]
        wx, wy = x1 - x0, y1 - y0
        return r - abs((cx - x0) * wy - (cy - y0) * wx) / mpmath.sqrt(wx * wx + wy * wy)
    (ax, ay, ra), (bx, by, rb) = first[1][:3], second[1][:3]
    apart = mpmath.hypot(bx - ax, by - ay)
    return min(ra + rb - apart, apart - abs(ra - rb))


def answer(texts):
    """The points two curves truly meet at, each the double nearest it, with how far it moves when
    the curves do; None where the case owes something to a tolerance, and is not kept."""
    first, second = (truth.parsed(text) for text in texts)
    found = truth.crossings(first, second)
    if not found or any(touch for _, _, touch in found):
        return None
    size = max([1.0] + [abs(float(value)) for curve in (first, second) for value in curve[1][:2]])
    clear = CLEAR * TOLERANCE * size
    points = []
    for point, (along, beyond), _ in found:
        size = max(size, abs(float(point[0])), abs(float(point[1])))
        other = beyond if first[0] != "arc" else along
        sides = (placed(first, point, along), placed(second, point, other))
        if 0 in sides:
            return None
        if sides == (1, 1):
            points.append(point)
    clear = CLEAR * TOLERANCE * size
    if len(found) == 2 and parting(first, second) < clear:
        return None
    for curve, other in ((first, second), (second, first)):
        if any(off(other, end) < clear for end in ends(curve)):
            return None
    if not points:
        return None
    return [(truth.nearest(x), truth.nearest(y), slack(first, second, (x, y))) for x, y in points]


def forged():
    """Every case's lines of text, family by family."""
    out = ["# Stress cases of lines and arcs, each with where they truly meet, worked out at 120 "
           "digits and",
           "# rounded once, forged by tools/reference/stress.py. MIT."]
    for family, title, pose in FAMILIES:
        rng = random.Random(f"geoxy stress {family}")
        made = 0
        while made < EACH:
            texts = pose(rng)
            if rng.random() < 0.5:
                texts.reverse()
            points = answer(texts)
            if points is None:
                continue
            made += 1
            out.append(f"case stress {family}: {title}, {made}")
            out.append(f"tolerance {TOLERANCE!r} relative")
            out.append(f"first {texts[0]}")
            out.append(f"second {texts[1]}")
            out.extend(f"hit {x!r} {y!r} Cross {loose!r}" for x, y, loose in points)
            out.append("end")
    return out


def main():
    lines = forged()
    if len(sys.argv) == 3 and sys.argv[1] == "--check":
        kept = open(sys.argv[2], encoding="utf-8").read().splitlines()
        if kept != lines:
            sys.exit(f"{sys.argv[2]} is not what tools/reference/stress.py forges")
        return
    print("\n".join(lines))


if __name__ == "__main__":
    main()
