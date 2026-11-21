#!/usr/bin/env python3
"""Checks, at 120 significant digits, that every point a case file says two curves meet at lies on
both of them: for every kind of curve, and sharing nothing with geoxy or with the forge's oracle.

    python3 tools/reference/sound.py CASES...

Each case's curves are read as the doubles every library is given. For each hit and each end of a
shared stretch the case expects, the distance from the point to each curve is worked out in mpmath:
in closed form to a line, a ray, an endless line and a circular arc; to an ellipse, a Bézier curve
and each piece of a spline as the least over the places where the curve runs square to the way to
the point, each found between two of 96 places a piece where that squareness changes sign and
narrowed to 100 digits, and over the piece's ends; to a polyline as the least over its segments,
a bulge's arc built from its ends and its bulge. A shared stretch that never ends, as a ray's
along an endless line, has no point at that end to check.

What is worked out is never less than the true distance: a foot not found leaves a distance too
large, and a point is then told off, never let through. So a point passed lies on both curves:

    exact    within two steps of a double, at the case's size, of both curves: what a point that
             is a true meeting, rounded once to the nearest double, cannot miss;
    near     within the case's tolerance of both, or within the hit's own slack: an answer the
             tolerance rules give, where the curves come within it and do not truly meet;
    off      farther than both from a curve: the case is wrong, or this check missed a foot.

The exit status is 1 if any point is off. This tells that each answer written is a place of both
curves. It does not tell that no meeting is missing from a case: for lines and arcs truth.py does,
for the other kinds the forge's oracle, by exact arithmetic.

mpmath is needed: python3 -m venv venv && venv/bin/pip install mpmath.
"""

import sys

from mpmath import mp, mpf

mp.dps = 120
PLACES = 96
FINE = mpf(10) ** -100


def angle(word):
    """An angle from its word, in radians: degrees where a d follows it."""
    if word.endswith("d"):
        return mpf(float(word[:-1])) * mp.pi / 180
    return mpf(float(word))


def parsed(text):
    """A curve from its text: its kind and what it is made of, every number the double written."""
    words = text.split()
    kind, rest = words[0], words[1:]
    if kind in ("line", "ray", "xline", "quadratic", "cubic"):
        return kind, [mpf(float(word)) for word in rest]
    if kind == "arc":
        return kind, [mpf(float(word)) for word in rest[:3]] + [angle(rest[3]), angle(rest[4])]
    if kind == "ellipse":
        return kind, [mpf(float(word)) for word in rest[:6]] + [angle(rest[6]), angle(rest[7])]
    if kind == "polyline":
        count = int(rest[1])
        numbers = [mpf(float(word)) for word in rest[2:2 + 3 * count]]
        return kind, (rest[0] == "closed", [numbers[3 * i:3 * i + 3] for i in range(count)])
    if kind == "spline":
        degree, knots = int(rest[0]), int(rest[1])
        ties = [mpf(float(word)) for word in rest[2:2 + knots]]
        poles = int(rest[2 + knots])
        flat = [mpf(float(word)) for word in rest[3 + knots:3 + knots + 2 * poles]]
        after = 3 + knots + 2 * poles
        weights = [mpf(float(word)) for word in rest[after + 1:after + 1 + int(rest[after])]]
        return kind, (degree, ties, [flat[2 * i:2 * i + 2] for i in range(poles)], weights)
    raise ValueError(f"no such curve: {kind}")


def away(point, other):
    """The distance between two points."""
    return mp.hypot(point[0] - other[0], point[1] - other[1])


def from_run(point, start, way, stretch):
    """The distance from a point to a straight curve, by its foot kept within the stretch."""
    square = way[0] * way[0] + way[1] * way[1]
    if square == 0:
        return away(point, start)
    along = ((point[0] - start[0]) * way[0] + (point[1] - start[1]) * way[1]) / square
    along = min(max(along, stretch[0]), stretch[1])
    return away(point, (start[0] + along * way[0], start[1] + along * way[1]))


def from_arc(point, values):
    """The distance from a point to a circular arc: to its circle where it sweeps over the way to
    the point, to the nearer of its ends where it does not."""
    cx, cy, radius, start, sweep = values
    ends = [(cx + radius * mp.cos(turn), cy + radius * mp.sin(turn))
            for turn in (start, start + sweep)]
    least = min(away(point, end) for end in ends)
    if (point[0], point[1]) == (cx, cy):
        return least if sweep == 0 else min(least, abs(radius))
    toward = mp.atan2(point[1] - cy, point[0] - cx)
    turned = (toward - start) % (2 * mp.pi) if sweep > 0 else (start - toward) % (2 * mp.pi)
    if turned <= abs(sweep):
        least = min(least, abs(away(point, (cx, cy)) - abs(radius)))
    return least


def narrowed(square, low, high):
    """Where a function that changes sign between two places is nothing, by bisection."""
    below = square(low)
    while high - low > FINE:
        middle = (low + high) / 2
        value = square(middle)
        if value == 0:
            return middle
        if (value > 0) == (below > 0):
            low, below = middle, value
        else:
            high = middle
    return (low + high) / 2


def from_piece(point, at, slope, span=(mpf(0), mpf(1))):
    """The distance from a point to a piece of a curve given by its point and its slope at a
    parameter: the least over its ends and over where it runs square to the way to the point."""

    def square(along):
        here, way = at(along), slope(along)
        return (here[0] - point[0]) * way[0] + (here[1] - point[1]) * way[1]

    places = [span[0] + (span[1] - span[0]) * step / PLACES for step in range(PLACES + 1)]
    values = [square(place) for place in places]
    least = min(away(point, at(span[0])), away(point, at(span[1])))
    for index in range(PLACES):
        if values[index] == 0:
            least = min(least, away(point, at(places[index])))
        elif (values[index] > 0) != (values[index + 1] > 0) and values[index + 1] != 0:
            foot = narrowed(square, places[index], places[index + 1])
            least = min(least, away(point, at(foot)))
    return least


def from_ellipse(point, values):
    """The distance from a point to an ellipse's arc, a turn of its sweep at a time."""
    xx, xy, yx, yy, ox, oy, start, sweep = values

    def at(turn):
        cos, sin = mp.cos(turn), mp.sin(turn)
        return ox + xx * cos + yx * sin, oy + xy * cos + yy * sin

    def slope(turn):
        cos, sin = mp.cos(turn), mp.sin(turn)
        return -xx * sin + yx * cos, -xy * sin + yy * cos

    turns = max(1, int(mp.ceil(abs(sweep) / (2 * mp.pi))))
    cuts = [start + sweep * step / turns for step in range(turns + 1)]
    return min(from_piece(point, at, slope, (cuts[i], cuts[i + 1])) for i in range(min(turns, 2)))


def casteljau(poles, along):
    """A rational Bézier piece's point and slope at a parameter, its poles homogeneous."""
    rows = [list(pole) for pole in poles]
    degree = len(rows) - 1
    last = None
    for step in range(degree):
        if step == degree - 1:
            last = [[degree * (rows[1][k] - rows[0][k]) for k in range(3)]]
        rows = [[(1 - along) * rows[i][k] + along * rows[i + 1][k] for k in range(3)]
                for i in range(len(rows) - 1)]
    x, y, w = rows[0]
    dx, dy, dw = last[0]
    here = (x / w, y / w)
    return here, ((dx - here[0] * dw) / w, (dy - here[1] * dw) / w)


def from_bezier(point, poles):
    """The distance from a point to a rational Bézier piece."""
    return from_piece(point, lambda along: casteljau(poles, along)[0],
                      lambda along: casteljau(poles, along)[1])


def blossom(spline, span, values):
    """A spline's blossom over a span, by de Boor's triangle, each round at its own value."""
    degree, knots, poles, weights = spline
    cells = []
    for index in range(span - degree, span + 1):
        weight = weights[index] if weights else mpf(1)
        cells.append([poles[index][0] * weight, poles[index][1] * weight, weight])
    for round_ in range(1, degree + 1):
        for index in range(degree, round_ - 1, -1):
            low = knots[span - degree + index]
            high = knots[span + 1 + index - round_]
            part = (values[round_ - 1] - low) / (high - low)
            cells[index] = [(1 - part) * cells[index - 1][k] + part * cells[index][k]
                            for k in range(3)]
    return cells[degree]


def from_spline(point, spline):
    """The distance from a point to a spline: the least over its spans, each as a Bézier piece."""
    degree, knots, poles, _ = spline
    least = None
    for span in range(degree, len(poles)):
        if knots[span] < knots[span + 1]:
            piece = [blossom(spline, span, [knots[span]] * (degree - i) + [knots[span + 1]] * i)
                     for i in range(degree + 1)]
            distance = from_bezier(point, piece)
            least = distance if least is None else min(least, distance)
    return least


def from_polyline(point, polyline):
    """The distance from a point to a polyline: the least over its segments, a bulge's arc built
    from its ends and its bulge, the tangent of a quarter of its sweep."""
    closed, vertices = polyline
    if len(vertices) == 1:
        return away(point, vertices[0][:2])
    least = None
    for index in range(len(vertices) if closed else len(vertices) - 1):
        (x0, y0, bulge), (x1, y1, _) = vertices[index], vertices[(index + 1) % len(vertices)]
        if bulge == 0 or (x0, y0) == (x1, y1):
            distance = from_run(point, (x0, y0), (x1 - x0, y1 - y0), (mpf(0), mpf(1)))
        else:
            offset = (1 - bulge * bulge) / (4 * bulge)
            cx = (x0 + x1) / 2 - (y1 - y0) * offset
            cy = (y0 + y1) / 2 + (x1 - x0) * offset
            radius = away((x0, y0), (cx, cy))
            start = mp.atan2(y0 - cy, x0 - cx)
            distance = from_arc(point, (cx, cy, radius, start, 4 * mp.atan(bulge)))
        least = distance if least is None else min(least, distance)
    return least


def distance(point, curve):
    """The distance from a point to a curve, never less than the true one."""
    kind, values = curve
    if kind == "line":
        way = (values[2] - values[0], values[3] - values[1])
        return from_run(point, values[:2], way, (mpf(0), mpf(1)))
    if kind in ("ray", "xline"):
        stretch = (mpf(0) if kind == "ray" else -mp.inf, mp.inf)
        return from_run(point, values[:2], values[2:4], stretch)
    if kind == "arc":
        return from_arc(point, values)
    if kind == "ellipse":
        return from_ellipse(point, values)
    if kind in ("quadratic", "cubic"):
        poles = [(values[2 * i], values[2 * i + 1], mpf(1)) for i in range(len(values) // 2)]
        return from_bezier(point, poles)
    if kind == "spline":
        return from_spline(point, values)
    return from_polyline(point, values)


def sized(curve):
    """The largest number a curve is made of, for the size its doubles' steps are taken at."""
    kind, values = curve
    if kind == "polyline":
        numbers = [abs(number) for vertex in values[1] for number in vertex[:2]]
    elif kind == "spline":
        numbers = [abs(number) for pole in values[2] for number in pole]
    elif kind == "arc":
        numbers = [abs(values[0]) + abs(values[2]), abs(values[1]) + abs(values[2])]
    elif kind == "ellipse":
        numbers = [abs(number) for number in values[:6]]
    else:
        numbers = [abs(number) for number in values]
    return max(numbers + [mpf(1)])


def cases(path):
    """Each case of a file: its name, tolerance, two curves and the points it expects, each with
    the slack it is given."""
    name, tolerance, curves, points = None, mpf(0), {}, []
    with open(path, encoding="utf-8") as lines:
        for line in lines:
            words = line.split()
            if not words or words[0].startswith("#"):
                continue
            if words[0] == "case":
                name, tolerance, curves, points = line[5:].strip(), mpf(0), {}, []
            elif words[0] == "tolerance":
                tolerance = mpf(float(words[1]))
            elif words[0] in ("first", "second"):
                curves[words[0]] = parsed(line.split(None, 1)[1])
            elif words[0] == "hit":
                slack = mpf(float(words[4])) if len(words) > 4 else mpf(0)
                points.append(((mpf(float(words[1])), mpf(float(words[2]))), slack, words[3]))
            elif words[0] == "shared":
                for index in (1, 3):
                    end = (mpf(float(words[index])), mpf(float(words[index + 1])))
                    if mp.isfinite(end[0]) and mp.isfinite(end[1]):
                        points.append((end, mpf(0), "shared"))
            elif words[0] == "end":
                yield name, tolerance, curves, points


def main():
    """Checks every case file named, printing each point that is off and a count of the rest."""
    wrong = 0
    for path in sys.argv[1:]:
        tally = {"exact": 0, "near": 0, "off": 0}
        worst = mpf(0)
        for name, tolerance, curves, points in cases(path):
            size = max(sized(curves["first"]), sized(curves["second"]))
            step = size * mpf(2) ** -52
            for point, slack, kind in points:
                farthest = max(distance(point, curves["first"]),
                               distance(point, curves["second"]))
                if farthest <= 2 * step:
                    tally["exact"] += 1
                    worst = max(worst, farthest / step)
                elif farthest <= max(tolerance, slack):
                    tally["near"] += 1
                else:
                    tally["off"] += 1
                    where = f"{mp.nstr(point[0], 17)} {mp.nstr(point[1], 17)}"
                    print(f"OFF  {name}: {kind} {where} is {mp.nstr(farthest, 5)} from a curve; "
                          f"tolerance {mp.nstr(tolerance, 5)}, slack {mp.nstr(slack, 5)}")
        wrong += tally["off"]
        print(f"{path}: {tally['exact']} exact (the farthest {mp.nstr(worst, 3)} of a step), "
              f"{tally['near']} near, {tally['off']} off")
    return 1 if wrong else 0


if __name__ == "__main__":
    sys.exit(main())
