#!/usr/bin/env python3
"""Scores fillets and chamfers against the truth worked out at 120 significant digits.

    uv run --offline --no-project --with mpmath python -I tools/reference/fillet.py CASES

Each line of the case file is a family's name, two curves' carriers, then what was asked and each
library's answer, fields parted by "|", numbers as hexadecimal doubles:

    lines L sx sy ex ey | C cx cy r | r radius | geoxy ox oy ax ay bx by | occt ...
    chamfer L sx sy ex ey | L sx sy ex ey | d d1 d2 | geoxy ax ay bx by | occt ...

L is the line through two points, C a circle. For a fillet, every circle of the radius tangent to
both carriers is worked out in closed form, a line's offsets by the radius and a circle grown or
shrunk by it, crossed; the one whose center is nearest a library's is its truth, with the foot of
its center on each line and the point C + (O - C) R / (R +- r) on each circle. For a chamfer of two
lines, the truth is their corner and the points each distance from it toward the line's first
point, the end away from the corner. A library that wrote "fail" is counted, not scored.

Each coordinate's error is |value - truth| in steps of a double at the truth: 0 for the double
nearest, never more than 0.5 for it. Printed per family and library: how many answered, the largest
and the median error, and how many coordinates were the nearest double.

mpmath is needed: uv run --offline --no-project --with mpmath, as above.
"""

import math
import sys

from mpmath import mp, mpf, sqrt

mp.dps = 120


def number(text):
    """A double from its hexadecimal text, as an mpf."""
    return mpf(float.fromhex(text))


def carrier(text):
    """A carrier from its text: ("L", start, end) or ("C", center, radius)."""
    words = text.split()
    values = [number(word) for word in words[1:]]
    if words[0] == "L":
        return ("L", (values[0], values[1]), (values[2], values[3]))
    return ("C", (values[0], values[1]), values[2])


def line_offsets(line, radius):
    """A line's two offsets by a radius, as (unit normal, constant): n . x = k."""
    _, start, end = line
    way = (end[0] - start[0], end[1] - start[1])
    size = sqrt(way[0] ** 2 + way[1] ** 2)
    normal = (-way[1] / size, way[0] / size)
    base = normal[0] * start[0] + normal[1] * start[1]
    return [(normal, base + side * radius) for side in (1, -1)]


def circle_offsets(circle, radius):
    """A circle's offsets: (center, signed reach R + r and R - r)."""
    _, center, size = circle
    return [(center, size + radius), (center, size - radius)]


def lines_meet(one, other):
    """Where n1 . x = k1 and n2 . x = k2 meet; none for parallel lines."""
    (n1, k1), (n2, k2) = one, other
    det = n1[0] * n2[1] - n1[1] * n2[0]
    if det == 0:
        return []
    return [((k1 * n2[1] - k2 * n1[1]) / det, (n1[0] * k2 - n2[0] * k1) / det)]


def line_circle(line, circle):
    """Where n . x = k meets the circle of a center and a reach."""
    (normal, base), (center, reach) = line, circle
    gap = base - (normal[0] * center[0] + normal[1] * center[1])
    rest = reach * reach - gap * gap
    if rest < 0:
        return []
    foot = (center[0] + normal[0] * gap, center[1] + normal[1] * gap)
    along = (-normal[1], normal[0])
    root = sqrt(rest)
    return [(foot[0] + along[0] * root, foot[1] + along[1] * root),
            (foot[0] - along[0] * root, foot[1] - along[1] * root)]


def circles_meet(one, other):
    """Where two circles of centers and reaches meet."""
    (c1, r1), (c2, r2) = one, other
    way = (c2[0] - c1[0], c2[1] - c1[1])
    square = way[0] ** 2 + way[1] ** 2
    if square == 0:
        return []
    middle = (square + r1 * r1 - r2 * r2) / (2 * square)
    rest = r1 * r1 / square - middle * middle
    if rest < 0:
        return []
    root = sqrt(rest)
    base = (c1[0] + way[0] * middle, c1[1] + way[1] * middle)
    return [(base[0] - way[1] * root, base[1] + way[0] * root),
            (base[0] + way[1] * root, base[1] - way[0] * root)]


def touch(curve, center, reach):
    """Where a circle about a center touches a carrier: a line's foot, a circle's point."""
    if curve[0] == "L":
        _, start, end = curve
        way = (end[0] - start[0], end[1] - start[1])
        along = ((center[0] - start[0]) * way[0] + (center[1] - start[1]) * way[1]) / (
            way[0] ** 2 + way[1] ** 2)
        return (start[0] + way[0] * along, start[1] + way[1] * along)
    _, middle, size = curve
    scale = size / reach
    return (middle[0] + (center[0] - middle[0]) * scale, middle[1] + (center[1] - middle[1]) * scale)


def circles(first, second, radius):
    """Every circle of a radius tangent to both carriers: (center, reach on each)."""
    found = []
    one = line_offsets(first, radius) if first[0] == "L" else circle_offsets(first, radius)
    two = line_offsets(second, radius) if second[0] == "L" else circle_offsets(second, radius)
    for a in one:
        for b in two:
            if first[0] == "L" and second[0] == "L":
                points = lines_meet(a, b)
            elif first[0] == "L":
                points = line_circle(a, b)
            elif second[0] == "L":
                points = line_circle(b, a)
            else:
                points = circles_meet(a, b)
            reach_a = a[1] if first[0] == "C" else None
            reach_b = b[1] if second[0] == "C" else None
            for point in points:
                found.append((point, reach_a, reach_b))
    return found


def ulps(value, truth):
    """|value - truth| in steps of a double at the truth."""
    nearest = float(truth)
    step = math.ulp(nearest) if nearest != 0 else math.ulp(1e-300)
    return float(abs(mpf(value) - truth) / step)


def scored_fillet(first, second, radius, answer):
    """The errors of a fillet's six coordinates against the truth nearest its center."""
    center = (number(answer[0]), number(answer[1]))
    best = min(circles(first, second, radius),
               key=lambda item: (item[0][0] - center[0]) ** 2 + (item[0][1] - center[1]) ** 2)
    point, reach_a, reach_b = best
    one = touch(first, point, reach_a)
    two = touch(second, point, reach_b)
    truth = [point[0], point[1], one[0], one[1], two[0], two[1]]
    return [ulps(float.fromhex(text), value) for text, value in zip(answer, truth)]


def scored_chamfer(first, second, distances, answer):
    """The errors of a chamfer's four coordinates."""
    lines = [line_offsets(curve, 0)[0] for curve in (first, second)]
    corner = lines_meet(lines[0], lines[1])[0]
    truth = []
    for curve, distance in zip((first, second), distances):
        _, start, end = curve
        far = start if (start[0] - corner[0]) ** 2 + (start[1] - corner[1]) ** 2 > (
            end[0] - corner[0]) ** 2 + (end[1] - corner[1]) ** 2 else end
        way = (far[0] - corner[0], far[1] - corner[1])
        size = sqrt(way[0] ** 2 + way[1] ** 2)
        truth += [corner[0] + way[0] / size * distance, corner[1] + way[1] / size * distance]
    return [ulps(float.fromhex(text), value) for text, value in zip(answer, truth)]


def main():
    """Reads the case file and prints the scores."""
    scores = {}
    for text in open(sys.argv[1], encoding="ascii"):
        fields = [field.strip() for field in text.split("|")]
        family, first_text = fields[0].split(" ", 1)
        first, second = carrier(first_text), carrier(fields[1])
        asked = fields[2].split()
        for field in fields[3:]:
            name, *answer = field.split()
            entry = scores.setdefault((family, name), {"answered": 0, "failed": 0, "errors": []})
            if answer == ["fail"]:
                entry["failed"] += 1
                continue
            entry["answered"] += 1
            if asked[0] == "r":
                entry["errors"] += scored_fillet(first, second, number(asked[1]), answer)
            else:
                entry["errors"] += scored_chamfer(first, second,
                                                  [number(asked[1]), number(asked[2])], answer)
    print(f"{'family':16} {'library':6} {'answered':>8} {'failed':>6} {'max ulps':>12}"
          f" {'median':>8} {'nearest':>8}")
    for (family, name), entry in scores.items():
        errors = sorted(entry["errors"])
        if not errors:
            print(f"{family:16} {name:6} {entry['answered']:8} {entry['failed']:6}")
            continue
        nearest = sum(1 for error in errors if error <= 0.5) / len(errors)
        print(f"{family:16} {name:6} {entry['answered']:8} {entry['failed']:6}"
              f" {errors[-1]:12.4g} {errors[len(errors) // 2]:8.3g} {nearest:8.2%}")


if __name__ == "__main__":
    main()
