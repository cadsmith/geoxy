#!/usr/bin/env python3
"""Checks the exact oracle on cases worked out by hand, before it is trusted with geoxy's suite.

    python3 -I tools/suites/check.py

Each case states its answer from algebra done on paper: the roots of a polynomial built from
them, where a line crosses a circle at a Pythagorean point, where a parabola meets a line. The
oracle must give exactly those points, said tangent exactly where they are, and the stretches
two curves share. The exit status is the number of cases that fail.
"""

import pathlib
import sys
from fractions import Fraction as F

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import oracle  # noqa: E402  pylint: disable=wrong-import-position

EXACT = F(1, 2 ** 100)


def same(found, wanted):
    """Whether two points agree to within 2^-100."""
    return abs(found[0] - wanted[0]) <= EXACT and abs(found[1] - wanted[1]) <= EXACT


def segment(start, end):
    """The oracle's piece for a segment between two points."""
    return oracle.line(start, oracle.minus(end, start), F(0), F(1))


def circle(center, radius):
    """The oracle's piece for a whole circle."""
    return oracle.conic(center, (radius, F(0)), (F(0), radius), (F(1), F(0)), (F(1), F(0)),
                        True, True)


def check_hits(name, found, wanted):
    """Whether the hits found are the ones wanted, each a point and whether it is a tangency,
    in any order; prints what differs."""
    left = list(found.hits)
    for point, tangent in wanted:
        match = next((hit for hit in left if same(hit[0], point) and hit[1] == tangent), None)
        if match is None:
            print(f"{name}: no hit at {point}, tangent {tangent}; found {found.hits}")
            return False
        left.remove(match)
    if left:
        print(f"{name}: hits not expected {left}")
        return False
    return True


def check_shared(name, found, wanted):
    """Whether the stretches found shared are the ones wanted, each by its two ends either
    way round; prints what differs."""
    if len(found.shared) != len(wanted):
        print(f"{name}: shared {found.shared}, wanted {wanted}")
        return False
    for (a, b), (p, q) in zip(found.shared, wanted):
        if not ((same(a, p) and same(b, q)) or (same(a, q) and same(b, p))):
            print(f"{name}: shared {(a, b)}, wanted {(p, q)}")
            return False
    return True


def narrowed(found, wanted):
    """Whether roots found are the ones wanted, each to within the 2^-140 it is narrowed to, and
    each said repeated or not as wanted."""
    return len(found) == len(wanted) and all(
        abs(root - want) <= F(1, 2 ** 140) and repeated == twice
        for (root, repeated), (want, twice) in zip(found, wanted))


def roots():
    """Polynomials whose roots are known: one simple and one double, and the square root of 2."""
    failed = 0
    # (x - 1)(x - 2)^2 = x^3 - 5x^2 + 8x - 4
    found = oracle.roots([F(-4), F(8), F(-5), F(1)], F(0), F(3))
    if not narrowed(found, [(F(1), False), (F(2), True)]):
        print(f"roots of (x - 1)(x - 2)^2: {found}")
        failed += 1
    found = oracle.roots([F(-2), F(0), F(1)], F(0), F(2))
    if len(found) != 1 or abs(found[0][0] ** 2 - 2) > F(1, 2 ** 130) or found[0][1]:
        print(f"roots of x^2 - 2: {found}")
        failed += 1
    # A root on the stretch's start, which bisection must keep and not drift from.
    found = oracle.roots([F(0), F(-1), F(1)], F(0), F(2))
    if not narrowed(found, [(F(0), False), (F(1), False)]):
        print(f"roots of x^2 - x: {found}")
        failed += 1
    return failed


def lines():
    """Segments that cross, and segments on one line that share a stretch or only an end."""
    half = F(1, 2)
    cases = [
        ("segments cross",
         oracle.meet_pieces(segment((F(0), F(0)), (F(2), F(2))),
                            segment((F(0), F(2)), (F(2), F(0)))),
         [((F(1), F(1)), False)], []),
        ("segments share",
         oracle.meet_pieces(segment((F(0), F(0)), (F(2), F(0))),
                            segment((F(1), F(0)), (F(3), F(0)))),
         [], [((F(1), F(0)), (F(2), F(0)))]),
        ("segments end to end",
         oracle.meet_pieces(segment((F(0), F(0)), (F(1), F(0))),
                            segment((F(1), F(0)), (F(3), F(0)))),
         [((F(1), F(0)), False)], []),
        ("segments apart",
         oracle.meet_pieces(segment((F(0), F(0)), (F(1), F(1))),
                            segment((F(0), half), (F(0), F(2)))),
         [], []),
    ]
    return sum(not (check_hits(name, found, hits) and check_shared(name, found, shared))
               for name, found, hits, shared in cases)


def conics():
    """A circle a line crosses at Pythagorean points, one it touches, an ellipse, and two arcs
    of one circle that share a stretch."""
    one, zero = F(1), F(0)
    three, four = F(3, 5), F(4, 5)
    ellipse = oracle.conic((zero, zero), (F(2), zero), (zero, one), (one, zero), (one, zero), True,
                           True)
    quarter = oracle.conic((zero, zero), (one, zero), (zero, one), (one, zero), (zero, one), True,
                           False)
    later = oracle.conic((zero, zero), (one, zero), (zero, one), (three, four), (-one, zero), True,
                         False)
    cases = [
        ("line crosses circle",
         oracle.meet_pieces(circle((zero, zero), one), segment((-one, three), (one, three))),
         [((-four, three), False), ((four, three), False)], []),
        ("line touches circle",
         oracle.meet_pieces(circle((zero, zero), one), segment((-one, one), (one, one))),
         [((zero, one), True)], []),
        ("line crosses ellipse",
         oracle.meet_pieces(ellipse, segment((-F(2), three), (F(2), three))),
         [((-F(8, 5), three), False), ((F(8, 5), three), False)], []),
        ("arcs share",
         oracle.meet_pieces(quarter, later),
         [], [((three, four), (zero, one))]),
        ("circles cross",
         oracle.meet_pieces(circle((zero, zero), one), circle((F(6, 5), zero), one)),
         [((F(3, 5), four), False), ((F(3, 5), -four), False)], []),
    ]
    return sum(not (check_hits(name, found, hits) and check_shared(name, found, shared))
               for name, found, hits, shared in cases)


def beziers():
    """The parabola y = x^2 as a quadratic, met by a line through two of its points, by its
    tangent at the vertex, and by the parabola turned over, which it touches there."""
    one, zero = F(1), F(0)
    parabola = oracle.bezier([(-one, one), (zero, -one), (one, one)])
    cases = [
        ("line crosses parabola",
         oracle.meet_pieces(parabola, oracle.bezier([(-one, -one), (F(2), F(2))])),
         [((zero, zero), False), ((one, one), False)], []),
        ("tangent touches parabola",
         oracle.meet_pieces(parabola, oracle.bezier([(-one, zero), (one, zero)])),
         [((zero, zero), True)], []),
        ("parabolas touch",
         oracle.meet_pieces(parabola, oracle.bezier([(-one, -one), (zero, one), (one, -one)])),
         [((zero, zero), True)], []),
        ("parabola and its half share",
         oracle.meet_pieces(parabola, oracle.bezier([(zero, zero), (F(1, 2), zero), (one, one)])),
         [], [((zero, zero), (one, one))]),
    ]
    return sum(not (check_hits(name, found, hits) and check_shared(name, found, shared))
               for name, found, hits, shared in cases)


def angles_and_places():
    """The angle a crossing is at, and a loop's double point met twice: once for each of the
    loop's places through it, whichever piece the oracle solves in."""
    failed = 0
    square = oracle.meet_pieces(segment((F(0), F(0)), (F(2), F(2))),
                                segment((F(0), F(2)), (F(2), F(0))))
    if square.sines != [1.0]:
        print(f"segments at a right angle: sines {square.sines}")
        failed += 1
    tangent = oracle.meet_pieces(circle((F(0), F(0)), F(1)), segment((-F(1), F(1)), (F(1), F(1))))
    if tangent.sines != [0.0]:
        print(f"line touching a circle: sines {tangent.sines}")
        failed += 1
    # x = s^2, y = s^3 - s for s in [-3/2, 3/2]: through (1, 0) at s = -1 and s = 1.
    loop = oracle.bezier([(F(9, 4), -F(15, 8)), (-F(3, 4), F(31, 8)), (-F(3, 4), -F(31, 8)),
                          (F(9, 4), F(15, 8))])
    upright = segment((F(1), -F(2)), (F(1), F(2)))
    straight = oracle.bezier([(F(1), -F(2)), (F(1), F(2))])
    for name, found in (("line first", oracle.meet_pieces(upright, loop)),
                        ("loop first", oracle.meet_pieces(loop, straight))):
        at = [p for p, _ in found.hits if same(p, (F(1), F(0)))]
        tangents = [t for _, t in found.hits]
        if len(at) != 2 or any(tangents):
            print(f"loop's double point, {name}: hits {found.hits}")
            failed += 1
    return failed


def rationals():
    """A rational quadratic, weights 1, 1/2, 1: y = x crosses it at its middle, (2/3, 2/3); the
    plain quadratic on its points meets it only at their ends, tangent to it at both."""
    one, zero = F(1), F(0)
    points = [(one, zero), (one, one), (zero, one)]
    conic = oracle.bezier(points, [one, F(1, 2), one])
    cases = [
        ("line crosses a rational quadratic",
         oracle.meet_pieces(conic, oracle.bezier([(zero, zero), (one, one)])),
         [((F(2, 3), F(2, 3)), False)], []),
        ("quadratics on one control polygon",
         oracle.meet_pieces(conic, oracle.bezier(points)),
         [((one, zero), True), ((zero, one), True)], []),
    ]
    return sum(not (check_hits(name, found, hits) and check_shared(name, found, shared))
               for name, found, hits, shared in cases)


def main():
    """Runs every case; the exit status is how many failed."""
    failed = roots() + lines() + conics() + beziers() + angles_and_places() + rationals()
    print(f"oracle: {failed} failed" if failed else "oracle: every case agrees")
    return failed


if __name__ == "__main__":
    sys.exit(main())
