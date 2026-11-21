#!/usr/bin/env python3
"""Checks the whole-number code integers.v proves, as the compiler builds it, against exact
rational arithmetic.

    python3 -I tests/proof/integers.py [compiler ...]

integers.v proves, step for step, a model of source/exact/grid.cpp and whole.cpp. This runs the
code itself: it builds the two files with a probe, with each compiler named (Clang and Homebrew's
GCC by default), once with the compiler's 128 bit products and once with MultiplyHalves in their
place, and asks OrientGrid, OnGrid, OrientWhole and CircleWhole
the signs of many determinants: random doubles across the whole range, subnormals and the largest
included; places on one grid, whose span lies at the grids' limits; and points exactly collinear
or cocircular, scaled by powers of two, where the sign is zero. Python's fractions give every
determinant exactly. Each sign must be the exact one, and each grid must answer exactly where its
span check lets it.

The exit status is 1 when a check fails.
"""

import fractions
import pathlib
import random
import shutil
import struct
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]

PROBE = r"""#include "exact/grid.hpp"
#include "exact/whole.hpp"

#include <bit>
#include <cstdint>
#include <cstdio>
#include <optional>

namespace
{
  auto Read(std::uint64_t bits) -> double { return std::bit_cast<double>(bits); }

  void Say(std::optional<int> sign)
  {
    if (sign) { std::printf("%d ", *sign); } else { std::printf("n "); }
  }
}

int main()
{
  char kind = 0;
  while (std::scanf(" %c", &kind) == 1)
  {
    unsigned long long raw[8] = {};
    const int count = kind == 'o' ? 6 : 8;
    for (int index = 0; index < count; ++index)
    {
      if (std::scanf("%llx", &raw[index]) != 1) { return 1; }
    }
    const Geoxy::Exact::Triangle corners = {.First  = {Read(raw[0]), Read(raw[1])},
                                           .Second = {Read(raw[2]), Read(raw[3])},
                                           .Third  = {Read(raw[4]), Read(raw[5])}};
    if (kind == 'o')
    {
      Say(Geoxy::Exact::OrientGrid(corners));
      std::printf("%d\n", Geoxy::Exact::OrientWhole(corners));
    }
    else
    {
      const Geoxy::Point point = {Read(raw[6]), Read(raw[7])};
      Say(Geoxy::Exact::OnGrid(corners, point));
      std::printf("%d\n", Geoxy::Exact::CircleWhole(corners, point));
    }
  }
  return 0;
}
"""


def bits(value):
    return struct.unpack(">Q", struct.pack(">d", value))[0]


def exact(value):
    return fractions.Fraction(value)


def sign(value):
    return (value > 0) - (value < 0)


def orient(points):
    ax, ay, bx, by, cx, cy = map(exact, points)
    return sign((ax - cx) * (by - cy) - (ay - cy) * (bx - cx))


def circle(points):
    ax, ay, bx, by, cx, cy, dx, dy = map(exact, points)
    ax, ay, bx, by, cx, cy = ax - dx, ay - dy, bx - dx, by - dy, cx - dx, cy - dy
    return sign((ax * ax + ay * ay) * (bx * cy - by * cx)
                + (bx * bx + by * by) * (cx * ay - cy * ax)
                + (cx * cx + cy * cy) * (ax * by - ay * bx))


def span(points):
    """The grids' reading: None where a coordinate is subnormal, else the highest place minus the
    lowest set one, or None too when all are zero."""
    lowest, highest = None, None
    for value in points:
        if value == 0:
            continue
        mantissa, exponent = (bits(value) & ((1 << 52) - 1)), (bits(value) >> 52) & 0x7FF
        if exponent == 0:
            return "subnormal"
        significand = mantissa | (1 << 52)
        place = exponent - 1075
        low = place + (significand & -significand).bit_length() - 1
        lowest = low if lowest is None else min(lowest, low)
        highest = place + 52 if highest is None else max(highest, place + 52)
    return None if lowest is None else highest - lowest


def any_double(draw):
    """A finite double of any bits."""
    while True:
        value = struct.unpack(">d", struct.pack(">Q", draw.getrandbits(64)))[0]
        if value == value and abs(value) != float("inf"):
            return value


def banded(draw, exponent, width):
    """A double near 2^exponent whose set bits lie within width places of it."""
    whole = draw.getrandbits(width) | 1
    return draw.choice((1, -1)) * whole * 2.0 ** (exponent - width)


def cases(draw):
    """(kind, points): random, banded at the grids' limits, and degenerate."""
    edges = [0.0, -0.0, 5e-324, -5e-324, 2.2250738585072014e-308, 1.7976931348623157e308,
             -1.7976931348623157e308, 1.0, -1.0, 0.5, 3.0]
    for _ in range(3000):
        yield "o", [any_double(draw) for _ in range(6)]
        yield "c", [any_double(draw) for _ in range(8)]
        yield "o", [draw.choice(edges) for _ in range(6)]
        yield "c", [draw.choice(edges) for _ in range(8)]
    for width in (40, 52, 53, 60, 61, 62, 70, 72, 73, 100, 124, 125, 126):
        for _ in range(300):
            exponent = draw.randint(-1000, 1000)
            yield "o", [banded(draw, exponent - draw.randint(0, max(0, width - 53)),
                               min(53, width)) for _ in range(6)]
            yield "c", [banded(draw, exponent - draw.randint(0, max(0, width - 53)),
                               min(53, width)) for _ in range(8)]
    for _ in range(3000):
        scale = 2.0 ** draw.randint(-1060, 960)
        step = draw.choice((scale, scale * 2.0 ** draw.randint(-60, 60)))
        x, y = draw.randint(-2 ** 20, 2 ** 20), draw.randint(-2 ** 20, 2 ** 20)
        dx, dy = draw.randint(-2 ** 10, 2 ** 10), draw.randint(-2 ** 10, 2 ** 10)
        s, t = draw.randint(-2 ** 10, 2 ** 10), draw.randint(-2 ** 10, 2 ** 10)
        points = [x * step, y * step, (x + s * dx) * step, (y + s * dy) * step,
                  (x + t * dx) * step, (y + t * dy) * step]
        yield "o", points
        a, b = draw.choice(((3, 4), (5, 12), (8, 15), (20, 21), (119, 120)))
        images = [(a, b), (-a, b), (a, -b), (-a, -b), (b, a), (-b, a), (b, -a), (-b, -a)]
        cx, cy = draw.randint(-2 ** 30, 2 ** 30), draw.randint(-2 ** 30, 2 ** 30)
        chosen = [draw.choice(images) for _ in range(4)]
        yield "c", [coordinate for u, v in chosen for coordinate in ((cx + u) * scale,
                                                                     (cy + v) * scale)]


def compilers():
    named = sys.argv[1:]
    if named:
        return named
    found = [shutil.which("clang++")]
    found += sorted(str(path) for path in pathlib.Path("/opt/homebrew/bin").glob("g++-[0-9]*"))[-1:]
    return [compiler for compiler in found if compiler]


def run(compiler, folder, flags):
    target = pathlib.Path(folder) / f"probe{len(flags)}-{pathlib.Path(compiler).name}"
    source = pathlib.Path(folder) / "probe.cpp"
    source.write_text(PROBE)
    command = [compiler, "-std=c++23", "-O2", "-ffp-contract=off", "-DGEOXY_STATIC", *flags,
               "-I", str(ROOT / "source"), "-I", str(ROOT / "include"), "-o", str(target),
               str(source), str(ROOT / "source/exact/grid.cpp"),
               str(ROOT / "source/exact/whole.cpp")]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return target


def main():
    failed = False

    def report(name, held, told=""):
        nonlocal failed
        failed = failed or not held
        print(f"{name:58} {'held' if held else 'BROKEN':7} {told}", flush=True)

    all_cases = list(cases(random.Random(20261113)))
    lines = "".join(f"{kind} {' '.join(f'{bits(v):x}' for v in points)}\n"
                    for kind, points in all_cases)
    with tempfile.TemporaryDirectory() as folder:
        builds = (("128 bit products", []), ("MultiplyHalves", ["-U__SIZEOF_INT128__"]))
        for compiler, (kind, flags) in ((compiler, build) for compiler in compilers()
                                         for build in builds):
            name = f"{pathlib.Path(compiler).name} {kind}"
            probe = run(compiler, folder, flags)
            out = subprocess.run([str(probe)], input=lines, capture_output=True, text=True,
                                 check=True).stdout.split("\n")
            wrong = {"grid": [], "whole": [], "grid refused": []}
            answered = {"o": 0, "c": 0}
            for (kind, points), line in zip(all_cases, out):
                grid, whole = line.split()
                expected = orient(points) if kind == "o" else circle(points)
                if int(whole) != expected:
                    wrong["whole"].append(points)
                limit = 125 if kind == "o" else 61
                reach = span(points)
                fits = reach != "subnormal" and (reach is None or reach < limit)
                if grid == "n":
                    if fits:
                        wrong["grid refused"].append(points)
                else:
                    answered[kind] += 1
                    if int(grid) != expected or not fits:
                        wrong["grid"].append(points)
            report(f"{name}: OrientWhole and CircleWhole exact, {len(all_cases)} cases",
                   not wrong["whole"], str(wrong["whole"][:1]))
            report(f"{name}: OrientGrid, OnGrid exact, {answered['o']}+{answered['c']} answered",
                   not wrong["grid"], str(wrong["grid"][:1]))
            report(f"{name}: the grids answer wherever their span fits",
                   not wrong["grid refused"], str(wrong["grid refused"][:1]))
            zeros = sum(1 for (kind, points) in all_cases
                        if (orient(points) if kind == "o" else circle(points)) == 0)
            report(f"{name}: the cases hold exact zeros", zeros > 1000, str(zeros))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
