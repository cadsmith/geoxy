/***************************************************************************************************
 *  @file      tests/geoxy/axis.cpp
 *  @brief     Tests that two lines' or circles' carriers are told to cross on an axis exactly where
 *             they do, and nowhere a step of a double off it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "axis.hpp"

#include "abut.hpp"
#include "real.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/tolerance.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

namespace Geoxy
{
  namespace
  {
    /**
     *  A circle about a center.
     *
     *  @param  x      The center's x.
     *  @param  y      The center's y.
     *  @param  radius The radius.
     *  @return        The circle, as a whole arc.
     */
    auto Round(double x, double y, double radius) -> Arc
    {
      return {.Center = {.X = x, .Y = y},
              .Radius = radius,
              .Start  = Radians(0.0),
              .Sweep  = Radians(2.0 * Pi)};
    }

    /**
     *  Whether two curves made ready meet as the curves themselves, either first, to the bit: each
     *  mark's point, a zero's sign too, and kind.
     *
     *  @param  first     One curve.
     *  @param  second    The other.
     *  @param  tolerance The tolerance.
     *  @return           Whether they do.
     */
    auto Bits(const Curve &first, const Curve &second, const Tolerance &tolerance) -> bool
    {
      const auto once = [&tolerance](const Curve &one, const Curve &other) -> bool
      {
        const Prepared ahead{one};
        const Prepared behind{other};
        Intersections  meeting;
        Marking        marking;
        Intersect({one, other}, tolerance, meeting);
        Intersect({ahead, behind}, tolerance, marking);
        const auto same = [](const Hit &hit, const Mark &mark) -> bool
        {
          return std::bit_cast<std::uint64_t>(hit.At.X) == std::bit_cast<std::uint64_t>(mark.At.X)
              && std::bit_cast<std::uint64_t>(hit.At.Y) == std::bit_cast<std::uint64_t>(mark.At.Y)
              && hit.Kind                               == mark.Kind;
        };
        return std::ranges::equal(meeting.Hits, marking.Marks, same)
            && meeting.Overlaps.size() == marking.Overlaps.size();
      };
      return once(first, second) && once(second, first);
    }

    /**
     *  Pythagorean triples: a point on a circle about the origin of the last's radius.
     */
    constexpr std::array<std::array<double, 3>, 5> Triples = {
        {{3, 4, 5}, {5, 12, 13}, {8, 15, 17}, {7, 24, 25}, {20, 21, 29}}};
  } // namespace

  TEST_CASE(
      "axis: lines and arcs crossing on an axis, made ready, meet as Intersect answers, to the bit")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    for (const double scale : {1.0, 0x1p-20, 0x1p20, 0x1p93})
    {
      for (const std::array<double, 3> &triple : Triples)
      {
        const double radius = triple[2] * scale;
        const Point  point  = {.X = triple[0] * scale, .Y = triple[1] * scale};
        const Arc    circle = Round(0, 0, radius);
        const Arc    half   = {
            .Center = {}, .Radius = radius, .Start = Degrees(-90), .Sweep = Degrees(180)};
        const Tolerance near = {.Length = tolerance.Length * scale, .Radians = 1e-9};
        for (const Point axis : {Point{.X = 0, .Y = -radius}, Point{.X = radius, .Y = 0},
                                 Point{.X = -radius, .Y = 0}, Point{.X = 0, .Y = radius}})
        {
          const Point beyond = point + ((point - axis) * 0.5);
          const Point before = axis + ((axis - point) * 0.5);
          CHECK(Bits(Line{.Start = before, .End = beyond}, circle, near));
          CHECK(Bits(Line{.Start = before, .End = beyond}, half, near));
          /**
           *  A second circle through the point on the axis and the other, its center the first's
           *  mirrored across their chord.
           */
          const Arc other = Round(point.X + axis.X, point.Y + axis.Y, radius);
          CHECK(Bits(circle, other, near));
          CHECK(Bits(half, other, near));
          /**
           *  A level or upright line on the axis's crossing, at +0 and -0.
           */
          for (const double zero : {0.0, -0.0})
          {
            const bool  level = Zero(axis.Y);
            const Point from =
                level ? Point{.X = zero, .Y = -2 * radius} : Point{.X = -2 * radius, .Y = zero};
            const Point until =
                level ? Point{.X = zero, .Y = 2 * radius} : Point{.X = 2 * radius, .Y = zero};
            CHECK(Bits(Line{.Start = from, .End = until}, circle, near));
            CHECK(Bits(Line{.Start = from, .End = until}, half, near));
          }
        }
      }
    }
  }

  TEST_CASE("axis: crossings on an axis and a hair off it, made ready, meet as Intersect answers")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    /**
     *  A circle a hair off the origin, too fine to be asked of: its crossing within a bound of the
     *  axis is not told on it.
     */
    const Arc  hair  = Round(1e-35, 0, 5);
    const Line slant = {.Start = {.X = -1, .Y = -6}, .End = {.X = 6, .Y = 1}};
    CHECK(Bits(slant, hair, tolerance));
    CHECK(Bits(hair, Round(-5, 5, 5), tolerance));
    /**
     *  A level line a hair above the axis, its crossing 1.25e-41 off it; and level and upright
     *  lines on the axes through a circle's crossings with them at the origin.
     */
    const Arc offset = Round(4, 3, 5);
    CHECK(Bits(Line{.Start = {.X = -3, .Y = 1e-20}, .End = {.X = 9, .Y = 1e-20}}, Round(4, 0, 4),
               tolerance));
    CHECK(Bits(Line{.Start = {.X = -3, .Y = 0}, .End = {.X = 9, .Y = 0}}, offset, tolerance));
    CHECK(Bits(Line{.Start = {.X = 0, .Y = -3}, .End = {.X = 0, .Y = 9}}, offset, tolerance));
    /**
     *  Two circles crossing shallowly at (0, 1), as the grazing pass tells them.
     */
    CHECK(Bits(Round(0, -144, 145), Round(17, -143, 145), tolerance));
    /**
     *  The pass for ends, on a wall at +0 and one a hair off it, ending on an arc.
     */
    const Arc swing = {
        .Center = {.X = 4, .Y = 3}, .Radius = 5, .Start = Degrees(180), .Sweep = Degrees(60)};
    const Line wall  = {.Start = {.X = 0, .Y = 0}, .End = {.X = 9, .Y = 0}};
    const Line above = {.Start = {.X = -2, .Y = 1e-20}, .End = {.X = 9, .Y = 1e-20}};
    const Arc  lower = {
        .Center = {.X = 4, .Y = 0}, .Radius = 4, .Start = Degrees(90), .Sweep = Degrees(180)};
    CHECK(Bits(wall, swing, tolerance));
    CHECK(Bits(above, lower, tolerance));
    const Curve    wallCurve  = wall;
    const Curve    swingCurve = swing;
    const Curve    aboveCurve = above;
    const Curve    lowerCurve = lower;
    const Prepared wallReady{wallCurve};
    const Prepared swingReady{swingCurve};
    const Prepared aboveReady{aboveCurve};
    const Prepared lowerReady{lowerCurve};
    CHECK(Abutting::Rimmed({wallReady, swingReady}, 1e-9).has_value());
    CHECK(!Abutting::Rimmed({aboveReady, lowerReady}, 1e-9).has_value());
  }

  TEST_CASE("axis: numbers held, and not")
  {
    CHECK(Axis::Held(0.0));
    CHECK(Axis::Held(-0.0));
    CHECK(Axis::Held(1.0));
    CHECK(Axis::Held(-0x1p100));
    CHECK(Axis::Held(0x1p-100));
    CHECK(!Axis::Held(0x1p101));
    CHECK(!Axis::Held(-0x1p-101));
    CHECK(!Axis::Held(std::numeric_limits<double>::infinity()));
    CHECK(!Axis::Held(std::numeric_limits<double>::quiet_NaN()));
  }

  TEST_CASE("axis: a point on a line's carrier exactly, and a step of a double off it")
  {
    const Line slant = {.Start = {.X = -5, .Y = -5}, .End = {.X = 5, .Y = 5}};
    CHECK(Axis::On(slant, {.X = 0, .Y = 0}));
    CHECK(Axis::On(slant, {.X = -0.0, .Y = 7}) == false);
    CHECK(Axis::On(slant, {.X = 0x1p-60, .Y = 0x1p-60}));
    CHECK(!Axis::On(slant, {.X = 0, .Y = std::nextafter(0.0, 1.0)}));
    const Line steep = {.Start = {.X = 0.1, .Y = 0.3}, .End = {.X = 0.35, .Y = 2.2}};
    CHECK(Axis::On(steep, steep.Start));
    CHECK(Axis::On(steep, steep.End));
    CHECK(!Axis::On(steep, {.X = 0, .Y = 0}));
    const Line level = {.Start = {.X = 1, .Y = 0}, .End = {.X = 6, .Y = 0}};
    CHECK(Axis::On(level, {.X = 2.5, .Y = -0.0}));
    CHECK(!Axis::On(level, {.X = 2.5, .Y = 0x1p-90}));
    /**
     *  A line of no length, or a number not held: never told.
     */
    CHECK(!Axis::On(Line{.Start = {.X = 1, .Y = 1}, .End = {.X = 1, .Y = 1}}, {.X = 1, .Y = 1}));
    CHECK(!Axis::On(Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 0x1p120, .Y = 0x1p120}},
                    {.X = 1, .Y = 1}));
    CHECK(!Axis::On(slant, {.X = 0x1p-120, .Y = 0x1p-120}));
  }

  TEST_CASE("axis: a point on an arc's circle exactly, and a step of a double off it")
  {
    const Arc circle = Round(0, 0, 5);
    CHECK(Axis::On(circle, {.X = 0, .Y = -5}));
    CHECK(Axis::On(circle, {.X = 3, .Y = 4}));
    CHECK(Axis::On(circle, {.X = -4, .Y = -0.0}) == false);
    CHECK(!Axis::On(circle, {.X = 0, .Y = std::nextafter(5.0, 6.0)}));
    const Arc moved = Round(-5, 5, 5);
    CHECK(Axis::On(moved, {.X = 0, .Y = 5}));
    CHECK(Axis::On(moved, {.X = -5, .Y = 0}));
    const Arc fine = Round(0.125, 0.375, 0.5);
    CHECK(Axis::On(fine, {.X = 0.125, .Y = 0.875}));
    CHECK(!Axis::On(Round(0.1, 0.3, 0.5), {.X = 0.1, .Y = 0.8}));
    CHECK(!Axis::On(fine, {.X = 0, .Y = 0.875}));
    /**
     *  A number not held: never told.
     */
    CHECK(!Axis::On(Round(0, 0, 0x1p-120), {.X = 0, .Y = 0x1p-120}));
    CHECK(!Axis::On(circle, {.X = 0x1p101, .Y = 0}));
  }

  TEST_CASE("axis: numerators that Intersect's exact way takes to zero")
  {
    CHECK(Axis::Fits(1.0));
    CHECK(Axis::Fits(-0x1p279));
    CHECK(!Axis::Fits(0.0));
    CHECK(!Axis::Fits(0x1p-281));
    CHECK(!Axis::Fits(0x1p281));
  }
} // namespace Geoxy
