/***************************************************************************************************
 *  @file      tests/geoxy/arc.cpp
 *  @brief     Tests of arcs and circles.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/arc.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/point.hpp"

#include <doctest/doctest.h>

namespace Geoxy
{
  TEST_CASE("arc: a circle is the arc of a whole turn from the x axis")
  {
    const Arc circle = Circle({.X = 5e6, .Y = 1}, 2);
    CHECK(circle.Center == Point{.X = 5e6, .Y = 1});
    CHECK(circle.Radius == 2.0);
    CHECK(circle.Start  == Radians(0));
    CHECK(circle.Sweep  == Turn);
  }

  TEST_CASE("arc: equal by center, radius, start and sweep")
  {
    const Arc arc = {
        .Center = {.X = 1, .Y = 2}, .Radius = 3, .Start = Degrees(30), .Sweep = Degrees(-90)};
    Arc other = arc;
    CHECK(arc == other);
    other.Center.X = 0;
    CHECK_FALSE(arc == other);
    other        = arc;
    other.Radius = 4;
    CHECK_FALSE(arc == other);
    other       = arc;
    other.Start = Degrees(31);
    CHECK_FALSE(arc == other);
    other       = arc;
    other.Sweep = Degrees(90);
    CHECK_FALSE(arc == other);
  }
} // namespace Geoxy
