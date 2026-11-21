/***************************************************************************************************
 *  @file      tests/geoxy/segment.cpp
 *  @brief     Tests of a polyline's segments and of bulges.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/segment.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/length.hpp"
#include "geoxy/line.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <limits>
#include <variant>

namespace Geoxy
{
  /**
   *  Measured in AutoCAD's console: from 0,0 to 10,0 a bulge of 0.5 has center 5,3.75 and radius
   *  6.25, and 1 is the half circle under the chord.
   */
  TEST_CASE("segment: a bulge's arc as AutoCAD draws it")
  {
    const Arc arc = std::get<Arc>(Bent({.X = 0, .Y = 0}, {.X = 10, .Y = 0}, 0.5));
    CHECK(arc.Center            == Point{.X = 5, .Y = 3.75});
    CHECK(arc.Radius            == 6.25);
    CHECK(arc.Sweep.InRadians() == doctest::Approx(4 * std::atan(0.5)).epsilon(1e-15));
    const Arc half = std::get<Arc>(Bent({.X = 0, .Y = 0}, {.X = 10, .Y = 0}, 1));
    CHECK(half.Center            == Point{.X = 5, .Y = 0});
    CHECK(half.Radius            == 5.0);
    CHECK(half.Start.InRadians() == Pi);
    const Arc under = std::get<Arc>(Bent({.X = 0, .Y = 0}, {.X = 10, .Y = 0}, -0.5));
    CHECK(under.Center == Point{.X = 5, .Y = -3.75});
  }

  TEST_CASE("segment: a bulge under 2^-60, coincident ends or no bulge make a line")
  {
    /**
     *  An arc however slight its bulge: one of a thousand millionth rises half a thousandth from a
     *  chord of a million, and is an arc. Only under 2^-60, where it rises less than a step of a
     *  double, is it the chord.
     */
    CHECK(std::holds_alternative<Arc>(Bent({.X = 0, .Y = 0}, {.X = 1e6, .Y = 0}, 1e-6)));
    CHECK(std::holds_alternative<Arc>(Bent({.X = 0, .Y = 0}, {.X = 1e6, .Y = 0}, 9.99999e-7)));
    CHECK(std::holds_alternative<Arc>(Bent({.X = 0, .Y = 0}, {.X = 1e6, .Y = 0}, 1e-9)));
    CHECK(std::holds_alternative<Arc>(Bent({.X = 0, .Y = 0}, {.X = 1e6, .Y = 0}, 0x1p-60)));
    CHECK(std::holds_alternative<Line>(Bent({.X = 0, .Y = 0}, {.X = 1e6, .Y = 0}, 0x1p-61)));
    CHECK(std::holds_alternative<Line>(Bent({.X = 0, .Y = 0}, {.X = 1, .Y = 0}, 0)));
    CHECK(std::holds_alternative<Line>(Bent({.X = 1, .Y = 1}, {.X = 1, .Y = 1}, 1)));
    CHECK(std::holds_alternative<Line>(
        Bent({.X = 0, .Y = 0}, {.X = 1, .Y = 0}, std::numeric_limits<double>::quiet_NaN())));
  }

  TEST_CASE("segment: how many a polyline has, open and closed")
  {
    const Vertex a = {.Position = {.X = 0, .Y = 0}, .Bulge = 0};
    const Vertex b = {.Position = {.X = 1, .Y = 0}, .Bulge = 1};
    const Vertex c = {.Position = {.X = 1, .Y = 1}, .Bulge = 0};
    CHECK(Count(Polyline{}) == 0);
    CHECK(Count(Polyline{.Vertices = {a}, .Closed = true})        == 0);
    CHECK(Count(Polyline{.Vertices = {a, b, c}, .Closed = false}) == 2);
    CHECK(Count(Polyline{.Vertices = {a, b, c}, .Closed = true})  == 3);
    CHECK(Count(Polyline{.Vertices = {a, b}, .Closed = true})     == 2);
  }

  TEST_CASE("segment: the nth runs from its vertex to the next, the last back to the first")
  {
    const Polyline polyline = {
        .Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                     {.Position = {.X = 10, .Y = 0}, .Bulge = 1},
                     {.Position = {.X = 10, .Y = 10}, .Bulge = 0}},
        .Closed   = true,
    };
    CHECK(std::get<Line>(Nth(polyline, 0))
          == Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}});
    CHECK(std::get<Arc>(Nth(polyline, 1)).Center == Point{.X = 10, .Y = 5});
    CHECK(std::get<Line>(Nth(polyline, 2))
          == Line{.Start = {.X = 10, .Y = 10}, .End = {.X = 0, .Y = 0}});
  }

  TEST_CASE("segment: a nearly straight arc is read as the arc it is, wherever a polyline is read")
  {
    /**
     *  A chord of a thousand under a bulge of a millionth or a thousand millionth: the arc's middle
     *  lies half the bulge times the chord to its right, its length is the chord's by the angle its
     *  bulge is the tangent of over that tangent, and the point of it nearest one below its middle
     *  is its middle.
     */
    for (const double bulge : {1e-6, 1e-9})
    {
      CAPTURE(bulge);
      const Polyline bowed  = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = bulge},
                                            {.Position = {.X = 1000, .Y = 0}}},
                               .Closed   = false};
      const Point    middle = {.X = 500, .Y = -500 * bulge};
      CHECK(Distance(At(bowed, 0.5), middle) < 1e-12);
      CHECK(std::abs(Length(bowed) - (1000 * (1 + (bulge * bulge)) * std::atan(bulge) / bulge))
            < 1e-12);
      CHECK(Distance(Nearest(Curve{bowed}, Point{.X = 500, .Y = -10}).At, middle) < 1e-12);
    }
  }
} // namespace Geoxy
