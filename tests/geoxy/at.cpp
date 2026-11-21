/***************************************************************************************************
 *  @file      tests/geoxy/at.cpp
 *  @brief     Tests of the points of curves.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/at.hpp"

#include "curves.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <limits>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether two points are within a distance.
     *
     *  @param  lhs      One point.
     *  @param  rhs      The other.
     *  @param  distance The distance.
     *  @return          Whether they are.
     */
    auto Near(Point lhs, Point rhs, double distance) -> bool
    {
      return Distance(lhs, rhs) <= distance;
    }
  } // namespace

  TEST_CASE("at: a line, a ray and an xline")
  {
    const Line line = {.Start = {.X = 5e6, .Y = 1}, .End = {.X = 5e6 + 10, .Y = 2}};
    CHECK(At(line, 0)   == line.Start);
    CHECK(At(line, 1)   == line.End);
    CHECK(At(line, 0.5) == Point{.X = 5e6 + 5, .Y = 1.5});
    CHECK(At(Ray{.Start = {.X = 1, .Y = 1}, .Way = {.X = 2, .Y = 0}}, 3) == Point{.X = 7, .Y = 1});
    CHECK(At(Xline{.Through = {.X = 1, .Y = 1}, .Way = {.X = 0, .Y = 2}}, -1)
          == Point{.X = 1, .Y = -1});
  }

  TEST_CASE("at: an arc, exact at its quarter turns, and run on round its circle")
  {
    const Arc circle = Circle({.X = 5e6, .Y = -3e6}, 2);
    CHECK(At(circle, 0)    == Point{.X = 5e6 + 2, .Y = -3e6});
    CHECK(At(circle, 0.25) == Point{.X = 5e6, .Y = -3e6 + 2});
    CHECK(At(circle, 0.5)  == Point{.X = 5e6 - 2, .Y = -3e6});
    const Arc arc = {
        .Center = {.X = 0, .Y = 0}, .Radius = 1, .Start = Degrees(90), .Sweep = Degrees(-90)};
    CHECK(At(arc, 1) == Point{.X = 1, .Y = 0});
    CHECK(At(arc, 2) == Point{.X = 0, .Y = -1});
  }

  TEST_CASE("at: an ellipse, the unit circle mapped")
  {
    const Ellipse ellipse = {
        .Map = Translation({.X = 10, .Y = 0}) * Scaling(2, 1), .Start = {}, .Sweep = Turn};
    CHECK(At(ellipse, 0)    == Point{.X = 12, .Y = 0});
    CHECK(At(ellipse, 0.25) == Point{.X = 10, .Y = 1});
    CHECK(Near(At(ellipse, 1.0 / 8), {.X = 10 + std::sqrt(2.0), .Y = std::sqrt(0.5)}, 1e-15));
  }

  TEST_CASE("at: a polyline, exact at its vertices, its bulges along their chords")
  {
    const Polyline polyline = {
        .Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                     {.Position = {.X = 10, .Y = 0}, .Bulge = 0},
                     {.Position = {.X = 10, .Y = 10}, .Bulge = -1}},
        .Closed   = true,
    };
    CHECK(At(polyline, 0) == Point{.X = 0, .Y = 0});
    CHECK(At(polyline, 1) == Point{.X = 10, .Y = 0});
    CHECK(At(polyline, 2) == Point{.X = 10, .Y = 10});
    CHECK(At(polyline, 3) == Point{.X = 0, .Y = 0});
    CHECK(Near(At(polyline, 0.5), {.X = 5, .Y = -5}, 1e-14));
    CHECK(At(polyline, 1.5) == Point{.X = 10, .Y = 5});
    CHECK(Near(At(polyline, 2.5), {.X = 10, .Y = 0}, 1e-13));
    CHECK(Near(At(polyline, -0.5), {.X = 5, .Y = 5}, 1e-14));
    CHECK(At(Polyline{}, 1) == Point{});
    CHECK(
        At(Polyline{.Vertices = {{.Position = {.X = 3, .Y = 4}, .Bulge = 1}}, .Closed = false}, 0.5)
        == Point{.X = 3, .Y = 4});
    CHECK(std::isnan(At(polyline, std::numeric_limits<double>::quiet_NaN()).X));
    const Polyline doubled = {.Vertices = {{.Position = {.X = 1, .Y = 1}, .Bulge = 1},
                                           {.Position = {.X = 1, .Y = 1}, .Bulge = 0}},
                              .Closed   = false};
    CHECK(At(doubled, 0.5) == Point{.X = 1, .Y = 1});
  }

  /**
   *  At 5e6 a chord of 1000 bulged 1e-6 has its center 2.5e8 away: the point worked along the chord
   *  stays on the arc to the rounding of the vertices, where one worked from the center does not.
   */
  TEST_CASE("at: a small bulge far from the origin keeps its precision")
  {
    const Polyline polyline = {.Vertices = {{.Position = {.X = 5e6, .Y = 5e6}, .Bulge = 1e-6},
                                            {.Position = {.X = 5e6 + 1000, .Y = 5e6}, .Bulge = 0}},
                               .Closed   = false};
    const Point    middle   = At(polyline, 0.5);
    CHECK(middle.X == doctest::Approx(5e6 + 500).epsilon(1e-15));
    CHECK(std::abs(middle.Y - (5e6 - 0.0005)) <= 4 * 0x1.0p-30);
  }

  TEST_CASE("at: a spline equals its Bézier curve, its basis sums and its circle")
  {
    for (const double along : {0.0, 0.1, 0.25, 0.5, 0.9, 1.0})
    {
      CHECK(Near(At(Tests::Clamped(), along), At(Tests::Wave, along), 1e-14));
      CHECK(std::abs(Distance(At(Tests::Quarter(), along), {}) - 1.0) <= 1e-15);
    }
    CHECK(At(Tests::Clamped(), 0) == Tests::Wave.Start);
    CHECK(At(Tests::Clamped(), 1) == Tests::Wave.End);
    for (const double along : {0.0, 0.1, 0.3, 0.7, 1.1, 1.5, 1.99})
    {
      CHECK(Near(At(Tests::Uneven(), along), Tests::Summed(Tests::Uneven(), along), 1e-13));
    }
    CHECK(Near(At(Tests::Uneven(), 2), {.X = 8, .Y = 1}, 1e-14));
  }

  TEST_CASE("at: a spline whose last knots repeat ends on its last span that is not empty")
  {
    const Spline spline = {.Degree  = 1,
                           .Knots   = {0, 0, 1, 1, 1},
                           .Poles   = {{.X = 0, .Y = 0}, {.X = 2, .Y = 2}, {.X = 5, .Y = 0}},
                           .Weights = {}};
    CHECK(At(spline, 1)   == Point{.X = 2, .Y = 2});
    CHECK(At(spline, 0.5) == Point{.X = 1, .Y = 1});
    CHECK(Near(At(Tests::Clamped(), 1.5), At(Tests::Wave, 1.5), 1e-12));
    CHECK(Near(At(Tests::Clamped(), -0.5), At(Tests::Wave, -0.5), 1e-12));
  }

  TEST_CASE("at: Bézier curves, exact at their ends")
  {
    const Bezier::Quadratic quadratic = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 2, .Y = 0}};
    CHECK(At(quadratic, 0)   == quadratic.Start);
    CHECK(At(quadratic, 1)   == quadratic.End);
    CHECK(At(quadratic, 0.5) == Point{.X = 1, .Y = 1});
    CHECK(At(Tests::Wave, 0) == Tests::Wave.Start);
    CHECK(At(Tests::Wave, 1) == Tests::Wave.End);
  }

  TEST_CASE("at: any curve, by its kind")
  {
    CHECK(At(Curve{Line{.Start = {}, .End = {.X = 2, .Y = 0}}}, 0.5) == Point{.X = 1, .Y = 0});
    CHECK(At(Curve{Circle({}, 1)}, 0.25) == Point{.X = 0, .Y = 1});
  }
} // namespace Geoxy
