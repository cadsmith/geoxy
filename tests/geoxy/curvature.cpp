/***************************************************************************************************
 *  @file      tests/geoxy/curvature.cpp
 *  @brief     Tests of the unit tangent, the normal and the signed curvature of each curve: a
 *             circle's 1 / r exactly, an ellipse's against mpmath, and every kind's against the
 *             same curve spelled another way.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/curvature.hpp"

#include "curves.hpp"
#include "real.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/transform.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <optional>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether two directions are within a distance of each other.
     *
     *  @param  lhs      One, or none.
     *  @param  rhs      The other.
     *  @param  distance The distance.
     *  @return          Whether the first is one and they are.
     */
    auto Near(std::optional<Vector> lhs, Vector rhs, double distance) -> bool
    {
      return lhs && Length(*lhs - rhs) <= distance;
    }

    /**
     *  Whether a curvature is a value, within a part of it.
     *
     *  @param  curvature The curvature, or none.
     *  @param  value     The value.
     *  @param  part      The part of it.
     *  @return           Whether it is one and is.
     */
    auto Near(std::optional<double> curvature, double value, double part) -> bool
    {
      return curvature && std::abs(*curvature - value) <= part * std::abs(value);
    }
  } // namespace

  TEST_CASE("curvature: a circle's is 1 / r exactly, its tangent its radius turned")
  {
    for (const double radius : {0.1, 1.0, 3.0, 7.25, 1e-7, 4e6})
    {
      const Arc circle = Circle({.X = 5e6, .Y = -2}, radius);
      for (const double along : {0.0, 0.1, 0.25, 0.6})
      {
        const std::optional<double> curvature = Curvature(circle, along);
        REQUIRE(curvature);
        CHECK(Exactly(*curvature, 1.0 / radius));
        const Vector out = Direction(circle.Start + (circle.Sweep * along));
        CHECK(Tangent(circle, along) == Perpendicular(out));
        CHECK(Normal(circle, along)  == -out);
      }
    }
    const Arc backward = {
        .Center = {.X = 1, .Y = 1}, .Radius = 2, .Start = Degrees(90), .Sweep = Degrees(-90)};
    CHECK(Curvature(backward, 0.5) == -0.5);
    CHECK(Tangent(backward, 0)     == Vector{.X = 1, .Y = 0});
    CHECK(Normal(backward, 0)      == Vector{.X = 0, .Y = 1});
    const Arc inward = {
        .Center = {.X = 1, .Y = 1}, .Radius = -2, .Start = Degrees(0), .Sweep = Degrees(90)};
    CHECK(Curvature(inward, 0.5) == 0.5);
    CHECK(Tangent(inward, 0)     == Vector{.X = 0, .Y = -1});
    const Arc both = {
        .Center = {.X = 1, .Y = 1}, .Radius = -2, .Start = Degrees(0), .Sweep = Degrees(-90)};
    CHECK(Tangent(both, 0)   == Vector{.X = 0, .Y = 1});
    CHECK(Curvature(both, 0) == -0.5);
    CHECK_FALSE(Tangent(Circle({}, 0), 0.5));
    CHECK_FALSE(Curvature(Circle({}, 0), 0.5));
    CHECK_FALSE(Normal(Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = {}}, 0.5));
  }

  TEST_CASE("curvature: an ellipse's det M / |M t|^3, against mpmath")
  {
    const Ellipse wide = {.Map = Scaling(3, 1), .Start = Degrees(0), .Sweep = Degrees(360)};
    CHECK(Curvature(wide, 0) == 3.0);
    CHECK(Near(Curvature(wide, 0.25), 1.0 / 9, 1e-15));
    CHECK(Tangent(wide, 0)   == Vector{.X = 0, .Y = 1});
    CHECK(Normal(wide, 0.25) == Vector{.X = 0, .Y = -1});

    const Ellipse sheared = {
        .Map   = {.X = {.X = 2, .Y = 0.5}, .Y = {.X = -0.3, .Y = 1.2}, .Origin = {.X = -7, .Y = 3}},
        .Start = Radians(0.3),
        .Sweep = Radians(-2)};
    constexpr std::array<std::array<double, 2>, 5> Mpmath = {{
        {0.0, -1.0851380466292249326},
        {0.25, -1.2170967143835653079},
        {0.5, -0.5881670984330363741},
        {1.0, -0.29574343270163799987},
        {0.7, -0.35969699595559526835},
    }};
    for (const auto &[along, curvature] : Mpmath)
    {
      CHECK(Near(Curvature(sheared, along), curvature, 4e-16));
      CHECK(Near(Tangent(sheared, along), *Unit(Derivative(sheared, along)), 1e-15));
    }
    const Ellipse still = {.Map = Scaling(3, 1), .Start = Degrees(10), .Sweep = Degrees(0)};
    CHECK_FALSE(Tangent(still, 0.5));
    CHECK_FALSE(Curvature(still, 0.5));
    const Ellipse flat = {
        .Map = {.X = {.X = 2, .Y = 0}, .Y = {.X = 0, .Y = 0}}, .Start = {}, .Sweep = Turn};
    CHECK_FALSE(Tangent(flat, 0));
    CHECK_FALSE(Curvature(flat, 0));
    CHECK(Curvature(flat, 0.25) == 0.0);
  }

  TEST_CASE("curvature: lines, rays and xlines run straight, but where they stop")
  {
    const Line line = {.Start = {.X = 1, .Y = 1}, .End = {.X = 4, .Y = 5}};
    CHECK(Tangent(line, 0.3) == Vector{.X = 0.6, .Y = 0.8});
    CHECK(Normal(line, 0.3)  == Vector{.X = -0.8, .Y = 0.6});
    CHECK(Curvature(line, 2) == 0.0);
    CHECK_FALSE(Tangent(Line{.Start = {.X = 1, .Y = 1}, .End = {.X = 1, .Y = 1}}, 0));
    CHECK_FALSE(Curvature(Line{.Start = {.X = 1, .Y = 1}, .End = {.X = 1, .Y = 1}}, 0));
    CHECK(Tangent(Ray{.Start = {}, .Way = {.X = 0, .Y = -3}}, 9)       == Vector{.X = 0, .Y = -1});
    CHECK(Normal(Ray{.Start = {}, .Way = {.X = 0, .Y = -3}}, 9)        == Vector{.X = 1, .Y = 0});
    CHECK(Curvature(Ray{.Start = {}, .Way = {.X = 0, .Y = -3}}, 9)     == 0.0);
    CHECK(Tangent(Xline{.Through = {}, .Way = {.X = -2, .Y = 0}}, 9)   == Vector{.X = -1, .Y = 0});
    CHECK(Normal(Xline{.Through = {}, .Way = {.X = -2, .Y = 0}}, 9)    == Vector{.X = 0, .Y = -1});
    CHECK(Curvature(Xline{.Through = {}, .Way = {.X = -2, .Y = 0}}, 9) == 0.0);
  }

  TEST_CASE("curvature: a polyline's bulges by their circles, its lines flat")
  {
    const Polyline polyline = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                                            {.Position = {.X = 2, .Y = 0}, .Bulge = 0},
                                            {.Position = {.X = 2, .Y = 3}, .Bulge = -0.5}},
                               .Closed   = true};
    CHECK(Near(Curvature(polyline, 0.5), 1.0, 1e-15));
    CHECK(Near(Tangent(polyline, 0.5), {.X = 1, .Y = 0}, 1e-15));
    CHECK(Near(Normal(polyline, 0.5), {.X = 0, .Y = 1}, 1e-15));
    CHECK(Curvature(polyline, 1.5) == 0.0);
    CHECK(Tangent(polyline, 1.5)   == Vector{.X = 0, .Y = 1});
    const double chord = std::sqrt(13.0);
    CHECK(Near(Curvature(polyline, 2.25), -2.0 / (chord * 1.25), 1e-15));
    CHECK_FALSE(Tangent(Polyline{}, 0));
    CHECK_FALSE(Curvature(Polyline{}, 0));
  }

  TEST_CASE("curvature: Bézier curves and splines, a parabola's and the same curve's spelled twice")
  {
    const Bezier::Quadratic parabola = {
        .Start = {.X = -1, .Y = 1}, .Control = {.X = 0, .Y = -1}, .End = {.X = 1, .Y = 1}};
    CHECK(Curvature(parabola, 0.5) == 2.0);
    CHECK(Near(Curvature(parabola, 1), 0.17888543819998317571, 4e-16));
    CHECK(Tangent(parabola, 0.5) == Vector{.X = 1, .Y = 0});
    CHECK(Normal(parabola, 0.5)  == Vector{.X = 0, .Y = 1});

    const Spline clamped = Tests::Clamped();
    for (const double along : {0.0, 0.2, 0.5, 0.9, 1.0})
    {
      const std::optional<double> curvature = Curvature(Tests::Wave, along);
      REQUIRE(curvature);
      CHECK(Near(Curvature(clamped, along), *curvature, 1e-13));
      CHECK(Near(Tangent(clamped, along), *Tangent(Tests::Wave, along), 1e-15));
      CHECK(Near(Normal(clamped, along), *Normal(Tests::Wave, along), 1e-15));
    }
    const Bezier::Cubic stopped = {.Start = {.X = 0, .Y = 0},
                                   .Out   = {.X = 0, .Y = 0},
                                   .In    = {.X = 1, .Y = 1},
                                   .End   = {.X = 2, .Y = 0}};
    CHECK_FALSE(Tangent(stopped, 0));
    CHECK_FALSE(Curvature(stopped, 0));
    CHECK(Tangent(stopped, 0.5));

    for (const double along : {0.0, 0.3, 0.5, 1.0})
    {
      CHECK(Near(Curvature(Tests::Quarter(), along), 1.0, 1e-14));
      const Point point = At(Tests::Quarter(), along);
      CHECK(Near(Normal(Tests::Quarter(), along), Point{} - point, 1e-14));
    }
    const Spline uneven = Tests::Uneven();
    for (const double along : {0.1, 0.3, 0.7, 1.1, 1.9})
    {
      const std::optional<double> scaled = Curvature(Transformed(uneven, Scaling(2)), along);
      REQUIRE(scaled);
      CHECK(Near(Curvature(uneven, along), 2 * *scaled, 1e-12));
    }
    const Spline broken = {.Degree  = 1,
                           .Knots   = {0, 0, 1, 2, 2},
                           .Poles   = {{.X = 0, .Y = 0}, {.X = 3, .Y = 4}, {.X = 3, .Y = 0}},
                           .Weights = {}};
    CHECK(Curvature(broken, 0.5) == 0.0);
    CHECK(Tangent(broken, 0.5)   == Vector{.X = 0.6, .Y = 0.8});
    CHECK(Tangent(broken, 1.5)   == Vector{.X = 0, .Y = -1});
    CHECK_FALSE(Tangent(Spline{}, 0));
    CHECK_FALSE(Curvature(Spline{}, 0));
  }

  TEST_CASE("curvature: any curve, as its own kind")
  {
    const std::array<Curve, 3> curves = {
        Curve{Circle({}, 4)},
        Curve{Line{.Start = {}, .End = {.X = 0, .Y = 2}}},
        Curve{Tests::Wave},
    };
    for (const Curve &curve : curves)
    {
      CHECK(Tangent(curve, 0.5)
            == std::visit([](const auto &shape) -> std::optional<Vector>
                          { return Tangent(shape, 0.5); }, curve));
      CHECK(Normal(curve, 0.5)
            == std::visit([](const auto &shape) -> std::optional<Vector>
                          { return Normal(shape, 0.5); }, curve));
      CHECK(Curvature(curve, 0.5)
            == std::visit([](const auto &shape) -> std::optional<double>
                          { return Curvature(shape, 0.5); }, curve));
    }
  }
} // namespace Geoxy
