/***************************************************************************************************
 *  @file      tests/geoxy/derivative.cpp
 *  @brief     Tests of the derivatives of curves, against finite differences.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/derivative.hpp"

#include "curves.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <array>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether a curve's derivative matches its central finite difference at a value.
     *
     *  @param  curve The curve.
     *  @param  along The value.
     *  @return       Whether they agree to a part in 1e6 of the derivative's length, and one.
     */
    auto Agrees(const Curve &curve, double along) -> bool
    {
      constexpr double Step   = 1e-6;
      const Vector difference = (At(curve, along + Step) - At(curve, along - Step)) / (2 * Step);
      const Vector derivative = Derivative(curve, along);
      return Length(derivative - difference) <= 1e-6 * (Length(derivative) + 1.0);
    }
  } // namespace

  TEST_CASE("derivative: every kind agrees with its finite difference")
  {
    const Polyline polyline = {
        .Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0.5},
                     {.Position = {.X = 10, .Y = 0}, .Bulge = 0},
                     {.Position = {.X = 10, .Y = 10}, .Bulge = -2}},
        .Closed   = true,
    };
    const std::array<Curve, 11> curves = {
        Line{.Start = {.X = 1, .Y = 2}, .End = {.X = 4, .Y = -2}},
        Ray{.Start = {.X = 1, .Y = 2}, .Way = {.X = 3, .Y = 1}},
        Xline{.Through = {.X = 1, .Y = 2}, .Way = {.X = -1, .Y = 1}},
        Arc{.Center = {.X = 1, .Y = 2}, .Radius = 3, .Start = Degrees(30), .Sweep = Degrees(-250)},
        Ellipse{.Map   = Translation({.X = 1, .Y = 1}) * Rotation(Radians(0.4)) * Scaling(3, 1),
                .Start = Radians(0.2),
                .Sweep = Radians(4)},
        polyline,
        Tests::Clamped(),
        Tests::Quarter(),
        Tests::Uneven(),
        Bezier::Quadratic{
            .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 3, .Y = 1}},
        Tests::Wave,
    };
    for (const Curve &curve : curves)
    {
      for (const double along : {0.1, 0.35, 0.5, 0.8, 1.2, 1.7, 2.4})
      {
        CHECK(Agrees(curve, along));
      }
    }
  }

  TEST_CASE("derivative: the exact ones")
  {
    CHECK(Derivative(Line{.Start = {.X = 1, .Y = 2}, .End = {.X = 4, .Y = -2}}, 0.3)
          == Vector{.X = 3, .Y = -4});
    CHECK(Derivative(Ray{.Start = {}, .Way = {.X = 3, .Y = 1}}, 7) == Vector{.X = 3, .Y = 1});
    CHECK(Derivative(Xline{.Through = {}, .Way = {.X = -1, .Y = 1}}, -7)
          == Vector{.X = -1, .Y = 1});
    CHECK(Derivative(Circle({}, 2), 0) == Vector{.X = 0, .Y = 2 * 2 * Pi});
    CHECK(Derivative(Ellipse{.Map = Scaling(2, 1), .Start = {}, .Sweep = Turn}, 0.25)
          == Vector{.X = -2 * 2 * Pi, .Y = 0});
    CHECK(Derivative(Polyline{}, 0) == Vector{});
  }

  /**
   *  A spline of degree 1 whose doubled knot cuts it in two: the slope across the cut has no span,
   *  and counts for nothing.
   */
  TEST_CASE("derivative: a spline across a knot of full multiplicity")
  {
    const Spline spline = {
        .Degree  = 1,
        .Knots   = {0, 0, 1, 1, 2, 2},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 5, .Y = 1}, {.X = 6, .Y = 3}},
        .Weights = {},
    };
    CHECK(Derivative(spline, 0.5) == Vector{.X = 1, .Y = 1});
    CHECK(Derivative(spline, 1.5) == Vector{.X = 1, .Y = 2});
  }

  /**
   *  A rational spline's weight sums to nothing where its run-on past the domain passes the line at
   *  infinity: there it is at no point and has no tangent, where the slope was divided by nothing.
   *  Weights 1e-170 apart from nothing are levelled first, and their tangent is the spline's own.
   */
  TEST_CASE("derivative: none where a spline's weight sums to nothing")
  {
    const Spline beyond = {
        .Degree  = 1,
        .Knots   = {0, 0, 1, 1},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 0}},
        .Weights = {1, 3},
    };
    CHECK(Derivative(beyond, -0.5) == Vector{});
    CHECK(Derivative(beyond, 0.5)  == Vector{.X = 0.75, .Y = 0});
    const Spline light = {
        .Degree  = 1,
        .Knots   = {0, 0, 1, 1},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}},
        .Weights = {1e-170, 1e-170},
    };
    CHECK(Derivative(light, 0.5) == Vector{.X = 1, .Y = 1});
  }
} // namespace Geoxy
