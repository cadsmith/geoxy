/***************************************************************************************************
 *  @file      tests/geoxy/inertia.cpp
 *  @brief     Tests of second moments against Green's theorem integrated by mpmath at 50 digits,
 *             and against closed forms.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/inertia.hpp"

#include "curves.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/area.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/fit.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>

namespace Geoxy
{
  namespace
  {
    /**
     *  Moments checked against a reference: the area, the centroid and the three moments, each
     *  within a slack relative to its own size.
     *
     *  @param moments   The moments.
     *  @param reference The area, the centroid's coordinates, and the moments about x and y and
     *                   their product, as mpmath gives them.
     *  @param slack     How far each may be off, relative.
     */
    void Matched(const std::optional<Moments> &moments, const std::array<double, 6> &reference,
                 double slack)
    {
      REQUIRE(moments.has_value());
      const double                area   = std::abs(reference[0]);
      const double                size   = std::sqrt(area);
      const double                scale  = area * area;
      const std::array<double, 6> found  = {moments->Area,       moments->Centroid.X,
                                            moments->Centroid.Y, moments->AboutX,
                                            moments->AboutY,     moments->Product};
      const std::array<double, 6> slacks = {
          area, size + std::abs(reference[1]), size + std::abs(reference[2]), scale, scale, scale};
      for (std::size_t index = 0; index < found.size(); ++index)
      {
        CHECK(std::abs(found.at(index) - reference.at(index)) <= slack * slacks.at(index));
      }
    }
  } // namespace

  TEST_CASE("inertia: against mpmath")
  {
    for (const double offset : {0.0, 5e6})
    {
      const Polyline bulged = {
          .Vertices = {{.Position = {.X = offset, .Y = offset}, .Bulge = 0.5},
                       {.Position = {.X = offset + 4, .Y = offset}, .Bulge = 0},
                       {.Position = {.X = offset + 4, .Y = offset + 3}, .Bulge = 1e-4},
                       {.Position = {.X = offset, .Y = offset + 3}, .Bulge = -0.25}},
          .Closed   = true};
      Matched(Inertia(bulged, {}),
              {14.03683545527437, offset + 2.100029491803371, offset + 1.1201043100023407,
               16.99690847148106, 15.624270961821772, 0.5334105950897857},
              1e-14);
    }
    const Bezier::Cubic arch = {.Start = {.X = 0, .Y = 0},
                                .Out   = {.X = 1, .Y = 3},
                                .In    = {.X = 3, .Y = 3},
                                .End   = {.X = 4, .Y = 0}};
    Matched(Inertia(arch, {}),
            {-6.3, 2.0, 0.9183673469387755, -2.2080174927113703, -5.343506493506493, 0.0}, 1e-14);
    Matched(Inertia(Tests::Quarter(), {}),
            {0.28539816339744833, 0.5839794646280362, 0.5839794646280362, 0.015686296744689375,
             0.015686296744689375, -0.013996577438006037},
            1e-13);
    const Ellipse ellipse = {
        .Map   = Translation({.X = 5e6, .Y = -5e6}) * Rotation(Degrees(30)) * Scaling(3, 1),
        .Start = {},
        .Sweep = Turn};
    Matched(Inertia(ellipse, {}),
            {9.42477796076938, 5e6, -5e6, 7.0685834705770345, 16.493361431346415, 8.16209713905398},
            1e-13);
    const Arc arc = {
        .Center = {.X = 1, .Y = 1}, .Radius = 2, .Start = Degrees(20), .Sweep = Degrees(300)};
    Matched(Inertia(arc, {}),
            {12.204026319534854, 0.9462031230662602, 1.0094858408812712, 12.456784839851176,
             11.337499676378998, 0.2036932415798266},
            1e-14);
  }

  TEST_CASE("inertia: closed forms, and a flat bulge")
  {
    const double radius = 3;
    Matched(Inertia(Circle({.X = -7, .Y = 2}, radius), {}),
            {std::numbers::pi * 9, -7, 2, std::numbers::pi * 81 / 4, std::numbers::pi * 81 / 4, 0},
            1e-14);
    const Polyline box = {.Vertices = {{.Position = {.X = 1e8, .Y = 1e8}},
                                       {.Position = {.X = 1e8 + 4, .Y = 1e8}},
                                       {.Position = {.X = 1e8 + 4, .Y = 1e8 + 2}},
                                       {.Position = {.X = 1e8, .Y = 1e8 + 2}}},
                          .Closed   = false};
    Matched(Inertia(box, {}), {8, 1e8 + 2, 1e8 + 1, 4.0 * 8 / 12, 2.0 * 64 / 12, 0}, 1e-15);

    const double   bulge = 1e-5;
    const Polyline lens  = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = bulge},
                                         {.Position = {.X = 2, .Y = 0}, .Bulge = bulge}},
                            .Closed   = true};
    const std::optional<Moments> thin = Inertia(lens, {});
    REQUIRE(thin.has_value());
    const double area = 2 * (4.0 * bulge / 3);
    CHECK(std::abs(thin->Area - area) <= 1e-8 * area);
    CHECK(thin->AboutY > 0);
    CHECK(std::abs(thin->AboutY - (area * 0.2)) <= 1e-6 * area);

    const Principal axes = Axes(*Inertia(box, {}));
    CHECK(std::abs(axes.Least - (8.0 / 3)) <= 1e-14);
    CHECK(std::abs(axes.Most - (32.0 / 3)) <= 1e-14);
    CHECK(std::abs(axes.Axis.X - 1) <= 1e-15);
    Polyline backward = box;
    std::ranges::reverse(backward.Vertices);
    const std::optional<Moments> clockwise = Inertia(backward, {});
    REQUIRE(clockwise.has_value());
    CHECK(clockwise->Area < 0);
    CHECK(std::abs(Axes(*clockwise).Least - (8.0 / 3)) <= 1e-14);
  }

  TEST_CASE("inertia: splines span by span, to an accuracy")
  {
    const std::array<Point, 5> loop   = {Point{.X = 0, .Y = 0}, Point{.X = 4, .Y = 1},
                                         Point{.X = 5, .Y = 4},
                                         Point{.X = 1, .Y = 5},
                                         Point{.X = -1, .Y = 2}};
    const auto                 spline = Fit(loop, {.Closed = true});
    REQUIRE(spline.has_value());
    const std::optional<Moments> fine   = Inertia(*spline, {});
    const std::optional<Moments> coarse = Inertia(*spline, {.Error = 1e-3});
    REQUIRE(fine.has_value());
    REQUIRE(coarse.has_value());
    CHECK(std::abs(fine->Area - coarse->Area) <= 1e-3);
    CHECK(std::abs(fine->AboutX - coarse->AboutX) <= 1e-1);
    Matched(Inertia(Tests::Uneven(), {}), {fine->Area, 0, 0, 0, 0, 0}, 1e30);
  }

  TEST_CASE("inertia: a curve that closes on itself")
  {
    const Bezier::Cubic          drop    = {.Start = {.X = 1, .Y = 1},
                                            .Out   = {.X = 4, .Y = 4},
                                            .In    = {.X = -2, .Y = 4},
                                            .End   = {.X = 1, .Y = 1}};
    const std::optional<Moments> moments = Inertia(drop, {});
    REQUIRE(moments.has_value());
    CHECK(std::abs(moments->Area - Area(drop)) <= 1e-14 * std::abs(Area(drop)));
    CHECK(std::abs(moments->Centroid.X - 1) <= 1e-14);
    CHECK(std::abs(moments->Product) <= 1e-14);
  }

  TEST_CASE("inertia: what encloses nothing")
  {
    CHECK_FALSE(Inertia(Line{.Start = {}, .End = {.X = 1, .Y = 1}}, {}).has_value());
    CHECK_FALSE(Inertia(Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}, {}).has_value());
    CHECK_FALSE(Inertia(Xline{.Through = {}, .Way = {.X = 1, .Y = 0}}, {}).has_value());
    CHECK_FALSE(Inertia(Spline{}, {}).has_value());
    CHECK_FALSE(Inertia(Polyline{}, {}).has_value());
    const Polyline huge = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                        {.Position = {.X = 1e300, .Y = 0}},
                                        {.Position = {.X = 1e300, .Y = 1e300}}},
                           .Closed   = true};
    CHECK_FALSE(Inertia(huge, {}).has_value());
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(Inertia(Bezier::Cubic{.Start = {},
                                      .Out   = {.X = nan, .Y = 1},
                                      .In    = {.X = 2, .Y = 1},
                                      .End   = {.X = 3, .Y = 0}},
                        {})
                    .has_value());
  }
} // namespace Geoxy
