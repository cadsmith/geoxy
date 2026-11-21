/***************************************************************************************************
 *  @file      tests/geoxy/domain.cpp
 *  @brief     Tests of the ranges of curves' parameters.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/domain.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <limits>

namespace Geoxy
{
  TEST_CASE("domain: each kind's range")
  {
    constexpr double Infinity = std::numeric_limits<double>::infinity();
    const Range      unit     = {.Start = 0, .End = 1};
    CHECK(Domain(Curve{Line{}})    == unit);
    CHECK(Domain(Curve{Arc{}})     == unit);
    CHECK(Domain(Curve{Ellipse{}}) == unit);
    CHECK(Domain(Curve{Bezier::Quadratic{}}) == unit);
    CHECK(Domain(Curve{Bezier::Cubic{}})     == unit);
    CHECK(Domain(Curve{Ray{}})   == Range{.Start = 0, .End = Infinity});
    CHECK(Domain(Curve{Xline{}}) == Range{.Start = -Infinity, .End = Infinity});
    const Polyline polyline = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                            {.Position = {.X = 1, .Y = 0}, .Bulge = 0},
                                            {.Position = {.X = 1, .Y = 1}, .Bulge = 0}},
                               .Closed   = true};
    CHECK(Domain(Curve{polyline}) == Range{.Start = 0, .End = 3});
    const Spline spline = {
        .Degree  = 2,
        .Knots   = {0, 0, 0, 0.5, 2, 2, 2},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = 0}, {.X = 3, .Y = 1}},
        .Weights = {}};
    CHECK(Domain(Curve{spline}) == Range{.Start = 0, .End = 2});
    CHECK_FALSE(unit            == Range{.Start = 0, .End = 2});
    CHECK_FALSE(unit            == Range{.Start = 1, .End = 1});
  }
} // namespace Geoxy
