/***************************************************************************************************
 *  @file      tests/geoxy/ends.cpp
 *  @brief     Tests of where a curve starts and ends: both for a bounded curve, the start alone for
 *             a ray, none for an xline or a curve of nothing.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/ends.hpp"

#include "curves.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/line.hpp"
#include "geoxy/place.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

namespace Geoxy
{
  TEST_CASE("ends: both of a bounded curve, at its domain's ends")
  {
    const Line line = {.Start = {.X = 1, .Y = 2}, .End = {.X = 3, .Y = 4}};
    CHECK(Ends(line)
          == Endpoints{.Start = Place{.At = line.Start, .Along = 0},
                       .End   = Place{.At = line.End, .Along = 1}});
    const Arc arc = {
        .Center = {.X = 0, .Y = 0}, .Radius = 2, .Start = Degrees(90), .Sweep = Degrees(90)};
    CHECK(Ends(arc).End             == Place{.At = {.X = -2, .Y = 0}, .Along = 1});
    CHECK(Ends(Tests::Uneven()).End == Place{.At = {.X = 8, .Y = 1}, .Along = 2});
    CHECK(Ends(Tests::Wave).Start   == Place{.At = {.X = 0, .Y = 0}, .Along = 0});

    const Polyline closed = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                          {.Position = {.X = 2, .Y = 0}, .Bulge = 0.5},
                                          {.Position = {.X = 2, .Y = 3}, .Bulge = 0}},
                             .Closed   = true};
    CHECK(Ends(closed).End == Place{.At = {.X = 0, .Y = 0}, .Along = 3});
    const Polyline point = {.Vertices = {{.Position = {.X = 5, .Y = 6}, .Bulge = 0}},
                            .Closed   = false};
    CHECK(Ends(point).Start == Ends(point).End);
  }

  TEST_CASE("ends: a ray's start alone, none of an xline or of a curve of nothing")
  {
    const Endpoints ray = Ends(Ray{.Start = {.X = 1, .Y = 1}, .Way = {.X = 1, .Y = 0}});
    CHECK(ray.Start == Place{.At = {.X = 1, .Y = 1}, .Along = 0});
    CHECK_FALSE(ray.End);
    CHECK(Ends(Xline{.Through = {}, .Way = {.X = 1, .Y = 0}}) == Endpoints{});
    CHECK(Ends(Polyline{}) == Endpoints{});
    CHECK(Ends(Spline{})   == Endpoints{});
    CHECK_FALSE(ray        == Endpoints{});
    CHECK_FALSE(ray        == Endpoints{.Start = ray.Start, .End = ray.Start});
  }
} // namespace Geoxy
