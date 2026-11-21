/***************************************************************************************************
 *  @file      tests/geoxy/polyline.cpp
 *  @brief     Tests of polylines.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/polyline.hpp"

#include "geoxy/point.hpp"

#include <doctest/doctest.h>

namespace Geoxy
{
  TEST_CASE("polyline: equal vertex for vertex, and by whether it closes")
  {
    const Polyline polyline = {
        .Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0.5},
                     {.Position = {.X = 10, .Y = 0}, .Bulge = 0}},
        .Closed   = false,
    };
    Polyline other = polyline;
    CHECK(polyline == other);
    other.Closed = true;
    CHECK_FALSE(polyline == other);
    other                   = polyline;
    other.Vertices[0].Bulge = -0.5;
    CHECK_FALSE(polyline == other);
    other                        = polyline;
    other.Vertices[1].Position.Y = 1;
    CHECK_FALSE(polyline == other);
    other = polyline;
    other.Vertices.pop_back();
    CHECK_FALSE(polyline == other);
  }
} // namespace Geoxy
