/***************************************************************************************************
 *  @file      tests/geoxy/bezier/cubic.cpp
 *  @brief     Tests of cubic Bézier curves.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/bezier/cubic.hpp"

#include "geoxy/point.hpp"

#include <doctest/doctest.h>

namespace Geoxy::Bezier
{
  TEST_CASE("cubic: equal point for point")
  {
    const Cubic cubic = {
        .Start = {.X = 0, .Y = 0},
        .Out   = {.X = 1, .Y = 2},
        .In    = {.X = 3, .Y = 2},
        .End   = {.X = 4, .Y = 0},
    };
    Cubic other = cubic;
    CHECK(cubic == other);
    other.Start.Y = 1;
    CHECK_FALSE(cubic == other);
    other       = cubic;
    other.Out.X = 0;
    CHECK_FALSE(cubic == other);
    other      = cubic;
    other.In.Y = 0;
    CHECK_FALSE(cubic == other);
    other       = cubic;
    other.End.X = 5;
    CHECK_FALSE(cubic == other);
  }
} // namespace Geoxy::Bezier
