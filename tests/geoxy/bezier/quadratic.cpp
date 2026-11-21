/***************************************************************************************************
 *  @file      tests/geoxy/bezier/quadratic.cpp
 *  @brief     Tests of quadratic Bézier curves.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/bezier/quadratic.hpp"

#include "geoxy/point.hpp"

#include <doctest/doctest.h>

namespace Geoxy::Bezier
{
  TEST_CASE("quadratic: equal point for point")
  {
    const Quadratic quadratic = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 2, .Y = 0}};
    Quadratic other = quadratic;
    CHECK(quadratic == other);
    other.Start.X = 1;
    CHECK_FALSE(quadratic == other);
    other           = quadratic;
    other.Control.Y = 1;
    CHECK_FALSE(quadratic == other);
    other       = quadratic;
    other.End.X = 3;
    CHECK_FALSE(quadratic == other);
  }
} // namespace Geoxy::Bezier
