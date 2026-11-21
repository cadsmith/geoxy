/***************************************************************************************************
 *  @file      tests/geoxy/orient.cpp
 *  @brief     Tests of orientation.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/orient.hpp"

#include "geoxy/point.hpp"

#include <doctest/doctest.h>

namespace Geoxy
{
  TEST_CASE("orient: which way three points turn")
  {
    CHECK(Orient({.X = 0, .Y = 0}, {.X = 1, .Y = 0}, {.X = 0, .Y = 1})
          == Orientation::Counterclockwise);
    CHECK(Orient({.X = 0, .Y = 0}, {.X = 0, .Y = 1}, {.X = 1, .Y = 0}) == Orientation::Clockwise);
    CHECK(Orient({.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = 2}) == Orientation::Collinear);
    CHECK(Orient({.X = 1, .Y = 1}, {.X = 1, .Y = 1}, {.X = 5, .Y = 3}) == Orientation::Collinear);
  }

  /**
   *  Kettner et al.'s corner: a point a step of a double above the diagonal through two far ones.
   */
  TEST_CASE("orient: exact where the plain determinant is not")
  {
    const Point q     = {.X = 12, .Y = 12};
    const Point r     = {.X = 24, .Y = 24};
    const Point above = {.X = 0.5, .Y = 0.5 + 0x1.0p-53};
    const Point below = {.X = 0.5 + 0x1.0p-53, .Y = 0.5};
    CHECK(Orient(above, q, r) == Orientation::Counterclockwise);
    CHECK(Orient(q, r, above) == Orientation::Counterclockwise);
    CHECK(Orient(r, q, above) == Orientation::Clockwise);
    CHECK(Orient(below, q, r) == Orientation::Clockwise);
    CHECK(Orient({.X = 0.5, .Y = 0.5}, q, r) == Orientation::Collinear);
  }
} // namespace Geoxy
