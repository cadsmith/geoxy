/***************************************************************************************************
 *  @file      tests/geoxy/line.cpp
 *  @brief     Tests of lines, rays and xlines.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/line.hpp"

#include "geoxy/ray.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

namespace Geoxy
{
  TEST_CASE("line: equal end for end")
  {
    const Line line = {.Start = {.X = 1, .Y = 2}, .End = {.X = 3, .Y = 4}};
    CHECK(line       == Line{.Start = {.X = 1, .Y = 2}, .End = {.X = 3, .Y = 4}});
    CHECK_FALSE(line == Line{.Start = {.X = 0, .Y = 2}, .End = {.X = 3, .Y = 4}});
    CHECK_FALSE(line == Line{.Start = {.X = 1, .Y = 2}, .End = {.X = 3, .Y = 0}});
    CHECK(Line{}     == Line{.Start = {}, .End = {}});
  }

  TEST_CASE("ray: equal by start and way")
  {
    const Ray ray = {.Start = {.X = 1, .Y = 2}, .Way = {.X = 0, .Y = 1}};
    CHECK(ray       == Ray{.Start = {.X = 1, .Y = 2}, .Way = {.X = 0, .Y = 1}});
    CHECK_FALSE(ray == Ray{.Start = {.X = 0, .Y = 2}, .Way = {.X = 0, .Y = 1}});
    CHECK_FALSE(ray == Ray{.Start = {.X = 1, .Y = 2}, .Way = {.X = 0, .Y = 2}});
  }

  TEST_CASE("xline: equal by point and way")
  {
    const Xline xline = {.Through = {.X = 1, .Y = 2}, .Way = {.X = 1, .Y = 1}};
    CHECK(xline       == Xline{.Through = {.X = 1, .Y = 2}, .Way = {.X = 1, .Y = 1}});
    CHECK_FALSE(xline == Xline{.Through = {.X = 1, .Y = 0}, .Way = {.X = 1, .Y = 1}});
    CHECK_FALSE(xline == Xline{.Through = {.X = 1, .Y = 2}, .Way = {.X = -1, .Y = -1}});
  }
} // namespace Geoxy
