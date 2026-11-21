/***************************************************************************************************
 *  @file      tests/geoxy/point.cpp
 *  @brief     Tests of points.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/point.hpp"

#include "draws.hpp"

#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <limits>

namespace Geoxy
{
  static_assert(Point{.X = 5, .Y = 7} - Point{.X = 2, .Y = 3}  == Vector{.X = 3, .Y = 4});
  static_assert(Point{.X = 2, .Y = 3} + Vector{.X = 3, .Y = 4} == Point{.X = 5, .Y = 7});
  static_assert(Point{.X = 5, .Y = 7} - Vector{.X = 3, .Y = 4} == Point{.X = 2, .Y = 3});

  TEST_CASE("point: the arithmetic at run time, and equality coordinate by coordinate")
  {
    const Point start = {.X = 2, .Y = 3};
    const Point end   = {.X = 5, .Y = 7};
    CHECK(end - start == Vector{.X = 3, .Y = 4});
    CHECK(start + Vector{.X = 3, .Y = 4} == end);
    CHECK(end - Vector{.X = 3, .Y = 4}   == start);
    CHECK_FALSE(start == Point{.X = 0, .Y = 3});
    CHECK_FALSE(start == Point{.X = 2, .Y = 0});
  }

  TEST_CASE("point: halfway between")
  {
    CHECK(Midpoint({.X = 0, .Y = 0}, {.X = 2, .Y = 4}) == Point{.X = 1, .Y = 2});
    CHECK(Midpoint({.X = 5e6, .Y = -3e6}, {.X = 5e6 + 2, .Y = -3e6 + 4})
          == Point{.X = 5e6 + 1, .Y = -3e6 + 2});
    constexpr double Huge = std::numeric_limits<double>::max();
    CHECK(Midpoint({.X = Huge, .Y = -Huge}, {.X = Huge, .Y = -Huge})
          == Point{.X = Huge, .Y = -Huge});

    /**
     *  Below the normal doubles a half rounds: a point with an odd last step is still its own
     *  midpoint, and two a step apart meet on the even one.
     */
    const Point odd = {.X = 0x1.6c78c51295bfcp-15, .Y = -0x1.bdp-1066};
    CHECK(Midpoint(odd, odd) == odd);
    CHECK(Midpoint({.X = 0x1p-1074, .Y = Huge}, {.X = 0x2p-1074, .Y = -Huge})
          == Point{.X = 0x2p-1074, .Y = 0});
    CHECK(Midpoint({.X = 0x1p-1074, .Y = 1}, {.X = Huge, .Y = 0x1p-1074})
          == Point{.X = Huge / 2, .Y = 0.5});
    CHECK(Midpoint({.X = Huge, .Y = 0x1p-1074}, {.X = 0x1p-1074, .Y = Huge})
          == Point{.X = Huge / 2, .Y = Huge / 2});
  }

  TEST_CASE("point: a fraction of the way is exact at both ends")
  {
    Draws draws{0x6C657270U};
    for (const double offset : {0.0, 5e6, 1e8})
    {
      for (int index = 0; index < 1000; ++index)
      {
        const Point start = {.X = draws.Near(offset, 1e3), .Y = draws.Near(-offset, 1e3)};
        const Point end   = {.X = draws.Near(offset, 1e3), .Y = draws.Near(-offset, 1e3)};
        CHECK(Lerp(start, end, 0.0) == start);
        CHECK(Lerp(start, end, 1.0) == end);
      }
    }
    CHECK(Lerp({.X = 0, .Y = 0}, {.X = 4, .Y = 8}, 0.25) == Point{.X = 1, .Y = 2});
    CHECK(Lerp({.X = 0, .Y = 0}, {.X = 4, .Y = 8}, 0.75) == Point{.X = 3, .Y = 6});
    CHECK(Lerp({.X = 0, .Y = 0}, {.X = 4, .Y = 8}, 2.0)  == Point{.X = 8, .Y = 16});
  }

  TEST_CASE("point: the distance between two points at every scale")
  {
    CHECK(Distance({.X = 0, .Y = 0}, {.X = 3, .Y = 4}) == 5.0);
    CHECK(Distance({.X = 5e6, .Y = -3e6}, {.X = 5e6 + 3, .Y = -3e6 + 4}) == 5.0);
    CHECK(Distance({.X = 1e8, .Y = 1e8}, {.X = 1e8 + 3, .Y = 1e8 - 4})   == 5.0);
    CHECK(Distance({.X = 1, .Y = 1}, {.X = 1, .Y = 1}) == 0.0);
  }
} // namespace Geoxy
