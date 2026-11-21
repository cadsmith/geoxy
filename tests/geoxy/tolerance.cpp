/***************************************************************************************************
 *  @file      tests/geoxy/tolerance.cpp
 *  @brief     Tests of tolerances and of the step of a double.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/tolerance.hpp"

#include "geoxy/box.hpp"
#include "geoxy/point.hpp"

#include <doctest/doctest.h>

#include <limits>

namespace Geoxy
{
  TEST_CASE("tolerance: of a drawing, 1e-10 of its larger side")
  {
    const Tolerance tolerance =
        Tolerance::Of({.Min = {.X = 0, .Y = 0}, .Max = {.X = 1000, .Y = 500}});
    CHECK(tolerance.Length  == doctest::Approx(1e-7).epsilon(1e-15));
    CHECK(tolerance.Radians == doctest::Approx(1e-10).epsilon(1e-15));
  }

  /**
   *  A survey drawing 100 units across at 5e6: 1e-10 of it, 1e-8, is under 16 steps of a double
   *  there.
   */
  TEST_CASE("tolerance: of a drawing at survey coordinates, no less than 16 steps of a double")
  {
    const Tolerance tolerance =
        Tolerance::Of({.Min = {.X = 5e6, .Y = 3e6}, .Max = {.X = 5e6 + 100, .Y = 3e6 + 100}});
    CHECK(tolerance.Length  == 16 * 0x1.0p-30);
    CHECK(tolerance.Radians == doctest::Approx(16 * 0x1.0p-30 / 100).epsilon(1e-15));
  }

  TEST_CASE("tolerance: of no extent, the steps alone")
  {
    const Tolerance point = Tolerance::Of({.Min = {.X = 1, .Y = 1}, .Max = {.X = 1, .Y = 1}});
    CHECK(point.Length  == 16 * 0x1.0p-52);
    CHECK(point.Radians == 1e-10);
    const Tolerance empty = Tolerance::Of({});
    CHECK(empty.Length  == 16 * 0x1.0p-1074);
    CHECK(empty.Radians == 1e-10);
  }

  /**
   *  A box with a side or a corner past every double, as the extents of a curve whose numbers run
   *  past them are: its length was past every double too, and its angle not a number.
   */
  TEST_CASE("tolerance: of a box past every double, the largest double, at the ratio's angle")
  {
    constexpr double Most    = std::numeric_limits<double>::max();
    constexpr double Endless = std::numeric_limits<double>::infinity();
    for (const Box box : {Box{.Min = {.X = 0, .Y = 0}, .Max = {.X = Endless, .Y = 1}},
                          Box{.Min = {.X = -Endless, .Y = 0}, .Max = {.X = 1, .Y = Endless}},
                          Box{.Min = {.X = -Most, .Y = 0}, .Max = {.X = Most, .Y = 1}},
                          Box{.Min = {.X = 0, .Y = Endless}, .Max = {.X = 1, .Y = Endless}},
                          Box{.Min = {.X = Endless, .Y = 0}, .Max = {.X = Endless, .Y = 1}}})
    {
      const Tolerance tolerance = Tolerance::Of(box);
      CHECK(tolerance.Length  == Most);
      CHECK(tolerance.Radians == 1e-10);
    }
  }

  TEST_CASE("tolerance: the same point within a length, its edge included")
  {
    CHECK(Same({.X = 0, .Y = 0}, {.X = 3, .Y = 4}, {.Length = 5}));
    CHECK_FALSE(Same({.X = 0, .Y = 0}, {.X = 3, .Y = 4}, {.Length = 4.9}));
    CHECK(Same({.X = 5e6, .Y = 1}, {.X = 5e6, .Y = 1}, {}));
    CHECK_FALSE(Same({.X = 5e6, .Y = 1}, {.X = 5e6 + 0x1.0p-30, .Y = 1}, {}));
  }

  TEST_CASE("tolerance: the step of a double")
  {
    CHECK(Grain(1.0)  == 0x1.0p-52);
    CHECK(Grain(-1.0) == 0x1.0p-52);
    CHECK(Grain(5e6)  == 0x1.0p-30);
    CHECK(Grain(0.0)  == 0x1.0p-1074);
    CHECK(Grain(Point{.X = 1, .Y = -5e6}) == 0x1.0p-30);
    CHECK(Grain(Point{.X = -5e6, .Y = 1}) == 0x1.0p-30);
  }
} // namespace Geoxy
