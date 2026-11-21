/***************************************************************************************************
 *  @file      tests/geoxy/vector.cpp
 *  @brief     Tests of vectors.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/vector.hpp"

#include "draws.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <functional>
#include <limits>
#include <optional>

namespace Geoxy
{
  static_assert(Vector{.X = 1, .Y = 2} + Vector{.X = 3, .Y = 4} == Vector{.X = 4, .Y = 6});
  static_assert(Vector{.X = 1, .Y = 2} - Vector{.X = 3, .Y = 5} == Vector{.X = -2, .Y = -3});
  static_assert(-Vector{.X = 1, .Y = -2}     == Vector{.X = -1, .Y = 2});
  static_assert(Vector{.X = 1, .Y = 2} * 3.0 == Vector{.X = 3, .Y = 6});
  static_assert(3.0 * Vector{.X = 1, .Y = 2} == Vector{.X = 3, .Y = 6});
  static_assert(Vector{.X = 3, .Y = 6} / 3.0 == Vector{.X = 1, .Y = 2});
  static_assert(Perpendicular(Vector{.X = 1, .Y = 2}) == Vector{.X = -2, .Y = 1});
  static_assert(std::equal_to<double>{}(Largest(Vector{.X = -3, .Y = 2}), 3.0));
  static_assert(std::equal_to<double>{}(Largest(Vector{.X = 1, .Y = -4}), 4.0));
  static_assert(Vector{.X = 0.0, .Y = 1} == Vector{.X = -0.0, .Y = 1});

  TEST_CASE("vector: the arithmetic at run time, and equality component by component")
  {
    const Vector one = {.X = 1, .Y = 2};
    const Vector two = {.X = 3, .Y = 4};
    CHECK(one + two == Vector{.X = 4, .Y = 6});
    CHECK(two - one == Vector{.X = 2, .Y = 2});
    CHECK(-one      == Vector{.X = -1, .Y = -2});
    CHECK(one * 3.0 == Vector{.X = 3, .Y = 6});
    CHECK(3.0 * one == Vector{.X = 3, .Y = 6});
    CHECK(two / 2.0 == Vector{.X = 1.5, .Y = 2});
    CHECK(Perpendicular(one)         == Vector{.X = -2, .Y = 1});
    CHECK(Largest({.X = -3, .Y = 2}) == 3.0);
    CHECK(Largest({.X = 1, .Y = -4}) == 4.0);
    CHECK_FALSE(one == Vector{.X = 0, .Y = 2});
    CHECK_FALSE(one == Vector{.X = 1, .Y = 0});
  }

  TEST_CASE("vector: dot and cross products")
  {
    CHECK(Dot({.X = 1, .Y = 2}, {.X = 3, .Y = 4})   == 11.0);
    CHECK(Cross({.X = 1, .Y = 2}, {.X = 3, .Y = 4}) == -2.0);
    CHECK(Cross({.X = 1, .Y = 0}, {.X = 0, .Y = 1}) == 1.0);
    CHECK(Cross({.X = 2, .Y = 2}, {.X = 1, .Y = 1}) == 0.0);
  }

  TEST_CASE("vector: a length at every scale")
  {
    CHECK(Length({.X = 3, .Y = 4})           == 5.0);
    CHECK(Length({.X = 3e6, .Y = 4e6})       == 5e6);
    CHECK(Length({.X = 3e300, .Y = 4e300})   == doctest::Approx(5e300).epsilon(1e-15));
    CHECK(Length({.X = 3e-300, .Y = 4e-300}) == doctest::Approx(5e-300).epsilon(1e-15));
    CHECK(Length({.X = std::numeric_limits<double>::max(), .Y = 0})
          == std::numeric_limits<double>::max());
    CHECK(Length({.X = 0x1.0p-1074, .Y = 0}) == 0x1.0p-1074);
    CHECK(Length(Vector{}) == 0.0);
    CHECK(std::isinf(Length({.X = std::numeric_limits<double>::infinity(), .Y = 0})));
  }

  TEST_CASE("vector: a unit vector has length 1 at every scale")
  {
    CHECK(Unit({.X = 1e300, .Y = 0})        == Vector{.X = 1, .Y = 0});
    CHECK(Unit({.X = 0, .Y = -0x1.0p-1074}) == Vector{.X = 0, .Y = -1});
    CHECK(Unit({.X = 3, .Y = 4})            == Vector{.X = 0.6, .Y = 0.8});
    Draws draws{0x756E6974U};
    for (const double scale : {1e-300, 1e-20, 1.0, 5e6, 1e8, 1e20, 1e300})
    {
      for (int index = 0; index < 1000; ++index)
      {
        const std::optional<Vector> unit =
            Unit({.X = draws.Near(0, scale), .Y = draws.Near(0, scale)});
        REQUIRE(unit);
        CHECK(Length(*unit) == doctest::Approx(1.0).epsilon(1e-15));
      }
    }
  }

  TEST_CASE("vector: no unit vector for no displacement or a non-finite one")
  {
    CHECK_FALSE(Unit({}));
    CHECK_FALSE(Unit({.X = std::numeric_limits<double>::infinity(), .Y = 0}));
    CHECK_FALSE(Unit({.X = 1, .Y = std::numeric_limits<double>::infinity()}));
    CHECK_FALSE(Unit({.X = 1, .Y = std::numeric_limits<double>::quiet_NaN()}));
  }
} // namespace Geoxy
