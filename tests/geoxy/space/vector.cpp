/***************************************************************************************************
 *  @file      tests/geoxy/space/vector.cpp
 *  @brief     Tests of vectors in space.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/space/vector.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <limits>

namespace Geoxy::Space
{
  static_assert(East == Vector{.X = 1, .Y = 0, .Z = 0});
  static_assert(Up != North);

  TEST_CASE("space vector: equality component by component")
  {
    const Vector vector = {.X = 1, .Y = 2, .Z = 3};
    CHECK(vector       == Vector{.X = 1, .Y = 2, .Z = 3});
    CHECK_FALSE(vector == Vector{.X = 0, .Y = 2, .Z = 3});
    CHECK_FALSE(vector == Vector{.X = 1, .Y = 0, .Z = 3});
    CHECK_FALSE(vector == Vector{.X = 1, .Y = 2, .Z = 0});
  }

  TEST_CASE("space vector: dot and cross products, right-handed")
  {
    CHECK(Dot({.X = 1, .Y = 2, .Z = 3}, {.X = 4, .Y = 5, .Z = 6}) == 32.0);
    CHECK(Cross(East, North) == Up);
    CHECK(Cross(North, Up)   == East);
    CHECK(Cross(Up, East)    == North);
    CHECK(Cross({.X = 1, .Y = 2, .Z = 3}, {.X = 4, .Y = 5, .Z = 6})
          == Vector{.X = -3, .Y = 6, .Z = -3});
  }

  TEST_CASE("space vector: a length at every scale")
  {
    CHECK(Length({.X = 2, .Y = 3, .Z = 6})             == 7.0);
    CHECK(Length({.X = 2e300, .Y = 3e300, .Z = 6e300}) == doctest::Approx(7e300).epsilon(1e-15));
    CHECK(Length({.X = 2e-300, .Y = 3e-300, .Z = 6e-300})
          == doctest::Approx(7e-300).epsilon(1e-15));
    CHECK(Length({}) == 0.0);
    CHECK(std::isinf(Length({.X = 0, .Y = 0, .Z = std::numeric_limits<double>::infinity()})));
  }

  TEST_CASE("space vector: a unit vector, from an extrusion of 1e300 too")
  {
    CHECK(Unit({.X = 0, .Y = 0, .Z = 1e300}) == Up);
    CHECK(Unit({.X = 2, .Y = 3, .Z = 6})     == Vector{.X = 2.0 / 7, .Y = 3.0 / 7, .Z = 6.0 / 7});
    CHECK(Length(*Unit({.X = 1e-300, .Y = 2e-300, .Z = 3e-300}))
          == doctest::Approx(1.0).epsilon(1e-15));
    CHECK_FALSE(Unit({}));
    CHECK_FALSE(Unit({.X = std::numeric_limits<double>::infinity(), .Y = 0, .Z = 0}));
    CHECK_FALSE(Unit({.X = 0, .Y = std::numeric_limits<double>::quiet_NaN(), .Z = 1}));
    CHECK_FALSE(Unit({.X = 0, .Y = 0, .Z = std::numeric_limits<double>::infinity()}));
  }
} // namespace Geoxy::Space
