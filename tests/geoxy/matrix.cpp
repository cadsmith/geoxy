/***************************************************************************************************
 *  @file      tests/geoxy/matrix.cpp
 *  @brief     Tests of affine maps.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/matrix.hpp"

#include "draws.hpp"
#include "real.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/box.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <optional>

namespace Geoxy
{
  static_assert(Matrix{} == Matrix{.X = {.X = 1, .Y = 0}, .Y = {.X = 0, .Y = 1}, .Origin = {}});
  static_assert(Translation({.X = 2, .Y = 3}).Origin == Point{.X = 2, .Y = 3});
  static_assert(Scaling(2, 3) == Matrix{.X = {.X = 2, .Y = 0}, .Y = {.X = 0, .Y = 3}});
  static_assert(Scaling(2)    == Scaling(2, 2));

  namespace
  {
    auto Near(Point actual, Point expected) -> bool
    {
      return actual.X == doctest::Approx(expected.X).epsilon(1e-14)
          && actual.Y == doctest::Approx(expected.Y).epsilon(1e-14);
    }
  } // namespace

  TEST_CASE("matrix: the simple maps at run time, and equality number by number")
  {
    CHECK(Translation({.X = 2, .Y = 3}) == Matrix{.Origin = {.X = 2, .Y = 3}});
    CHECK(Scaling(2)     == Scaling(2, 2));
    CHECK_FALSE(Matrix{} == Scaling(2, 1));
    CHECK_FALSE(Matrix{} == Scaling(1, 2));
    CHECK_FALSE(Matrix{} == Translation({.X = 1, .Y = 0}));
  }

  TEST_CASE("matrix: points and displacements mapped")
  {
    const Matrix map = {
        .X = {.X = 2, .Y = 1}, .Y = {.X = -1, .Y = 3}, .Origin = {.X = 10, .Y = 20}};
    CHECK(map * Point{.X = 1, .Y = 2}  == Point{.X = 10, .Y = 27});
    CHECK(map * Vector{.X = 1, .Y = 2} == Vector{.X = 0, .Y = 7});
    CHECK(Translation({.X = 5e6, .Y = 1}) * Vector{.X = 1, .Y = 2} == Vector{.X = 1, .Y = 2});
  }

  TEST_CASE("matrix: a turn of a quarter is exact, about the origin and about a point")
  {
    CHECK(Rotation(Degrees(90)) * Point{.X = 1, .Y = 0}  == Point{.X = 0, .Y = 1});
    CHECK(Rotation(Degrees(180)) * Point{.X = 1, .Y = 2} == Point{.X = -1, .Y = -2});
    const Point about = {.X = 5e6, .Y = 5e6};
    CHECK(Rotation(Degrees(90), about) * Point{.X = 5e6 + 1, .Y = 5e6}
          == Point{.X = 5e6, .Y = 5e6 + 1});
    CHECK(Near(Rotation(Degrees(30)) * Point{.X = 2, .Y = 0}, {.X = std::sqrt(3.0), .Y = 1}));
  }

  TEST_CASE("matrix: a mirror across a line")
  {
    CHECK(Mirror({.X = 0, .Y = 0}, {.X = 1, .Y = 0})
          == Matrix{.X = {.X = 1, .Y = 0}, .Y = {.X = 0, .Y = -1}});
    CHECK(*Mirror({.X = 0, .Y = 5e6}, {.X = 1, .Y = 5e6}) * Point{.X = 2, .Y = 5e6 + 3}
          == Point{.X = 2, .Y = 5e6 - 3});
    CHECK(Near(*Mirror({.X = 0, .Y = 0}, {.X = 1, .Y = 1}) * Point{.X = 1, .Y = 0},
               {.X = 0, .Y = 1}));
    CHECK(Determinant(*Mirror({.X = 1, .Y = 2}, {.X = 3, .Y = 7})) == doctest::Approx(-1.0));
    CHECK_FALSE(Mirror({.X = 1, .Y = 1}, {.X = 1, .Y = 1}));
  }

  TEST_CASE("matrix: the right one first")
  {
    const Matrix turn = Rotation(Degrees(90));
    const Matrix move = Translation({.X = 10, .Y = 0});
    CHECK((move * turn) * Point{.X = 1, .Y = 0} == Point{.X = 10, .Y = 1});
    CHECK((turn * move) * Point{.X = 1, .Y = 0} == Point{.X = 0, .Y = 11});
    CHECK((move * turn) * Point{.X = 3, .Y = 4} == move * (turn * Point{.X = 3, .Y = 4}));
  }

  TEST_CASE("matrix: a box mapped is the box of its corners mapped")
  {
    const Box box = {.Min = {.X = 0, .Y = 0}, .Max = {.X = 10, .Y = 10}};
    CHECK(Translation({.X = 1, .Y = 2}) * box
          == Box{.Min = {.X = 1, .Y = 2}, .Max = {.X = 11, .Y = 12}});
    const Box turned = Rotation(Degrees(45)) * box;
    CHECK(turned.Max.Y == doctest::Approx(10 * std::sqrt(2.0)));
    CHECK(turned.Min.X == doctest::Approx(-5 * std::sqrt(2.0)));
    CHECK(IsEmpty(Rotation(Degrees(45)) * Box{}));
  }

  TEST_CASE("matrix: the determinant scales areas")
  {
    CHECK(Determinant({})             == 1.0);
    CHECK(Determinant(Scaling(2, 3))  == 6.0);
    CHECK(Determinant(Scaling(1, -1)) == -1.0);
  }

  TEST_CASE("matrix: an inverse undoes the map")
  {
    Draws draws{0x696E76U};
    for (int index = 0; index < 1000; ++index)
    {
      const Matrix moved  = Translation({.X = draws.Near(5e6, 1e3), .Y = draws.Near(-3e6, 1e3)});
      const Matrix turned = Rotation(Radians(draws.Near(0, 4)));
      const double alongX = draws.Near(2, 1);
      const double alongY = draws.Near(-2, 1);
      const Matrix map    = moved * turned * Scaling(alongX, alongY);
      const std::optional<Matrix> inverse = Inverse(map);
      REQUIRE(inverse);
      const Point point = {.X = draws.Near(5e6, 1e3), .Y = draws.Near(-3e6, 1e3)};
      const Point back  = *inverse * (map * point);
      CHECK(Distance(back, point) <= 1e-7);
    }
    CHECK(Inverse(Scaling(2, 4))        == Scaling(0.5, 0.25));
    CHECK(Inverse(Scaling(1e-300))->X.X == doctest::Approx(1e300).epsilon(1e-15));
  }

  TEST_CASE("matrix: no inverse for a flat map, a nearly flat one or one past a double")
  {
    CHECK_FALSE(Inverse(Scaling(1, 0)));
    CHECK_FALSE(Inverse(Scaling(0)));
    /**
     *  Axes a ten million millionth from parallel have an inverse all the same, and it undoes the
     *  map; axes whose determinant doubles would round to 0 have theirs too.
     */
    const Matrix near = {.X = {.X = 1, .Y = 0}, .Y = {.X = 1, .Y = 1e-13}};
    REQUIRE(Inverse(near));
    CHECK(Inverse(near)->Y.Y == doctest::Approx(1e13).epsilon(1e-15));
    CHECK(Inverse(near)->Y.X == doctest::Approx(-1e13).epsilon(1e-15));
    CHECK(Inverse(near)->X   == Vector{.X = 1, .Y = 0});
    const Matrix nearer = {.X = {.X = 1 + 0x1p-30, .Y = 1}, .Y = {.X = 1, .Y = 1 - 0x1p-30}};
    REQUIRE(Inverse(nearer));
    CHECK(Inverse(nearer)->X.X == doctest::Approx(-(1 - 0x1p-30) * 0x1p60).epsilon(1e-15));
    CHECK_FALSE(Inverse({.X = {.X = 1, .Y = 2}, .Y = {.X = 2, .Y = 4}}));
    CHECK_FALSE(Inverse({.X = {.X = Infinite, .Y = 0}, .Y = {.X = 0, .Y = 1}}));
    CHECK_FALSE(Inverse({.X = {.X = 0x1.0p-1030, .Y = 0}, .Y = {.X = 0, .Y = 1}}));
    CHECK_FALSE(Inverse({.X = {.X = 1, .Y = 0}, .Y = {.X = 0, .Y = 0x1.0p-1030}}));
    CHECK_FALSE(Inverse(
        {.X = {.X = 1e-300, .Y = 0}, .Y = {.X = 0, .Y = 1e-300}, .Origin = {.X = 1e300, .Y = 0}}));
  }

  TEST_CASE("matrix: the stretch is the largest singular value, not the longest axis")
  {
    CHECK(Stretch({})             == 1.0);
    CHECK(Stretch(Scaling(2, -3)) == 3.0);
    CHECK(Stretch(Rotation(Radians(0.3))) == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(Stretch({.X = {.X = 1, .Y = 0}, .Y = {.X = 1, .Y = 1}})
          == doctest::Approx((1 + std::sqrt(5.0)) / 2).epsilon(1e-15));
    CHECK(Stretch({.X = {.X = 1, .Y = 0}, .Y = {.X = 0, .Y = 0}}) == 1.0);
    CHECK(Stretch(Scaling(1e300, 2e300)) == doctest::Approx(2e300).epsilon(1e-15));
    CHECK(Stretch(Scaling(0))            == 0.0);
    CHECK(std::isinf(Stretch(Scaling(Infinite))));
  }
} // namespace Geoxy
