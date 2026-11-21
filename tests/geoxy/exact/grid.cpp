/***************************************************************************************************
 *  @file      tests/geoxy/exact/grid.cpp
 *  @brief     Tests of the circle test's determinant in integers, against Shewchuk's own.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/grid.hpp"

#include "../draws.hpp"
#include "exact/triangle.hpp"
#include "exact/two.hpp"
#include "original.hpp"

#include "geoxy/point.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>

namespace Geoxy::Exact
{
  namespace
  {
    auto GridOriginal(const Triangle &circle, Point point) -> int
    {
      std::array<double, 2> a = {circle.First.X, circle.First.Y};
      std::array<double, 2> b = {circle.Second.X, circle.Second.Y};
      std::array<double, 2> c = {circle.Third.X, circle.Third.Y};
      std::array<double, 2> d = {point.X, point.Y};
      return Sign(incircle(a.data(), b.data(), c.data(), d.data()));
    }

    /**
     *  OnGrid decides, and gives the original's sign.
     */
    auto GridAgrees(const Triangle &circle, Point point) -> bool
    {
      const std::optional<int> sign = OnGrid(circle, point);
      return sign.has_value() && *sign == GridOriginal(circle, point);
    }

    /**
     *  OrientGrid decides, and gives the original's sign.
     */
    auto Turns(const Triangle &corners) -> bool
    {
      std::array<double, 2>    a    = {corners.First.X, corners.First.Y};
      std::array<double, 2>    b    = {corners.Second.X, corners.Second.Y};
      std::array<double, 2>    c    = {corners.Third.X, corners.Third.Y};
      const std::optional<int> sign = OrientGrid(corners);
      return sign.has_value() && *sign == Sign(orient2d(a.data(), b.data(), c.data()));
    }

    struct GridFixture
    {
      GridFixture() { exactinit(); }
    };
  } // namespace

  TEST_CASE("grid: a product from halves is the product")
  {
    Draws draws{0x6861'6C66U};

    constexpr std::uint64_t Top = std::numeric_limits<std::uint64_t>::max();
    for (const std::array<std::uint64_t, 2> edge :
         {std::array<std::uint64_t, 2>{0, 0}, {Top, Top}, {Top, 1}, {0xFFFF'FFFFU, 0x1'0000'0000U}})
    {
      const Product halves = MultiplyHalves(edge[0], edge[1]);
      const Product whole  = Multiply(edge[0], edge[1]);
      CHECK(halves.High == whole.High);
      CHECK(halves.Low  == whole.Low);
    }
    int same = 0;
    for (int index = 0; index < 100'000; ++index)
    {
      const std::uint64_t first  = draws.Next();
      const std::uint64_t second = draws.Next() >> (draws.Next() % 64U);
      const Product       halves = MultiplyHalves(first, second);
      const Product       whole  = Multiply(first, second);
      same += static_cast<int>(halves.High == whole.High && halves.Low == whole.Low);
    }
    CHECK(same == 100'000);
  }

  TEST_CASE_FIXTURE(GridFixture, "grid: four zero places are on their circle")
  {
    CHECK(OnGrid({}, {}) == std::optional<int>{0});
  }

  TEST_CASE_FIXTURE(GridFixture,
                    "grid: a place not normal, or places spanning 61 bits, are left out")
  {
    const Triangle circle = {
        .First = {.X = 0, .Y = 1}, .Second = {.X = -1, .Y = -1}, .Third = {.X = 1, .Y = -1}};
    CHECK_FALSE(OnGrid(circle, {.X = std::numeric_limits<double>::denorm_min(), .Y = 0}));
    CHECK_FALSE(OnGrid(circle, {.X = std::numeric_limits<double>::infinity(), .Y = 0}));
    CHECK_FALSE(OnGrid(circle, {.X = std::numeric_limits<double>::quiet_NaN(), .Y = 0}));
    CHECK_FALSE(OnGrid(circle, {.X = 0x1.0p-61, .Y = 0}));
    CHECK(GridAgrees(circle, {.X = 0x1.0p-60, .Y = 0}));
  }

  /**
   *  The corners of a rectangle lie on one circle, exactly: snapped to the millimetre at survey
   *  coordinates, every determinant is zero, and a millimetre off it, not.
   */
  TEST_CASE_FIXTURE(GridFixture, "grid: a rectangle's corners in millimetres at survey scale")
  {
    constexpr int Count = 20'000;

    Draws draws{0x7265'6374U};

    const auto millimetres = [](double value) -> double
    {
      return std::round(value * 1000) / 1000;
    };

    int zero   = 0;
    int agreed = 0;
    for (int index = 0; index < Count; ++index)
    {
      const double   left   = millimetres(draws.Near(5e5, 1e3));
      const double   bottom = millimetres(draws.Near(4.5e6, 1e3));
      const double   right  = millimetres(left + draws.Near(50, 49));
      const double   top    = millimetres(bottom + draws.Near(50, 49));
      const Triangle circle = {.First  = {.X = left, .Y = bottom},
                               .Second = {.X = right, .Y = bottom},
                               .Third  = {.X = right, .Y = top}};
      zero   += static_cast<int>(OnGrid(circle, {.X = left, .Y = top}) == std::optional<int>{0});
      agreed += static_cast<int>(GridAgrees(circle, {.X = left, .Y = top}));
      agreed += static_cast<int>(GridAgrees(circle, {.X = left, .Y = millimetres(top + 0.001)}));
    }
    CHECK(zero   == Count);
    CHECK(agreed == 2 * Count);
  }

  /**
   *  Places of whole numbers up to 2^60 in magnitude, either sign, spread and on a lattice circle:
   *  products near 2^250, carries through every word, both signs of every term.
   */
  TEST_CASE_FIXTURE(GridFixture, "grid: whole numbers up to 2^60 agree with the original")
  {
    constexpr int Count = 50'000;

    Draws draws{0x6269'6773U};

    const auto whole = [&draws](unsigned bits) -> double
    {
      const auto magnitude = static_cast<double>(draws.Next() >> (64U - bits));
      return (draws.Next() & 1U) == 0 ? magnitude : -magnitude;
    };
    const auto place = [&whole, &draws] -> Point
    {
      const auto bits = static_cast<unsigned>(1 + (draws.Next() % 60U));
      return {.X = whole(bits), .Y = whole(bits)};
    };

    int agreed = 0;
    for (int index = 0; index < Count; ++index)
    {
      const Triangle circle  = {.First = place(), .Second = place(), .Third = place()};
      agreed                += static_cast<int>(GridAgrees(circle, place()));
    }
    CHECK(agreed == Count);
  }

  /**
   *  Places on a circle by the eight images of a Pythagorean triple scaled near 2^57, about a
   *  centre whose low bits set the grid: zero, and a unit off it, not.
   */
  TEST_CASE_FIXTURE(GridFixture, "grid: a lattice circle near the top of the grid is exact")
  {
    constexpr double Scale = 0x1.0p57;

    const Point    centre = {.X = 0x1.0p-2, .Y = -0x1.8p-2};
    const Triangle circle = {.First  = {.X = centre.X + (3 * Scale), .Y = centre.Y + (4 * Scale)},
                             .Second = {.X = centre.X - (4 * Scale), .Y = centre.Y + (3 * Scale)},
                             .Third  = {.X = centre.X - (3 * Scale), .Y = centre.Y - (4 * Scale)}};
    CHECK(OnGrid(circle, {.X = 5 * Scale, .Y = centre.Y}) == std::optional<int>{});
    const Triangle whole = {.First  = {.X = 3 * Scale, .Y = 4 * Scale},
                            .Second = {.X = -4 * Scale, .Y = 3 * Scale},
                            .Third  = {.X = -3 * Scale, .Y = -4 * Scale}};
    CHECK(OnGrid(whole, {.X = 5 * Scale, .Y = 0})          == std::optional<int>{0});
    CHECK(OnGrid(whole, {.X = 5 * Scale, .Y = 1})          == std::optional<int>{-1});
    CHECK(OnGrid(whole, {.X = (5 * Scale) - 4096, .Y = 0}) == std::optional<int>{1});
    CHECK(GridAgrees(whole, {.X = 5 * Scale, .Y = 1}));
  }

  TEST_CASE_FIXTURE(GridFixture, "grid: three zero corners are collinear")
  {
    CHECK(OrientGrid({}) == std::optional<int>{0});
  }

  TEST_CASE_FIXTURE(GridFixture,
                    "grid: a corner not normal, or corners spanning 125 bits, are left out of two "
                    "words")
  {
    const auto with = [](double x) -> Triangle
    {
      return {.First = {.X = 1, .Y = 0}, .Second = {.X = 0, .Y = 1}, .Third = {.X = x, .Y = 0}};
    };
    CHECK_FALSE(OrientGrid(with(std::numeric_limits<double>::denorm_min())));
    CHECK_FALSE(OrientGrid(with(std::numeric_limits<double>::infinity())));
    CHECK_FALSE(OrientGrid(with(std::numeric_limits<double>::quiet_NaN())));
    CHECK_FALSE(OrientGrid(with(0x1.0p-125)));
    CHECK(OrientGrid(with(0x1.0p-124)) == std::optional<int>{1});
  }

  /**
   *  Corners on a diagonal from 2^60 down to 2^-60: their units shift past a word, and they turn
   *  exactly as the original says, collinear and nudged off either way.
   */
  TEST_CASE_FIXTURE(GridFixture, "grid: corners 120 bits apart turn as the original says")
  {
    const Point far  = {.X = 0x1.0p60, .Y = 0x1.0p60};
    const Point near = {.X = 0x1.0p-60, .Y = 0x1.0p-60};
    CHECK(OrientGrid({.First = far, .Second = near, .Third = {}}) == std::optional<int>{0});
    CHECK(OrientGrid({.First = {.X = -far.X, .Y = -far.Y}, .Second = near, .Third = {}})
          == std::optional<int>{0});
    CHECK(Turns({.First = far, .Second = {.X = near.X, .Y = 0x1.8p-60}, .Third = {}}));
    CHECK(Turns({.First = far, .Second = {.X = 0x1.8p-60, .Y = near.Y}, .Third = {}}));
    CHECK(Turns({.First = {.X = 0x1.0p60, .Y = -0x1.0p-60}, .Second = near, .Third = far}));
  }

  TEST_CASE_FIXTURE(GridFixture, "grid: places spanning up to 124 bits turn as the original says")
  {
    constexpr int Count = 50'000;

    Draws draws{0x7475'726EU};

    const auto place = [&draws](int top) -> double
    {
      const auto   bits      = static_cast<unsigned>(1 + (draws.Next() % 53U));
      const auto   magnitude = static_cast<double>(draws.Next() >> (64U - bits));
      const int    shift     = top - static_cast<int>(draws.Next() % 72U);
      const double value     = std::ldexp(magnitude, shift - static_cast<int>(bits));
      return (draws.Next() & 1U) == 0 ? value : -value;
    };

    int agreed = 0;
    for (int index = 0; index < Count; ++index)
    {
      const int      top      = static_cast<int>(draws.Next() % 600U) - 300;
      const Triangle corners  = {.First  = {.X = place(top), .Y = place(top)},
                                 .Second = {.X = place(top), .Y = place(top)},
                                 .Third  = {.X = place(top), .Y = place(top)}};
      agreed                 += static_cast<int>(Turns(corners));
    }
    CHECK(agreed == Count);
  }
} // namespace Geoxy::Exact
