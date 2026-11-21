/***************************************************************************************************
 *  @file      tests/geoxy/box.cpp
 *  @brief     Tests of boxes.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/box.hpp"

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <limits>

namespace Geoxy
{
  namespace
  {
    constexpr Box TenSquare = {.Min = {.X = 0, .Y = 0}, .Max = {.X = 10, .Y = 10}};
  } // namespace

  static_assert(IsEmpty(Box{}));
  static_assert(!IsEmpty(TenSquare));
  static_assert(Grown(Box{}, {.X = 1, .Y = 2})
                == Box{.Min = {.X = 1, .Y = 2}, .Max = {.X = 1, .Y = 2}});
  static_assert(!IsEmpty(Grown(Box{}, {.X = 1, .Y = 2})));
  static_assert(Grown(TenSquare, {.X = -1, .Y = 20})
                == Box{.Min = {.X = -1, .Y = 0}, .Max = {.X = 10, .Y = 20}});
  static_assert(Joined(TenSquare, Box{}) == TenSquare);
  static_assert(Joined(Box{}, TenSquare) == TenSquare);
  static_assert(Joined(TenSquare, Box{.Min = {.X = 5, .Y = -5}, .Max = {.X = 20, .Y = 5}})
                == Box{.Min = {.X = 0, .Y = -5}, .Max = {.X = 20, .Y = 10}});
  static_assert(Contains(TenSquare, {.X = 10, .Y = 0}));
  static_assert(!Contains(TenSquare, {.X = 10.5, .Y = 0}));
  static_assert(!Contains(TenSquare, {.X = 5, .Y = -1}));
  static_assert(!Contains(Box{}, {.X = 0, .Y = 0}));
  static_assert(Overlaps(TenSquare, Box{.Min = {.X = 10, .Y = 10}, .Max = {.X = 20, .Y = 20}}));
  static_assert(!Overlaps(TenSquare, Box{.Min = {.X = 11, .Y = 0}, .Max = {.X = 20, .Y = 10}}));
  static_assert(!Overlaps(TenSquare, Box{.Min = {.X = 0, .Y = 11}, .Max = {.X = 10, .Y = 20}}));
  static_assert(!Overlaps(TenSquare, Box{}));
  static_assert(Size(TenSquare) == Vector{.X = 10, .Y = 10});

  TEST_CASE("box: the tests at run time")
  {
    const Box square = TenSquare;
    CHECK(IsEmpty(Box{}));
    CHECK_FALSE(IsEmpty(square));
    CHECK(IsEmpty(Box{.Min = {.X = 0, .Y = 2}, .Max = {.X = 1, .Y = 1}}));
    CHECK(Grown(Box{}, {.X = 1, .Y = 2}) == Box{.Min = {.X = 1, .Y = 2}, .Max = {.X = 1, .Y = 2}});
    CHECK(Joined(square, Box{})          == square);
    CHECK(Contains(square, {.X = 10, .Y = 0}));
    CHECK_FALSE(Contains(square, {.X = -1, .Y = 0}));
    CHECK_FALSE(Contains(square, {.X = 11, .Y = 0}));
    CHECK_FALSE(Contains(square, {.X = 5, .Y = -1}));
    CHECK_FALSE(Contains(square, {.X = 5, .Y = 11}));
    CHECK(Overlaps(square, Box{.Min = {.X = 10, .Y = 10}, .Max = {.X = 20, .Y = 20}}));
    CHECK_FALSE(Overlaps(square, Box{.Min = {.X = 11, .Y = 0}, .Max = {.X = 20, .Y = 10}}));
    CHECK_FALSE(Overlaps(square, Box{.Min = {.X = -20, .Y = 0}, .Max = {.X = -11, .Y = 10}}));
    CHECK_FALSE(Overlaps(square, Box{.Min = {.X = 0, .Y = 11}, .Max = {.X = 10, .Y = 20}}));
    CHECK_FALSE(Overlaps(square, Box{.Min = {.X = 0, .Y = -20}, .Max = {.X = 10, .Y = -11}}));
    CHECK(Size(square) == Vector{.X = 10, .Y = 10});
    CHECK_FALSE(square == Box{.Min = {.X = 0, .Y = 0}, .Max = {.X = 10, .Y = 11}});
    CHECK_FALSE(square == Box{.Min = {.X = 1, .Y = 0}, .Max = {.X = 10, .Y = 10}});
  }

  TEST_CASE("box: a box with a NaN corner is empty")
  {
    CHECK(IsEmpty(Box{.Min = {.X = std::numeric_limits<double>::quiet_NaN(), .Y = 0},
                      .Max = {.X = 1, .Y = 1}}));
    CHECK(IsEmpty(Box{.Min = {.X = 0, .Y = 2}, .Max = {.X = 1, .Y = 1}}));
  }

  TEST_CASE("box: the middle of a box at every scale")
  {
    CHECK(Center(TenSquare) == Point{.X = 5, .Y = 5});
    CHECK(Center({.Min = {.X = 5e6, .Y = 1e8}, .Max = {.X = 5e6 + 2, .Y = 1e8 + 4}})
          == Point{.X = 5e6 + 1, .Y = 1e8 + 2});
    constexpr double Huge = std::numeric_limits<double>::max();
    CHECK(Center({.Min = {.X = -Huge, .Y = -Huge}, .Max = {.X = Huge, .Y = Huge}})
          == Point{.X = 0, .Y = 0});
  }

  TEST_CASE("box: padded on every side, and past nothing to empty")
  {
    CHECK(Padded(TenSquare, 1)  == Box{.Min = {.X = -1, .Y = -1}, .Max = {.X = 11, .Y = 11}});
    CHECK(Padded(TenSquare, -5) == Box{.Min = {.X = 5, .Y = 5}, .Max = {.X = 5, .Y = 5}});
    CHECK(IsEmpty(Padded(TenSquare, -6)));
  }
} // namespace Geoxy
