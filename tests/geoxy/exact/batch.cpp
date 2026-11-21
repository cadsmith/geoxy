/***************************************************************************************************
 *  @file      tests/geoxy/exact/batch.cpp
 *  @brief     Tests of the orientation test of many points against one line.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/batch.hpp"

#include "../draws.hpp"
#include "exact/orient.hpp"

#include "geoxy/point.hpp"

#include <doctest/doctest.h>

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace Geoxy::Exact
{
  namespace
  {
    /**
     *  How many of the points get a sign other than Orient's.
     */
    auto Wrong(const Line &line, std::span<const Point> points) -> std::size_t
    {
      std::vector<std::int8_t> signs(points.size(), 7);
      OrientSigns(line, points, signs);
      std::size_t wrong = 0;
      for (std::size_t index = 0; index < points.size(); ++index)
      {
        wrong += static_cast<std::size_t>(signs[index]
                                          != Exact::Orient(line.First, line.Second, points[index]));
      }
      return wrong;
    }

    /**
     *  Kettner et al.'s grid near the diagonal through (12, 12) and (24, 24), every point of which
     *  the first stage leaves open, each followed by a point far off the line that it decides: so
     *  every block mixes the two.
     */
    auto Grid() -> std::vector<Point>
    {
      std::vector<Point> points;
      for (int i = 0; i < 64; ++i)
      {
        for (int j = 0; j < 64; ++j)
        {
          points.push_back({.X = 0.5 + (i * 0x1.0p-53), .Y = 0.5 + (j * 0x1.0p-53)});
          points.push_back({.X = 0.5 + j, .Y = 0.5 - i - 1});
        }
      }
      return points;
    }

    constexpr Line Diagonal = {.First = {.X = 12, .Y = 12}, .Second = {.X = 24, .Y = 24}};
  } // namespace

  TEST_CASE("batch: every point near a line gets Orient's sign")
  {
    CHECK(Wrong(Diagonal, Grid()) == 0);
  }

  TEST_CASE("batch: random points against random lines get Orient's sign")
  {
    Draws draws{0x6261'7463U};

    std::size_t wrong = 0;
    for (int line = 0; line < 100; ++line)
    {
      const Line         random = {.First  = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)},
                                   .Second = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)}};
      std::vector<Point> points;
      points.reserve(600);
      for (int index = 0; index < 600; ++index)
      {
        points.push_back({.X = draws.Near(0, 1), .Y = draws.Near(0, 1)});
      }
      wrong += Wrong(random, points);
    }
    CHECK(wrong == 0);
  }

  /**
   *  Lengths around the block of 256 and none at all, from offsets that break alignment: every
   *  block and its remainder, and nothing written past the signs asked for.
   */
  TEST_CASE("batch: lengths across the block, unaligned, writing nothing past the end")
  {
    const std::vector<Point> grid = Grid();
    for (const std::size_t length : {0UZ, 1UZ, 255UZ, 256UZ, 257UZ, 513UZ})
    {
      for (const std::size_t offset : {0UZ, 3UZ})
      {
        const std::span<const Point> points = std::span{grid}.subspan(offset * 97, length);
        std::vector<std::int8_t>     signs(offset + length + 1, 7);
        OrientSigns(Diagonal, points, std::span{signs}.subspan(offset, length));

        std::size_t wrong = 0;
        for (std::size_t index = 0; index < length; ++index)
        {
          wrong += static_cast<std::size_t>(
              signs[offset + index]
              != Exact::Orient(Diagonal.First, Diagonal.Second, points[index]));
        }
        std::size_t untouched = 0;
        for (std::size_t index = 0; index < offset; ++index)
        {
          untouched += static_cast<std::size_t>(signs[index] == 7);
        }
        CHECK(wrong        == 0);
        CHECK(untouched    == offset);
        CHECK(signs.back() == 7);
      }
    }
  }
} // namespace Geoxy::Exact
