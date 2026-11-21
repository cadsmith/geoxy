/***************************************************************************************************
 *  @file      source/exact/batch.cpp
 *  @brief     The orientation test of many points against one line, a block at a time.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/batch.hpp"

#include "exact/orient.hpp"
#include "exact/triangle.hpp"
#include "exact/two.hpp"

#include "geoxy/point.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

namespace Geoxy::Exact
{
  namespace
  {
    /**
     *  How many points a block holds: enough that the loop runs long in the vector unit, few enough
     *  that its flags stay in the first level cache.
     */
    constexpr std::size_t Block = 256;

    /**
     *  The corners a point makes with the line.
     *
     *  @param  line  The line.
     *  @param  point The point.
     *  @return       The line's two points, then the point.
     */
    auto CornersOf(const Line &line, Point point) -> Triangle
    {
      return {.First = line.First, .Second = line.Second, .Third = point};
    }
  } // namespace

  void OrientSigns(const Line &line, std::span<const Point> points, std::span<std::int8_t> signs)
  {
    std::array<std::uint8_t, Block> flags = {};
    const std::span<std::uint8_t>   open  = flags;
    for (std::size_t start = 0; start < points.size(); start += Block)
    {
      const std::span<const Point> block =
          points.subspan(start, std::min(Block, points.size() - start));
      const std::span<std::int8_t> out = signs.subspan(start, block.size());

      std::uint8_t anyOpen = 0;
      for (std::size_t index = 0; index < block.size(); ++index)
      {
        const Approximation approximation = Approximate(CornersOf(line, block[index]));
        out[index]   = static_cast<std::int8_t>(Sign(approximation.Determinant));
        open[index]  = static_cast<std::uint8_t>(!Decided(approximation));
        anyOpen     |= open[index];
      }
      if (anyOpen == 0) [[likely]]
      {
        continue;
      }
      for (std::size_t index = 0; index < block.size(); ++index)
      {
        if (open[index] != 0)
        {
          out[index] = static_cast<std::int8_t>(OrientAdapt(CornersOf(line, block[index])));
        }
      }
    }
  }
} // namespace Geoxy::Exact
