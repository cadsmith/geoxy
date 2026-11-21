/***************************************************************************************************
 *  @file      tests/bench/exact/batch.cpp
 *  @brief     Benchmarks of the orientation test of many points against one line.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/batch.hpp"

#include "../../geoxy/draws.hpp"
#include "../../geoxy/exact/original.hpp"
#include "exact/orient.hpp"
#include "exact/two.hpp"

#include "geoxy/point.hpp"

#include <array>
#include <benchmark/benchmark.h>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace Geoxy::Exact
{
  namespace
  {
    constexpr std::size_t BatchCount = 65'536;

    constexpr Line Through = {.First = {.X = -0.3, .Y = 0.2}, .Second = {.X = 0.7, .Y = -0.4}};

    auto Points() -> std::vector<Point>
    {
      Draws draws{21};

      std::vector<Point> points;
      points.reserve(BatchCount);
      for (std::size_t index = 0; index < BatchCount; ++index)
      {
        points.push_back({.X = draws.Near(0, 1), .Y = draws.Near(0, 1)});
      }
      return points;
    }

    void Finish(benchmark::State &state)
    {
      state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(BatchCount));
    }

    /**
     *  Shewchuk's orient2d in a loop, one call a point.
     */
    void LineOriginal(benchmark::State &state)
    {
      exactinit();
      const std::vector<Point> points = Points();
      std::vector<std::int8_t> signs(BatchCount);
      std::array<double, 2>    first  = {Through.First.X, Through.First.Y};
      std::array<double, 2>    second = {Through.Second.X, Through.Second.Y};
      for ([[maybe_unused]] auto iteration : state)
      {
        for (std::size_t index = 0; index < BatchCount; ++index)
        {
          std::array<double, 2> third       = {points[index].X, points[index].Y};
          const double          determinant = orient2d(first.data(), second.data(), third.data());
          signs[index]                      = static_cast<std::int8_t>(Sign(determinant));
        }
        benchmark::ClobberMemory();
      }
      Finish(state);
    }

    /**
     *  Orient in a loop, its first stage inline.
     */
    void LineOrient(benchmark::State &state)
    {
      const std::vector<Point> points = Points();
      std::vector<std::int8_t> signs(BatchCount);
      for ([[maybe_unused]] auto iteration : state)
      {
        for (std::size_t index = 0; index < BatchCount; ++index)
        {
          signs[index] =
              static_cast<std::int8_t>(Orient(Through.First, Through.Second, points[index]));
        }
        benchmark::ClobberMemory();
      }
      Finish(state);
    }

    /**
     *  OrientSigns, its first stage over blocks in the vector unit.
     */
    void LineOrientSigns(benchmark::State &state)
    {
      const std::vector<Point> points = Points();
      std::vector<std::int8_t> signs(BatchCount);
      for ([[maybe_unused]] auto iteration : state)
      {
        OrientSigns(Through, points, signs);
        benchmark::ClobberMemory();
      }
      Finish(state);
    }
  } // namespace

  BENCHMARK(LineOriginal);
  BENCHMARK(LineOrient);
  BENCHMARK(LineOrientSigns);
} // namespace Geoxy::Exact
