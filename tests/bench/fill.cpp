/***************************************************************************************************
 *  @file      tests/bench/fill.cpp
 *  @brief     Benchmarks of fill triangles, by earcut's measure: time per vertex, star polygons of
 *             15 to 10,000 vertices with holes; and of splines fit through points.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/fit.hpp"
#include "geoxy/point.hpp"
#include "geoxy/triangulate.hpp"

#include <array>
#include <benchmark/benchmark.h>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <random>
#include <vector>

/**
 *  earcut.hpp's README gives 119,680 map polygons, 1.9 M vertices, in 276 ms on an M1 Pro: about
 *  145 ns a vertex. The Vertices counter here is a rate, so its inverse is the time per vertex.
 */

namespace Geoxy::Bench
{
  namespace
  {
    /**
     *  A star of points about a center, its radius wavering.
     *
     *  @param  center The center.
     *  @param  count  How many points.
     *  @param  radius Its radius, and the seed.
     *  @return        The points, counterclockwise.
     */
    auto Star(Point center, std::size_t count, std::array<double, 2> radius) -> std::vector<Point>
    {
      std::mt19937_64                        draw(static_cast<std::uint64_t>(radius[1]));
      std::uniform_real_distribution<double> waver(0.7, 1.0);
      std::vector<Point>                     points;
      for (std::size_t index = 0; index < count; ++index)
      {
        const double angle =
            2.0 * std::numbers::pi * static_cast<double>(index) / static_cast<double>(count);
        const double reach = radius[0] * waver(draw);
        points.push_back(
            {.X = center.X + (reach * std::cos(angle)), .Y = center.Y + (reach * std::sin(angle))});
      }
      return points;
    }

    /**
     *  A region of a star outline of a count of vertices, and a hole for every hundred.
     *
     *  @param  count The outline's vertices.
     *  @return       The loops.
     */
    auto Region(std::size_t count) -> std::vector<std::vector<Point>>
    {
      std::vector<std::vector<Point>> loops = {Star({}, count, {100, 1})};
      const std::size_t               holes = count / 100;
      for (std::size_t hole = 0; hole < holes; ++hole)
      {
        const double angle =
            2.0 * std::numbers::pi * static_cast<double>(hole) / static_cast<double>(holes);
        loops.push_back(
            Star({.X = 40 * std::cos(angle), .Y = 40 * std::sin(angle)}, 16,
                 {60.0 / static_cast<double>(holes + 2), 2.0 + static_cast<double>(hole)}));
      }
      return loops;
    }

    void ByFill(benchmark::State &state, Fill fill)
    {
      const std::vector<std::vector<Point>> loops =
          Region(static_cast<std::size_t>(state.range(0)));
      std::size_t vertices = 0;
      for (const std::vector<Point> &loop : loops)
      {
        vertices += loop.size();
      }
      for ([[maybe_unused]] auto iteration : state)
      {
        Mesh mesh = Triangulate(loops, fill);
        benchmark::DoNotOptimize(mesh);
      }
      state.counters["Vertices"] = benchmark::Counter(
          static_cast<double>(vertices * static_cast<std::size_t>(state.iterations())),
          benchmark::Counter::kIsRate);
    }

    void Ears(benchmark::State &state)
    {
      ByFill(state, Fill::Ears);
    }

    void Delaunay(benchmark::State &state)
    {
      ByFill(state, Fill::Delaunay);
    }

    void Fitted(benchmark::State &state)
    {
      const std::vector<Point> points =
          Star({}, static_cast<std::size_t>(state.range(0)), {100, 3});
      for ([[maybe_unused]] auto iteration : state)
      {
        auto spline = Fit(points, {.Closed = true});
        benchmark::DoNotOptimize(spline);
      }
    }
  } // namespace

  BENCHMARK(Ears)->Arg(15)->Arg(100)->Arg(1000)->Arg(10'000);
  BENCHMARK(Delaunay)->Arg(15)->Arg(100)->Arg(1000)->Arg(10'000);
  BENCHMARK(Fitted)->Arg(8)->Arg(100)->Arg(10'000);
} // namespace Geoxy::Bench
