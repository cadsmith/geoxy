/***************************************************************************************************
 *  @file      tests/bench/region.cpp
 *  @brief     Benchmarks of regions: the union and intersection of two smooth stars of 8 to 512
 *             sides and of two grids of rounded rectangles, a point's winding about a star, and the
 *             loop around a point picked in a plan, made once or picked again.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/region.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/arrangement.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/line.hpp"
#include "geoxy/loop.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"

#include <benchmark/benchmark.h>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The tolerance every case is met within.
     */
    constexpr Tolerance BenchFine = {.Length = 1e-9, .Radians = 1e-9};

    /**
     *  A smooth star: a polygon of some sides about a center, its radius waving with the angle, as
     *  the probe against other engines draws it.
     *
     *  @param  sides  How many sides.
     *  @param  center The center.
     *  @param  waves  How many times the radius waves in a turn.
     *  @return        The region.
     */
    auto Star(std::size_t sides, Point center, double waves) -> Region
    {
      Loop               loop = {};
      std::vector<Point> points;
      for (std::size_t index = 0; index < sides; ++index)
      {
        const double angle =
            2 * std::numbers::pi * static_cast<double>(index) / static_cast<double>(sides);
        const double radius = 1 + (0.25 * std::sin(waves * angle));
        points.push_back({.X = center.X + (radius * std::cos(angle)),
                          .Y = center.Y + (radius * std::sin(angle))});
      }
      for (std::size_t index = 0; index < sides; ++index)
      {
        loop.Edges.emplace_back(Line{.Start = points[index], .End = points[(index + 1) % sides]});
      }
      return {.Loops = {loop}};
    }

    /**
     *  A grid of rounded rectangles, 0.8 a side with corners of radius 0.2, a unit apart.
     *
     *  @param  count  How many along each side.
     *  @param  corner The least corner of the first.
     *  @return        The region.
     */
    auto Pills(std::size_t count, Point corner) -> Region
    {
      Region region = {};
      for (std::size_t row = 0; row < count; ++row)
      {
        for (std::size_t column = 0; column < count; ++column)
        {
          const double left   = corner.X + static_cast<double>(column);
          const double bottom = corner.Y + static_cast<double>(row);
          const double right  = left + 0.8;
          const double top    = bottom + 0.8;
          const double r      = 0.2;
          region.Loops.push_back(
              {.Edges = {
                   Line{.Start = {.X = left + r, .Y = bottom},
                        .End   = {.X = right - r, .Y = bottom}},
                   Arc{.Center = {.X = right - r, .Y = bottom + r},
                       .Radius = r,
                       .Start  = Degrees(-90),
                       .Sweep  = Degrees(90)},
                   Line{.Start = {.X = right, .Y = bottom + r}, .End = {.X = right, .Y = top - r}},
                   Arc{.Center = {.X = right - r, .Y = top - r},
                       .Radius = r,
                       .Start  = Degrees(0),
                       .Sweep  = Degrees(90)},
                   Line{.Start = {.X = right - r, .Y = top}, .End = {.X = left + r, .Y = top}},
                   Arc{.Center = {.X = left + r, .Y = top - r},
                       .Radius = r,
                       .Start  = Degrees(90),
                       .Sweep  = Degrees(90)},
                   Line{.Start = {.X = left, .Y = top - r}, .End = {.X = left, .Y = bottom + r}},
                   Arc{.Center = {.X = left + r, .Y = bottom + r},
                       .Radius = r,
                       .Start  = Degrees(180),
                       .Sweep  = Degrees(90)},
               }});
        }
      }
      return region;
    }

    /**
     *  A plan of rooms: a grid of walls, each wall a line across the whole plan.
     *
     *  @param  rooms How many rooms along each side.
     *  @return       The walls.
     */
    auto RegionPlan(std::size_t rooms) -> std::vector<Curve>
    {
      std::vector<Curve> walls;
      const auto         side = static_cast<double>(rooms) * 4;
      for (std::size_t index = 0; index <= rooms; ++index)
      {
        const auto line = static_cast<double>(index) * 4;
        walls.emplace_back(Line{.Start = {.X = line, .Y = 0}, .End = {.X = line, .Y = side}});
        walls.emplace_back(Line{.Start = {.X = 0, .Y = line}, .End = {.X = side, .Y = line}});
        for (std::size_t other = 0; other < rooms; ++other)
        {
          walls.emplace_back(
              Circle({.X = line + 2, .Y = (static_cast<double>(other) * 4) + 2}, 0.5));
        }
      }
      return walls;
    }

    /**
     *  The union of two smooth stars of some sides.
     *
     *  @param state The benchmark's state; its range, the sides.
     */
    void RegionStars(benchmark::State &state)
    {
      const auto   sides  = static_cast<std::size_t>(state.range(0));
      const Region first  = Star(sides, {}, 5);
      const Region second = Star(sides, {.X = 0.3, .Y = 0.2}, 7);
      for (auto pass : state)
      {
        benchmark::DoNotOptimize(Union(first, second, BenchFine));
      }
    }

    /**
     *  The intersection of two smooth stars of some sides.
     *
     *  @param state The benchmark's state; its range, the sides.
     */
    void RegionStarsShared(benchmark::State &state)
    {
      const auto   sides  = static_cast<std::size_t>(state.range(0));
      const Region first  = Star(sides, {}, 5);
      const Region second = Star(sides, {.X = 0.3, .Y = 0.2}, 7);
      for (auto pass : state)
      {
        benchmark::DoNotOptimize(Intersection(first, second, BenchFine));
      }
    }

    /**
     *  The union of two grids of rounded rectangles, the second half a unit across from the first.
     *
     *  @param state The benchmark's state; its range, the rectangles along a side.
     */
    void RegionPills(benchmark::State &state)
    {
      const auto   count  = static_cast<std::size_t>(state.range(0));
      const Region first  = Pills(count, {});
      const Region second = Pills(count, {.X = 0.5, .Y = 0.5});
      for (auto pass : state)
      {
        benchmark::DoNotOptimize(Union(first, second, BenchFine));
      }
    }

    /**
     *  How many times a star of 512 sides winds about a point inside it.
     *
     *  @param state The benchmark's state.
     */
    void LoopWinding(benchmark::State &state)
    {
      const Loop loop = Star(512, {}, 5).Loops.front();
      for (auto pass : state)
      {
        benchmark::DoNotOptimize(Winding(loop, {.X = 0.1, .Y = 0.05}, BenchFine));
      }
    }

    /**
     *  A plan's walls and columns arranged.
     *
     *  @param state The benchmark's state; its range, the rooms along a side.
     */
    void RegionArranged(benchmark::State &state)
    {
      const std::vector<Curve> walls = RegionPlan(static_cast<std::size_t>(state.range(0)));
      for (auto pass : state)
      {
        benchmark::DoNotOptimize(Arrangement{walls, {.Within = BenchFine, .Gap = 0}});
      }
    }

    /**
     *  The loop around a point picked in a plan arranged once.
     *
     *  @param state The benchmark's state; its range, the rooms along a side.
     */
    void RegionPicked(benchmark::State &state)
    {
      const Arrangement plan{RegionPlan(static_cast<std::size_t>(state.range(0))),
                             {.Within = BenchFine, .Gap = 0}};
      for (auto pass : state)
      {
        benchmark::DoNotOptimize(Boundary(plan, {.X = 5, .Y = 7}));
      }
    }
  } // namespace

  BENCHMARK(RegionStars)->Arg(8)->Arg(64)->Arg(512);
  BENCHMARK(RegionStarsShared)->Arg(512);
  BENCHMARK(RegionPills)->Arg(4)->Arg(16);
  BENCHMARK(LoopWinding);
  BENCHMARK(RegionArranged)->Arg(8)->Arg(32);
  BENCHMARK(RegionPicked)->Arg(8)->Arg(32);
} // namespace Geoxy
