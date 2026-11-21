/***************************************************************************************************
 *  @file      tests/bench/offset.cpp
 *  @brief     Benchmarks of Offset, on profiles of lines and arcs and on curves that are neither
 *             lines nor arcs.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/offset.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/tolerance.hpp"

#include <array>
#include <benchmark/benchmark.h>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <numbers>
#include <vector>

/**
 *  Each benchmark offsets one closed polyline every step of its distance from one step to its
 *  count, both ways, round joins; its time is that of one offset. Its Loops counter is how many
 *  loops an offset gives, on average.
 */

namespace Geoxy::Bench
{
  namespace
  {
    /**
     *  A profile: a closed polyline, and the offsets asked of it.
     */
    struct Profile
    {
      Polyline    Shape = {};  /**< The polyline.         */
      std::size_t Count = 0;   /**< How many steps.       */
      double      Step  = 0.0; /**< How far each step is. */
    };

    /**
     *  A closed polyline from its vertices.
     *
     *  @param  vertices The vertices, each its x, y and bulge.
     *  @return          The polyline.
     */
    auto Closed(std::initializer_list<std::array<double, 3>> vertices) -> Polyline
    {
      Polyline polyline = {.Closed = true};
      for (const std::array<double, 3> &vertex : vertices)
      {
        polyline.Vertices.push_back(
            {.Position = {.X = vertex[0], .Y = vertex[1]}, .Bulge = vertex[2]});
      }
      return polyline;
    }

    auto Square() -> Profile
    {
      return {.Shape = Closed({{0, 0, 0}, {40, 0, 0}, {40, 40, 0}, {0, 40, 0}}),
              .Count = 30,
              .Step  = 1};
    }

    auto Diamond() -> Profile
    {
      return {.Shape = Closed({{0, 20, 0}, {20, 0, 0}, {40, 20, 0}, {20, 40, 0}}),
              .Count = 30,
              .Step  = 1};
    }

    auto Round() -> Profile
    {
      return {.Shape = Closed({{-40, 0, 1}, {40, 0, 1}}), .Count = 30, .Step = 1};
    }

    auto Rounded() -> Profile
    {
      const double bulge = std::tan(std::numbers::pi / 8.0);
      return {.Shape = Closed({{5, 0, 0},
                               {35, 0, bulge},
                               {40, 5, 0},
                               {40, 15, bulge},
                               {35, 20, 0},
                               {5, 20, bulge},
                               {0, 15, 0},
                               {0, 5, bulge}}),
              .Count = 30,
              .Step  = 0.5};
    }

    /**
     *  A plate with a round outer corner and a notch cut into one side by a half circle.
     *
     *  @return The profile.
     */
    auto Plate() -> Profile
    {
      const double bulge = std::tan(std::numbers::pi / 8.0);
      return {.Shape = Closed({{0, 0, 0},
                               {8, 0, -1},
                               {14, 0, 0},
                               {24, 0, 0},
                               {24, 10, bulge},
                               {18, 16, 0},
                               {0, 16, 0}}),
              .Count = 40,
              .Step  = 0.1};
    }

    /**
     *  A wheel of teeth: points by turns on two circles, each joined to the next by an arc bulging
     *  out.
     *
     *  @return The profile.
     */
    auto Gear() -> Profile
    {
      constexpr std::size_t Teeth    = 12;
      Polyline              polyline = {.Closed = true};
      for (std::size_t index = 0; index < 2 * Teeth; ++index)
      {
        const double angle  = static_cast<double>(index) * std::numbers::pi / Teeth;
        const double radius = index % 2 == 0 ? 20.0 : 15.0;
        polyline.Vertices.push_back(
            {.Position = {.X = radius * std::cos(angle), .Y = radius * std::sin(angle)},
             .Bulge    = 0.3});
      }
      return {.Shape = polyline, .Count = 40, .Step = 0.1};
    }

    /**
     *  A flower: a circle's points, each joined to the next by a petal that bulges out past a half
     *  circle.
     *
     *  @param  count How many petals.
     *  @return       The profile.
     */
    auto Flower(std::size_t count) -> Profile
    {
      Polyline polyline = {.Closed = true};
      for (std::size_t index = 0; index < count; ++index)
      {
        const double angle =
            static_cast<double>(index) * 2.0 * std::numbers::pi / static_cast<double>(count);
        polyline.Vertices.push_back(
            {.Position = {.X = 40 * std::cos(angle), .Y = 40 * std::sin(angle)}, .Bulge = 1.25});
      }
      return {.Shape = polyline, .Count = 30, .Step = 1};
    }

    /**
     *  Offsets a profile every step both ways, round after round.
     *
     *  @param state   The benchmark's state.
     *  @param profile The profile.
     */
    void Offsets(benchmark::State &state, const Profile &profile)
    {
      constexpr Tolerance Fine  = {.Length = 1e-9, .Radians = 1e-10};
      std::size_t         loops = 0;
      std::size_t         step  = 0;
      for ([[maybe_unused]] auto iteration : state)
      {
        const std::size_t multiple = (step / 2) + 1;
        const double      distance =
            profile.Step * static_cast<double>(multiple) * (step % 2 == 0 ? 1.0 : -1.0);
        auto found  = Offset(profile.Shape, {.Distance = distance, .Gaps = Gap::Fillet}, Fine);
        loops      += found->size();
        benchmark::DoNotOptimize(found);
        step = (step + 1) % (2 * profile.Count);
      }
      state.counters["Loops"] =
          benchmark::Counter(static_cast<double>(loops), benchmark::Counter::kAvgIterations);
    }

    void BySquare(benchmark::State &state)
    {
      Offsets(state, Square());
    }

    void ByDiamond(benchmark::State &state)
    {
      Offsets(state, Diamond());
    }

    void ByCircle(benchmark::State &state)
    {
      Offsets(state, Round());
    }

    void ByRounded(benchmark::State &state)
    {
      Offsets(state, Rounded());
    }

    void ByPlate(benchmark::State &state)
    {
      Offsets(state, Plate());
    }

    void ByGear(benchmark::State &state)
    {
      Offsets(state, Gear());
    }

    void ByFlower(benchmark::State &state)
    {
      Offsets(state, Flower(static_cast<std::size_t>(state.range(0))));
    }

    /**
     *  An ellipse offset out and in, to a millionth of its size.
     *
     *  @param state The benchmark's state.
     */
    void ByEllipse(benchmark::State &state)
    {
      const Ellipse ellipse = {
          .Map = {.X = {.X = 10}, .Y = {.Y = 5}}, .Start = Radians(0), .Sweep = Radians(2 * Pi)};
      constexpr Tolerance Fine = {.Length = 1e-9, .Radians = 1e-10};
      double              side = 1.0;
      for ([[maybe_unused]] auto iteration : state)
      {
        auto found = Offset(ellipse, {.Distance = 2 * side, .Approximation = {1e-5}}, Fine);
        benchmark::DoNotOptimize(found);
        side = -side;
      }
    }
  } // namespace

  BENCHMARK(BySquare);
  BENCHMARK(ByDiamond);
  BENCHMARK(ByCircle);
  BENCHMARK(ByRounded);
  BENCHMARK(ByPlate);
  BENCHMARK(ByGear);
  BENCHMARK(ByFlower)->Arg(10)->Arg(25)->Arg(50)->Arg(100);
  BENCHMARK(ByEllipse);
} // namespace Geoxy::Bench
