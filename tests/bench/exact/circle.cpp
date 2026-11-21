/***************************************************************************************************
 *  @file      tests/bench/exact/circle.cpp
 *  @brief     Benchmarks of the circle test, beside Shewchuk's own.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/circle.hpp"

#include "../../geoxy/draws.hpp"
#include "../../geoxy/exact/original.hpp"
#include "exact/triangle.hpp"

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <benchmark/benchmark.h>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <vector>

namespace Geoxy::Exact
{
  namespace
  {
    struct Case
    {
      Triangle Circle;
      Point    Point;
    };

    constexpr std::size_t CircleCount = 65'536;

    /**
     *  Places spread over the plane: the first stage decides them all.
     */
    auto CircleSpread() -> std::vector<Case>
    {
      Draws draws{11};

      const auto next = [&draws] -> Geoxy::Point
      {
        return {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)};
      };

      std::vector<Case> all;
      all.reserve(CircleCount);
      for (std::size_t index = 0; index < CircleCount; ++index)
      {
        all.push_back(
            {.Circle = {.First = next(), .Second = next(), .Third = next()}, .Point = next()});
      }
      return all;
    }

    /**
     *  Places on the unit circle, rounded: the later stages decide them, the compensated second
     *  stage most.
     */
    auto CircleNear() -> std::vector<Case>
    {
      Draws draws{12};

      const auto onCircle = [&draws] -> Geoxy::Point
      {
        const double angle = draws.Near(0, std::numbers::pi);
        return {.X = std::cos(angle), .Y = std::sin(angle)};
      };

      std::vector<Case> all;
      all.reserve(CircleCount);
      for (std::size_t index = 0; index < CircleCount; ++index)
      {
        all.push_back({.Circle = {.First = onCircle(), .Second = onCircle(), .Third = onCircle()},
                       .Point  = onCircle()});
      }
      return all;
    }

    /**
     *  Integer places on a circle about a survey-scale centre, from the eight images of (3, 4):
     *  every difference exact and every determinant zero, so the exact sum decides each, in
     *  integers.
     */
    auto Lattice() -> std::vector<Case>
    {
      static constexpr std::array<Vector, 8> Around{
          Vector{.X =  3, .Y =  4},
          Vector{.X = -3, .Y =  4},
          Vector{.X =  3, .Y = -4},
          Vector{.X = -3, .Y = -4},
          Vector{.X =  4, .Y =  3},
          Vector{.X = -4, .Y =  3},
          Vector{.X =  4, .Y = -3},
          Vector{.X = -4, .Y = -3},
      };

      Draws draws{13};

      std::vector<Case> all;
      all.reserve(CircleCount);
      for (std::size_t index = 0; index < CircleCount; ++index)
      {
        const Geoxy::Point centre = {.X = std::round(draws.Near(5e6, 1e3)),
                                     .Y = std::round(draws.Near(-3e6, 1e3))};
        const auto         next   = [&draws, centre] -> Geoxy::Point
        {
          return centre + Around.at(draws.Next() % Around.size());
        };
        all.push_back(
            {.Circle = {.First = next(), .Second = next(), .Third = next()}, .Point = next()});
      }
      return all;
    }

    auto CircleReference(const Case &one) -> double
    {
      std::array<double, 2> a = {one.Circle.First.X, one.Circle.First.Y};
      std::array<double, 2> b = {one.Circle.Second.X, one.Circle.Second.Y};
      std::array<double, 2> c = {one.Circle.Third.X, one.Circle.Third.Y};
      std::array<double, 2> d = {one.Point.X, one.Point.Y};
      return incircle(a.data(), b.data(), c.data(), d.data());
    }

    template <typename TestT>
    void CircleRun(benchmark::State &state, const std::vector<Case> &all, TestT test)
    {
      for ([[maybe_unused]] auto iteration : state)
      {
        for (const Case &one : all)
        {
          benchmark::DoNotOptimize(test(one));
        }
      }
      state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(all.size()));
    }

    /**
     *  Shewchuk's incircle over a set of cases.
     */
    template <typename MakeT>
    void CircleByOriginal(benchmark::State &state, MakeT make)
    {
      exactinit();
      CircleRun(state, make(), CircleReference);
    }

    /**
     *  InCircle over a set of cases.
     */
    template <typename MakeT>
    void ByInCircle(benchmark::State &state, MakeT make)
    {
      CircleRun(state, make(),
                [](const Case &one) -> int { return InCircle(one.Circle, one.Point); });
    }
  } // namespace

  BENCHMARK_CAPTURE(CircleByOriginal, CircleSpread, CircleSpread);
  BENCHMARK_CAPTURE(ByInCircle, CircleSpread, CircleSpread);
  BENCHMARK_CAPTURE(CircleByOriginal, CircleNear, CircleNear);
  BENCHMARK_CAPTURE(ByInCircle, CircleNear, CircleNear);
  BENCHMARK_CAPTURE(CircleByOriginal, CircleLattice, Lattice);
  BENCHMARK_CAPTURE(ByInCircle, CircleLattice, Lattice);
} // namespace Geoxy::Exact
