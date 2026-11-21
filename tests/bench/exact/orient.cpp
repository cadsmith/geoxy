/***************************************************************************************************
 *  @file      tests/bench/exact/orient.cpp
 *  @brief     Benchmarks of the orientation test, beside Shewchuk's own.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/orient.hpp"

#include "../../geoxy/draws.hpp"
#include "../../geoxy/exact/original.hpp"

#include "geoxy/point.hpp"

#include <array>
#include <benchmark/benchmark.h>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

/**
 *  Each run goes over 65,536 cases, 3 MiB, held in the second level cache, so what is measured is
 *  the test and not memory; and too many for a branch predictor to learn by heart, as it learns a
 *  few thousand repeated, which would flatter Shewchuk's sign branches.
 */

namespace Geoxy::Exact
{
  namespace
  {
    struct Corners
    {
      Point First;
      Point Second;
      Point Third;
    };

    constexpr std::size_t OrientCount = 65'536;

    /**
     *  Corners spread over the plane: the first stage decides them all.
     */
    auto OrientSpread() -> std::vector<Corners>
    {
      Draws draws{1};

      std::vector<Corners> all;
      all.reserve(OrientCount);
      for (std::size_t index = 0; index < OrientCount; ++index)
      {
        const Point first  = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)};
        const Point second = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)};
        const Point third  = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)};
        all.push_back({.First = first, .Second = second, .Third = third});
      }
      return all;
    }

    /**
     *  Corners a step off a line, at survey coordinates.
     */
    auto OrientNear() -> std::vector<Corners>
    {
      Draws draws{2};

      std::vector<Corners> all;
      all.reserve(OrientCount);
      for (std::size_t index = 0; index < OrientCount; ++index)
      {
        const Point a = {.X = draws.Near(5e6, 1e3), .Y = draws.Near(-3e6, 1e3)};
        const Point b = {.X = draws.Near(5e6, 1e3), .Y = draws.Near(-3e6, 1e3)};

        const double along = draws.Near(0.5, 2);
        const Point  c     = {.X = std::nextafter(a.X + (along * (b.X - a.X)), HUGE_VAL),
                              .Y = a.Y + (along * (b.Y - a.Y))};
        all.push_back({.First = a, .Second = b, .Third = c});
      }
      return all;
    }

    /**
     *  Integers just off a line: exact differences, so Orient decides by one fused comparison.
     */
    auto Integers() -> std::vector<Corners>
    {
      Draws draws{3};

      std::vector<Corners> all;
      all.reserve(OrientCount);
      for (std::size_t index = 0; index < OrientCount; ++index)
      {
        const double base   = 0x1.0p51 + static_cast<double>(draws.Next() >> 13U) - 0x1.0p21;
        const double step   = 1.0 + static_cast<double>(draws.Next() >> 44U);
        const Point  origin = {.X = static_cast<double>(draws.Next() >> 24U),
                               .Y = static_cast<double>(draws.Next() >> 24U)};
        all.push_back({.First  = {.X = origin.X + base + step, .Y = origin.Y + base},
                       .Second = {.X = origin.X + base, .Y = origin.Y + base - step},
                       .Third  = origin});
      }
      return all;
    }

    /**
     *  Corners exactly on a diagonal whose differences round: only the exact stage decides them.
     */
    auto Diagonal() -> std::vector<Corners>
    {
      Draws draws{4};

      std::vector<Corners> all;
      all.reserve(OrientCount);
      for (std::size_t index = 0; index < OrientCount; ++index)
      {
        const double x = draws.Near(0, 1);
        const double y = draws.Near(0.25, 1e-3);
        const double z = draws.Near(5e6, 1);
        all.push_back(
            {.First = {.X = x, .Y = x}, .Second = {.X = y, .Y = y}, .Third = {.X = z, .Y = z}});
      }
      return all;
    }

    auto OrientReference(const Corners &one) -> double
    {
      std::array<double, 2> first  = {one.First.X, one.First.Y};
      std::array<double, 2> second = {one.Second.X, one.Second.Y};
      std::array<double, 2> third  = {one.Third.X, one.Third.Y};
      return orient2d(first.data(), second.data(), third.data());
    }

    template <typename TestT>
    void OrientRun(benchmark::State &state, const std::vector<Corners> &all, TestT test)
    {
      for ([[maybe_unused]] auto iteration : state)
      {
        for (const Corners &one : all)
        {
          benchmark::DoNotOptimize(test(one));
        }
      }
      state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(all.size()));
    }

    /**
     *  Shewchuk's orient2d over a set of corners.
     */
    template <typename MakeT>
    void OrientByOriginal(benchmark::State &state, MakeT make)
    {
      exactinit();
      OrientRun(state, make(), OrientReference);
    }

    /**
     *  Orient over a set of corners.
     */
    template <typename MakeT>
    void ByOrient(benchmark::State &state, MakeT make)
    {
      OrientRun(state, make(),
                [](const Corners &one) -> int { return Orient(one.First, one.Second, one.Third); });
    }

    /**
     *  The plain determinant, for scale: what the first stage costs at least.
     */
    void Plain(benchmark::State &state)
    {
      OrientRun(state, OrientSpread(),
                [](const Corners &one) -> double
                {
                  const Point &a = one.First;
                  const Point &b = one.Second;
                  const Point &c = one.Third;
                  return ((a.X - c.X) * (b.Y - c.Y)) - ((a.Y - c.Y) * (b.X - c.X));
                });
    }
  } // namespace

  BENCHMARK(Plain);
  BENCHMARK_CAPTURE(OrientByOriginal, Spread, OrientSpread);
  BENCHMARK_CAPTURE(ByOrient, Spread, OrientSpread);
  BENCHMARK_CAPTURE(OrientByOriginal, Near, OrientNear);
  BENCHMARK_CAPTURE(ByOrient, Near, OrientNear);
  BENCHMARK_CAPTURE(OrientByOriginal, Integers, Integers);
  BENCHMARK_CAPTURE(ByOrient, Integers, Integers);
  BENCHMARK_CAPTURE(OrientByOriginal, Diagonal, Diagonal);
  BENCHMARK_CAPTURE(ByOrient, Diagonal, Diagonal);
} // namespace Geoxy::Exact
