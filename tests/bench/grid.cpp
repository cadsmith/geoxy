/***************************************************************************************************
 *  @file      tests/bench/grid.cpp
 *  @brief     Benchmarks of meeting a whole grid: floor plans of ten and a hundred thousand curves
 *             made ready and gridded, their pairs found, met kind by kind, met one by one, and met
 *             by brute force; and pairs of lines scattered at random, by what they meet in.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/grid.hpp"

#include "../geoxy/draws.hpp"
#include "gridwork.hpp"
#include "plain.hpp"
#include "pooled.hpp"
#include "sane.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/workers.hpp"

#include <algorithm>
#include <array>
#include <benchmark/benchmark.h>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <numbers>
#include <span>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A floor plan of about some curves, as the audit's probe draws it: rooms of 10 by 10, each
     *  with double walls, a door's swing and leaf, a column, a line and an arc of furniture; and a
     *  few tiling lines across the whole plan, one each way for every 2,000 curves.
     *
     *  @param  count About how many curves.
     *  @return       The curves.
     */
    auto GridPlan(std::size_t count) -> std::vector<Curve>
    {
      Draws      draws{7};
      const auto unit = [&draws]() -> double
      {
        return draws.Near(0.5, 0.5);
      };
      const std::size_t  rooms = std::max<std::size_t>(1, count / 10);
      const auto         side  = static_cast<std::size_t>(std::ceil(std::sqrt(rooms)));
      std::vector<Curve> curves;
      for (std::size_t room = 0; room < rooms; ++room)
      {
        const std::size_t row = room / side;
        const double      x   = (10.0 * static_cast<double>(room % side)) + (0.05 * unit());
        const double      y   = (10.0 * static_cast<double>(row)) + (0.05 * unit());
        curves.emplace_back(Line{.Start = {.X = x, .Y = y}, .End = {.X = x + 10, .Y = y}});
        curves.emplace_back(
            Line{.Start = {.X = x, .Y = y + 0.2}, .End = {.X = x + 3, .Y = y + 0.2}});
        curves.emplace_back(
            Line{.Start = {.X = x + 4, .Y = y + 0.2}, .End = {.X = x + 10, .Y = y + 0.2}});
        curves.emplace_back(Line{.Start = {.X = x, .Y = y}, .End = {.X = x, .Y = y + 10}});
        curves.emplace_back(
            Line{.Start = {.X = x + 0.2, .Y = y + 0.2}, .End = {.X = x + 0.2, .Y = y + 10}});
        curves.emplace_back(Arc{.Center = {.X = x + 3, .Y = y + 0.2},
                                .Radius = 1.0,
                                .Start  = Degrees(0),
                                .Sweep  = Degrees(90)});
        curves.emplace_back(
            Line{.Start = {.X = x + 3, .Y = y + 0.2}, .End = {.X = x + 3, .Y = y + 1.2}});
        curves.emplace_back(
            Circle({.X = x + 0.1 + (0.4 * unit()), .Y = y + 0.1 + (0.4 * unit())}, 0.35));
        const Point desk = {.X = x + 2 + (6 * unit()), .Y = y + 2 + (6 * unit())};
        curves.emplace_back(
            Line{.Start = desk,
                 .End   = {.X = desk.X + (3 * unit()) - 1.5, .Y = desk.Y + (3 * unit()) - 1.5}});
        curves.emplace_back(Arc{.Center = {.X = x + 2 + (6 * unit()), .Y = y + 2 + (6 * unit())},
                                .Radius = 0.5 + (2 * unit()),
                                .Start  = Degrees(360 * unit()),
                                .Sweep  = Degrees(30 + (300 * unit()))});
      }
      const double size = 10.0 * static_cast<double>(side);
      for (std::size_t axis = 0; axis < std::max<std::size_t>(1, count / 2000); ++axis)
      {
        const double where = size * unit();
        curves.emplace_back(
            Line{.Start = {.X = -5, .Y = where}, .End = {.X = size + 5, .Y = where}});
        curves.emplace_back(
            Line{.Start = {.X = where, .Y = -5}, .End = {.X = where, .Y = size + 5}});
      }
      return curves;
    }

    /**
     *  The tolerance of a plan: Tolerance::Of its extents.
     *
     *  @param  curves The plan's curves.
     *  @return        The tolerance.
     */
    auto Within(const std::vector<Curve> &curves) -> Tolerance
    {
      Box all;
      for (const Curve &curve : curves)
      {
        all = Joined(all, Extents(curve));
      }
      return Tolerance::Of(all);
    }

    /**
     *  Lines scattered at random over a square, of lengths from half to one and a half and of any
     *  heading, the square as wide as makes each come near about two others.
     *
     *  @param  count How many.
     *  @return       The lines.
     */
    auto Hatch(std::size_t count) -> std::vector<Curve>
    {
      Draws              draws{3};
      const double       side = 0.6 * std::sqrt(static_cast<double>(count));
      std::vector<Curve> curves;
      for (std::size_t index = 0; index < count; ++index)
      {
        const Point  start   = {.X = draws.Near(0, side), .Y = draws.Near(0, side)};
        const double heading = draws.Near(0, std::numbers::pi);
        const double length  = draws.Near(1, 0.5);
        curves.emplace_back(Line{.Start = start,
                                 .End   = {.X = start.X + (length * std::cos(heading)),
                                           .Y = start.Y + (length * std::sin(heading))}});
      }
      return curves;
    }

    /**
     *  The least of a benchmark's repetitions: what it takes when nothing else runs, as near as a
     *  loaded machine tells it.
     *
     *  @param  values The repetitions' times.
     *  @return        The least.
     */
    auto Least(const std::vector<double> &values) -> double
    {
      return std::ranges::min(values);
    }

    /**
     *  A counter of the time each of some items took, an iteration doing them all.
     *
     *  @param  items How many.
     *  @return       The counter.
     */
    auto Each(std::size_t items) -> benchmark::Counter
    {
      return {static_cast<double>(items),
              benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert};
    }

    /**
     *  A run of pairs met kind by kind, as IntersectAll meets them.
     *
     *  @param  grid The grid.
     *  @param  run  The pairs, their kind and the tolerance.
     *  @return      How many meet.
     */
    auto Batched(const Grid &grid, const Gridwork::Run &run) -> std::size_t
    {
      std::size_t met = 0;
      Gridwork::Met(grid, run, [&met](const Grid::Pair &, const Marking &) -> void { ++met; });
      return met;
    }

    /**
     *  A run of pairs met one by one through Intersect, as IntersectAll met them before.
     *
     *  @param  grid    The grid.
     *  @param  run     The pairs and the tolerance.
     *  @param  marking A marking to meet them into.
     *  @return         How many meet.
     */
    auto OneByOne(const Grid &grid, const Gridwork::Run &run, Marking &marking) -> std::size_t
    {
      std::size_t met = 0;
      for (const Grid::Pair &pair : run.Pairs)
      {
        Intersect({grid[pair[0]], grid[pair[1]]}, run.Within, marking);
        met += static_cast<std::size_t>(
            !marking.Marks.empty() || !marking.Overlaps.empty() || marking.Spent);
      }
      return met;
    }

    /**
     *  A plan made into a grid: its curves copied, made ready, boxed and gridded.
     */
    void GridMade(benchmark::State &state)
    {
      const std::vector<Curve> curves    = GridPlan(static_cast<std::size_t>(state.range(0)));
      const Tolerance          tolerance = Within(curves);
      for ([[maybe_unused]] auto iteration : state)
      {
        const Grid grid{curves, tolerance};
        benchmark::DoNotOptimize(grid.Size());
      }
      state.counters["curves"] = static_cast<double>(curves.size());
    }

    /**
     *  The pairs of a grid that may meet, found in its tiling and sorted by kind.
     */
    void GridPairs(benchmark::State &state)
    {
      const std::vector<Curve> curves    = GridPlan(static_cast<std::size_t>(state.range(0)));
      const Tolerance          tolerance = Within(curves);
      const Grid               grid{curves, tolerance};
      std::size_t              count = 0;
      for ([[maybe_unused]] auto iteration : state)
      {
        const Gridwork::Pairs pairs = Gridwork::Candidates(grid, tolerance);
        count                       = 0;
        for (const auto &kind : pairs)
        {
          count += kind.size();
        }
        benchmark::DoNotOptimize(count);
      }
      state.counters["pairs"] = static_cast<double>(count);
    }

    /**
     *  Those pairs met, each into one marking.
     */
    void GridNarrow(benchmark::State &state)
    {
      const std::vector<Curve> curves    = GridPlan(static_cast<std::size_t>(state.range(0)));
      const Tolerance          tolerance = Within(curves);
      const Grid               grid{curves, tolerance};
      const Gridwork::Pairs    pairs = Gridwork::Candidates(grid, tolerance);
      Marking                  marking;
      std::size_t              marks = 0;
      for ([[maybe_unused]] auto iteration : state)
      {
        marks = 0;
        for (const auto &kind : pairs)
        {
          for (const Grid::Pair &pair : kind)
          {
            Intersect({grid[pair[0]], grid[pair[1]]}, tolerance, marking);
            marks += marking.Marks.size();
          }
        }
        benchmark::DoNotOptimize(marks);
      }
      state.counters["marks"] = static_cast<double>(marks);
    }

    /**
     *  Those pairs met kind by kind, as IntersectAll meets them.
     */
    void GridBatched(benchmark::State &state)
    {
      const std::vector<Curve> curves    = GridPlan(static_cast<std::size_t>(state.range(0)));
      const Tolerance          tolerance = Within(curves);
      const Grid               grid{curves, tolerance};
      const Gridwork::Pairs    pairs = Gridwork::Candidates(grid, tolerance);
      std::size_t              marks = 0;
      const Grid::Sink         sink  = [&marks](const Grid::Pair &, const Marking &marking) -> void
      {
        marks += marking.Marks.size();
      };
      for ([[maybe_unused]] auto iteration : state)
      {
        marks = 0;
        for (std::size_t kind = 0; kind < Gridwork::Kinds; ++kind)
        {
          Gridwork::Met(grid,
                        {.Of     = static_cast<Gridwork::Kind>(kind),
                         .Pairs  = pairs.at(kind),
                         .Within = tolerance},
                        sink);
        }
        benchmark::DoNotOptimize(marks);
      }
      state.counters["marks"] = static_cast<double>(marks);
    }

    /**
     *  The pairs of one kind met as IntersectAll meets them, or one by one through Intersect: the
     *  plan's size, the kind, and 1 for one by one.
     */
    void GridKind(benchmark::State &state)
    {
      const std::vector<Curve> curves    = GridPlan(static_cast<std::size_t>(state.range(0)));
      const Tolerance          tolerance = Within(curves);
      const Grid               grid{curves, tolerance};
      const Gridwork::Pairs    pairs = Gridwork::Candidates(grid, tolerance);
      const auto               kind  = static_cast<Gridwork::Kind>(state.range(1));
      const Gridwork::Run      run   = {
          .Of = kind, .Pairs = pairs.at(static_cast<std::size_t>(kind)), .Within = tolerance};
      Marking     marking;
      std::size_t met = 0;
      for ([[maybe_unused]] auto iteration : state)
      {
        met = state.range(2) == 0 ? Batched(grid, run) : OneByOne(grid, run, marking);
        benchmark::DoNotOptimize(met);
      }
      state.counters["pair"] = Each(run.Pairs.size());
      state.counters["met"]  = static_cast<double>(met);
    }

    /**
     *  The pairs of one kind of a plan that meet where an end of one lies on the other, a corner or
     *  a joint of two walls, a door's swing on its wall, or the rest of that kind: the plan's size,
     *  the kind, 1 for those at an end and 0 for the rest, and 1 for one by one through Intersect.
     */
    void GridEnded(benchmark::State &state)
    {
      const std::vector<Curve> curves    = GridPlan(static_cast<std::size_t>(state.range(0)));
      const Tolerance          tolerance = Within(curves);
      const Grid               grid{curves, tolerance};
      const Gridwork::Pairs    pairs = Gridwork::Candidates(grid, tolerance);
      const auto               kind  = static_cast<Gridwork::Kind>(state.range(1));
      std::vector<Grid::Pair>  chosen;
      Marking                  marking;
      for (const Grid::Pair &pair : pairs.at(static_cast<std::size_t>(kind)))
      {
        Intersect({grid[pair[0]], grid[pair[1]]}, tolerance, marking);
        const bool ended = std::ranges::any_of(marking.Marks, [](const Mark &mark) -> bool
                                               { return mark.Kind == Contact::End; });
        if (ended == (state.range(2) == 1))
        {
          chosen.push_back(pair);
        }
      }
      const Gridwork::Run run = {.Of = kind, .Pairs = chosen, .Within = tolerance};
      std::size_t         met = 0;
      for ([[maybe_unused]] auto iteration : state)
      {
        met = state.range(3) == 0 ? Batched(grid, run) : OneByOne(grid, run, marking);
        benchmark::DoNotOptimize(met);
      }
      state.counters["pair"]  = Each(run.Pairs.size());
      state.counters["pairs"] = static_cast<double>(run.Pairs.size());
    }

    /**
     *  Pairs of lines scattered at random, those their exact pass tells miss, 0, those it tells
     *  cross, 1, and those it leaves to Intersect, 2: met as IntersectAll meets them, or one by one
     *  through Intersect, 1 for one by one.
     */
    void GridLines(benchmark::State &state)
    {
      const std::vector<Curve>               curves    = Hatch(40'000);
      const Tolerance                        tolerance = Within(curves);
      const Grid                             grid{curves, tolerance};
      const Gridwork::Pairs                  pairs = Gridwork::Candidates(grid, tolerance);
      std::array<std::vector<Grid::Pair>, 3> told;
      Marking                                marking;
      for (const Grid::Pair &pair : pairs.at(static_cast<std::size_t>(Gridwork::Kind::Lines)))
      {
        const Plain plain = Plainly({grid[pair[0]], grid[pair[1]]}, Sane(tolerance), marking);
        told.at(plain.Told ? plain.Count : 2).push_back(pair);
      }
      const Gridwork::Run run = {.Of     = Gridwork::Kind::Lines,
                                 .Pairs  = told.at(static_cast<std::size_t>(state.range(0))),
                                 .Within = tolerance};
      std::size_t         met = 0;
      for ([[maybe_unused]] auto iteration : state)
      {
        met = state.range(1) == 0 ? Batched(grid, run) : OneByOne(grid, run, marking);
        benchmark::DoNotOptimize(met);
      }
      state.counters["pair"]  = Each(run.Pairs.size());
      state.counters["pairs"] = static_cast<double>(run.Pairs.size());
    }

    /**
     *  Every pair of a grid met through IntersectAll, its pairs found and met.
     */
    void GridIntersectAll(benchmark::State &state)
    {
      const std::vector<Curve> curves    = GridPlan(static_cast<std::size_t>(state.range(0)));
      const Tolerance          tolerance = Within(curves);
      const Grid               grid{curves, tolerance};
      std::size_t              met   = 0;
      std::size_t              marks = 0;
      for ([[maybe_unused]] auto iteration : state)
      {
        met   = 0;
        marks = 0;
        IntersectAll(grid, tolerance,
                     [&met, &marks](const Grid::Pair &, const Marking &marking) -> void
                     {
                       ++met;
                       marks += marking.Marks.size();
                     });
        benchmark::DoNotOptimize(marks);
      }
      state.counters["met"]   = static_cast<double>(met);
      state.counters["marks"] = static_cast<double>(marks);
    }

    /**
     *  Every pair of a grid met through IntersectAll on a pool kept from one call to the next,
     *  timed on the wall's clock, and the processor time of every thread together as its work: the
     *  plan's size, and the pool's threads, the calling thread among them; 0 for IntersectAll given
     *  no pool, which makes one for the call.
     */
    void GridPooled(benchmark::State &state)
    {
      const std::vector<Curve> curves    = GridPlan(static_cast<std::size_t>(state.range(0)));
      const Tolerance          tolerance = Within(curves);
      const Grid               grid{curves, tolerance};
      const auto               threads = static_cast<std::size_t>(state.range(1));
      Workers                  pool{std::max<std::size_t>(1, threads)};
      std::size_t              met   = 0;
      std::size_t              marks = 0;
      const Grid::Sink sink = [&met, &marks](const Grid::Pair &, const Marking &marking) -> void
      {
        ++met;
        marks += marking.Marks.size();
      };
      const std::clock_t began = std::clock();
      for ([[maybe_unused]] auto iteration : state)
      {
        met   = 0;
        marks = 0;
        if (threads == 0)
        {
          IntersectAll(grid, tolerance, sink);
        }
        else
        {
          IntersectAll({grid, pool}, tolerance, sink);
        }
        benchmark::DoNotOptimize(marks);
      }
      state.counters["met"]   = static_cast<double>(met);
      state.counters["marks"] = static_cast<double>(marks);
      state.counters["work"] =
          benchmark::Counter(static_cast<double>(std::clock() - began) / CLOCKS_PER_SEC,
                             benchmark::Counter::kAvgIterations);
    }

    /**
     *  A pool woken for a job that does nothing, and met again once each thread is done with it:
     *  what a call on a pool pays before any work, by the pool's threads.
     */
    void WorkersWoken(benchmark::State &state)
    {
      Workers           pool{static_cast<std::size_t>(state.range(0))};
      const Pooled::Job job = [](std::size_t slot) -> void
      {
        benchmark::DoNotOptimize(slot);
      };
      for ([[maybe_unused]] auto iteration : state)
      {
        benchmark::DoNotOptimize(Pooled::Run(pool, job));
      }
    }

    /**
     *  Every pair of a grid met by brute force, each pair made ready told far apart by Intersect.
     */
    void GridBrute(benchmark::State &state)
    {
      const std::vector<Curve> curves    = GridPlan(static_cast<std::size_t>(state.range(0)));
      const Tolerance          tolerance = Within(curves);
      const Grid               grid{curves, tolerance};
      Marking                  marking;
      std::size_t              marks = 0;
      for ([[maybe_unused]] auto iteration : state)
      {
        marks = 0;
        for (std::size_t first = 0; first < grid.Size(); ++first)
        {
          for (std::size_t second = first + 1; second < grid.Size(); ++second)
          {
            Intersect({grid[first], grid[second]}, tolerance, marking);
            marks += marking.Marks.size();
          }
        }
        benchmark::DoNotOptimize(marks);
      }
      state.counters["marks"] = static_cast<double>(marks);
    }

    /**
     *  One curve moved, and moved back: its entries taken out of the tiling and put in again.
     */
    void GridMoved(benchmark::State &state)
    {
      const std::vector<Curve> curves    = GridPlan(static_cast<std::size_t>(state.range(0)));
      const Tolerance          tolerance = Within(curves);
      Grid                     grid{curves, tolerance};
      const std::size_t        index = curves.size() / 2;
      const Curve              moved = Line{.Start = {.X = 1, .Y = 1}, .End = {.X = 25, .Y = 3}};
      for ([[maybe_unused]] auto iteration : state)
      {
        grid.Moved(index, moved);
        grid.Moved(index, curves[index]);
      }
    }

    /**
     *  One curve met with a grid: a line across a few rooms.
     */
    void GridIntersectWith(benchmark::State &state)
    {
      const std::vector<Curve> curves    = GridPlan(static_cast<std::size_t>(state.range(0)));
      const Tolerance          tolerance = Within(curves);
      const Grid               grid{curves, tolerance};
      const Curve              across = Line{.Start = {.X = 1, .Y = 1}, .End = {.X = 25, .Y = 3}};
      const Prepared           ready{across};
      std::size_t              met = 0;
      for ([[maybe_unused]] auto iteration : state)
      {
        met = 0;
        IntersectWith({ready, grid}, tolerance,
                      [&met](const Grid::Pair &, const Marking &) -> void { ++met; });
        benchmark::DoNotOptimize(met);
      }
      state.counters["met"] = static_cast<double>(met);
    }
  } // namespace

  BENCHMARK(GridMade)->Arg(10'000)->Arg(100'000)->Unit(benchmark::kMillisecond);
  BENCHMARK(GridPairs)->Arg(10'000)->Arg(100'000)->Unit(benchmark::kMillisecond);
  BENCHMARK(GridNarrow)
      ->Arg(10'000)
      ->Arg(100'000)
      ->Unit(benchmark::kMillisecond)
      ->ComputeStatistics("min", Least);
  BENCHMARK(GridBatched)
      ->Arg(10'000)
      ->Arg(100'000)
      ->Unit(benchmark::kMillisecond)
      ->ComputeStatistics("min", Least);
  BENCHMARK(GridKind)
      ->ArgsProduct({{10'000, 100'000}, {0, 1, 2, 3}, {0, 1}})
      ->Unit(benchmark::kMicrosecond)
      ->ComputeStatistics("min", Least);
  BENCHMARK(GridEnded)
      ->ArgsProduct({{10'000, 100'000}, {0, 2}, {0, 1}, {0, 1}})
      ->Unit(benchmark::kMicrosecond)
      ->ComputeStatistics("min", Least);
  BENCHMARK(GridLines)
      ->ArgsProduct({{0, 1, 2}, {0, 1}})
      ->Unit(benchmark::kMicrosecond)
      ->ComputeStatistics("min", Least);
  BENCHMARK(GridIntersectAll)
      ->Arg(10'000)
      ->Arg(100'000)
      ->Unit(benchmark::kMillisecond)
      ->ComputeStatistics("min", Least);
  BENCHMARK(GridPooled)
      ->ArgsProduct({{10'000, 100'000}, {1, 2, 4, 8, 10, 0}})
      ->Unit(benchmark::kMillisecond)
      ->UseRealTime()
      ->ComputeStatistics("min", Least);
  BENCHMARK(WorkersWoken)->Arg(2)->Arg(4)->Arg(10)->Unit(benchmark::kMicrosecond)->UseRealTime();
  BENCHMARK(GridBrute)->Arg(10'000)->Unit(benchmark::kMillisecond)->Iterations(1);
  BENCHMARK(GridMoved)->Arg(10'000)->Arg(100'000)->Unit(benchmark::kMicrosecond);
  BENCHMARK(GridIntersectWith)->Arg(10'000)->Arg(100'000)->Unit(benchmark::kMicrosecond);
} // namespace Geoxy
