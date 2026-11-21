/***************************************************************************************************
 *  @file      tests/geoxy/scratch.cpp
 *  @brief     Tests for the pool each thread keeps: lists drawing on it, meetings in many threads
 *             at once each the same as alone.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "scratch.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/tolerance.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <memory_resource>
#include <thread>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether two meetings are the same, hit for hit and stretch for stretch.
     *
     *  @param  one   One.
     *  @param  other The other.
     *  @return       Whether they are.
     */
    auto ScratchIdentical(const Intersections &one, const Intersections &other) -> bool
    {
      const auto hit = [](const Hit &lhs, const Hit &rhs) -> bool
      {
        return lhs.At     == rhs.At
            && lhs.First  == rhs.First
            && lhs.Second == rhs.Second
            && lhs.Kind   == rhs.Kind;
      };
      const auto overlap = [](const Overlap &lhs, const Overlap &rhs) -> bool
      {
        return lhs.First == rhs.First && lhs.Second == rhs.Second && lhs.Along == rhs.Along;
      };
      return std::ranges::equal(one.Hits, other.Hits, hit)
          && std::ranges::equal(one.Overlaps, other.Overlaps, overlap)
          && one.Spent == other.Spent;
    }
  } // namespace

  TEST_CASE("scratch: lists draw on the thread's pool, and give back to it")
  {
    List<double> values;
    for (int index = 0; index < 1000; ++index)
    {
      values.push_back(index);
    }
    CHECK(values.size() == 1000);
    CHECK(values.back() == 999);
    const List<double> copied = values;
    CHECK(copied            == values);
    CHECK(Scratch<double>{} == Scratch<double>{Scratch<int>{}});
  }

  TEST_CASE("scratch: the arena as a memory resource, shared by every thread, is itself alone")
  {
    Scratchpad &shared = Scratchpad::Shared();
    CHECK(shared.is_equal(Scratchpad::Shared()));
    CHECK_FALSE(shared.is_equal(*std::pmr::new_delete_resource()));
    const std::size_t        before = Arena::Local().Kept();
    std::pmr::vector<double> values(&shared);
    values.assign(1000, 1.0);
    CHECK(values.size() == 1000);
    CHECK(Arena::Local().Kept() >= before);
  }

  TEST_CASE("scratch: the arena grows past its first stretch, and starts over in one it keeps")
  {
    /**
     *  A list grown far past the first stretch, beside one held throughout, so that blocks are
     *  given back out of turn; then all given back, and the same taken again in one stretch.
     */
    for (int round = 0; round < 3; ++round)
    {
      List<int>    held(3, round);
      List<double> grown;
      for (int index = 0; index < 200000; ++index)
      {
        grown.push_back(index);
      }
      const List<double> copied(grown.begin(), grown.begin() + 1000);
      CHECK(grown.back()  == 199999);
      CHECK(copied.back() == 999);
      CHECK(held.front()  == round);
    }
    List<char> none;
    none.reserve(0);
    CHECK(none.empty());
    const Scratch<char> one;
    char               *room = Scratch<char>{one}.allocate(0);
    Scratch<char>{one}.deallocate(room, 0);
    CHECK(Scratch<double>{} == Scratch<double>{one});
  }

  TEST_CASE("scratch: blocks given back out of turn keep the arena within its ceiling")
  {
    /**
     *  A block held throughout, and each round one taken, one taken after it, and the first given
     *  back out of turn: its room stays taken until the arena starts over. Without a ceiling the
     *  arena would hold every such block, 18 MiB over these rounds.
     */
    constexpr std::size_t Ceiling = std::size_t{1} << 20U;
    Arena                 arena{Ceiling};
    void                 *held = arena.Take(64);
    for (std::size_t round = 0; round < 1000; ++round)
    {
      const std::size_t size  = 4096 * (1 + (round % 8));
      void             *first = arena.Take(size);
      void             *after = arena.Take(1024);
      arena.Give(first, size);
      arena.Give(after, 1024);
      CHECK(arena.Kept() <= Ceiling);
    }
    const std::size_t kept = arena.Kept();
    void             *past = arena.Take(Ceiling / 2);
    CHECK(arena.Kept() == kept);
    arena.Give(past, Ceiling / 2);
    arena.Give(held, 64);
    CHECK(arena.Kept() <= Ceiling / 16);
    CHECK(Arena::Ceiling == std::size_t{1} << 26U);
  }

  TEST_CASE("scratch: a block past the ceiling is the system's, and the arena starts over small")
  {
    constexpr std::size_t Ceiling = std::size_t{1} << 20U;
    constexpr std::size_t First   = std::size_t{1} << 16U;
    Arena                 arena{Ceiling};
    void                 *whole = arena.Take(2 * Ceiling);
    CHECK(arena.Kept() == 0);
    arena.Give(whole, 2 * Ceiling);
    CHECK(arena.Kept() == 0);
    void *small = arena.Take(100);
    CHECK(arena.Kept() == First);
    void *large = arena.Take(200000);
    CHECK(arena.Kept() > First);
    arena.Give(large, 200000);
    arena.Give(small, 100);
    CHECK(arena.Kept() == First);

    Arena other{Ceiling};
    void *wide = other.Take(100000);
    CHECK(other.Kept() > Ceiling / 16);
    other.Give(wide, 100000);
    CHECK(other.Kept() == 0);
    void *again = other.Take(16);
    CHECK(other.Kept() == First);
    other.Give(again, 16);
  }

  TEST_CASE("scratch: meetings in many threads at once are each as alone")
  {
    const Curve                wave   = Bezier::Cubic{.Start = {.X = 0, .Y = 0},
                                                      .Out   = {.X = 1, .Y = 3},
                                                      .In    = {.X = 2, .Y = -3},
                                                      .End   = {.X = 3, .Y = 0}};
    const std::array<Curve, 3> others = {
        Curve{Line{.Start = {.X = -1, .Y = 0.5}, .End = {.X = 4, .Y = 0.5}}},
        Curve{Arc{.Center = {.X = 1.5, .Y = 0}, .Radius = 1.25, .Start = {}, .Sweep = Turn}},
        Curve{Bezier::Cubic{.Start = {.X = 0, .Y = 1},
                            .Out   = {.X = 1, .Y = -2},
                            .In    = {.X = 2, .Y = 4},
                            .End   = {.X = 3, .Y = -1}}}};
    const Tolerance            tolerance = {.Length = 1e-9, .Radians = 1e-9};
    std::vector<Intersections> alone;
    alone.reserve(others.size());
    for (const Curve &other : others)
    {
      alone.push_back(Intersect(wave, other, tolerance));
    }
    std::atomic<int> differ = 0;
    {
      std::vector<std::jthread> threads;
      threads.reserve(8);
      for (int thread = 0; thread < 8; ++thread)
      {
        threads.emplace_back(
            [&]() -> void
            {
              for (int round = 0; round < 200; ++round)
              {
                const std::size_t index = static_cast<std::size_t>(round) % others.size();
                if (!ScratchIdentical(Intersect(wave, others.at(index), tolerance),
                                      alone.at(index)))
                {
                  ++differ;
                }
              }
            });
      }
    }
    CHECK(differ == 0);
  }

  TEST_CASE("scratch: one meeting refilled pair after pair is as each alone, and keeps its room")
  {
    const Curve                wave   = Bezier::Cubic{.Start = {.X = 0, .Y = 0},
                                                      .Out   = {.X = 1, .Y = 3},
                                                      .In    = {.X = 2, .Y = -3},
                                                      .End   = {.X = 3, .Y = 0}};
    const std::array<Curve, 3> others = {
        Curve{Line{.Start = {.X = -1, .Y = 0.5}, .End = {.X = 4, .Y = 0.5}}},
        Curve{Arc{.Center = {.X = 1.5, .Y = 0}, .Radius = 1.25, .Start = {}, .Sweep = Turn}},
        Curve{wave}};
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    Intersections   kept;
    kept.Hits.reserve(16);
    kept.Overlaps.reserve(16);
    const Hit     *room    = kept.Hits.data();
    const Overlap *shelves = kept.Overlaps.data();
    for (int round = 0; round < 30; ++round)
    {
      for (const Curve &other : others)
      {
        Intersect({wave, other}, tolerance, kept);
        CHECK(ScratchIdentical(kept, Intersect(wave, other, tolerance)));
        CHECK(kept.Hits.data()     == room);
        CHECK(kept.Overlaps.data() == shelves);
      }
    }
  }
} // namespace Geoxy
