/***************************************************************************************************
 *  @file      tests/geoxy/boxes.cpp
 *  @brief     Tests of boxes kept by tags, against asking every box.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/boxes.hpp"

#include "draws.hpp"

#include "geoxy/box.hpp"
#include "geoxy/point.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Every box, by its tag: what the boxes are asked against.
     */
    using Ledger = std::map<std::uint64_t, Box>;

    /**
     *  Whether a box is kept in the tree, not apart: not empty, with finite sides.
     *
     *  @param  box The box.
     *  @return     Whether it is.
     */
    auto Treed(const Box &box) -> bool
    {
      return !IsEmpty(box) && std::isfinite((box.Max.X - box.Min.X) + (box.Max.Y - box.Min.Y));
    }

    /**
     *  How far a point lies from a box.
     *
     *  @param  box   The box.
     *  @param  point The point.
     *  @return       The distance; 0 within it.
     */
    auto Gap(const Box &box, Point point) -> double
    {
      const double across = std::max({box.Min.X - point.X, 0.0, point.X - box.Max.X});
      const double upward = std::max({box.Min.Y - point.Y, 0.0, point.Y - box.Max.Y});
      return std::hypot(across, upward);
    }

    /**
     *  A tag as packed from a kind and a slot, sparse over 64 bits.
     *
     *  @param  draws The draws.
     *  @return       The tag.
     */
    auto Tagged(Draws &draws) -> std::uint64_t
    {
      return ((draws.Next() % 8) << 40U) | (draws.Next() % 4096);
    }

    /**
     *  A box about a place: mostly finite and small, now and then empty, a point, or with a side
     *  that is not finite.
     *
     *  @param  draws The draws.
     *  @param  about The place, across and upward.
     *  @return       The box.
     */
    auto Drawn(Draws &draws, double about) -> Box
    {
      const double x    = draws.Near(about, 100.0);
      const double y    = draws.Near(about, 100.0);
      const double wide = std::abs(draws.Near(0.0, 8.0));
      const double high = std::abs(draws.Near(0.0, 8.0));
      switch (draws.Next() % 32)
      {
      case 0:
        return {};
      case 1:
        return {.Min = {.X = x, .Y = y}, .Max = {.X = x, .Y = y}};
      case 2:
        return {.Min = {.X = x, .Y = y},
                .Max = {.X = std::numeric_limits<double>::infinity(), .Y = y}};
      case 3:
        return {.Min = {.X = std::numeric_limits<double>::quiet_NaN(), .Y = y},
                .Max = {.X = x, .Y = y + high}};
      default:
        return {.Min = {.X = x, .Y = y}, .Max = {.X = x + wide, .Y = y + high}};
      }
    }

    /**
     *  The tags whose boxes overlap a box, by asking every one.
     *
     *  @param  ledger The boxes.
     *  @param  box    The box.
     *  @return        Their tags, in order.
     */
    auto Overlapping(const Ledger &ledger, const Box &box) -> std::vector<std::uint64_t>
    {
      std::vector<std::uint64_t> overlapping;
      for (const auto &[tag, kept] : ledger)
      {
        if (Overlaps(kept, box))
        {
          overlapping.push_back(tag);
        }
      }
      return overlapping;
    }

    /**
     *  The tags whose boxes lie wholly within a box, by asking every one.
     *
     *  @param  ledger The boxes.
     *  @param  box    The box.
     *  @return        Their tags, in order.
     */
    auto Within(const Ledger &ledger, const Box &box) -> std::vector<std::uint64_t>
    {
      std::vector<std::uint64_t> within;
      for (const auto &[tag, kept] : ledger)
      {
        if (!IsEmpty(kept) && Contains(box, kept.Min) && Contains(box, kept.Max))
        {
          within.push_back(tag);
        }
      }
      return within;
    }

    /**
     *  The box of every finite box, by asking every one.
     *
     *  @param  ledger The boxes.
     *  @return        The box.
     */
    auto Extents(const Ledger &ledger) -> Box
    {
      Box all;
      for (const auto &[tag, kept] : ledger)
      {
        if (Treed(kept))
        {
          all = Joined(all, kept);
        }
      }
      return all;
    }

    /**
     *  The tags nearest a point, by asking every one, nearest first.
     *
     *  @param  ledger The boxes.
     *  @param  point  The point.
     *  @param  most   How many at most, and how far at most.
     *  @return        The nearest, the lower tag first of those as near.
     */
    auto Nearest(const Ledger &ledger, Point point, const std::pair<std::size_t, double> &most)
        -> std::vector<Boxes::Close>
    {
      std::vector<Boxes::Close> all;
      for (const auto &[tag, kept] : ledger)
      {
        const double distance = Gap(kept, point);
        if (distance <= most.second)
        {
          all.push_back({.Tag = tag, .Distance = distance});
        }
      }
      std::ranges::sort(all, {}, [](const Boxes::Close &one) -> std::tuple<double, std::uint64_t>
                        { return {one.Distance, one.Tag}; });
      all.resize(std::min(all.size(), most.first));
      return all;
    }

    /**
     *  How many of the boxes' answers differ from asking every box, over puts, moves and erases
     *  drawn about a place.
     *
     *  @param  seed  The draws' seed.
     *  @param  about The place.
     *  @return       How many differ.
     */
    auto Wrong(std::uint64_t seed, double about) -> int
    {
      Draws                     draws{seed};
      Ledger                    ledger;
      std::vector<Boxes::Entry> entries;
      for (int index = 0; index < 300; ++index)
      {
        const Boxes::Entry entry = {.Tag = Tagged(draws), .Bounds = Drawn(draws, about)};
        entries.push_back(entry);
        ledger[entry.Tag] = entry.Bounds;
      }
      Boxes boxes{entries};
      int   wrong = 0;
      for (int round = 0; round < Rounds(4000); ++round)
      {
        const std::uint64_t tag = Tagged(draws);
        if (draws.Next() % 3 == 0)
        {
          boxes.Erase(tag);
          ledger.erase(tag);
        }
        else
        {
          const Box box = Drawn(draws, about);
          boxes.Insert(tag, box);
          ledger[tag] = box;
        }
        const double x     = draws.Near(about, 120.0);
        const double y     = draws.Near(about, 120.0);
        const double wide  = std::abs(draws.Near(0.0, 40.0));
        const Box    query = {.Min = {.X = x, .Y = y}, .Max = {.X = x + wide, .Y = y + wide}};
        const Point  point = {.X = x, .Y = y};
        const std::pair<std::size_t, double> most = {
            1 + (draws.Next() % 6), draws.Next() % 2 == 0 ? Boxes::Anywhere : wide};
        const auto distance = [&boxes, point](std::uint64_t which) -> double
        {
          return Gap(*boxes.Bounds(which), point);
        };
        wrong += static_cast<int>(boxes.Size() != ledger.size());
        wrong += static_cast<int>(boxes.Overlapping(query) != Overlapping(ledger, query));
        wrong += static_cast<int>(boxes.Within(query) != Within(ledger, query));
        wrong += static_cast<int>(boxes.Extents() != Extents(ledger));
        std::vector<Boxes::Close> found =
            boxes.Nearest(point, {.Count = most.first, .Within = most.second}, distance);
        std::ranges::sort(found, {},
                          [](const Boxes::Close &one) -> std::tuple<double, std::uint64_t>
                          { return {one.Distance, one.Tag}; });
        wrong += static_cast<int>(found != Nearest(ledger, point, most));
      }
      return wrong;
    }
  } // namespace

  TEST_CASE("boxes: every answer is asking every box's, near the origin and far from it")
  {
    for (const auto &[seed, about] :
         std::array{std::pair{0xB0C5U, 0.0}, std::pair{0xB0C6U, 5e6}, std::pair{0xB0C7U, 1e8}})
    {
      CHECK(Wrong(seed, about) == 0);
    }
  }

  TEST_CASE("boxes: a repeated tag keeps its last box, and a tag not kept has none")
  {
    const std::array entries = {
        Boxes::Entry{.Tag = 7, .Bounds = {.Min = {.X = 0, .Y = 0}, .Max = {.X = 1, .Y = 1}}},
        Boxes::Entry{.Tag = 7, .Bounds = {.Min = {.X = 5, .Y = 5}, .Max = {.X = 6, .Y = 6}}}};
    const Boxes boxes{entries};
    CHECK(boxes.Size() == 1);
    CHECK(boxes.Holds(7));
    CHECK_FALSE(boxes.Holds(8));
    CHECK(boxes.Bounds(7) == entries[1].Bounds);
    CHECK_FALSE(boxes.Bounds(8));
  }

  TEST_CASE("boxes: extents shrink when the farthest box goes, and none are left empty")
  {
    Boxes boxes;
    CHECK(IsEmpty(boxes.Extents()));
    boxes.Insert(1, {.Min = {.X = 0, .Y = 0}, .Max = {.X = 1, .Y = 1}});
    boxes.Insert(2, {.Min = {.X = 9, .Y = 9}, .Max = {.X = 10, .Y = 10}});
    CHECK(boxes.Extents() == Box{.Min = {.X = 0, .Y = 0}, .Max = {.X = 10, .Y = 10}});
    boxes.Insert(2, {.Min = {.X = 2, .Y = 2}, .Max = {.X = 3, .Y = 3}});
    CHECK(boxes.Extents() == Box{.Min = {.X = 0, .Y = 0}, .Max = {.X = 3, .Y = 3}});
    boxes.Erase(2);
    boxes.Erase(2);
    CHECK(boxes.Extents() == Box{.Min = {.X = 0, .Y = 0}, .Max = {.X = 1, .Y = 1}});
    boxes.Erase(1);
    CHECK(boxes.Size() == 0);
    CHECK(IsEmpty(boxes.Extents()));
    CHECK(boxes.Overlapping({.Min = {.X = -1e9, .Y = -1e9}, .Max = {.X = 1e9, .Y = 1e9}}).empty());
    CHECK(boxes.Nearest({}, {.Count = 0}, [](std::uint64_t) -> double { return 0.0; }).empty());
  }

  TEST_CASE("boxes: asked for more than are kept, Nearest gives every one, nearest first")
  {
    Boxes boxes;
    boxes.Insert(5, {.Min = {.X = 4, .Y = 0}, .Max = {.X = 5, .Y = 1}});
    boxes.Insert(6, {.Min = {.X = 1, .Y = 0}, .Max = {.X = 2, .Y = 1}});
    const auto distance = [&boxes](std::uint64_t tag) -> double
    {
      return boxes.Bounds(tag)->Min.X;
    };
    const std::vector<Boxes::Close> nearest = boxes.Nearest({}, {.Count = 9}, distance);
    CHECK(nearest
          == std::vector<Boxes::Close>{{.Tag = 6, .Distance = 1}, {.Tag = 5, .Distance = 4}});
  }

  TEST_CASE("boxes: a visitor stops a query, and the walk outward goes on until it says stop")
  {
    Boxes boxes;
    for (std::uint64_t tag = 0; tag < 40; ++tag)
    {
      const auto corner = static_cast<double>(tag % 8);
      boxes.Insert(tag,
                   {.Min = {.X = corner, .Y = corner}, .Max = {.X = corner + 1, .Y = corner + 1}});
    }
    boxes.Insert(99, {.Min = {.X = std::numeric_limits<double>::quiet_NaN(), .Y = 0},
                      .Max = {.X = 1, .Y = 1}});
    const Box  all  = {.Min = {.X = -1, .Y = -1}, .Max = {.X = 20, .Y = 20}};
    int        seen = 0;
    const auto stop = [&seen](std::uint64_t) -> bool
    {
      return ++seen < 3;
    };
    boxes.Overlapping(all, stop);
    CHECK(seen == 3);
    seen = 0;
    boxes.Within(all, stop);
    CHECK(seen == 3);
    const auto distance = [&boxes](std::uint64_t tag) -> double
    {
      return boxes.Bounds(tag)->Min.X;
    };
    std::vector<Boxes::Close> walked;
    boxes.Nearest({}, distance,
                  [&walked](const Boxes::Close &close) -> bool
                  {
                    walked.push_back(close);
                    return walked.size() < 7;
                  });
    REQUIRE(walked.size() == 7);
    CHECK(walked[0]       == Boxes::Close{.Tag = 0, .Distance = 0});
    CHECK(walked[4]       == Boxes::Close{.Tag = 32, .Distance = 0});
    CHECK(walked[5]       == Boxes::Close{.Tag = 1, .Distance = 1});
    CHECK(boxes.Nearest({}, {.Count = 50, .Within = 1}, distance).size() == 10);

    Boxes endless;
    for (const std::uint64_t tag : {std::uint64_t{1}, std::uint64_t{2}})
    {
      endless.Insert(tag, {.Min = {.X = 0, .Y = 0},
                           .Max = {.X = std::numeric_limits<double>::infinity(), .Y = 1}});
    }
    int told = 0;
    endless.Overlapping({.Min = {.X = 5, .Y = 0}, .Max = {.X = 6, .Y = 1}},
                        [&told](std::uint64_t) -> bool { return ++told < 1; });
    CHECK(told == 1);
  }

  TEST_CASE("boxes: copies keep trees of their own, and a move leaves none behind")
  {
    Boxes one;
    one.Insert(3, {.Min = {.X = 0, .Y = 0}, .Max = {.X = 1, .Y = 1}});
    Boxes two{one};
    two.Insert(4, {.Min = {.X = 2, .Y = 2}, .Max = {.X = 3, .Y = 3}});
    CHECK(one.Size() == 1);
    CHECK(two.Size() == 2);
    Boxes three;
    three = two;
    three.Erase(3);
    CHECK(two.Holds(3));
    Boxes four{std::move(three)};
    CHECK(four.Size() == 1);
    Boxes five;
    five = std::move(four);
    CHECK(five.Holds(4));
    CHECK(Boxes::Close{.Tag = 1, .Distance = 2}       == Boxes::Close{.Tag = 1, .Distance = 2});
    CHECK_FALSE(Boxes::Close{.Tag = 1, .Distance = 2} == Boxes::Close{.Tag = 2, .Distance = 2});
    CHECK_FALSE(Boxes::Close{.Tag = 1, .Distance = 2} == Boxes::Close{.Tag = 1, .Distance = 3});
  }
} // namespace Geoxy
