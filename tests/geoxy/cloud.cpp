/***************************************************************************************************
 *  @file      tests/geoxy/cloud.cpp
 *  @brief     Tests for points kept in a tree: every answer against asking every point.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/cloud.hpp"

#include "draws.hpp"

#include "geoxy/box.hpp"
#include "geoxy/point.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Where points are drawn: about a place, within a spread.
     */
    struct Patch
    {
      double About  = 0.0; /**< The middle, across and upward. */
      double Spread = 1.0; /**< How far from it at most.       */
    };

    /**
     *  A point somewhere in a patch, on a coarse lattice so that many are as far from another.
     *
     *  @param  draws The draws.
     *  @param  patch The patch.
     *  @return       The point.
     */
    auto Dotted(Draws &draws, const Patch &patch) -> Point
    {
      const double step = patch.Spread / 64;
      return {.X = patch.About + (step * static_cast<double>(draws.Next() % 128)) - patch.Spread,
              .Y = patch.About + (step * static_cast<double>(draws.Next() % 128)) - patch.Spread};
    }

    /**
     *  The points nearest a point, by asking every one.
     *
     *  @param  cloud The cloud, and one past its greatest index.
     *  @param  point The point.
     *  @param  most  How many at most, and how far at most.
     *  @return       The nearest, as Cloud::Nearest orders them.
     */
    auto Asked(const std::pair<const Cloud &, std::size_t> &cloud, Point point,
               const std::pair<std::size_t, double> &most) -> std::vector<Cloud::Close>
    {
      std::vector<Cloud::Close> all;
      for (std::size_t index = 0; index < cloud.second; ++index)
      {
        if (cloud.first.Holds(index) && Distance(cloud.first[index], point) <= most.second)
        {
          all.push_back({.Index = index, .Distance = Distance(cloud.first[index], point)});
        }
      }
      std::ranges::sort(all, {}, [](const Cloud::Close &one) -> std::tuple<double, std::size_t>
                        { return {one.Distance, one.Index}; });
      all.resize(std::min(all.size(), most.first));
      return all;
    }

    /**
     *  How many of a cloud's answers differ from asking every point, over boxes and points drawn in
     *  its patch.
     *
     *  @param  cloud The cloud, and one past its greatest index.
     *  @param  patch The patch.
     *  @param  draws The draws.
     *  @return       The count, 0 for a cloud that answers as asking every point does.
     */
    auto Wrong(const std::pair<const Cloud &, std::size_t> &cloud, const Patch &patch, Draws &draws)
        -> int
    {
      int wrong = 0;
      for (int round = 0; round < 60; ++round)
      {
        const Point corner = Dotted(draws, patch);
        const Box   box    = {
            .Min = corner,
            .Max = {.X = corner.X + (patch.Spread / 4), .Y = corner.Y + (patch.Spread / 8)}};
        std::vector<std::size_t> held;
        for (std::size_t index = 0; index < cloud.second; ++index)
        {
          if (cloud.first.Holds(index) && Contains(box, cloud.first[index]))
          {
            held.push_back(index);
          }
        }
        wrong               += cloud.first.Within(box) == held ? 0 : 1;
        const double within  = round % 2 == 0 ? Cloud::Anywhere : patch.Spread * (1 + round) / 256;
        for (const std::size_t count : {std::size_t{1}, std::size_t{9}, cloud.second + 1})
        {
          const auto got   = cloud.first.Nearest(corner, count, within);
          const auto want  = Asked(cloud, corner, {count, within});
          wrong           += got == want ? 0 : 1;
        }
      }
      return wrong;
    }

    /**
     *  The patches the tests draw in: small, of unit size, and far from the origin.
     */
    constexpr std::array Patches{
        Patch{.About =   0, .Spread = 1e-6},
        Patch{.About =   0, .Spread =   10},
        Patch{.About = 5e6, .Spread =  100},
    };
  } // namespace

  TEST_CASE("cloud: a cloud of no point answers nothing")
  {
    const Cloud none;
    CHECK(none.Size() == 0);
    CHECK(none.Within(Box{.Min = {.X = -1, .Y = -1}, .Max = {.X = 1, .Y = 1}}).empty());
    CHECK(none.Nearest(Point{}, 3, Cloud::Anywhere).empty());
    const Cloud packed{std::span<const Point>{}};
    CHECK(packed.Nearest(Point{}, 3, Cloud::Anywhere).empty());
  }

  TEST_CASE("cloud: points packed at once are found as asking every point finds them")
  {
    for (const Patch &patch : Patches)
    {
      for (const std::size_t count : {std::size_t{9}, std::size_t{700}, std::size_t{5000}})
      {
        Draws              draws{61};
        std::vector<Point> points;
        points.reserve(count);
        for (std::size_t index = 0; index < count; ++index)
        {
          points.push_back(Dotted(draws, patch));
        }
        const Cloud cloud{points};
        CHECK(cloud.Size() == count);
        CHECK(cloud.Nearest(points.front(), 0, Cloud::Anywhere).empty());
        CHECK(Wrong({cloud, count}, patch, draws) == 0);
      }
    }
  }

  TEST_CASE("cloud: points put in and taken out one at a time are found as asking every point does")
  {
    for (const Patch &patch : Patches)
    {
      Draws       draws{71};
      Cloud       cloud;
      std::size_t top = 0;
      for (int index = 0; index < 6000; ++index)
      {
        top = std::max(top, cloud.Insert(Dotted(draws, patch)) + 1);
      }
      CHECK(cloud.Size() == 6000);
      CHECK(Wrong({cloud, top}, patch, draws) == 0);
      for (std::size_t index = 0; index < top; ++index)
      {
        if (index % 5 != 0)
        {
          cloud.Erase(index);
        }
      }
      cloud.Erase(1);
      cloud.Erase(top + 3);
      CHECK(cloud.Size() == 1200);
      CHECK(Wrong({cloud, top}, patch, draws) == 0);
      for (int index = 0; index < 500; ++index)
      {
        CHECK(cloud.Insert(Dotted(draws, patch)) < top);
      }
      CHECK(Wrong({cloud, top}, patch, draws) == 0);
      for (std::size_t index = 0; index < top; ++index)
      {
        cloud.Erase(index);
      }
      CHECK(cloud.Size() == 0);
      CHECK(cloud.Within(Box{.Min = {.X = -1e9, .Y = -1e9}, .Max = {.X = 1e9, .Y = 1e9}}).empty());
      for (int index = 0; index < 900; ++index)
      {
        CHECK(cloud.Insert(Dotted(draws, patch)) < top);
      }
      CHECK(Wrong({cloud, top}, patch, draws) == 0);
    }
  }

  TEST_CASE("cloud: points that are not finite are asked every time, and never the nearest")
  {
    constexpr double NaN    = std::numeric_limits<double>::quiet_NaN();
    constexpr double Far    = std::numeric_limits<double>::infinity();
    const std::array points = {Point{.X = 1, .Y = 1}, Point{.X = NaN, .Y = 0},
                               Point{.X = 2, .Y = 2}, Point{.X = Far, .Y = 0}};
    Cloud            cloud{points};
    CHECK(cloud.Within(Box{.Min = {.X = 0, .Y = 0}, .Max = {.X = 3, .Y = 3}})
          == std::vector<std::size_t>{0, 2});
    CHECK(cloud.Within(Box{.Min = {.X = 0, .Y = 0}, .Max = {.X = Far, .Y = 3}})
          == std::vector<std::size_t>{0, 2, 3});
    CHECK(cloud.Nearest(Point{}, 4, 1e300).size()           == 2);
    CHECK(cloud.Nearest(Point{}, 4, Cloud::Anywhere).size() == 3);
    const std::size_t more = cloud.Insert(Point{.X = 0, .Y = NaN});
    CHECK(cloud.Nearest(Point{.X = NaN, .Y = 0}, 4, Cloud::Anywhere).empty());
    cloud.Erase(1);
    cloud.Erase(more);
    CHECK(cloud.Size() == 3);
    CHECK(cloud.Nearest(Point{}, 4, Cloud::Anywhere).size() == 3);
  }

  TEST_CASE("cloud: a cloud copied or taken over holds the same points at the same indices")
  {
    Draws draws{81};
    Cloud cloud;
    for (int index = 0; index < 200; ++index)
    {
      cloud.Insert(Dotted(draws, Patches[1]));
    }
    cloud.Erase(3);
    Cloud copied{cloud};
    copied.Erase(4);
    CHECK(cloud.Holds(4));
    CHECK(copied.Size() == 198);
    Cloud assigned;
    assigned = cloud;
    CHECK(assigned.Size() == 199);
    Cloud taken{std::move(assigned)};
    Cloud given;
    given = std::move(taken);
    CHECK(given.Insert(Point{}) == 3);
    CHECK(Cloud::Close{.Index = 1, .Distance = 2}       == Cloud::Close{.Index = 1, .Distance = 2});
    CHECK_FALSE(Cloud::Close{.Index = 1, .Distance = 2} == Cloud::Close{.Index = 3, .Distance = 2});
    CHECK_FALSE(Cloud::Close{.Index = 1, .Distance = 2} == Cloud::Close{.Index = 1, .Distance = 4});
    CHECK(Wrong({given, 200}, Patches[1], draws)        == 0);
  }
} // namespace Geoxy
