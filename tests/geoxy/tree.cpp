/***************************************************************************************************
 *  @file      tests/geoxy/tree.cpp
 *  @brief     Tests for curves kept in a tree: every answer against asking every curve.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/tree.hpp"

#include "draws.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/box.hpp"
#include "geoxy/contain.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
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
     *  Where curves are drawn: about a place, within a spread, each of a size.
     */
    struct Field
    {
      double About  = 0.0; /**< The middle, across and upward. */
      double Spread = 1.0; /**< How far from it at most.       */
      double Size   = 1.0; /**< How large a curve is at most.  */
    };

    /**
     *  A curve of some kind somewhere in a field: a line, an arc, or a polyline with an arc.
     *
     *  @param  draws The draws.
     *  @param  field The field.
     *  @return       The curve.
     */
    auto Drawn(Draws &draws, const Field &field) -> Curve
    {
      const Point  start = {.X = draws.Near(field.About, field.Spread),
                            .Y = draws.Near(field.About, field.Spread)};
      const Point  end   = {.X = start.X + draws.Near(0, field.Size),
                            .Y = start.Y + draws.Near(0, field.Size)};
      const Point  last  = {.X = end.X + draws.Near(0, field.Size),
                            .Y = end.Y + draws.Near(0, field.Size)};
      const double kind  = draws.Near(0, 1);
      if (kind < -0.3)
      {
        return Line{.Start = start, .End = end};
      }
      if (kind < 0.3)
      {
        return Arc{.Center = start,
                   .Radius = field.Size * (1.5 + draws.Near(0, 1)) / 4,
                   .Start  = Degrees(draws.Near(0, 180)),
                   .Sweep  = Degrees(draws.Near(0, 300))};
      }
      return Polyline{.Vertices = {{.Position = start, .Bulge = draws.Near(0, 1)},
                                   {.Position = end, .Bulge = 0},
                                   {.Position = last, .Bulge = 0}},
                      .Closed   = false};
    }

    /**
     *  A box somewhere in a field, a few curves wide.
     *
     *  @param  draws The draws.
     *  @param  field The field.
     *  @return       The box.
     */
    auto Boxed(Draws &draws, const Field &field) -> Box
    {
      const Point corner = {.X = draws.Near(field.About, field.Spread),
                            .Y = draws.Near(field.About, field.Spread)};
      return {.Min = corner,
              .Max = {.X = corner.X + (field.Size * (2 + draws.Near(0, 2))),
                      .Y = corner.Y + (field.Size * (2 + draws.Near(0, 2)))}};
    }

    /**
     *  The indices a tree holds, in order.
     *
     *  @param  tree The tree.
     *  @param  top  One past the greatest index it may hold.
     *  @return      Them.
     */
    auto Held(const Tree &tree, std::size_t top) -> std::vector<std::size_t>
    {
      std::vector<std::size_t> held;
      for (std::size_t index = 0; index < top; ++index)
      {
        if (tree.Holds(index))
        {
          held.push_back(index);
        }
      }
      return held;
    }

    /**
     *  The curves nearest a point, by asking every one.
     *
     *  @param  tree  The tree, and one past its greatest index.
     *  @param  point The point.
     *  @param  most  How many at most, and how far at most.
     *  @return       The nearest, as Tree::Nearest orders them.
     */
    auto Asked(const std::pair<const Tree &, std::size_t> &tree, Point point,
               const std::pair<std::size_t, double> &most) -> std::vector<Tree::Close>
    {
      std::vector<Tree::Close> all;
      for (const std::size_t index : Held(tree.first, tree.second))
      {
        const Place  place = Nearest(tree.first[index], point);
        const double far   = Distance(place.At, point);
        if (far <= most.second)
        {
          all.push_back({.Index = index, .Nearest = place, .Distance = far});
        }
      }
      std::ranges::sort(all, {}, [](const Tree::Close &one) -> std::tuple<double, std::size_t>
                        { return {one.Distance, one.Index}; });
      all.resize(std::min(all.size(), most.first));
      return all;
    }

    /**
     *  Where a ray strikes the curves, by meeting every one.
     *
     *  @param  tree      The tree, and one past its greatest index.
     *  @param  ray       The ray.
     *  @param  tolerance The tolerance.
     *  @return           The strikes, as Tree::Cast orders them.
     */
    auto Met(const std::pair<const Tree &, std::size_t> &tree, const Ray &ray,
             const Tolerance &tolerance) -> std::vector<Tree::Strike>
    {
      std::vector<Tree::Strike> all;
      const Curve               aimed = ray;
      for (const std::size_t index : Held(tree.first, tree.second))
      {
        const Intersections met = Intersect(aimed, tree.first[index], tolerance);
        for (const Hit &hit : met.Hits)
        {
          all.push_back({.Index = index, .At = hit.At, .Along = hit.First.Along});
        }
        for (const Overlap &overlap : met.Overlaps)
        {
          all.push_back(
              {.Index = index, .At = At(aimed, overlap.First.Start), .Along = overlap.First.Start});
          all.push_back(
              {.Index = index, .At = At(aimed, overlap.First.End), .Along = overlap.First.End});
        }
      }
      std::ranges::sort(all, {}, [](const Tree::Strike &one) -> std::tuple<double, std::size_t>
                        { return {one.Along, one.Index}; });
      return all;
    }

    /**
     *  How many of a tree's answers differ from asking every curve, over boxes, points and rays
     *  drawn in its field.
     *
     *  @param  tree  The tree, and one past its greatest index.
     *  @param  field The field.
     *  @param  draws The draws.
     *  @return       The count, 0 for a tree that answers as asking every curve does.
     */
    auto Wrong(const std::pair<const Tree &, std::size_t> &tree, const Field &field, Draws &draws)
        -> int
    {
      int                            wrong     = 0;
      const Tolerance                tolerance = {.Length = 1e-9 * field.Size, .Radians = 1e-9};
      const std::vector<std::size_t> held      = Held(tree.first, tree.second);
      for (int round = 0; round < 40; ++round)
      {
        const Box                box    = Boxed(draws, field);
        std::vector<std::size_t> inside = held;
        std::vector<std::size_t> across = held;
        std::erase_if(inside, [&tree, &box](std::size_t index) -> bool
                      { return !Inside(tree.first[index], box); });
        std::erase_if(across, [&tree, &box](std::size_t index) -> bool
                      { return !Intersects(tree.first[index], box); });
        const std::vector<std::size_t> near  = tree.first.Near(box);
        wrong                               += tree.first.Window(box) == inside ? 0 : 1;
        wrong                               += tree.first.Crossing(box) == across ? 0 : 1;
        wrong                               += std::ranges::includes(near, across) ? 0 : 1;

        const Point  point  = box.Min;
        const double within = round % 2 == 0 ? Tree::Anywhere : field.Size * (1 + round) / 8;
        for (const std::size_t count : {std::size_t{1}, std::size_t{7}, held.size() + 1})
        {
          wrong += tree.first.Nearest(point, count, within) == Asked(tree, point, {count, within})
                     ? 0
                     : 1;
        }
        const Ray ray  = {.Start = point,
                          .Way   = {.X = draws.Near(0, field.Size), .Y = draws.Near(0, field.Size)}};
        wrong         += tree.first.Cast(ray, tolerance) == Met(tree, ray, tolerance) ? 0 : 1;
      }
      return wrong;
    }

    /**
     *  The fields the tests draw in: small, of unit size, far from the origin, and of curves of
     *  very unlike sizes.
     */
    constexpr std::array Fields{
        Field{.About =   0, .Spread = 1e-3, .Size = 1e-5},
        Field{.About =   0, .Spread =   10, .Size =  0.5},
        Field{.About = 5e6, .Spread =  100, .Size =    3},
        Field{.About =   0, .Spread = 1000, .Size =  400},
    };
  } // namespace

  TEST_CASE("tree: a tree of no curve answers nothing, and one of few answers as asking each")
  {
    const Tree none;
    CHECK(none.Size() == 0);
    CHECK(none.Near(Box{.Min = {.X = -1, .Y = -1}, .Max = {.X = 1, .Y = 1}}).empty());
    CHECK(none.Nearest(Point{}, 3, Tree::Anywhere).empty());
    CHECK(none.Cast(Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}, Tolerance{}).empty());

    const Tree packed{std::span<const Curve>{}};
    CHECK(packed.Size() == 0);
    CHECK(packed.Nearest(Point{}, 1, Tree::Anywhere).empty());

    Draws              draws{11};
    std::vector<Curve> few;
    few.reserve(5);
    for (int index = 0; index < 5; ++index)
    {
      few.push_back(Drawn(draws, Fields[1]));
    }
    const Tree tree{few};
    CHECK(tree.Size() == 5);
    CHECK(tree.Nearest(Point{}, 0, Tree::Anywhere).empty());
    CHECK(Wrong({tree, few.size()}, Fields[1], draws) == 0);
  }

  TEST_CASE(
      "tree: curves packed at once are found as asking every curve finds them, at every scale")
  {
    for (const Field &field : Fields)
    {
      Draws              draws{21};
      std::vector<Curve> curves;
      curves.reserve(900);
      for (int index = 0; index < 900; ++index)
      {
        const Field sized = {.About  = field.About,
                             .Spread = field.Spread,
                             .Size   = index % 97 == 0 ? field.Spread : field.Size};
        curves.push_back(Drawn(draws, sized));
      }
      const Tree tree{curves};
      CHECK(tree.Size() == curves.size());
      CHECK(Wrong({tree, curves.size()}, field, draws) == 0);
    }
  }

  TEST_CASE("tree: curves put in and taken out one at a time are found as asking every curve does")
  {
    for (const Field &field : Fields)
    {
      Draws       draws{31};
      Tree        tree;
      std::size_t top = 0;
      for (int index = 0; index < 700; ++index)
      {
        top = std::max(top, tree.Insert(Drawn(draws, field)) + 1);
      }
      CHECK(tree.Size() == 700);
      CHECK(Wrong({tree, top}, field, draws) == 0);
      for (std::size_t index = 0; index < top; index += 3)
      {
        tree.Erase(index);
      }
      tree.Erase(0);
      tree.Erase(top + 5);
      CHECK(tree.Size() == 700 - 234);
      CHECK_FALSE(tree.Holds(0));
      CHECK(Wrong({tree, top}, field, draws) == 0);
      const std::size_t again = tree.Insert(Drawn(draws, field));
      CHECK(again < top);
      for (int index = 0; index < 300; ++index)
      {
        top = std::max(top, tree.Insert(Drawn(draws, field)) + 1);
      }
      CHECK(Wrong({tree, top}, field, draws) == 0);
      for (std::size_t index = 0; index < top; ++index)
      {
        tree.Erase(index);
      }
      CHECK(tree.Size() == 0);
      CHECK(tree.Nearest(Point{}, 2, Tree::Anywhere).empty());
      CHECK(tree.Insert(Drawn(draws, field)) < top);
      CHECK(Wrong({tree, top}, field, draws) == 0);
    }
  }

  TEST_CASE("tree: curves with no box are asked every time")
  {
    constexpr double   NaN = std::numeric_limits<double>::quiet_NaN();
    Draws              draws{41};
    const Field       &field = Fields[1];
    std::vector<Curve> curves;
    curves.reserve(64);
    for (int index = 0; index < 60; ++index)
    {
      curves.push_back(Drawn(draws, field));
    }
    curves.emplace_back(Ray{.Start = {.X = 1, .Y = 2}, .Way = {.X = 1, .Y = 1}});
    curves.emplace_back(Xline{.Through = {.X = -3, .Y = 0}, .Way = {.X = 0, .Y = 1}});
    curves.emplace_back(Line{.Start = {.X = NaN, .Y = 0}, .End = {.X = 1, .Y = 1}});
    Tree tree{curves};
    CHECK(Wrong({tree, curves.size()}, field, draws) == 0);
    const std::size_t ray = tree.Insert(Ray{.Start = {.X = 0, .Y = 5}, .Way = {.X = -1, .Y = 0}});
    CHECK(Wrong({tree, ray + 1}, field, draws) == 0);
    tree.Erase(60);
    tree.Erase(ray);
    CHECK(tree.Size() == curves.size() - 1);
    CHECK(Wrong({tree, ray + 1}, field, draws) == 0);
    CHECK(tree.Nearest(Point{.X = NaN, .Y = 0}, 3, Tree::Anywhere).empty());
  }

  /**
   *  Far from the origin Intersect takes a length under 16 steps of a double as that many steps, so
   *  a line a few steps off a ray's shares a stretch with it whatever finer tolerance is asked: the
   *  tree must not pass its box by.
   */
  TEST_CASE("tree: a ray strikes what Intersect finds within the steps of a double, far away")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    for (const double base : {1e8, 5e6, 1e3, 1.0})
    {
      const Ray   ray   = {.Start = {.X = 0, .Y = base}, .Way = {.X = 1, .Y = 0}};
      const Curve aimed = ray;
      double      above = base;
      double      below = base;
      for (int steps = 0; steps <= 40; ++steps)
      {
        for (const double level : {above, below})
        {
          Tree        tree;
          const Curve line = Line{.Start = {.X = 10, .Y = level}, .End = {.X = 20, .Y = level}};
          const std::size_t   index = tree.Insert(line);
          const Intersections found = Intersect(aimed, line, tolerance);
          CHECK(tree.Cast(ray, tolerance).size()
                == found.Hits.size() + (2 * found.Overlaps.size()));
          CHECK(Tree{std::span{&line, 1}}.Cast(ray, tolerance).size()
                == tree.Cast(ray, tolerance).size());
          CHECK(Tree{tree}.Cast(ray, tolerance).size() == tree.Cast(ray, tolerance).size());
          CHECK(index == 0);
        }
        above = std::nextafter(above, 2 * base);
        below = std::nextafter(below, 0.0);
      }
    }
  }

  TEST_CASE("tree: a ray along a line strikes it at both ends of what they share")
  {
    const std::vector<Curve> curves = {Line{.Start = {.X = 2, .Y = 0}, .End = {.X = 5, .Y = 0}},
                                       Line{.Start = {.X = 3, .Y = -1}, .End = {.X = 3, .Y = 1}},
                                       Line{.Start = {.X = -4, .Y = -1}, .End = {.X = -4, .Y = 1}},
                                       Line{.Start = {.X = 1, .Y = 3}, .End = {.X = 9, .Y = 3}}};
    const Tree               tree{curves};
    const Ray                ray     = {.Start = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 0}};
    const auto               strikes = tree.Cast(ray, Tolerance{.Length = 1e-9, .Radians = 1e-9});
    REQUIRE(strikes.size() == 3);
    CHECK(strikes[0].Index == 0);
    CHECK(strikes[0].At    == Point{.X = 2, .Y = 0});
    CHECK(strikes[1].Index == 1);
    CHECK(strikes[1].At    == Point{.X = 3, .Y = 0});
    CHECK(strikes[2].Index == 0);
    CHECK(strikes[2].At    == Point{.X = 5, .Y = 0});
  }

  TEST_CASE("tree: a tree copied or taken over holds the same curves at the same indices")
  {
    Draws draws{51};
    Tree  tree;
    for (int index = 0; index < 100; ++index)
    {
      tree.Insert(Drawn(draws, Fields[1]));
    }
    tree.Erase(7);
    const Point point   = {.X = 1, .Y = 1};
    const auto  nearest = tree.Nearest(point, 4, Tree::Anywhere);

    Tree copied{tree};
    CHECK((copied.Nearest(point, 4, Tree::Anywhere) == nearest));
    copied.Erase(nearest.front().Index);
    CHECK((tree.Nearest(point, 4, Tree::Anywhere) == nearest));
    CHECK(copied.Insert(Line{}) == nearest.front().Index);

    Tree assigned;
    assigned = tree;
    CHECK(assigned.Size() == 99);
    CHECK((assigned.Nearest(point, 4, Tree::Anywhere) == nearest));

    Tree taken{std::move(assigned)};
    CHECK((taken.Nearest(point, 4, Tree::Anywhere) == nearest));
    Tree given;
    given = std::move(taken);
    CHECK((given.Nearest(point, 4, Tree::Anywhere) == nearest));
    CHECK(given.Insert(Line{}) == 7);

    const Tree::Close  close  = nearest.front();
    const Tree::Strike strike = {.Index = 1, .At = {.X = 2, .Y = 3}, .Along = 4};
    CHECK(close == close);
    CHECK_FALSE(close
                == Tree::Close{.Index    = close.Index + 1,
                               .Nearest  = close.Nearest,
                               .Distance = close.Distance});
    CHECK_FALSE(close
                == Tree::Close{.Index = close.Index, .Nearest = {}, .Distance = close.Distance});
    CHECK_FALSE(close
                == Tree::Close{.Index = close.Index, .Nearest = close.Nearest, .Distance = -1});
    CHECK(strike       == strike);
    CHECK_FALSE(strike == Tree::Strike{.Index = 2, .At = strike.At, .Along = 4});
    CHECK_FALSE(strike == Tree::Strike{.Index = 1, .At = {}, .Along = 4});
    CHECK_FALSE(strike == Tree::Strike{.Index = 1, .At = strike.At, .Along = 5});
  }
} // namespace Geoxy
