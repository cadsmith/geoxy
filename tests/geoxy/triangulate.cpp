/***************************************************************************************************
 *  @file      tests/geoxy/triangulate.cpp
 *  @brief     Tests of fills cut into triangles: every triangle counterclockwise and inside, no
 *             side used twice the same way, their areas summing exactly to the region's, and the
 *             Delaunay fill's circles empty.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/triangulate.hpp"

#include "exact/circle.hpp"
#include "expansion.hpp"

#include "geoxy/orient.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <random>
#include <set>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Twice the signed area of a triangle, exactly.
     *
     *  @param  first  The first corner.
     *  @param  second The second.
     *  @param  third  The third.
     *  @return        The area, as an expansion.
     */
    auto Twice(Point first, Point second, Point third) -> Expansion::Number
    {
      using namespace Expansion;
      const Number alongX = Difference(Of(second.X), Of(first.X));
      const Number alongY = Difference(Of(second.Y), Of(first.Y));
      const Number upX    = Difference(Of(third.X), Of(first.X));
      const Number upY    = Difference(Of(third.Y), Of(first.Y));
      return Difference(Product(alongX, upY), Product(alongY, upX));
    }

    /**
     *  Twice the area a loop encloses, taken positive, exactly.
     *
     *  @param  loop The loop.
     *  @return      The area, as an expansion.
     */
    auto Enclosed(const std::vector<Point> &loop) -> Expansion::Number
    {
      Expansion::Number sum = Expansion::Of(0.0);
      for (std::size_t index = 1; index + 1 < loop.size(); ++index)
      {
        sum = Expansion::Sum(sum, Twice(loop[0], loop[index], loop[index + 1]));
      }
      return Expansion::Sign(sum) < 0 ? Expansion::Difference(Expansion::Of(0.0), sum) : sum;
    }

    /**
     *  Whether a point is within loops by the even-odd rule, by a ray to its right.
     *
     *  @param  loops The loops.
     *  @param  point The point, off every side.
     *  @return       Whether it is.
     */
    auto Within(std::span<const std::vector<Point>> loops, Point point) -> bool
    {
      bool inside = false;
      for (const std::vector<Point> &loop : loops)
      {
        for (std::size_t index = 0; index < loop.size(); ++index)
        {
          const Point from = loop[index];
          const Point next = loop[(index + 1) % loop.size()];
          if ((from.Y > point.Y) != (next.Y > point.Y)
              && point.X < from.X + ((point.Y - from.Y) * (next.X - from.X) / (next.Y - from.Y)))
          {
            inside = !inside;
          }
        }
      }
      return inside;
    }

    /**
     *  Every loop's points, in order.
     *
     *  @param  loops The loops.
     *  @return       The points.
     */
    auto Flat(std::span<const std::vector<Point>> loops) -> std::vector<Point>
    {
      std::vector<Point> points;
      for (const std::vector<Point> &loop : loops)
      {
        points.insert(points.end(), loop.begin(), loop.end());
      }
      return points;
    }

    /**
     *  A triangle's corners as points.
     *
     *  @param  points   The loops' points.
     *  @param  triangle The triangle.
     *  @return          Its corners.
     */
    auto Cornered(std::span<const Point> points, const std::array<std::uint32_t, 3> &triangle)
        -> std::array<Point, 3>
    {
      return {points[triangle.front()], points[triangle[1]], points[triangle.back()]};
    }

    /**
     *  How many triangles turn clockwise or not at all, or have their centroid outside the region.
     *
     *  @param  loops The region's loops.
     *  @param  mesh  The fill.
     *  @return       The count.
     */
    auto Misplaced(std::span<const std::vector<Point>> loops, const Mesh &mesh) -> std::size_t
    {
      const std::vector<Point> points = Flat(loops);
      std::size_t              wrong  = 0;
      for (const auto &triangle : mesh.Triangles)
      {
        const auto [one, two, three]  = Cornered(points, triangle);
        const Point middle            = {.X = (one.X + two.X + three.X) / 3,
                                         .Y = (one.Y + two.Y + three.Y) / 3};
        const bool  turned            = Orient(one, two, three) == Orientation::Counterclockwise;
        wrong                        += turned && Within(loops, middle) ? 0U : 1U;
      }
      return wrong;
    }

    /**
     *  How many sides are used twice the same way round.
     *
     *  @param  mesh The fill.
     *  @return      The count.
     */
    auto Repeated(const Mesh &mesh) -> std::size_t
    {
      std::set<std::pair<std::uint32_t, std::uint32_t>> sides;
      std::size_t                                       again = 0;
      for (const auto &[one, two, three] : mesh.Triangles)
      {
        for (const auto &side : {std::pair{one, two}, std::pair{two, three}, std::pair{three, one}})
        {
          again += sides.insert(side).second ? 0U : 1U;
        }
      }
      return again;
    }

    /**
     *  How far the triangles' areas sum from the region's, exactly: its sign.
     *
     *  @param  loops The region's loops.
     *  @param  mesh  The fill.
     *  @return       The sign of the difference.
     */
    auto Unsummed(std::span<const std::vector<Point>> loops, const Mesh &mesh) -> int
    {
      const std::vector<Point> points = Flat(loops);
      Expansion::Number        sum    = Expansion::Of(0.0);
      for (const auto &triangle : mesh.Triangles)
      {
        const auto [one, two, three] = Cornered(points, triangle);
        sum                          = Expansion::Sum(sum, Twice(one, two, three));
      }
      Expansion::Number region = Enclosed(loops.front());
      for (const std::vector<Point> &hole : loops.subspan(1))
      {
        region = Expansion::Difference(region, Enclosed(hole));
      }
      return Expansion::Sign(Expansion::Difference(sum, region));
    }

    /**
     *  A fill checked: every triangle counterclockwise, its centroid within the region, no side
     *  used twice the same way, and the areas summing to the region's exactly.
     *
     *  @param loops The region's loops.
     *  @param mesh  The fill.
     */
    void Checked(std::span<const std::vector<Point>> loops, const Mesh &mesh)
    {
      CHECK(mesh.Whole);
      CHECK(Misplaced(loops, mesh) == 0);
      CHECK(Repeated(mesh)         == 0);
      CHECK(Unsummed(loops, mesh)  == 0);
    }

    /**
     *  How many sides two triangles share are not locally Delaunay: the far corner of one inside
     *  the other's circle.
     *
     *  @param  loops The region's loops.
     *  @param  mesh  The fill.
     *  @return       How many.
     */
    auto Undelaunay(std::span<const std::vector<Point>> loops, const Mesh &mesh) -> std::size_t
    {
      const std::vector<Point> points = Flat(loops);
      std::size_t              bad    = 0;
      for (const auto &first : mesh.Triangles)
      {
        const auto [one, two, three] = Cornered(points, first);
        for (const auto &second : mesh.Triangles)
        {
          for (const std::uint32_t far : second)
          {
            const bool opposite =
                std::ranges::count(first, far) == 0
                && std::ranges::count_if(second, [&first](std::uint32_t corner) -> bool
                                         { return std::ranges::count(first, corner) > 0; })
                       == 2;
            const int inside =
                Exact::InCircle({.First = one, .Second = two, .Third = three}, points[far]);
            bad += opposite && inside > 0 ? 1 : 0;
          }
        }
      }
      return bad;
    }

    /**
     *  A star of points about a center, its radius wavering, counterclockwise or not.
     *
     *  @param  center The center.
     *  @param  count  How many points.
     *  @param  shape  Its radius, how much it wavers, and the seed.
     *  @return        The points.
     */
    auto Star(Point center, std::size_t count, std::array<double, 3> shape) -> std::vector<Point>
    {
      std::mt19937_64                        draw(static_cast<std::uint64_t>(shape[2]));
      std::uniform_real_distribution<double> waver(1.0 - shape[1], 1.0);
      std::vector<Point>                     points;
      for (std::size_t index = 0; index < count; ++index)
      {
        const double angle =
            2.0 * std::numbers::pi * static_cast<double>(index) / static_cast<double>(count);
        const double radius = shape[0] * waver(draw);
        points.push_back({.X = center.X + (radius * std::cos(angle)),
                          .Y = center.Y + (radius * std::sin(angle))});
      }
      return points;
    }

    /**
     *  A square, counterclockwise.
     *
     *  @param  low  Its lower left corner.
     *  @param  side Its side.
     *  @return      Its corners.
     */
    auto Square(Point low, double side) -> std::vector<Point>
    {
      return {low, low + Vector{.X = side, .Y = 0}, low + Vector{.X = side, .Y = side},
              low + Vector{.X = 0, .Y = side}};
    }
  } // namespace

  TEST_CASE("triangulate: a square, a concave polygon and its reverse")
  {
    const std::vector<std::vector<Point>> square = {Square({}, 2)};
    const Mesh                            two    = Triangulate(square, Fill::Ears);
    CHECK(two.Triangles.size() == 2);
    Checked(square, two);

    std::vector<Point> comb = {{.X = 0, .Y = 0}, {.X = 10, .Y = 0}, {.X = 10, .Y = 5}};
    for (int tooth = 4; tooth >= 0; --tooth)
    {
      comb.push_back({.X = (2.0 * tooth) + 1.5, .Y = 5});
      comb.push_back({.X = (2.0 * tooth) + 1, .Y = 1});
      comb.push_back({.X = (2.0 * tooth) + 0.5, .Y = 5});
    }
    comb.push_back({.X = 0, .Y = 5});
    const std::vector<std::vector<Point>> combed = {comb};
    Checked(combed, Triangulate(combed, Fill::Ears));
    std::ranges::reverse(comb);
    const std::vector<std::vector<Point>> reversed = {comb};
    Checked(reversed, Triangulate(reversed, Fill::Ears));
    Checked(reversed, Triangulate(reversed, Fill::Delaunay));
    CHECK(Undelaunay(reversed, Triangulate(reversed, Fill::Delaunay)) == 0);
  }

  TEST_CASE("triangulate: holes bridged")
  {
    for (const Vector offset : {Vector{}, Vector{.X = 5e6, .Y = 5e6}, Vector{.X = -1e8, .Y = 1e8}})
    {
      const Point                           base  = Point{} + offset;
      const std::vector<std::vector<Point>> holed = {
          Square(base, 10), Square(base + Vector{.X = 2, .Y = 2}, 2),
          Square(base + Vector{.X = 6, .Y = 2}, 2),
          Star(base + Vector{.X = 5, .Y = 7}, 9, {1.5, 0.5, 3})};
      for (const Fill fill : {Fill::Ears, Fill::Delaunay})
      {
        const Mesh mesh = Triangulate(holed, fill);
        Checked(holed, mesh);
        CHECK(mesh.Triangles.size() == 4 + 4 + 4 + 9 + (2 * 3) - 2);
      }
      CHECK(Undelaunay(holed, Triangulate(holed, Fill::Delaunay)) == 0);
    }
  }

  TEST_CASE("triangulate: large regions along the Z curve")
  {
    for (const Vector offset : {Vector{}, Vector{.X = 5e6, .Y = -5e6}})
    {
      std::vector<std::vector<Point>> region = {Star(Point{} + offset, 4000, {100, 0.3, 1})};
      for (int hole = 0; hole < 6; ++hole)
      {
        const double angle = hole * std::numbers::pi / 3;
        region.push_back(Star(Point{.X = 40 * std::cos(angle), .Y = 40 * std::sin(angle)} + offset,
                              50 + (10 * static_cast<std::size_t>(hole)), {12, 0.4, 7.0 + hole}));
      }
      const Mesh  ears   = Triangulate(region, Fill::Ears);
      std::size_t points = 0;
      for (const std::vector<Point> &loop : region)
      {
        points += loop.size();
      }
      Checked(region, ears);
      CHECK(ears.Triangles.size() == points + (2 * (region.size() - 1)) - 2);
      const Mesh delaunay = Triangulate(region, Fill::Delaunay);
      Checked(region, delaunay);
      CHECK(delaunay.Triangles.size() == ears.Triangles.size());
    }
  }

  TEST_CASE("triangulate: corners that are none, and loops that touch")
  {
    const std::vector<std::vector<Point>> repeated = {{{.X = 0, .Y = 0},
                                                       {.X = 1, .Y = 0},
                                                       {.X = 2, .Y = 0},
                                                       {.X = 2, .Y = 0},
                                                       {.X = 2, .Y = 2},
                                                       {.X = 0, .Y = 2},
                                                       {.X = 0, .Y = 1},
                                                       {.X = 0, .Y = 0}}};
    const Mesh                            square   = Triangulate(repeated, Fill::Ears);
    CHECK(square.Triangles.size() == 2);
    Checked(repeated, square);

    const std::vector<std::vector<Point>> flat = {
        {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = 2}}};
    CHECK(Triangulate(flat, Fill::Ears).Triangles.empty());
    CHECK(Triangulate(flat, Fill::Ears).Whole);
    CHECK(Triangulate({}, Fill::Ears).Triangles.empty());
    const std::vector<std::vector<Point>> empty = {{}};
    CHECK(Triangulate(empty, Fill::Ears).Triangles.empty());

    const std::vector<std::vector<Point>> flatHole = {
        Square({}, 4), {{.X = 1, .Y = 1}, {.X = 2, .Y = 2}, {.X = 3, .Y = 3}}};
    CHECK(Triangulate(flatHole, Fill::Ears).Triangles.size() == 2);

    const std::vector<std::vector<Point>> touching = {
        Square({}, 4), {{.X = 0, .Y = 2}, {.X = 1, .Y = 1}, {.X = 2, .Y = 2}, {.X = 1, .Y = 3}}};
    Checked(touching, Triangulate(touching, Fill::Ears));

    const std::vector<std::vector<Point>> notch = {
        {{.X = 0, .Y = 0}, {.X = 4, .Y = 0}, {.X = 4, .Y = 4}, {.X = 0, .Y = 4}, {.X = 1, .Y = 2}},
        {{.X = 1, .Y = 2}, {.X = 2, .Y = 1}, {.X = 3, .Y = 2}, {.X = 2, .Y = 3}}};
    Checked(notch, Triangulate(notch, Fill::Delaunay));

    const std::vector<std::vector<Point>> apart = {Square({}, 4), Square({.X = -3, .Y = 1}, 1)};
    const Mesh                            lost  = Triangulate(apart, Fill::Ears);
    CHECK_FALSE(lost.Whole);
    CHECK(lost.Triangles.size() == 2);

    const std::vector<std::vector<Point>> crossed = {{{.X = 0, .Y = 0},
                                                      {.X = 4, .Y =  0},
                                                      {.X = 4, .Y =  4},
                                                      {.X = 3, .Y = -1},
                                                      {.X = 2, .Y = 4}}};
    CHECK_FALSE(Triangulate(crossed, Fill::Ears).Whole);
  }
} // namespace Geoxy
