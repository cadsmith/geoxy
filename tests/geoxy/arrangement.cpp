/***************************************************************************************************
 *  @file      tests/geoxy/arrangement.cpp
 *  @brief     Tests of the loop around a picked point: the rooms of a plan, a column as an island,
 *             squares nested, doors and short walls closed by the gap, loose lines pruned, picks on
 *             a wall and outside everything, at every scale.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/arrangement.hpp"

#include "loops.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/error.hpp"
#include "geoxy/line.hpp"
#include "geoxy/loop.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/region.hpp"
#include "geoxy/tolerance.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <limits>
#include <numbers>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A tolerance fine enough for every case drawn about unit size.
     */
    constexpr Closing ArrangedFine = {.Within = {.Length = 1e-9, .Radians = 1e-9}, .Gap = 0.0};

    /**
     *  The places each case is drawn at: the origin, survey coordinates and far past them.
     */
    constexpr std::array<double, 3> ArrangedOffsets = {0.0, 5e6, 1e8};

    /**
     *  A line between two points, both moved by a place.
     *
     *  @param  offset The place.
     *  @param  start  Where it starts, from the place.
     *  @param  end    Where it ends, from the place.
     *  @return        The line.
     */
    auto Wall(double offset, Point start, Point end) -> Curve
    {
      return Line{.Start = {.X = start.X + offset, .Y = start.Y + offset},
                  .End   = {.X = end.X + offset, .Y = end.Y + offset}};
    }

    /**
     *  The loops of a region's boundary, or 99 where it has none.
     *
     *  @param  region The region, or why there is none.
     *  @return        The count.
     */
    auto BoundaryLoops(const std::expected<Region, Error> &region) -> std::size_t
    {
      return region ? region->Loops.size() : 99;
    }

    /**
     *  A region's area, or not a number where there is none.
     *
     *  @param  region The region, or why there is none.
     *  @return        The area.
     */
    auto AreaOf(const std::expected<Region, Error> &region) -> double
    {
      return region ? Area(*region, Accuracy{}) : std::numeric_limits<double>::quiet_NaN();
    }

    /**
     *  A square's sides, each a curve of its own, counterclockwise.
     *
     *  @param  least Its least corner.
     *  @param  side  Its side.
     *  @return       The curves.
     */
    auto Sides(Point least, double side) -> std::vector<Curve>
    {
      return Tests::Square(least, side).Edges;
    }
  } // namespace

  TEST_CASE("arrangement: the rooms of a plan, a pick on a wall and one outside, at every scale")
  {
    for (const double offset : ArrangedOffsets)
    {
      const std::vector<Curve> walls = {
          Wall(offset, {.X = 0, .Y = 0}, {.X = 10, .Y = 0}),
          Wall(offset, {.X = 10, .Y = 0}, {.X = 10, .Y = 6}),
          Wall(offset, {.X = 10, .Y = 6}, {.X = 0, .Y = 6}),
          Wall(offset, {.X = 0, .Y = 6}, {.X = 0, .Y = 0}),
          Wall(offset, {.X = 4, .Y = 0}, {.X = 4, .Y = 6}),
          Wall(offset, {.X = 10, .Y = 3}, {.X = 4, .Y = 3}),
      };
      const Tolerance tolerance = Tolerance::Of(
          {.Min = {.X = offset, .Y = offset}, .Max = {.X = offset + 10, .Y = offset + 6}});
      const Arrangement plan{walls, {.Within = tolerance, .Gap = 0}};
      const auto        pick = [&plan, offset](double x, double y) -> std::expected<Region, Error>
      {
        return Boundary(plan, {.X = offset + x, .Y = offset + y});
      };
      CHECK(AreaOf(pick(2, 3))         == 24.0);
      CHECK(AreaOf(pick(7, 1))         == 18.0);
      CHECK(AreaOf(pick(7, 5))         == 18.0);
      CHECK(BoundaryLoops(pick(2, 3))  == 1);
      CHECK(BoundaryLoops(pick(-1, 1)) == 0);
      CHECK(BoundaryLoops(pick(4, 1))  == 0);
      CHECK(BoundaryLoops(pick(11, 7)) == 0);
      CHECK(pick(7, 5)->Loops.front().Edges.size() == 4);
    }
  }

  TEST_CASE("arrangement: a column is an island of its room, and squares nest, in any order")
  {
    std::vector<Curve> room = Sides({}, 10);
    room.emplace_back(Circle({.X = 5, .Y = 5}, 1));
    const Arrangement columned{room, ArrangedFine};
    CHECK(BoundaryLoops(Boundary(columned, {.X = 2, .Y = 2})) == 2);
    CHECK(AreaOf(Boundary(columned, {.X = 2, .Y = 2}))
          == doctest::Approx(100 - std::numbers::pi).epsilon(1e-15));
    CHECK(AreaOf(Boundary(columned, {.X = 5, .Y = 5}))
          == doctest::Approx(std::numbers::pi).epsilon(1e-15));
    CHECK(Boundary(columned, {.X = 5, .Y = 5})->Loops.front().Edges.size() == 1);
    const std::array<std::vector<Curve>, 3> squares = {Sides({}, 10), Sides({.X = 2, .Y = 2}, 6),
                                                       Sides({.X = 4, .Y = 4}, 2)};
    for (const auto &order :
         {std::array<std::size_t, 3>{0, 1, 2}, std::array<std::size_t, 3>{2, 1, 0}})
    {
      std::vector<Curve> curves;
      for (const std::size_t index : order)
      {
        curves.insert(curves.end(), squares.at(index).begin(), squares.at(index).end());
      }
      const Arrangement nested{curves, ArrangedFine};
      CHECK(AreaOf(Boundary(nested, {.X = 1, .Y = 1}))        == 64.0);
      CHECK(BoundaryLoops(Boundary(nested, {.X = 1, .Y = 1})) == 2);
      CHECK(AreaOf(Boundary(nested, {.X = 3, .Y = 3}))        == 32.0);
      CHECK(AreaOf(Boundary(nested, {.X = 5, .Y = 5}))        == 4.0);
    }
  }

  TEST_CASE("arrangement: a door's gap and a short wall closed by the gap, not past it")
  {
    const std::vector<Curve> door = {
        Wall(0, {.X = 0, .Y = 0}, {.X = 4, .Y = 0}),
        Wall(0, {.X = 4.5, .Y = 0}, {.X = 10, .Y = 0}),
        Wall(0, {.X = 10, .Y = 0}, {.X = 10, .Y = 10}),
        Wall(0, {.X = 10, .Y = 10}, {.X = 0, .Y = 10}),
        Wall(0, {.X = 0, .Y = 10}, {.X = 0, .Y = 0}),
        Wall(0, {.X = 4.2, .Y = 8}, {.X = 4.2, .Y = 7}),
    };
    const Arrangement shut{door, {.Within = ArrangedFine.Within, .Gap = 1}};
    CHECK(AreaOf(Boundary(shut, {.X = 5, .Y = 5})) == 100.0);
    const Arrangement open{door, {.Within = ArrangedFine.Within, .Gap = 0.1}};
    CHECK(BoundaryLoops(Boundary(open, {.X = 5, .Y = 5})) == 0);
    for (const bool back : {false, true})
    {
      std::vector<Curve> rooms = Sides({}, 10);
      rooms.push_back(back ? Wall(0, {.X = 5, .Y = 9.7}, {.X = 5, .Y = 0})
                           : Wall(0, {.X = 5, .Y = 0}, {.X = 5, .Y = 9.7}));
      rooms.push_back(Wall(0, {.X = 7, .Y = 9.6}, {.X = 7, .Y = 9.2}));
      const Arrangement walled{rooms, {.Within = ArrangedFine.Within, .Gap = 0.5}};
      CHECK(AreaOf(Boundary(walled, {.X = 2, .Y = 5})) == 50.0);
      CHECK(AreaOf(Boundary(walled, {.X = 8, .Y = 5})) == 50.0);
    }
    const std::vector<Curve> three = {
        Wall(0, {.X = 0, .Y = 0}, {.X = 4, .Y = 0}),
        Wall(0, {.X = 4.3, .Y = 0}, {.X = 10, .Y = 0}),
        Wall(0, {.X = 4.1, .Y = -2}, {.X = 4.1, .Y = -0.2}),
        Wall(0, {.X = 10, .Y = 0}, {.X = 10, .Y = 10}),
        Wall(0, {.X = 10, .Y = 10}, {.X = 0, .Y = 10}),
        Wall(0, {.X = 0, .Y = 10}, {.X = 0, .Y = 0}),
    };
    const Arrangement crowded{three, {.Within = ArrangedFine.Within, .Gap = 0.5}};
    CHECK(AreaOf(Boundary(crowded, {.X = 5, .Y = 5})) == doctest::Approx(100.03).epsilon(1e-15));
    const Arrangement closed{Sides({}, 10), {.Within = ArrangedFine.Within, .Gap = 0.5}};
    CHECK(AreaOf(Boundary(closed, {.X = 5, .Y = 5})) == 100.0);
    std::vector<Curve> corner = Sides({}, 10);
    corner.push_back(Wall(0, {.X = 5, .Y = 5}, {.X = 9.8, .Y = 9.7}));
    const Arrangement cornered{corner, {.Within = ArrangedFine.Within, .Gap = 0.5}};
    CHECK(BoundaryLoops(Boundary(cornered, {.X = 9, .Y = 2})) == 1);
  }

  TEST_CASE("arrangement: loose lines are pruned, arcs bound rooms, and polylines are cut")
  {
    std::vector<Curve> room = Sides({}, 10);
    room.push_back(Wall(0, {.X = 2, .Y = 2}, {.X = 3, .Y = 3}));
    room.push_back(Wall(0, {.X = 10, .Y = 5}, {.X = 12, .Y = 5}));
    const Arrangement loose{room, ArrangedFine};
    CHECK(AreaOf(Boundary(loose, {.X = 5, .Y = 5}))        == 100.0);
    CHECK(BoundaryLoops(Boundary(loose, {.X = 5, .Y = 5})) == 1);
    const std::vector<Curve> swing = {
        Wall(0, {.X = 0, .Y = 0}, {.X = 1, .Y = 0}),
        Arc{.Center = {}, .Radius = 1, .Start = Degrees(0), .Sweep = Degrees(90)},
        Wall(0, {.X = 0, .Y = 1}, {.X = 0, .Y = 0}),
    };
    const Arrangement door{swing, ArrangedFine};
    CHECK(AreaOf(Boundary(door, {.X = 0.3, .Y = 0.3}))
          == doctest::Approx(std::numbers::pi / 4).epsilon(1e-15));
    const Polyline    outline = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                              {.Position = {.X = 4, .Y = 0}, .Bulge = 1},
                                              {.Position = {.X = 4, .Y = 4}},
                                              {.Position = {.X = 0, .Y = 4}}},
                                 .Closed   = true};
    const Arrangement bulged{
        std::vector<Curve>{outline, Wall(0, {.X = 2, .Y = 0}, {.X = 2, .Y = 4})}, ArrangedFine};
    CHECK(AreaOf(Boundary(bulged, {.X = 3, .Y = 2}))
          == doctest::Approx(8 + (2 * std::numbers::pi)).epsilon(1e-15));
    CHECK(AreaOf(Boundary(bulged, {.X = 1, .Y = 2})) == 8.0);
  }

  TEST_CASE("arrangement: a pick says why it cannot answer, and an arrangement moved from has none")
  {
    const Arrangement oval{std::vector<Curve>{Ellipse{.Map = {}, .Start = {}, .Sweep = Turn}},
                           ArrangedFine};
    CHECK(Boundary(oval, {}).error() == Error::Degenerate);
    const Arrangement nan{
        std::vector<Curve>{Wall(std::numeric_limits<double>::quiet_NaN(), {}, {.X = 1, .Y = 0})},
        {.Within = ArrangedFine.Within, .Gap = 1}};
    CHECK(Boundary(nan, {}).error() == Error::NotFinite);
    Arrangement room{Sides({}, 10), ArrangedFine};
    CHECK(Boundary(room, {.X = std::numeric_limits<double>::infinity(), .Y = 0}).error()
          == Error::NotFinite);
    CHECK(Boundary(room, {.X = 0, .Y = std::numeric_limits<double>::quiet_NaN()}).error()
          == Error::NotFinite);
    const Arrangement copy  = room;
    const Arrangement taken = std::move(room);
    CHECK(AreaOf(Boundary(taken, {.X = 5, .Y = 5}))       == 100.0);
    CHECK(AreaOf(Boundary(copy, {.X = 5, .Y = 5}))        == 100.0);
    CHECK(BoundaryLoops(Boundary(room, {.X = 5, .Y = 5})) == 0);
    const Arrangement none{std::vector<Curve>{}, ArrangedFine};
    CHECK(BoundaryLoops(Boundary(none, {})) == 0);
    const Arrangement disc{std::vector<Curve>{Circle({}, 1)}, ArrangedFine};
    CHECK(BoundaryLoops(Boundary(disc, {.X = 0.9, .Y = 0.9})) == 0);
    CHECK(BoundaryLoops(Boundary(disc, {.X = 0.5, .Y = 0.5})) == 1);
  }

  TEST_CASE("arrangement: faces no plane's are mended as a boolean's are")
  {
    const auto ring = [](const std::vector<Point> &points) -> std::vector<Curve>
    {
      std::vector<Curve> lines;
      lines.reserve(points.size());
      for (std::size_t index = 0; index < points.size(); ++index)
      {
        lines.emplace_back(
            Line{.Start = points[index], .End = points[(index + 1) % points.size()]});
      }
      return lines;
    };
    const auto both = [](std::vector<Curve>        first,
                         const std::vector<Curve> &second) -> std::vector<Curve>
    {
      first.insert(first.end(), second.begin(), second.end());
      return first;
    };
    /**
     *  Two quadrilaterals whose crossings merged make two lines cross twice, mended by the split; a
     *  pick inside their intersection found nothing, and finds it.
     */
    const Arrangement quads{
        both(ring({{.X = 0x1.417d8d205e980p+11, .Y = 0x1.eaec0a23c08aap+12},
                   {.X = 0x1.29369128024fdp+12, .Y = 0x1.3e32cfc222974p+12},
                   {.X = 0x1.8b4b3b57ef5c7p+12, .Y = 0x1.8bb107f2339a4p+12},
                   {.X = 0x1.02d370c01c58ap+12, .Y = 0x1.1c352129e8c6dp+13}}),
             ring({{.X = 0x1.449b6a8da6934p+12, .Y = 0x1.70a80f0b025e9p+7},
                   {.X = 0x1.441e9cfa9de44p+11, .Y = 0x1.29c795f2b5056p+13},
                   {.X = 0x1.376a7e9af9df5p+13, .Y = 0x1.97cf9e535f909p+11},
                   {.X = 0x1.2ddfdc03b4f26p+11, .Y = 0x1.ce9574f1eceb6p+12}})),
        {.Within = {.Length = 9.344619977064311e-07, .Radians = 1e-10}, .Gap = 0.0}};
    const auto quad = Boundary(quads, {.X = 4476, .Y = 7012});
    REQUIRE(quad.has_value());
    CHECK(std::abs(Area(*quad, Accuracy{0}) - 4433981.640161665) < 1e-8);
    /**
     *  Two triangles whose faces stay no plane's once split, mended by the two pieces between two
     *  merged vertices made one; the face is their exact intersection.
     */
    const Arrangement triangles{
        both(ring({{.X = -0x1.3000617251c3fp+0, .Y = 0x1.a7640e41caee3p-1},
                   {.X = -0x1.3d189a0e23de9p+0, .Y = 0x1.3ea6ff1bba511p-3},
                   {.X = 0x1.dd6a3dbfeebbdp-1, .Y = -0x1.0fb25c4a70e91p+0}}),
             ring({{.X = -0x1.5e8b290a7ae33p+0, .Y = -0x1.2677ce00a54c4p-1},
                   {.X = 0x1.f0f743bf601e2p-2, .Y = 0x1.5f6c4675a0b16p-2},
                   {.X = -0x1.0aaf3acc1f4d9p+0, .Y = 0x1.547ab6e816ceep+0}})),
        {.Within = {.Length = 1e-2, .Radians = 1e-9}, .Gap = 0.0}};
    const auto triangle = Boundary(triangles, {.X = -0.89, .Y = 0.19});
    REQUIRE(triangle.has_value());
    CHECK(std::abs(Area(*triangle, Accuracy{0}) - 0.444582325792361) < 1e-15);
    /**
     *  Lenses of arcs whose faces are no plane's: one however mended, the first kept; one a plane's
     *  traced banded.
     */
    const double      thin  = 0x1.e9168e1c0714p-7;
    const double      thick = 0x1.ef23c9a2e69d8p-6;
    const double      round = 0x1.15783cb4e3668p-3;
    const Arrangement lenses{
        std::vector<Curve>{
            Polyline{
                .Vertices = {{.Position = {.X = 0x1.4c08p+11, .Y = -0x1.4244p+13}, .Bulge = thin},
                             {.Position = {.X = 0x1.2016p+13, .Y = -0x1.55ccp+12}, .Bulge = thick}},
                .Closed   = true},
            Polyline{
                .Vertices = {{.Position = {.X = 0x1.86ap+9, .Y = -0x1.6472p+13}, .Bulge = round},
                             {.Position = {.X = 0x1.81bep+13, .Y = 0x1.388p+8}, .Bulge = round}},
                .Closed   = true}},
        {.Within = {.Length = 100, .Radians = 1e-9}, .Gap = 0.0}};
    CHECK(Boundary(lenses, {.X = 6000, .Y = -6000}).has_value());
    const double      cap = 0x1.a856f3b3e3828p-4;
    const Arrangement caps{
        std::vector<Curve>{
            Polyline{
                .Vertices = {{.Position = {.X = 0x1.2fd7c83739794p+22, .Y = 0x1.d682b76bb622cp+21},
                              .Bulge    = cap},
                             {.Position = {.X = 0x1.7718a44a24eeap+22, .Y = 0x1.2fd7c83739794p+22},
                              .Bulge    = cap}},
                .Closed   = true},
            Polyline{.Vertices =
                         {{.Position = {.X = 0x1.7718a44a24eeap+22, .Y = 0x1.2fd7c83739794p+22},
                           .Bulge    = cap},
                          {.Position = {.X = 0x1.41386df77b4e3p+22, .Y = 0x1.5c3537d0b5a19p+21}}},
                     .Closed = true}},
        {.Within = {.Length = 1000, .Radians = 1e-9}, .Gap = 0.0}};
    CHECK(Boundary(caps, {.X = 5.5e6, .Y = 4.5e6}).has_value());
  }
} // namespace Geoxy
