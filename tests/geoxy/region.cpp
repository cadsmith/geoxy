/***************************************************************************************************
 *  @file      tests/geoxy/region.cpp
 *  @brief     Tests of regions: where a point lies, areas, and the booleans of regions of lines and
 *             arcs against exact areas and the invariants every boolean keeps, at every scale.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/region.hpp"

#include "draws.hpp"
#include "loops.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/error.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/loop.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <numbers>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A tolerance fine enough for every case drawn about unit size.
     */
    constexpr Tolerance RegionFine = {.Length = 1e-9, .Radians = 1e-9};

    /**
     *  The places each case is drawn at: the origin, survey coordinates and far past them.
     */
    constexpr std::array<double, 3> Offsets = {0.0, 5e6, 1e8};

    /**
     *  A region's area, or not a number where there is none.
     *
     *  @param  region The region, or why there is none.
     *  @return        The area.
     */
    auto Measured(const std::expected<Region, Error> &region) -> double
    {
      return region ? Area(*region, Accuracy{}) : std::numeric_limits<double>::quiet_NaN();
    }

    /**
     *  How many loops a region has, or none where there is no region.
     *
     *  @param  region The region, or why there is none.
     *  @return        The count; 99 for none.
     */
    auto Loops(const std::expected<Region, Error> &region) -> std::size_t
    {
      return region ? region->Loops.size() : 99;
    }

    /**
     *  Whether no two edges of a region cross or share a stretch, nor any edge itself.
     *
     *  @param  region    The region.
     *  @param  tolerance The tolerance its curves are met within.
     *  @return           Whether none does.
     */
    auto Simple(const Region &region, const Tolerance &tolerance) -> bool
    {
      std::vector<Curve> edges;
      for (const Loop &loop : region.Loops)
      {
        edges.insert(edges.end(), loop.Edges.begin(), loop.Edges.end());
      }
      for (std::size_t first = 0; first < edges.size(); ++first)
      {
        for (std::size_t second = first + 1; second < edges.size(); ++second)
        {
          const Intersections meeting = Intersect(edges[first], edges[second], tolerance);
          const bool          crossed = std::ranges::any_of(meeting.Hits, [](const Hit &hit) -> bool
                                                            { return hit.Kind == Contact::Cross; });
          if (crossed || !meeting.Overlaps.empty())
          {
            return false;
          }
        }
      }
      return true;
    }

    /**
     *  The area to which two areas of a size agree: a few steps of a double at the largest.
     *
     *  @param  size The largest area.
     *  @return      The slack.
     */
    auto Slack(double size) -> double
    {
      return 1e-12 * std::max(1.0, std::abs(size));
    }

    /**
     *  Every boolean of two regions checked: the areas add up, either order gives the same, no
     *  result crosses itself, and every point of a grid lies where the operands say.
     *
     *  @param one       One region.
     *  @param two       The other.
     *  @param tolerance The tolerance.
     */
    void Invariants(const Region &one, const Region &two, const Tolerance &tolerance)
    {
      const auto joined = Union(one, two, tolerance);
      const auto shared = Intersection(one, two, tolerance);
      const auto taken  = Difference(one, two, tolerance);
      const auto other  = Difference(two, one, tolerance);
      const auto turned = Union(two, one, tolerance);
      const auto both   = Intersection(two, one, tolerance);
      const bool given  = joined && shared && taken && other && turned && both;
      REQUIRE(given);
      const double                first  = Measured(Filled(one, tolerance));
      const double                second = Measured(Filled(two, tolerance));
      const double                united = Measured(joined);
      const double                common = Measured(shared);
      const std::array<double, 5> misses = {
          std::abs(united + common - first - second),  std::abs(Measured(taken) + common - first),
          std::abs(Measured(other) + common - second), std::abs(Measured(turned) - united),
          std::abs(Measured(both) - common),
      };
      CHECK(std::ranges::max(misses) <= Slack(first + second));
      const bool simple =
          Simple(*joined, tolerance) && Simple(*shared, tolerance) && Simple(*taken, tolerance);
      CHECK(simple);
    }

    /**
     *  Where every point of a grid over a box lies against the booleans of two regions, as the
     *  regions say: a point inside both, inside one, inside neither.
     *
     *  @param one       One region.
     *  @param two       The other.
     *  @param tolerance The tolerance.
     */
    void Agrees(const Region &one, const Region &two, const Tolerance &tolerance)
    {
      const Region joined = *Union(one, two, tolerance);
      const Region shared = *Intersection(one, two, tolerance);
      const Region taken  = *Difference(one, two, tolerance);
      const auto   placed = [](bool inside) -> Location
      {
        return inside ? Location::Inside : Location::Outside;
      };
      for (int index = 0; index < 23 * 23; ++index)
      {
        const int      row    = index / 23;
        const int      column = index % 23;
        const Point    point  = {.X = -0.4 + (column * 0.2093), .Y = -0.4 + (row * 0.2131)};
        const Location first  = Locate(one, point, tolerance);
        const Location second = Locate(two, point, tolerance);
        const bool     inOne  = first == Location::Inside;
        const bool     inTwo  = second == Location::Inside;
        const bool     clear  = first != Location::On && second != Location::On;
        CHECK((!clear
               || (Locate(joined, point, tolerance) == placed(inOne || inTwo)
                   && Locate(shared, point, tolerance) == placed(inOne && inTwo)
                   && Locate(taken, point, tolerance) == placed(inOne && !inTwo))));
      }
    }

    /**
     *  A random polygon: points about a center at random radii, in order of angle.
     *
     *  @param  draws  The draws.
     *  @param  center The center.
     *  @param  count  How many points.
     *  @return        The loop.
     */
    auto Star(Draws &draws, Point center, int count) -> Loop
    {
      std::vector<Point> points;
      for (int index = 0; index < count; ++index)
      {
        const double angle  = 2 * std::numbers::pi * index / count;
        const double radius = draws.Near(1.2, 0.8);
        points.push_back({.X = center.X + (radius * std::cos(angle)),
                          .Y = center.Y + (radius * std::sin(angle))});
      }
      return Tests::Polygon(points);
    }

    /**
     *  A smooth star of many sides waving a few times about a center, as a polyline its sides
     *  bulged alike: its sides run in long stretches that rise or fall along both axes.
     *
     *  @param  center Its center.
     *  @param  count  How many sides; its radius waves three to seven times round it, by the count.
     *  @param  bulge  Each side's bulge.
     *  @return        Its loop.
     */
    auto Wavy(Point center, int count, double bulge) -> Loop
    {
      Polyline polyline = {.Closed = true};
      for (int index = 0; index < count; ++index)
      {
        const double angle  = 2 * std::numbers::pi * index / count;
        const double radius = 1 + (0.25 * std::sin(((count % 5) + 3) * angle));
        polyline.Vertices.push_back({.Position = {.X = center.X + (radius * std::cos(angle)),
                                                  .Y = center.Y + (radius * std::sin(angle))},
                                     .Bulge    = bulge});
      }
      return {.Edges = {polyline}};
    }

    /**
     *  A random rounded rectangle: its corners quarter circles.
     *
     *  @param  draws  The draws.
     *  @param  center The center.
     *  @return        The loop.
     */
    auto Pill(Draws &draws, Point center) -> Loop
    {
      const double width  = draws.Near(1.5, 1.0);
      const double height = draws.Near(1.5, 1.0);
      const double radius = draws.Near(0.4, 0.3);
      const double left   = center.X - width;
      const double right  = center.X + width;
      const double bottom = center.Y - height;
      const double top    = center.Y + height;
      return {
          .Edges = {
              Line{.Start = {.X = left + radius, .Y = bottom},
                   .End   = {.X = right - radius, .Y = bottom}},
              Arc{.Center = {.X = right - radius, .Y = bottom + radius},
                  .Radius = radius,
                  .Start  = Degrees(-90),
                  .Sweep  = Degrees(90)},
              Line{.Start = {.X = right, .Y = bottom + radius},
                   .End   = {.X = right, .Y = top - radius}},
              Arc{.Center = {.X = right - radius, .Y = top - radius},
                  .Radius = radius,
                  .Start  = Degrees(0),
                  .Sweep  = Degrees(90)},
              Line{.Start = {.X = right - radius, .Y = top}, .End = {.X = left + radius, .Y = top}},
              Arc{.Center = {.X = left + radius, .Y = top - radius},
                  .Radius = radius,
                  .Start  = Degrees(90),
                  .Sweep  = Degrees(90)},
              Line{.Start = {.X = left, .Y = top - radius},
                   .End   = {.X = left, .Y = bottom + radius}},
              Arc{.Center = {.X = left + radius, .Y = bottom + radius},
                  .Radius = radius,
                  .Start  = Degrees(180),
                  .Sweep  = Degrees(90)},
          }};
    }
  } // namespace

  TEST_CASE("region: a point lies inside, outside, on, or in a hole")
  {
    const Region holed = {
        .Loops = {Tests::Square({}, 4), Tests::Reversed(Tests::Square({.X = 1, .Y = 1}, 2))}};
    CHECK(Locate(holed, {.X = 0.5, .Y = 0.5}, RegionFine) == Location::Inside);
    CHECK(Locate(holed, {.X = 2, .Y = 2}, RegionFine)     == Location::Outside);
    CHECK(Locate(holed, {.X = 5, .Y = 2}, RegionFine)     == Location::Outside);
    CHECK(Locate(holed, {.X = 1, .Y = 2}, RegionFine)     == Location::On);
    CHECK(Locate(holed, {.X = 0, .Y = 2}, RegionFine)     == Location::On);
    CHECK(Locate(holed, {.X = std::numeric_limits<double>::quiet_NaN(), .Y = 2}, RegionFine)
          == Location::On);
    CHECK(Locate({}, {}, RegionFine) == Location::Outside);
    CHECK(Area(holed, Accuracy{})    == 12.0);
    CHECK(holed
          == Region{.Loops = {Tests::Square({}, 4),
                              Tests::Reversed(Tests::Square({.X = 1, .Y = 1}, 2))}});
    CHECK_FALSE(holed == Region{.Loops = {Tests::Square({}, 4)}});
    CHECK(*Filled({.Loops = {Tests::Square({}, 2)}}, RegionFine)
          == Region{.Loops = {Tests::Square({}, 2)}});
    CHECK(Area(Region{}, Accuracy{}) == 0.0);
  }

  TEST_CASE("region: two squares that overlap, at every scale")
  {
    for (const double offset : Offsets)
    {
      const Tolerance tolerance = Tolerance::Of(
          {.Min = {.X = offset, .Y = offset}, .Max = {.X = offset + 3, .Y = offset + 3}});
      const Vector shift = {.X = offset, .Y = offset};
      const Region first = Tests::Moved(Region{.Loops = {Tests::Square({}, 2)}}, shift);
      const Region second =
          Tests::Moved(Region{.Loops = {Tests::Square({.X = 1, .Y = 1}, 2)}}, shift);
      CHECK(Measured(Union(first, second, tolerance))        == 7.0);
      CHECK(Measured(Intersection(first, second, tolerance)) == 1.0);
      CHECK(Measured(Difference(first, second, tolerance))   == 3.0);
      CHECK(Loops(Union(first, second, tolerance))           == 1);
      CHECK(Union(first, second, tolerance)->Loops.front().Edges.size() == 8);
      Invariants(first, second, tolerance);
    }
  }

  TEST_CASE("region: squares apart, nested, the same, and sharing an edge or a corner")
  {
    const Region square = {.Loops = {Tests::Square({}, 2)}};
    const Region apart  = {.Loops = {Tests::Square({.X = 5, .Y = 0}, 2)}};
    CHECK(Measured(Union(square, apart, RegionFine))      == 8.0);
    CHECK(Loops(Union(square, apart, RegionFine))         == 2);
    CHECK(Loops(Intersection(square, apart, RegionFine))  == 0);
    CHECK(Measured(Difference(square, apart, RegionFine)) == 4.0);
    const Region inner = {.Loops = {Tests::Square({.X = 0.5, .Y = 0.5}, 1)}};
    CHECK(Measured(Union(square, inner, RegionFine))        == 4.0);
    CHECK(Loops(Union(square, inner, RegionFine))           == 1);
    CHECK(Measured(Intersection(square, inner, RegionFine)) == 1.0);
    CHECK(Measured(Difference(square, inner, RegionFine))   == 3.0);
    CHECK(Loops(Difference(square, inner, RegionFine))      == 2);
    CHECK(Loops(Difference(inner, square, RegionFine))      == 0);
    Invariants(square, inner, RegionFine);
    CHECK(Measured(Union(square, square, RegionFine))        == 4.0);
    CHECK(Loops(Union(square, square, RegionFine))           == 1);
    CHECK(Measured(Intersection(square, square, RegionFine)) == 4.0);
    CHECK(Loops(Difference(square, square, RegionFine))      == 0);
    const Region back = {.Loops = {Tests::Reversed(Tests::Square({}, 2))}};
    CHECK(Measured(Union(square, back, RegionFine))   == 4.0);
    CHECK(Loops(Difference(square, back, RegionFine)) == 0);
    const Region beside = {.Loops = {Tests::Square({.X = 2, .Y = 0}, 2)}};
    CHECK(Measured(Union(square, beside, RegionFine))     == 8.0);
    CHECK(Loops(Union(square, beside, RegionFine))        == 1);
    CHECK(Loops(Intersection(square, beside, RegionFine)) == 0);
    const Region step = {.Loops = {Tests::Square({.X = 2, .Y = 1}, 2)}};
    CHECK(Measured(Union(square, step, RegionFine)) == 8.0);
    CHECK(Loops(Union(square, step, RegionFine))    == 1);
    Invariants(square, step, RegionFine);
    const Region corner = {.Loops = {Tests::Square({.X = 2, .Y = 2}, 2)}};
    CHECK(Measured(Union(square, corner, RegionFine))     == 8.0);
    CHECK(Loops(Intersection(square, corner, RegionFine)) == 0);
    Invariants(square, corner, RegionFine);
  }

  TEST_CASE("region: discs crossing, cocircular, and a disc and its half")
  {
    const double lens   = (2 * std::acos(0.5)) - (0.5 * std::numbers::sqrt3);
    const Region disc   = {.Loops = {Tests::Disc({}, 1)}};
    const Region beside = {.Loops = {Tests::Disc({.X = 1, .Y = 0}, 1)}};
    CHECK(Measured(Intersection(disc, beside, RegionFine)) == doctest::Approx(lens).epsilon(1e-14));
    CHECK(Measured(Union(disc, beside, RegionFine))
          == doctest::Approx((2 * std::numbers::pi) - lens).epsilon(1e-14));
    CHECK(Intersection(disc, beside, RegionFine)->Loops.front().Edges.size() == 2);
    Invariants(disc, beside, RegionFine);
    Agrees(disc, beside, RegionFine);
    const Region halves = {
        .Loops = {
            {.Edges = {
                 Arc{.Center = {}, .Radius = 1, .Start = Degrees(0), .Sweep = Degrees(180)},
                 Arc{.Center = {}, .Radius = 1, .Start = Degrees(180), .Sweep = Degrees(180)}}}}};
    CHECK(Measured(Union(disc, halves, RegionFine))
          == doctest::Approx(std::numbers::pi).epsilon(1e-15));
    CHECK(Measured(Intersection(disc, halves, RegionFine))
          == doctest::Approx(std::numbers::pi).epsilon(1e-15));
    CHECK(Loops(Difference(disc, halves, RegionFine)) == 0);
    const Region half = {
        .Loops = {
            {.Edges = {Arc{.Center = {}, .Radius = 1, .Start = Degrees(0), .Sweep = Degrees(180)},
                       Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 1, .Y = 0}}}}}};
    CHECK(Measured(Intersection(disc, half, RegionFine))
          == doctest::Approx(std::numbers::pi / 2).epsilon(1e-15));
    CHECK(Measured(Difference(disc, half, RegionFine))
          == doctest::Approx(std::numbers::pi / 2).epsilon(1e-15));
    Invariants(disc, half, RegionFine);
    const Region square = {.Loops = {Tests::Square({.X = -2, .Y = -2}, 4)}};
    const auto   ring   = Difference(square, disc, RegionFine);
    CHECK(Measured(ring) == doctest::Approx(16 - std::numbers::pi).epsilon(1e-15));
    REQUIRE(Loops(ring)  == 2);
    CHECK(ring->Loops.back().Edges.size() == 1);
    const Region left = {.Loops = {Tests::Square({.X = -1.5, .Y = -0.5}, 1)}};
    const auto   seam = Union(disc, left, RegionFine);
    REQUIRE(Loops(seam) == 1);
    CHECK(seam->Loops.front().Edges.size() == 4);
    Invariants(disc, left, RegionFine);
  }

  TEST_CASE("region: discs far from the origin, and a disc turned far round before it starts")
  {
    const double lens = (2 * std::acos(0.5)) - (0.5 * std::numbers::sqrt3);
    for (const double offset : Offsets)
    {
      const Vector shift = {.X = offset, .Y = offset};
      const Region disc  = Tests::Moved(Region{.Loops = {Tests::Disc({}, 1)}}, shift);
      const Region beside =
          Tests::Moved(Region{.Loops = {Tests::Disc({.X = 1, .Y = 0}, 1)}}, shift);
      CHECK(Measured(Intersection(disc, beside, RegionFine))
            == doctest::Approx(lens).epsilon(1e-6));
      CHECK(Loops(Union(disc, beside, RegionFine)) == 1);
    }
    const Region disc   = {.Loops = {Tests::Disc({}, 1)}};
    const Region turned = {
        .Loops = {
            {.Edges = {Arc{.Center = {}, .Radius = 1, .Start = Degrees(1e7), .Sweep = Turn}}}}};
    const Region square = {.Loops = {Tests::Square({}, 2)}};
    CHECK(Measured(Union(turned, square, RegionFine))
          == doctest::Approx(Measured(Union(disc, square, RegionFine))).epsilon(1e-9));
    CHECK(Measured(Intersection(turned, square, RegionFine))
          == doctest::Approx(std::numbers::pi / 4).epsilon(1e-9));
  }

  TEST_CASE("region: a disc run clockwise, and arcs leaving a point along one tangent")
  {
    const Region square    = {.Loops = {Tests::Square({}, 2)}};
    const Region disc      = {.Loops = {Tests::Disc({.X = 2, .Y = 1}, 0.5)}};
    const Region clockwise = {.Loops = {Tests::Reversed(Tests::Disc({.X = 2, .Y = 1}, 0.5))}};
    CHECK(Measured(Union(square, clockwise, RegionFine))
          == doctest::Approx(Measured(Union(square, disc, RegionFine))).epsilon(1e-15));
    const Region inverted = {.Loops = {Tests::Disc({.X = 2, .Y = 1}, -0.5)}};
    CHECK(Measured(Union(square, inverted, RegionFine))
          == doctest::Approx(Measured(Union(square, disc, RegionFine))).epsilon(1e-15));
    CHECK(Measured(Intersection(square, clockwise, RegionFine))
          == doctest::Approx(std::numbers::pi / 8).epsilon(1e-15));
    const auto quarter = [](double radius) -> Arc
    {
      return {.Center = {.X = 0, .Y = radius},
              .Radius = radius,
              .Start  = Degrees(-90),
              .Sweep  = Degrees(90)};
    };
    const Loop inner = {.Edges = {quarter(1),
                                  Line{.Start = {.X = 1, .Y = 1}, .End = {.X = 2, .Y = 2}},
                                  Tests::Reversed(quarter(2))}};
    const Loop outer = {.Edges = {quarter(2),
                                  Line{.Start = {.X = 2, .Y = 2}, .End = {.X = 3, .Y = 3}},
                                  Tests::Reversed(quarter(3))}};
    const auto fan   = Filled({.Loops = {inner, outer}}, RegionFine);
    CHECK(Measured(fan)
          == doctest::Approx(-Area(inner, Accuracy{}) - Area(outer, Accuracy{})).epsilon(1e-14));
    CHECK(Loops(fan) == 1);
  }

  TEST_CASE("region: what one region fills, its loops crossing, winding twice, or running back")
  {
    const Region bow = {
        .Loops = {Tests::Polygon(
            {{.X = 0, .Y = 0}, {.X = 2, .Y = 2}, {.X = 2, .Y = 0}, {.X = 0, .Y = 2}})}};
    CHECK(Area(bow, Accuracy{})             == 0.0);
    CHECK(Measured(Filled(bow, RegionFine)) == 2.0);
    CHECK(Loops(Filled(bow, RegionFine))    == 2);
    const Region twice = {
        .Loops = {{.Edges = {Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = Degrees(360)},
                             Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = Degrees(360)}}}}};
    CHECK(Measured(Filled(twice, RegionFine)) == doctest::Approx(std::numbers::pi).epsilon(1e-15));
    const Region back = {.Loops = {Tests::Reversed(Tests::Square({}, 2))}};
    CHECK(Measured(Filled(back, RegionFine)) == 4.0);
    const Region crescent = {
        .Loops = {{.Edges = {Arc{.Center = {.X = 0, .Y = 2},
                                 .Radius = 2,
                                 .Start  = Degrees(-90),
                                 .Sweep  = Degrees(180)},
                             Line{.Start = {.X = 0, .Y = 4}, .End = {.X = 0, .Y = 2}},
                             Arc{.Center = {.X = 0, .Y = 1},
                                 .Radius = 1,
                                 .Start  = Degrees(90),
                                 .Sweep  = Degrees(-180)}}}}};
    CHECK(Measured(Filled(crescent, RegionFine))
          == doctest::Approx(1.5 * std::numbers::pi).epsilon(1e-15));
    CHECK(Filled(crescent, RegionFine)->Loops.front().Edges.size() == 3);
    const Region slot  = {.Loops = {Tests::Slot({}, 3, 1)}};
    const Region shelf = {.Loops = {Tests::Rectangle({.X = 1, .Y = -1}, {.X = 6, .Y = -0.5})}};
    const double cut   = 1 + (((std::numbers::pi / 3) - (0.5 * std::sqrt(0.75))) / 2);
    CHECK(Measured(Union(slot, shelf, RegionFine))
          == doctest::Approx(6 + std::numbers::pi + 2.5 - cut).epsilon(1e-15));
    Invariants(slot, shelf, RegionFine);
    const Region dropped = {.Loops = {Tests::Square({}, 10)}};
    const Region spikes  = {
        .Loops = {
            Tests::Polygon({{.X = 5, .Y = -2}, {.X = 5 + 1e-11, .Y = 3}, {.X = 4, .Y = -2}}),
            Tests::Polygon({{.X = 5.5, .Y = -3}, {.X = 5 + 2e-11, .Y = 3}, {.X = 5.3, .Y = -3}})}};
    Invariants(dropped, spikes, RegionFine);
    /**
     *  A square with a spike out and back along one line, whose tip one edge alone reaches; and a
     *  loop running back and forth along one line, its boxes all one.
     */
    const Region spiked = {.Loops = {Tests::Polygon({{.X = 0, .Y = 0},
                                                     {.X = 2, .Y = 0},
                                                     {.X = 2, .Y = 1},
                                                     {.X = 3, .Y = 1},
                                                     {.X = 2, .Y = 1},
                                                     {.X = 2, .Y = 2},
                                                     {.X = 0, .Y = 2}})}};
    CHECK(Measured(Filled(spiked, RegionFine))         == 4.0);
    CHECK(Measured(Union(spiked, dropped, RegionFine)) == 100.0);
    const Region shuttle = {.Loops = {Tests::Polygon({{.X = 0, .Y = 0},
                                                      {.X = 1, .Y = 1},
                                                      {.X = 0, .Y = 0},
                                                      {.X = 1, .Y = 1},
                                                      {.X = 0, .Y = 0},
                                                      {.X = 1, .Y = 1}})}};
    CHECK(Measured(Filled(shuttle, RegionFine)) == 0.0);
  }

  TEST_CASE("region: bulged polygons met at a coarse tolerance, two cuts at one place made one")
  {
    /**
     *  At a fifth of their size, a mark lies where a curve's end is, its place there the end's, yet
     *  farther than the tolerance from the end's vertex: the two are made one at that place.
     */
    const Polyline three = {
        .Vertices = {{.Position = {.X = 1.125, .Y = -0.0625}, .Bulge = 0.875},
                     {.Position = {.X = 0.125, .Y = 1.28125}, .Bulge = -0.125},
                     {.Position = {.X = -0.75, .Y = -0.859375}, .Bulge = -0.125}},
        .Closed   = true};
    const Polyline four = {
        .Vertices = {{.Position = {.X = 1.6875, .Y = 0.296875}, .Bulge = -1},
                     {.Position = {.X = -0.078125, .Y = 1.265625}, .Bulge = -0.625},
                     {.Position = {.X = -0.703125, .Y = -0.484375}, .Bulge = 0.625},
                     {.Position = {.X = 0.296875, .Y = -1.046875}, .Bulge = -0.25}},
        .Closed   = true};
    const Region    triangle = {.Loops = {{.Edges = {three}}}};
    const Region    quad     = {.Loops = {{.Edges = {four}}}};
    const Tolerance coarse   = {.Length = 0.2, .Radians = 1e-9};
    const double    one      = Measured(Filled(triangle, coarse));
    const double    other    = Measured(Filled(quad, coarse));
    const double    united   = Measured(Union(triangle, quad, coarse));
    CHECK(united >= std::max(one, other));
    CHECK(united <= one + other);
    CHECK(Measured(Intersection(triangle, quad, coarse)) >= 0.0);
    CHECK(Measured(Difference(triangle, quad, coarse)) >= 0.0);
    CHECK(Measured(Difference(quad, triangle, coarse)) >= 0.0);
  }

  TEST_CASE("region: a line's pieces that come round to where they began are made one")
  {
    /**
     *  At a tolerance of 1, a side cut twice has its pieces' ends merged into two vertices, the
     *  later piece coming back to where the earlier began: a loop of two pieces of one line, which
     *  threw where it was taken for a circle's. Its faces are no plane's, split or not; the two
     *  pieces, drawn as one segment, are made one, their tallies cancelling, and the area is 309.1,
     *  within 4 of the 313.0 at a tolerance of 1e-9, where it was 77.7.
     */
    const Region jagged = {.Loops = {Tests::Polygon({{.X = 49.75, .Y = 154},
                                                     {.X =  21.5, .Y =     71},
                                                     {.X =    21, .Y =  70.75},
                                                     {.X = 45.25, .Y = 156.75},
                                                     {.X = 15.25, .Y =  53.75},
                                                     {.X = 23.25, .Y =  84.75},
                                                     {.X = 19.25, .Y =   71.5},
                                                     {.X = 39.75, .Y =    152},
                                                     {.X = 29, .Y = 114}})}};
    const auto   filled = Filled(jagged, {.Length = 1, .Radians = 1e-9});
    REQUIRE(filled.has_value());
    CHECK(std::ranges::none_of(filled->Loops,
                               [](const Loop &loop) -> bool { return loop.Edges.size() == 2; }));
    CHECK(Measured(filled) >= Measured(Filled(jagged, RegionFine)) - 4.0);
    CHECK(Measured(filled) <= Measured(Filled(jagged, RegionFine)));
  }

  TEST_CASE("region: a loop of one line's pieces out and back is two sides")
  {
    /**
     *  Random polygons about 1,000 at a tolerance of 1e-2, cut down: a side's pieces merged so that
     *  the later comes back to where the earlier began, in faces a plane's, the loop between them
     *  kept. It is two sides, the second the first turned round, as no circle's is.
     */
    const auto ring = [](const std::vector<Point> &points) -> Region
    {
      Polyline polyline = {.Closed = true};
      for (const Point point : points)
      {
        polyline.Vertices.push_back({.Position = point});
      }
      return {.Loops = {{.Edges = {polyline}}}};
    };
    const Region first  = ring({{.X = 0x1.f3b924660a970p+9, .Y = 0x1.f48561bc104ddp+9},
                                {.X = 0x1.f3ad748ed8004p+9, .Y = 0x1.f47208f494fe7p+9},
                                {.X = 0x1.f393913f611d6p+9, .Y = 0x1.f493d59dd8d23p+9},
                                {.X = 0x1.f39a7eda5a004p+9, .Y = 0x1.f48894baea111p+9},
                                {.X = 0x1.f3a6fb9a5c47ap+9, .Y = 0x1.f46a9e3aca951p+9},
                                {.X = 0x1.f392890683492p+9, .Y = 0x1.f4817259c08c5p+9},
                                {.X = 0x1.f3b650155ba9ap+9, .Y = 0x1.f4560932cfdaap+9},
                                {.X = 0x1.f3510a43c9df7p+9, .Y = 0x1.f40b5129e8974p+9},
                                {.X = 0x1.f48bea9ca4040p+9, .Y = 0x1.f421cb9356f08p+9}});
    const Region second = ring({{.X = 0x1.f490af371c565p+9, .Y = 0x1.f423bb368dfbap+9},
                                {.X = 0x1.f44d1a56ff2e3p+9, .Y = 0x1.f440d6eea2c5bp+9},
                                {.X = 0x1.f362503be9ebcp+9, .Y = 0x1.f4a54edbdc75fp+9},
                                {.X = 0x1.f38bf27a79690p+9, .Y = 0x1.f491910fed87cp+9},
                                {.X = 0x1.f3a70709169d3p+9, .Y = 0x1.f48592c29cb10p+9},
                                {.X = 0x1.f36c4759b5271p+9, .Y = 0x1.f49cdfa742f71p+9},
                                {.X = 0x1.f39fc22192330p+9, .Y = 0x1.f48750102d301p+9},
                                {.X = 0x1.f389ce230b557p+9, .Y = 0x1.f48f49ab66086p+9},
                                {.X = 0x1.f394bdda34337p+9, .Y = 0x1.f482102ac8bc9p+9},
                                {.X = 0x1.f38da76e13e79p+9, .Y = 0x1.f482c33a86bfap+9},
                                {.X = 0x1.f3a1fd90bb4ecp+9, .Y = 0x1.f47c5b5d37944p+9},
                                {.X = 0x1.f3774d5fd87c5p+9, .Y = 0x1.f48773575f9dfp+9},
                                {.X = 0x1.f3d2360e1bed8p+9, .Y = 0x1.f40d3a23e0a74p+9}});
    const auto   joined = Union(first, second, {.Length = 1e-2, .Radians = 1e-9});
    REQUIRE(joined.has_value());
    const auto paired = std::ranges::find_if(joined->Loops, [](const Loop &loop) -> bool
                                             { return loop.Edges.size() == 2; });
    REQUIRE(paired != joined->Loops.end());
    const Line &out  = std::get<Line>(paired->Edges.front());
    const Line &back = std::get<Line>(paired->Edges.back());
    CHECK(out.Start == back.End);
    CHECK(out.End   == back.Start);
  }

  TEST_CASE("region: polylines with bulges, and a region with a hole met by a square")
  {
    const double   bulge   = std::tan(std::numbers::pi / 8);
    const Polyline rounded = {.Vertices = {{.Position = {.X = 0.5, .Y = 0}},
                                           {.Position = {.X = 1.5, .Y = 0}, .Bulge = bulge},
                                           {.Position = {.X = 2, .Y = 0.5}},
                                           {.Position = {.X = 2, .Y = 1.5}, .Bulge = bulge},
                                           {.Position = {.X = 1.5, .Y = 2}},
                                           {.Position = {.X = 0.5, .Y = 2}, .Bulge = bulge},
                                           {.Position = {.X = 0, .Y = 1.5}},
                                           {.Position = {.X = 0, .Y = 0.5}, .Bulge = bulge}},
                              .Closed   = true};
    const Region   first   = {.Loops = {{.Edges = {rounded}}}};
    const Region   second  = {.Loops = {Tests::Square({.X = 1, .Y = 1}, 2)}};
    CHECK(Measured(Filled(first, RegionFine))
          == doctest::Approx(3 + (std::numbers::pi / 4)).epsilon(1e-14));
    /**
     *  A bulged side between a vertex and its own repeat is a line of no length, as Nth makes it.
     */
    Polyline repeated = rounded;
    repeated.Vertices.insert(repeated.Vertices.begin() + 2, rounded.Vertices[1]);
    CHECK(Measured(Filled({.Loops = {{.Edges = {repeated}}}}, RegionFine))
          == Measured(Filled(first, RegionFine)));
    Invariants(first, second, RegionFine);
    Agrees(first, second, RegionFine);
    const Region holed = {
        .Loops = {Tests::Square({}, 4), Tests::Reversed(Tests::Square({.X = 1, .Y = 1}, 2))}};
    const Region patch = {.Loops = {Tests::Square({.X = 2, .Y = 2}, 3)}};
    CHECK(Measured(Union(holed, patch, RegionFine)) == 18.0);
    Invariants(holed, patch, RegionFine);
    Agrees(holed, patch, RegionFine);
  }

  TEST_CASE("region: a circle and a bar combined, and a shared edge")
  {
    const Polyline circle = {.Vertices = {{.Position = {.X = 0, .Y = 1}, .Bulge = 1},
                                          {.Position = {.X = 10, .Y = 1}, .Bulge = 1}},
                             .Closed   = true};
    const Region   round  = {.Loops = {{.Edges = {circle}}}};
    const Region   bar    = {.Loops = {Tests::Rectangle({.X = 3, .Y = -10}, {.X = 6, .Y = 10})}};
    const double   strip =
        std::sqrt(24.0) + (2 * std::sqrt(21.0)) + (25 * (std::asin(0.2) + std::asin(0.4)));
    CHECK(Measured(Intersection(round, bar, RegionFine)) == doctest::Approx(strip).epsilon(1e-14));
    CHECK(Measured(Intersection(round, bar, RegionFine))
          == doctest::Approx(29.386000046924).epsilon(1e-12));
    CHECK(Measured(Union(round, bar, RegionFine))
          == doctest::Approx(109.15381629282).epsilon(1e-12));
    CHECK(Measured(Difference(round, bar, RegionFine))
          == doctest::Approx(29.336980664548 + 19.816835628274).epsilon(1e-12));
    CHECK(Loops(Difference(round, bar, RegionFine)) == 2);
    Invariants(round, bar, RegionFine);
    const double   corner = -0.414214;
    const Polyline small  = {.Vertices = {{.Position = {.X = -0.105, .Y = 0.235}},
                                          {.Position = {.X = -0.095, .Y = 0.235}},
                                          {.Position = {.X = -0.095, .Y = 0}, .Bulge = -1},
                                          {.Position = {.X = -0.105, .Y = 0}}},
                             .Closed   = true};
    const Polyline slab   = {.Vertices = {{.Position = {.X = -0.25, .Y = 0.235}, .Bulge = corner},
                                          {.Position = {.X = -0.255, .Y = 0.24}},
                                          {.Position = {.X = -0.255, .Y = 0.29}, .Bulge = corner},
                                          {.Position = {.X = -0.25, .Y = 0.295}},
                                          {.Position = {.X = 0.25, .Y = 0.295}, .Bulge = corner},
                                          {.Position = {.X = 0.255, .Y = 0.29}},
                                          {.Position = {.X = 0.255, .Y = 0.24}, .Bulge = corner},
                                          {.Position = {.X = 0.25, .Y = 0.235}}},
                             .Closed   = true};
    const Region   one    = {.Loops = {{.Edges = {small}}}};
    const Region   two    = {.Loops = {{.Edges = {slab}}}};
    CHECK(Measured(Union(one, two, RegionFine))
          == doctest::Approx(0.032967809756574).epsilon(1e-12));
    CHECK(Loops(Union(one, two, RegionFine)) == 1);
    CHECK(Measured(Difference(one, two, RegionFine))
          == doctest::Approx(0.0023892699081699).epsilon(1e-12));
    CHECK(Measured(Difference(two, one, RegionFine))
          == doctest::Approx(0.030578539848405).epsilon(1e-12));
    CHECK(Loops(Intersection(one, two, RegionFine)) == 0);
    Invariants(one, two, RegionFine);
  }

  TEST_CASE("region: random polygons and rounded rectangles keep every invariant")
  {
    Draws draws{2026};
    for (int round = 0; round < Rounds(48); ++round)
    {
      const Point  center = {.X = draws.Near(1.5, 0.6), .Y = draws.Near(1.5, 0.6)};
      const Region first  = {.Loops = {Star(draws, {.X = 1.5, .Y = 1.5}, 3 + (round % 9))}};
      const Region second = {
          .Loops = {round % 2 == 0 ? Star(draws, center, 4 + (round % 7)) : Pill(draws, center)}};
      Invariants(first, second, RegionFine);
      Agrees(first, second, RegionFine);
    }
  }

  TEST_CASE("region: smooth stars of many sides, their sides in long runs, keep every invariant")
  {
    for (const int count : {40, 96})
    {
      for (const double bulge : {0.0, 0.15})
      {
        const Region first  = {.Loops = {Wavy({.X = 1.5, .Y = 1.5}, count, bulge)}};
        const Region second = {.Loops = {Wavy({.X = 1.8, .Y = 1.7}, count + 7, bulge)}};
        Invariants(first, second, RegionFine);
        Agrees(first, second, RegionFine);
      }
    }
    /**
     *  A star against itself shifted along its own sides, and against a smaller one inside it.
     */
    const Region star  = {.Loops = {Wavy({.X = 1.5, .Y = 1.5}, 64, 0.0)}};
    const Region small = {.Loops = {Wavy({.X = 1.5, .Y = 1.5}, 48, 0.0)}};
    Invariants(star, star, RegionFine);
    Invariants(star, small, RegionFine);
    Agrees(star, small, RegionFine);
  }

  TEST_CASE("region: booleans say why they cannot answer")
  {
    const Region square = {.Loops = {Tests::Square({}, 2)}};
    const Region open   = {.Loops = {{.Edges = {Line{.Start = {}, .End = {.X = 1, .Y = 0}}}}}};
    CHECK(Union(square, open, RegionFine).error() == Error::Degenerate);
    const Region nan = {
        .Loops = {{.Edges = {Line{.Start = {.X = std::numeric_limits<double>::quiet_NaN(), .Y = 0},
                                  .End   = {.X = 1, .Y = 0}}}}}};
    CHECK(Union(square, nan, RegionFine).error() == Error::NotFinite);
    const Region rayed = {.Loops = {{.Edges = {Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}}}}};
    CHECK(Filled(rayed, RegionFine).error() == Error::Degenerate);
    const Region endless = {.Loops = {{.Edges = {Xline{.Through = {}, .Way = {.X = 1, .Y = 0}}}}}};
    CHECK(Filled(endless, RegionFine).error() == Error::Degenerate);
    const Region gaps = {.Loops = {open.Loops.front(), Tests::Square({}, 2)}};
    CHECK(Filled(gaps, RegionFine).error() == Error::Degenerate);
    const Region oval = {.Loops = {{.Edges = {Ellipse{.Map = {}, .Start = {}, .Sweep = Turn}}}}};
    CHECK(Filled(oval, RegionFine).error() == Error::Degenerate);
    const Region coil = {
        .Loops = {{.Edges = {Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = Degrees(720)}}}}};
    CHECK(Filled(coil, RegionFine).error() == Error::Degenerate);
    const Region both = {
        .Loops = {{.Edges = {Ellipse{.Map = {}, .Start = {}, .Sweep = Turn}}}, nan.Loops.front()}};
    CHECK(Filled(both, RegionFine).error() == Error::Degenerate);
    const Polyline broken = {.Vertices = {{.Position = {}},
                                          {.Position = {.X = 1, .Y = 0},
                                           .Bulge    = std::numeric_limits<double>::infinity()}},
                             .Closed   = true};
    CHECK(Filled({.Loops = {{.Edges = {broken}}}}, RegionFine).error() == Error::NotFinite);
    CHECK(Loops(Filled({}, RegionFine)) == 0);
  }

  TEST_CASE("region: crossings merged by distance never cross two lines twice")
  {
    /**
     *  Two random polygons' three almost concurrent lines, their crossings 0.90e-6, 0.89e-6 and
     *  1.65e-6 apart against a tolerance of 0.93e-6, cut down to two quadrilaterals. Two of the
     *  crossings merged, and the faces traced about them gave an empty intersection; the areas are
     *  the exact ones, worked out in rationals, to a few steps of a double.
     */
    const auto ring = [](const std::vector<Point> &points) -> Region
    {
      Polyline polyline = {.Closed = true};
      for (const Point point : points)
      {
        polyline.Vertices.push_back({.Position = point});
      }
      return {.Loops = {{.Edges = {polyline}}}};
    };
    const Region    first  = ring({{.X = 0x1.417d8d205e980p+11, .Y = 0x1.eaec0a23c08aap+12},
                                   {.X = 0x1.29369128024fdp+12, .Y = 0x1.3e32cfc222974p+12},
                                   {.X = 0x1.8b4b3b57ef5c7p+12, .Y = 0x1.8bb107f2339a4p+12},
                                   {.X = 0x1.02d370c01c58ap+12, .Y = 0x1.1c352129e8c6dp+13}});
    const Region    second = ring({{.X = 0x1.449b6a8da6934p+12, .Y = 0x1.70a80f0b025e9p+7},
                                   {.X = 0x1.441e9cfa9de44p+11, .Y = 0x1.29c795f2b5056p+13},
                                   {.X = 0x1.376a7e9af9df5p+13, .Y = 0x1.97cf9e535f909p+11},
                                   {.X = 0x1.2ddfdc03b4f26p+11, .Y = 0x1.ce9574f1eceb6p+12}});
    const Tolerance wide   = {.Length = 9.344619977064311e-07, .Radians = 1e-10};
    const auto      met    = Intersection(first, second, wide);
    const auto      joined = Union(first, second, wide);
    REQUIRE(met.has_value());
    REQUIRE(joined.has_value());
    CHECK(std::abs(Area(*met, Accuracy{0}) - 4433981.640161665) < 1e-8);
    CHECK(std::abs(Area(*joined, Accuracy{0}) - 12586554.781327523) < 1e-8);
  }

  TEST_CASE("region: arcs that cross within twice the length are ordered by where they part")
  {
    /**
     *  Bulged polygons of about 0.003 met within 1e-5, arcs crossing twice within the length. Their
     *  faces are no plane's whether or not the clusters are split, and the union was empty. Traced
     *  banded, the halves of two arcs that stay within twice the length until they cross ordered as
     *  they part, they are a plane's; the union lies within 4e-9 of the 2.78365e-6 at a tolerance a
     *  tenth as wide, against the tolerance times the perimeter, 8.7e-8.
     */
    const double bulge = 0x1.3333333333333p-3;
    const auto   ring  = [](const std::vector<Vertex> &vertices) -> Region
    {
      return {.Loops = {{.Edges = {Polyline{.Vertices = vertices, .Closed = true}}}}};
    };
    const Region first = ring(
        {{.Position = {.X = 0x1.312cfffeccccdp+22, .Y = 0x1.312d000041893p+22}, .Bulge = bulge},
         {.Position = {.X = 0x1.312cffff2f1aap+22, .Y = 0x1.312cffff58106p+22}, .Bulge = bulge},
         {.Position = {.X = 0x1.312d00007ae14p+22, .Y = 0x1.312cffff645a2p+22}, .Bulge = bulge},
         {.Position = {.X = 0x1.312d0000bc6a8p+22, .Y = 0x1.312d00005a1cbp+22}, .Bulge = bulge},
         {.Position = {.X = 0x1.312d000004189p+22, .Y = 0x1.312d00013f7cfp+22}, .Bulge = bulge}});
    const Region second = ring(
        {{.Position = {.X = 0x1.312cfffeccccdp+22, .Y = 0x1.312d000041893p+22}, .Bulge = bulge},
         {.Position = {.X = 0x1.312cffff2f1aap+22, .Y = 0x1.312cffff58106p+22}, .Bulge = bulge},
         {.Position = {.X = 0x1.312cfffe41893p+22, .Y = 0x1.312d000039581p+22}}});
    const auto joined = Union(first, second, {.Length = 0x1.4f8b588e368f1p-17, .Radians = 1e-9});
    REQUIRE(joined.has_value());
    CHECK(std::abs(Area(*joined, Accuracy{0}) - 2.7836470399661931e-06) < 4e-9);
  }

  TEST_CASE("region: pieces between two merged vertices that cross there are made one")
  {
    const auto ring = [](const std::vector<Vertex> &vertices) -> Region
    {
      return {.Loops = {{.Edges = {Polyline{.Vertices = vertices, .Closed = true}}}}};
    };
    /**
     *  Two triangles from random ones of 16 vertices at 1e-2 of their size, faces no plane's
     *  whether or not split, gave an empty intersection. With the two pieces between two merged
     *  vertices made one, the intersection is the exact one and the union that of the triangles
     *  drawn through the merged vertices, both worked out in rationals.
     */
    const Region first =
        ring({{.Position = {.X = -0x1.3000617251c3fp+0, .Y = 0x1.a7640e41caee3p-1}},
              {.Position = {.X = -0x1.3d189a0e23de9p+0, .Y = 0x1.3ea6ff1bba511p-3}},
              {.Position = {.X = 0x1.dd6a3dbfeebbdp-1, .Y = -0x1.0fb25c4a70e91p+0}}});
    const Region second =
        ring({{.Position = {.X = -0x1.5e8b290a7ae33p+0, .Y = -0x1.2677ce00a54c4p-1}},
              {.Position = {.X = 0x1.f0f743bf601e2p-2, .Y = 0x1.5f6c4675a0b16p-2}},
              {.Position = {.X = -0x1.0aaf3acc1f4d9p+0, .Y = 0x1.547ab6e816ceep+0}}});
    const Tolerance wide   = {.Length = 1e-2, .Radians = 1e-9};
    const auto      met    = Intersection(first, second, wide);
    const auto      joined = Union(first, second, wide);
    REQUIRE(met.has_value());
    REQUIRE(joined.has_value());
    CHECK(std::abs(Area(*met, Accuracy{0}) - 0.444582325792361) < 1e-15);
    CHECK(std::abs(Area(*joined, Accuracy{0}) - 1.9296039990211726) < 1e-15);
    /**
     *  Two triangles of arcs at 1e-2 of their size, one passing 40 from the other's corner, gave an
     *  empty intersection; both now lie within the tolerance times the boundary, 1.5e7, of the
     *  areas at a millionth of it, the union to the bit.
     */
    const double bulge = 0x1.cea51c5fc53ecp-3;
    const Region arced = ring(
        {{.Position = {.X = -0x1.1385a1656fc39p+14, .Y = -0x1.9aecb3a80a638p+11}, .Bulge = bulge},
         {.Position = {.X = 0x1.46f7d1fbec5e3p+13, .Y = -0x1.7e038f4d9268dp+13}, .Bulge = bulge},
         {.Position = {.X = 0x1.91311da837cefp+11, .Y = 0x1.1c2f4def97181p+13}, .Bulge = bulge}});
    const Region other = ring(
        {{.Position = {.X = 0x1.80332862d750dp+12, .Y = 0x1.4de1d815d09bcp+12}, .Bulge = bulge},
         {.Position = {.X = -0x1.e45dfddfb65b0p+12, .Y = 0x1.5e7681d5bee7ep+14}, .Bulge = bulge},
         {.Position = {.X = -0x1.9aa04808ab65ap+13, .Y = -0x1.2a24c39918b47p+12}, .Bulge = bulge}});
    const Tolerance coarse  = {.Length = 100, .Radians = 1e-9};
    const auto      crossed = Intersection(arced, other, coarse);
    const auto      covered = Union(arced, other, coarse);
    REQUIRE(crossed.has_value());
    REQUIRE(covered.has_value());
    CHECK(std::abs(Area(*crossed, Accuracy{0}) - 163247372.399934) < 1e4);
    CHECK(std::abs(Area(*covered, Accuracy{0}) - 619728864.954121) < 1e-3);
    /**
     *  Two lenses of two arcs each at 1e-2 of their size, still no plane's banded or not, pieces
     *  made one or not: the union is made from the first arrangement.
     */
    const Region lens = ring(
        {{.Position = {.X = 0x1.4c08p+11, .Y = -0x1.4244p+13}, .Bulge = 0x1.e9168e1c0714p-7},
         {.Position = {.X = 0x1.2016p+13, .Y = -0x1.55ccp+12}, .Bulge = 0x1.ef23c9a2e69d8p-6}});
    const Region bowed =
        ring({{.Position = {.X = 0x1.86ap+9, .Y = -0x1.6472p+13}, .Bulge = 0x1.15783cb4e3668p-3},
              {.Position = {.X = 0x1.81bep+13, .Y = 0x1.388p+8}, .Bulge = 0x1.15783cb4e3668p-3}});
    CHECK(Union(lens, bowed, coarse).has_value());
  }

  TEST_CASE("region: faces no plane's once split are traced and mended after ferrying")
  {
    /**
     *  A triangle against a triangle and a pentagon, each side its own line, met at 1e-2 of their
     *  size: the windings ferried disagree, as do those of the split arrangement, so the faces are
     *  traced and mended, as they were before windings were ferried, to the bit.
     */
    const auto ring = [](const std::vector<Point> &points) -> Loop
    {
      Loop loop = {};
      for (std::size_t index = 0; index < points.size(); ++index)
      {
        loop.Edges.emplace_back(
            Line{.Start = points[index], .End = points[(index + 1) % points.size()]});
      }
      return loop;
    };
    const Region first = {
        .Loops = {ring({{.X = 0x1.9a9f4b0bb5ce4p+2, .Y = -0x1.3c6c1b77322e5p+2},
                        {.X = 0x1.4de1b99cf745ep+2, .Y = 0x1.c8ed3236fb61dp+1},
                        {.X = -0x1.c075a9c776e39p+0, .Y = 0x1.7d5d8b6f37e8bp+2}})}};
    const Region second = {
        .Loops = {ring({{.X = -0x1.7283d5f8db38p-4, .Y = 0x1.75f8f25c9ec54p+1},
                        {.X = -0x1.5d071680918d1p+1, .Y = 0x1.e4fd85e778c1p-2},
                        {.X = -0x1.e07f24f37d879p+1, .Y = -0x1.69874772307afp+1}}),
                  ring({{.X = -0x1.11c3959b6015p+1, .Y = 0x1.75cdbeb49e183p+1},
                        {.X = -0x1.f4d15140445dap+0, .Y =  0x1.0b9c538bf305dp+0},
                        {.X = -0x1.82cdd47c6ea92p+1, .Y =  0x1.4cb2eb8152f47p+1},
                        {.X = -0x1.8caea5652d959p+1, .Y = -0x1.e36c187856567p-1},
                        {.X = 0x1.f1210ca6c4467p-3, .Y = -0x1.4cce5dc4aa902p-4}})}};
    const auto joined = Union(first, second, {.Length = 0x1.1eb851eb851ecp-4, .Radians = 1e-9});
    REQUIRE(joined.has_value());
    CHECK(std::bit_cast<std::uint64_t>(Area(*joined, Accuracy{}))
          == std::bit_cast<std::uint64_t>(0x1.200cd3fdac7fap+5));
  }
} // namespace Geoxy
