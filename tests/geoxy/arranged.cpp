/***************************************************************************************************
 *  @file      tests/geoxy/arranged.cpp
 *  @brief     Tests that an arrangement made for a boolean, its runs of lines that meet nothing
 *             else one edge each and its areas worked out plainly, traces, winds and stitches as
 *             one made for picking does, to the bit, and that lone joints merge as any vertex.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "arranged.hpp"

#include "real.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/line.hpp"
#include "geoxy/loop.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/region.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <tuple>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A tolerance fine enough for every case drawn about unit size.
     */
    constexpr Tolerance ArrangedTight = {.Length = 1e-9, .Radians = 1e-9};

    /**
     *  A closed run of lines through points, the last back to the first.
     *
     *  @param  points The points.
     *  @return        The lines.
     */
    auto Ringed(const std::vector<Point> &points) -> std::vector<Curve>
    {
      std::vector<Curve> lines;
      lines.reserve(points.size());
      for (std::size_t index = 0; index < points.size(); ++index)
      {
        lines.emplace_back(
            Line{.Start = points[index], .End = points[(index + 1) % points.size()]});
      }
      return lines;
    }

    /**
     *  A star of many sides about a center, its radius waving.
     *
     *  @param  center Its center.
     *  @param  count  How many sides.
     *  @param  size   Its radius, about which it waves by a quarter.
     *  @return        Its sides.
     */
    auto Waved(Point center, int count, double size) -> std::vector<Curve>
    {
      std::vector<Point> points;
      for (int index = 0; index < count; ++index)
      {
        const double angle  = 2 * std::numbers::pi * index / count;
        const double radius = size * (1 + (0.25 * std::sin(5 * angle)));
        points.push_back({.X = center.X + (radius * std::cos(angle)),
                          .Y = center.Y + (radius * std::sin(angle))});
      }
      return Ringed(points);
    }

    /**
     *  Two runs of curves arranged, the first weighing in the first operand and the second in the
     *  second, and traced.
     *
     *  @param  runs      The two runs.
     *  @param  chained   Whether runs of lines may be one edge.
     *  @param  tolerance The tolerance.
     *  @return           The arrangement.
     */
    auto Laid(const std::array<std::span<const Curve>, 2> &runs, bool chained,
              const Tolerance &tolerance = ArrangedTight) -> Arranged
    {
      const Inputs inputs = {
          .Runs = {runs[0], runs[1]}, .Tallies = {Tally{1, 0}, Tally{0, 1}}, .Chained = chained};
      Arranged arranged = Arrange(inputs, tolerance);
      Trace(arranged);
      return arranged;
    }

    /**
     *  Whether two arrangements of the same runs, one chained, trace, wind and stitch alike: the
     *  same cycles, each with the same first half's vertex, the same outer cycles and windings, and
     *  the same loops for the union of the two operands.
     *
     *  @param  runs The two runs.
     *  @return      Whether they do.
     */
    auto Alike(const std::array<std::span<const Curve>, 2> &runs) -> bool
    {
      const Arranged           one   = Laid(runs, false);
      const Arranged           other = Laid(runs, true);
      const std::vector<Tally> wound = Windings(one);
      const std::vector<Tally> winds = Windings(other);
      bool                     same  = one.Cycles.size()  == other.Cycles.size()
                                    && one.Outer          == other.Outer
                                    && one.Samples.size() == other.Samples.size()
                                    && wound              == winds;
      std::vector<bool>        kept(wound.size());
      for (std::size_t index = 0; same && index < wound.size(); ++index)
      {
        kept[index] = wound[index][0] != 0 || wound[index][1] != 0;
        same        = same
                   && one.Points[one.Edges[one.Cycles[index].First / 2].Ends.at(one.Cycles[index].First
                                                                                % 2)]
                          == other.Points[other.Edges[other.Cycles[index].First / 2].Ends.at(
                              other.Cycles[index].First % 2)];
      }
      for (std::size_t index = 0; same && index < one.Samples.size(); ++index)
      {
        same = one.Points[one.Samples[index]] == other.Points[other.Samples[index]];
      }
      return same && Stitched(one, kept).Loops == Stitched(other, kept).Loops;
    }

    /**
     *  A region's loops of lines, each as the ends of its sides from its least side on, the loops
     *  in order: what is the same however the loops are ordered and where each starts.
     *
     *  @param  region The region, of lines.
     *  @return        The loops.
     */
    auto Canonical(const Region &region) -> std::vector<std::vector<std::array<double, 4>>>
    {
      std::vector<std::vector<std::array<double, 4>>> loops;
      for (const Loop &loop : region.Loops)
      {
        std::vector<std::array<double, 4>> sides;
        for (const Curve &edge : loop.Edges)
        {
          const Line &line = std::get<Line>(edge);
          sides.push_back({line.Start.X, line.Start.Y, line.End.X, line.End.Y});
        }
        std::ranges::rotate(sides, std::ranges::min_element(sides));
        loops.push_back(sides);
      }
      std::ranges::sort(loops);
      return loops;
    }

    /**
     *  Whether a boolean of two runs of lines read from the windings ferried along them is the one
     *  the traced faces stitch, the arrangement made as a boolean makes it.
     *
     *  @param  runs      The two runs.
     *  @param  keeps     The rule.
     *  @param  tolerance The tolerance.
     *  @return           Whether it is read, and the same.
     */
    auto Ferries(const std::array<std::span<const Curve>, 2> &runs, Keeps keeps,
                 const Tolerance &tolerance = ArrangedTight) -> bool
    {
      const Inputs inputs   = {.Runs    = {runs[0], runs[1]},
                               .Tallies = {Tally{1, 0}, Tally{0, 1}},
                               .Chained = true,
                               .Spare   = true};
      Arranged     arranged = Arrange(inputs, tolerance);
      Region       ferried  = {};
      const bool   read     = Ferried(arranged, keeps, ferried) == Ferry::Done;
      Trace(arranged);
      const std::vector<Tally> windings = Windings(arranged);
      std::vector<bool>        kept(windings.size());
      for (std::size_t index = 0; index < windings.size(); ++index)
      {
        kept[index] = keeps(windings[index][0] != 0, windings[index][1] != 0);
      }
      return read && Planar(arranged) && Canonical(ferried) == Canonical(Stitched(arranged, kept));
    }

    /**
     *  An arrangement drawn by hand: loops of lines through points, each line one edge between the
     *  vertices at its ends, weighing in the first operand, every vertex a lone joint.
     *
     *  @param  loops The loops' points.
     *  @return       The arrangement, not traced.
     */
    auto Drawn(const std::vector<std::vector<Point>> &loops) -> Arranged
    {
      Arranged arranged;
      arranged.Within = ArrangedTight;
      for (const std::vector<Point> &points : loops)
      {
        const std::size_t first = arranged.Points.size();
        arranged.Runs.push_back(arranged.Sources.size());
        for (std::size_t index = 0; index < points.size(); ++index)
        {
          const std::size_t next = (index + 1) % points.size();
          arranged.Points.push_back(points[index]);
          arranged.Sources.emplace_back(Line{.Start = points[index], .End = points[next]});
          arranged.Edges.push_back({.Source  = arranged.Sources.size() - 1,
                                    .Stretch = {.Start = 0.0, .End = 1.0},
                                    .Ends    = {first + index, first + next},
                                    .Weighs  = {1, 0}});
          arranged.Tallies.push_back({1, 0});
          arranged.Lone.push_back({Arranged::Joint, Arranged::Joint});
        }
      }
      arranged.Runs.push_back(arranged.Sources.size());
      return arranged;
    }

    /**
     *  The rule that keeps what the first operand fills.
     *
     *  @param  first Whether the first fills it.
     *  @return       Whether kept.
     */
    auto Firstly(bool first, bool /*second*/) -> bool
    {
      return first;
    }
  } // namespace

  TEST_CASE("arranged: a joint beside a line no longer than the length merges as any vertex")
  {
    /**
     *  A stub and a long line joined at its end, then another, meeting nothing else: the stub's
     *  ends merge where it is no longer than the length, by its length where its coordinates' gaps
     *  are shorter; chained, the two long lines are one edge, their joint no vertex, where the stub
     *  is told at its joint, and so leaves the first of them in no pair.
     */
    for (const auto &[stub, points, strung] :
         {std::tuple{Point{.X = 4e-10, .Y = 3e-10}, 3UZ, 1UZ},
          std::tuple{Point{.X = 9.6e-10, .Y = 7.2e-10}, 4UZ, 0UZ}})
    {
      const std::vector<Curve> walls = {Line{.Start = {}, .End = stub},
                                        Line{.Start = stub, .End = {.X = 1, .Y = -1}},
                                        Line{.Start = {.X = 1, .Y = -1}, .End = {.X = 2, .Y = 3}}};
      for (const bool chained : {false, true})
      {
        const Arranged arranged =
            Arrange({.Runs = {walls}, .Tallies = {Tally{}}, .Chained = chained}, ArrangedTight);
        const std::size_t vertices = points - (chained ? strung : 0);
        std::size_t       lines    = 0;
        for (const Arranged::Edge &edge : arranged.Edges)
        {
          lines += edge.Lines;
        }
        CHECK(arranged.Points.size() == vertices);
        CHECK(arranged.Edges.size()  == vertices - 1);
        CHECK(lines == points - 1);
      }
    }
    /**
     *  Two lines in a row, apart: no joint.
     */
    const std::vector<Curve> apart = {Line{.Start = {}, .End = {.X = 1, .Y = 0}},
                                      Line{.Start = {.X = 0, .Y = 5}, .End = {.X = 1, .Y = 5}}};
    for (const bool chained : {false, true})
    {
      const Arranged loose =
          Arrange({.Runs = {apart}, .Tallies = {Tally{}}, .Chained = chained}, ArrangedTight);
      CHECK(loose.Points.size() == 4);
      CHECK(loose.Edges.size()  == 2);
    }
  }

  TEST_CASE("arranged: lines in a row that meet nothing else are one edge, and trace alike")
  {
    const std::vector<Curve> star  = Waved({.X = 1.5, .Y = 1.5}, 96, 1.0);
    const std::vector<Curve> other = Waved({.X = 1.8, .Y = 1.7}, 64, 1.0);
    const std::vector<Curve> apart = Waved({.X = 9.0, .Y = 1.5}, 40, 0.5);
    const std::vector<Curve> small = Waved({.X = 1.5, .Y = 1.5}, 48, 0.3);
    CHECK(Alike({star, other}));
    CHECK(Alike({other, star}));
    CHECK(Alike({star, apart}));
    CHECK(Alike({star, small}));
    /**
     *  A square and one inside it on part of its side: a stretch shared beside runs of lines.
     */
    const std::vector<Curve> square = Waved({.X = 2, .Y = 2}, 4, 2.0);
    const std::vector<Curve> inside =
        Ringed({{.X = 1, .Y = 0}, {.X = 3, .Y = 0}, {.X = 3, .Y = 2}, {.X = 1, .Y = 2}});
    const std::vector<Curve> outer = Ringed({{.X = 0, .Y = 0},
                                             {.X = 2, .Y = 0},
                                             {.X = 4, .Y = 0},
                                             {.X = 4, .Y = 2},
                                             {.X = 4, .Y = 4},
                                             {.X = 2, .Y = 4},
                                             {.X = 0, .Y = 4},
                                             {.X = 0, .Y = 2}});
    CHECK(Alike({outer, inside}));
    CHECK(Alike({inside, outer}));
    CHECK(Alike({square, inside}));
    const Arranged chained = Laid({star, other}, true);
    std::size_t    lines   = 0;
    for (const Arranged::Edge &edge : chained.Edges)
    {
      lines += edge.Lines;
    }
    CHECK(chained.Edges.size() < Laid({star, other}, false).Edges.size());
    CHECK(lines > chained.Edges.size());
  }

  TEST_CASE("arranged: a sliver's faces told apart exactly where their plain areas lie too near")
  {
    /**
     *  A sliver along a diagonal, its sides 2^-41 off each way: its area is 2^-40, under what the
     *  plain areas' slacks leave, so its outer cycle is taken from the exact areas, as unchained.
     */
    const double       away = 0x1p-41;
    std::vector<Point> points;
    for (int index = 0; index <= 64; ++index)
    {
      points.push_back({.X = index / 64.0, .Y = index / 64.0});
    }
    for (int index = 64; index >= 0; --index)
    {
      points.push_back({.X = (index / 64.0) - away, .Y = (index / 64.0) + away});
    }
    const std::vector<Curve> sliver = Ringed(points);
    const std::vector<Curve> none   = {};
    const Tolerance          fine   = {.Length = 1e-15, .Radians = 1e-9};
    const Arranged           one    = Laid({sliver, none}, false, fine);
    const Arranged           other  = Laid({sliver, none}, true, fine);
    REQUIRE(one.Cycles.size()   == 2);
    REQUIRE(other.Cycles.size() == 2);
    for (std::size_t index = 0; index < 2; ++index)
    {
      CHECK(one.Cycles[index].Area == other.Cycles[index].Area);
    }
    CHECK(one.Outer == other.Outer);
    CHECK(std::abs(one.Cycles[0].Area) == 0x1p-40);
    /**
     *  Run the other way, its outer cycle first, beside a square apart, whose cycles' plain areas
     *  lie far apart: only the sliver's are worked out again.
     */
    std::ranges::reverse(points);
    const std::vector<Curve> turned = Ringed(points);
    const std::vector<Curve> square =
        Ringed({{.X = 3, .Y = 0}, {.X = 4, .Y = 0}, {.X = 4, .Y = 1}, {.X = 3, .Y = 1}});
    const Arranged both = Laid({turned, square}, false, fine);
    const Arranged each = Laid({turned, square}, true, fine);
    REQUIRE(both.Cycles.size() == each.Cycles.size());
    CHECK(both.Outer           == each.Outer);
    CHECK(Windings(both)       == Windings(each));
  }

  TEST_CASE("arranged: lines through a merged vertex that meet elsewhere twist it, and are split")
  {
    /**
     *  The smallest case found, from two random polygons of 4,000 vertices: three lines whose
     *  crossings lie 0.90e-6, 0.89e-6 and 1.65e-6 apart, against a length of 0.93e-6. The first two
     *  merge; the first line and the third pass through that vertex and meet at the third crossing,
     *  so the halves about the two are no plane's.
     */
    const std::vector<Curve> first =
        Ringed({{.X = 0x1.417d8d205e980p+11, .Y = 0x1.eaec0a23c08aap+12},
                {.X = 0x1.29369128024fdp+12, .Y = 0x1.3e32cfc222974p+12},
                {.X = 0x1.8b4b3b57ef5c7p+12, .Y = 0x1.8bb107f2339a4p+12},
                {.X = 0x1.02d370c01c58ap+12, .Y = 0x1.1c352129e8c6dp+13}});
    const std::vector<Curve> second =
        Ringed({{.X = 0x1.449b6a8da6934p+12, .Y = 0x1.70a80f0b025e9p+7},
                {.X = 0x1.441e9cfa9de44p+11, .Y = 0x1.29c795f2b5056p+13},
                {.X = 0x1.376a7e9af9df5p+13, .Y = 0x1.97cf9e535f909p+11},
                {.X = 0x1.2ddfdc03b4f26p+11, .Y = 0x1.ce9574f1eceb6p+12}});
    const Tolerance wide = {.Length = 9.344619977064311e-07, .Radians = 1e-10};
    for (const bool chained : {false, true})
    {
      Inputs inputs = {
          .Runs = {first, second}, .Tallies = {Tally{1, 0}, Tally{0, 1}}, .Chained = chained};
      Arranged twisted = Arrange(inputs, wide);
      Trace(twisted);
      CHECK_FALSE(twisted.Twisted);
      CHECK_FALSE(Planar(twisted));
      inputs.Untwisted = true;
      Arranged split   = Arrange(inputs, wide);
      Trace(split);
      CHECK(split.Twisted);
      CHECK(Planar(split));
      CHECK(split.Points.size() == twisted.Points.size() + 1);
    }
    /**
     *  A square whose two lower corners a triangle's sides pass 1e-10 inside, crossing both sides
     *  at each within the length: each corner's cluster has its lines meet two by two in it, and is
     *  kept, as the lines through one point.
     */
    const std::vector<Curve> square =
        Ringed({{.X = 0, .Y = 0}, {.X = 4, .Y = 0}, {.X = 4, .Y = 4}, {.X = 0, .Y = 4}});
    const std::vector<Curve> triangle =
        Ringed({{.X = -1, .Y = 1 + 1e-10}, {.X = 2, .Y = -2 + 1e-10}, {.X = 5, .Y = 1 + 1e-10}});
    const Inputs inputs = {.Runs      = {square, triangle},
                           .Tallies   = {Tally{1, 0}, Tally{0, 1}},
                           .Chained   = true,
                           .Untwisted = true};
    Arranged     kept   = Arrange(inputs, ArrangedTight);
    Trace(kept);
    CHECK_FALSE(kept.Twisted);
    CHECK(Planar(kept));
    /**
     *  Two squares apart, nothing merged across a gap: nothing to tell.
     */
    const std::vector<Curve> apart =
        Ringed({{.X = 6, .Y = 0}, {.X = 7, .Y = 0}, {.X = 7, .Y = 1}, {.X = 6, .Y = 1}});
    const Inputs lone = {
        .Runs = {square, apart}, .Tallies = {Tally{1, 0}, Tally{0, 1}}, .Untwisted = true};
    Arranged both = Arrange(lone, ArrangedTight);
    Trace(both);
    CHECK_FALSE(both.Twisted);
    CHECK(Planar(both));
  }

  TEST_CASE("arranged: pieces between two merged vertices that cross there are made one")
  {
    const auto ring = [](const std::vector<Vertex> &vertices) -> std::vector<Curve>
    {
      return {Polyline{.Vertices = vertices, .Closed = true}};
    };
    /**
     *  Each from random polygons of up to 25 vertices: the first arrangement and the one split are
     *  no plane's, and the pieces made one leave one that is. Two triangles at 1e-2 of their size:
     *  a side passes 0.0025 from the second's corner, its foot the corner's own place, and the two
     *  cuts one point to Untwisted; it meets the side before the corner 0.0037 on, and the two
     *  pieces between are drawn as one segment. Two triangles of arcs at 1e-2, one passing 40 from
     *  the other's corner. Five arcs and three of their chords at 1e-12, the chords' crossings
     *  2e-12 from their ends: three pairs of pieces between the same two vertices, the two short
     *  ones made one, the long one a face of its own and too far apart.
     */
    const std::vector<std::tuple<std::vector<Curve>, std::vector<Curve>, double, std::size_t>>
        cases = {
            {ring({{.Position = {.X = -0x1.3000617251c3fp+0, .Y = 0x1.a7640e41caee3p-1}},
                   {.Position = {.X = -0x1.3d189a0e23de9p+0, .Y = 0x1.3ea6ff1bba511p-3}},
                   {.Position = {.X = 0x1.dd6a3dbfeebbdp-1, .Y = -0x1.0fb25c4a70e91p+0}}}),
             ring({{.Position = {.X = -0x1.5e8b290a7ae33p+0, .Y = -0x1.2677ce00a54c4p-1}},
                   {.Position = {.X = 0x1.f0f743bf601e2p-2, .Y = 0x1.5f6c4675a0b16p-2}},
                   {.Position = {.X = -0x1.0aaf3acc1f4d9p+0, .Y = 0x1.547ab6e816ceep+0}}}),
             0x1.47ae147ae147bp-7, 1},
            {ring({{.Position = {.X = -0x1.1385a1656fc39p+14, .Y = -0x1.9aecb3a80a638p+11},
                    .Bulge    = 0x1.cea51c5fc53ecp-3},
                   {.Position = {.X = 0x1.46f7d1fbec5e3p+13, .Y = -0x1.7e038f4d9268dp+13},
                    .Bulge    = 0x1.cea51c5fc53ecp-3},
                   {.Position = {.X = 0x1.91311da837cefp+11, .Y = 0x1.1c2f4def97181p+13},
                    .Bulge    = 0x1.cea51c5fc53ecp-3}}),
             ring({{.Position = {.X = 0x1.80332862d750dp+12, .Y = 0x1.4de1d815d09bcp+12},
                    .Bulge    = 0x1.cea51c5fc53ecp-3},
                   {.Position = {.X = -0x1.e45dfddfb65b0p+12, .Y = 0x1.5e7681d5bee7ep+14},
                    .Bulge    = 0x1.cea51c5fc53ecp-3},
                   {.Position = {.X = -0x1.9aa04808ab65ap+13, .Y = -0x1.2a24c39918b47p+12},
                    .Bulge    = 0x1.cea51c5fc53ecp-3}}),
             100, 1},
            {ring({{.Position = {.X = -0x1.7d1b4368dd940p+2, .Y = 0x1.7bf7ae176f726p+6},
                    .Bulge    = -0x1.958bcff769946p-5},
                   {.Position = {.X = -0x1.265eb24fad9afp+6, .Y = 0x1.2e9c0dd59ddddp+4},
                    .Bulge    = 0x1.a6628a22e4640p-6},
                   {.Position = {.X = -0x1.4b78781d36562p+5, .Y = -0x1.0506c070fdea3p+6},
                    .Bulge    = 0x1.360e0b07be5acp-4},
                   {.Position = {.X = 0x1.f011adb102ea9p+5, .Y = -0x1.2bf5e3534abfbp+6},
                    .Bulge    = -0x1.6eeaac2e4c3e8p-4},
                   {.Position = {.X = 0x1.9306862f3614bp+6, .Y = 0x1.3eed538a6a092p+5},
                    .Bulge    = -0x1.874cc9ac1f886p-5}}),
             ring({{.Position = {.X = -0x1.7c6dc8baa720ep+6, .Y = 0x1.84df623cb2588p+3}},
                   {.Position = {.X = -0x1.4b78781d36562p+5, .Y = -0x1.0506c070fdea3p+6}},
                   {.Position = {.X = -0x1.265eb24fad9afp+6, .Y = 0x1.2e9c0dd59ddddp+4}}}),
             1e-12, 2},
    };
    for (const auto &[first, second, length, melded] : cases)
    {
      const Tolerance tolerance = {.Length = length, .Radians = 1e-9};
      Inputs          inputs    = {
          .Runs = {first, second}, .Tallies = {Tally{1, 0}, Tally{0, 1}}, .Chained = true};
      Arranged merged = Arrange(inputs, tolerance);
      Trace(merged);
      CHECK_FALSE(Planar(merged));
      inputs.Untwisted = true;
      Arranged split   = Arrange(inputs, tolerance);
      Trace(split);
      CHECK_FALSE(Planar(split));
      const std::size_t edges = merged.Edges.size();
      CHECK(Melded(merged));
      CHECK(Planar(merged));
      CHECK(merged.Edges.size() == edges - melded);
    }
    /**
     *  A chord and an arc 2e-3 off it, against a length of 1e-3, then one 1e-2 off: each a face of
     *  its own, and none made one.
     */
    for (const double bulge : {2e-3, 1e-2})
    {
      const std::vector<Curve> lens =
          ring({{.Position = {.X = -1, .Y = 0}}, {.Position = {.X = 1, .Y = 0}, .Bulge = bulge}});
      Arranged arranged = Laid({lens, {}}, false, {.Length = 1e-3, .Radians = 1e-9});
      CHECK(Planar(arranged));
      CHECK_FALSE(Melded(arranged));
      CHECK(arranged.Edges.size() == 2);
    }
  }

  TEST_CASE("arranged: windings ferried along the runs read the loops the traced faces stitch")
  {
    const Keeps either = [](bool first, bool second) -> bool
    {
      return first || second;
    };
    const Keeps both = [](bool first, bool second) -> bool
    {
      return first && second;
    };
    const std::vector<Curve> star  = Waved({.X = 1.5, .Y = 1.5}, 96, 1.0);
    const std::vector<Curve> other = Waved({.X = 1.8, .Y = 1.7}, 64, 1.0);
    CHECK(Ferries({star, other}, either));
    CHECK(Ferries({star, other}, both));
    /**
     *  A clockwise hexagon whose farthest point west is a lone joint of a run of lines, its outer
     *  side the run's left, met by a square at its east end.
     */
    const std::vector<Curve> hexagon = Ringed({{.X = 10, .Y = 0},
                                               {.X = 8, .Y = -2},
                                               {.X = 2, .Y = -2},
                                               {.X = 0, .Y =  0},
                                               {.X = 2, .Y =  2},
                                               {.X = 8, .Y = 2}});
    const std::vector<Curve> square =
        Ringed({{.X = 9, .Y = -1}, {.X = 12, .Y = -1}, {.X = 12, .Y = 1}, {.X = 9, .Y = 1}});
    CHECK(Ferries({hexagon, square}, either));
    CHECK(Ferries({hexagon, square}, both));
    /**
     *  A bowtie filled, its two lobes meeting at its crossing; and a square with a triangle inside
     *  it smaller than the length, of no edge.
     */
    const std::vector<Curve> bowtie =
        Ringed({{.X = 0, .Y = 0}, {.X = 2, .Y = 2}, {.X = 2, .Y = 0}, {.X = 0, .Y = 2}});
    CHECK(Ferries({bowtie, {}}, Firstly));
    std::vector<Curve> speck =
        Ringed({{.X = 1, .Y = 1}, {.X = 1 + 1e-10, .Y = 1}, {.X = 1, .Y = 1 + 1e-10}});
    CHECK(Ferries({square, speck}, either));
  }

  TEST_CASE("arranged: windings ferried are unread where no outer side is told, or a run breaks")
  {
    Region region = {};
    /**
     *  A run whose second piece starts where the first does not end.
     */
    Arranged broken         = Drawn({{{.X = 0, .Y = 0}, {.X = 4, .Y = 0}, {.X = 4, .Y = 4}}});
    broken.Edges[1].Ends[0] = 0;
    CHECK(Ferried(broken, Firstly, region) == Ferry::Unfit);
    /**
     *  A loop folded flat: west its two halves leave the same way, and both sectors hold the way.
     */
    const Arranged flat = Drawn({{{.X = 0, .Y = 0}, {.X = 4, .Y = 0}, {.X = 8, .Y = 0}}});
    CHECK(Ferried(flat, Firstly, region) == Ferry::Unfit);
    /**
     *  A square one of whose lines runs against its edge: the outer side told west and east
     *  disagree.
     */
    Arranged flipped =
        Drawn({{{.X = 0, .Y = 0}, {.X = 4, .Y = 0}, {.X = 4, .Y = 1}, {.X = 0, .Y = 1}}});
    flipped.Sources[1] = Line{.Start = {.X = 4, .Y = 1}, .End = {.X = 4, .Y = 0}};
    CHECK(Ferried(flipped, Firstly, region) == Ferry::Unfit);
    /**
     *  Two squares apart: two components.
     */
    const Arranged apart =
        Drawn({{{.X = 0, .Y = 0}, {.X = 1, .Y = 0}, {.X = 1, .Y = 1}, {.X = 0, .Y = 1}},
               {{.X = 5, .Y = 0}, {.X = 6, .Y = 0}, {.X = 6, .Y = 1}, {.X = 5, .Y = 1}}});
    CHECK(Ferried(apart, Firstly, region)      == Ferry::Unfit);
    CHECK(Ferried(Arranged{}, Firstly, region) == Ferry::Unfit);
  }

  TEST_CASE("arranged: windings ferried across crowded vertices, ordered as Trace orders them")
  {
    /**
     *  A diamond whose west vertex is crowded, then every one: told from the ordered halves.
     */
    Region   region = {};
    Arranged diamond =
        Drawn({{{.X = 0, .Y = 1}, {.X = 1, .Y = 0}, {.X = 2, .Y = 1}, {.X = 1, .Y = 2}}});
    diamond.Lone[0] = {Arranged::Crowded, Arranged::Crowded};
    REQUIRE(Ferried(diamond, Firstly, region) == Ferry::Done);
    CHECK(Area(region, Accuracy{0})           == doctest::Approx(2.0));
    std::ranges::fill(diamond.Lone, std::array{Arranged::Crowded, Arranged::Crowded});
    region = {};
    REQUIRE(Ferried(diamond, Firstly, region) == Ferry::Done);
    CHECK(Area(region, Accuracy{0})           == doctest::Approx(2.0));
    /**
     *  Two lines crossing within twice the angle: their halves not told by sign but ordered as
     *  Fanned orders them, by number, which twists the faces here.
     */
    const std::vector<Curve> one =
        Ringed({{.X = 0, .Y = 0}, {.X = 10, .Y = 0.1}, {.X = 0, .Y = 1}});
    const std::vector<Curve> other =
        Ringed({{.X = 10, .Y = 0.05}, {.X = 0, .Y = 0.05}, {.X = 5, .Y = -1}});
    const Inputs   inputs = {.Runs    = {one, other},
                             .Tallies = {Tally{1, 0}, Tally{0, 1}},
                             .Chained = true,
                             .Spare   = true};
    const Arranged near   = Arrange(inputs, {.Length = 1e-9, .Radians = 0.05});
    CHECK(Ferried(near, Firstly, region) == Ferry::Twisted);
  }
} // namespace Geoxy

namespace Geoxy
{
  namespace
  {
    /**
     *  A level line 100 long and parallel slanted lines across it, each crossing it at the place
     *  given, in the order given: each crossing's vertex is made in that order.
     *
     *  @param  places Where along the level line each slanted one crosses it, from 0 to 1.
     *  @return        The level line's pieces' stretches, in the edges' order.
     */
    auto Speared(const std::vector<double> &places) -> std::vector<Range>
    {
      const std::vector<Curve> level = {Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 100, .Y = 0}}};
      std::vector<Curve>       slants;
      for (const double place : places)
      {
        const double middle = 100.0 * place;
        slants.emplace_back(
            Line{.Start = {.X = middle - 0.5, .Y = -1.0}, .End = {.X = middle + 0.5, .Y = 1.0}});
      }
      const Inputs       inputs = {.Runs = {level, slants}, .Tallies = {Tally{1, 0}, Tally{0, 1}}};
      const Arranged     arranged = Arrange(inputs, ArrangedTight);
      std::vector<Range> stretches;
      for (const Arranged::Edge &edge : arranged.Edges)
      {
        if (edge.Source == 0)
        {
          stretches.push_back(edge.Stretch);
        }
      }
      return stretches;
    }

    /**
     *  Whether pieces' stretches run one after another from 0 to 1, each ending where the next
     *  starts and longer than nothing, as many as there are crossings and one more.
     *
     *  @param  stretches The stretches.
     *  @param  crossings How many crossings.
     *  @return           Whether they do.
     */
    auto Successive(const std::vector<Range> &stretches, std::size_t crossings) -> bool
    {
      bool after = stretches.size() == crossings + 1
                && Zero(stretches.front().Start)
                && Exactly(stretches.back().End, 1.0);
      for (std::size_t index = 0; after && index < stretches.size(); ++index)
      {
        after = stretches[index].Start < stretches[index].End
             && (index == 0 || Exactly(stretches[index - 1].End, stretches[index].Start));
      }
      return after;
    }
  } // namespace

  TEST_CASE("arranged: a curve's many cuts are put in order along it, spread or crowded")
  {
    /**
     *  Eighty crossings spread along the line, made in a scrambled order, some two a hair apart and
     *  made the later one first: dealt into bins of a few each.
     */
    std::vector<double> spread;
    spread.reserve(96);
    for (int index = 0; index < 80; ++index)
    {
      spread.push_back(0.01 + (0.98 * static_cast<double>((index * 37) % 80) / 80.0));
    }
    for (int index = 0; index < 8; ++index)
    {
      spread.push_back(0.1 + (0.1 * index) - 1e-4);
      spread.push_back(0.1 + (0.1 * index) - 2e-4);
    }
    CHECK(Successive(Speared(spread), spread.size()));
    /**
     *  Seventy crossings within a thousandth of the line: in one bin, sorted whole.
     */
    std::vector<double> crowded;
    crowded.reserve(70);
    for (int index = 0; index < 70; ++index)
    {
      crowded.push_back(0.5 + (1e-5 * static_cast<double>((index * 29) % 70)));
    }
    CHECK(Successive(Speared(crowded), crowded.size()));
  }
} // namespace Geoxy

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether no two of an arrangement's vertices lie within a length of each other, as Merge
     *  leaves them: swept in order along x, each against those after it within the length.
     *
     *  @param  arranged The arrangement.
     *  @param  length   The length.
     *  @return          Whether none do.
     */
    auto Apart(const Arranged &arranged, double length) -> bool
    {
      std::vector<Point> points(arranged.Points.begin(), arranged.Points.end());
      std::ranges::sort(points, [](Point lhs, Point rhs) -> bool { return lhs.X < rhs.X; });
      for (std::size_t index = 0; index < points.size(); ++index)
      {
        for (std::size_t other = index + 1;
             other < points.size() && points[other].X - points[index].X <= length; ++other)
        {
          if (Distance(points[index], points[other]) <= length)
          {
            return false;
          }
        }
      }
      return true;
    }

    /**
     *  How many of an arrangement's vertices lie within a distance of a point.
     *
     *  @param  arranged The arrangement.
     *  @param  point    The point.
     *  @param  within   The distance.
     *  @return          How many.
     */
    auto Around(const Arranged &arranged, Point point, double within) -> std::size_t
    {
      return static_cast<std::size_t>(
          std::ranges::count_if(arranged.Points, [point, within](Point vertex) -> bool
                                { return Distance(vertex, point) <= within; }));
    }
  } // namespace

  TEST_CASE("arranged: thousands of vertices merged as one by one, the lone ones never looked for")
  {
    /**
     *  A grid of 260 level and 260 upright lines, a level and an upright one a thousandth beside
     *  two of them; apart, three lines nearly through one point, their crossings 1e-12 apart, ten
     *  short lines through another, and a square's loop: tens of thousands of vertices, in buckets,
     *  a few near others, the loop's joints twins.
     */
    std::vector<Curve> lines;
    for (int index = 0; index < 260; ++index)
    {
      const double place = (static_cast<double>(index) + 0.37) / 260.0;
      lines.emplace_back(Line{.Start = {.X = -0.1, .Y = place}, .End = {.X = 1.1, .Y = place}});
      lines.emplace_back(Line{.Start = {.X = place, .Y = -0.1}, .End = {.X = place, .Y = 1.1}});
    }
    const double beside = (10.37 / 260.0) + 1e-3;
    lines.emplace_back(Line{.Start = {.X = -0.1, .Y = beside}, .End = {.X = 1.1, .Y = beside}});
    lines.emplace_back(Line{.Start = {.X = beside, .Y = -0.1}, .End = {.X = beside, .Y = 1.1}});
    const Point knot = {.X = 3.0, .Y = 0.5};
    for (const double turn : {0.1, 1.2, 2.3})
    {
      const Vector way = {.X = std::cos(turn), .Y = std::sin(turn)};
      const Vector off = Vector{.X = -way.Y, .Y = way.X} * (1e-12 * turn);
      lines.emplace_back(
          Line{.Start = knot + off - (way * 0.05), .End = knot + off + (way * 0.05)});
    }
    const Point star = {.X = 5.0, .Y = 5.0};
    for (int index = 0; index < 10; ++index)
    {
      const double angle = std::numbers::pi * static_cast<double>(index) / 10.0;
      const Vector way   = {.X = 0.1 * std::cos(angle), .Y = 0.1 * std::sin(angle)};
      lines.emplace_back(Line{.Start = star - way, .End = star + way});
    }
    const std::vector<Curve> square =
        Ringed({{.X = 3, .Y = 3}, {.X = 3.5, .Y = 3}, {.X = 3.5, .Y = 3.5}, {.X = 3, .Y = 3.5}});
    const Inputs   inputs   = {.Runs = {lines, square}, .Tallies = {Tally{1, 0}, Tally{0, 1}}};
    const Arranged arranged = Arrange(inputs, ArrangedTight);
    CHECK(arranged.Points.size() > 68000);
    CHECK(Apart(arranged, arranged.Within.Length));
    CHECK(Around(arranged, knot, 1e-9) == 1);
    CHECK(Around(arranged, star, 1e-9) == 1);
    CHECK(Around(arranged, {.X = 10.37 / 260.0, .Y = beside}, 1e-6) == 1);
    CHECK(Around(arranged, {.X = 3.25, .Y = 3.25}, 0.36)            == 4);
  }
} // namespace Geoxy
