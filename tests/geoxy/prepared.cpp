/***************************************************************************************************
 *  @file      tests/geoxy/prepared.cpp
 *  @brief     Tests that curves made ready meet as the curves themselves do, to the bit, with each
 *             hit's places and with its point alone, on random pairs and at every turn of the
 *             rules.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/prepared.hpp"

#include "abut.hpp"
#include "draws.hpp"
#include "expansion.hpp"
#include "polylines.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
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
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A tight tolerance.
     */
    constexpr Tolerance Tight = {.Length = 1e-9, .Radians = 1e-9};

    /**
     *  Whether two curves made ready and met one way round answer as Intersect does for the curves,
     *  to the bit: the same meeting, and a marking of each hit's point and kind, in order, with the
     *  same stretches.
     *
     *  @param  curves    The two curves, the first first.
     *  @param  tolerance The tolerance.
     *  @return           Whether they do.
     */
    auto Matched(const Curves &curves, const Tolerance &tolerance) -> bool
    {
      const Curve   &first  = curves[0];
      const Curve   &second = curves[1];
      const Prepared one{first};
      const Prepared other{second};
      Intersections  meeting;
      Intersections  readied;
      Marking        marking;
      Intersect(curves, tolerance, meeting);
      Intersect({one, other}, tolerance, readied);
      Intersect({one, other}, tolerance, marking);
      const auto hit = [](const Hit &lhs, const Hit &rhs) -> bool
      {
        return lhs.At     == rhs.At
            && lhs.First  == rhs.First
            && lhs.Second == rhs.Second
            && lhs.Kind   == rhs.Kind;
      };
      const auto mark = [](const Hit &lhs, const Mark &rhs) -> bool
      {
        return rhs == Mark{.At = lhs.At, .Kind = lhs.Kind};
      };
      const auto overlap = [](const Overlap &lhs, const Overlap &rhs) -> bool
      {
        return lhs.First == rhs.First && lhs.Second == rhs.Second && lhs.Along == rhs.Along;
      };
      return std::ranges::equal(meeting.Hits, readied.Hits, hit)
          && std::ranges::equal(meeting.Overlaps, readied.Overlaps, overlap)
          && meeting.Spent == readied.Spent
          && std::ranges::equal(meeting.Hits, marking.Marks, mark)
          && std::ranges::equal(meeting.Overlaps, marking.Overlaps, overlap)
          && meeting.Spent == marking.Spent;
    }

    /**
     *  Whether two curves made ready meet as the curves themselves do, either first.
     *
     *  @param  first     One curve.
     *  @param  second    The other.
     *  @param  tolerance The tolerance.
     *  @return           Whether they do.
     */
    auto Agreed(const Curve &first, const Curve &second, const Tolerance &tolerance) -> bool
    {
      return Matched({first, second}, tolerance) && Matched({second, first}, tolerance);
    }

    /**
     *  Whether two curves made ready meet as the curves themselves do, either first, bit for bit:
     *  each hit's point and places, and each mark's point, the same doubles to the sign of a zero.
     *
     *  @param  first     One curve.
     *  @param  second    The other.
     *  @param  tolerance The tolerance.
     *  @return           Whether they do.
     */
    auto Bitwise(const Curve &first, const Curve &second, const Tolerance &tolerance) -> bool
    {
      const auto bits = [](double value) -> std::uint64_t
      {
        return std::bit_cast<std::uint64_t>(value);
      };
      const auto same = [&bits](Point lhs, Point rhs) -> bool
      {
        return bits(lhs.X) == bits(rhs.X) && bits(lhs.Y) == bits(rhs.Y);
      };
      const auto once = [&](const Curves &curves) -> bool
      {
        const Prepared one{curves[0]};
        const Prepared other{curves[1]};
        Intersections  meeting;
        Intersections  readied;
        Marking        marking;
        Intersect(curves, tolerance, meeting);
        Intersect({one, other}, tolerance, readied);
        Intersect({one, other}, tolerance, marking);
        const auto hit = [&](const Hit &lhs, const Hit &rhs) -> bool
        {
          return same(lhs.At, rhs.At)
              && same(lhs.First.At, rhs.First.At)
              && same(lhs.Second.At, rhs.Second.At)
              && bits(lhs.First.Along) == bits(rhs.First.Along)
              && bits(lhs.Second.Along) == bits(rhs.Second.Along)
              && lhs.Kind == rhs.Kind;
        };
        const auto mark = [&](const Hit &lhs, const Mark &rhs) -> bool
        {
          return same(lhs.At, rhs.At) && lhs.Kind == rhs.Kind;
        };
        return std::ranges::equal(meeting.Hits, readied.Hits, hit)
            && std::ranges::equal(meeting.Hits, marking.Marks, mark)
            && meeting.Overlaps.size() == readied.Overlaps.size()
            && meeting.Overlaps.size() == marking.Overlaps.size()
            && meeting.Spent == readied.Spent
            && meeting.Spent == marking.Spent;
      };
      return once({first, second}) && once({second, first}) && Agreed(first, second, tolerance);
    }

    /**
     *  The points two curves made ready meet at, in order along the first.
     *
     *  @param  first     The first curve.
     *  @param  second    The second.
     *  @param  tolerance The tolerance.
     *  @return           Each mark's point.
     */
    auto Marked(const Curve &first, const Curve &second, const Tolerance &tolerance)
        -> std::vector<Point>
    {
      const Prepared one{first};
      const Prepared other{second};
      Marking        marking;
      Intersect({one, other}, tolerance, marking);
      std::vector<Point> points;
      points.reserve(marking.Marks.size());
      for (const Mark &mark : marking.Marks)
      {
        points.push_back(mark.At);
      }
      return points;
    }

    /**
     *  The marks two curves made ready meet in, in order along the first.
     *
     *  @param  first     The first curve.
     *  @param  second    The second.
     *  @param  tolerance The tolerance.
     *  @return           The marks.
     */
    auto Marks(const Curve &first, const Curve &second, const Tolerance &tolerance)
        -> std::vector<Mark>
    {
      const Prepared one{first};
      const Prepared other{second};
      Marking        marking;
      Intersect({one, other}, tolerance, marking);
      return marking.Marks;
    }

    /**
     *  Where two lines' carriers cross, worked out in expansions and rounded once: the oracle the
     *  exact pass is held to.
     *
     *  @param  first  One line.
     *  @param  second The other, not parallel to it.
     *  @return        The double nearest the crossing.
     */
    auto Precisely(const Line &first, const Line &second) -> Point
    {
      using namespace Expansion;
      const std::array<Number, 2> way   = {Difference(Of(first.End.X), Of(first.Start.X)),
                                           Difference(Of(first.End.Y), Of(first.Start.Y))};
      const std::array<Number, 2> other = {Difference(Of(second.End.X), Of(second.Start.X)),
                                           Difference(Of(second.End.Y), Of(second.Start.Y))};
      const std::array<Number, 2> gap   = {Difference(Of(second.Start.X), Of(first.Start.X)),
                                           Difference(Of(second.Start.Y), Of(first.Start.Y))};
      const Number turn  = Difference(Product(way[0], other[1]), Product(way[1], other[0]));
      const Number reach = Difference(Product(gap[0], other[1]), Product(gap[1], other[0]));
      const auto   axis  = [&turn, &reach](double start, const Number &step) -> double
      {
        return Rounded({.Top    = Sum(Product(Of(start), turn), Product(step, reach)),
                        .Bottom = {},
                        .Root   = {},
                        .Below  = turn},
                       1);
      };
      return {.X = axis(first.Start.X, way[0]), .Y = axis(first.Start.Y, way[1])};
    }

    /**
     *  Where a line's carrier crosses an arc's circle, worked out in expansions and rounded once: P
     *  + d (-b +- sqrt(b^2 - a c)) / a, the root added first.
     *
     *  @param  line The line, cutting the circle.
     *  @param  arc  The arc whose circle.
     *  @return      The doubles nearest the two crossings.
     */
    auto Precisely(const Line &line, const Arc &arc) -> std::array<Point, 2>
    {
      using namespace Expansion;
      const std::array<double, 2> start  = {line.Start.X, line.Start.Y};
      const std::array<Number, 2> way    = {Difference(Of(line.End.X), Of(line.Start.X)),
                                            Difference(Of(line.End.Y), Of(line.Start.Y))};
      const std::array<Number, 2> gap    = {Difference(Of(line.Start.X), Of(arc.Center.X)),
                                            Difference(Of(line.Start.Y), Of(arc.Center.Y))};
      const Number                square = Sum(Product(way[0], way[0]), Product(way[1], way[1]));
      const Number                along  = Sum(Product(way[0], gap[0]), Product(way[1], gap[1]));
      const Number rest = Difference(Sum(Product(gap[0], gap[0]), Product(gap[1], gap[1])),
                                     Product(Of(arc.Radius), Of(arc.Radius)));
      const Number root = Difference(Product(along, along), Product(square, rest));
      const auto   axis = [&](std::size_t index, int side) -> double
      {
        return Rounded(
            {.Top = Difference(Product(Of(start.at(index)), square), Product(way.at(index), along)),
             .Bottom = way.at(index),
             .Root   = root,
             .Below  = square},
            side);
      };
      return {Point{.X = axis(0, 1), .Y = axis(1, 1)}, Point{.X = axis(0, -1), .Y = axis(1, -1)}};
    }

    /**
     *  Where two arcs' circles cross, worked out in expansions and rounded once: O1 + (k g +-
     *  sqrt(4 |g|^2 r1^2 - k^2) g') / (2 |g|^2), g = O2 - O1, g' it turned a quarter, k = |g|^2 +
     *  r1^2 - r2^2.
     *
     *  @param  first  One arc, its circle crossing the other's.
     *  @param  second The other.
     *  @return        The doubles nearest the two crossings.
     */
    auto Precisely(const Arc &first, const Arc &second) -> std::array<Point, 2>
    {
      using namespace Expansion;
      const std::array<double, 2> centre = {first.Center.X, first.Center.Y};
      const std::array<Number, 2> gap    = {Difference(Of(second.Center.X), Of(first.Center.X)),
                                            Difference(Of(second.Center.Y), Of(first.Center.Y))};
      const std::array<Number, 2> turned = {Scaled(gap[1], -1.0), gap[0]};
      const Number                apart  = Sum(Product(gap[0], gap[0]), Product(gap[1], gap[1]));
      const Number                inner  = Product(Of(first.Radius), Of(first.Radius));
      const Number                chord =
          Difference(Sum(apart, inner), Product(Of(second.Radius), Of(second.Radius)));
      const Number big   = Difference(Scaled(Product(apart, inner), 4.0), Product(chord, chord));
      const Number twice = Scaled(apart, 2.0);
      const auto   axis  = [&](std::size_t index, int side) -> double
      {
        return Rounded(
            {.Top    = Sum(Product(Of(centre.at(index)), twice), Product(chord, gap.at(index))),
             .Bottom = turned.at(index),
             .Root   = big,
             .Below  = twice},
            side);
      };
      return {Point{.X = axis(0, 1), .Y = axis(1, 1)}, Point{.X = axis(0, -1), .Y = axis(1, -1)}};
    }

    /**
     *  An arc of a circle, a whole circle one time in four, else of a sweep from about 0.5 to 6
     *  radians either way.
     *
     *  @param  draws  The draws.
     *  @param  centre Its center.
     *  @param  radius Its radius.
     *  @return        The arc.
     */
    auto Bent(Draws &draws, Point centre, double radius) -> Arc
    {
      if (draws.Next() % 4 == 0)
      {
        return Circle(centre, radius);
      }
      const double sweepSize = draws.Near(3.25, 2.75);
      const double sweep     = sweepSize * (draws.Next() % 2 == 0 ? 1.0 : -1.0);
      return {.Center = centre,
              .Radius = radius,
              .Start  = Radians(draws.Near(Pi, Pi)),
              .Sweep  = Radians(sweep)};
    }

    /**
     *  Where and at what size the exact kernels' lab drew a pair.
     */
    struct Placing
    {
      Point  Offset = {};  /**< Where about.             */
      double Scale  = 1.0; /**< The size of its curves.  */
      double Size   = 1.0; /**< The size of its numbers. */
    };

    /**
     *  Where and at what size the lab drew a pair of a kind: about the origin, far from it at
     *  1.2e8, about the origin again for narrow angles, on a CAD drawing of 1000 units, at a size
     *  from 2^-30 to 2^30, and about the origin for crossings halfway between two doubles.
     *
     *  @param  draws The draws.
     *  @param  kind  The kind, 0 to 5.
     *  @return       Where and at what size.
     */
    auto Pictured(Draws &draws, int kind) -> Placing
    {
      if (kind == 1)
      {
        return {.Offset = {.X = draws.Near(1.2345e8, 1e6), .Y = draws.Near(-6.789e7, 1e6)},
                .Scale  = 1.0,
                .Size   = 1e8};
      }
      if (kind == 3)
      {
        return {.Offset = {.X = 500, .Y = 500}, .Scale = 500.0, .Size = 1e3};
      }
      if (kind == 4)
      {
        const double scale = std::ldexp(1.0, static_cast<int>(draws.Next() % 61) - 30);
        return {.Offset = {}, .Scale = scale, .Size = scale};
      }
      return {};
    }

    /**
     *  A line from one point to another.
     *
     *  @param  start Where it starts.
     *  @param  end   Where it ends.
     *  @return       The line, as a curve.
     */
    auto Stroke(Point start, Point end) -> Curve
    {
      return Line{.Start = start, .End = end};
    }

    /**
     *  An arc about the origin, of radius 5.
     *
     *  @param  start Where it starts, in degrees.
     *  @param  sweep How far it runs, in degrees.
     *  @return       The arc, as a curve.
     */
    auto Bow(double start, double sweep) -> Curve
    {
      return Arc{.Center = {}, .Radius = 5, .Start = Degrees(start), .Sweep = Degrees(sweep)};
    }

    /**
     *  A whole circle.
     *
     *  @param  center Its center.
     *  @param  radius Its radius.
     *  @return        The circle, as a curve.
     */
    auto Ring(Point center, double radius) -> Curve
    {
      return Circle(center, radius);
    }

    /**
     *  What the pass for ends tells of two curves made ready, within a tolerance whose length is
     *  the one they are met within: at or above their grain.
     *
     *  @param  first     The first curve.
     *  @param  second    The second.
     *  @param  tolerance The tolerance.
     *  @return           What it tells; nothing where it leaves them to Intersect.
     */
    auto Abutted(const Curve &first, const Curve &second, const Tolerance &tolerance)
        -> std::optional<Abutment>
    {
      const Prepared one{first};
      const Prepared other{second};
      return Abutting::Met({one, other}, tolerance);
    }

    /**
     *  Whether the pass for ends tells two curves, either first, a single mark at a point, of a
     *  kind.
     *
     *  @param  curves    The two curves.
     *  @param  tolerance The tolerance, its length at or above their grain.
     *  @param  mark      The mark.
     *  @return           Whether it does.
     */
    auto Abuts(const Curves &curves, const Tolerance &tolerance, const Mark &mark) -> bool
    {
      const auto single = [&mark](const std::optional<Abutment> &told) -> bool
      {
        return told && told->Count == 1 && told->Marks[0] == mark;
      };
      return single(Abutted(curves[0], curves[1], tolerance))
          && single(Abutted(curves[1], curves[0], tolerance));
    }

    /**
     *  Whether the pass for ends leaves two curves to Intersect, either first.
     *
     *  @param  one       One curve.
     *  @param  other     The other.
     *  @param  tolerance The tolerance.
     *  @return           Whether it does.
     */
    auto Left(const Curve &one, const Curve &other, const Tolerance &tolerance) -> bool
    {
      return !Abutted(one, other, tolerance) && !Abutted(other, one, tolerance);
    }
  } // namespace

  TEST_CASE("prepared: random lines and arcs made ready meet as the curves themselves, to the bit")
  {
    Draws      draws{0x5EED};
    const auto pick = [&draws](std::uint64_t count) -> std::uint64_t
    {
      return draws.Next() % count;
    };
    const auto either = [&pick]() -> double
    {
      return pick(2) == 0 ? 1.0 : -1.0;
    };
    /**
     *  How many tolerances a crossing is set from an end: none, within the tolerance, about the
     *  four and the sixty-four that decide, and far.
     */
    const std::array<double, 12> nears  = {0, 0.3, 1, 3.9, 4.1, 30, 63.5, 64.5, 130, 1e3, 1e5, 1e7};
    int                          differ = 0;
    const int                    rounds = Rounds(20000);
    for (int round = 0; round < rounds; ++round)
    {
      const double    size      = std::ldexp(1.0, static_cast<int>(pick(10)) - 3);
      const double    length    = std::ldexp(size, -static_cast<int>(pick(26)) - 8);
      const Tolerance tolerance = {.Length  = length,
                                   .Radians = pick(8) == 0 ? draws.Near(0.15, 0.15) : 1e-9};
      const auto      near      = [&]() -> double
      {
        return nears.at(pick(nears.size())) * length * either();
      };
      const auto apart = [&]() -> double
      {
        return pick(3) == 0 ? near() : size * draws.Near(1.6, 1.5);
      };
      /**
       *  An arc, a point of it within its sweep or some tolerances from an end, and through that
       *  point a line at a slight or a plain angle, another line, and another arc, each of their
       *  own ends set as near.
       */
      const double baseX = draws.Near(0, 3);
      const Point  base  = {.X = baseX * size * (pick(3) == 0 ? 1e4 : 1.0),
                            .Y = draws.Near(0, 3) * size};
      const double radius =
          pick(6) == 0 ? length * draws.Near(4500, 4450) : size * draws.Near(1.6, 1.4);
      const double start = draws.Near(0, Pi);
      const double sweep = pick(5) == 0 ? 2 * Pi : draws.Near(3.2, 3.0) * either();
      const Arc    arc   = {.Center = base,
                            .Radius = radius,
                            .Start  = pick(2) == 0 ? Radians(start) : Degrees(start * 180 / Pi),
                            .Sweep  = pick(2) == 0 ? Radians(sweep) : Degrees(sweep * 180 / Pi)};
      const double turn = pick(3) == 0 ? start + (sweep * draws.Near(0.5, 0.45))
                                       : (pick(2) == 0 ? start : start + sweep) + (near() / radius);
      const Point  spot = base + (Vector{.X = std::cos(turn), .Y = std::sin(turn)} * radius);
      const auto   lean = [&]() -> double
      {
        return pick(3) == 0 ? draws.Near(0.16, 0.14) : draws.Near(Pi / 2, (Pi / 2) - 0.3);
      };
      const double way     = turn + (Pi / 2) + lean();
      const double other   = way + lean();
      const auto   through = [&](double heading) -> Curve
      {
        const Vector along = {.X = std::cos(heading), .Y = std::sin(heading)};
        return Line{.Start = spot - (along * apart()), .End = spot + (along * apart())};
      };
      const double reach =
          pick(6) == 0 ? length * draws.Near(4500, 4450) : size * draws.Near(1.6, 1.4);
      const double swing   = pick(3) == 0 ? draws.Near(0, 0.3) : draws.Near(Pi / 2, 1.1);
      const double bearing = turn + swing + (pick(2) == 0 ? 0.0 : Pi);
      const double round2  = pick(5) == 0 ? 2 * Pi : draws.Near(3.2, 3.0) * either();
      const double start2  = pick(3) == 0
                               ? bearing - (round2 * draws.Near(0.5, 0.45))
                               : (pick(2) == 0 ? bearing : bearing - round2) + (near() / reach);
      const Curve  second =
          Arc{.Center = spot - (Vector{.X = std::cos(bearing), .Y = std::sin(bearing)} * reach),
              .Radius = reach,
              .Start  = Radians(start2),
              .Sweep  = Radians(round2)};
      const Curve line  = through(way);
      differ           += Agreed(line, arc, tolerance) ? 0 : 1;
      differ           += Agreed(line, through(other), tolerance) ? 0 : 1;
      differ           += Agreed(arc, second, tolerance) ? 0 : 1;
    }
    CHECK(differ == 0);
  }

  TEST_CASE("prepared: lines told across, apart and far, and left to Intersect near every edge")
  {
    const Curve rising = Stroke({.X = 0, .Y = 0}, {.X = 10, .Y = 10});
    const Curve across = Stroke({.X = 0, .Y = 10}, {.X = 10, .Y = 0});
    REQUIRE(Marked(rising, across, Tight).size() == 1);
    CHECK(Marked(rising, across, Tight).front()  == Point{.X = 5, .Y = 5});
    CHECK(Agreed(rising, across, Tight));

    /**
     *  Their carriers crossing past an end of one, of the other, of both; and boxes far apart, each
     *  way.
     */
    const Curve short1 = Stroke({.X = 0, .Y = 8}, {.X = 3, .Y = 6});
    CHECK(Marked(rising, short1, Tight).empty());
    CHECK(Marked(short1, rising, Tight).empty());
    CHECK(Agreed(rising, short1, Tight));
    for (const Vector move : {Vector{.X = 50, .Y = 0}, Vector{.X = -50, .Y = 0},
                              Vector{.X = 0, .Y = 50}, Vector{.X = 0, .Y = -50}})
    {
      const Curve far = Stroke(Point{.X = 0, .Y = 10} + move, Point{.X = 10, .Y = 0} + move);
      CHECK(Marked(rising, far, Tight).empty());
      CHECK(Agreed(rising, far, Tight));
    }

    /**
     *  Left to Intersect: at a slight angle; a crossing within sixty-four tolerances of an end of
     *  either; a line of no length; lines sharing a stretch; and a short line nearly square to a
     *  long one, its shadow on it about a tolerance long, either first.
     */
    const Curve level   = Stroke({.X = 0, .Y = 0}, {.X = 10, .Y = 0});
    const Curve slight  = Stroke({.X = 0, .Y = -0.5}, {.X = 10, .Y = 0.6});
    const Curve ending  = Stroke({.X = 5, .Y = -4}, {.X = 5, .Y = 3e-8});
    const Curve nothing = Stroke({.X = 5, .Y = 0}, {.X = 5, .Y = 0});
    const Curve along   = Stroke({.X = 4, .Y = 0}, {.X = 14, .Y = 0});
    const Curve leaning = Stroke({.X = 5, .Y = -1e-4}, {.X = 5 + 1.5e-9, .Y = 1e-4});
    for (const Curve *other : {&slight, &ending, &nothing, &along, &leaning})
    {
      CHECK(Agreed(level, *other, Tight));
    }
    CHECK(Marked(level, slight, Tight).size()  == 1);
    CHECK(Marked(ending, level, Tight).size()  == 1);
    CHECK(Marked(level, leaning, Tight).size() == 1);

    /**
     *  A crossing that lies halfway between two doubles is not certain in double-double:
     *  Intersect's exact way rounds it, and the mark is its hit.
     */
    const Curve upright = Stroke({.X = 1, .Y = 0}, {.X = 1 + 0x1p-52, .Y = 2});
    const Curve middle  = Stroke({.X = 0, .Y = 1}, {.X = 3, .Y = 1});
    CHECK(Agreed(upright, middle, Tight));
    CHECK(Marked(upright, middle, Tight).size() == 1);

    /**
     *  An angle below nothing, or past what the fast ways allow for, leaves all to Intersect.
     */
    CHECK(Agreed(rising, across, {.Length = 1e-9, .Radians = -1}));
    CHECK(Agreed(rising, across, {.Length = 1e-9, .Radians = 0.2}));
    CHECK(Marked(rising, across, {.Length = 1e-9, .Radians = 0.2}).size() == 1);
  }

  TEST_CASE("prepared: a line and an arc told across, apart and far, and left to Intersect near "
            "every edge")
  {
    const Curve upper = Bow(0, 180);
    const Curve chord = Stroke({.X = -10, .Y = 3}, {.X = 10, .Y = 3});
    const Curve back  = Stroke({.X = 10, .Y = 3}, {.X = -10, .Y = 3});
    const Point left  = {.X = -4, .Y = 3};
    const Point right = {.X = 4, .Y = 3};

    /**
     *  Two crossings, in order along the first: along the line either way, round the arc the way it
     *  sweeps.
     */
    CHECK(Marked(chord, upper, Tight)          == std::vector{left, right});
    CHECK(Marked(back, upper, Tight)           == std::vector{right, left});
    CHECK(Marked(upper, chord, Tight)          == std::vector{right, left});
    CHECK(Marked(Bow(180, -180), chord, Tight) == std::vector{left, right});
    CHECK(Marked(Bow(200, -220), chord, Tight) == std::vector{left, right});
    CHECK(Marked(Bow(-20, 220), back, Tight)   == std::vector{right, left});
    for (const Curve &bow : {upper, Bow(180, -180), Bow(200, -220), Bow(-20, 220)})
    {
      CHECK(Agreed(chord, bow, Tight));
      CHECK(Agreed(back, bow, Tight));
    }

    /**
     *  One crossing off the arc, on its near side or its far side; a whole circle; a line stopping
     *  between the crossings; a line well outside the circle.
     */
    const Curve upright = Stroke({.X = 3, .Y = -9}, {.X = 3, .Y = 9});
    const Curve stopped = Stroke({.X = -10, .Y = 3}, {.X = 0, .Y = 3});
    const Curve outside = Stroke({.X = 3, .Y = 8}, {.X = 8, .Y = 3});
    CHECK(Marked(upright, upper, Tight)        == std::vector{Point{.X = 3, .Y = 4}});
    CHECK(Marked(upright, Bow(0, 270), Tight)  == std::vector{Point{.X = 3, .Y = 4}});
    CHECK(Marked(upright, Bow(0, -270), Tight) == std::vector{Point{.X = 3, .Y = -4}});
    CHECK(Marked(upright, Ring({}, 5), Tight).size() == 2);
    CHECK(Marked(stopped, upper, Tight) == std::vector{left});
    CHECK(Marked(outside, upper, Tight).empty());
    for (const Curve &bow : {upper, Bow(0, 270), Bow(0, -270), Ring({}, 5)})
    {
      CHECK(Agreed(upright, bow, Tight));
      CHECK(Agreed(stopped, bow, Tight));
      CHECK(Agreed(outside, bow, Tight));
    }

    /**
     *  Left to Intersect: a line a few tolerances outside the circle, or inside it; cutting it at a
     *  slight angle; of no length, or of a length whose square is lost below every double; ending
     *  within sixty-four tolerances of a crossing; a crossing as near an end of the arc, either
     *  end, or a whole circle's seam; an arc too small for the tolerance, or of no radius, or
     *  sweeping more than a turn.
     */
    const double past    = (5 + 6e-9) / std::numbers::sqrt2;
    const Curve  grazing = Stroke({.X = past - 3, .Y = past + 3}, {.X = past + 3, .Y = past - 3});
    const Curve  kissing = Stroke({.X = -10, .Y = 5 - 2e-9}, {.X = 10, .Y = 5 - 2e-9});
    const Curve  shallow = Stroke({.X = -10, .Y = 4.97}, {.X = 10, .Y = 4.97});
    const Curve  nothing = Stroke({.X = 4, .Y = 3}, {.X = 4, .Y = 3});
    const Curve  mote    = Stroke({.X = 0, .Y = 0}, {.X = 1e-200, .Y = 1e-200});
    const Curve  ending  = Stroke({.X = -10, .Y = 3}, {.X = 4 + 3e-8, .Y = 3});
    for (const Curve *line : {&grazing, &kissing, &shallow, &nothing, &mote, &ending})
    {
      CHECK(Agreed(*line, upper, Tight));
    }
    const double there = std::atan2(3.0, 4.0) * 180 / Pi;
    const Curve  seam  = Stroke({.X = 5, .Y = -8}, {.X = 5 - 1e-8, .Y = 8});
    CHECK(Agreed(chord, Bow(there + 1e-8, 120), Tight));
    CHECK(Agreed(chord, Bow(-30, 30 + there - 1e-8), Tight));
    CHECK(Agreed(seam, Ring({}, 5), Tight));
    const Curve speck = Arc{.Center = {.X = 4, .Y = 3}, .Radius = 1e-6, .Start = {}, .Sweep = Turn};
    const Curve point = Arc{.Center = {.X = 4, .Y = 3}, .Radius = 0, .Start = {}, .Sweep = Turn};
    const Curve twice = Arc{.Center = {}, .Radius = 5, .Start = {}, .Sweep = Degrees(540)};
    const Curve still = Arc{.Center = {}, .Radius = 5, .Start = {}, .Sweep = Radians(0)};
    for (const Curve *odd : {&speck, &point, &twice, &still})
    {
      CHECK(Agreed(chord, *odd, Tight));
    }
    CHECK(Marked(chord, twice, Tight).size() == Marked(twice, chord, Tight).size());

    /**
     *  One crossing certain in double-double and the other halfway between two doubles, whichever
     *  comes first: Intersect's exact way rounds it.
     */
    const Curve ring   = Ring({.X = 12 + 0x1p-49, .Y = 0}, 5);
    const Curve level  = Stroke({.X = 0, .Y = 3}, {.X = 30, .Y = 3});
    const Curve turned = Stroke({.X = 30, .Y = 3}, {.X = 0, .Y = 3});
    CHECK(Agreed(level, ring, Tight));
    CHECK(Agreed(turned, ring, Tight));
    CHECK(Marked(level, ring, Tight).size() == 2);

    /**
     *  A line that crosses the circle once has that crossing worked out alone, of the carrier's
     *  two, whichever it is along the line: the other, halfway between two doubles, is not looked
     *  at. Where the one wanted is the one halfway, Intersect's exact way rounds it.
     */
    const Point sure    = {.X = 8 + 0x1p-49, .Y = 3};
    const Curve nearer  = Stroke({.X = 0, .Y = 3}, {.X = 12, .Y = 3});
    const Curve coming  = Stroke({.X = 12, .Y = 3}, {.X = 0, .Y = 3});
    const Curve farther = Stroke({.X = 12, .Y = 3}, {.X = 30, .Y = 3});
    CHECK(Marked(nearer, ring, Tight)         == std::vector{sure});
    CHECK(Marked(coming, ring, Tight)         == std::vector{sure});
    CHECK(Marked(ring, nearer, Tight)         == std::vector{sure});
    CHECK(Marked(farther, ring, Tight).size() == 1);
    for (const Curve *line : {&nearer, &coming, &farther})
    {
      CHECK(Agreed(*line, ring, Tight));
    }

    /**
     *  Where the exact pass is not certain of a crossing the ways before it take over: a level line
     *  meeting first a crossing halfway between two doubles, 2 + 2^-53, then one that is not; and a
     *  line all but upright whose crossings lie 1e-290 or so from the axis, below every size the
     *  pass certifies, met twice or once.
     */
    const Curve small  = Ring({.X = 1.75 + 0x1p-52, .Y = 0}, 0.3125);
    const Curve inward = Stroke({.X = 4, .Y = 0.1875}, {.X = 0, .Y = 0.1875});
    CHECK(Agreed(inward, small, Tight));
    CHECK(Marked(inward, small, Tight).size() == 2);
    const Curve wheel   = Ring({.X = 0.3125, .Y = 0.375}, 0.4);
    const Curve slender = Stroke({.X = 2e-290, .Y = 0}, {.X = 0, .Y = 1});
    const Curve brief   = Stroke({.X = 2e-290, .Y = 0}, {.X = 0, .Y = 0.4});
    CHECK(Agreed(slender, wheel, Tight));
    CHECK(Agreed(brief, wheel, Tight));
    CHECK(Marked(slender, wheel, Tight).size() == 2);
    CHECK(Marked(brief, wheel, Tight).size()   == 1);
  }

  TEST_CASE("prepared: two arcs told across, apart and far, and left to Intersect near every edge")
  {
    const Curve ring  = Ring({}, 5);
    const Curve other = Ring({.X = 6, .Y = 0}, 5);
    const Point above = {.X = 3, .Y = 4};
    const Point below = {.X = 3, .Y = -4};

    /**
     *  Two crossings, in order round the first; arcs holding one of them, or none.
     */
    CHECK(Marked(ring, other, Tight) == std::vector{above, below});
    CHECK(Marked(Arc{.Center = {}, .Radius = 5, .Start = {}, .Sweep = -Turn}, other, Tight)
          == std::vector{below, above});
    CHECK(Marked(Bow(0, 180), other, Tight)   == std::vector{above});
    CHECK(Marked(Bow(180, 150), other, Tight) == std::vector{below});
    CHECK(Marked(Bow(90, 180), other, Tight).empty());
    for (const Curve &bow : {ring, Bow(0, 180), Bow(180, 150), Bow(90, 180), Bow(100, -330)})
    {
      CHECK(Agreed(bow, other, Tight));
      CHECK(Agreed(
          bow,
          Arc{.Center = {.X = 6, .Y = 0}, .Radius = 5, .Start = Degrees(90), .Sweep = Degrees(100)},
          Tight));
    }

    /**
     *  Circles well apart with their boxes overlapping, one well within the other, and a small one
     *  crossing a large one, worked from the small one, either first.
     */
    const Curve aside  = Ring({.X = 8, .Y = 8}, 5);
    const Curve within = Ring({.X = 1, .Y = 0}, 2);
    const Curve small  = Ring({.X = 5, .Y = 0}, 0.25);
    CHECK(Marked(ring, aside, Tight).empty());
    CHECK(Marked(ring, within, Tight).empty());
    CHECK(Marked(ring, small, Tight).size() == 2);
    CHECK(Marked(small, ring, Tight).size() == 2);
    for (const Curve *round : {&aside, &within, &small})
    {
      CHECK(Agreed(ring, *round, Tight));
    }

    /**
     *  Circles crossing where a coordinate is exactly zero, which no bound tells, left to
     *  Intersect's exact way; and circles six tolerances apart, their boxes overlapping, left to
     *  Intersect.
     */
    const Curve left  = Ring({.X = -3, .Y = 0}, 5);
    const Curve right = Ring({.X = 3, .Y = 0}, 5);
    CHECK(Agreed(left, right, Tight));
    CHECK(Marked(left, right, Tight) == std::vector{Point{.X = 0, .Y = 4}, Point{.X = 0, .Y = -4}});
    const double diagonal = (10 + 6e-9) / std::numbers::sqrt2;
    CHECK(Agreed(ring, Ring({.X = diagonal, .Y = diagonal}, 5), Tight));

    /**
     *  Left to Intersect: circles touching, or a few tolerances from it, apart or within; crossing
     *  at a slight angle; one too small for the tolerance; a crossing near an end of either arc.
     */
    const std::array<Curve, 6> edges = {
        Ring({.X = 10, .Y = 0}, 5),       Ring({.X = 10 + 6e-9, .Y = 0}, 5),
        Ring({.X = 3 - 6e-9, .Y = 0}, 2), Ring({.X = 9.99, .Y = 0}, 5),
        Ring({.X = 5, .Y = 0}, 1e-6),     Ring({.X = 5, .Y = 0}, 5e-6)};
    for (const Curve &edge : edges)
    {
      CHECK(Agreed(ring, edge, Tight));
    }
    const double there = std::atan2(4.0, 3.0) * 180 / Pi;
    CHECK(Agreed(Bow(there + 1e-8, 120), other, Tight));
    CHECK(Agreed(ring,
                 Arc{.Center = {.X = 6, .Y = 0},
                     .Radius = 5,
                     .Start  = Degrees(180 - there + 1e-8),
                     .Sweep  = Degrees(120)},
                 Tight));
  }

  /**
   *  Lines as the exact kernels' lab drew them: about the origin, far from it at 1.2e8, at angles
   *  from 1e-3 down to 1e-9 radians, on a CAD grid of 0.001, of sizes from 2^-30 to 2^30, and
   *  crossing exactly halfway between two doubles. Each pair meets ready as Intersect meets it, and
   *  every crossing is the double nearest the true one, worked out in expansions.
   */
  TEST_CASE("prepared: lines of every regime cross at the exact point, as Intersect answers")
  {
    Draws      draws{0x1AB5};
    const auto pick = [&draws](std::uint64_t count) -> std::uint64_t
    {
      return draws.Next() % count;
    };
    int       differ  = 0;
    int       wrong   = 0;
    int       crossed = 0;
    const int rounds  = Rounds(6000);
    for (int round = 0; round < rounds; ++round)
    {
      const auto    kind    = static_cast<int>(pick(6));
      const Placing place   = Pictured(draws, kind);
      const double  scale   = place.Scale;
      const Point   offset  = place.Offset;
      const Vector  centre  = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)};
      const double  heading = draws.Near(Pi, Pi);
      const double  turn    = [&draws, &pick, kind]() -> double
      {
        if (kind != 2)
        {
          return draws.Near(1.55, 1.25);
        }
        const double tiny = std::pow(10.0, draws.Near(-6, 3));
        return tiny * (pick(2) == 0 ? 1 : -1);
      }();
      const auto stroke = [&](double angle) -> Line
      {
        const Vector along = {.X = std::cos(angle), .Y = std::sin(angle)};
        Point        start = offset + ((centre - (along * draws.Near(0.55, 0.45))) * scale);
        Point        end   = offset + ((centre + (along * draws.Near(0.55, 0.45))) * scale);
        if (kind == 3)
        {
          for (double *value : {&start.X, &start.Y, &end.X, &end.Y})
          {
            *value = std::round(*value * 1000) / 1000;
          }
        }
        return {.Start = start, .End = end};
      };
      Line first  = stroke(heading);
      Line second = stroke(heading + turn);
      if (kind == 5)
      {
        const double where = draws.Near(0, 1e3);
        const double level = draws.Near(0, 1e3);
        const double step  = std::nextafter(where, Pi * 1e3) - where;
        first = {.Start = {.X = where, .Y = level - 1}, .End = {.X = where + step, .Y = level + 1}};
        second = {.Start = {.X = where - 3, .Y = level}, .End = {.X = where + 4, .Y = level}};
      }
      const Tolerance tolerance  = {.Length = 1e-12 * place.Size, .Radians = 1e-9};
      differ                    += Agreed(first, second, tolerance) ? 0 : 1;
      for (const Mark &mark : Marks(first, second, tolerance))
      {
        if (mark.Kind == Contact::Cross)
        {
          ++crossed;
          wrong += mark.At == Precisely(first, second) ? 0 : 1;
        }
      }
    }
    CHECK(differ == 0);
    CHECK(wrong  == 0);
    CHECK(crossed > rounds / 2);
  }

  /**
   *  Lines and arcs, and arcs and arcs, as the exact kernels' lab drew them: about the origin, far
   *  from it at 1.2e8, of sizes from 2^-30 to 2^30, a circle of radius 1e6 with a short line across
   *  its rim, lines and circles grazing to 1e-12 of the radius, and level and upright lines. Each
   *  pair meets ready as Intersect meets it, and every crossing is the double nearest the true one,
   *  worked out in expansions.
   */
  TEST_CASE(
      "prepared: lines and arcs of every regime cross at the exact points, as Intersect answers")
  {
    Draws      draws{0xC40D};
    const auto pick = [&draws](std::uint64_t count) -> std::uint64_t
    {
      return draws.Next() % count;
    };
    const auto towards = [](double angle) -> Vector
    {
      return {.X = std::cos(angle), .Y = std::sin(angle)};
    };
    const std::array<int, 6> placings = {0, 1, 4, 0, 0, 0};
    int                      differ   = 0;
    int                      wrong    = 0;
    int                      crossed  = 0;
    const int                rounds   = Rounds(3000);
    for (int round = 0; round < rounds; ++round)
    {
      const auto    kind  = static_cast<std::size_t>(pick(6));
      const Placing place = Pictured(draws, placings.at(kind));
      const bool    huge  = kind == 3;
      const bool    graze = kind == 4;
      const double  scale = huge ? 1e6 : place.Scale;
      const Point   centre =
          place.Offset + (Vector{.X = draws.Near(0, 1), .Y = draws.Near(0, 1)} * place.Scale);
      const double radius = scale * draws.Near(1.1, 0.4);
      const double off    = graze ? radius * (1 - std::pow(10.0, draws.Near(-7.5, 4.5)))
                                  : radius * draws.Near(0, 0.9);
      const Vector normal = towards(draws.Near(Pi, Pi));
      const Vector along  = Perpendicular(normal);
      const double half   = std::sqrt(std::max(0.0, (radius * radius) - (off * off)));
      const Point  middle = centre + (normal * off);
      Line         line   = {.Start = middle - (along * (half * draws.Near(1.25, 0.25))),
                             .End   = middle + (along * (half * draws.Near(1.25, 0.25)))};
      if (huge)
      {
        const Point  rim  = centre + (normal * radius);
        const Vector tilt = towards(draws.Near(Pi, Pi));
        line              = {.Start = rim - tilt, .End = rim + tilt};
      }
      if (kind == 5)
      {
        const double level = centre.Y + off;
        line               = {.Start = {.X = centre.X - (2 * radius), .Y = level},
                              .End   = {.X = centre.X + (2 * radius), .Y = level}};
        if (pick(2) == 0)
        {
          line = {.Start = {.X = centre.X + off, .Y = centre.Y + (2 * radius)},
                  .End   = {.X = centre.X + off, .Y = centre.Y - (2 * radius)}};
        }
      }
      const Arc       arc        = Bent(draws, centre, radius);
      const Tolerance tolerance  = {.Length  = 1e-12 * place.Size * (huge ? 1e6 : 1.0),
                                    .Radians = 1e-9};
      differ                    += Agreed(line, arc, tolerance) ? 0 : 1;
      for (const Mark &mark : Marks(line, arc, tolerance))
      {
        if (mark.Kind == Contact::Cross)
        {
          const std::array<Point, 2> exact = Precisely(line, arc);
          ++crossed;
          wrong += mark.At == exact[0] || mark.At == exact[1] ? 0 : 1;
        }
      }

      const double other   = scale * draws.Near(1.1, 0.4);
      const double spread  = std::abs(radius - other);
      const double apart   = graze ? (radius + other) * (1 - std::pow(10.0, draws.Near(-7.5, 4.5)))
                                   : spread + ((radius + other - spread) * draws.Near(0.5, 0.45));
      const Arc    second  = Bent(draws, centre + (towards(draws.Near(Pi, Pi)) * apart), other);
      differ              += Agreed(arc, second, tolerance) ? 0 : 1;
      for (const Mark &mark : Marks(arc, second, tolerance))
      {
        if (mark.Kind == Contact::Cross)
        {
          const std::array<Point, 2> exact = Precisely(arc, second);
          ++crossed;
          wrong += mark.At == exact[0] || mark.At == exact[1] ? 0 : 1;
        }
      }
    }
    CHECK(differ == 0);
    CHECK(wrong  == 0);
    CHECK(crossed > rounds);
  }

  /**
   *  A level or an upright line holds its coordinate fixed, and a crossing with it keeps it
   *  exactly, either first; a crossing whose coordinates are both zero, or one far below every size
   *  the pass certifies, is left to the ways before it, which tell them exactly; and a line no
   *  longer than two Wide tolerances is told as before.
   */
  TEST_CASE("prepared: level and upright lines, crossings at zero, and short lines")
  {
    const Curve level   = Stroke({.X = -4, .Y = 1}, {.X = 6, .Y = 1});
    const Curve upright = Stroke({.X = 2, .Y = -3}, {.X = 2, .Y = 5});
    const Curve leaning = Stroke({.X = -1, .Y = -2}, {.X = 5, .Y = 4});
    const Curve steep   = Stroke({.X = 0, .Y = 0}, {.X = 3, .Y = 7});
    const Point there   = {.X = 2, .Y = 1};
    CHECK(Marked(level, upright, Tight)   == std::vector{there});
    CHECK(Marked(upright, level, Tight)   == std::vector{there});
    CHECK(Marked(level, leaning, Tight)   == std::vector{there});
    CHECK(Marked(leaning, upright, Tight) == std::vector{there});
    CHECK(Marked(level, steep, Tight)     == std::vector{Point{.X = 3.0 / 7.0, .Y = 1}});
    for (const Curve *other : {&upright, &leaning, &steep})
    {
      CHECK(Agreed(level, *other, Tight));
      CHECK(Agreed(upright, *other, Tight));
    }

    const Curve rising  = Stroke({.X = -1, .Y = -1}, {.X = 1, .Y = 1});
    const Curve falling = Stroke({.X = -1, .Y = 1}, {.X = 1, .Y = -1});
    CHECK(Marked(rising, falling, Tight) == std::vector{Point{}});
    CHECK(Agreed(rising, falling, Tight));
    const Curve slender = Stroke({.X = 2e-300, .Y = -1}, {.X = 0, .Y = 1});
    const Curve axis    = Stroke({.X = -1, .Y = 0.5}, {.X = 1, .Y = 0.5});
    CHECK(Marked(slender, axis, Tight) == std::vector{Point{.X = 2e-300 / 4, .Y = 0.5}});
    CHECK(Agreed(slender, axis, Tight));

    const Tolerance coarse = {.Length = 0.01, .Radians = 1e-9};
    const Curve     brief  = Stroke({.X = 0, .Y = 0}, {.X = 1, .Y = 0});
    const Curve     beyond = Stroke({.X = 1.02, .Y = -1}, {.X = 3, .Y = 1});
    CHECK(Marked(brief, beyond, coarse).empty());
    CHECK(Agreed(brief, beyond, coarse));
  }

  /**
   *  Curves so large or so small that Intersect scales them to a size near 1, and rays of ways so
   *  long or so short, are met ready as Intersect meets them: no mark from what was kept at their
   *  own size.
   */
  TEST_CASE("prepared: curves of any size and rays of any way, as Intersect answers")
  {
    for (const int power : {-300, -150, 150, 300, 1000})
    {
      CAPTURE(power);
      const double    scale = std::ldexp(1.0, power);
      const Tolerance near  = {.Length = 0x1p-30 * scale, .Radians = 1e-9};
      const Curve     line  = Line{.Start = {.X = -3 * scale, .Y = -1 * scale},
                                   .End   = {.X = 4 * scale, .Y = 2 * scale}};
      const Curve     arc   = Arc{.Center = {.X = 0.5 * scale, .Y = 0.25 * scale},
                                  .Radius = 2 * scale,
                                  .Start  = Degrees(10),
                                  .Sweep  = Degrees(250)};
      const Curve     cross = Line{.Start = {.X = -2 * scale, .Y = 2 * scale},
                                   .End   = {.X = 3 * scale, .Y = -2 * scale}};
      CHECK(Agreed(line, arc, near));
      CHECK(Agreed(line, cross, near));
      CHECK(Agreed(arc, arc, near));
      CHECK(Intersect(line, arc, near).Hits.size() == 2);

      const Curve ray = Ray{.Start = {.X = -2, .Y = 3}, .Way = {.X = scale, .Y = -scale}};
      const Curve bar = Line{.Start = {.X = -3, .Y = -1}, .End = {.X = 4, .Y = 2}};
      CHECK(Agreed(ray, bar, Tight));
      CHECK(Intersect(ray, bar, Tight).Hits.size() == 1);

      const Curve endless =
          Xline{.Through = {.X = 0, .Y = 0.5}, .Way = {.X = scale, .Y = 0.125 * scale}};
      CHECK(Agreed(endless, bar, Tight));
      CHECK(Agreed(bar, bar, Tight));
    }

    /**
     *  A spline with one weight far heavier than the rest is met ready as Intersect meets it,
     *  through the spline that stands for it.
     */
    const Curve stiff = Spline{.Degree  = 2,
                               .Knots   = {0, 0, 0, 1, 1, 1},
                               .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = 0}},
                               .Weights = {1, 0x1p60, 1}};
    const Curve level = Line{.Start = {.X = -0.9584593461557238, .Y = 0.4705701413978771},
                             .End   = {.X = 2.97593582170314, .Y = 0.6647300940468062}};
    CHECK(Agreed(stiff, level, Tight));
    CHECK(Intersect(stiff, level, Tight).Hits.size() == 2);
  }

  TEST_CASE("prepared: curves of any other kind, and of numbers that are not finite, as Intersect "
            "answers")
  {
    const Curve  line  = Stroke({.X = 0, .Y = 0}, {.X = 10, .Y = 10});
    const Curve  ray   = Ray{.Start = {.X = 0, .Y = 10}, .Way = {.X = 1, .Y = -1}};
    const Curve  cubic = Bezier::Cubic{.Start = {.X = 0, .Y = 10},
                                       .Out   = {.X = 3, .Y = 9},
                                       .In    = {.X = 7, .Y = 2},
                                       .End   = {.X = 10, .Y = 0}};
    const Curve  far   = Bezier::Cubic{.Start = {.X = 100, .Y = 110},
                                       .Out   = {.X = 103, .Y = 109},
                                       .In    = {.X = 107, .Y = 102},
                                       .End   = {.X = 110, .Y = 100}};
    const double none  = std::numeric_limits<double>::quiet_NaN();
    const Curve  lost  = Stroke({.X = 0, .Y = 10}, {.X = none, .Y = 0});
    const Curve  huge  = Arc{.Center = {.X = 1.7e308, .Y = -1.7e308},
                             .Radius = 1e308,
                             .Start  = {},
                             .Sweep  = Degrees(90)};
    for (const Curve *other : {&ray, &cubic, &far, &lost})
    {
      CHECK(Agreed(line, *other, Tight));
      CHECK(Agreed(Bow(0, 180), *other, Tight));
    }
    CHECK(Marked(line, ray, Tight).size()   == 1);
    CHECK(Marked(line, cubic, Tight).size() == 1);
    CHECK(Marked(line, far, Tight).empty());
    CHECK(Marked(cubic, far, Tight).empty());

    /**
     *  Two pairs the fast ways once told wrongly, found by meeting 240 million pairs both ways. A
     *  short line all but square to a long one, its shadow on it a tolerance long, which
     *  Intersect's thorough way calls a touch: now left to Intersect. And a small circle across a
     *  large arc at its end, whose rough crossings, worked from the large one, were sixty-five
     *  tolerances off: now worked from the small one.
     */
    const Tolerance coarse = {.Length = 2.38419e-07, .Radians = 1e-9};
    const Curve     longer = Stroke({.X = -3.2337450380471235, .Y = -1.7163711986068417},
                                    {.X = -2.0708627097325829, .Y = -0.23947298569964121});
    const Curve     square = Stroke({.X = -2.3275710125169584, .Y = -0.56547386555563661},
                                    {.X = -2.3275333705787689, .Y = -0.565503137064716});
    CHECK(Agreed(longer, square, coarse));
    const Tolerance fine  = {.Length = 4.65661e-10, .Radians = 1e-9};
    const Curve     large = Arc{.Center = {.X = -5.8434853691939237, .Y = -0.87642324881616229},
                                .Radius = 8.1016280913568526,
                                .Start  = Radians(7.6558189267473722),
                                .Sweep  = Radians(-5.3016552256332954)};
    const Curve     small = Arc{.Center = {.X = -11.560557733992361, .Y = 4.8639172783266069},
                                .Radius = 2.9907305153504877e-06,
                                .Start  = Radians(1.0826886953509725),
                                .Sweep  = Radians(-2.0186664875319513)};
    CHECK(Agreed(large, small, fine));
    CHECK(Marked(large, small, fine).size() == 1);

    /**
     *  An arc whose box runs past every double is far from nothing; the same curve met with itself
     *  shares all of it.
     */
    const Prepared one{line};
    const Prepared other{huge};
    Intersections  meeting;
    Intersect({one, other}, Tight, meeting);
    CHECK(meeting.Hits.empty());
    Marking marking;
    Intersect({one, one}, Tight, marking);
    CHECK(marking.Marks.empty());
    CHECK(marking.Overlaps.size() == 1);
    CHECK(&one.Shape()            == &line);

    /**
     *  Marks are one where their points and kinds are.
     */
    const Mark here = {.At = {.X = 1, .Y = 2}, .Kind = Contact::Cross};
    CHECK(here       == Mark{.At = {.X = 1, .Y = 2}, .Kind = Contact::Cross});
    CHECK_FALSE(here == Mark{.At = {.X = 1, .Y = 3}, .Kind = Contact::Cross});
    CHECK_FALSE(here == Mark{.At = {.X = 1, .Y = 2}, .Kind = Contact::End});
  }

  /**
   *  Arcs of one circle, or of two within a hair of each other, meeting at ends they share, along
   *  stretches, at ends a few tolerances apart or nowhere, as Intersect answers: whole circles too,
   *  and circles farther apart, which other passes or Intersect tell.
   */
  TEST_CASE("prepared: arcs of one circle, or of circles a hair apart, as Intersect answers")
  {
    const Point     centre    = {.X = 3.0, .Y = -2.0};
    const Tolerance tolerance = {.Length = 0x1p-30, .Radians = 1e-9};
    const auto      arc       = [&centre](double start, double sweep) -> Curve
    {
      return Arc{.Center = centre, .Radius = 5.0, .Start = Degrees(start), .Sweep = Degrees(sweep)};
    };
    const std::array<std::array<double, 4>, 9> pairs = {{
        {0.0, 90.0, 90.0, 90.0},
        {0.0, 120.0, 60.0, 140.0},
        {0.0, 180.0, 180.0, 180.0},
        {0.0, 300.0, 200.0, 260.0},
        {0.0, 90.0, 200.0, 60.0},
        {0.0, 90.0, 90.0, -45.0},
        {0.0, 360.0, 30.0, 60.0},
        {30.0, 360.0, 30.0, 60.0},
        {0.0, 360.0, 0.0, -360.0},
    }};
    for (const std::array<double, 4> &pair : pairs)
    {
      CHECK(Agreed(arc(pair[0], pair[1]), arc(pair[2], pair[3]), tolerance));
    }
    /**
     *  Ends a quarter, two and twenty tolerances apart, along the circle.
     */
    for (const double shy : {0.25, 2.0, 20.0})
    {
      const double turn = (shy * tolerance.Length / 5.0) * (180.0 / Pi);
      CHECK(Agreed(arc(0.0, 90.0), arc(90.0 + turn, 60.0), tolerance));
    }
    /**
     *  The second circle moved, or its radius changed, by a hundredth, a fifth, a half and two
     *  tolerances.
     */
    for (const double hair : {0.01, 0.2, 0.5, 2.0})
    {
      const double off = hair * tolerance.Length;
      for (const Curve &moved : {Curve{Arc{.Center = {.X = centre.X + off, .Y = centre.Y},
                                           .Radius = 5.0,
                                           .Start  = Degrees(90.0),
                                           .Sweep  = Degrees(90.0)}},
                                 Curve{Arc{.Center = centre,
                                           .Radius = 5.0 + off,
                                           .Start  = Degrees(60.0),
                                           .Sweep  = Degrees(140.0)}}})
      {
        CHECK(Agreed(arc(0.0, 90.0), moved, tolerance));
        CHECK(Agreed(arc(0.0, 120.0), moved, tolerance));
      }
    }
    /**
     *  Arcs of one circle a rounding apart, touching at two ends, far from the origin.
     */
    const Tolerance coarse = {.Length = 1e-8 * 5000001.5, .Radians = 1e-9};
    CHECK(Agreed(Arc{.Center = {.X = 5000000.25, .Y = 3000000.5},
                     .Radius = 1.2499999998,
                     .Start  = Radians(3.141592653589793),
                     .Sweep  = Radians(3.141592653589793)},
                 Arc{.Center = {.X = 5000000.25, .Y = 3000000.4999999995},
                     .Radius = 1.2500000003,
                     .Start  = Radians(3.5e-10),
                     .Sweep  = Radians(3.1415926530)},
                 coarse));
    /**
     *  At random, arcs of one circle, ends at the other's ends or anywhere, in degrees or radians.
     */
    Draws       draws{0x5EED'2048};
    std::size_t differ = 0;
    for (int round = 0; round < Rounds(2000); ++round)
    {
      const double radiusFraction = draws.Near(1.0, 0.5);
      const double radius = std::ldexp(radiusFraction, static_cast<int>(draws.Next() % 41) - 20);
      const Point middle = {.X = draws.Near(0.0, 2.0) * radius, .Y = draws.Near(0.0, 2.0) * radius};
      const double start = draws.Near(0.0, 180.0);
      const double sweep = draws.Near(0.0, 300.0);
      const double joint =
          std::array{start, start + sweep, draws.Near(0.0, 180.0)}.at(draws.Next() % 3);
      const double next = draws.Near(0.0, 300.0);
      const Curve  first =
          Arc{.Center = middle, .Radius = radius, .Start = Degrees(start), .Sweep = Degrees(sweep)};
      const Curve second = draws.Next() % 2 == 0 ? Curve{Arc{.Center = middle,
                                                             .Radius = radius,
                                                             .Start  = Degrees(joint),
                                                             .Sweep  = Degrees(next)}}
                                                 : Curve{Arc{.Center = middle,
                                                             .Radius = radius,
                                                             .Start  = Radians(joint * Pi / 180.0),
                                                             .Sweep  = Radians(next * Pi / 180.0)}};
      differ += Agreed(first, second, {.Length = 1e-9 * radius, .Radians = 1e-9}) ? 0U : 1U;
    }
    CHECK(differ == 0);
  }

  /**
   *  Polylines keep their legs, arcs and grid once made ready, and meet as their curves do, with
   *  every other kind and with each other, of a few legs or many, both ways round; a copy of one
   *  made ready keeps what it keeps.
   */
  TEST_CASE("prepared: polylines made ready, their legs kept, meet as the curves themselves")
  {
    Draws       draws{0x5EED'2046};
    std::size_t differ = 0;
    for (int round = 0; round < Rounds(600); ++round)
    {
      const std::size_t count = round % 4 == 0 ? 25 + (draws.Next() % 40) : 2 + (draws.Next() % 9);
      const Curve       strung =
          round % 2 == 0 ? Curve{Zigzagged(draws, count)} : Curve{Strung(draws, count)};
      const Curve line  = Stroke({.X = draws.Near(0.5, 0.8), .Y = draws.Near(0.5, 0.8)},
                                 {.X = draws.Near(0.5, 0.8), .Y = draws.Near(0.5, 0.8)});
      const Curve arc   = Arc{.Center = {.X = draws.Near(0.5, 0.5), .Y = draws.Near(0.5, 0.5)},
                              .Radius = draws.Near(0.6, 0.4),
                              .Start  = Radians(draws.Near(0, 3)),
                              .Sweep  = Radians(draws.Near(0, 4))};
      const Curve other = Strung(draws, 1 + (draws.Next() % 30));
      for (const Curve *with : {&line, &arc, &other})
      {
        differ += Agreed(strung, *with, Tight) ? 0U : 1U;
      }
    }
    CHECK(differ == 0);
    const Curve    strung = Strung(draws, 30);
    const Curve    line   = Stroke({.X = -1, .Y = 0.5}, {.X = 2, .Y = 0.5});
    const Prepared kept{strung};
    const Prepared across{line};
    Prepared       copy{line};
    copy = kept;
    Marking marking;
    Marking copied;
    Intersect({kept, across}, Tight, marking);
    Intersect({copy, across}, Tight, copied);
    CHECK(marking.Marks == copied.Marks);
    CHECK(Agreed(strung, line, Tight));
  }

  /**
   *  A level line and an upright one, at a corner, at a joint, across each other, or beside each
   *  other, are told by the pass for ends, at every scale and far from the origin, as Intersect
   *  answers; an end near the other's carrier but not on it, a line or a tolerance too short, a
   *  crossing on an axis, or an angle past what the pass allows for, are left to Intersect.
   */
  TEST_CASE(
      "prepared: level and upright lines at a corner, a joint or across, as Intersect answers")
  {
    for (const int power : {-40, 0, 40})
    {
      for (const double offset : {3.0, 1.2345e8})
      {
        CAPTURE(power);
        CAPTURE(offset);
        const double    scale  = std::ldexp(1.0, power);
        const double    length = (offset > 10.0 ? 1e-4 : 1e-9) * scale;
        const Tolerance near   = {.Length = length, .Radians = 1e-9};
        const auto      spot   = [offset, scale](double x, double y) -> Point
        {
          return {.X = (offset + x) * scale, .Y = (offset + y) * scale};
        };
        const auto end = [&near](const Curve &first, const Curve &second, Point point) -> bool
        {
          return Agreed(first, second, near)
              && Abuts({first, second}, near, {.At = point, .Kind = Contact::End});
        };
        CHECK(end(Stroke(spot(0, 0), spot(10, 0)), Stroke(spot(0, 0), spot(0, 10)), spot(0, 0)));
        CHECK(end(Stroke(spot(10, 0), spot(0, 0)), Stroke(spot(0, -10), spot(0, 0)), spot(0, 0)));
        CHECK(end(Stroke(spot(-5, 0), spot(5, 0)), Stroke(spot(2, 0), spot(2, 7)), spot(2, 0)));
        CHECK(end(Stroke(spot(2, 3), spot(9, 3)), Stroke(spot(2, -4), spot(2, 6)), spot(2, 3)));
        const Curve level   = Stroke(spot(-5, 1), spot(5, 1));
        const Curve upright = Stroke(spot(2, -4), spot(2, 6));
        CHECK(Agreed(level, upright, near));
        CHECK(Abuts({level, upright}, near, {.At = spot(2, 1), .Kind = Contact::Cross}));

        /**
         *  An end past four tolerances of the other's carrier: across it where the carriers cross
         *  on both, beside it where they do not; within four tolerances, left to Intersect.
         */
        const Point above    = spot(2, 1) + Vector{.X = 0, .Y = 10 * length};
        const Point beyond   = spot(2, 1) + Vector{.X = 10 * length, .Y = 0};
        const Curve reaching = Stroke(spot(2, -4), above);
        CHECK(Agreed(level, reaching, near));
        CHECK(Abuts({level, reaching}, near, {.At = spot(2, 1), .Kind = Contact::Cross}));
        for (const auto &[one, beside] :
             {std::pair{level, Stroke(above, spot(2, 6))},
              std::pair{upright, Stroke(beyond, spot(9, 1))},
              std::pair{upright, Stroke(spot(-5, 1), spot(2, 1) - Vector{.X = 10 * length})}})
        {
          CHECK(Agreed(one, beside, near));
          const std::optional<Abutment> told = Abutted(one, beside, near);
          CHECK((told && told->Count == 0));
        }
        for (const double gap : {0.25, 0.5, 1.0, 2.0, 3.0, 3.99})
        {
          CAPTURE(gap);
          const Curve short1 = Stroke(spot(2, 1) + Vector{.X = 0, .Y = gap * length}, spot(2, 6));
          const Curve short2 = Stroke(spot(-5, 1), spot(2, 1) - Vector{.X = gap * length, .Y = 0});
          CHECK(Agreed(level, short1, near));
          CHECK(Agreed(short2, upright, near));
          CHECK(Left(level, short1, near));
          CHECK(Left(short2, upright, near));
        }
        const Curve brief = Stroke(spot(2, 1), spot(2, 1) + Vector{.X = 0, .Y = 3 * length});
        CHECK(Agreed(level, brief, near));
        CHECK(Left(level, brief, near));
        CHECK(
            Left(Stroke(spot(2, 1), spot(2, 1) + Vector{.X = 3 * length, .Y = 0}), upright, near));
        CHECK(Agreed(level, upright, {.Length = 0, .Radians = 1e-9}));
        CHECK(Agreed(level, upright, {.Length = length, .Radians = 0.1}));
        CHECK(Left(level, upright, {.Length = length, .Radians = 0.1}));

        /**
         *  What Squarely leaves, a tolerance too fine for it, lines parallel, a slanted line, the
         *  slanted pass tells as Intersect does.
         */
        const Tolerance fine  = {.Length = 0x1p-48 * offset * scale, .Radians = 1e-9};
        const Curve     high  = Stroke(spot(-5, 2), spot(5, 2));
        const Curve     side  = Stroke(spot(3, -4), spot(3, 6));
        const Curve     slant = Stroke(spot(0, -4), spot(4, 6));
        CHECK(Agreed(level, upright, fine));
        CHECK(Agreed(level, high, near));
        CHECK(Agreed(upright, side, near));
        CHECK(Agreed(level, slant, near));
        CHECK(!Left(level, upright, fine));
        CHECK(!Left(level, high, near));
        CHECK(!Left(upright, side, near));
        CHECK(!Left(level, slant, near));
      }
    }

    /**
     *  A corner or a joint on an axis, whose zero Intersect's exact way gives a sign to, of either
     *  sign: left by Squarely, told by the slanted pass as Intersect's quick crossing gives it.
     */
    for (const double zero : {0.0, -0.0})
    {
      const Curve level   = Stroke({.X = zero, .Y = 2}, {.X = 10, .Y = 2});
      const Curve upright = Stroke({.X = zero, .Y = 2}, {.X = zero, .Y = 9});
      const Curve floor   = Stroke({.X = 1, .Y = zero}, {.X = 6, .Y = zero});
      const Curve wall    = Stroke({.X = 3, .Y = zero}, {.X = 3, .Y = 5});
      CHECK(Agreed(level, upright, Tight));
      CHECK(Agreed(floor, wall, Tight));
      CHECK(!Left(level, upright, Tight));
      CHECK(!Left(floor, wall, Tight));
      CHECK(Marks(level, upright, Tight)
            == std::vector{Mark{.At = {.X = 0, .Y = 2}, .Kind = Contact::End}});
    }
  }

  /**
   *  A level or upright line and an arc where an end of either lies at a crossing of the line's
   *  carrier and the arc's circle, a door's swing on its wall and its leaf, a wall ending on an
   *  arc, are told by the pass for ends, either first, at every scale, as Intersect answers; and
   *  left to Intersect where a test is not clear of what Intersect compares it with.
   */
  /**
   *  A line along an axis, at +0 or -0, crossed by a slanted line or cutting a circle: the exact
   *  passes once gave the zero the line holds, where Intersect's quick crossing adds a zero of
   *  either sign to it, or its exact way gives the sign of the ways' cross product.
   */
  TEST_CASE(
      "prepared: a line along an axis keeps Intersect's signed zero where it crosses, to the bit")
  {
    const auto bits = [](const Curve &first, const Curve &second) -> bool
    {
      const Prepared one{first};
      const Prepared other{second};
      Intersections  meeting;
      Marking        marking;
      Intersect({first, second}, Tight, meeting);
      Intersect({one, other}, Tight, marking);
      const auto same = [](const Hit &hit, const Mark &mark) -> bool
      {
        return std::bit_cast<std::uint64_t>(hit.At.X) == std::bit_cast<std::uint64_t>(mark.At.X)
            && std::bit_cast<std::uint64_t>(hit.At.Y) == std::bit_cast<std::uint64_t>(mark.At.Y)
            && hit.Kind                               == mark.Kind;
      };
      return std::ranges::equal(meeting.Hits, marking.Marks, same);
    };
    const Curve circle = Arc{.Center = {.X = 0.3, .Y = 0.2},
                             .Radius = 5.0,
                             .Start  = Radians(0.1),
                             .Sweep  = Radians(6.0)};
    for (const double zero : {0.0, -0.0})
    {
      const Curve floor   = Line{.Start = {.X = 1, .Y = zero}, .End = {.X = 6, .Y = zero}};
      const Curve slant   = Line{.Start = {.X = 2, .Y = -1}, .End = {.X = 4, .Y = 3}};
      const Curve level   = Line{.Start = {.X = -7, .Y = zero}, .End = {.X = 7, .Y = zero}};
      const Curve upright = Line{.Start = {.X = zero, .Y = 7}, .End = {.X = zero, .Y = -7}};
      CHECK(bits(floor, slant));
      CHECK(bits(slant, floor));
      CHECK(bits(level, circle));
      CHECK(bits(circle, level));
      CHECK(bits(upright, circle));
      CHECK(bits(circle, upright));
    }
  }

  TEST_CASE("prepared: a level or upright line and an arc meeting at an end, as Intersect answers")
  {
    for (const int power : {-40, 0, 40})
    {
      for (const double offset : {3.0, 1.2345e8})
      {
        CAPTURE(power);
        CAPTURE(offset);
        const double    scale  = std::ldexp(1.0, power);
        const double    length = (offset > 10.0 ? 1e-3 : 1e-9) * scale;
        const Tolerance near   = {.Length = length, .Radians = 1e-9};
        const auto      spot   = [offset, scale](double x, double y) -> Point
        {
          return {.X = (offset + x) * scale, .Y = (offset + y) * scale};
        };
        const auto bow = [&spot, scale](double start, double sweep) -> Curve
        {
          return Arc{.Center = spot(0, 0),
                     .Radius = 5 * scale,
                     .Start  = Degrees(start),
                     .Sweep  = Degrees(sweep)};
        };
        const auto told = [&near](const Curve &first, const Curve &second,
                                  const std::vector<Mark> &marks) -> bool
        {
          const std::optional<Abutment> abutted = Abutted(first, second, near);
          return Agreed(first, second, near)
              && abutted
              && std::ranges::equal(std::span{abutted->Marks}.first(abutted->Count), marks);
        };
        const Mark start = {.At = spot(5, 0), .Kind = Contact::End};
        for (const Curve &swing : {bow(0, 90), bow(90, -90)})
        {
          CHECK(told(swing, Stroke(spot(5, 0), spot(12, 0)), {start}));
          CHECK(told(Stroke(spot(12, 0), spot(5, 0)), swing, {start}));
          CHECK(told(Stroke(spot(0, 0), spot(0, 5)), swing,
                     {{.At = spot(0, 5), .Kind = Contact::End}}));
          CHECK(told(Stroke(spot(2, 0), spot(9, 0)), swing, {start}));
          CHECK(told(Stroke(spot(-8, 0), spot(5, 0)), swing, {start}));
          CHECK(told(swing, Stroke(spot(4, 3), spot(4, 10)),
                     {{.At = spot(4, 3), .Kind = Contact::End}}));
        }
        const Mark through = {.At = spot(-4, 3), .Kind = Contact::Cross};
        const Mark ending  = {.At = spot(4, 3), .Kind = Contact::End};
        CHECK(told(Stroke(spot(-10, 3), spot(4, 3)), bow(0, 180), {through, ending}));
        CHECK(told(Stroke(spot(4, 3), spot(-10, 3)), bow(0, 180), {ending, through}));
        CHECK(told(bow(0, 180), Stroke(spot(-10, 3), spot(4, 3)), {ending, through}));
        CHECK(told(bow(180, -180), Stroke(spot(-10, 3), spot(4, 3)), {through, ending}));
        CHECK(told(Stroke(spot(-8, 0), spot(5, 0)), bow(0, 270),
                   {{.At = spot(-5, 0), .Kind = Contact::Cross}, start}));
        CHECK(told(Stroke(spot(4, 3), spot(10, 3)), bow(90, 90), {}));

        /**
         *  Left to Intersect: a carrier near tangent; a whole circle; an arc no wider than 1024
         *  tolerances; a slanted line, a short one; an arc's end beyond the line's, or near the
         *  line but not on it; a line's end near the circle but not on it, or just past the
         *  crossing, where Intersect does not place it on the line.
         */
        const Curve wall  = Stroke(spot(-10, 0), spot(5, 0));
        const Curve swing = bow(0, 90);
        for (const Curve &other :
             {Stroke(spot(-3, 4.95), spot(3, 4.95)), Ring(spot(0, 0), 5 * scale),
              Stroke(spot(0, 0), spot(3, 4)),
              Stroke(spot(5, 0), spot(5, 0) + Vector{.X = 3 * length}),
              Stroke(spot(6, 0), spot(9, 0)),
              Stroke(spot(-10, 0) + Vector{.Y = 2 * length}, spot(9, 0) + Vector{.Y = 2 * length}),
              Stroke(spot(-10, 0), spot(5, 0) - Vector{.X = 2 * length}),
              Stroke(spot(-10, 0), spot(5, 0) - Vector{.X = 1e-12 * scale * offset})})
        {
          CHECK(Agreed(swing, other, near));
          CHECK(Left(swing, other, near));
        }
        CHECK(Agreed(Ring(spot(0, 0), 5 * scale), wall, near));
        CHECK(Left(Ring(spot(0, 0), 5 * scale), wall, near));
        CHECK(Agreed(wall, swing, {.Length = scale, .Radians = 1e-9}));
        CHECK(Left(wall, swing, {.Length = scale, .Radians = 1e-9}));
        CHECK(Left(wall, swing, {.Length = 0x1p-43 * offset * scale, .Radians = 1e-9}));
      }
    }

    /**
     *  A swing whose end on its wall is halfway between two doubles, which only Intersect's exact
     *  way rounds; one swept round all but a hair, whose wall crosses its circle where the turns
     *  from both its ends are a hair; a wall on an axis; and curves Intersect does not take as
     *  plain.
     */
    const double odd = 63 + 0x1p-47;
    const Curve  tie =
        Arc{.Center = {.X = odd, .Y = 7}, .Radius = 1, .Start = Degrees(0), .Sweep = Degrees(90)};
    const Curve hair = Arc{.Center = {.X = 3, .Y = 3},
                           .Radius = 5,
                           .Start  = Radians(0),
                           .Sweep  = Radians((2 * Pi) - 0x1p-49)};
    const Curve axis =
        Arc{.Center = {.X = 3, .Y = 0}, .Radius = 5, .Start = Degrees(0), .Sweep = Degrees(90)};
    const double huge = std::numeric_limits<double>::infinity();
    const Curve  lost = Stroke({.X = -huge, .Y = 3}, {.X = 8, .Y = 3});
    const Curve  coiled =
        Arc{.Center = {.X = 3, .Y = 3}, .Radius = 5, .Start = Degrees(0), .Sweep = Degrees(720)};
    for (const auto &[first, second] :
         {std::pair{tie, Stroke({.X = 64, .Y = 7}, {.X = 70, .Y = 7})},
          std::pair{hair, Stroke({.X = -5, .Y = 3}, {.X = 8, .Y = 3})},
          std::pair{coiled, Stroke({.X = 8, .Y = 3}, {.X = 12, .Y = 3})},
          std::pair{Bow(0, 90), lost}, std::pair{Bow(0, 90), Bow(90, 90)}})
    {
      CHECK(Agreed(first, second, Tight));
      CHECK(Left(first, second, Tight));
    }
    /**
     *  A wall on an axis at +0 crosses there at +0, as every way of Intersect's gives it; at -0 it
     *  is left to Intersect, whose quick crossing adds a zero of either sign to it.
     */
    const Curve floor = Stroke({.X = 8, .Y = 0}, {.X = 12, .Y = 0});
    const Curve below = Stroke({.X = 8, .Y = -0.0}, {.X = 12, .Y = -0.0});
    CHECK(Agreed(axis, floor, Tight));
    CHECK(Agreed(axis, below, Tight));
    CHECK(!Left(axis, floor, Tight));
    CHECK(Left(axis, below, Tight));
  }

  /**
   *  Walls of a plan at random, level and upright, their corners, joints and crossings on a grid,
   *  some ends a fraction of the tolerance or a few tolerances off another wall; doors' swings
   *  whose ends lie on walls; and walls ending on arcs, at every scale and far from the origin:
   *  made ready, every pair meets as Intersect answers, to the bit, and the pass for ends tells
   *  most pairs that meet at an end.
   */
  TEST_CASE("prepared: random walls, corners, joints and door swings meet as Intersect answers")
  {
    Draws     draws{0xAB07};
    int       wrong  = 0;
    int       ended  = 0;
    int       told   = 0;
    const int rounds = Rounds(1600);
    for (int round = 0; round < rounds; ++round)
    {
      const double         scale     = std::ldexp(1.0, static_cast<int>(draws.Next() % 81) - 40);
      const bool           far       = draws.Next() % 3 == 0;
      const Point          offset    = {.X = scale * (far ? 1e7 : 50) * draws.Near(0, 1),
                                        .Y = scale * (far ? -3e6 : 50) * draws.Near(0, 1)};
      const double         unit      = scale * (draws.Next() % 2 == 0 ? 1e-6 : 1e-3);
      const double         length    = unit * draws.Near(1, 0.5);
      const Tolerance      tolerance = {.Length = length, .Radians = 1e-9};
      constexpr std::array Fractions = {0.0, 0.0, 0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 5.0, 70.0};
      const auto           spot      = [&offset, scale](double x, double y) -> Point
      {
        return offset + (Vector{.X = x, .Y = y} * scale);
      };
      const auto grid = [&draws](double spread) -> double
      {
        return std::round(draws.Near(0, spread));
      };
      const auto nudge = [&draws, &Fractions, length]() -> double
      {
        const double fraction = Fractions.at(draws.Next() % Fractions.size());
        return fraction * length * (draws.Next() % 2 == 0 ? 1.0 : -1.0);
      };
      const double y                     = grid(5);
      const double x                     = grid(5);
      Point        start                 = spot(grid(8), y);
      Point        bottom                = spot(x, grid(8));
      Point        top                   = spot(x, grid(8));
      start.X                            = draws.Next() % 3 == 0 ? bottom.X + nudge() : start.X;
      top.Y                              = draws.Next() % 3 == 0 ? start.Y + nudge() : top.Y;
      const Point              centre    = spot(grid(5), grid(5));
      const double             turn      = 90.0 * static_cast<double>(draws.Next() % 4);
      const double             sweepSize = draws.Next() % 2 == 0 ? 90.0 : 270.0;
      const double             sweep     = sweepSize * (draws.Next() % 2 == 0 ? 1 : -1);
      const double             radius    = std::round(draws.Near(5, 4)) + 1;
      const Point              cut       = centre + (Direction(Degrees(turn)) * (radius * scale));
      const Point              stop      = cut + (Direction(Degrees(turn)) * (scale * grid(6)));
      const std::vector<Curve> curves    = {
          Stroke(start, spot(grid(8), y)),
          Stroke(bottom, top),
          Arc{.Center = centre,
              .Radius = radius * scale,
              .Start  = Degrees(turn),
              .Sweep  = Degrees(sweep)},
          Stroke(cut, stop + Vector{.X = nudge(), .Y = 0}),
          Stroke(cut + Vector{.X = 0, .Y = nudge()}, centre),
          Arc{.Center = centre,
              .Radius = 5 * scale,
              .Start  = Degrees(360 * draws.Near(0.5, 0.5)),
              .Sweep  = Degrees(360 * draws.Near(0, 0.99))},
          Stroke(centre + (Vector{.X = -10, .Y = 3} * scale),
                 centre + (Vector{.X = 4, .Y = 3} * scale) + Vector{.X = nudge(), .Y = 0}),
          Stroke(centre + (Vector{.X = 4, .Y = -10} * scale),
                 centre + (Vector{.X = 4, .Y = 3} * scale)),
      };
      for (std::size_t first = 0; first < curves.size(); ++first)
      {
        for (std::size_t second = first + 1; second < curves.size(); ++second)
        {
          CAPTURE(round);
          const Curve        &one     = curves.at(first);
          const Curve        &other   = curves.at(second);
          const Intersections meeting = Intersect(one, other, tolerance);
          const bool          end     = std::ranges::any_of(meeting.Hits, [](const Hit &hit) -> bool
                                                            { return hit.Kind == Contact::End; });
          ended += end ? 1 : 0;
          told  += end && Abutted(one, other, tolerance) ? 1 : 0;
          wrong += Agreed(one, other, tolerance) ? 0 : 1;
        }
      }
    }
    CHECK(wrong == 0);
    CHECK(told > ended / 3);
  }

  /**
   *  Prepared curves are met within the tolerance Intersect takes: one whose length is not a
   *  number, which failed every test of a distance and so found nothing, is taken as 0, and one
   *  past every double as the largest double.
   */
  TEST_CASE(
      "prepared: a tolerance not a number, below 0 or past every double, as Intersect takes it")
  {
    constexpr double None    = std::numeric_limits<double>::quiet_NaN();
    constexpr double Endless = std::numeric_limits<double>::infinity();
    const Curve      line    = Stroke({.X = 0, .Y = 0}, {.X = 10, .Y = 10});
    const Curve      across  = Stroke({.X = 0, .Y = 10}, {.X = 10, .Y = 0});
    const Curve      cubic   = Bezier::Cubic{.Start = {.X = 0, .Y = 10},
                                             .Out   = {.X = 3, .Y = 9},
                                             .In    = {.X = 7, .Y = 2},
                                             .End   = {.X = 10, .Y = 0}};
    for (const Tolerance tolerance :
         {Tolerance{.Length = None, .Radians = None}, Tolerance{.Length = -1, .Radians = -1},
          Tolerance{.Length = Endless, .Radians = 1e-9},
          Tolerance{.Length = 1e-9, .Radians = Endless}})
    {
      CAPTURE(tolerance.Length);
      CAPTURE(tolerance.Radians);
      CHECK(Agreed(line, across, tolerance));
      CHECK(Agreed(line, cubic, tolerance));
      CHECK(Agreed(Bow(0, 180), across, tolerance));
      CHECK(Marked(line, across, tolerance).size() == 1);
      CHECK(Marked(line, cubic, tolerance).size()  == 1);
    }
  }

  /**
   *  Rays and endless lines through a point of an arc, of a line or of one another, at a slight or
   *  a plain angle, a ray's start some tolerances before or past that point or far, of ways of any
   *  length, one in eight level or upright, made ready and met as Intersect meets the curves, to
   *  the bit, either first: the exact passes tell each as a line whose end, or both ends, lie at an
   *  infinity of places.
   */
  TEST_CASE("prepared: random rays and endless lines made ready meet as the curves themselves")
  {
    Draws      draws{0x5EED7};
    const auto pick = [&draws](std::uint64_t count) -> std::uint64_t
    {
      return draws.Next() % count;
    };
    const auto either = [&pick]() -> double
    {
      return pick(2) == 0 ? 1.0 : -1.0;
    };
    const std::array<double, 13> nears  = {0,    0.3,  1,   3.9, 4.1, 8,  30,
                                           63.5, 64.5, 130, 1e3, 1e5, 1e7};
    const std::array<Vector, 4>  axes   = {Vector{.X = 1, .Y = 0}, Vector{.X = 0, .Y = 1},
                                           Vector{.X = -1, .Y = 0}, Vector{.X = 0, .Y = -1}};
    int                          differ = 0;
    const int                    rounds = Rounds(20000);
    for (int round = 0; round < rounds; ++round)
    {
      const double    size      = std::ldexp(1.0, static_cast<int>(pick(10)) - 3);
      const double    length    = std::ldexp(size, -static_cast<int>(pick(26)) - 8);
      const Tolerance tolerance = {.Length  = length,
                                   .Radians = pick(8) == 0 ? draws.Near(0.15, 0.15) : 1e-9};
      const auto      near      = [&]() -> double
      {
        return nears.at(pick(nears.size())) * length * either();
      };
      const auto apart = [&]() -> double
      {
        return pick(3) == 0 ? near() : size * draws.Near(1.6, 1.5);
      };
      const auto pace = [&]() -> double
      {
        return pick(4) == 0 ? std::ldexp(either(), static_cast<int>(pick(201)) - 100)
                            : draws.Near(1.0, 0.9);
      };
      const Point  base = {.X = draws.Near(0, 3) * size, .Y = draws.Near(0, 3) * size};
      const double radius =
          pick(6) == 0 ? length * draws.Near(4500, 4450) : size * draws.Near(1.6, 1.4);
      const double start = draws.Near(0, Pi);
      const double sweep = pick(5) == 0 ? 2 * Pi : draws.Near(3.2, 3.0) * either();
      const Arc    arc   = {
          .Center = base, .Radius = radius, .Start = Radians(start), .Sweep = Radians(sweep)};
      const double turn = pick(3) == 0 ? start + (sweep * draws.Near(0.5, 0.45))
                                       : (pick(2) == 0 ? start : start + sweep) + (near() / radius);
      const Point  spot = base + (Vector{.X = std::cos(turn), .Y = std::sin(turn)} * radius);
      const auto   lean = [&]() -> double
      {
        return pick(3) == 0 ? draws.Near(0.16, 0.14) : draws.Near(Pi / 2, (Pi / 2) - 0.3);
      };
      const auto endless = [&](double heading) -> Curve
      {
        const Vector unit = pick(8) == 0 ? axes.at(pick(axes.size()))
                                         : Vector{.X = std::cos(heading), .Y = std::sin(heading)};
        if (pick(2) == 0)
        {
          return Ray{.Start = spot - (unit * apart()), .Way = unit * pace()};
        }
        return Xline{.Through = spot + (unit * apart()), .Way = unit * pace()};
      };
      const double way    = turn + (Pi / 2) + lean();
      const double other  = way + lean();
      const Vector across = {.X = std::cos(other), .Y = std::sin(other)};
      const Curve  first  = endless(way);
      const Curve line = Line{.Start = spot - (across * apart()), .End = spot + (across * apart())};
      differ += Bitwise(first, arc, tolerance) ? 0 : 1;
      differ += Bitwise(first, line, tolerance) ? 0 : 1;
      differ += Bitwise(first, endless(other), tolerance) ? 0 : 1;
    }
    CHECK(differ == 0);
  }

  /**
   *  Rays and endless lines told across, behind their starts, along and beside lines, through arcs
   *  in order along either, outside circles, and left to Intersect at a start on the other curve,
   *  near a circle, at a crossing on an axis or halfway between two doubles; and those Prepared
   *  keeps no way for, or meets at a size or an angle the exact passes leave, as Intersect answers.
   */
  TEST_CASE(
      "prepared: rays and endless lines told across, behind their starts, and left to Intersect")
  {
    const Curve across  = Stroke({.X = 0, .Y = 10}, {.X = 10, .Y = 0});
    const Point middle  = {.X = 5, .Y = 5};
    const Curve toward  = Ray{.Start = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 1}};
    const Curve behind  = Ray{.Start = {.X = 6, .Y = 6}, .Way = {.X = 1, .Y = 1}};
    const Curve back    = Ray{.Start = {.X = 6, .Y = 6}, .Way = {.X = -2, .Y = -2}};
    const Curve through = Xline{.Through = {.X = 2, .Y = 2}, .Way = {.X = 3, .Y = 3}};
    const Curve sloping = Xline{.Through = {.X = 0, .Y = 10}, .Way = {.X = 1, .Y = -1}};
    CHECK(Marked(toward, across, Tight) == std::vector{middle});
    CHECK(Marked(across, back, Tight)   == std::vector{middle});
    CHECK(Marked(behind, across, Tight).empty());
    CHECK(Marked(through, across, Tight) == std::vector{middle});
    CHECK(Marked(through, behind, Tight).empty());
    CHECK(Marked(toward, sloping, Tight)  == std::vector{middle});
    CHECK(Marked(sloping, through, Tight) == std::vector{middle});
    for (const Curve *ray : {&toward, &behind, &back, &through})
    {
      CHECK(Bitwise(*ray, across, Tight));
      CHECK(Bitwise(*ray, sloping, Tight));
      CHECK(Bitwise(*ray, toward, Tight));
      CHECK(Bitwise(*ray, through, Tight));
    }

    /**
     *  Left to Intersect: a start on the other line, a ray along it, one beside it.
     */
    const Curve onto   = Ray{.Start = middle, .Way = {.X = 1, .Y = 0}};
    const Curve along  = Ray{.Start = {.X = 2, .Y = 8}, .Way = {.X = 1, .Y = -1}};
    const Curve beside = Ray{.Start = {.X = 0, .Y = 11}, .Way = {.X = 1, .Y = -1}};
    CHECK(Marks(onto, across, Tight) == std::vector{Mark{.At = middle, .Kind = Contact::End}});
    for (const Curve *ray : {&onto, &along, &beside})
    {
      CHECK(Bitwise(*ray, across, Tight));
      CHECK(Bitwise(*ray, sloping, Tight));
    }

    /**
     *  Through an arc of radius 5: in order along the endless line either way, round the arc the
     *  way it sweeps; a ray from within the circle, either way; past it; level, upright and
     *  slanted; outside a whole circle, or a few tolerances from it.
     */
    const Curve upper = Bow(0, 180);
    const Point left  = {.X = -4, .Y = 3};
    const Point right = {.X = 4, .Y = 3};
    const Curve chord = Xline{.Through = {.X = -10, .Y = 3}, .Way = {.X = 1, .Y = 0}};
    const Curve back2 = Xline{.Through = {.X = 10, .Y = 3}, .Way = {.X = -0.5, .Y = 0}};
    CHECK(Marked(chord, upper, Tight) == std::vector{left, right});
    CHECK(Marked(back2, upper, Tight) == std::vector{right, left});
    CHECK(Marked(upper, chord, Tight) == std::vector{right, left});
    const Curve outward = Ray{.Start = {.X = 0, .Y = 3}, .Way = {.X = 1, .Y = 0}};
    const Curve inward  = Ray{.Start = {.X = 0, .Y = 3}, .Way = {.X = -3, .Y = 0}};
    const Curve past    = Ray{.Start = {.X = 5, .Y = 3}, .Way = {.X = 1, .Y = 0}};
    const Curve upright = Xline{.Through = {.X = 3, .Y = 0}, .Way = {.X = 0, .Y = 2}};
    const Curve slanted = Ray{.Start = {.X = -10, .Y = -9}, .Way = {.X = 1, .Y = 1}};
    const Curve outside = Ray{.Start = {.X = -10, .Y = 6}, .Way = {.X = 1, .Y = 0}};
    const Curve kissing = Ray{.Start = {.X = -10, .Y = 5 - 2e-9}, .Way = {.X = 1, .Y = 0}};
    CHECK(Marked(outward, upper, Tight) == std::vector{right});
    CHECK(Marked(inward, upper, Tight)  == std::vector{left});
    CHECK(Marked(past, upper, Tight).empty());
    CHECK(Marked(upright, upper, Tight)        == std::vector{Point{.X = 3, .Y = 4}});
    CHECK(Marked(upright, Bow(0, -180), Tight) == std::vector{Point{.X = 3, .Y = -4}});
    CHECK(Marked(slanted, Ring({}, 5), Tight).size() == 2);
    CHECK(Marked(outside, Ring({}, 5), Tight).empty());
    for (const Curve *ray :
         {&chord, &back2, &outward, &inward, &past, &upright, &slanted, &outside, &kissing})
    {
      for (const Curve &bow : {upper, Bow(0, -180), Bow(200, -220), Ring({}, 5)})
      {
        CHECK(Bitwise(*ray, bow, Tight));
      }
    }

    /**
     *  A crossing on an axis: a level or upright one at either zero, a way of a negative zero, a
     *  slanted one through a point of the circle at zero, and a line's alike; and a crossing
     *  halfway between two doubles.
     */
    const Curve circle = Arc{.Center = {.X = 0.3, .Y = 0.2},
                             .Radius = 5.0,
                             .Start  = Radians(0.1),
                             .Sweep  = Radians(6.0)};
    const Curve slant  = Stroke({.X = 2, .Y = -1}, {.X = 4, .Y = 3});
    for (const double zero : {0.0, -0.0})
    {
      const Curve level  = Xline{.Through = {.X = -7, .Y = zero}, .Way = {.X = 1, .Y = 0}};
      const Curve rising = Ray{.Start = {.X = zero, .Y = 7}, .Way = {.X = 0, .Y = -1}};
      const Curve flat   = Ray{.Start = {.X = 1, .Y = zero}, .Way = {.X = 1, .Y = zero}};
      const Curve tilted = Ray{.Start = {.X = -7, .Y = 0.5}, .Way = {.X = 1, .Y = -zero}};
      for (const Curve *ray : {&level, &rising, &flat, &tilted})
      {
        CHECK(Bitwise(*ray, circle, Tight));
        CHECK(Bitwise(*ray, slant, Tight));
      }
    }
    const Curve wheel = Arc{.Center = {.X = 3, .Y = 0}, .Radius = 5, .Start = {}, .Sweep = Turn};
    const Curve steep = Xline{.Through = {.X = 0, .Y = 4}, .Way = {.X = 1, .Y = -2}};
    const Curve bar   = Stroke({.X = -1, .Y = 6}, {.X = 5, .Y = -6});
    CHECK(Bitwise(steep, wheel, Tight));
    CHECK(Bitwise(bar, wheel, Tight));
    CHECK(Marked(steep, wheel, Tight).size() == 2);
    const Curve halfway = Ray{.Start = {.X = 1, .Y = 0}, .Way = {.X = 0x1p-52, .Y = 2}};
    const Curve level   = Stroke({.X = 0, .Y = 1}, {.X = 3, .Y = 1});
    CHECK(Bitwise(halfway, level, Tight));
    CHECK(Marked(halfway, level, Tight).size() == 1);

    /**
     *  Rays Prepared keeps no way for: a start not a number, a way of nothing, of ways past 2^96 or
     *  under 2^-96; curves at sizes from 2^-1000 to 2^1000, or below every normal double; at an
     *  angle past what the exact passes allow for; with a polyline, an arc of more than a turn, and
     *  an arc whose size Shrunk scales.
     */
    const double none = std::numeric_limits<double>::quiet_NaN();
    for (const Curve &ray :
         {Curve{Ray{.Start = {.X = none, .Y = 0}, .Way = {.X = 1, .Y = 1}}},
          Curve{Ray{.Start = {.X = 0, .Y = 0}, .Way = {.X = 0, .Y = 0}}},
          Curve{Ray{.Start = {.X = 0, .Y = 0}, .Way = {.X = 0x1p100, .Y = 0x1p100}}},
          Curve{Xline{.Through = {.X = 0, .Y = 0}, .Way = {.X = 0x1p-100, .Y = 0x1p-101}}},
          Curve{Ray{.Start = {.X = 0x1p-1074, .Y = 0}, .Way = {.X = 1, .Y = 1}}}})
    {
      CHECK(Bitwise(ray, across, Tight));
      CHECK(Bitwise(ray, upper, Tight));
    }
    for (const int power : {-1000, -300, 300, 1000})
    {
      CAPTURE(power);
      const double    scale = std::ldexp(1.0, power);
      const Tolerance near  = {.Length = 0x1p-30 * scale, .Radians = 1e-9};
      const Curve ray  = Ray{.Start = {.X = -2 * scale, .Y = 3 * scale}, .Way = {.X = 1, .Y = -1}};
      const Curve line = Line{.Start = {.X = -3 * scale, .Y = -1 * scale},
                              .End   = {.X = 4 * scale, .Y = 2 * scale}};
      const Curve arc  = Arc{.Center = {.X = 0.5 * scale, .Y = 0.25 * scale},
                             .Radius = 2 * scale,
                             .Start  = Degrees(10),
                             .Sweep  = Degrees(250)};
      const Curve endless = Xline{.Through = {.X = 0, .Y = 0.5 * scale}, .Way = {.X = 3, .Y = 1}};
      CHECK(Bitwise(ray, line, near));
      CHECK(Bitwise(ray, arc, near));
      CHECK(Bitwise(endless, arc, near));
      CHECK(Bitwise(endless, ray, near));
    }
    const Curve polyline = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 8}, .Bulge = 0},
                                                 {.Position = {.X = 8, .Y = 0}, .Bulge = 0.5}}};
    const Curve twice    = Arc{.Center = {}, .Radius = 5, .Start = {}, .Sweep = Degrees(540)};
    const Curve vast     = Arc{.Center = {}, .Radius = 0x1p97, .Start = {}, .Sweep = Degrees(90)};
    CHECK(Bitwise(toward, across, {.Length = 1e-9, .Radians = 0.2}));
    CHECK(Bitwise(chord, upper, {.Length = 1e-9, .Radians = 0.2}));
    for (const Curve *other : {&polyline, &twice, &vast})
    {
      CHECK(Bitwise(toward, *other, Tight));
      CHECK(Bitwise(chord, *other, Tight));
    }
  }
} // namespace Geoxy
