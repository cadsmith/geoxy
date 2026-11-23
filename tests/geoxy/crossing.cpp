/***************************************************************************************************
 *  @file      tests/geoxy/crossing.cpp
 *  @brief     Tests for where two lines or conics cross, worked out exactly and rounded once.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "crossing.hpp"

#include "draws.hpp"
#include "exact/whole.hpp"
#include "expansion.hpp"
#include "tangency.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>
#include <variant>

namespace Geoxy::Crossings
{
  namespace
  {
    /**
     *  Where two curves cross near a point, the places taken from the point's nearest on each.
     *
     *  @param  first  One curve.
     *  @param  second The other.
     *  @param  near   Near where.
     *  @return        The crossing, if any.
     */
    auto Near(const Curve &first, const Curve &second, Point near) -> std::optional<Crossed>
    {
      return Crossing({&first, &second}, {Nearest(first, near).Along, Nearest(second, near).Along},
                      near);
    }

    /**
     *  A crossing found, or one at no place: every check on its place then fails.
     *
     *  @param  found The crossing, if any.
     *  @return       It, or one at no place.
     */
    auto Found(const std::optional<Crossed> &found) -> Crossed
    {
      constexpr double Nothing = std::numeric_limits<double>::quiet_NaN();
      return found.value_or(
          Crossed{.At = {.X = Nothing, .Y = Nothing}, .Along = {Nothing, Nothing}});
    }

    /**
     *  Where two lines' carriers cross, worked out in expansions and rounded once, as the exact way
     *  does: the oracle the quick way is held to.
     *
     *  @param  first  One line.
     *  @param  second The other.
     *  @return        The double nearest the crossing.
     */
    auto ExactCrossing(const Line &first, const Line &second) -> Point
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
      const auto   axis  = [&](double start, const Number &step) -> double
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
     *  Splitmix64: numbers that look random, the same on every run.
     */
    struct Mixer
    {
      std::uint64_t State = 0x9E3779B97F4A7C15; /**< Where the sequence is. */

      /**
       *  A double from -1 to 1.
       *
       *  @return The double.
       */
      auto Unit() -> double
      {
        State             += 0x9E3779B97F4A7C15;
        std::uint64_t mix  = State;
        mix                = (mix ^ (mix >> 30U)) * 0xBF58476D1CE4E5B9;
        mix                = (mix ^ (mix >> 27U)) * 0x94D049BB133111EB;
        return (static_cast<double>((mix ^ (mix >> 31U)) >> 11U) * 0x1p-52) - 1.0;
      }
    };

    /**
     *  The crossing of a square root's two signs nearer a point, each worked out in expansions and
     *  rounded once: the oracle the quick conic ways are held to.
     *
     *  @param  axes Each coordinate's quotient, its root's sign to choose.
     *  @param  near The point the nearer is chosen by.
     *  @return      The double nearest the nearer crossing.
     */
    auto Nearer(const std::array<Expansion::Quotient, 2> &axes, Point near) -> Point
    {
      const Point plus = {.X = Expansion::Rounded(axes[0], 1), .Y = Expansion::Rounded(axes[1], 1)};
      const Point minus = {.X = Expansion::Rounded(axes[0], -1),
                           .Y = Expansion::Rounded(axes[1], -1)};
      return Distance(plus, near) <= Distance(minus, near) ? plus : minus;
    }

    /**
     *  Where a line's carrier crosses a circle nearer a point, exactly.
     *
     *  @param  line The line.
     *  @param  arc  The arc whose circle.
     *  @param  near The point.
     *  @return      The double nearest the crossing.
     */
    auto ExactCrossing(const Line &line, const Arc &arc, Point near) -> Point
    {
      using namespace Expansion;
      const std::array<Number, 2> way = {Difference(Of(line.End.X), Of(line.Start.X)),
                                         Difference(Of(line.End.Y), Of(line.Start.Y))};
      const std::array<Number, 2> gap = {Difference(Of(line.Start.X), Of(arc.Center.X)),
                                         Difference(Of(line.Start.Y), Of(arc.Center.Y))};
      const auto                  dot = [](const std::array<Number, 2> &one,
                                           const std::array<Number, 2> &other) -> Number
      {
        return Sum(Product(one[0], other[0]), Product(one[1], other[1]));
      };
      const Number a    = dot(way, way);
      const Number b    = dot(way, gap);
      const Number c    = Difference(dot(gap, gap), Product(Of(arc.Radius), Of(arc.Radius)));
      const Number d    = Difference(Product(b, b), Product(a, c));
      const auto   axis = [&](double start, const Number &step) -> Quotient
      {
        return {.Top    = Difference(Product(Of(start), a), Product(step, b)),
                .Bottom = step,
                .Root   = d,
                .Below  = a};
      };
      return Nearer({axis(line.Start.X, way[0]), axis(line.Start.Y, way[1])}, near);
    }

    /**
     *  Where two circles cross nearer a point, exactly.
     *
     *  @param  first  One arc's circle.
     *  @param  second The other's.
     *  @param  near   The point.
     *  @return        The double nearest the crossing.
     */
    auto ExactCrossing(const Arc &first, const Arc &second, Point near) -> Point
    {
      using namespace Expansion;
      const std::array<Number, 2> gap   = {Difference(Of(second.Center.X), Of(first.Center.X)),
                                           Difference(Of(second.Center.Y), Of(first.Center.Y))};
      const Number                apart = Sum(Product(gap[0], gap[0]), Product(gap[1], gap[1]));
      const Number                inner = Product(Of(first.Radius), Of(first.Radius));
      const Number                chord =
          Difference(Sum(apart, inner), Product(Of(second.Radius), Of(second.Radius)));
      const Number big = Difference(Product(Product(Of(4.0), apart), inner), Product(chord, chord));
      const Number twice = Sum(apart, apart);
      return Nearer(
          {Quotient{.Top    = Sum(Product(Of(first.Center.X), twice), Product(chord, gap[0])),
                    .Bottom = Difference(Number{}, gap[1]),
                    .Root   = big,
                    .Below  = twice},
           Quotient{.Top    = Sum(Product(Of(first.Center.Y), twice), Product(chord, gap[1])),
                    .Bottom = gap[0],
                    .Root   = big,
                    .Below  = twice}},
          near);
    }
  } // namespace

  TEST_CASE("crossing: lines, rays and endless lines, exactly")
  {
    /**
     *  0.1 is not 1/10: the crossing of these is worked from the doubles as they are.
     */
    const Line    one     = {.Start = {.X = 0, .Y = 0}, .End = {.X = 0.3, .Y = 0.1}};
    const Line    other   = {.Start = {.X = 0, .Y = 0.1}, .End = {.X = 0.3, .Y = 0}};
    const Crossed crossed = Found(Near(one, other, {.X = 0.15, .Y = 0.05}));
    CHECK(crossed.At.X == 0.15);
    CHECK(crossed.At.Y == 0.05);
    CHECK(!Near(one, Line{.Start = {.X = 0, .Y = 1}, .End = {.X = 0.3, .Y = 1.1}}, {}));

    /**
     *  Past a line's end either way, behind a ray's start: not on them.
     */
    const Line stub = {.Start = {.X = 0, .Y = 0}, .End = {.X = 1, .Y = 0}};
    CHECK(!Near(stub, Line{.Start = {.X = 2, .Y = -1}, .End = {.X = 2, .Y = 1}}, {.X = 1, .Y = 0}));
    CHECK(!Near(stub, Line{.Start = {.X = -1, .Y = -1}, .End = {.X = -1, .Y = 1}}, {}));
    CHECK(!Near(Ray{.Start = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 0}},
                Line{.Start = {.X = -1, .Y = -1}, .End = {.X = -1, .Y = 1}}, {}));
    const Crossed endless =
        Found(Near(Xline{.Through = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 1}},
                   Ray{.Start = {.X = 2, .Y = 0}, .Way = {.X = 0, .Y = 1}}, {.X = 2, .Y = 2}));
    CHECK(endless.At == Point{.X = 2, .Y = 2});
  }

  TEST_CASE("crossing: polylines' straight segments, and bulges by their exact circles")
  {
    const Curve   path     = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                                   {.Position = {.X = 2, .Y = 0}, .Bulge = 0},
                                                   {.Position = {.X = 2, .Y = 2}, .Bulge = 1},
                                                   {.Position = {.X = 0, .Y = 2}, .Bulge = 0}},
                                      .Closed   = false};
    const Curve   cut      = Line{.Start = {.X = 1, .Y = -1}, .End = {.X = 1, .Y = 3}};
    const Crossed straight = Found(Crossing({&path, &cut}, {0.5, 0.25}, {.X = 1, .Y = 0}));
    CHECK(straight.At       == Point{.X = 1, .Y = 0});
    CHECK(straight.Along[0] == 0.5);
    const Crossed bulged = Found(Crossing({&path, &cut}, {2.5, 0.75}, {.X = 1, .Y = 3}));
    CHECK(bulged.At       == Point{.X = 1, .Y = 3});
    CHECK(bulged.Along[0] == 2.5);
    CHECK(bulged.Along[1] == 1.0);
    /**
     *  The bulge's circle crossed off its arc, below its chord: no crossing of the polyline's.
     */
    const Curve across = Line{.Start = {.X = 0, .Y = 3.5}, .End = {.X = 2.5, .Y = 1}};
    CHECK(Crossing({&path, &across}, {2.5, 0.25}, {.X = 0.6, .Y = 2.9}));
    CHECK(!Crossing({&path, &across}, {2.5, 0.75}, {.X = 1.9, .Y = 1.6}));
    /**
     *  A line of no length, or one missing the bulge's circle: nothing.
     */
    const Curve dot  = Line{.Start = {.X = 1, .Y = 3}, .End = {.X = 1, .Y = 3}};
    const Curve miss = Line{.Start = {.X = 0, .Y = 4}, .End = {.X = 2, .Y = 4}};
    CHECK(!Crossing({&path, &dot}, {2.5, 0.0}, {.X = 1, .Y = 3}));
    CHECK(!Crossing({&path, &miss}, {2.5, 0.5}, {.X = 1, .Y = 4}));
    const Curve beyond = Line{.Start = {.X = 3, .Y = -1}, .End = {.X = 3, .Y = 1}};
    CHECK(!Crossing({&path, &beyond}, {0.9, 0.5}, {.X = 3, .Y = 0}));
    CHECK(Crossing({&cut, &path}, {0.25, 0.5}, {.X = 1, .Y = 0}));
    const Curve none = Polyline{};
    CHECK(!Crossing({&none, &cut}, {0.0, 0.5}, Point{}));
    CHECK(!Crossing({&path, &cut}, {std::numeric_limits<double>::quiet_NaN(), 0.25}, Point{}));
    const Curve before = Line{.Start = {.X = -1, .Y = -1}, .End = {.X = -1, .Y = 1}};
    CHECK(!Crossing({&path, &before}, {0.1, 0.5}, {.X = -1, .Y = 0}));
  }

  TEST_CASE("crossing: a point placed on an ellipse by its angle, on a bulge by its chord's")
  {
    const Curve oval =
        Ellipse{.Map   = {.X = {.X = 2, .Y = 0}, .Y = {.X = 0, .Y = 1}, .Origin = {.X = 1, .Y = 1}},
                .Start = {},
                .Sweep = Degrees(90)};
    const Curve   cut = Line{.Start = {.X = 1, .Y = -2}, .End = {.X = 1, .Y = 4}};
    const Crossed top = Found(Placed({&oval, &cut}, {0.5, 0.5}, {.X = 1, .Y = 2}));
    CHECK(top.Along[0] == 1.0);
    CHECK(top.On[0]    == Point{.X = 1, .Y = 2});

    /**
     *  At its center, no angle, and past its sweep: Nearest places it, not on it.
     */
    CHECK(!Placed({&oval, &cut}, {0.5, 0.5}, {.X = 1, .Y = 1}));
    CHECK(!Placed({&oval, &cut}, {0.5, 0.5}, {.X = 1, .Y = 0}));

    /**
     *  At a bulge's first vertex: the bulge's start.
     */
    const Curve   bulge = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                                                {.Position = {.X = 2, .Y = 0}, .Bulge = 0}},
                                   .Closed   = false};
    const Curve   down  = Line{.Start = {.X = 0, .Y = 1}, .End = {.X = 0, .Y = -1}};
    const Crossed start = Found(Placed({&bulge, &down}, {0.25, 0.5}, Point{}));
    CHECK(start.Along[0] == 0.0);
  }

  TEST_CASE("crossing: a point on a bulge placed where At puts it, at its ends and past them")
  {
    /**
     *  Points of a polyline's bulges, inside them, at their ends and a hair either side: each
     *  placed with the point At gives at its place, to the bit, whether worked out from the bow
     *  Located made or by At where the place falls at an end or rounds onto the next segment.
     */
    Draws       draws{0x5EED'4701};
    std::size_t differ = 0;
    std::size_t placed = 0;
    for (int round = 0; round < Rounds(2000); ++round)
    {
      Polyline polyline = {.Vertices = {}, .Closed = round % 3 == 0};
      for (int vertex = 0; vertex < 5; ++vertex)
      {
        polyline.Vertices.push_back({.Position = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)},
                                     .Bulge    = vertex % 2 == 0 ? draws.Near(0, 1.5) : 0.0});
      }
      const Curve       curve  = polyline;
      const std::size_t index  = 2 * (draws.Next() % 2);
      const Curve       single = Polyline{
          .Vertices = {polyline.Vertices[index], polyline.Vertices[index + 1]}, .Closed = false};
      const double part =
          std::array{0.0, 1.0, draws.Near(0.5, 0.5), 1.0 - 0x1p-45, 1.0 + 0x1p-50, -0x1p-50}.at(
              draws.Next() % 6);
      const std::optional<Place> place =
          Located(curve, static_cast<double>(index) + 0.5, At(single, part));
      placed += place.has_value() ? 1U : 0U;
      differ += !place
                     || std::bit_cast<std::array<std::uint64_t, 2>>(place->At)
                            == std::bit_cast<std::array<std::uint64_t, 2>>(At(curve, place->Along))
                  ? 0U
                  : 1U;
    }
    CHECK(differ == 0);
    CHECK(placed > static_cast<std::size_t>(Rounds(2000)) / 2);

    /**
     *  A point a hair short of a bulge's end, far along a polyline of many: its place rounds onto
     *  the next segment's start, whose point is that vertex, as At gives it.
     */
    Polyline chain = {};
    for (int vertex = 0; vertex <= 9000; ++vertex)
    {
      chain.Vertices.push_back(
          {.Position = {.X = static_cast<double>(vertex), .Y = 0}, .Bulge = 0.5});
    }
    const Curve many  = chain;
    const Curve piece = Polyline{.Vertices = {chain.Vertices[8192], chain.Vertices[8193]}};
    const std::optional<Place> there = Located(many, 8192.5, At(piece, 1.0 - 0x1p-45));
    CHECK(there.value_or(Place{}).Along == 8193.0);
    CHECK(there.value_or(Place{}).At    == Point{.X = 8193, .Y = 0});
  }

  TEST_CASE("crossing: nothing for two bulges or a bulge and an arc, told before their carriers")
  {
    /**
     *  The exact way crosses a bulge's circle with a line alone: two bulges, or a bulge and an arc,
     *  have no crossing it works out, wherever they cross, nor a part with no carrier.
     */
    const Curve one   = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                                              {.Position = {.X = 2, .Y = 0}, .Bulge = 0}},
                                 .Closed   = false};
    const Curve other = Polyline{.Vertices = {{.Position = {.X = 1, .Y = 1}, .Bulge = -1},
                                              {.Position = {.X = 1, .Y = -2}, .Bulge = 0}},
                                 .Closed   = false};
    const Curve round =
        Arc{.Center = {.X = 1, .Y = 0}, .Radius = 1.5, .Start = {}, .Sweep = Degrees(360)};
    const Point near = {.X = 1.7, .Y = -0.7};
    CHECK(!Crossing({&one, &other}, {0.5, 0.5}, near));
    CHECK(!Crossing({&round, &one}, {0.5, 0.5}, near));
    CHECK(!Crossing({&one, &round}, {0.5, 0.5}, near));
    CHECK(!Crossing({&one, &other}, {0.5, std::numeric_limits<double>::infinity()}, near));
  }

  TEST_CASE("crossing: aimed as crossed, an arc's turn kept where its angle places it")
  {
    /**
     *  Two arcs, an arc and a line: the point and places Crossing gives, and on each arc placed by
     *  its angle the direction whose point is At's and whose quarter turn is Derivative's, to the
     *  bit.
     */
    const Curve round =
        Arc{.Center = {.X = 0, .Y = 0}, .Radius = 2, .Start = Degrees(10), .Sweep = Degrees(160)};
    const Curve across = Arc{
        .Center = {.X = 0.3, .Y = 2.5}, .Radius = 2, .Start = Degrees(190), .Sweep = Degrees(160)};
    const Curve cut   = Line{.Start = {.X = -3, .Y = 1.1}, .End = {.X = 3, .Y = 1.3}};
    std::size_t aimed = 0;
    for (const auto &[one, other] : {std::pair{&round, &across}, std::pair{&cut, &round}})
    {
      const std::array<const Curve *, 2> curves = {one, other};
      const Roots                        roots  = Quickly(curves);
      for (const Point near : {Point{.X = 1.6, .Y = 1.2}, Point{.X = -1.4, .Y = 1.3}})
      {
        const Approach                approach = {.Along = {0.5, 0.5}, .At = near};
        const std::optional<Crossed>  crossed  = Crossing(curves, approach, roots);
        const std::optional<Bearings> bearings = Aiming(curves, approach, roots);
        CHECK(crossed.has_value() == bearings.has_value());
        const Crossed  made = crossed.value_or(Crossed{});
        const Bearings kept = bearings.value_or(Bearings{});
        CHECK(made.At    == kept.At);
        CHECK(made.Along == kept.Along);
        for (std::size_t side = 0; side < 2; ++side)
        {
          const auto *arc = std::get_if<Arc>(curves.at(side));
          CHECK((arc != nullptr) == kept.Ways.at(side).has_value());
          if (arc != nullptr && kept.Ways.at(side))
          {
            ++aimed;
            const Vector way   = *kept.Ways.at(side);
            const double along = kept.Along.at(side);
            CHECK(arc->Center + (way * arc->Radius) == At(*arc, along));
            CHECK(Perpendicular(way) * (arc->Radius * arc->Sweep.InRadians())
                  == Derivative(*arc, along));
          }
        }
      }
    }
    CHECK(aimed >= 4);

    /**
     *  A crossing off an arc's sweep, and none at all: nothing, as Crossing gives.
     */
    const Curve brief =
        Arc{.Center = {.X = 0, .Y = 0}, .Radius = 2, .Start = Degrees(10), .Sweep = Degrees(20)};
    const Curve apart =
        Arc{.Center = {.X = 9, .Y = 0}, .Radius = 1, .Start = {}, .Sweep = Degrees(90)};
    const Approach far = {.Along = {0.5, 0.5}, .At = {.X = -1.4, .Y = 1.3}};
    CHECK(!Aiming({&brief, &across}, far, Quickly({&brief, &across})));
    CHECK(!Aiming({&round, &apart}, far, Quickly({&round, &apart})));
  }

  TEST_CASE("crossing: no place on a polyline's segment of no length, bulged or not")
  {
    /**
     *  A bulge from a vertex to itself is no arc: Nth makes it a line of no length, which places no
     *  point, as a straight segment of no length does.
     */
    const Curve doubled = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                {.Position = {.X = 1, .Y = 0.5}, .Bulge = 0.5},
                                                {.Position = {.X = 1, .Y = 0.5}},
                                                {.Position = {.X = 2, .Y = 0}}}};
    CHECK_FALSE(Located(doubled, 1.5, {.X = 1, .Y = 0.5}));
    CHECK_FALSE(Located(doubled, 1.0, {.X = 1.25, .Y = 0.5}));
  }

  TEST_CASE("crossing: a line's place at its ends, and a hair past one, as Nearest takes it")
  {
    const Curve   ground = Line{.Start = {.X = 1, .Y = 0}, .End = {.X = 3, .Y = 0}};
    const Curve   atEnd  = Line{.Start = {.X = 3, .Y = -1}, .End = {.X = 3, .Y = 1}};
    const Crossed ended  = Found(Crossing({&ground, &atEnd}, {1.0, 0.5}, Point{.X = 3, .Y = 0}));
    CHECK(ended.At       == Point{.X = 3, .Y = 0});
    CHECK(ended.Along[0] == 1.0);
    CHECK(ended.On[0]    == Point{.X = 3, .Y = 0});

    /**
     *  A crossing a step of a double before the start, within the slack a rounded point is given:
     *  placed at the start, its nearest point.
     */
    const double  hair   = std::nextafter(1.0, 0.0);
    const Curve   before = Line{.Start = {.X = hair, .Y = -1}, .End = {.X = hair, .Y = 1}};
    const Crossed early = Found(Crossing({&ground, &before}, {0.0, 0.5}, Point{.X = hair, .Y = 0}));
    CHECK(early.At       == Point{.X = hair, .Y = 0});
    CHECK(early.Along[0] == 0.0);
    CHECK(early.On[0]    == Point{.X = 1, .Y = 0});
  }

  TEST_CASE("crossing: a line and a conic, two circles, by their exact roots")
  {
    const Arc circle = Circle({}, 5);
    /**
     *  3-4-5: the crossing is a double, found so.
     */
    const Crossed pythagorean = Found(
        Near(Line{.Start = {.X = -8, .Y = 4}, .End = {.X = 8, .Y = 4}}, circle, {.X = 3, .Y = 4}));
    CHECK(pythagorean.At == Point{.X = 3, .Y = 4});
    /**
     *  sqrt(24): the double nearest.
     */
    const Crossed rooted = Found(Near(Line{.Start = {.X = -8, .Y = 1}, .End = {.X = 8, .Y = 1}},
                                      circle, {.X = 4.9, .Y = 1}));
    CHECK(rooted.At.X == std::sqrt(24.0));
    CHECK(
        !Near(Line{.Start = {.X = -8, .Y = 6}, .End = {.X = 8, .Y = 6}}, circle, {.X = 0, .Y = 6}));
    CHECK(!Near(Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 0, .Y = 0}}, circle, {}));

    /**
     *  An ellipse, its axes 4 and 2 turned a quarter: 4 (x/2)^2 + (y/4)^2 = ... crossing x = 1.
     */
    const Ellipse oval = {.Map = Rotation(Degrees(90)) * Scaling(4, 2), .Start = {}, .Sweep = Turn};
    const Crossed oblong = Found(
        Near(Line{.Start = {.X = 1, .Y = -5}, .End = {.X = 1, .Y = 5}}, oval, {.X = 1, .Y = 3}));
    CHECK(oblong.At.X == 1.0);
    CHECK(oblong.At.Y == std::sqrt(12.0));

    /**
     *  Two circles 6 apart, radii 5: at (3, +-4).
     */
    const Crossed two = Found(Near(circle, Circle({.X = 6, .Y = 0}, 5), {.X = 3, .Y = -3.9}));
    CHECK(two.At == Point{.X = 3, .Y = -4});
    CHECK(!Near(circle, Circle({}, 3), {.X = 4, .Y = 0}));
    CHECK(!Near(circle, Circle({.X = 20, .Y = 0}, 5), {.X = 10, .Y = 0}));
    CHECK(!Near(oval, circle, {.X = 2, .Y = 2}));

    /**
     *  A touch asked for: found where the carriers touch exactly, as y = 5 does the circle; nowhere
     *  where they cross, as y = 4 does, though no double-double roots rule it out first.
     */
    const Curve    round   = circle;
    const Curve    tangent = Line{.Start = {.X = -8, .Y = 5}, .End = {.X = 8, .Y = 5}};
    const Curve    across  = Line{.Start = {.X = -8, .Y = 4}, .End = {.X = 8, .Y = 4}};
    const Approach touch   = {.Along = {0.5, 0.25}, .At = {.X = 0, .Y = 5}, .Touch = true};
    const Crossed  kissed  = Found(Crossing({&tangent, &round}, touch, Roots{}));
    CHECK(kissed.At == Point{.X = 0, .Y = 5});
    CHECK(!Crossing({&across, &round}, touch, Roots{}));
  }

  TEST_CASE("crossing: a touch decided by the circles' own sums, or left to the exact crossing")
  {
    const Curve round   = Arc{.Center = {}, .Radius = 5, .Start = {}, .Sweep = Turn};
    const Curve top     = Line{.Start = {.X = -8, .Y = 5}, .End = {.X = 8, .Y = 5}};
    const Curve across  = Line{.Start = {.X = -8, .Y = 4}, .End = {.X = 8, .Y = 4}};
    const Curve ray     = Ray{.Start = {.X = 9, .Y = -5}, .Way = {.X = -1, .Y = 0}};
    const Curve endless = Xline{.Through = {.X = 5, .Y = 9}, .Way = {.X = 0, .Y = 1}};

    /**
     *  A line, a ray, an endless line touching: decided, at the foot; one crossing: decided,
     *  nowhere.
     */
    CHECK(Tangency({&top, &round}).At     == Point{.X = 0, .Y = 5});
    CHECK(Tangency({&round, &ray}).At     == Point{.X = 0, .Y = -5});
    CHECK(Tangency({&endless, &round}).At == Point{.X = 5, .Y = 0});
    const Curve inexact = Line{.Start = {.X = 0.1, .Y = 5}, .End = {.X = 0.7, .Y = 5}};
    const Curve offset = Arc{.Center = {.X = 0.4, .Y = 0}, .Radius = 5, .Start = {}, .Sweep = Turn};
    CHECK(Tangency({&inexact, &offset}).At == Point{.X = 0.4, .Y = 5});
    CHECK(Tangency({&across, &round}).Decided);
    CHECK(!Tangency({&across, &round}).At);

    /**
     *  Circles touching outside and inside, or crossing: decided.
     */
    const Curve outside = Arc{.Center = {.X = 8, .Y = 0}, .Radius = 3, .Start = {}, .Sweep = Turn};
    const Curve inside  = Arc{.Center = {.X = 2, .Y = 0}, .Radius = 3, .Start = {}, .Sweep = Turn};
    const Curve crossed = Arc{.Center = {.X = 6, .Y = 0}, .Radius = 3, .Start = {}, .Sweep = Turn};
    CHECK(Tangency({&round, &outside}).At == Point{.X = 5, .Y = 0});
    CHECK(Tangency({&inside, &round}).At  == Point{.X = 5, .Y = 0});
    CHECK(Tangency({&round, &crossed}).Decided);
    CHECK(!Tangency({&round, &crossed}).At);

    /**
     *  Left undecided: a line of no length, circles on one center, other curves; a touch a hair
     *  from a halfway point between doubles, or past the range the bound holds over.
     */
    const Curve dot  = Line{.Start = {.X = 5, .Y = 0}, .End = {.X = 5, .Y = 0}};
    const Curve same = Arc{.Center = {}, .Radius = 3, .Start = {}, .Sweep = Turn};
    const Curve oval = Ellipse{.Map = Scaling(2, 1), .Start = {}, .Sweep = Turn};
    CHECK(!Tangency({&dot, &round}).Decided);
    CHECK(!Tangency({&round, &same}).Decided);
    CHECK(!Tangency({&top, &across}).Decided);
    CHECK(!Tangency({&round, &oval}).Decided);
    const double tiny = 0x1p-53;
    const Curve  left = Arc{.Center = {.X = 1, .Y = 0}, .Radius = tiny, .Start = {}, .Sweep = Turn};
    const Curve  right =
        Arc{.Center = {.X = 1 + (2 * tiny), .Y = 0}, .Radius = tiny, .Start = {}, .Sweep = Turn};
    CHECK(!Tangency({&left, &right}).Decided);
    const double huge = 0x1p950;
    const Curve  far  = Arc{.Center = {.X = huge, .Y = 0}, .Radius = 1, .Start = {}, .Sweep = Turn};
    const Curve  next = Arc{.Center = {.X = huge, .Y = 2}, .Radius = 1, .Start = {}, .Sweep = Turn};
    CHECK(!Tangency({&far, &next}).Decided);
    const double small = 0x1p-950;
    const Curve  low = Arc{.Center = {.X = small, .Y = 0}, .Radius = 1, .Start = {}, .Sweep = Turn};
    const Curve  above =
        Arc{.Center = {.X = small, .Y = 2}, .Radius = 1, .Start = {}, .Sweep = Turn};
    CHECK(!Tangency({&low, &above}).Decided);
    const Curve under = Arc{.Center = {.X = 0, .Y = 1}, .Radius = tiny, .Start = {}, .Sweep = Turn};
    const Curve over =
        Arc{.Center = {.X = 0, .Y = 1 + (2 * tiny)}, .Radius = tiny, .Start = {}, .Sweep = Turn};
    CHECK(!Tangency({&under, &over}).Decided);

    /**
     *  Undecided, the exact crossing decides: a touch where the point sits on a halfway point, and
     *  none where a line crosses a circle spelled as an ellipse.
     */
    const Approach touch = {.Along = {0.0, 0.0}, .At = {.X = 1, .Y = 0}, .Touch = true};
    CHECK(Crossing({&left, &right}, touch, Roots{}).has_value());
    const Curve circled = Ellipse{.Map = Scaling(5, 5), .Start = {}, .Sweep = Turn};
    CHECK(!Crossing({&across, &circled}, touch, Roots{}));
  }

  TEST_CASE("crossing: by the curves' own short sums, the exact crossing, where hard")
  {
    /**
     *  Crossings the double-double roots leave open: on an axis, at small whole numbers, and near
     *  the origin on circles a million wide. Where the short sums decide, the same double as the
     *  exact way gives.
     */
    Mixer mixer;
    int   differ  = 0;
    int   decided = 0;
    for (int round = 0; round < 20000; ++round)
    {
      const bool   wide  = round % 3 == 0;
      const double scale = wide ? 1e6 : 1.0;
      const auto   whole = [&mixer](double reach) -> double
      {
        return std::round(mixer.Unit() * reach);
      };
      const Point  centre = wide ? Point{.X = -scale * 0.7, .Y = scale * 0.7}
                                 : Point{.X = whole(8.0), .Y = whole(8.0)};
      const double radius =
          wide ? Distance(centre, {.X = mixer.Unit(), .Y = 0}) : 1.0 + std::abs(whole(6.0));
      const Arc   circle = {.Center = centre, .Radius = radius, .Start = {}, .Sweep = Turn};
      const Curve ring   = circle;
      const Point near   = {.X = mixer.Unit() * radius * 0.1, .Y = mixer.Unit() * radius * 0.1};
      if (round % 2 == 0)
      {
        const Line     line     = {.Start = {.X = whole(9.0) - 10.0, .Y = whole(9.0)},
                                   .End   = {.X = whole(9.0) + 10.0, .Y = whole(9.0)}};
        const Curve    straight = line;
        const Touching got      = Shortly({&straight, &ring}, near);
        if (got.Decided)
        {
          ++decided;
          differ += *got.At == ExactCrossing(line, circle, near) ? 0 : 1;
        }
      }
      else
      {
        const Arc      other  = {.Center = {.X = centre.X + whole(9.0), .Y = centre.Y + whole(9.0)},
                                 .Radius = 1.0 + std::abs(whole(9.0)),
                                 .Start  = {},
                                 .Sweep  = Turn};
        const Curve    second = other;
        const Touching got    = Shortly({&ring, &second}, near);
        if (got.Decided)
        {
          ++decided;
          differ += *got.At == ExactCrossing(circle, other, near) ? 0 : 1;
        }
      }
    }
    CHECK(differ == 0);
    CHECK(decided > 2000);
    MESSAGE("decided by short sums: " << decided);
    const Curve dot   = Line{.Start = {.X = 1, .Y = 1}, .End = {.X = 1, .Y = 1}};
    const Curve oval  = Ellipse{.Map = Scaling(2, 1), .Start = {}, .Sweep = Turn};
    const Curve round = Arc{.Center = {}, .Radius = 5, .Start = {}, .Sweep = Turn};
    CHECK(!Shortly({&dot, &round}, {}).Decided);
    CHECK(!Shortly({&oval, &round}, {}).Decided);
    CHECK(!Shortly({&round, &dot}, {}).Decided);
    CHECK(!Shortly({&dot, &dot}, {}).Decided);

    /**
     *  Undecided: a crossing on a halfway point between doubles, x = 1 + 2^-53; circles past the
     *  range the bounds hold over.
     */
    const Curve axis = Line{.Start = {.X = -3, .Y = 0}, .End = {.X = 3, .Y = 0}};
    const Curve speck =
        Arc{.Center = {.X = 1, .Y = 0}, .Radius = 0x1p-53, .Start = {}, .Sweep = Turn};
    CHECK(!Shortly({&axis, &speck}, {.X = 1, .Y = 0}).Decided);
    const Curve vast =
        Arc{.Center = {.X = 0x1p120, .Y = 0}, .Radius = 0x1p119, .Start = {}, .Sweep = Turn};
    const Curve vaster =
        Arc{.Center = {.X = 0x1p120, .Y = 0x1p119}, .Radius = 0x1p119, .Start = {}, .Sweep = Turn};
    CHECK(!Shortly({&vast, &vaster}, {}).Decided);
  }

  TEST_CASE("crossing: on an arc, its ends within a few steps, past them not")
  {
    const Curve   quarter = Arc{.Center = {}, .Radius = 5, .Start = {}, .Sweep = Degrees(90)};
    const Curve   across  = Line{.Start = {.X = -8, .Y = 0}, .End = {.X = 8, .Y = 0}};
    const Crossed start = Found(Crossing({&quarter, &across}, {0.0, 13.0 / 16}, {.X = 5, .Y = 0}));
    CHECK(start.At == Point{.X = 5, .Y = 0});
    CHECK(!Crossing({&quarter, &across}, {0.0, 3.0 / 16}, {.X = -5, .Y = 0}));
  }

  TEST_CASE("crossing: two lines in double-double where it settles, exactly where it does not")
  {
    /**
     *  Random lines at every size and slant, their carriers' crossing the same double as the exact
     *  way gives, however it was reached.
     */
    Mixer mixer;
    int   differ = 0;
    for (int round = 0; round < 20000; ++round)
    {
      const double size        = std::ldexp(1.0, static_cast<int>(40 * mixer.Unit()));
      const double farFraction = mixer.Unit();
      const double far         = std::ldexp(farFraction, static_cast<int>(30 * (mixer.Unit() + 1)));
      const double slantFraction = mixer.Unit();
      const double slant = std::ldexp(slantFraction, -static_cast<int>(20 * (mixer.Unit() + 1)));
      const Line   one   = {
          .Start = {.X = far + (size * mixer.Unit()), .Y = far + (size * mixer.Unit())},
          .End   = {.X = far + (size * mixer.Unit()), .Y = far + (size * mixer.Unit())}};
      const Vector way    = one.End - one.Start;
      const Point  middle = one.Start + (way * 0.5);
      const Line other = {.Start = middle - (way * 0.7) + (Vector{.X = -way.Y, .Y = way.X} * slant),
                          .End = middle + (way * 0.6) - (Vector{.X = -way.Y, .Y = way.X} * slant)};
      const Curve                  first   = one;
      const Curve                  second  = other;
      const std::optional<Crossed> crossed = Crossing({&first, &second}, {0.5, 0.5}, middle);
      if (crossed && !(crossed->At == ExactCrossing(one, other)))
      {
        ++differ;
      }
    }
    CHECK(differ == 0);

    /**
     *  A crossing on the halfway point between two doubles, and a hair either side: the quick way
     *  cannot settle them, and the exact way rounds each to the nearest, a tie to the even.
     */
    const Line slope = {.Start = {.X = 1, .Y = 0}, .End = {.X = 2, .Y = 1}};
    for (const double hair : {0.0, 0x1p-110, -0x1p-110})
    {
      const double height = 0x1p-53 + hair;
      const Line   level  = {.Start = {.X = 0, .Y = height}, .End = {.X = 3, .Y = height}};
      const Curve  first  = slope;
      const Curve  second = level;
      const std::optional<Crossed> crossed =
          Crossing({&first, &second}, {0.0, 0.33}, {.X = 1, .Y = 0});
      REQUIRE(crossed.has_value());
      CHECK(crossed.value_or(Crossed{}).At == ExactCrossing(slope, level));
    }
    const Curve tie    = slope;
    const Curve halved = Line{.Start = {.X = 0, .Y = 0x1p-53}, .End = {.X = 3, .Y = 0x1p-53}};
    CHECK(Crossing({&tie, &halved}, {0.0, 0.33}, {.X = 1, .Y = 0}).value_or(Crossed{}).At.X == 1.0);
  }

  /**
   *  Two lines sharing an end exactly cross there, whichever is first, met plain or made ready: an
   *  end 2^-38 off the axis beside numbers of 2^17, where the exact numerator's parts, summed as
   *  they come, cancel to a sliver of the largest and its three largest put the point a step off;
   *  and an end three steps of the least double off the axis beside ends at 1, whose products fall
   *  below the least double, as do the halfway points the point is rounded by.
   */
  TEST_CASE("crossing: lines sharing an end exactly cross at it, met plain or ready, either first")
  {
    struct Joint
    {
      Point Start  = {}; /**< Where the first line starts.          */
      Point Shared = {}; /**< Where it ends, and the second starts. */
      Point End    = {}; /**< Where the second ends.                */
    };

    const std::array joints = {
        Joint{.Start  = {.X = 0x1.08e9ebb34f082p+17, .Y = 0x1.4268d0b61bc2p+15},
              .Shared = {.X = 0x1.517f718f5ec06p+16, .Y = 0x1.cff2e210c47bep-39},
              .End    = {.X = 0x1.2e1b87aedf062p+16, .Y = -0x1.0abe8bb04a16ep+15}},
        Joint{.Start  = {.X = 1, .Y = -1},
              .Shared = {.X = 1, .Y = 0x3p-1074},
              .End    = {.X = 1.5, .Y = 1}}};
    for (const Joint &joint : joints)
    {
      const Curve into = Line{.Start = joint.Start, .End = joint.Shared};
      const Curve out  = Line{.Start = joint.Shared, .End = joint.End};
      CHECK(Found(Crossing({&into, &out}, {1.0, 0.0}, joint.Shared)).At == joint.Shared);
      CHECK(Found(Crossing({&out, &into}, {0.0, 1.0}, joint.Shared)).At == joint.Shared);
      const Prepared inward{into};
      const Prepared outward{out};
      for (const double length : {1e-10, 0.0})
      {
        const Tolerance tolerance = {.Length = length, .Radians = 1e-12};
        for (const bool swapped : {false, true})
        {
          const Curve         &first  = swapped ? out : into;
          const Curve         &second = swapped ? into : out;
          const Prepared      &one    = swapped ? outward : inward;
          const Prepared      &other  = swapped ? inward : outward;
          const Intersections  plain  = Intersect(first, second, tolerance);
          Intersections        ready;
          Marking              marked;
          const Prepared::Pair pair = {one, other};
          Intersect(pair, tolerance, ready);
          Intersect(pair, tolerance, marked);
          for (const Intersections *meeting : std::array<const Intersections *, 2>{&plain, &ready})
          {
            REQUIRE(meeting->Hits.size() == 1);
            const Hit &hit = meeting->Hits[0];
            CHECK(hit.At           == joint.Shared);
            CHECK(hit.First.At     == joint.Shared);
            CHECK(hit.Second.At    == joint.Shared);
            CHECK(hit.First.Along  == (swapped ? 0.0 : 1.0));
            CHECK(hit.Second.Along == (swapped ? 1.0 : 0.0));
          }
          REQUIRE(marked.Marks.size() == 1);
          CHECK(marked.Marks[0].At    == joint.Shared);
        }
      }
    }

    /**
     *  Lines whose numbers part too far are crossed in whole numbers, and never touch; lines of no
     *  numbers, or of one not finite, are left to the sums as before, which find no crossing.
     */
    const Joint   &last    = joints.back();
    const Curve    into    = Line{.Start = last.Start, .End = last.Shared};
    const Curve    out     = Line{.Start = last.Shared, .End = last.End};
    const Approach touch   = {.Along = {1.0, 0.0}, .At = last.Shared, .Touch = true};
    const Curve    nothing = Line{};
    const Curve    endless =
        Line{.Start = {}, .End = {.X = std::numeric_limits<double>::infinity(), .Y = 1}};
    CHECK(!Crossing({&into, &out}, touch, Roots{}).has_value());
    CHECK(!Crossing({&nothing, &nothing}, {0.0, 0.0}, Point{}).has_value());
    CHECK(!Crossing({&endless, &out}, {0.5, 0.5}, last.Shared).has_value());
  }

  /**
   *  Lines some five hundred binades below 1, whose products fall near the subnormals, where their
   *  roundings are steps of the least double and no share of a bound of their sizes: the crossing
   *  in double-double settles none the whole numbers put at another double, where it settled one in
   *  two hundred a step off before its bounds took those roundings in, and the exact crossing is
   *  the whole numbers' double.
   */
  TEST_CASE("crossing: lines whose products near the subnormals settle only where exact")
  {
    Mixer mixer;
    int   crossed = 0;
    int   wrong   = 0;
    for (int round = 0; round < Rounds(4000); ++round)
    {
      const int  power = -505 - static_cast<int>(20.0 * std::abs(mixer.Unit()));
      const auto place = [&mixer, power] -> Point
      {
        return {.X = std::ldexp(mixer.Unit(), power), .Y = std::ldexp(mixer.Unit(), power)};
      };
      const Line                 one    = {.Start = place(), .End = place()};
      const Line                 other  = {.Start = place(), .End = place()};
      const Curve                first  = one;
      const Curve                second = other;
      const std::optional<Point> whole =
          Exact::CrossWhole({Exact::Straightway{.Through = one.Start, .Toward = one.End},
                             Exact::Straightway{.Through = other.Start, .Toward = other.End}});
      const Point                exact = whole.value_or(Point{});
      const Roots                roots = Quickly({&first, &second});
      const std::optional<Point> spot =
          Spotted({&first, &second}, {.Along = {0.5, 0.5}, .At = exact}, 0);
      wrong   += static_cast<int>(roots.At[0].has_value() && !(*roots.At[0] == exact));
      wrong   += static_cast<int>(spot.has_value() && !(*spot == exact));
      crossed += static_cast<int>(spot.has_value());
    }
    CHECK(wrong == 0);
    CHECK(crossed > Rounds(4000) / 2);
  }

  TEST_CASE("crossing: lines and circles, circles and circles in double-double, as exactly")
  {
    /**
     *  Random lines and circles, circles and circles, at every size, crossing at every angle: the
     *  crossing nearer the hit the same double as the exact way gives, however reached. The hit is
     *  a foot from a line's start, or a point anywhere on the first circle, never on the line of
     *  centers, where both crossings lie as near and rounding alone would choose.
     */
    Mixer mixer;
    int   differ  = 0;
    int   crossed = 0;
    for (int round = 0; round < 20000; ++round)
    {
      const double size        = std::ldexp(1.0, static_cast<int>(30 * mixer.Unit()));
      const double farFraction = mixer.Unit();
      const double far         = std::ldexp(farFraction, static_cast<int>(25 * (mixer.Unit() + 1)));
      const Point  centre = {.X = far + (size * mixer.Unit()), .Y = far + (size * mixer.Unit())};
      const Arc    circle = {.Center = centre,
                             .Radius = size * (0.5 + std::abs(mixer.Unit())),
                             .Start  = Radians(Pi * mixer.Unit()),
                             .Sweep  = Turn};
      const Curve  ring   = circle;
      if (round % 2 == 0)
      {
        const Line line = {
            .Start = centre + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size * 3.0),
            .End   = centre + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size * 3.0)};
        const Curve straight = line;
        for (const Place &foot : Perpendiculars(ring, line.Start))
        {
          const Point                  near = foot.At;
          const std::optional<Crossed> got =
              Crossing({&straight, &ring}, {Nearest(straight, near).Along, foot.Along}, near);
          if (got)
          {
            ++crossed;
            differ += got->At == ExactCrossing(line, circle, near) ? 0 : 1;
          }
        }
      }
      else
      {
        const Arc other = {.Center = centre + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size),
                           .Radius = size * (0.5 + std::abs(mixer.Unit())),
                           .Start  = Radians(Pi * mixer.Unit()),
                           .Sweep  = Turn};
        const Curve                  beside = other;
        const double                 along  = (mixer.Unit() + 1.0) / 2.0;
        const Point                  near   = At(ring, along);
        const std::optional<Crossed> got =
            Crossing({&ring, &beside}, {along, Nearest(beside, near).Along}, near);
        if (got)
        {
          ++crossed;
          differ += got->At == ExactCrossing(circle, other, near) ? 0 : 1;
        }
      }
    }
    MESSAGE("crossings checked: " << crossed);
    CHECK(differ == 0);
  }

  TEST_CASE("crossing: one of a line's two crossings with a circle, worked alone, as of the two")
  {
    /**
     *  Either of the two, the same point to the bit as when both are worked; nothing for a line
     *  that misses the circle, grazes it, or has no length.
     */
    const Arc circle  = {.Center = {.X = 1.25, .Y = -0.5}, .Radius = 5, .Start = {}, .Sweep = Turn};
    const Curve round = circle;
    const Line  chord = {.Start = {.X = -9, .Y = 2.5}, .End = {.X = 11, .Y = 3.75}};
    const Curve line  = chord;
    const Roots both  = Quickly({&line, &round});
    REQUIRE(both.Twice);
    REQUIRE(both.At[0].has_value());
    REQUIRE(both.At[1].has_value());
    CHECK(Quickly(chord, circle, true)  == both.At[0]);
    CHECK(Quickly(chord, circle, false) == both.At[1]);
    const Line outside = {.Start = {.X = -9, .Y = 7}, .End = {.X = 11, .Y = 7}};
    const Line grazing = {.Start = {.X = -9, .Y = 4.5}, .End = {.X = 11, .Y = 4.5}};
    const Line nothing = {.Start = {.X = 1, .Y = 1}, .End = {.X = 1, .Y = 1}};
    CHECK(!Quickly(outside, circle, true));
    CHECK(!Quickly(grazing, circle, false));
    CHECK(!Quickly(nothing, circle, true));
  }
} // namespace Geoxy::Crossings
