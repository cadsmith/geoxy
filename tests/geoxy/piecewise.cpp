/***************************************************************************************************
 *  @file      tests/geoxy/piecewise.cpp
 *  @brief     Tests that a polyline met segment by segment answers as the thorough way does, to the
 *             bit, crossings on its arcs too, from eight vertices to five hundred and twelve; and
 *             that whatever ties one segment's hits to another's is left to the thorough way.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "piecewise.hpp"

#include "draws.hpp"
#include "polylines.hpp"
#include "thorough.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <utility>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  Two curves met the thorough way.
     *
     *  @param  one       One curve.
     *  @param  other     The other.
     *  @param  tolerance The tolerance.
     *  @return           What they meet in.
     */
    auto PiecewiseThoroughly(const Curve &one, const Curve &other, const Tolerance &tolerance)
        -> Intersections
    {
      Intersections meeting;
      Thorough({one, other}, tolerance, meeting);
      return meeting;
    }

    /**
     *  Whether two curves are met piecewise, at the tolerance Intersect hands its direct ways.
     *
     *  @param  first     The first curve.
     *  @param  second    The second.
     *  @param  tolerance The tolerance asked for.
     *  @return           Whether they are.
     */
    auto MetPiecewise(const Curve &first, const Curve &second, const Tolerance &tolerance) -> bool
    {
      const std::optional<double> length = Reachable({first, second}, tolerance);
      Intersections               meeting;
      return length.has_value()
          && Piecewise({.Shapes = {&first, &second}},
                       {.Length = *length, .Radians = tolerance.Radians}, meeting);
    }

    /**
     *  Whether two meetings are the same to the bit, hit for hit and stretch for stretch.
     *
     *  @param  one   One.
     *  @param  other The other.
     *  @return       Whether they are.
     */
    auto PiecewiseIdentical(const Intersections &one, const Intersections &other) -> bool
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

    /**
     *  Whether two curves meet as they meet thoroughly, to the bit.
     *
     *  @param  first     The first curve.
     *  @param  second    The second.
     *  @param  tolerance The tolerance.
     *  @return           Whether they do.
     */
    auto AsThorough(const Curve &first, const Curve &second, const Tolerance &tolerance) -> bool
    {
      return PiecewiseIdentical(Intersect(first, second, tolerance),
                                PiecewiseThoroughly(first, second, tolerance));
    }

    /**
     *  A polyline with a tail of nine more segments in steps from its last vertex, so that the
     *  swift kernel, asked first only of a polyline of a few segments, leaves it to the legs.
     *
     *  @param  polyline The polyline.
     *  @param  step     Each step of the tail.
     *  @return          It with its tail.
     */
    auto Tailed(Polyline polyline, Vector step) -> Polyline
    {
      for (int index = 0; index < 9; ++index)
      {
        polyline.Vertices.push_back({.Position = polyline.Vertices.back().Position + step});
      }
      return polyline;
    }

    /**
     *  A polyline with its coordinates swapped: a zigzag across turned to zigzag up.
     *
     *  @param  polyline The polyline.
     *  @return          It turned.
     */
    auto Upright(Polyline polyline) -> Polyline
    {
      for (Vertex &vertex : polyline.Vertices)
      {
        vertex.Position = {.X = vertex.Position.Y, .Y = vertex.Position.X};
        vertex.Bulge    = -vertex.Bulge;
      }
      return polyline;
    }

    /**
     *  A curve to meet a polyline with, drawn at random: a line, one through two of its vertices, a
     *  ray, an endless line, an arc, a polyline, the same polyline moved a hair, a zigzag, a line
     *  across, or a polyline through some of its vertices.
     *
     *  @param  draws    The draws.
     *  @param  polyline The polyline.
     *  @return          The curve.
     */
    auto Other(Draws &draws, const Polyline &polyline) -> Curve
    {
      const std::size_t count = polyline.Vertices.size();
      const auto        point = [&draws]() -> Point
      {
        return {.X = draws.Near(0.5, 1), .Y = draws.Near(0.5, 1)};
      };
      const auto vertex = [&draws, &polyline, count]() -> Point
      {
        return polyline.Vertices[draws.Next() % count].Position;
      };
      switch (draws.Next() % 10)
      {
      case 0:
        return Line{.Start = vertex(), .End = vertex()};
      case 1:
        return Ray{.Start = point(), .Way = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)}};
      case 2:
        return Xline{.Through = point(), .Way = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)}};
      case 3:
        return Arc{.Center = point(),
                   .Radius = draws.Near(0.6, 0.5),
                   .Start  = Degrees(draws.Near(0, 180)),
                   .Sweep  = Degrees(draws.Near(0, 300))};
      case 4:
        return Strung(draws, 2 + (draws.Next() % 10));
      case 5:
      {
        Polyline moved = polyline;
        for (Vertex &corner : moved.Vertices)
        {
          corner.Position.X += draws.Near(0, 1e-9);
        }
        return moved;
      }
      case 6:
        return Upright(Zigzagged(draws, 2 + (draws.Next() % 40)));
      case 7:
      {
        Polyline shared = {.Vertices = {{.Position = vertex(), .Bulge = Bulged(draws)},
                                        {.Position = vertex(), .Bulge = Bulged(draws)},
                                        {.Position = point(), .Bulge = 0}}};
        return shared;
      }
      default:
      {
        const double high = draws.Near(0.5, 0.3);
        return Line{.Start = {.X = -0.2, .Y = high},
                    .End   = {.X = 1.2, .Y = high + draws.Near(0, 0.2)}};
      }
      }
    }
  } // namespace

  namespace
  {
    /**
     *  Whether two curves met piecewise, where they are, meet as they meet thoroughly, to the bit.
     *
     *  @param  first     The first curve.
     *  @param  second    The second.
     *  @param  tolerance The tolerance.
     *  @return           Whether they are met piecewise, and whether they then meet so.
     */
    auto Matched(const Curve &first, const Curve &second, const Tolerance &tolerance)
        -> std::pair<bool, bool>
    {
      const std::optional<double> length = Reachable({first, second}, tolerance);
      Intersections               meeting;
      const bool decided = length.has_value()
                        && Piecewise({.Shapes = {&first, &second}},
                                     {.Length = *length, .Radians = tolerance.Radians}, meeting);
      return {decided,
              !decided
                  || PiecewiseIdentical(meeting, PiecewiseThoroughly(first, second, tolerance))};
    }

    /**
     *  A point of a grid of whole numbers from 0 to 8.
     *
     *  @param  draws The draws.
     *  @return       The point.
     */
    auto Lattice(Draws &draws) -> Point
    {
      return {.X = static_cast<double>(draws.Next() % 9),
              .Y = static_cast<double>(draws.Next() % 9)};
    }

    /**
     *  A curve moved: scaled by a power of two and moved by an offset, each coordinate exactly
     *  where the offset is a whole number and the coordinates are.
     *
     *  @param  curve The curve, a line or a polyline.
     *  @param  power The power of two.
     *  @param  shift The offset, along x and against y.
     *  @return       The curve moved.
     */
    auto Moved(Curve curve, int power, double shift) -> Curve
    {
      const auto move = [power, shift](Point point) -> Point
      {
        return {.X = std::ldexp(point.X, power) + shift, .Y = std::ldexp(point.Y, power) - shift};
      };
      if (auto *polyline = std::get_if<Polyline>(&curve))
      {
        for (Vertex &vertex : polyline->Vertices)
        {
          vertex.Position = move(vertex.Position);
        }
      }
      if (auto *line = std::get_if<Line>(&curve))
      {
        *line = {.Start = move(line->Start), .End = move(line->End)};
      }
      return curve;
    }

    /**
     *  A curve to meet a polyline of the grid with at its vertices and on its segments: a line from
     *  one of its vertices, a line on the grid, a polyline through some of its vertices, one that
     *  starts or ends halfway along one of its segments, or a line through two of its vertices and
     *  on past them.
     *
     *  @param  draws    The draws.
     *  @param  polyline The polyline.
     *  @return          The curve.
     */
    auto Abutter(Draws &draws, const Polyline &polyline) -> Curve
    {
      const std::span<const Vertex> vertices = polyline.Vertices;
      const auto                    vertex   = [&draws, vertices]() -> Point
      {
        return vertices[draws.Next() % vertices.size()].Position;
      };
      switch (draws.Next() % 5)
      {
      case 0:
        return Line{.Start = vertex(), .End = Lattice(draws)};
      case 1:
        return Line{.Start = Lattice(draws), .End = Lattice(draws)};
      case 2:
      {
        Polyline through = {.Closed = draws.Next() % 4 == 0};
        for (std::size_t count = 2 + (draws.Next() % 5); count > 0; --count)
        {
          through.Vertices.push_back({.Position = draws.Next() % 2 == 0 ? vertex() : Lattice(draws),
                                      .Bulge    = draws.Next() % 8 == 0 ? -0.2 : 0.0});
        }
        return through;
      }
      case 3:
      {
        const std::size_t index = draws.Next() % (vertices.size() - 1);
        const Point       half  = Midpoint(vertices[index].Position, vertices[index + 1].Position);
        Polyline          tee   = {.Vertices = {{.Position = half}, {.Position = Lattice(draws)}}};
        if (draws.Next() % 2 == 0)
        {
          tee.Vertices.push_back({.Position = Lattice(draws)});
        }
        if (draws.Next() % 2 == 0)
        {
          std::ranges::reverse(tee.Vertices);
        }
        return tee;
      }
      default:
      {
        const Point one   = vertex();
        const Point other = vertex();
        return Line{.Start = one + (one - other), .End = other + (other - one)};
      }
      }
    }
  } // namespace

  TEST_CASE("piecewise: a line across zigzags of 8 to 512 vertices, as met thoroughly")
  {
    /**
     *  A line crossing every segment of a zigzag, a third of them arcs: met segment by segment at
     *  every size, either way round, every crossing the thorough way's to the bit.
     */
    Draws           draws{0x5EED'2015};
    const Tolerance tolerance = Tolerance::Of({.Min = {.X = 0, .Y = 0}, .Max = {.X = 1, .Y = 1}});
    for (const std::size_t count : {8UZ, 16UZ, 32UZ, 64UZ, 128UZ, 256UZ, 512UZ})
    {
      CAPTURE(count);
      for (int round = 0; round < Rounds(4); ++round)
      {
        const Curve zigzag = Zigzagged(draws, count);
        const Curve line   = Line{.Start = {.X = -0.1, .Y = draws.Near(0.5, 0.1)},
                                  .End   = {.X = 1.1, .Y = draws.Near(0.5, 0.1)}};
        for (const auto &[one, other] : {std::pair{&line, &zigzag}, std::pair{&zigzag, &line}})
        {
          CHECK(MetPiecewise(*one, *other, tolerance));
          CHECK(AsThorough(*one, *other, tolerance));
        }
      }
    }
  }

  TEST_CASE("piecewise: zigzags across zigzags of 8 to 64 vertices, as met thoroughly")
  {
    Draws           draws{0x5EED'2016};
    const Tolerance tolerance = Tolerance::Of({.Min = {.X = 0, .Y = 0}, .Max = {.X = 1, .Y = 1}});
    for (const std::size_t count : {8UZ, 16UZ, 32UZ, 64UZ})
    {
      CAPTURE(count);
      const Curve across  = Zigzagged(draws, count);
      const Curve upright = Upright(Zigzagged(draws, count));
      CHECK(MetPiecewise(across, upright, tolerance));
      CHECK(AsThorough(across, upright, tolerance));
      CHECK(AsThorough(upright, across, tolerance));
    }
  }

  TEST_CASE("piecewise: polylines of every kind and any curve met with them, as met thoroughly")
  {
    /**
     *  Vertices on lines and on each other, segments of no length, arcs past a half turn, closed
     *  polylines, curves through vertices and along segments: met piecewise where nothing ties one
     *  segment's hits to another's, and thoroughly where something does, the same either way.
     */
    Draws       draws{0x5EED'2017};
    std::size_t decided = 0;
    std::size_t left    = 0;
    std::size_t differ  = 0;
    for (int round = 0; round < Rounds(3000); ++round)
    {
      const Polyline polyline = round % 4 == 0 ? Zigzagged(draws, 2 + (draws.Next() % 40))
                                               : Strung(draws, 1 + (draws.Next() % 12));
      const Curve    strung   = polyline;
      const Curve    other    = Other(draws, polyline);
      Tolerance tolerance     = Tolerance::Of({.Min = {.X = -1, .Y = -1}, .Max = {.X = 2, .Y = 2}});
      if (round % 7 == 0)
      {
        tolerance.Length = 1e-3;
      }
      for (const auto &[one, two] : {std::pair{&strung, &other}, std::pair{&other, &strung}})
      {
        const bool piecewise  = MetPiecewise(*one, *two, tolerance);
        decided              += piecewise ? 1 : 0;
        left                 += piecewise ? 0 : 1;
        differ               += AsThorough(*one, *two, tolerance) ? 0U : 1U;
      }
    }
    MESSAGE("met piecewise: " << decided << ", left to the thorough way: " << left);
    CHECK(differ == 0);
    CHECK(decided > left / 4);
    CHECK(left > decided / 4);
  }

  TEST_CASE("piecewise: legwork kept by a polyline made ready meets as legwork worked out anew")
  {
    /**
     *  Each polyline's legs, arcs and grid kept once, by either curve, both or neither: every pair
     *  decided alike, to the bit.
     */
    Draws       draws{0x5EED'2045};
    std::size_t differ = 0;
    for (int round = 0; round < Rounds(1500); ++round)
    {
      const Polyline polyline = round % 3 == 0 ? Zigzagged(draws, 2 + (draws.Next() % 60))
                                               : Strung(draws, 1 + (draws.Next() % 12));
      const Curve    strung   = polyline;
      const Curve    other =
          round % 5 == 0 ? Curve{Strung(draws, 2 + (draws.Next() % 30))} : Other(draws, polyline);
      const Tolerance tolerance =
          Tolerance::Of({.Min = {.X = -1, .Y = -1}, .Max = {.X = 2, .Y = 2}});
      for (const auto &[one, two] : {std::pair{&strung, &other}, std::pair{&other, &strung}})
      {
        const std::optional<double> length = Reachable({*one, *two}, tolerance);
        if (!length)
        {
          continue;
        }
        const Tolerance usable = {.Length = *length, .Radians = tolerance.Radians};
        const std::shared_ptr<const Legwork> first  = Legworked(*one);
        const std::shared_ptr<const Legwork> second = Legworked(*two);
        Intersections                        anew;
        const bool decided = Piecewise({.Shapes = {one, two}}, usable, anew);
        for (const std::array<const Legwork *, 2> &kept :
             {std::array{first.get(), second.get()},
              std::array<const Legwork *, 2>{first.get(), nullptr},
              std::array<const Legwork *, 2>{nullptr, second.get()}})
        {
          Intersections meeting;
          const bool    same =
              Piecewise({.Shapes = {one, two}, .Legs = kept}, usable, meeting) == decided;
          differ += same && PiecewiseIdentical(meeting, anew) ? 0U : 1U;
        }
      }
    }
    CHECK(differ == 0);
    CHECK(Legworked(Curve{Line{}}) == nullptr);

    /**
     *  An arc keeps its one leg, and one swept a whole turn keeps it unsound, as Legged leaves it:
     *  met piecewise by neither.
     */
    const Curve quarter =
        Arc{.Center = {.X = 0, .Y = 0}, .Radius = 1, .Start = {}, .Sweep = Degrees(90)};
    const Curve whole   = Arc{.Center = {.X = 0, .Y = 0}, .Radius = 1, .Start = {}, .Sweep = Turn};
    const Curve through = Polyline{.Vertices = {{.Position = {.X = -2, .Y = 0.5}},
                                                {.Position = {.X = 0, .Y = 0.4}, .Bulge = 0.2},
                                                {.Position = {.X = 2, .Y = 0.5}}},
                                   .Closed   = false};
    const std::shared_ptr<const Legwork> round = Legworked(quarter);
    const std::shared_ptr<const Legwork> ring  = Legworked(whole);
    CHECK(round != nullptr);
    CHECK(ring != nullptr);
    const Tolerance usable = {.Length = 1e-9, .Radians = 1e-9};
    Intersections   kept;
    Intersections   anew;
    CHECK(Piecewise({.Shapes = {&through, &quarter}, .Legs = {nullptr, round.get()}}, usable, kept)
          == Piecewise({.Shapes = {&through, &quarter}}, usable, anew));
    CHECK(PiecewiseIdentical(kept, anew));
    Intersections circled;
    CHECK(
        !Piecewise({.Shapes = {&through, &whole}, .Legs = {nullptr, ring.get()}}, usable, circled));
  }

  TEST_CASE("piecewise: two arcs crossing twice or once, met by the direct way's own steps")
  {
    /**
     *  Two half circles crossing twice, each order along the first, and an arc of its own across a
     *  bulge: met piecewise, the same as met thoroughly.
     */
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Curve     upper     = Polyline{
        .Vertices = {{.Position = {.X = 2, .Y = 0}, .Bulge = 1}, {.Position = {.X = -2, .Y = 0}}},
        .Closed   = false};
    const Curve lower = Polyline{.Vertices = {{.Position = {.X = -2, .Y = 2.5}, .Bulge = 1},
                                              {.Position = {.X = 2.1, .Y = 2.4}}},
                                 .Closed   = false};
    const Curve round = Arc{
        .Center = {.X = 0.2, .Y = 3}, .Radius = 2, .Start = Degrees(200), .Sweep = Degrees(140)};
    for (const auto &[one, two] : {std::pair{&upper, &lower}, std::pair{&lower, &upper},
                                   std::pair{&upper, &round}, std::pair{&round, &upper}})
    {
      CHECK(AsThorough(*one, *two, tolerance));
    }
    CHECK(MetPiecewise(upper, lower, tolerance));
    CHECK(MetPiecewise(lower, upper, tolerance));
    CHECK(Intersect(upper, lower, tolerance).Hits.size() == 2);

    /**
     *  A bulge of 1e17, its sweep a whole turn as a double, across a gentle one: met piecewise at
     *  the arc's own grain, its one crossing the double nearest the true one, which 120 digits put
     *  at (140000.0000000000022, -3.919972000000000124e-7). The thorough way, held to the vertices'
     *  finer grain, finds none.
     */
    const Curve whole = Polyline{
        .Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1e17}, {.Position = {.X = 1, .Y = 0}}},
        .Closed   = false};
    const Curve gentle = Polyline{.Vertices = {{.Position = {.X = 1e5, .Y = -4e5}, .Bulge = 0.1},
                                               {.Position = {.X = 1e5, .Y = 4e5}}},
                                  .Closed   = false};
    const Point truth  = {.X = 0x1.117p+17, .Y = -0x1.a4e75eaf66af4p-22};
    for (const auto &[one, two] : {std::pair{&whole, &gentle}, std::pair{&gentle, &whole}})
    {
      const Intersections met   = Intersect(*one, *two, tolerance);
      const bool          found = met.Hits.size() == 1 && met.Hits.front().At == truth;
      CHECK(found);
    }
  }

  TEST_CASE("piecewise: polylines' arcs met from their circles alone, as met thoroughly")
  {
    /**
     *  Polylines of a few bulges, half turns and less, past a half turn, flat ones of 1e-4 and ones
     *  swept a whole turn as a double, crossing at every angle, an end now and then near the
     *  other's circle or a crossing near a flat bulge's end, every third pair an arc of its own
     *  across one: met from their circles where all is clear, by the direct way's steps for two
     *  arcs where not, the same as met thoroughly, either first. An arc of its own crossing a bulge
     *  next to an axis, where Newton's method alone polishes them, is left to those steps.
     */
    constexpr std::array<double, 7> Bulges = {0.25, -0.5, 1.0, -2.0, 1e-4, -0.75, 0x1p21};
    Draws                           draws{0xB0E5'0005};
    std::size_t                     decided = 0;
    std::size_t                     differ  = 0;
    for (int round = 0; round < Rounds(4000); ++round)
    {
      const auto bulged = [&draws, &Bulges]() -> Curve
      {
        Polyline   polyline;
        const auto count = 2 + static_cast<std::size_t>(draws.Next() % 3);
        for (std::size_t index = 0; index < count; ++index)
        {
          polyline.Vertices.push_back(
              {.Position = {.X = draws.Near(0.0, 1.0), .Y = draws.Near(0.0, 1.0)},
               .Bulge    = Bulges.at(static_cast<std::size_t>(draws.Next() % Bulges.size()))});
        }
        return polyline;
      };
      const auto arced = [&draws]() -> Curve
      {
        return Arc{.Center = {.X = draws.Near(0.0, 1.0), .Y = draws.Near(0.0, 1.0)},
                   .Radius = draws.Near(0.75, 0.5),
                   .Start  = Degrees(draws.Near(0.0, 180.0)),
                   .Sweep  = Degrees(draws.Near(0.0, 300.0))};
      };
      const Curve     one       = bulged();
      const Curve     other     = round % 3 == 0 ? arced() : bulged();
      const Tolerance tolerance = {.Length = round % 5 == 0 ? 1e-4 : 1e-9, .Radians = 1e-9};
      for (const auto &[first, second] : {std::pair{&one, &other}, std::pair{&other, &one}})
      {
        decided += MetPiecewise(*first, *second, tolerance) ? 1U : 0U;
        differ  += AsThorough(*first, *second, tolerance) ? 0U : 1U;
      }
    }
    MESSAGE("met piecewise: " << decided);
    CHECK(differ == 0);
    CHECK(decided > 0);
    const Curve half = Polyline{
        .Vertices = {{.Position = {.X = 0, .Y = -1}, .Bulge = 1}, {.Position = {.X = 0, .Y = 1}}},
        .Closed   = false};
    const Curve beside = Arc{
        .Center = {.X = 1, .Y = 1.001}, .Radius = 1, .Start = Degrees(-150), .Sweep = Degrees(90)};
    const Tolerance fine = {.Length = 1e-9, .Radians = 1e-9};
    CHECK(AsThorough(half, beside, fine));
    CHECK(AsThorough(beside, half, fine));
    CHECK(Intersect(half, beside, fine).Hits.size() == 1);
  }

  TEST_CASE("piecewise: a polyline sized by its bulges' circles, the thorough way at their grain")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Curve     whole     = Polyline{
        .Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1e17}, {.Position = {.X = 1, .Y = 0}}},
        .Closed   = false};
    const Point truth = {.X = 0x1.117p+17, .Y = -0x1.a4e75eaf66af4p-22};

    /**
     *  The bulge of 1e17 across the gentle one, the gentle one going on through its end: a vertex
     *  on the other curve's end leaves the pair to the thorough way, which, held to the vertices'
     *  grain of 1e-9, found the end alone. Held to the bulge's own grain, 128, it finds the
     *  crossing too, the double nearest the true one, which 400 digits put at
     *  (140000.0000000000022, -3.919972000000000124e-7), and the end, (1, 0), once: the bulge's
     *  start lies within the grain of the gentle one's vertex there, and the hit found at it is put
     *  exactly at the end, where the hit found at the end is put too.
     */
    const Curve tail = Polyline{.Vertices = {{.Position = {.X = 1e5, .Y = -4e5}, .Bulge = 0.1},
                                             {.Position = {.X = 1e5, .Y = 4e5}},
                                             {.Position = {.X = 1, .Y = 0}},
                                             {.Position = {.X = 1, .Y = -1}}},
                                .Closed   = false};
    for (const auto &[one, two] : {std::pair{&whole, &tail}, std::pair{&tail, &whole}})
    {
      CHECK_FALSE(MetPiecewise(*one, *two, tolerance));
      const Intersections met     = Intersect(*one, *two, tolerance);
      const auto          crossed = [&truth](const Hit &hit) -> bool
      {
        return hit.At == truth && hit.Kind == Contact::Cross;
      };
      const auto ended = [](const Hit &hit) -> bool
      {
        return hit.At == Point{.X = 1, .Y = 0} && hit.Kind == Contact::End;
      };
      REQUIRE(met.Hits.size() == 2);
      CHECK(crossed(met.Hits[0]));
      CHECK(ended(met.Hits[1]));
    }

    /**
     *  A straight segment on from the bulge, across the gentle one 3.9e-7 from the bulge's
     *  crossing, nearer than the bulge's grain keeps two hits apart: met the thorough way, each
     *  crossing the double nearest the true one, (140000.0000000000022, 0) the segment's.
     */
    const Curve onward = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1e17},
                                               {.Position = {.X = 1, .Y = 0}},
                                               {.Position = {.X = 2e5, .Y = 0}}},
                                  .Closed   = false};
    const Curve gentle = Polyline{.Vertices = {{.Position = {.X = 1e5, .Y = -4e5}, .Bulge = 0.1},
                                               {.Position = {.X = 1e5, .Y = 4e5}}},
                                  .Closed   = false};
    for (const auto &[one, two] : {std::pair{&onward, &gentle}, std::pair{&gentle, &onward}})
    {
      CHECK_FALSE(MetPiecewise(*one, *two, tolerance));
      const Intersections met = Intersect(*one, *two, tolerance);
      REQUIRE(met.Hits.size() == 2);
      CHECK(met.Hits[0].At    == truth);
      CHECK(met.Hits[1].At    == Point{.X = 0x1.117p+17, .Y = 0});
    }

    /**
     *  A polyline's size: its vertices', and each bulged segment's circle's, its center's grown by
     *  its radius, as Nth's arc has them; the closing segment's of a closed polyline, not the last
     *  vertex's bulge of an open one, nor a bulge under Straight, nor a slight one, under 1e-4,
     *  whose arc is worked out from its ends, nor a segment of no length, nor a circle past every
     *  double.
     */
    const auto circled = [](const Polyline &polyline, std::size_t index) -> double
    {
      const Arc arc = std::get<Arc>(Nth(polyline, index));
      return std::max(std::abs(arc.Center.X), std::abs(arc.Center.Y)) + arc.Radius;
    };
    const Polyline open = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1e-3},
                                        {.Position = {.X = 1, .Y = 0}, .Bulge = 0x1p-61},
                                        {.Position = {.X = 2, .Y = 0}, .Bulge = 1e-4},
                                        {.Position = {.X = 2, .Y = 0}, .Bulge = 0},
                                        {.Position = {.X = 3, .Y = 1}, .Bulge = 1e-5}},
                           .Closed   = false};
    CHECK(std::equal_to<double>{}(Sized(open), circled(open, 0)));
    Polyline shut = open;
    shut.Closed   = true;
    CHECK(std::equal_to<double>{}(Sized(shut), Sized(open)));
    shut.Vertices[4].Bulge = 2e-4;
    CHECK(std::equal_to<double>{}(Sized(shut), circled(shut, 4)));
    const Polyline straight = {.Vertices = {{.Position = {.X = -5, .Y = 0}, .Bulge = 0x1p-61},
                                            {.Position = {.X = 1, .Y = 4}, .Bulge = 1e-3}},
                               .Closed   = false};
    CHECK(std::equal_to<double>{}(Sized(straight), 5.0));
    const Polyline vast = {.Vertices = {{.Position = {.X = -5, .Y = 0}, .Bulge = 1e300},
                                        {.Position = {.X = 1, .Y = 4}}},
                           .Closed   = false};
    CHECK(std::equal_to<double>{}(Sized(vast), 5.0));
  }

  TEST_CASE("piecewise: what ties one segment's hits to another's goes the thorough way")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const auto      left      = [&tolerance](const Curve &one, const Curve &other) -> void
    {
      CHECK_FALSE(MetPiecewise(one, other, tolerance));
      CHECK(AsThorough(one, other, tolerance));
      CHECK(AsThorough(other, one, tolerance));
    };
    const Curve wedge = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                              {.Position = {.X = 1, .Y = 1}},
                                              {.Position = {.X = 2, .Y = 0}}}};
    SUBCASE("a touch at a vertex, or a hit near one but not on it")
    {
      left(wedge, Line{.Start = {.X = 0, .Y = 1}, .End = {.X = 2, .Y = 1}});
      left(wedge, Line{.Start = {.X = 1 + 1e-8, .Y = 0.5}, .End = {.X = 1 + 1e-8, .Y = 1.5}});
    }
    SUBCASE("a stretch along two segments")
    {
      const Curve along = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                {.Position = {.X = 1, .Y = 0}},
                                                {.Position = {.X = 2, .Y = 0}},
                                                {.Position = {.X = 3, .Y = 1}}}};
      left(along, Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 2.5, .Y = 0}});
      CHECK(Intersect(along, Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 2.5, .Y = 0}}, tolerance)
                .Overlaps.size()
            == 1);
    }
    SUBCASE("a segment so short and so square to a line across it that its shadow seems shared")
    {
      /**
       *  Six hundred tolerances long, across a line a five hundred and twelfth of a radian from
       *  square: its shadow on the line is longer than the tolerance and both ends of the shadow
       *  lie within it of the segment's carrier, though its vertices lie three hundred away. The
       *  thorough way crosses them, the shadow not a stretch they share.
       */
      const Curve tick = Polyline{.Vertices = {{.Position = {.X = 0.5, .Y = -3e-7}},
                                               {.Position = {.X = 0.5 + 1.17e-9, .Y = 3e-7}},
                                               {.Position = {.X = 0.75, .Y = 1}}}};
      const Curve flat = Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 1, .Y = 0}};
      CHECK_FALSE(MetPiecewise(flat, tick, tolerance));
      CHECK(AsThorough(flat, tick, tolerance));
      const Intersections met = Intersect(flat, tick, tolerance);
      CHECK(met.Overlaps.empty());
      CHECK(std::ranges::count_if(met.Hits,
                                  [](const Hit &hit) -> bool { return hit.Kind == Contact::Cross; })
            == 1);
    }
    SUBCASE("a closed polyline's seam")
    {
      const Curve square = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                 {.Position = {.X = 2, .Y = 0}},
                                                 {.Position = {.X = 2, .Y = 2}},
                                                 {.Position = {.X = 0, .Y = 2}}},
                                    .Closed   = true};
      left(square, Line{.Start = {.X = -1, .Y = -1}, .End = {.X = 1, .Y = 1}});
    }
    SUBCASE("a segment of no length on the other curve")
    {
      const Curve doubled = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                  {.Position = {.X = 1, .Y = 0.5}},
                                                  {.Position = {.X = 1, .Y = 0.5}},
                                                  {.Position = {.X = 2, .Y = 0}}}};
      left(doubled, Line{.Start = {.X = 0, .Y = 0.5}, .End = {.X = 2, .Y = 0.5}});
    }
    SUBCASE("a shallow crossing")
    {
      left(Tailed(Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                        {.Position = {.X = 10, .Y = 1}},
                                        {.Position = {.X = 10, .Y = 5}}}},
                  {.X = 0, .Y = 1}),
           Line{.Start = {.X = 0, .Y = 0.5}, .End = {.X = 9, .Y = 0.5}});
    }
    SUBCASE("an arc touching a line, or another arc")
    {
      const Curve dome = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = -1},
                                               {.Position = {.X = 2, .Y = 0}}}};
      left(dome, Line{.Start = {.X = -1, .Y = 1}, .End = {.X = 3, .Y = 1}});
      left(dome, Arc{.Center = {.X = 1, .Y = 3},
                     .Radius = 2,
                     .Start  = Degrees(200),
                     .Sweep  = Degrees(140)});
    }
    SUBCASE("arcs of one circle, near each other or far apart")
    {
      const double sweep = 20.0;
      const Curve  top   = Polyline{
          .Vertices = {{.Position = {.X = 10 * std::cos(80.0 * std::numbers::pi / 180),
                                     .Y = 10 * std::sin(80.0 * std::numbers::pi / 180)},
                        .Bulge    = std::tan(sweep * std::numbers::pi / 720)},
                       {.Position = {.X = 10 * std::cos(100.0 * std::numbers::pi / 180),
                                     .Y = 10 * std::sin(100.0 * std::numbers::pi / 180)}}}};
      left(top, Arc{.Center = {}, .Radius = 10, .Start = Degrees(260), .Sweep = Degrees(20)});
      left(top, Arc{.Center = {}, .Radius = 10, .Start = Degrees(85), .Sweep = Degrees(40)});
    }
    SUBCASE("two crossings near each other")
    {
      const Curve cross = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                {.Position = {.X = 2, .Y = 2}},
                                                {.Position = {.X = 2, .Y = 0}},
                                                {.Position = {.X = 0, .Y = 2}}}};
      left(cross, Line{.Start = {.X = 0, .Y = 1 + 1e-12}, .End = {.X = 2, .Y = 1 + 1e-12}});
    }
    SUBCASE("hits on short segments with long ones between, on both polylines")
    {
      /**
       *  Where the thorough way polished each crossing to moves the place it looks at between them
       *  along a segment five thousand times faster: no look is certain.
       */
      const Curve comb = Polyline{.Vertices = {{.Position = {.X = 0.25, .Y = -0.01}},
                                               {.Position = {.X = 0.25, .Y = 0.01}},
                                               {.Position = {.X = 0.75, .Y = 100}},
                                               {.Position = {.X = 1.25, .Y = 0.01}},
                                               {.Position = {.X = 1.25, .Y = -0.01}}}};
      const Curve rake = Polyline{.Vertices = {{.Position = {.X = 0.24, .Y = 0}},
                                               {.Position = {.X = 0.26, .Y = 0}},
                                               {.Position = {.X = 0.75, .Y = -100}},
                                               {.Position = {.X = 1.24, .Y = 0}},
                                               {.Position = {.X = 1.26, .Y = 0}}}};
      CHECK(Intersect(comb, rake, tolerance).Hits.size() == 2);
      left(comb, rake);
    }
  }

  TEST_CASE("piecewise: arcs too tight and pairs the direct ways leave, at coarse tolerances")
  {
    /**
     *  A shallow crossing first, at a sine of 0.29; a steep one; an arc past a half turn of a
     *  radius of a fifth, its chord three tenths from the line; and a tail away from it.
     */
    const Curve zig =
        Tailed(Polyline{.Vertices = {{.Position = {.X = -0.5, .Y = 0.7}},
                                     {.Position = {.X = 1.5, .Y = 1.3}},
                                     {.Position = {.X = 2.5, .Y = 4}},
                                     {.Position = {.X = 3.5, .Y = 0.7}, .Bulge = -7.87},
                                     {.Position = {.X = 3.6, .Y = 0.7}},
                                     {.Position = {.X = 5, .Y = 0.7}}}},
               {.X = 0, .Y = -1});
    const Curve across = Line{.Start = {.X = -1, .Y = 1}, .End = {.X = 6, .Y = 1}};
    CHECK(MetPiecewise(zig, across, {.Length = 1e-9, .Radians = 1e-9}));
    /**
     *  A tolerance of a thousandth makes the arc's radius under the clearance; an angle of half a
     *  radian makes the shallow crossing a touch, and one of a fifth leaves it short of twice its
     *  sine.
     */
    for (const Tolerance tolerance :
         {Tolerance{.Length = 1e-3, .Radians = 1e-9}, Tolerance{.Length = 1e-9, .Radians = 0.5},
          Tolerance{.Length = 1e-9, .Radians = 0.2}})
    {
      CAPTURE(tolerance.Length);
      CAPTURE(tolerance.Radians);
      CHECK_FALSE(MetPiecewise(zig, across, tolerance));
      CHECK(AsThorough(zig, across, tolerance));
    }
  }

  TEST_CASE("piecewise: closed polylines, crossings on the axes and arcs either side")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const auto      met       = [&tolerance](const Curve &one, const Curve &other) -> void
    {
      CHECK(MetPiecewise(one, other, tolerance));
      CHECK(AsThorough(one, other, tolerance));
    };
    const Curve square = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                               {.Position = {.X = 2, .Y = 0}},
                                               {.Position = {.X = 2, .Y = 2}},
                                               {.Position = {.X = 0, .Y = 2}}},
                                  .Closed   = true};
    /**
     *  Two crossings half the square apart either way round, which the thorough way's look between
     *  may take either way; two more than half apart, looked at across the seam.
     */
    met(square, Line{.Start = {.X = 1, .Y = -1}, .End = {.X = 1, .Y = 3}});
    met(Line{.Start = {.X = -0.3, .Y = 0.6}, .End = {.X = 1.8, .Y = -0.1}}, square);
    met(square, Line{.Start = {.X = -0.3, .Y = 0.6}, .End = {.X = 1.8, .Y = -0.1}});
    const Curve zigzag = Polyline{.Vertices = {{.Position = {.X = -1, .Y = 0.5}},
                                               {.Position = {.X = 0.5, .Y = 2.5}},
                                               {.Position = {.X = 1, .Y = -0.5}},
                                               {.Position = {.X = 1.5, .Y = 2.5}},
                                               {.Position = {.X = 3, .Y = 1.5}}}};
    met(square, zigzag);
    met(zigzag, square);
    /**
     *  Crossings out of order along each polyline: the look between two of them along one lands on
     *  a third.
     */
    const Curve bar = Polyline{
        .Vertices = {{.Position = {.X = 0, .Y = 0.25}}, {.Position = {.X = 4, .Y = 0.25}}}};
    const Curve back = Polyline{.Vertices = {{.Position = {.X = 0.8, .Y = -1}},
                                             {.Position = {.X = 1.2, .Y = 1.5}},
                                             {.Position = {.X = 2.8, .Y = 1.5}},
                                             {.Position = {.X = 3.2, .Y = -1}},
                                             {.Position = {.X = 2.2, .Y = -1}},
                                             {.Position = {.X = 1.8, .Y = 1.5}}}};
    met(bar, back);
    /**
     *  Arcs on the first polyline, on the second, on both, and as the other curve; rays and endless
     *  lines.
     */
    const Curve bowed = Polyline{.Vertices = {{.Position = {.X = -1, .Y = 0.3}, .Bulge = 0.4},
                                              {.Position = {.X = 0.5, .Y = 2.7}, .Bulge = -0.3},
                                              {.Position = {.X = 1, .Y = -0.7}, .Bulge = 0.2},
                                              {.Position = {.X = 1.5, .Y = 2.6}}}};
    met(bowed, square);
    met(square, bowed);
    met(bowed, Upright(std::get<Polyline>(bowed)));
    met(bowed, Arc{.Center = {.X = 0.3, .Y = 1},
                   .Radius = 1.7,
                   .Start  = Degrees(-60),
                   .Sweep  = Degrees(120)});
    met(Ray{.Start = {.X = -2, .Y = 1}, .Way = {.X = 1, .Y = 0.1}}, bowed);
    met(bowed, Xline{.Through = {.X = 0, .Y = 1}, .Way = {.X = 1, .Y = -0.2}});
  }

  TEST_CASE("piecewise: hits at vertices and ends lying exactly on the other curve")
  {
    /**
     *  Polylines of the grid, a sixth of their segments bulged, met with lines and polylines
     *  through their vertices, along their segments and ending on them: scaled by 2^±20 and 2^±40,
     *  or moved 1e8 off, every coordinate still exact. A hit where both curves pass through a
     *  vertex, or where one ends on the other, is met piecewise where nothing else ties it, the
     *  thorough way's to the bit; a touch there, a bulge beside it or a closed polyline's seam is
     *  left to the thorough way.
     */
    Draws       draws{0x5EED'2036};
    std::size_t decided = 0;
    std::size_t differ  = 0;
    for (int round = 0; round < Rounds(4000); ++round)
    {
      Polyline polyline = {.Closed = false};
      for (std::size_t count = 2 + (draws.Next() % 6); count > 0; --count)
      {
        polyline.Vertices.push_back(
            {.Position = Lattice(draws), .Bulge = draws.Next() % 6 == 0 ? 0.3 : 0.0});
      }
      polyline.Closed                = polyline.Vertices.size() >= 3 && draws.Next() % 4 == 0;
      const Curve          abutting  = Abutter(draws, polyline);
      constexpr std::array Powers    = {0, 0, 20, -20, 40, -40};
      const int            power     = Powers.at(draws.Next() % Powers.size());
      const double         shift     = power == 0 && draws.Next() % 5 == 0 ? 1e8 : 0.0;
      const Curve          strung    = Moved(polyline, power, shift);
      const Curve          other     = Moved(abutting, power, shift);
      const double         size      = std::ldexp(8.0, power) + shift;
      const Tolerance      tolerance = {
          .Length  = std::max(std::ldexp(1e-9, power), 64.0 * (std::nextafter(size, 1e300) - size)),
          .Radians = 1e-9};
      for (const auto &[one, two] : {std::pair{&strung, &other}, std::pair{&other, &strung}})
      {
        const auto [piecewise, same]  = Matched(*one, *two, tolerance);
        decided                      += piecewise ? 1 : 0;
        differ                       += same ? 0 : 1;
      }
    }
    MESSAGE("met piecewise: " << decided << " of " << 2 * Rounds(4000));
    CHECK(differ == 0);
    CHECK(decided > static_cast<std::size_t>(Rounds(4000)) / 2);
  }

  TEST_CASE("piecewise: a crossing at a vertex, an end on a segment, a vertex shared")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const auto      met = [&tolerance](const Curve &one, const Curve &other, Contact kind) -> void
    {
      CHECK(MetPiecewise(one, other, tolerance));
      CHECK(MetPiecewise(other, one, tolerance));
      CHECK(AsThorough(one, other, tolerance));
      CHECK(AsThorough(other, one, tolerance));
      CHECK(std::ranges::contains(Intersect(one, other, tolerance).Hits, kind, &Hit::Kind));
    };
    /**
     *  A zigzag across a line at two of its vertices, from below to above and back.
     */
    const Curve zigzag = Polyline{.Vertices = {{.Position = {.X = 0, .Y = -1}},
                                               {.Position = {.X = 1, .Y = 0}},
                                               {.Position = {.X = 1.5, .Y = 1}},
                                               {.Position = {.X = 3, .Y = 0}},
                                               {.Position = {.X = 4, .Y = -0.5}}}};
    met(zigzag, Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 5, .Y = 0}}, Contact::Cross);
    /**
     *  A line ending on a segment, and a polyline beginning on one.
     */
    met(zigzag, Line{.Start = {.X = 2.25, .Y = 0.5}, .End = {.X = 2.25, .Y = 3}}, Contact::End);
    met(zigzag,
        Polyline{
            .Vertices = {{.Position = {.X = 3.5, .Y = 2}}, {.Position = {.X = 3.5, .Y = -0.25}}}},
        Contact::End);
    /**
     *  Two polylines crossing at a vertex both have, and two that go on from each other's end.
     */
    met(zigzag,
        Polyline{.Vertices = {{.Position = {.X = 1, .Y = 3}},
                              {.Position = {.X = 1.5, .Y = 1}},
                              {.Position = {.X = 2, .Y = -1}}}},
        Contact::Cross);
    met(zigzag,
        Polyline{.Vertices = {{.Position = {.X = 4, .Y = -0.5}}, {.Position = {.X = 5, .Y = 2}}}},
        Contact::End);
    /**
     *  A vertex on the other's vertex where the other turns back, which may be a touch; a closed
     *  polyline's seam on the other; a bulge beside the vertex; and a vertex the polyline passes
     *  twice, or ends on itself at: each left to the thorough way.
     */
    const auto left = [&tolerance](const Curve &one, const Curve &other) -> void
    {
      CHECK_FALSE(MetPiecewise(one, other, tolerance));
      CHECK(AsThorough(one, other, tolerance));
      CHECK(AsThorough(other, one, tolerance));
    };
    left(zigzag, Polyline{.Vertices = {{.Position = {.X = 0.5, .Y = 2}},
                                       {.Position = {.X = 1.5, .Y = 1}},
                                       {.Position = {.X = 2.5, .Y = 2}}}});
    left(Polyline{.Vertices = {{.Position = {.X = 1, .Y = 0}},
                               {.Position = {.X = 1, .Y = 2}},
                               {.Position = {.X = 2, .Y = -1}}},
                  .Closed   = true},
         Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 5, .Y = 0}});
    left(Polyline{.Vertices = {{.Position = {.X = 0, .Y = -1}, .Bulge = 0.2},
                               {.Position = {.X = 1, .Y = 0}},
                               {.Position = {.X = 2, .Y = 1}}}},
         Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 5, .Y = 0}});
    left(Polyline{.Vertices = {{.Position = {.X = 0, .Y = -1}},
                               {.Position = {.X = 1, .Y = 0}},
                               {.Position = {.X = 2, .Y = 1}},
                               {.Position = {.X = 2, .Y = -1}},
                               {.Position = {.X = 1, .Y = 0}},
                               {.Position = {.X = 0, .Y = 1}}}},
         Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 5, .Y = 0}});
    left(Polyline{.Vertices = {{.Position = {.X = 0, .Y = -1}},
                               {.Position = {.X = 2, .Y = 1}},
                               {.Position = {.X = 2, .Y = -1}},
                               {.Position = {.X = 1, .Y = 0}}}},
         Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 5, .Y = 0}});
  }

  TEST_CASE("piecewise: joints left, bulges crossed on an axis, grids of many legs")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const auto      left      = [&tolerance](const Curve &one, const Curve &other) -> void
    {
      CHECK_FALSE(MetPiecewise(one, other, tolerance));
      CHECK(AsThorough(one, other, tolerance));
      CHECK(AsThorough(other, one, tolerance));
    };
    const Curve flat = Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 1, .Y = 0}};
    /**
     *  A vertex on the line's carrier a little past either end of it; a polyline ending on itself
     *  where a line starts, and one crossing itself there.
     */
    left(flat, Polyline{.Vertices = {{.Position = {.X = 1 + 2e-9, .Y = 0}},
                                     {.Position = {.X = 1 + 2e-9, .Y = 1}}}});
    left(flat, Polyline{.Vertices = {{.Position = {.X = -2e-9, .Y = 0}},
                                     {.Position = {.X = -2e-9, .Y = 1}}}});
    const Curve upward = Line{.Start = {.X = 1, .Y = 0}, .End = {.X = 1, .Y = 3}};
    left(Polyline{.Vertices = {{.Position = {.X = 0, .Y = -1}},
                               {.Position = {.X = 2, .Y = 1}},
                               {.Position = {.X = 2, .Y = -1}},
                               {.Position = {.X = 1, .Y = 0}}}},
         upward);
    left(Polyline{.Vertices = {{.Position = {.X = 0, .Y = -1}},
                               {.Position = {.X = 2, .Y = 1}},
                               {.Position = {.X = 2, .Y = -1}},
                               {.Position = {.X = 0, .Y = 1}}}},
         upward);
    /**
     *  Thirty vertices at one point, its legs filed in a grid of one cell.
     */
    Polyline point;
    point.Vertices.resize(30);
    left(point, Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 1, .Y = 0}});
    /**
     *  A bulge crossed by a line along the x axis, each crossing on it, left by the swift kernel to
     *  the direct way; and one crossed shallowly.
     */
    const Curve bow = Polyline{.Vertices = {{.Position = {.X = -2, .Y = -1}, .Bulge = -1.5},
                                            {.Position = {.X = 2, .Y = -1}}}};
    const Curve axis =
        Polyline{.Vertices = {{.Position = {.X = -3, .Y = 0}}, {.Position = {.X = 3, .Y = 0}}}};
    CHECK(MetPiecewise(bow, axis, tolerance));
    CHECK(AsThorough(bow, axis, tolerance));
    CHECK(AsThorough(axis, bow, tolerance));
    CHECK(Intersect(bow, axis, tolerance).Hits.size() == 2);
    const Curve cup = Polyline{
        .Vertices = {{.Position = {.X = -1, .Y = 0}, .Bulge = 1}, {.Position = {.X = 1, .Y = 0}}}};
    left(cup, Polyline{.Vertices = {{.Position = {.X = -2, .Y = -0.995}},
                                    {.Position = {.X = 2, .Y = -0.995}}}});
    /**
     *  Crossings out of order along a polyline of thirty legs, its grid asked about a look along
     *  the other that lands on a third crossing.
     */
    Polyline back = {.Vertices = {{.Position = {.X = 0.8, .Y = -1}},
                                  {.Position = {.X = 1.2, .Y = 1.5}},
                                  {.Position = {.X = 2.8, .Y = 1.5}},
                                  {.Position = {.X = 3.2, .Y = -1}},
                                  {.Position = {.X = 2.2, .Y = -1}},
                                  {.Position = {.X = 1.8, .Y = 1.5}}}};
    for (int step = 1; step <= 25; ++step)
    {
      back.Vertices.push_back({.Position = {.X = 1.8, .Y = 1.5 + step}});
    }
    const Curve bar = Polyline{
        .Vertices = {{.Position = {.X = 0, .Y = 0.25}}, {.Position = {.X = 4, .Y = 0.25}}}};
    CHECK(MetPiecewise(back, bar, tolerance));
    CHECK(AsThorough(back, bar, tolerance));
    CHECK(Intersect(back, bar, tolerance).Hits.size() == 3);
  }

  TEST_CASE("piecewise: curves it does not meet")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Curve     strand    = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                      {.Position = {.X = 1, .Y = 2}},
                                                      {.Position = {.X = 2, .Y = 0}}}};
    const double    far       = std::numeric_limits<double>::infinity();
    for (const Curve &other :
         {Curve{Ellipse{.Map = {.X = {.X = 1, .Y = 0}, .Y = {.X = 0, .Y = 0.5}}}},
          Curve{Line{.Start = {.X = 0, .Y = 1}, .End = {.X = 0, .Y = 1}}},
          Curve{Line{.Start = {.X = 0, .Y = 1}, .End = {.X = far, .Y = 1}}},
          Curve{Arc{.Center = {.X = 1, .Y = 1}, .Radius = 0.5, .Sweep = Degrees(400)}},
          Curve{Arc{.Center = {.X = 1, .Y = 1}, .Radius = 0.5, .Sweep = Degrees(0)}},
          Curve{Arc{.Center = {.X = 1, .Y = 1}, .Radius = -0.5, .Sweep = Degrees(90)}},
          Curve{Polyline{.Vertices = {{.Position = {.X = 1, .Y = 1}}}}},
          Curve{Polyline{
              .Vertices = {{.Position = {.X = 0, .Y = 0}}, {.Position = {.X = 0, .Y = 0}}}}},
          Curve{Polyline{.Vertices = {{.Position = {.X = 0, .Y = 1}, .Bulge = 1e300},
                                      {.Position = {.X = 1, .Y = 1}}}}},
          Curve{Polyline{
              .Vertices = {{.Position = {.X = 0, .Y = 1}}, {.Position = {.X = far, .Y = 1}}}}}})
    {
      Intersections                        meeting;
      const std::shared_ptr<const Legwork> kept = Legworked(other);
      CHECK_FALSE(Piecewise({.Shapes = {&strand, &other}}, tolerance, meeting));
      CHECK_FALSE(Piecewise({.Shapes = {&other, &strand}}, tolerance, meeting));
      CHECK_FALSE(Piecewise({.Shapes = {&strand, &other}, .Legs = {nullptr, kept.get()}}, tolerance,
                            meeting));
      CHECK_FALSE(Piecewise({.Shapes = {&other, &strand}, .Legs = {kept.get(), nullptr}}, tolerance,
                            meeting));
      CHECK(meeting.Hits.empty());
    }
    Intersections meeting;
    const Curve   plain = Line{};
    CHECK_FALSE(Piecewise({.Shapes = {&plain, &plain}}, tolerance, meeting));
  }
} // namespace Geoxy
