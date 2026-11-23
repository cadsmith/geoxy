/***************************************************************************************************
 *  @file      tests/geoxy/slanted.cpp
 *  @brief     Tests that ready lines at any angle, an end of one near the other, meet as the curves
 *             themselves do, to the bit, signed zeros too, and that the slanted pass tells them.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "abut.hpp"
#include "draws.hpp"
#include "plain.hpp"
#include "real.hpp"
#include "sane.hpp"
#include "thorough.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/point.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <utility>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether two points are the same doubles, a zero's sign too.
     *
     *  @param  lhs One point.
     *  @param  rhs The other.
     *  @return     Whether they are.
     */
    auto Bitwise(Point lhs, Point rhs) -> bool
    {
      return std::bit_cast<std::uint64_t>(lhs.X) == std::bit_cast<std::uint64_t>(rhs.X)
          && std::bit_cast<std::uint64_t>(lhs.Y) == std::bit_cast<std::uint64_t>(rhs.Y);
    }

    /**
     *  Whether marks are a meeting's hits, to the bit.
     *
     *  @param  meeting The meeting.
     *  @param  marks   The marks.
     *  @return         Whether they are.
     */
    auto Marked(const Intersections &meeting, std::span<const Mark> marks) -> bool
    {
      if (meeting.Hits.size() != marks.size())
      {
        return false;
      }
      for (std::size_t index = 0; index < marks.size(); ++index)
      {
        const Hit &hit = meeting.Hits[index];
        if (!Bitwise(hit.At, marks[index].At) || hit.Kind != marks[index].Kind)
        {
          return false;
        }
      }
      return true;
    }

    /**
     *  Whether stretches are a meeting's, to the bit.
     *
     *  @param  meeting   The meeting.
     *  @param  stretches The stretches.
     *  @return           Whether they are.
     */
    auto Shared(const Intersections &meeting, std::span<const Overlap> stretches) -> bool
    {
      const auto bits = [](double lhs, double rhs) -> bool
      {
        return std::bit_cast<std::uint64_t>(lhs) == std::bit_cast<std::uint64_t>(rhs);
      };
      const auto same = [&bits](const Overlap &lhs, const Overlap &rhs) -> bool
      {
        return bits(lhs.First.Start, rhs.First.Start)
            && bits(lhs.First.End, rhs.First.End)
            && bits(lhs.Second.Start, rhs.Second.Start)
            && bits(lhs.Second.End, rhs.Second.End)
            && lhs.Along == rhs.Along;
      };
      return std::ranges::equal(meeting.Overlaps, stretches, same);
    }

    /**
     *  How two lines met: whether made ready, through Plainly, and as the slanted pass tells them,
     *  they met as the curves themselves, either first; and how many orders the pass told.
     */
    struct Sloped
    {
      bool Same = true; /**< Whether every way agrees to the bit. */
      int  Told = 0;    /**< How many orders the pass told.       */
    };

    /**
     *  Two lines met one way round, the curves themselves, made ready, through Plainly, and by the
     *  slanted pass asked directly.
     *
     *  @param  curves    The two lines, the first first.
     *  @param  tolerance The tolerance.
     *  @return           How they met.
     */
    auto SlopeOnce(const Curves &curves, const Tolerance &tolerance) -> Sloped
    {
      const Prepared one{curves[0].get()};
      const Prepared other{curves[1].get()};
      Intersections  meeting;
      Marking        marking;
      Marking        plainly;
      Intersect(curves, tolerance, meeting);
      Intersect({one, other}, tolerance, marking);
      const Tolerance sane  = Sane(tolerance);
      const Plain     plain = Plainly({one, other}, sane, plainly);
      const bool      clean = meeting.Overlaps.empty() && !meeting.Spent;
      const Mark      mark  = {.At = plain.At, .Kind = plain.Kind};
      Sloped met = {.Same = Marked(meeting, marking.Marks)
                         && Shared(meeting, marking.Overlaps)
                         && marking.Spent == meeting.Spent
                         && (plain.Told ? clean && Marked(meeting, std::span{&mark, plain.Count})
                                        : Marked(meeting, plainly.Marks)
                                              && Shared(meeting, plainly.Overlaps)),
                    .Told = 0};
      const std::optional<double> length = Reachable(curves, sane);
      if (!length || *length > Coarse || *length < Fine)
      {
        return met;
      }
      const std::optional<Abutment> told =
          Abutting::Slanted({one, other}, {.Length = *length, .Radians = sane.Radians});
      if (told)
      {
        met.Told = 1;
        met.Same = met.Same
                && !meeting.Spent
                && Marked(meeting, std::span{told->Marks}.first(told->Count))
                && Shared(meeting, std::span{&told->Stretch, told->Stretches});
      }
      return met;
    }

    /**
     *  Two lines met either first.
     *
     *  @param  first     One line.
     *  @param  second    The other.
     *  @param  tolerance The tolerance.
     *  @return           How they met.
     */
    auto SlopeBoth(const Line &first, const Line &second, const Tolerance &tolerance) -> Sloped
    {
      const Curve  one   = first;
      const Curve  other = second;
      const Sloped ahead = SlopeOnce({one, other}, tolerance);
      const Sloped back  = SlopeOnce({other, one}, tolerance);
      return {.Same = ahead.Same && back.Same, .Told = ahead.Told + back.Told};
    }

    /**
     *  A way of a length at a bearing.
     *
     *  @param  bearing The bearing, in radians.
     *  @param  reach   How far.
     *  @return         The way.
     */
    auto Heading(double bearing, double reach) -> Vector
    {
      return {.X = reach * std::cos(bearing), .Y = reach * std::sin(bearing)};
    }

    /**
     *  A point a way from another, rounded to a grid where the grid is not zero.
     *
     *  @param  from The point.
     *  @param  way  The way.
     *  @param  grid The grid's step, or zero for none.
     *  @return      The point.
     */
    auto Toward(Point from, Vector way, double grid) -> Point
    {
      const Point point = from + way;
      if (!(grid > 0.0))
      {
        return point;
      }
      return {.X = std::round(point.X / grid) * grid, .Y = std::round(point.Y / grid) * grid};
    }

    /**
     *  The tolerances two lines are met within, of the lines' size.
     */
    constexpr std::array Lengths = {1e-12, 1e-9, 1e-6, 1e-3};
  } // namespace

  TEST_CASE("slanted: polygons' joints at every angle meet as Intersect answers, to the bit")
  {
    Draws draws{0x51A7ED1};
    int   told  = 0;
    int   pairs = 0;
    for (int round = 0; round < Rounds(4000); ++round)
    {
      const double scale = std::ldexp(1.0, static_cast<int>(draws.Next() % 61) - 30);
      const double grid  = draws.Next() % 3 == 0 ? std::ldexp(scale, -12) : 0.0;
      const Point  vertex =
          Toward({.X = draws.Near(0.0, 8.0) * scale, .Y = draws.Near(0.0, 8.0) * scale}, {}, grid);
      const double bearing = draws.Near(0.0, Pi);
      const double turn    = [&draws]() -> double
      {
        if (draws.Next() % 2 == 0)
        {
          return draws.Near(0.0, Pi);
        }
        const double fraction = draws.Near(0.0, 1.0);
        return std::ldexp(fraction, -static_cast<int>(draws.Next() % 30));
      }();
      const Point from = Toward(vertex, Heading(bearing + Pi, draws.Near(2.5, 2.49) * scale), grid);
      const Point onward =
          Toward(vertex, Heading(bearing + turn, draws.Near(2.5, 2.49) * scale), grid);
      const Tolerance tolerance = {.Length  = Lengths.at(draws.Next() % Lengths.size()) * scale,
                                   .Radians = 1e-9};
      const Line      incoming  = {.Start = from, .End = vertex};
      const Line      outgoing  = {.Start = vertex, .End = onward};
      const Sloped    met = draws.Next() % 2 == 0
                              ? SlopeBoth(incoming, outgoing, tolerance)
                              : SlopeBoth({.Start = vertex, .End = from}, outgoing, tolerance);
      CHECK(met.Same);
      told  += met.Told;
      pairs += 2;
    }
    CHECK(told > pairs * 3 / 4);
  }

  TEST_CASE(
      "slanted: ends near the other line, tees and crossings near an end, as Intersect answers")
  {
    Draws draws{0x51A7ED2};
    int   told  = 0;
    int   pairs = 0;
    for (int round = 0; round < Rounds(4000); ++round)
    {
      const double scale  = std::ldexp(1.0, static_cast<int>(draws.Next() % 41) - 20);
      const double grid   = draws.Next() % 3 == 0 ? std::ldexp(scale, -10) : 0.0;
      const double length = Lengths.at(draws.Next() % Lengths.size()) * scale;
      const Point  start =
          Toward({.X = draws.Near(0.0, 8.0) * scale, .Y = draws.Near(0.0, 8.0) * scale}, {}, grid);
      const double bearing = draws.Near(0.0, Pi);
      const Point  end     = Toward(start, Heading(bearing, draws.Near(3.0, 2.9) * scale), grid);
      const double place   = draws.Next() % 4 == 0 ? 1.0 : draws.Near(0.5, 0.5);
      const double asideFraction = draws.Near(0.0, 6.0);
      const double aside         = length * asideFraction * static_cast<double>(draws.Next() % 2);
      const Point  meeting =
          Toward(Lerp(start, end, place), Heading(draws.Near(0.0, Pi), aside), grid);
      const double turn  = draws.Near(0.0, Pi);
      const double past  = draws.Next() % 2 == 0 ? 0.0 : draws.Near(0.0, 1.0) * scale;
      const Point  begin = Toward(meeting, Heading(bearing + turn + Pi, past), grid);
      const Point  finish =
          Toward(meeting, Heading(bearing + turn, draws.Near(2.0, 1.9) * scale), grid);
      const Sloped met = SlopeBoth({.Start = start, .End = end}, {.Start = begin, .End = finish},
                                   {.Length = length, .Radians = 1e-9});
      CHECK(met.Same);
      told  += met.Told;
      pairs += 2;
    }
    CHECK(told > pairs / 3);
  }

  TEST_CASE("slanted: corners, joints, ends end to end and short lines, by hand")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Point     origin    = {.X = 3.0, .Y = 1.0};
    const Point     slant     = {.X = 4.0, .Y = 2.5};
    const Point     onward    = {.X = 5.0, .Y = 4.1};
    const Point     back      = {.X = 6.0, .Y = 2.0};
    /**
     *  A corner at a steep angle, Endly's: an end at the exact crossing.
     */
    const Sloped corner =
        SlopeBoth({.Start = origin, .End = slant}, {.Start = slant, .End = back}, tolerance);
    CHECK(corner.Same);
    CHECK(corner.Told == 2);
    /**
     *  A joint at a narrow angle, Butt's: halfway between the ends.
     */
    const Sloped joint =
        SlopeBoth({.Start = origin, .End = slant}, {.Start = slant, .End = onward}, tolerance);
    CHECK(joint.Same);
    CHECK(joint.Told == 2);
    /**
     *  Two lines on one carrier, end to end: no angle, Butt's.
     */
    const Sloped along =
        SlopeBoth({.Start = {.X = 1.0, .Y = 1.0}, .End = {.X = 2.0, .Y = 3.0}},
                  {.Start = {.X = 2.0, .Y = 3.0}, .End = {.X = 3.0, .Y = 5.0}}, tolerance);
    CHECK(along.Same);
    CHECK(along.Told == 2);
    /**
     *  A tee: an end on the other's middle.
     */
    const Sloped tee =
        SlopeBoth({.Start = origin, .End = back},
                  {.Start = {.X = 4.5, .Y = 1.5}, .End = {.X = 4.0, .Y = 5.0}}, tolerance);
    CHECK(tee.Same);
    CHECK(tee.Told == 2);
    /**
     *  Ends half a tolerance and three tolerances from the other: Endly keeps the first, and gives
     *  the second up to Intersect's thorough way.
     */
    const Point  near = {.X = slant.X + 5e-10, .Y = slant.Y};
    const Sloped half =
        SlopeBoth({.Start = origin, .End = slant}, {.Start = near, .End = back}, tolerance);
    CHECK(half.Same);
    const Point  far = {.X = slant.X + 3e-9, .Y = slant.Y};
    const Sloped three =
        SlopeBoth({.Start = origin, .End = slant}, {.Start = far, .End = back}, tolerance);
    CHECK(three.Same);
    /**
     *  A line shorter than the tolerance at another's end, and one of no length.
     */
    const Sloped tiny =
        SlopeBoth({.Start = origin, .End = slant},
                  {.Start = slant, .End = {.X = slant.X + 1e-10, .Y = slant.Y + 2e-10}}, tolerance);
    CHECK(tiny.Same);
    const Sloped none =
        SlopeBoth({.Start = origin, .End = slant}, {.Start = slant, .End = slant}, tolerance);
    CHECK(none.Same);
    /**
     *  A crossing ten tolerances from an end, every end far: Intersect's quick crossing.
     */
    const Sloped cross =
        SlopeBoth({.Start = origin, .End = slant},
                  {.Start = {.X = 4.0, .Y = 2.0}, .End = {.X = 3.0, .Y = 3.000000005}}, tolerance);
    CHECK(cross.Same);
    /**
     *  Beside: an end on the other's carrier past its end, and lines folded back on each other.
     */
    const Sloped beside =
        SlopeBoth({.Start = origin, .End = slant},
                  {.Start = {.X = 5.0, .Y = 4.0}, .End = {.X = 6.0, .Y = 9.0}}, tolerance);
    CHECK(beside.Same);
    const Sloped folded = SlopeBoth({.Start = origin, .End = slant},
                                    {.Start = slant, .End = {.X = 3.0, .Y = 1.0000001}}, tolerance);
    CHECK(folded.Same);
  }

  TEST_CASE("slanted: lines on one carrier sharing a stretch, or not, as Intersect answers")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Point     origin    = {.X = 1.0, .Y = 2.0};
    const Vector    way       = {.X = 3.0, .Y = 4.0};
    const auto      place     = [origin, way](double along) -> Point
    {
      return origin + (way * along);
    };
    const Line base = {.Start = place(0.0), .End = place(1.0)};
    /**
     *  A stretch shared, one line within the other, the other way round, end to end, and lines
     *  sharing sixteen and three tolerances.
     */
    for (const auto &[from, until] :
         {std::pair{0.5, 1.5}, std::pair{0.25, 0.75}, std::pair{0.75, 0.25}, std::pair{1.0, 2.0},
          std::pair{1.0 - (16e-9 / 5.0), 2.0}, std::pair{1.0 - (3e-9 / 5.0), 2.0},
          std::pair{-0.5, 0.0}, std::pair{1.5, 0.5}})
    {
      CHECK(SlopeBoth(base, {.Start = place(from), .End = place(until)}, tolerance).Same);
    }
    /**
     *  A stretch along the first whose second's end lies two tolerances off it, and a short line
     *  lying along the first with an end a tolerance and a half past the stretch's.
     */
    const Vector across = {.X = -4.0e-9 * 2.0 / 5.0, .Y = 3.0e-9 * 2.0 / 5.0};
    CHECK(SlopeBoth(base, {.Start = place(0.5), .End = place(1.5) + across}, tolerance).Same);
    CHECK(SlopeBoth(base, {.Start = place(0.3) + across, .End = place(0.6) + across}, tolerance)
              .Same);
    CHECK(SlopeBoth(base, {.Start = place(0.3), .End = place(0.6) + (across * 0.75)}, tolerance)
              .Same);
    CHECK(SlopeBoth({.Start = place(0.2), .End = place(0.3)},
                    {.Start = place(0.25), .End = place(0.3 + 3e-10)}, tolerance)
              .Same);
    CHECK(SlopeBoth(base, {.Start = place(0.5), .End = place(1.0)}, tolerance).Told == 2);
  }

  TEST_CASE("slanted: lines whose ways square to nothing, as Intersect answers")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Line      wall      = {.Start = {.X = 0, .Y = 0}, .End = {.X = 3, .Y = 4}};
    const Line      speck     = {.Start = {.X = 0, .Y = 0}, .End = {.X = 1e-170, .Y = 0}};
    const Line      mote      = {.Start = {.X = 0, .Y = 0}, .End = {.X = 0, .Y = 1e-170}};
    const Line      far       = {.Start = {.X = 0, .Y = 0}, .End = {.X = -1e-170, .Y = 2e-170}};
    CHECK(SlopeBoth(wall, speck, tolerance).Same);
    CHECK(SlopeBoth(speck, mote, tolerance).Same);
    CHECK(SlopeBoth(mote, far, tolerance).Same);
    CHECK(SlopeBoth({.Start = {.X = 1, .Y = 1}, .End = {.X = 3, .Y = 4}}, speck, tolerance).Same);
  }

  TEST_CASE("slanted: joints and crossings on an axis keep Intersect's signed zeros")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    for (const double zero : {0.0, -0.0})
    {
      const Line level   = {.Start = {.X = zero, .Y = 2}, .End = {.X = 10, .Y = 2}};
      const Line upright = {.Start = {.X = zero, .Y = 2}, .End = {.X = zero, .Y = 9}};
      const Line floor   = {.Start = {.X = 1, .Y = zero}, .End = {.X = 6, .Y = zero}};
      const Line wall    = {.Start = {.X = 3, .Y = zero}, .End = {.X = 3, .Y = 5}};
      const Line slant   = {.Start = {.X = -3, .Y = 1}, .End = {.X = 3, .Y = 2}};
      CHECK(SlopeBoth(level, upright, tolerance).Same);
      CHECK(SlopeBoth(floor, wall, tolerance).Same);
      CHECK(SlopeBoth(upright, slant, tolerance).Same);
      CHECK(
          SlopeBoth(floor, {.Start = {.X = 2, .Y = -1}, .End = {.X = 4, .Y = 3}}, tolerance).Same);
    }
    for (const double sign : {1.0, -1.0})
    {
      const Point onX = {.X = 2.0, .Y = 0.0};
      const Point onY = {.X = 0.0, .Y = 3.0};
      CHECK(SlopeBoth({.Start = {.X = 1.0, .Y = -1.0 * sign}, .End = onX},
                      {.Start = onX, .End = {.X = 3.0, .Y = -1.0 * sign}}, tolerance)
                .Same);
      CHECK(SlopeBoth({.Start = {.X = -1.0 * sign, .Y = 2.0}, .End = onY},
                      {.Start = onY, .End = {.X = 1.0 * sign, .Y = 5.0}}, tolerance)
                .Same);
      CHECK(SlopeBoth({.Start = {.X = -5.0, .Y = -5.0 * sign}, .End = {.X = 5.0, .Y = 5.0 * sign}},
                      {.Start = {.X = -5.0, .Y = 5.0 * sign}, .End = {.X = 5.0, .Y = -5.0 * sign}},
                      tolerance)
                .Same);
    }
  }

  namespace
  {
    /**
     *  Two lines, one of no length, as Dotted tells them either first, against Intersect.
     *
     *  @param  point     The line of no length.
     *  @param  line      The other.
     *  @param  tolerance The tolerance.
     *  @return           Whether every answer Dotted gave is Intersect's, and how many orders it
     *                    told.
     */
    auto DotBoth(const Line &point, const Line &line, const Tolerance &tolerance) -> Sloped
    {
      Sloped met;
      for (const auto &[first, second] : {std::pair{&point, &line}, std::pair{&line, &point}})
      {
        const Curve                   one   = *first;
        const Curve                   other = *second;
        const Prepared                readied{one};
        const Prepared                across{other};
        const std::optional<Abutment> dotted = Abutting::Dotted({readied, across}, Sane(tolerance));
        if (dotted)
        {
          const Intersections meeting  = Intersect(one, other, tolerance);
          met.Told                    += 1;
          met.Same                     = met.Same
                                      && meeting.Overlaps.empty()
                                      && Marked(meeting, std::span{dotted->Marks}.first(dotted->Count));
        }
      }
      return met;
    }
  } // namespace

  TEST_CASE("slanted: a line of no length on, near or far from a line, as Intersect answers")
  {
    /**
     *  A point on a line, at an end, past an end, a quarter, a half, the whole and twice the
     *  tolerance off it, at every scale and slant, either first: as Intersect answers, made ready
     *  and through Plainly; and as Dotted tells it, where it does.
     */
    Draws       draws{0x5EED'2049};
    std::size_t told  = 0;
    const auto  check = [&told](const Line &point, const Line &line,
                                const Tolerance &tolerance) -> void
    {
      const Sloped met   = SlopeBoth(point, line, tolerance);
      const Sloped dots  = DotBoth(point, line, tolerance);
      told              += static_cast<std::size_t>(dots.Told);
      CHECK(met.Same);
      CHECK(dots.Same);
    };
    for (int round = 0; round < Rounds(3000); ++round)
    {
      const double scale = std::ldexp(1.0, static_cast<int>(draws.Next() % 41) - 20);
      const Point  start = {.X = draws.Near(0.0, 2.0) * scale, .Y = draws.Near(0.0, 2.0) * scale};
      const Vector way   = round % 4 == 0 ? Vector{.X = draws.Near(0.0, 1.0) * scale, .Y = 0.0}
                                          : Vector{.X = draws.Near(0.0, 1.0) * scale,
                                                   .Y = draws.Near(0.0, 1.0) * scale};
      const Line   line  = {.Start = start, .End = start + way};
      const double tolerance = 1e-9 * scale;
      const double off       = std::array{0.0, 0.25, 0.5, 1.0, 2.0}.at(draws.Next() % 5);
      const double along =
          std::array{0.0, 1.0, 1.0 + 1e-10, draws.Near(0.5, 0.5)}.at(draws.Next() % 4);
      const double length = Length(way);
      const Vector normal = {.X = -way.Y / length, .Y = way.X / length};
      const Point  spot   = start + (way * along) + (normal * (off * tolerance));
      check({.Start = spot, .End = spot}, line, {.Length = tolerance, .Radians = 1e-9});
    }
    /**
     *  A point exactly the tolerance off, a line too short, both of no length, an angle too wide,
     *  and an arc.
     */
    const Tolerance tolerance = {.Length = 0x1p-20, .Radians = 1e-9};
    const Line      level     = {.Start = {.X = 0.0, .Y = 0.0}, .End = {.X = 1.0, .Y = 0.0}};
    check({.Start = {.X = 0.5, .Y = 0x1p-20}, .End = {.X = 0.5, .Y = 0x1p-20}}, level, tolerance);
    check({.Start = {.X = 0.5, .Y = 0x1p-22}, .End = {.X = 0.5, .Y = 0x1p-22}},
          {.Start = {.X = 0.5, .Y = 0.0}, .End = {.X = 0.5 + 0x1p-16, .Y = 0.0}}, tolerance);
    check({.Start = {.X = 0.5, .Y = 0.0}, .End = {.X = 0.5, .Y = 0.0}},
          {.Start = {.X = 0.5, .Y = 0.0}, .End = {.X = 0.5, .Y = 0.0}}, tolerance);
    check({.Start = {.X = 0.5, .Y = 0.0}, .End = {.X = 0.5, .Y = 0.0}}, level,
          {.Length = 0x1p-20, .Radians = 0.1});
    const Curve    dot = Line{.Start = {.X = 0.5, .Y = 0.0}, .End = {.X = 0.5, .Y = 0.0}};
    const Curve    arc = Arc{.Center = {.X = 0.5, .Y = 1.0}, .Radius = 1.0, .Sweep = Degrees(90.0)};
    const Prepared point{dot};
    const Prepared round{arc};
    CHECK_FALSE(Abutting::Dotted({point, round}, tolerance).has_value());
    CHECK_FALSE(Abutting::Dotted({round, point}, tolerance).has_value());
    CHECK(told > 0);
  }

  namespace
  {
    /**
     *  Whether two lines told joined meet in the joint alone, made ready and as curves, and cut
     *  idly there: every mark and hit at the joint, no stretch, nothing spent, and Nearest placing
     *  the joint at each line's own end. Where they are not told, nothing is asked but that
     *  Seamless does not tell them either.
     *
     *  @param  first  One line.
     *  @param  second The other.
     *  @param  sane   The tolerance, as Sane makes it.
     *  @return        Whether they were told, and whether what was told holds.
     */
    auto Joining(const Line &first, const Line &second, const Tolerance &sane)
        -> std::pair<bool, bool>
    {
      const bool jointed  = Abutting::Jointed({&first, &second}, sane);
      const bool seamless = Seamless(Abutting::Limbed(first), Abutting::Limbed(second), sane);
      if (!jointed)
      {
        return {false, !seamless};
      }
      const bool     ahead = first.End == second.Start;
      const Point    joint = ahead ? first.End : first.Start;
      const Curve    one   = first;
      const Curve    other = second;
      const Prepared readied{one};
      const Prepared across{other};
      Marking        marking;
      Intersections  meeting;
      Intersect({readied, across}, sane, marking);
      Intersect(Curves{one, other}, sane, meeting);
      const auto joining = [joint](Point point) -> bool
      {
        return point == joint;
      };
      const bool held = marking.Overlaps.empty()
                     && meeting.Overlaps.empty()
                     && !marking.Spent
                     && !meeting.Spent
                     && std::ranges::all_of(marking.Marks, [&joining](const Mark &mark) -> bool
                                            { return joining(mark.At); })
                     && std::ranges::all_of(meeting.Hits, [&joining](const Hit &hit) -> bool
                                            { return joining(hit.At); })
                     && Exactly(Nearest(first, joint).Along, ahead ? 1.0 : 0.0)
                     && Exactly(Nearest(second, joint).Along, ahead ? 0.0 : 1.0);
      return {true, held};
    }
  } // namespace

  TEST_CASE("slanted: lines joined end to end told from the lines alone, as Intersect answers")
  {
    Draws draws{0x51A7ED7};
    int   told  = 0;
    int   pairs = 0;
    for (int round = 0; round < Rounds(8000); ++round)
    {
      const double scale = std::ldexp(1.0, static_cast<int>(draws.Next() % 61) - 30);
      const double grid  = draws.Next() % 3 == 0 ? std::ldexp(scale, -12) : 0.0;
      const Point  vertex =
          Toward({.X = draws.Near(0.0, 8.0) * scale, .Y = draws.Near(0.0, 8.0) * scale}, {}, grid);
      const double bearing = draws.Next() % 4 == 0
                               ? Pi / 2.0 * static_cast<double>(draws.Next() % 4)
                               : draws.Near(0.0, Pi);
      const double turn    = [&draws]() -> double
      {
        if (draws.Next() % 2 == 0)
        {
          return draws.Near(0.0, Pi);
        }
        const double fraction = draws.Near(0.0, 1.0);
        return std::ldexp(fraction, -static_cast<int>(draws.Next() % 40));
      }();
      const double far = draws.Next() % 8 == 0 ? std::ldexp(1.0, -40) : 1.0;
      const Point from = Toward(vertex, Heading(bearing + Pi, draws.Near(2.5, 2.49) * scale), grid);
      const Point onward =
          Toward(vertex, Heading(bearing + turn, draws.Near(2.5, 2.49) * scale * far), grid);
      const Tolerance tolerance = {.Length  = Lengths.at(draws.Next() % Lengths.size()) * scale,
                                   .Radians = 1e-9};
      const Tolerance sane      = Sane(tolerance);
      const Line      incoming  = {.Start = from, .End = vertex};
      const Line      outgoing  = {.Start = vertex, .End = onward};
      for (const auto &[first, second] :
           {std::pair{incoming, outgoing}, std::pair{outgoing, incoming},
            std::pair{Line{.Start = vertex, .End = from}, outgoing}})
      {
        const auto [joined, held] = Joining(first, second, sane);
        CHECK(held);
        told  += static_cast<int>(joined);
        pairs += 1;
      }
    }
    CHECK(told > pairs / 2);
  }

  TEST_CASE("slanted: joints told or left, by hand")
  {
    const Tolerance sane    = Sane({.Length = 1e-9, .Radians = 1e-9});
    const Line      level   = {.Start = {.X = 1.0, .Y = 2.0}, .End = {.X = 3.0, .Y = 2.0}};
    const Line      upright = {.Start = {.X = 3.0, .Y = 2.0}, .End = {.X = 3.0, .Y = 5.0}};
    const Line      slant   = {.Start = {.X = 3.0, .Y = 2.0}, .End = {.X = 4.0, .Y = 2.5}};
    const Line      steep   = {.Start = {.X = 3.0, .Y = 2.0}, .End = {.X = 2.0, .Y = 3.0}};
    const Line      back    = {.Start = {.X = 3.0, .Y = 2.0}, .End = {.X = 1.5, .Y = 2.0 + 1e-12}};
    const Line      dot     = {.Start = {.X = 3.0, .Y = 2.0}, .End = {.X = 3.0, .Y = 2.0}};
    const Line      stub =
        Line{.Start = {.X = 3.0, .Y = 2.0}, .End = {.X = 3.0 - 1e-9, .Y = 2.0 + 2e-9}};
    const Line apart = {.Start = {.X = 7.0, .Y = 2.0}, .End = {.X = 8.0, .Y = 2.0}};
    /**
     *  Lines along each other folded back, the far end apart from the other line by much, or by
     *  twice the tolerance, inside its box: Remote asks Away.
     */
    const Line rising = {.Start = {.X = 0.0, .Y = 0.0}, .End = {.X = 1.0, .Y = 1e-3}};
    const Line wide   = {.Start = {.X = 1.0, .Y = 1e-3}, .End = {.X = 0.5, .Y = 4e-4}};
    const Line narrow = {.Start = {.X = 1.0, .Y = 1e-3}, .End = {.X = 0.5, .Y = 5e-4 + 2e-9}};
    /**
     *  A line shorter than the tolerance's half along each axis but not across, a fold of two such,
     *  and ways too short to have a product; joints whose quick crossing does not settle, told by
     *  the end the lines share, where the exact crossing is: three steps of the least double off an
     *  axis, whose halfway point is itself, and 2^-38 off one beside numbers of 2^17.
     */
    const Line shortly =
        Line{.Start = {.X = 1.0 - 4e-10, .Y = 1.0 - 4e-10}, .End = {.X = 1.0, .Y = 1.0}};
    const Line onward = {.Start = {.X = 1.0, .Y = 1.0}, .End = {.X = 3.0, .Y = 3.0 + 1e-6}};
    const Line folded =
        Line{.Start = {.X = 1.0, .Y = 1.0}, .End = {.X = 1.0 - 3e-10, .Y = 1.0 - 3e-10 + 1e-12}};
    const Line ending = {.Start = {.X = -1.0, .Y = 1.0}, .End = {.X = 0x3p-1074, .Y = 1.0}};
    const Line going =
        Line{.Start = {.X = 0x3p-1074, .Y = 1.0}, .End = {.X = 1.0, .Y = 1.0 + 1e-3}};
    const Line  least     = {.Start = {.X = -1e-170, .Y = -1e-170}, .End = {.X = 0.0, .Y = 0.0}};
    const Line  leasts    = {.Start = {.X = 0.0, .Y = 0.0}, .End = {.X = 1e-170, .Y = 2e-170}};
    const Line  endless   = Line{.Start = {.X = 3.0, .Y = 2.0},
                                 .End   = {.X = std::numeric_limits<double>::infinity(), .Y = 2.5}};
    const Line  returning = {.Start = level.End, .End = level.Start};
    const Line  climbing  = {.Start = {.X = 1.0, .Y = -1.0}, .End = {.X = 1.0, .Y = 0x3p-1074}};
    const Line  onwards   = {.Start = {.X = 1.0, .Y = 0x3p-1074}, .End = {.X = 1.5, .Y = 1.0}};
    const Point sliver    = {.X = 0x1.517f718f5ec06p+16, .Y = 0x1.cff2e210c47bep-39};
    const Line  dropping  = {.Start = {.X = 0x1.08e9ebb34f082p+17, .Y = 0x1.4268d0b61bc2p+15},
                             .End   = sliver};
    const Line  leaving   = {.Start = sliver,
                             .End   = {.X = 0x1.2e1b87aedf062p+16, .Y = -0x1.0abe8bb04a16ep+15}};
    for (const auto &[first, second] :
         {std::pair{level, upright}, std::pair{level, slant}, std::pair{level, steep},
          std::pair{rising, wide}, std::pair{ending, going}, std::pair{climbing, onwards},
          std::pair{dropping, leaving}, std::pair{leaving, dropping}})
    {
      const auto [joined, held] = Joining(first, second, sane);
      CHECK(joined);
      CHECK(held);
    }
    for (const auto &[first, second] :
         {std::pair{level, back}, std::pair{level, dot}, std::pair{level, stub},
          std::pair{level, apart}, std::pair{rising, narrow}, std::pair{shortly, onward},
          std::pair{shortly, folded}, std::pair{least, leasts}, std::pair{level, endless},
          std::pair{level, returning}})
    {
      const auto [joined, held] = Joining(first, second, sane);
      CHECK_FALSE(joined);
      CHECK(held);
    }
    CHECK_FALSE(Abutting::Jointed({&level, &level}, sane));
    CHECK(Joining(level, slant, {.Length = 1e-9, .Radians = 0.1}) == std::pair{false, true});
    const Line far     = {.Start = {.X = 1e308, .Y = 1e308}, .End = {.X = 1.0, .Y = 1.0}};
    const Line beyond  = {.Start = {.X = 3.0, .Y = 3.0 + 1e-6}, .End = {.X = 1e308, .Y = 1e308}};
    const Line huge    = {.Start = {.X = 1e300, .Y = 1.0}, .End = {.X = 2e300, .Y = 3.0}};
    const Line huger   = {.Start = {.X = 2e300, .Y = 3.0}, .End = {.X = 3e300, .Y = 2.0}};
    const Line small   = {.Start = {.X = 1e-40, .Y = 0.0}, .End = {.X = 2e-40, .Y = 1e-40}};
    const Line smaller = {.Start = {.X = 2e-40, .Y = 1e-40}, .End = {.X = 3e-40, .Y = 3e-40}};
    CHECK_FALSE(Abutting::Jointed({&far, &onward}, sane));
    CHECK_FALSE(Abutting::Jointed({&onward, &beyond}, sane));
    CHECK_FALSE(Abutting::Jointed({&huge, &huger}, sane));
    CHECK(Joining(small, smaller, Sane({.Length = 1e-60, .Radians = 1e-9}))
          == std::pair{false, true});
  }
} // namespace Geoxy
