/***************************************************************************************************
 *  @file      tests/geoxy/direct.cpp
 *  @brief     Tests that Intersect's direct ways answer as the thorough way does, to the bit, on
 *             random pairs near every edge the direct ways leave to it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "draws.hpp"
#include "thorough.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Splitmix64: numbers that look random, the same on every run and every machine.
     */
    struct Mixer
    {
      std::uint64_t State = 0x2545F4914F6CDD1D; /**< Where the sequence is. */

      /**
       *  The next number.
       *
       *  @return The number.
       */
      auto Next() -> std::uint64_t
      {
        State             += 0x9E3779B97F4A7C15;
        std::uint64_t mix  = State;
        mix                = (mix ^ (mix >> 30U)) * 0xBF58476D1CE4E5B9;
        mix                = (mix ^ (mix >> 27U)) * 0x94D049BB133111EB;
        return mix ^ (mix >> 31U);
      }

      /**
       *  A double from -1 to 1.
       *
       *  @return The double.
       */
      auto Unit() -> double { return (static_cast<double>(Next() >> 11U) * 0x1p-52) - 1.0; }

      /**
       *  A choice among some.
       *
       *  @param  count How many.
       *  @return       Which, from 0.
       */
      auto Pick(std::uint64_t count) -> std::uint64_t { return Next() % count; }
    };

    /**
     *  Whether two meetings are the same to the bit, hit for hit and stretch for stretch.
     *
     *  @param  one   One.
     *  @param  other The other.
     *  @return       Whether they are.
     */
    auto Identical(const Intersections &one, const Intersections &other) -> bool
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
     *  A second line made to lie near an edge of what the direct way decides: crossing the first
     *  plainly, crossing it at a shallow angle, an end near it, parallel near it, along it.
     *
     *  @param  mixer Where the numbers come from.
     *  @param  first The first line.
     *  @param  near  The tolerance's length.
     *  @return       The second line.
     */
    auto Second(Mixer &mixer, const Line &first, double near) -> Line
    {
      const Vector way    = first.End - first.Start;
      const Vector across = Vector{.X = -way.Y, .Y = way.X} / Length(way);
      const Point  middle = first.Start + (way * (0.5 + (0.4 * mixer.Unit())));
      const double scale  = Length(way);
      switch (mixer.Pick(5))
      {
      case 0:
        return {.Start = middle + (across * scale * mixer.Unit()) + (way * mixer.Unit()),
                .End   = middle + (across * scale * mixer.Unit()) + (way * mixer.Unit())};
      case 1:
      {
        const double slope = std::ldexp(mixer.Unit(), -static_cast<int>(mixer.Pick(40)));
        return {.Start = middle - (way * 0.3) - (across * slope * scale),
                .End   = middle + (way * 0.3) + (across * slope * scale)};
      }
      case 2:
      {
        const double off = near * (1.0 + (6.0 * mixer.Unit()));
        return {.Start = middle + (across * off),
                .End   = middle + (across * scale * (1.0 + mixer.Unit()))};
      }
      case 3:
      {
        const double off = near * 8.0 * mixer.Unit();
        return {.Start = middle + (across * off), .End = middle + (across * off) + (way * 0.3)};
      }
      default:
        return {.Start = middle, .End = middle + (way * (0.2 + mixer.Unit()))};
      }
    }

    /**
     *  A line as itself, half the time, else as the ray from its start through its end or the
     *  endless line through both.
     *
     *  @param  mixer Where the choice comes from.
     *  @param  line  The line.
     *  @return       The curve.
     */
    auto Kinded(Mixer &mixer, const Line &line) -> Curve
    {
      switch (mixer.Pick(4))
      {
      case 0:
        return Ray{.Start = line.Start, .Way = line.End - line.Start};
      case 1:
        return Xline{.Through = line.Start, .Way = line.End - line.Start};
      default:
        return line;
      }
    }

    /**
     *  Where an arc goes, roughly, and the tolerance it is met within.
     */
    struct Spot
    {
      Point  Centre = {};  /**< Roughly where its center goes. */
      double Size   = 0.0; /**< Roughly its radius.            */
      double Near   = 0.0; /**< The tolerance's length.        */
    };

    /**
     *  An arc near a line or another arc at an edge of what the direct ways decide: crossing it
     *  plainly, near tangent to it within a few tolerances either way, an end near it, or a circle
     *  near the same; at times a whole circle.
     *
     *  @param  mixer Where the numbers come from.
     *  @param  spot  Where it goes and the tolerance.
     *  @return       The arc.
     */
    auto Round(Mixer &mixer, const Spot &spot) -> Arc
    {
      const double radius = spot.Size * (0.5 + std::abs(mixer.Unit()));
      const double offset = mixer.Pick(3) == 0 ? spot.Near * 6.0 * mixer.Unit() : 0.0;
      const Angle  start  = Radians(Pi * mixer.Unit());
      const Angle  sweep  = mixer.Pick(4) == 0 ? Turn : Radians(2.0 * Pi * mixer.Unit());
      return {.Center = spot.Centre, .Radius = radius + offset, .Start = start, .Sweep = sweep};
    }

    /**
     *  What the thorough way finds, as a meeting of its own.
     *
     *  @param  one       One curve.
     *  @param  other     The other.
     *  @param  tolerance The tolerance.
     *  @return           The meeting.
     */
    auto Thoroughly(const Curve &one, const Curve &other, const Tolerance &tolerance)
        -> Intersections
    {
      Intersections meeting;
      Thorough({one, other}, tolerance, meeting);
      return meeting;
    }

    /**
     *  Whether two meetings are one within a reach: as many hits, each of the kind of the other's
     *  and lying within the reach of it, and the same stretches; a touch's two points within twice
     *  the tolerance of each other, as two hits that near are one.
     *
     *  @param  pair   The two meetings, the one met in closed form first.
     *  @param  reach  How far apart two hits may lie.
     *  @param  length The tolerance's length.
     *  @return        Whether they are.
     */
    auto Close(std::pair<const Intersections &, const Intersections &> pair, double reach,
               double length) -> bool
    {
      const auto &[met, careful] = pair;
      const auto hit             = [reach, length](const Hit &lhs, const Hit &rhs) -> bool
      {
        const bool held =
            lhs.Kind != Contact::Touch || Distance(lhs.First.At, lhs.Second.At) <= 2.0 * length;
        /**
         *  The thorough way samples a band, and so may miss where the curves turn parallel by more
         *  than the tolerance's angle and call its touch a crossing; the closed form finds that
         *  place.
         */
        const bool kind =
            lhs.Kind == rhs.Kind || (lhs.Kind == Contact::Touch && rhs.Kind == Contact::Cross);
        return kind && Distance(lhs.At, rhs.At) <= reach && held;
      };
      /**
       *  Stretches of arcs on one circle met in closed form: the same, to a hair of their places.
       */
      const auto overlap = [](const Overlap &lhs, const Overlap &rhs) -> bool
      {
        constexpr double Hair = 1e-9;
        const auto       near = [](Range one, Range other) -> bool
        {
          return std::abs(one.Start - other.Start) <= Hair && std::abs(one.End - other.End) <= Hair;
        };
        return near(lhs.First, rhs.First) && near(lhs.Second, rhs.Second) && lhs.Along == rhs.Along;
      };
      /**
       *  The thorough way may keep, where arcs of one circle meet end to end, a stretch a rounding
       *  long; the closed form makes it the hit at an end it is.
       */
      const auto sliver = [](const Overlap &one) -> bool
      {
        constexpr double Hair = 1e-9;
        return one.First.End - one.First.Start <= Hair && one.Second.End - one.Second.Start <= Hair;
      };
      std::vector<Overlap> stretches;
      std::vector<double>  slivers;
      for (const Overlap &one : careful.Overlaps)
      {
        if (sliver(one))
        {
          slivers.push_back(one.First.Start);
        }
        else
        {
          stretches.push_back(one);
        }
      }
      std::vector<Hit> hits;
      for (const Hit &one : met.Hits)
      {
        const bool meets =
            one.Kind == Contact::End
            && std::ranges::any_of(slivers, [&one](double along) -> bool
                                   { return std::abs(one.First.Along - along) <= 1e-9; });
        if (!meets)
        {
          hits.push_back(one);
        }
      }
      std::vector<Overlap> kept;
      std::ranges::copy_if(met.Overlaps, std::back_inserter(kept),
                           [&sliver](const Overlap &one) -> bool { return !sliver(one); });
      return std::ranges::equal(hits, careful.Hits, hit)
          && std::ranges::equal(kept, stretches, overlap)
          && met.Spent == careful.Spent;
    }

    /**
     *  Whether a meeting of two straight curves finds the crossing the thorough way loses: the
     *  thorough way works a crossing out in doubles first, and where the two part by no more than
     *  twice the tolerance's angle its two places may land too far apart to be within the
     *  tolerance, so that it finds nothing. One hit where it found nothing, a touch or a crossing,
     *  its places on the curves and within the tolerance of each other, is the crossing it lost;
     *  its point is checked against rationals where two such lines are met on their own.
     *
     *  @param  curves    The two curves.
     *  @param  tolerance The tolerance.
     *  @param  pair      The meeting, and the thorough way's.
     *  @return           Whether it does.
     */
    auto Regained(const Curves &curves, const Tolerance &tolerance,
                  std::pair<const Intersections &, const Intersections &> pair) -> bool
    {
      const auto &[met, careful] = pair;
      const auto way             = [](const Curve &curve) -> std::optional<Vector>
      {
        const bool straight = std::holds_alternative<Line>(curve)
                           || std::holds_alternative<Ray>(curve)
                           || std::holds_alternative<Xline>(curve);
        return straight ? std::optional<Vector>{Derivative(curve, 0.0)} : std::nullopt;
      };
      const std::optional<Vector> one   = way(curves[0]);
      const std::optional<Vector> other = way(curves[1]);
      if (!one
          || !other
          || !careful.Hits.empty()
          || !careful.Overlaps.empty()
          || met.Hits.size() != 1
          || !met.Overlaps.empty())
      {
        return false;
      }
      const double sine = std::abs(Cross(*one, *other)) / (Length(*one) * Length(*other));
      const Hit   &hit  = met.Hits[0];
      return sine                                  <= 2.0 * std::sin(tolerance.Radians)
          && hit.Kind                              != Contact::End
          && hit.First.At                          == At(curves[0], hit.First.Along)
          && hit.Second.At                         == At(curves[1], hit.Second.Along)
          && Distance(hit.First.At, hit.Second.At) <= tolerance.Length;
    }

    /**
     *  Whether two curves are met as the thorough way meets them: to the bit; or, where they are a
     *  band met in closed form, within the band's reach, as the thorough way finds its touch
     *  somewhere along it: four times its half length, the band of a circle of the radius they part
     *  at; or, where two straight curves cross so shallowly that the thorough way loses the
     *  crossing, at it.
     *
     *  @param  curves    The two curves.
     *  @param  tolerance The tolerance.
     *  @return           Whether they are.
     */
    auto Agrees(const Curves &curves, const Tolerance &tolerance) -> bool
    {
      const Intersections met     = Intersect(curves[0], curves[1], tolerance);
      const Intersections careful = Thoroughly(curves[0], curves[1], tolerance);
      if (Identical(met, careful) || Regained(curves, tolerance, {met, careful}))
      {
        return true;
      }
      const Tolerance usable = Usable(curves, tolerance);
      Intersections   band;
      if (!Banded(curves, usable, band) || !Identical(band, met))
      {
        return false;
      }
      /**
       *  How sharply the two part where they come nearest, as a radius: a line's from an arc is the
       *  arc's, two circles' the radius of the difference of their curvatures, from inside, or of
       *  their sum, from outside.
       */
      std::array<const Arc *, 2> arcs  = {};
      std::size_t                count = 0;
      for (const Curve &curve : curves)
      {
        if (const auto *arc = std::get_if<Arc>(&curve))
        {
          arcs.at(count) = arc;
          ++count;
        }
      }
      double radius = count == 1 ? arcs[0]->Radius : 0.0;
      if (count == 2)
      {
        const double one   = arcs[0]->Radius;
        const double other = arcs[1]->Radius;
        const double apart = Distance(arcs[0]->Center, arcs[1]->Center);
        const bool   outside =
            std::abs(apart - (one + other)) < std::abs(apart - std::abs(one - other));
        /**
         *  Arcs of one circle part nowhere: their ends, met in closed form, are where they lie.
         */
        const double parting = outside ? one + other : std::abs(one - other);
        radius               = parting > usable.Length ? one * other / parting : 0.0;
      }
      const double reach = (4.0 * std::sqrt(2.0 * radius * usable.Length)) + (4.0 * usable.Length);
      return Close({met, careful}, reach, usable.Length);
    }

    /**
     *  Whether a direct way decides two curves' meeting.
     *
     *  @param  one       One curve.
     *  @param  other     The other.
     *  @param  tolerance The tolerance.
     *  @return           Whether it does.
     */
    auto Decided(const Curve &one, const Curve &other, const Tolerance &tolerance) -> bool
    {
      Intersections meeting;
      return Direct({one, other}, tolerance, meeting);
    }
  } // namespace

  TEST_CASE("direct: lines, rays and endless lines met directly are as met thoroughly")
  {
    Mixer     mixer;
    int       differ = 0;
    int       direct = 0;
    const int rounds = Rounds(50000);
    for (int round = 0; round < rounds; ++round)
    {
      const double away  = std::ldexp(mixer.Unit(), static_cast<int>(mixer.Pick(28)));
      const double size  = std::ldexp(1.0, static_cast<int>(mixer.Pick(12)) - 4);
      const Point  base  = {.X = away, .Y = -away};
      const Line   first = {.Start = base + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size),
                            .End   = base + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size)};
      const double near  = std::ldexp(size, -static_cast<int>(mixer.Pick(30)) - 6);
      const Line   other = Second(mixer, first, near);
      const Tolerance tolerance = {.Length = near, .Radians = 1e-9};
      const Curve     lined     = Kinded(mixer, first);
      const Curve     crossed   = Kinded(mixer, other);
      for (const bool swapped : {false, true})
      {
        const Curve &one = swapped ? crossed : lined;
        const Curve &two = swapped ? lined : crossed;
        if (!Agrees({one, two}, tolerance))
        {
          ++differ;
          if (differ <= 5)
          {
            MESSAGE("differ: ("
                    << first.Start.X
                    << ", "
                    << first.Start.Y
                    << ")-("
                    << first.End.X
                    << ", "
                    << first.End.Y
                    << ") as "
                    << lined.index()
                    << " and ("
                    << other.Start.X
                    << ", "
                    << other.Start.Y
                    << ")-("
                    << other.End.X
                    << ", "
                    << other.End.Y
                    << ") as "
                    << crossed.index()
                    << " within "
                    << near);
          }
        }
      }
      direct += Decided(lined, crossed, tolerance) ? 1 : 0;
    }
    CHECK(differ == 0);
    MESSAGE("decided directly: " << direct << " of " << rounds);
    CHECK(direct > rounds / 4);
  }

  TEST_CASE("direct: an end on a line, and two lines meeting at their ends, as met thoroughly")
  {
    Mixer     mixer;
    int       differ = 0;
    int       direct = 0;
    const int rounds = Rounds(40000);
    for (int round = 0; round < rounds; ++round)
    {
      const double away  = std::ldexp(mixer.Unit(), static_cast<int>(mixer.Pick(28)));
      const double size  = std::ldexp(1.0, static_cast<int>(mixer.Pick(12)) - 4);
      const Point  base  = {.X = away, .Y = -away};
      const Line   first = {.Start = base + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size),
                            .End   = base + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size)};
      const Vector way   = Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size;
      const double place = mixer.Pick(3) == 0 ? static_cast<double>(mixer.Pick(2)) : mixer.Unit();
      const Point  foot  = Lerp(first.Start, first.End, place);
      const Line   other = mixer.Pick(2) == 0 ? Line{.Start = foot, .End = foot + way}
                                              : Line{.Start = foot - way, .End = foot};
      const double near  = std::ldexp(size, -static_cast<int>(mixer.Pick(30)) - 6);
      const Tolerance tolerance = {.Length = near, .Radians = 1e-9};
      const Curve     lined     = Kinded(mixer, first);
      const Curve     crossed   = Kinded(mixer, other);
      for (const bool swapped : {false, true})
      {
        const Curve &one = swapped ? crossed : lined;
        const Curve &two = swapped ? lined : crossed;
        if (!Agrees({one, two}, tolerance))
        {
          ++differ;
          if (differ <= 5)
          {
            MESSAGE("differ: round " << round << " swapped " << swapped);
          }
        }
      }
      direct += Decided(lined, crossed, tolerance) ? 1 : 0;
    }
    CHECK(differ == 0);
    MESSAGE("decided directly: " << direct << " of " << rounds);
    CHECK(direct > rounds / 10);
  }

  TEST_CASE("direct: lines and arcs, arcs and arcs met directly are as met thoroughly")
  {
    Mixer     mixer;
    int       differ = 0;
    int       direct = 0;
    const int rounds = Rounds(30000);
    for (int round = 0; round < rounds; ++round)
    {
      const double away   = std::ldexp(mixer.Unit(), static_cast<int>(mixer.Pick(24)));
      const double size   = std::ldexp(1.0, static_cast<int>(mixer.Pick(10)) - 3);
      const double near   = std::ldexp(size, -static_cast<int>(mixer.Pick(30)) - 6);
      const Point  centre = {.X = away + (size * mixer.Unit()), .Y = -away + (size * mixer.Unit())};
      const Arc    arc    = Round(mixer, {.Centre = centre, .Size = size, .Near = near});
      Curve        other;
      switch (mixer.Pick(3))
      {
      case 0:
      {
        const double reach = arc.Radius + (near * 8.0 * mixer.Unit());
        const Vector way   = {.X = mixer.Unit(), .Y = mixer.Unit()};
        const Vector side  = Perpendicular(way) / Length(way);
        const Point  touch = centre + (side * reach);
        other = Line{.Start = touch - (way * size * 3.0), .End = touch + (way * size * 3.0)};
        break;
      }
      case 1:
        other = Line{.Start = centre + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size * 2.0),
                     .End   = centre + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size * 2.0)};
        break;
      default:
      {
        const Point beside = centre + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size);
        other = Round(mixer,
                      {.Centre = mixer.Pick(4) == 0 ? centre : beside, .Size = size, .Near = near});
        break;
      }
      }
      const Tolerance tolerance = {.Length = near, .Radians = 1e-9};
      for (const bool swapped : {false, true})
      {
        const Curve         one     = swapped ? other : Curve{arc};
        const Curve         two     = swapped ? Curve{arc} : other;
        const Intersections met     = Intersect(one, two, tolerance);
        const Intersections careful = Thoroughly(one, two, tolerance);
        if (!Agrees({one, two}, tolerance))
        {
          ++differ;
          if (differ <= 5)
          {
            MESSAGE("differ at round "
                    << round
                    << (swapped ? ", swapped" : "")
                    << ": "
                    << met.Hits.size()
                    << " hits against "
                    << careful.Hits.size());
            const auto show = [](const Curve &curve) -> std::string
            {
              if (const auto *line = std::get_if<Line>(&curve))
              {
                return "line "
                     + std::to_string(line->Start.X)
                     + ","
                     + std::to_string(line->Start.Y)
                     + " "
                     + std::to_string(line->End.X)
                     + ","
                     + std::to_string(line->End.Y);
              }
              const Arc &circular = std::get<Arc>(curve);
              return "arc c "
                   + std::to_string(circular.Center.X)
                   + ","
                   + std::to_string(circular.Center.Y)
                   + " r "
                   + std::to_string(circular.Radius)
                   + " start "
                   + std::to_string(circular.Start.InRadians())
                   + " sweep "
                   + std::to_string(circular.Sweep.InRadians());
            };
            MESSAGE(show(one) << " | " << show(two) << " | L " << near);
            for (const Hit &hit : met.Hits)
            {
              MESSAGE("  direct "
                      << hit.At.X
                      << ","
                      << hit.At.Y
                      << " kind "
                      << static_cast<int>(hit.Kind)
                      << " along "
                      << hit.First.Along
                      << " "
                      << hit.Second.Along);
            }
            for (const Hit &hit : careful.Hits)
            {
              MESSAGE("  thorough "
                      << hit.At.X
                      << ","
                      << hit.At.Y
                      << " kind "
                      << static_cast<int>(hit.Kind)
                      << " along "
                      << hit.First.Along
                      << " "
                      << hit.Second.Along);
            }
          }
        }
        direct += Decided(one, two, tolerance) ? 1 : 0;
      }
    }
    MESSAGE("decided directly: " << direct << " of " << (2 * rounds));
    CHECK(differ == 0);
    CHECK(direct > rounds / 2);
  }

  TEST_CASE("direct: ends on circles, lines through arcs' ends and arcs through arcs' ends")
  {
    Mixer     mixer;
    int       differ = 0;
    int       direct = 0;
    const int rounds = Rounds(40000);
    for (int round = 0; round < rounds; ++round)
    {
      const double away   = std::ldexp(mixer.Unit(), static_cast<int>(mixer.Pick(24)));
      const double size   = std::ldexp(1.0, static_cast<int>(mixer.Pick(12)) - 4);
      const Point  centre = {.X = away, .Y = -away};
      const bool   whole  = mixer.Pick(4) == 0;
      const Arc   arc = {.Center = centre,
                         .Radius = size,
                         .Start  = Degrees(std::round(mixer.Unit() * 360.0)),
                         .Sweep = whole ? Turn : Degrees(std::round(mixer.Unit() * 600.0) - 300.0)};
      const Point start  = At(arc, 0.0);
      const Point end    = At(arc, 1.0);
      const Point inside = At(arc, mixer.Unit());
      const Vector way   = Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size * 3.0;
      Curve        other = arc;
      switch (mixer.Pick(4))
      {
      case 0:
        other = Kinded(mixer, Line{.Start = start, .End = start + way});
        break;
      case 1:
        other = Kinded(mixer, Line{.Start = inside, .End = inside + way});
        break;
      case 2:
        other = Kinded(mixer, Line{.Start = start, .End = end});
        break;
      default:
      {
        const Point middle = centre + (way * 0.5);
        other              = Arc{.Center = middle,
                                 .Radius = Distance(middle, mixer.Pick(2) == 0 ? start : inside),
                                 .Start  = Degrees(std::round(mixer.Unit() * 360.0)),
                                 .Sweep  = Degrees(std::round(mixer.Unit() * 720.0) - 360.0)};
        break;
      }
      }
      const double    near      = std::ldexp(size, -static_cast<int>(mixer.Pick(30)) - 8);
      const Tolerance tolerance = {.Length = near, .Radians = 1e-9};
      const Curve     first     = arc;
      for (const bool swapped : {false, true})
      {
        const Curve &one = swapped ? other : first;
        const Curve &two = swapped ? first : other;
        if (!Agrees({one, two}, tolerance))
        {
          ++differ;
          if (differ <= 5)
          {
            MESSAGE("differ: round " << round << " swapped " << swapped);
          }
        }
      }
      direct += Decided(first, other, tolerance) ? 1 : 0;
    }
    CHECK(differ == 0);
    MESSAGE("decided directly: " << direct << " of " << rounds);
    CHECK(direct > rounds / 4);
  }

  TEST_CASE("direct: lines and circles touching exactly, as met thoroughly")
  {
    Mixer     mixer;
    int       differ = 0;
    int       direct = 0;
    const int rounds = Rounds(40000);
    for (int round = 0; round < rounds; ++round)
    {
      /**
       *  Half the arcs start on a quarter turn, where the touches lie, so that a seam is often at
       *  the touch.
       */
      const auto started = [&mixer]() -> Angle
      {
        return Degrees(mixer.Pick(2) == 0 ? 90.0 * static_cast<double>(mixer.Pick(4))
                                          : std::round(mixer.Unit() * 360.0));
      };
      const double unit   = std::ldexp(1.0, static_cast<int>(mixer.Pick(12)) - 6);
      const double radius = unit * static_cast<double>(1 + mixer.Pick(64));
      const Point  centre = {.X = unit * static_cast<double>(mixer.Pick(1U << 20U)),
                             .Y = -unit * static_cast<double>(mixer.Pick(1U << 20U))};
      const Arc    arc    = {
          .Center = centre,
          .Radius = radius,
          .Start  = started(),
          .Sweep  = mixer.Pick(3) == 0 ? Turn : Degrees(std::round(mixer.Unit() * 600.0) - 300.0)};
      const double side  = mixer.Pick(2) == 0 ? 1.0 : -1.0;
      const double reach = unit * static_cast<double>(1 + mixer.Pick(256));
      const double from  = unit * (static_cast<double>(mixer.Pick(256)) - 128.0);
      Curve        other = arc;
      switch (mixer.Pick(3))
      {
      case 0:
      {
        const double level = centre.Y + (side * radius);
        other              = Kinded(mixer, Line{.Start = {.X = centre.X + from, .Y = level},
                                                .End   = {.X = centre.X + from + reach, .Y = level}});
        break;
      }
      case 1:
      {
        const double level = centre.X + (side * radius);
        other              = Kinded(mixer, Line{.Start = {.X = level, .Y = centre.Y + from},
                                                .End   = {.X = level, .Y = centre.Y + from - reach}});
        break;
      }
      default:
      {
        const double size  = unit * static_cast<double>(1 + mixer.Pick(64));
        const bool   inner = mixer.Pick(2) == 0 && size < radius;
        const double apart = inner ? radius - size : radius + size;
        other              = Arc{
            .Center = {.X = centre.X + (side * apart), .Y = centre.Y},
            .Radius = size,
            .Start  = started(),
            .Sweep = mixer.Pick(3) == 0 ? Turn : Degrees(std::round(mixer.Unit() * 600.0) - 300.0)};
        break;
      }
      }
      const double    near      = std::ldexp(unit, -static_cast<int>(mixer.Pick(30)) - 8);
      const Tolerance tolerance = {.Length = near, .Radians = 1e-9};
      const Curve     first     = arc;
      for (const bool swapped : {false, true})
      {
        const Curve &one = swapped ? other : first;
        const Curve &two = swapped ? first : other;
        if (!Agrees({one, two}, tolerance))
        {
          ++differ;
          if (differ <= 5)
          {
            MESSAGE("differ: round " << round << " swapped " << swapped);
          }
        }
      }
      direct += Decided(first, other, tolerance) ? 1 : 0;
    }
    CHECK(differ == 0);
    MESSAGE("decided directly: " << direct << " of " << rounds);
    CHECK(direct > rounds / 4);
  }

  TEST_CASE("direct: lines and circles a few tolerances from touching, bands in closed form")
  {
    Mixer     mixer;
    int       differ = 0;
    int       direct = 0;
    const int rounds = Rounds(40000);
    for (int round = 0; round < rounds; ++round)
    {
      /**
       *  Half the arcs start on a quarter turn, where the touches lie, so that a seam is often at
       *  the touch.
       */
      const auto started = [&mixer]() -> Angle
      {
        return Degrees(mixer.Pick(2) == 0 ? 90.0 * static_cast<double>(mixer.Pick(4))
                                          : std::round(mixer.Unit() * 360.0));
      };
      const double unit   = std::ldexp(1.0, static_cast<int>(mixer.Pick(12)) - 6);
      const double radius = unit * static_cast<double>(1 + mixer.Pick(64));
      const Point  centre = {.X = unit * static_cast<double>(mixer.Pick(1U << 20U)),
                             .Y = -unit * static_cast<double>(mixer.Pick(1U << 20U))};
      const Arc    arc    = {
          .Center = centre,
          .Radius = radius,
          .Start  = started(),
          .Sweep  = mixer.Pick(3) == 0 ? Turn : Degrees(std::round(mixer.Unit() * 600.0) - 300.0)};
      const double side  = mixer.Pick(2) == 0 ? 1.0 : -1.0;
      const double near  = std::ldexp(unit, -static_cast<int>(mixer.Pick(30)) - 8);
      const double shift = near * ((mixer.Unit() * 9.0) - 4.5);
      const double reach = unit * static_cast<double>(1 + mixer.Pick(256));
      const double from  = unit * (static_cast<double>(mixer.Pick(256)) - 128.0);
      Curve        other = arc;
      switch (mixer.Pick(3))
      {
      case 0:
      {
        const double level = centre.Y + (side * (radius + shift));
        other              = Kinded(mixer, Line{.Start = {.X = centre.X + from, .Y = level},
                                                .End   = {.X = centre.X + from + reach, .Y = level}});
        break;
      }
      case 1:
      {
        const double level = centre.X + (side * (radius + shift));
        other              = Kinded(mixer, Line{.Start = {.X = level, .Y = centre.Y + from},
                                                .End   = {.X = level, .Y = centre.Y + from - reach}});
        break;
      }
      default:
      {
        const double size  = unit * static_cast<double>(1 + mixer.Pick(64));
        const bool   inner = mixer.Pick(2) == 0 && size < radius;
        const double apart = (inner ? radius - size : radius + size) + shift;
        other              = Arc{
            .Center = {.X = centre.X + (side * apart), .Y = centre.Y},
            .Radius = size,
            .Start  = started(),
            .Sweep = mixer.Pick(3) == 0 ? Turn : Degrees(std::round(mixer.Unit() * 600.0) - 300.0)};
        break;
      }
      }
      const Tolerance tolerance = {.Length = near, .Radians = 1e-9};
      const Curve     first     = arc;
      for (const bool swapped : {false, true})
      {
        const Curve &one = swapped ? other : first;
        const Curve &two = swapped ? first : other;
        if (!Agrees({one, two}, tolerance))
        {
          ++differ;
          if (differ <= 5)
          {
            MESSAGE("differ: round " << round << " swapped " << swapped);
          }
        }
      }
      direct += Decided(first, other, tolerance) ? 1 : 0;
    }
    CHECK(differ == 0);
    MESSAGE("decided directly: " << direct << " of " << rounds);
    CHECK(direct > rounds / 4);
  }

  TEST_CASE("direct: arcs of one circle, as met thoroughly")
  {
    Mixer     mixer;
    int       differ = 0;
    int       direct = 0;
    const int rounds = Rounds(40000);
    for (int round = 0; round < rounds; ++round)
    {
      const double away   = std::ldexp(mixer.Unit(), static_cast<int>(mixer.Pick(24)));
      const double radius = std::ldexp(1.0 + mixer.Unit(), static_cast<int>(mixer.Pick(12)) - 4);
      const Point  centre = {.X = away, .Y = -away};
      const auto   angle  = [&mixer](double reach) -> Angle
      {
        const double turns = (mixer.Unit() * 2.0 * reach) - reach;
        switch (mixer.Pick(3))
        {
        case 0:
          return Degrees(std::round(turns * 360.0));
        case 1:
          return Degrees(45.0 * std::round(turns * 8.0));
        default:
          return Radians(turns * 2.0 * Pi);
        }
      };
      const Arc first = {
          .Center = centre, .Radius = radius, .Start = angle(1.0), .Sweep = angle(1.0)};
      Arc other = {.Center = centre, .Radius = radius, .Start = angle(1.0), .Sweep = angle(1.0)};
      if (mixer.Pick(4) == 0)
      {
        other.Start = first.Start + first.Sweep;
      }
      if (mixer.Pick(8) == 0)
      {
        other.Sweep = Turn;
      }
      const double    near      = std::ldexp(radius, -static_cast<int>(mixer.Pick(30)) - 8);
      const Tolerance tolerance = {.Length = near, .Radians = 1e-9};
      const Curve     one       = first;
      const Curve     two       = other;
      for (const bool swapped : {false, true})
      {
        const Curve &lhs = swapped ? two : one;
        const Curve &rhs = swapped ? one : two;
        if (!Agrees({lhs, rhs}, tolerance))
        {
          ++differ;
          if (differ <= 5)
          {
            MESSAGE("differ: round " << round << " swapped " << swapped);
          }
        }
      }
      direct += Decided(one, two, tolerance) ? 1 : 0;
    }
    CHECK(differ == 0);
    MESSAGE("decided directly: " << direct << " of " << rounds);
    CHECK(direct > rounds / 10);
  }

  TEST_CASE("direct: arcs of one circle whose ends lie within the tolerance of each other")
  {
    /**
     *  An arc no longer than the tolerance, or all but a whole circle, has both its ends within
     *  half the tolerance of an end of another arc of its circle: as many as four pairs of ends,
     *  more than a closed form holds. Each is met as the thorough way meets it, never by throwing.
     */
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Point     origin    = {.X = 0, .Y = 0};
    const Arc       tiny      = {.Center = {.X = 0x1.2p+3, .Y = 0x1.2ca78e23c4708p+3},
                                 .Radius = 0x1.da16bae200d8bp+5,
                                 .Start  = Radians(0x1.5e6f7f45e7aacp+1),
                                 .Sweep  = Radians(0x1.7f00f350a8ddcp-1003)};
    const Arc       before    = {
        .Center = origin, .Radius = 1, .Start = Radians(-1e-12), .Sweep = Radians(1e-12)};
    const Arc quarter = {.Center = origin, .Radius = 1, .Start = {}, .Sweep = Degrees(90)};
    const Arc nearly  = {.Center = origin,
                         .Radius = 1,
                         .Start  = Radians(1e-12),
                         .Sweep  = Radians((2.0 * Pi) - 2e-12)};
    const Arc gap     = {
        .Center = origin, .Radius = 1, .Start = Radians(-1e-12), .Sweep = Radians(2e-12)};
    const std::array<std::array<Curve, 2>, 3> pairs = {std::array<Curve, 2>{tiny, tiny},
                                                       std::array<Curve, 2>{before, quarter},
                                                       std::array<Curve, 2>{nearly, gap}};
    for (const std::array<Curve, 2> &pair : pairs)
    {
      for (const bool swapped : {false, true})
      {
        const Curve &lhs = pair.at(swapped ? 1 : 0);
        const Curve &rhs = pair.at(swapped ? 0 : 1);
        CHECK(Agrees({lhs, rhs}, tolerance));
        const Prepared one{lhs};
        const Prepared two{rhs};
        Intersections  ready;
        Intersect(Prepared::Pair{std::cref(one), std::cref(two)}, tolerance, ready);
        CHECK(Identical(ready, Intersect(lhs, rhs, tolerance)));
      }
    }
    CHECK(Intersect(tiny, tiny, tolerance).Overlaps.size()  == 1);
    CHECK(Intersect(before, quarter, tolerance).Hits.size() == 1);
  }

  TEST_CASE("direct: lines stopping short of circles, within the tolerance and past it")
  {
    Mixer     mixer;
    int       differ = 0;
    int       direct = 0;
    const int rounds = Rounds(40000);
    for (int round = 0; round < rounds; ++round)
    {
      const double away   = std::ldexp(mixer.Unit(), static_cast<int>(mixer.Pick(24)));
      const double radius = std::ldexp(1.0 + mixer.Unit(), static_cast<int>(mixer.Pick(12)) - 4);
      const Point  centre = {.X = away, .Y = -away};
      const double near   = std::ldexp(radius, -static_cast<int>(mixer.Pick(30)) - 8);
      const double angle  = mixer.Unit() * 2.0 * Pi;
      const Vector toward = {.X = std::cos(angle), .Y = std::sin(angle)};
      const double tilt   = (mixer.Unit() - 0.5) * 2.4;
      const Vector way    = Vector{.X = (toward.X * std::cos(tilt)) - (toward.Y * std::sin(tilt)),
                                   .Y = (toward.X * std::sin(tilt)) + (toward.Y * std::cos(tilt))}
                          * radius;
      const double gap    = near * ((mixer.Unit() * 6.0) - 0.5) * (mixer.Pick(2) == 0 ? 1.0 : -1.0);
      const Point  end    = centre + (toward * (radius + gap));
      const Point  from   = end + (way * (mixer.Pick(2) == 0 ? 1.0 : -1.0) * (0.2 + mixer.Unit()));
      const Line   line =
          mixer.Pick(2) == 0 ? Line{.Start = from, .End = end} : Line{.Start = end, .End = from};
      const Curve other =
          mixer.Pick(4) == 0 ? Curve{Ray{.Start = end, .Way = from - end}} : Curve{line};
      const Arc arc = {
          .Center = centre,
          .Radius = radius,
          .Start  = Degrees(std::round(mixer.Unit() * 360.0)),
          .Sweep  = mixer.Pick(3) == 0 ? Turn : Degrees(std::round(mixer.Unit() * 600.0) - 300.0)};
      const Tolerance tolerance = {.Length = near, .Radians = 1e-9};
      const Curve     first     = arc;
      for (const bool swapped : {false, true})
      {
        const Curve &one = swapped ? other : first;
        const Curve &two = swapped ? first : other;
        if (!Agrees({one, two}, tolerance))
        {
          ++differ;
          if (differ <= 5)
          {
            MESSAGE("differ: round " << round << " swapped " << swapped);
          }
        }
      }
      direct += Decided(first, other, tolerance) ? 1 : 0;
    }
    CHECK(differ == 0);
    MESSAGE("decided directly: " << direct << " of " << rounds);
    CHECK(direct > rounds / 4);
  }

  TEST_CASE("direct: bands off their curves, curves far apart, circles about one center")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const double    hair      = tolerance.Length / 2;
    const Curve     unit      = Arc{.Center = {}, .Radius = 1, .Start = Degrees(90), .Sweep = Turn};
    const Curve     lower =
        Arc{.Center = {}, .Radius = 1, .Start = Degrees(180), .Sweep = Degrees(180)};
    const Curve upper = Arc{.Center = {}, .Radius = 1, .Start = Degrees(30), .Sweep = Degrees(120)};
    const auto  banded = [&tolerance](const Curve &one,
                                      const Curve &other) -> std::optional<Intersections>
    {
      Intersections meeting;
      if (!Banded({one, other}, tolerance, meeting))
      {
        return std::nullopt;
      }
      return meeting;
    };

    /**
     *  A hair from tangent where the segment does not reach, or the arc does not: nothing.
     */
    const Curve beyond = Line{.Start = {.X = 10, .Y = 1 + hair}, .End = {.X = 20, .Y = 1 + hair}};
    const Curve over   = Line{.Start = {.X = -3, .Y = 1 + hair}, .End = {.X = 3, .Y = 1 + hair}};
    CHECK(banded(beyond, unit).value_or(Intersections{.Hits = {{}}}).Hits.empty());
    CHECK(banded(over, lower).value_or(Intersections{.Hits = {{}}}).Hits.empty());
    CHECK(banded(over, unit).value_or(Intersections{}).Hits.size() == 1);

    /**
     *  Two circles a hair apart where one arc does not reach the touch: nothing; far apart, or
     *  about one center, not a band.
     */
    const Curve beside =
        Arc{.Center = {.X = 2 + hair, .Y = 0}, .Radius = 1, .Start = {}, .Sweep = Turn};
    const Curve aside = Arc{.Center = {.X = 2 + hair, .Y = 0},
                            .Radius = 1,
                            .Start  = Degrees(30),
                            .Sweep  = Degrees(120)};
    CHECK(banded(upper, beside).value_or(Intersections{.Hits = {{}}}).Hits.empty());
    CHECK(banded(unit, aside).value_or(Intersections{.Hits = {{}}}).Hits.empty());
    CHECK(banded(unit, beside).value_or(Intersections{}).Hits.size() == 1);
    const Curve away   = Line{.Start = {.X = -3, .Y = 5}, .End = {.X = 3, .Y = 5}};
    const Curve apart  = Arc{.Center = {.X = 5, .Y = 0}, .Radius = 1, .Start = {}, .Sweep = Turn};
    const Curve around = Arc{.Center = {}, .Radius = 1 + hair, .Start = {}, .Sweep = Turn};
    CHECK(!banded(away, unit));
    CHECK(!banded(unit, apart));
    CHECK(!banded(unit, around));
    CHECK(banded(upper, lower).value_or(Intersections{.Hits = {{}}}).Hits.empty());
    const Curve near = Arc{.Center = {}, .Radius = 1, .Start = Degrees(30), .Sweep = Degrees(60)};
    /**
     *  Arcs on circles a hair apart, meeting end to end: a point, in closed form and by the
     *  thorough way; an arc about a tolerance long within another, on circles a hair apart, the
     *  thorough way's, and on one circle, a stretch shared.
     */
    const Tolerance coarse = {.Length = 1e-8, .Radians = 1e-9};
    const Curve     upper2 =
        Arc{.Center = {}, .Radius = 2, .Start = Radians(Pi), .Sweep = Radians(-Pi)};
    const Curve quarter = Arc{.Center = {.X = 5.360630966617862e-16, .Y = -4.440892098500626e-16},
                              .Radius = 1.9999999999999996,
                              .Start  = Radians(-1.5707963267948968),
                              .Sweep  = Radians(1.570796326794897)};
    CHECK(Agrees({upper2, quarter}, coarse));
    CHECK(Thoroughly(upper2, quarter, coarse).Hits.size() == 1);
    const Curve big   = Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = Degrees(90)};
    const Curve speck = Arc{.Center = {.X = 1e-10, .Y = 0},
                            .Radius = 1,
                            .Start  = Radians(Pi / 4),
                            .Sweep  = Radians(1e-9)};
    const Curve grain =
        Arc{.Center = {}, .Radius = 1, .Start = Radians(Pi / 4), .Sweep = Radians(1e-9)};
    CHECK(!banded(big, speck));
    CHECK(banded(big, grain).value_or(Intersections{}).Overlaps.size() == 1);
    const Intersections shared = banded(upper, near).value_or(Intersections{});
    CHECK(shared.Hits.empty());
    CHECK(shared.Overlaps.size() == 1);

    /**
     *  A whole circle and an arc of it starting at its seam: the arc's stretch shared, decided
     *  directly to the bit, the seam's second meeting a turn round the end of that stretch.
     */
    const Curve circle = Arc{.Center = {}, .Radius = 5, .Start = {}, .Sweep = Turn};
    const Curve half   = Arc{.Center = {}, .Radius = 5, .Start = {}, .Sweep = Degrees(180)};
    CHECK(Decided(circle, half, tolerance));
    CHECK(Identical(Intersect(circle, half, tolerance), Thoroughly(circle, half, tolerance)));
    CHECK(Intersect(circle, half, tolerance).Overlaps.size() == 1);
  }

  TEST_CASE("direct: lines through a whole circle's seam, as met thoroughly")
  {
    Mixer     mixer;
    int       differ = 0;
    int       direct = 0;
    const int rounds = Rounds(40000);
    for (int round = 0; round < rounds; ++round)
    {
      const double away   = std::ldexp(mixer.Unit(), static_cast<int>(mixer.Pick(24)));
      const double radius = std::ldexp(1.0 + mixer.Unit(), static_cast<int>(mixer.Pick(12)) - 4);
      const Arc    circle = {.Center = {.X = away, .Y = -away},
                             .Radius = radius,
                             .Start  = Degrees(std::round(mixer.Unit() * 720.0) - 360.0),
                             .Sweep  = mixer.Pick(2) == 0 ? Turn : -Turn};
      const Point  seam   = At(circle, 0.0);
      const double angle  = mixer.Unit() * 2.0 * Pi;
      const Vector way    = Vector{.X = std::cos(angle), .Y = std::sin(angle)} * radius;
      const double from   = 0.2 + (2.0 * mixer.Unit());
      const Line line = {.Start = seam - (way * from), .End = seam + (way * (0.2 + mixer.Unit()))};
      const double    near      = std::ldexp(radius, -static_cast<int>(mixer.Pick(30)) - 8);
      const Tolerance tolerance = {.Length = near, .Radians = 1e-9};
      const Curve     first     = circle;
      const Curve     other     = Kinded(mixer, line);
      for (const bool swapped : {false, true})
      {
        const Curve &one = swapped ? other : first;
        const Curve &two = swapped ? first : other;
        if (!Agrees({one, two}, tolerance))
        {
          ++differ;
          if (differ <= 5)
          {
            MESSAGE("differ: round " << round << " swapped " << swapped);
          }
        }
      }
      direct += Decided(first, other, tolerance) ? 1 : 0;
    }
    CHECK(differ == 0);
    MESSAGE("decided directly: " << direct << " of " << rounds);
    CHECK(direct > rounds / 4);
  }

  TEST_CASE("direct: arcs of two circles sharing a vertex, as met thoroughly")
  {
    Mixer     mixer;
    int       differ = 0;
    int       direct = 0;
    const int rounds = Rounds(40000);
    for (int round = 0; round < rounds; ++round)
    {
      const double away   = std::ldexp(mixer.Unit(), static_cast<int>(mixer.Pick(24)));
      const double radius = std::ldexp(1.0 + mixer.Unit(), static_cast<int>(mixer.Pick(12)) - 4);
      const Arc    first  = {.Center = {.X = away, .Y = -away},
                             .Radius = radius,
                             .Start  = Radians(mixer.Unit() * 2.0 * Pi),
                             .Sweep  = Radians((mixer.Unit() - 0.5) * 3.0 * Pi)};
      const Point  vertex = At(first, static_cast<double>(mixer.Pick(2)));
      const double turn   = mixer.Unit() * 2.0 * Pi;
      const double reach =
          radius * std::ldexp(0.5 + mixer.Unit(), static_cast<int>(mixer.Pick(8)) - 4);
      const Point     centre = vertex + (Vector{.X = std::cos(turn), .Y = std::sin(turn)} * reach);
      const Arc       second = {.Center = centre,
                                .Radius = Distance(centre, vertex),
                                .Start  = Heading(vertex - centre),
                                .Sweep  = Radians((mixer.Unit() - 0.5) * 3.0 * Pi)};
      const double    near   = std::ldexp(radius, -static_cast<int>(mixer.Pick(30)) - 8);
      const Tolerance tolerance = {.Length = near, .Radians = 1e-9};
      const Curve     one       = first;
      const Curve     two       = second;
      for (const bool swapped : {false, true})
      {
        const Curve &lhs = swapped ? two : one;
        const Curve &rhs = swapped ? one : two;
        if (!Agrees({lhs, rhs}, tolerance))
        {
          ++differ;
          if (differ <= 5)
          {
            MESSAGE("differ: round " << round << " swapped " << swapped);
          }
        }
      }
      direct += Decided(one, two, tolerance) ? 1 : 0;
    }
    CHECK(differ == 0);
    MESSAGE("decided directly: " << direct << " of " << rounds);
    CHECK(direct > rounds / 10);
  }

  TEST_CASE("direct: odd arcs, tangent circles and lost crossings go the thorough way")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Curve     across    = Line{.Start = {.X = -3, .Y = 0.3}, .End = {.X = 3, .Y = 0.3}};
    const Curve     unit      = Arc{.Center = {}, .Radius = 1, .Start = Degrees(90), .Sweep = Turn};

    /**
     *  An arc of no radius, of one past every double or of none at all, of no sweep or of more than
     *  a turn; a line of no length: each left to the thorough way, which answers the same.
     */
    const double               endless = std::numeric_limits<double>::infinity();
    const std::array<Curve, 5> odds    = {
        Curve{Arc{.Center = {}, .Radius = 0, .Start = {}, .Sweep = Turn}},
        Curve{Arc{.Center = {}, .Radius = endless, .Start = {}, .Sweep = Turn}},
        Curve{Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = Radians(0)}},
        Curve{Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = Radians(3 * Pi)}},
        Curve{Line{.Start = {.X = 0.5, .Y = 0.5}, .End = {.X = 0.5, .Y = 0.5}}}};
    for (const Curve &odd : odds)
    {
      for (const Curve *other : {&across, &unit})
      {
        CHECK(!Decided(odd, *other, tolerance));
        CHECK(!Decided(*other, odd, tolerance));
        CHECK(Identical(Intersect(odd, *other, tolerance), Thoroughly(odd, *other, tolerance)));
      }
    }

    /**
     *  An ellipse of no sweep, a point, on another's carrier, either first: no stretch shared.
     */
    const Curve whole = Ellipse{.Map = Scaling(2, 1), .Start = {}, .Sweep = Turn};
    const Curve point = Ellipse{.Map = Scaling(2, 1), .Start = Degrees(30), .Sweep = Radians(0)};
    CHECK(Intersect(whole, point, tolerance).Overlaps.empty());
    CHECK(Intersect(point, whole, tolerance).Overlaps.empty());

    /**
     *  Two circles touching from outside, their seams far from the touch: a touch, decided directly
     *  at the exact point; and a hair apart, where the carriers do not touch exactly, a touch in
     *  closed form where they come nearest.
     */
    const Curve beside =
        Arc{.Center = {.X = 2, .Y = 0}, .Radius = 1, .Start = Degrees(90), .Sweep = Turn};
    CHECK(Decided(unit, beside, tolerance));
    CHECK(Identical(Intersect(unit, beside, tolerance), Thoroughly(unit, beside, tolerance)));
    CHECK(Intersect(unit, beside, tolerance).Hits.front().Kind == Contact::Touch);
    const Curve hair = Arc{.Center = {.X = 2 + (tolerance.Length / 2), .Y = 0},
                           .Radius = 1,
                           .Start  = Degrees(90),
                           .Sweep  = Turn};
    CHECK(Decided(unit, hair, tolerance));
    CHECK(Agrees({unit, hair}, tolerance));
    const Intersections band = Intersect(unit, hair, tolerance);
    CHECK(band.Hits.size()       == 1);
    CHECK(band.Hits.front().Kind == Contact::Touch);

    /**
     *  Places at which no exact crossing lies, as no direct way gives: nothing, for the thorough
     *  way to decide.
     */
    const Curve    far     = Line{.Start = {.X = -3, .Y = 5}, .End = {.X = 3, .Y = 5}};
    const Approach nowhere = {.Along = {0.5, 0.25}, .At = {.X = 0, .Y = 5}};
    Intersections  lost;
    CHECK(!Exactly({&far, &unit}, std::span{&nowhere, 1}, lost));
    CHECK(lost.Hits.empty());
  }

  TEST_CASE("direct: lines and arcs whose boxes lie far apart meet in nothing, as met thoroughly")
  {
    /**
     *  The box a line or an arc cannot leave, and how far apart two are along x or along y.
     */
    const auto boxed = [](const Curve &curve) -> std::array<Point, 2>
    {
      if (const auto *line = std::get_if<Line>(&curve))
      {
        return {Point{.X = std::min(line->Start.X, line->End.X),
                      .Y = std::min(line->Start.Y, line->End.Y)},
                Point{.X = std::max(line->Start.X, line->End.X),
                      .Y = std::max(line->Start.Y, line->End.Y)}};
      }
      const Arc &arc = std::get<Arc>(curve);
      return {Point{.X = arc.Center.X - arc.Radius, .Y = arc.Center.Y - arc.Radius},
              Point{.X = arc.Center.X + arc.Radius, .Y = arc.Center.Y + arc.Radius}};
    };
    const auto apart = [&boxed](const Curve &one, const Curve &other) -> double
    {
      const std::array<Point, 2> box  = boxed(one);
      const std::array<Point, 2> next = boxed(other);
      return std::max(
          {box[0].X - next[1].X, next[0].X - box[1].X, box[0].Y - next[1].Y, next[0].Y - box[1].Y});
    };
    Mixer     mixer;
    int       differ = 0;
    int       remote = 0;
    int       missed = 0;
    int       uneven = 0;
    const int rounds = Rounds(40000);
    for (int round = 0; round < rounds; ++round)
    {
      const double away  = std::ldexp(mixer.Unit(), static_cast<int>(mixer.Pick(28)));
      const double size  = std::ldexp(1.0, static_cast<int>(mixer.Pick(12)) - 4);
      const double near  = std::ldexp(size, -static_cast<int>(mixer.Pick(30)) - 6);
      const Point  base  = {.X = away, .Y = -away};
      const auto   drawn = [&mixer, size](Point about) -> Curve
      {
        if (mixer.Pick(2) == 0)
        {
          return Line{.Start = about + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size),
                      .End   = about + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size)};
        }
        return Arc{.Center = about + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size),
                   .Radius = size * (0.05 + std::abs(mixer.Unit())),
                   .Start  = Degrees(180.0 * mixer.Unit()),
                   .Sweep  = mixer.Pick(4) == 0 ? Turn : Degrees(350.0 * mixer.Unit())};
      };
      /**
       *  The second drawn about the first, then moved off along x or along y until their boxes lie
       *  some tolerances apart: none, a few either side of the four that decide, or many.
       */
      const Curve                 first     = drawn(base);
      const Curve                 about     = drawn(base);
      const Tolerance             tolerance = {.Length = near, .Radians = 1e-9};
      const Tolerance             usable    = Usable({first, about}, tolerance);
      const std::array<double, 6> gaps      = {-1.0, 0.5, 3.0, 5.0, 8.0, 1000.0};
      const double                gap       = gaps.at(mixer.Pick(gaps.size())) * usable.Length;
      const bool                  level     = mixer.Pick(2) == 0;
      const double                sign      = mixer.Pick(2) == 0 ? 1.0 : -1.0;
      const std::array<Point, 2>  box       = boxed(first);
      const std::array<Point, 2>  next      = boxed(about);
      const auto                  across    = [level](Point point) -> double
      {
        return level ? point.X : point.Y;
      };
      const double reach =
          sign > 0 ? across(box[1]) - across(next[0]) : across(next[1]) - across(box[0]);
      const Vector move   = level ? Vector{.X = sign * (reach + gap), .Y = 0}
                                  : Vector{.X = 0, .Y = sign * (reach + gap)};
      Curve        second = about;
      if (auto *line = std::get_if<Line>(&second))
      {
        line->Start = line->Start + move;
        line->End   = line->End + move;
      }
      else
      {
        std::get<Arc>(second).Center = std::get<Arc>(second).Center + move;
      }
      const Tolerance moved = Usable({first, second}, tolerance);
      const bool      far   = apart(first, second) > 4.5 * moved.Length;
      const bool      close = apart(first, second) < 3.5 * moved.Length;
      for (const bool swapped : {false, true})
      {
        const Curve &one = swapped ? second : first;
        const Curve &two = swapped ? first : second;
        /**
         *  Far apart, not reachable, and nothing met; nearer, reachable within the very length the
         *  curves' grain gives, to the bit.
         */
        const std::optional<double> length = Reachable({one, two}, tolerance);
        differ += Identical(Intersect(one, two, tolerance), Thoroughly(one, two, tolerance)) || !far
                    ? 0
                    : 1;
        differ += Agrees({one, two}, tolerance) ? 0 : 1;
        missed += far && length ? 1 : 0;
        missed += close && !length ? 1 : 0;
        uneven += length && !std::equal_to<double>{}(*length, moved.Length) ? 1 : 0;
        remote += far ? 1 : 0;
      }
    }
    CHECK(differ == 0);
    CHECK(missed == 0);
    CHECK(uneven == 0);
    MESSAGE("far apart: " << remote << " of " << 2 * rounds);
    CHECK(remote > rounds / 2);

    /**
     *  Boxes apart by less than the tolerance are not far: a line ending a hair short of another
     *  still ends on it, and lines a hair apart still share their stretch.
     */
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const double    hair      = tolerance.Length / 2;
    const Curve     level     = Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const Curve upright = Line{.Start = {.X = 10 + hair, .Y = -5}, .End = {.X = 10 + hair, .Y = 5}};
    const Curve beside  = Line{.Start = {.X = 0, .Y = hair}, .End = {.X = 10, .Y = hair}};
    CHECK(Intersect(level, upright, tolerance).Hits.size()    == 1);
    CHECK(Intersect(level, beside, tolerance).Overlaps.size() == 1);
    CHECK(Identical(Intersect(level, upright, tolerance), Thoroughly(level, upright, tolerance)));
    CHECK(Identical(Intersect(level, beside, tolerance), Thoroughly(level, beside, tolerance)));
    CHECK(Reachable({level, upright}, tolerance).has_value());
    CHECK(Reachable({level, beside}, tolerance).has_value());

    /**
     *  A line or an arc with a number that is not finite has no box: far from another or not, it is
     *  left to the other ways, which find nothing, as the thorough way does.
     */
    const double none    = std::numeric_limits<double>::quiet_NaN();
    const double endless = std::numeric_limits<double>::infinity();
    const Curve  distant = Line{.Start = {.X = 100, .Y = 100}, .End = {.X = 110, .Y = 105}};
    const std::array<Curve, 4> odds = {
        Curve{Line{.Start = {.X = 0, .Y = 0}, .End = {.X = none, .Y = 1}}},
        Curve{Line{.Start = {.X = 0, .Y = 0}, .End = {.X = endless, .Y = 1}}},
        Curve{Arc{.Center = {.X = none, .Y = 0}, .Radius = 1, .Start = {}, .Sweep = Turn}},
        Curve{Arc{.Center = {.X = 0, .Y = 0}, .Radius = none, .Start = {}, .Sweep = Turn}}};
    for (const Curve &odd : odds)
    {
      CHECK(Reachable({odd, distant}, tolerance).has_value());
      CHECK(Reachable({distant, odd}, tolerance).has_value());
      CHECK(Intersect(odd, distant, tolerance).Hits.empty());
      CHECK(Identical(Intersect(odd, distant, tolerance), Thoroughly(odd, distant, tolerance)));
      CHECK(Identical(Intersect(distant, odd, tolerance), Thoroughly(distant, odd, tolerance)));
    }

    /**
     *  An arc of numbers that are finite, whose box runs past every double all the same: its size
     *  is the one the curves' grain is taken from, which leaves such sums out, and it is not told
     *  far from anything.
     */
    const Curve huge = Arc{.Center = {.X = 1.7e308, .Y = -1.7e308},
                           .Radius = 1e308,
                           .Start  = {},
                           .Sweep  = Degrees(90)};
    for (const Curves &pair : {Curves{huge, distant}, Curves{distant, huge}})
    {
      const std::optional<double> reach = Reachable(pair, tolerance);
      CHECK(reach.has_value());
      CHECK(std::equal_to<double>{}(reach.value_or(0.0), Usable(pair, tolerance).Length));
    }

    /**
     *  A line wholly to one side of another's carrier, their boxes overlapping: some tolerances
     *  from it, none of them a few either side of the four that decide, or many; either first.
     */
    int sided = 0;
    for (int round = 0; round < rounds; ++round)
    {
      const double away  = std::ldexp(mixer.Unit(), static_cast<int>(mixer.Pick(28)));
      const double size  = std::ldexp(1.0, static_cast<int>(mixer.Pick(12)) - 4);
      const double near  = std::ldexp(size, -static_cast<int>(mixer.Pick(30)) - 6);
      const Point  base  = {.X = away, .Y = -away};
      const Line   first = {.Start = base + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size),
                            .End   = base + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size)};
      const Vector way   = first.End - first.Start;
      if (!(Length(way) > size / 16))
      {
        continue;
      }
      const Curve                 one    = first;
      const Tolerance             loose  = {.Length = near, .Radians = 1e-9};
      const Tolerance             usable = Usable({one, one}, loose);
      const Vector                across = Vector{.X = -way.Y, .Y = way.X} / Length(way);
      const std::array<double, 5> gaps   = {-1.0, 0.5, 3.0, 7.0, 1000.0};
      const double                gap    = gaps.at(mixer.Pick(gaps.size())) * usable.Length;
      const double                side   = mixer.Pick(2) == 0 ? 1.0 : -1.0;
      const auto                  off    = [&](double height) -> Point
      {
        return first.Start + (way * (0.5 + mixer.Unit())) + (across * (side * height));
      };
      const Line  second = {.Start = off(gap), .End = off(gap + (size * std::abs(mixer.Unit())))};
      const Curve two    = second;
      const bool  far    = gap > 6.0 * Usable({one, two}, loose).Length;
      for (const bool swapped : {false, true})
      {
        const Curve &lhs  = swapped ? two : one;
        const Curve &rhs  = swapped ? one : two;
        differ           += Agrees({lhs, rhs}, loose) ? 0 : 1;
        missed           += far && !Decided(lhs, rhs, Usable({lhs, rhs}, loose)) ? 1 : 0;
        sided            += far ? 1 : 0;
      }
    }
    CHECK(differ == 0);
    CHECK(missed == 0);
    MESSAGE("beside: " << sided << " of " << 2 * rounds);
    CHECK(sided > rounds / 4);

    /**
     *  A line that keeps clear of an arc's whole circle meets no arc of it, their boxes
     *  overlapping: passing outside it, lodged inside it, or stopping short of it; a ray leaving it
     *  behind. So do arcs of circles that keep clear of each other, one inside the other or apart.
     *  Within a few tolerances the other ways decide, the same.
     */
    const Curve upper = Arc{.Center = {}, .Radius = 5, .Start = Degrees(30), .Sweep = Degrees(120)};
    const Curve passing  = Line{.Start = {.X = -6, .Y = 3}, .End = {.X = 2, .Y = 8}};
    const Curve lodged   = Line{.Start = {.X = -3, .Y = 1}, .End = {.X = 2, .Y = 3}};
    const Curve stopping = Line{.Start = {.X = -9, .Y = 4}, .End = {.X = -4.5, .Y = 4}};
    const Curve leaving  = Ray{.Start = {.X = 6, .Y = 4}, .Way = {.X = 1, .Y = 0.25}};
    const Curve nearly =
        Line{.Start = {.X = -9, .Y = 4}, .End = {.X = -3 - (6 * tolerance.Length), .Y = 4}};
    const Curve nested =
        Arc{.Center = {.X = 1, .Y = 0.5}, .Radius = 2, .Start = {}, .Sweep = Degrees(200)};
    const Curve aside =
        Arc{.Center = {.X = 6, .Y = 5}, .Radius = 2, .Start = Degrees(90), .Sweep = Degrees(250)};
    for (const Curve *clear : {&passing, &lodged, &stopping, &leaving, &nested, &aside})
    {
      CHECK(Decided(*clear, upper, tolerance));
      CHECK(Decided(upper, *clear, tolerance));
      CHECK(Intersect(*clear, upper, tolerance).Hits.empty());
      CHECK(Identical(Intersect(*clear, upper, tolerance), Thoroughly(*clear, upper, tolerance)));
      CHECK(Identical(Intersect(upper, *clear, tolerance), Thoroughly(upper, *clear, tolerance)));
    }
    CHECK(Agrees({nearly, upper}, tolerance));
    CHECK(Agrees({upper, nearly}, tolerance));

    /**
     *  An arc of a radius below nothing lies on the circle of that radius: far is far all the same.
     */
    const Curve inside = Arc{.Center = {}, .Radius = -2, .Start = {}, .Sweep = Degrees(90)};
    CHECK(!Reachable({inside, distant}, tolerance));
    CHECK(!Reachable({distant, inside}, tolerance));
    CHECK(Identical(Intersect(inside, distant, tolerance), Thoroughly(inside, distant, tolerance)));
  }

  TEST_CASE(
      "direct: curves of every kind whose boxes lie far apart meet in nothing, as met thoroughly")
  {
    Mixer     mixer;
    int       differ = 0;
    int       missed = 0;
    int       remote = 0;
    const int rounds = Rounds(5000);
    for (int round = 0; round < rounds; ++round)
    {
      const double away = std::ldexp(mixer.Unit(), static_cast<int>(mixer.Pick(24)));
      const double size = std::ldexp(1.0, static_cast<int>(mixer.Pick(8)) - 2);
      const double near = std::ldexp(size, -static_cast<int>(mixer.Pick(20)) - 10);
      const Point  base = {.X = away, .Y = -away};
      const auto   spot = [&mixer, size](Point about) -> Point
      {
        return about + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * size);
      };
      const auto strung = [&mixer, &spot](Point about) -> Curve
      {
        const std::array<double, 5> bulges = {0.0, 1e-7, 0.4, -1.0, 2.5};
        Polyline                    polyline;
        polyline.Closed = mixer.Pick(3) == 0;
        for (std::uint64_t index = 0; index < 2 + mixer.Pick(4); ++index)
        {
          polyline.Vertices.push_back(
              {.Position = spot(about), .Bulge = bulges.at(mixer.Pick(bulges.size()))});
        }
        return polyline;
      };
      const auto splined = [&mixer, &spot](Point about) -> Curve
      {
        Spline spline;
        spline.Degree     = 2 + static_cast<int>(mixer.Pick(2));
        const int  poles  = spline.Degree + 2 + static_cast<int>(mixer.Pick(3));
        const int  spans  = poles - spline.Degree;
        const bool weighs = mixer.Pick(2) == 0;
        for (int index = 0; index < poles; ++index)
        {
          spline.Poles.push_back(spot(about));
          if (weighs)
          {
            spline.Weights.push_back(0.25 + (2.0 * std::abs(mixer.Unit())));
          }
        }
        for (int index = -spline.Degree; index <= spans + spline.Degree; ++index)
        {
          spline.Knots.push_back(static_cast<double>(std::clamp(index, 0, spans)));
        }
        return spline;
      };
      const auto drawn = [&](Point about) -> Curve
      {
        switch (mixer.Pick(6))
        {
        case 0:
          return Bezier::Cubic{
              .Start = spot(about), .Out = spot(about), .In = spot(about), .End = spot(about)};
        case 1:
          return Bezier::Quadratic{
              .Start = spot(about), .Control = spot(about), .End = spot(about)};
        case 2:
          return strung(about);
        case 3:
          return Ellipse{.Map   = Translation(spot(about) - Point{})
                                * Rotation(Degrees(180.0 * mixer.Unit()))
                                * Scaling(size * (0.1 + std::abs(mixer.Unit())),
                                          size * (0.05 + (0.5 * std::abs(mixer.Unit())))),
                         .Start = Degrees(180.0 * mixer.Unit()),
                         .Sweep = mixer.Pick(3) == 0 ? Turn : Degrees(350.0 * mixer.Unit())};
        case 4:
          return splined(about);
        default:
          return Line{.Start = spot(about), .End = spot(about)};
        }
      };
      /**
       *  The second about the first, a little off it, or far off, whichever way.
       */
      const std::array<double, 3> shifts = {0.0, 2.0, 200.0};
      const double                shift  = size * shifts.at(mixer.Pick(shifts.size()));
      const Curve                 first  = drawn(base);
      const Curve     second = drawn(base + (Vector{.X = mixer.Unit(), .Y = mixer.Unit()} * shift));
      const Tolerance tolerance = {.Length = near, .Radians = 1e-9};
      for (const bool swapped : {false, true})
      {
        const Curve        &one       = swapped ? second : first;
        const Curve        &two       = swapped ? first : second;
        const Intersections thorough  = Thoroughly(one, two, tolerance);
        const bool          reached   = Reachable({one, two}, tolerance).has_value();
        const bool          nothing   = thorough.Hits.empty() && thorough.Overlaps.empty();
        differ                       += Agrees({one, two}, tolerance) ? 0 : 1;
        missed                       += !reached && !nothing ? 1 : 0;
        remote                       += reached ? 0 : 1;
      }
    }
    CHECK(differ == 0);
    CHECK(missed == 0);
    MESSAGE("far apart: " << remote << " of " << 2 * rounds);
    CHECK(remote > rounds / 4);

    /**
     *  A polyline's arcs bow out of its vertices' box, up to a half turn by their height and past
     *  it round their chord's middle: lines the vertices' box alone would tell far are met.
     */
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Curve     bowl      = Polyline{
        .Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1}, {.Position = {.X = 10, .Y = 0}}}};
    const Curve under   = Line{.Start = {.X = 0, .Y = -4.9}, .End = {.X = 10, .Y = -4.9}};
    const Curve balloon = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = -2.5},
                                                {.Position = {.X = 10, .Y = 0}}}};
    const Curve beside  = Line{.Start = {.X = -2, .Y = 1}, .End = {.X = -2, .Y = 12}};
    CHECK(Reachable({bowl, under}, tolerance).has_value());
    CHECK(Intersect(bowl, under, tolerance).Hits.size() == 2);
    CHECK(Reachable({beside, balloon}, tolerance).has_value());
    CHECK(Intersect(beside, balloon, tolerance).Hits.size() == 2);
    CHECK(Identical(Intersect(bowl, under, tolerance), Thoroughly(bowl, under, tolerance)));
    CHECK(Identical(Intersect(beside, balloon, tolerance), Thoroughly(beside, balloon, tolerance)));

    /**
     *  An ellipse turned reaches past its axes' ends along x; a Bézier curve's hull is more than
     *  the curve: a line in the box is not far, whether it meets the curve or not.
     */
    const Curve oval =
        Ellipse{.Map = Rotation(Degrees(45)) * Scaling(5, 1), .Start = {}, .Sweep = Turn};
    const Curve upright = Line{.Start = {.X = 3.5, .Y = -9}, .End = {.X = 3.5, .Y = 9}};
    const Curve hooked  = Bezier::Cubic{.Start = {.X = 0, .Y = 0},
                                        .Out   = {.X = 0, .Y = 9},
                                        .In    = {.X = 9, .Y = 9},
                                        .End   = {.X = 9, .Y = 0}};
    const Curve over    = Line{.Start = {.X = -1, .Y = 8}, .End = {.X = 10, .Y = 8}};
    CHECK(Reachable({oval, upright}, tolerance).has_value());
    CHECK(Intersect(oval, upright, tolerance).Hits.size() == 2);
    CHECK(Reachable({hooked, over}, tolerance).has_value());
    CHECK(Intersect(hooked, over, tolerance).Hits.empty());

    /**
     *  Far from a line: a curve of each kind is told so; one that runs on without end, a spline
     *  that is not sound, or of a weight not above nothing, and a curve of a number that is not
     *  finite are not, and are left to the ways that were theirs.
     */
    const double none    = std::numeric_limits<double>::quiet_NaN();
    const Curve  distant = Line{.Start = {.X = 500, .Y = 500}, .End = {.X = 510, .Y = 505}};
    const Spline sound   = {.Degree  = 2,
                            .Knots   = {0, 0, 0, 1, 1, 1},
                            .Poles   = {{.X = 0, .Y = 0}, {.X = 4, .Y = 6}, {.X = 8, .Y = 0}},
                            .Weights = {}};
    Spline       unsound = sound;
    unsound.Knots        = {0, 0, 0, 1, 1};
    Spline sunk          = sound;
    sunk.Weights         = {1, -1, 1};
    const Curve lone     = Polyline{.Vertices = {{.Position = {.X = 1, .Y = 1}}}};
    for (const Curve &clear : {bowl, balloon, oval, hooked, Curve{sound}, lone})
    {
      CHECK(!Reachable({clear, distant}, tolerance));
      CHECK(!Reachable({distant, clear}, tolerance));
      CHECK(Identical(Intersect(clear, distant, tolerance), Thoroughly(clear, distant, tolerance)));
      CHECK(Intersect(distant, clear, tolerance).Hits.empty());
    }
    const std::array<Curve, 8> unboxed = {
        Curve{Ray{.Start = {}, .Way = {.X = -1, .Y = 0}}},
        Curve{Xline{.Through = {}, .Way = {.X = -1, .Y = 1}}},
        Curve{unsound},
        Curve{sunk},
        Curve{Polyline{}},
        Curve{Polyline{.Vertices = {{.Position = {.X = none, .Y = 0}}, {.Position = {.X = 1}}}}},
        Curve{Bezier::Quadratic{.Start = {}, .Control = {.X = none, .Y = 1}, .End = {.X = 1}}},
        Curve{Ellipse{.Map = Translation({.X = none, .Y = 0}), .Start = {}, .Sweep = Turn}}};
    for (const Curve &odd : unboxed)
    {
      CHECK(Reachable({odd, distant}, tolerance).has_value());
      CHECK(Reachable({distant, odd}, tolerance).has_value());
    }
    for (const Curve &endless : std::span{unboxed}.first(2))
    {
      CHECK(Identical(Intersect(endless, distant, tolerance),
                      Thoroughly(endless, distant, tolerance)));
    }

    /**
     *  A piece on a line's carrier that stops eight tolerances short of the line is not far: the
     *  thorough way weighs it, and finds its pass along the carrier out of the line's reach.
     */
    const Curve lying = Bezier::Quadratic{
        .Start = {.X = 1, .Y = 0}, .Control = {.X = 2, .Y = 0}, .End = {.X = 3, .Y = 0}};
    const Curve onward =
        Line{.Start = {.X = 3 + (8 * tolerance.Length), .Y = 0}, .End = {.X = 8, .Y = 0}};
    CHECK(Reachable({lying, onward}, tolerance).has_value());
    CHECK(Intersect(lying, onward, tolerance).Overlaps.empty());
    CHECK(Intersect(onward, lying, tolerance).Overlaps.empty());
    CHECK(Identical(Intersect(lying, onward, tolerance), Thoroughly(lying, onward, tolerance)));
  }

  TEST_CASE(
      "direct: a small circle across a large arc near its end, either first, as met thoroughly")
  {
    /**
     *  The rough crossings of two circles were worked from the first: from a large one, across a
     *  small one, they were off by sixty-five tolerances, the difference of two squares all but
     *  equal, and a crossing that near an end of the large arc was taken for one past it. Intersect
     *  then missed it with the large arc first and found it with the small one first. They are
     *  worked from the smaller now.
     *
     *  A large arc, and a small whole circle through a point of it some tolerances within its end,
     *  at each of many places and leanings: 26 of these were missed.
     */
    int differ = 0;
    for (const double radius : {4.0, 12.0, 190.0})
    {
      const double    length    = std::ldexp(radius, -34);
      const Tolerance tolerance = {.Length = length, .Radians = 1e-9};
      const Point     center    = {.X = -7.25, .Y = 19.5};
      const Curve     large =
          Arc{.Center = center, .Radius = radius, .Start = Radians(1.0), .Sweep = Radians(-2.5)};
      for (int step = 0; step < 240; ++step)
      {
        const double turn = 1.0 - 2.5 + ((5.0 + (0.5 * step)) * length / radius);
        const Point  spot = center + (Vector{.X = std::cos(turn), .Y = std::sin(turn)} * radius);
        for (int lean = 0; lean < 5; ++lean)
        {
          const double reach   = length * (3500.0 + (900.0 * lean));
          const double bearing = turn + (Pi / 2) + (0.35 * (lean - 2));
          const Vector toward  = {.X = std::cos(bearing), .Y = std::sin(bearing)};
          const Curve  small   = Arc{.Center = spot - (toward * reach),
                                     .Radius = reach,
                                     .Start  = Radians(0.5 * lean),
                                     .Sweep  = Turn};
          differ +=
              Identical(Intersect(large, small, tolerance), Thoroughly(large, small, tolerance))
                  ? 0
                  : 1;
          differ +=
              Identical(Intersect(small, large, tolerance), Thoroughly(small, large, tolerance))
                  ? 0
                  : 1;
        }
      }
    }
    CHECK(differ == 0);

    /**
     *  Three pairs found so, with where they truly cross, worked out at 120 digits and rounded
     *  once: each crossing 64.5 tolerances within an end of one arc.
     */
    struct Found
    {
      Arc                Large  = {};  /**< The large arc.          */
      Arc                Small  = {};  /**< The small one.          */
      double             Length = 0.0; /**< The tolerance.          */
      std::vector<Point> Truth  = {};  /**< Where they truly cross. */
    };

    const std::array<Found, 3> found = {
        Found{.Large  = {.Center = {.X = -7.230793150903132, .Y = 19.383756956586229},
                         .Radius = 11.942374336237423,
                         .Start  = Radians(9.3121637949993037),
                         .Sweep  = Radians(-4.8835340357984052)},
              .Small  = {.Center = {.X = -10.57425786398862, .Y = 7.9189643391306745},
                         .Radius = 3.4087946024515486e-06,
                         .Start  = Radians(-0.37073004924304259),
                         .Sweep  = Radians(6.2831853071795862)},
              .Length = 4.65661e-10,
              .Truth  = {{.X = -10.57425795416571, .Y = 7.918960931529064}}},
        Found{.Large  = {.Center = {.X = -29340.83865212823, .Y = 2.6730989291456524},
                         .Radius = 4.3213465212166344,
                         .Start  = Radians(1.4487560160332347),
                         .Sweep  = Radians(-5.3834988459490329)},
              .Small  = {.Center = {.X = -29344.799361325175, .Y = 4.401335818265216},
                         .Radius = 1.8421005181072815e-06,
                         .Start  = Radians(-5.1847425344704225),
                         .Sweep  = Radians(1.8180019783940649)},
              .Length = 2.32831e-10,
              .Truth  = {{.X = -29344.79936311737, .Y = 4.401336244137444},
                         {.X = -29344.79936285628, .Y = 4.401336842499345}}},
        Found{.Large  = {.Center = {.X = 1.4119934500643563, .Y = -0.54828194473887637},
                         .Radius = 1.4577372242491575,
                         .Start  = Radians(-2.0979945655897625),
                         .Sweep  = Radians(-5.7012327278479367)},
              .Small  = {.Center = {.X = 2.8043660331803921, .Y = -0.97990360408053723},
                         .Radius = 2.9850843285732897e-07,
                         .Start  = Radians(-3.2497310765061558),
                         .Sweep  = Radians(0.66071188697346805)},
              .Length = 5.82077e-11,
              .Truth  = {{.X = 2.8043657360339047, .Y = -0.9799035755980966}}}};
    for (const Found &one : found)
    {
      const Curve     large     = one.Large;
      const Curve     small     = one.Small;
      const Tolerance tolerance = {.Length = one.Length, .Radians = 1e-9};
      for (const Curves &pair : {Curves{large, small}, Curves{small, large}})
      {
        const Intersections meeting = Intersect(pair[0], pair[1], tolerance);
        CHECK(meeting.Hits.size() == one.Truth.size());
        for (const Point truth : one.Truth)
        {
          CHECK(std::ranges::any_of(meeting.Hits,
                                    [truth](const Hit &hit) -> bool { return hit.At == truth; }));
        }
        CHECK(Identical(meeting, Thoroughly(pair[0], pair[1], tolerance)));
      }
    }
  }

  TEST_CASE("direct: lines no longer than the tolerance that share an end, as met thoroughly")
  {
    /**
     *  Two lines that share an end and are no longer than the tolerance have as many as four pairs
     *  of ends within half of it of each other, where the direct way kept room for two, and
     *  Intersect threw. A second pair now leaves them to the thorough way.
     */
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Curve     level     = Line{.Start = {}, .End = {.X = 3e-10, .Y = 0}};
    const Curve     upright   = Line{.Start = {}, .End = {.X = 0, .Y = 4e-10}};
    CHECK(Identical(Intersect(level, upright, tolerance), Thoroughly(level, upright, tolerance)));
    CHECK(Identical(Intersect(upright, level, tolerance), Thoroughly(upright, level, tolerance)));

    /**
     *  Two pairs found so, each at the tolerance it threw at.
     */
    const Point shared = {.X = 71.63196388764962, .Y = 87.45604640753515};
    const Curve longer =
        Line{.Start = shared, .End = {.X = 71.63197099681285, .Y = 87.45597351118526}};
    const Curve shorter =
        Line{.Start = shared, .End = {.X = 71.63197253636723, .Y = 87.45601082236702}};
    const Curve     third  = Line{.Start = {.X = -2.769875673679021, .Y = -1.8769958721565112},
                                  .End   = {.X = -2.769875682096332, .Y = -1.876995875170634}};
    const Curve     fourth = Line{.Start = {.X = -2.7698756778876765, .Y = -1.8769958736635726},
                                  .End   = {.X = -2.7698756852997723, .Y = -1.876995878663262}};
    const Tolerance loose  = {.Length = 1e-4, .Radians = 1e-9};
    const Tolerance fine   = {.Length = 1e-8, .Radians = 1e-9};
    CHECK(Identical(Intersect(longer, shorter, loose), Thoroughly(longer, shorter, loose)));
    CHECK(Identical(Intersect(shorter, longer, loose), Thoroughly(shorter, longer, loose)));
    CHECK(Identical(Intersect(third, fourth, fine), Thoroughly(third, fourth, fine)));
    CHECK(Identical(Intersect(fourth, third, fine), Thoroughly(fourth, third, fine)));
  }

  TEST_CASE("direct: an arc of a radius near the largest double meets a line in nothing")
  {
    /**
     *  Beside a radius near the largest double, a line is so short in the circle's frame that its
     *  length squared is nothing, and the thorough way divided a root by it: a division by nothing,
     *  which the sanitizers stop at. That root is no longer asked for.
     */
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Curve     line      = Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 10}};
    const Curve     huge      = Arc{.Center = {.X = 1.7e308, .Y = -1.7e308},
                                    .Radius = 1e308,
                                    .Start  = {},
                                    .Sweep  = Degrees(90)};
    CHECK(Intersect(line, huge, tolerance).Hits.empty());
    CHECK(Intersect(huge, line, tolerance).Hits.empty());
    CHECK(Thoroughly(line, huge, tolerance).Hits.empty());
  }

  TEST_CASE("direct: a line a step of a double from tangent to an arc touches it once")
  {
    /**
     *  A line whose ends lie one step of a double above and below the bottom of a circle of radius
     *  3: met as ready curves for their marks, the carriers are asked to touch, and no exact touch
     *  is found on them, so the band about where they come nearest is told instead: one touch, at
     *  the circle's bottom.
     */
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-10};
    const Curve     line = Line{.Start = {.X = 0x1.0f50bb4419421p+1, .Y = -0x1.7ffffffffffffp+1},
                                .End   = {.X = 0x1.7857a25df35eep+2, .Y = -0x1.8000000000001p+1}};
    const Curve     arc  = Arc{.Center = {.X = 4, .Y = 0},
                               .Radius = 3,
                               .Start  = Radians(-0x1.46dc09ec29433p+1),
                               .Sweep  = Radians(0x1.f730bd281f69bp+0)};
    const Prepared  straight(line);
    const Prepared  round(arc);
    for (const Prepared::Pair &pair :
         {Prepared::Pair{straight, round}, Prepared::Pair{round, straight}})
    {
      Marking marking;
      Intersect(pair, tolerance, marking);
      REQUIRE(marking.Marks.size() == 1);
      CHECK(marking.Marks[0].Kind  == Contact::Touch);
      CHECK(Distance(marking.Marks[0].At, Point{.X = 4, .Y = -3}) < 1e-14);
    }
  }
} // namespace Geoxy
