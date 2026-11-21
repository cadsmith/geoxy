/***************************************************************************************************
 *  @file      source/met/linear.cpp
 *  @brief     Two straight curves met directly, lines, rays and endless lines: crossing, lying
 *             along each other, or meeting end to end.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "met/linear.hpp"

#include "met/judge.hpp"
#include "met/parts.hpp"
#include "met/solve.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "straight.hpp"
#include "thorough.hpp"
#include "trig.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <variant>
#include <vector>

namespace Geoxy::Met
{
  auto Straightened(const Curve &curve) -> std::optional<Straight>
  {
    if (const auto *line = std::get_if<Line>(&curve))
    {
      return Straight{.Run    = Running(*line),
                      .Ends   = {line->Start, line->End},
                      .Count  = 2,
                      .Shape  = &curve,
                      .Places = {0.0, 1.0}};
    }
    if (const auto *ray = std::get_if<Ray>(&curve))
    {
      return Straight{.Run    = Running(*ray),
                      .Ends   = {ray->Start, {}},
                      .Count  = 1,
                      .Shape  = &curve,
                      .Places = {0.0, 0.0}};
    }
    if (const auto *xline = std::get_if<Xline>(&curve))
    {
      return Straight{
          .Run = Running(*xline), .Ends = {}, .Count = 0, .Shape = &curve, .Places = {}};
    }
    return std::nullopt;
  }

  auto Away(Point point, const Run &run) -> double
  {
    /**
     *  A way shorter than the square root of the least double squares to nothing: the run lies
     *  within that of its start, far finer than any tolerance the curves' numbers allow.
     */
    const double square = Dot(run.Way, run.Way);
    const double along  = square > 0.0 ? std::clamp(Dot(point - run.Start, run.Way) / square,
                                                    run.Along.Start, run.Along.End)
                                       : run.Along.Start;
    return Distance(point, run.Start + (run.Way * along));
  }

  namespace
  {
    /**
     *  Two straight curves, the first lying along the second over a stretch, as Aligned finds it.
     */
    struct Alongside
    {
      std::array<const Straight *, 2> Runs = {}; /**< Each as a run with its ends.          */
      Range                           Span = {}; /**< The first's stretch along the second. */
    };

    /**
     *  Two straight curves met directly where they share one stretch and nothing else: lying on one
     *  line within the tolerance over a stretch more than four tolerances long, endless where both
     *  run on, as the thorough way finds it, and each end of each either an end of that stretch or
     *  farther than four tolerances from the other. There the stretch's ends are the curves' own,
     *  or at infinity, where it rests as found, and every end within the tolerance of the other is
     *  on it; so the thorough way's answer is that stretch alone, worked out by the same steps.
     *
     *  @param  lying     The two straight curves, each of some length, and the stretch of the first
     *                    that lies along the second.
     *  @param  tolerance The tolerance, no finer than their grain.
     *  @param  into      Where what they meet in goes, empty on entry.
     *  @return           Whether it was decided so.
     */
    auto Collinearly(const Alongside &lying, const Tolerance &tolerance, Intersections &into)
        -> bool
    {
      const Straight &first  = *lying.Runs[0];
      const Straight &second = *lying.Runs[1];
      const Run      &one    = first.Run;
      const Run      &other  = second.Run;
      const Range    &span   = lying.Span;
      const double    clear  = 4.0 * tolerance.Length;
      if (!((span.End - span.Start) * Length(one.Way) > clear))
      {
        return false;
      }
      const double square = Dot(one.Way, one.Way);
      const double base   = Dot(other.Start - one.Start, one.Way) / square;
      const double rate   = Dot(other.Way, one.Way) / square;
      const auto   bounds = [&span](double along) -> bool
      {
        return Exactly(along, span.Start) || Exactly(along, span.End);
      };
      for (std::size_t index = 0; index < first.Count; ++index)
      {
        if (!bounds(first.Places.at(index)) && !(Away(first.Ends.at(index), other) > clear))
        {
          return false;
        }
      }
      for (std::size_t index = 0; index < second.Count; ++index)
      {
        /**
         *  An end of the second bounds the stretch only where it lies within the tolerance of the
         *  first, as the thorough way takes its ends there; half of it, here, to be sure.
         */
        const double along = base + (rate * second.Places.at(index));
        const double away  = Away(second.Ends.at(index), one);
        const bool   bound = bounds(along) && away <= tolerance.Length / 2;
        if (!bound && !(away > clear))
        {
          return false;
        }
      }
      const Found found = Against(one, other, tolerance);
      into.Overlaps.assign(found.Shared.begin(), found.Shared.end());
      return true;
    }

    /**
     *  Two straight curves met: each as a curve and as a run with its ends, the tolerance, and the
     *  cross product of their ways.
     */
    struct Straights
    {
      std::array<const Straight *, 2> Runs      = {};  /**< Each as a run with its ends. */
      std::array<const Curve *, 2>    Curves    = {};  /**< Each as a curve.             */
      Geoxy::Tolerance                Tolerance = {};  /**< The tolerance.               */
      double                          Cross     = 0.0; /**< The cross product of ways.   */
    };

    /**
     *  Two straight curves met directly where an end of one or both lies on the other: crossing at
     *  a sine above Steep, every end either farther than four tolerances from the other curve or
     *  within half the tolerance of where their carriers cross, and that crossing, worked out
     *  exactly, on both. There the thorough way finds that crossing and the ends at it, merges them
     *  into one hit at an end, and polishes it to the same exact crossing.
     *
     *  @param  straights The two curves, as curves and runs, the tolerance and how they cross.
     *  @param  into      Where what they meet in goes, empty on entry.
     *  @return           Whether it was decided so; else the thorough way must decide.
     */
    auto Endly(const Straights &straights, Intersections &into) -> bool
    {
      const Run   &one     = straights.Runs[0]->Run;
      const Run   &other   = straights.Runs[1]->Run;
      const double lengths = Length(one.Way) * Length(other.Way);
      /**
       *  Two ways each shorter than the square root of the least double have lengths whose product
       *  is nothing, and part at no angle a double tells: they are left to the thorough way.
       */
      if (!(lengths > 0.0))
      {
        return false;
      }
      const double sine   = straights.Cross / lengths;
      const double length = straights.Tolerance.Length;
      if (!(std::abs(sine) > Steep))
      {
        return false;
      }
      const Vector gap   = other.Start - one.Start;
      const Point  point = one.Start + (one.Way * (Cross(gap, other.Way) / straights.Cross));
      const auto   ended = [&point, length](const Straight &ends, const Run &run) -> bool
      {
        return std::ranges::all_of(
            std::span{ends.Ends}.first(ends.Count), [&run, &point, length](Point end) -> bool
            { return Away(end, run) > 4.0 * length || Distance(end, point) <= length / 2.0; });
      };
      if (!ended(*straights.Runs[0], other) || !ended(*straights.Runs[1], one))
      {
        return false;
      }
      const std::array near = {Approach{
          .Along = {Cross(gap, other.Way) / straights.Cross, Cross(gap, one.Way) / straights.Cross},
          .At    = point}};
      if (!Exactly(straights.Curves, near, into))
      {
        return false;
      }
      into.Hits.front().Kind = Contact::End;
      return true;
    }

    /**
     *  Two straight curves met directly where they cross so shallowly that they stay within the
     *  tolerance of each other a long way either side: within twice the tolerance's angle, neither
     *  parallel to the other exactly nor lying along it, and every end farther than twice the
     *  tolerance from the other, so farther than the tolerance whatever the rounding, which keeps
     *  the second from lying along the first too, as such a stretch ends at an end within the
     *  tolerance of the other. Where curves so placed come within the tolerance at all, the stretch
     *  they stay within it over ends where the carriers part, on both sides, as an end within it
     *  would lie within the tolerance of the other; it is the band about a shallow crossing, one
     *  hit, at the crossing of the carriers worked out exactly, which then lies on both. That hit
     *  is a touch where they part by no more than the angle and a crossing where they part by more,
     *  as the thorough way tells them.
     *
     *  The thorough way works that crossing out in doubles first, and along lines so near parallel
     *  its two places land too far apart to be within the tolerance: it dropped the hit and found
     *  nothing, where two endless lines must meet. Where the exact crossing lies off either curve,
     *  or so far out that the doubles there are coarser than the tolerance, the thorough way
     *  decides.
     *
     *  @param  straights The two curves, as curves and runs, the tolerance and how they cross.
     *  @param  span      The stretch of the first that lies along the second, if any.
     *  @param  into      Where what they meet in goes, empty on entry.
     *  @return           Whether it was decided so; else the thorough way must decide.
     */
    auto Shallowly(const Straights &straights, const std::optional<Range> &span,
                   Intersections &into) -> bool
    {
      const Run   &one    = straights.Runs[0]->Run;
      const Run   &other  = straights.Runs[1]->Run;
      const double cross  = straights.Cross;
      const double length = straights.Tolerance.Length;
      const auto   clear  = [length](const Straight &ends, const Run &run) -> bool
      {
        return std::ranges::all_of(std::span{ends.Ends}.first(ends.Count),
                                   [&run, length](Point end) -> bool
                                   { return Away(end, run) > 2.0 * length; });
      };
      if (Zero(cross)
          || span
          || !clear(*straights.Runs[0], other)
          || !clear(*straights.Runs[1], one))
      {
        return false;
      }
      const Vector     gap   = other.Start - one.Start;
      const double     place = Cross(gap, other.Way) / cross;
      const std::array near  = {Approach{.Along = {place, Cross(gap, one.Way) / cross},
                                         .At    = one.Start + (one.Way * place)}};
      /**
       *  Every hit's two places lie within the tolerance of each other; a crossing so far out that
       *  the doubles there are coarser has none such, and is left to the thorough way.
       */
      if (!Exactly(straights.Curves, near, into)
          || !(Distance(into.Hits.front().First.At, into.Hits.front().Second.At) <= length))
      {
        into.Hits.clear();
        return false;
      }
      /**
       *  A crossing where each passes to a side of the other beyond the angle, as Judged asks of
       *  one both ways round; a touch otherwise.
       */
      const double radians   = straights.Tolerance.Radians;
      const int    sides     = Passing({one.Way, one.Way}, other.Way, radians)
                             * Passing({other.Way, other.Way}, one.Way, radians);
      into.Hits.front().Kind = sides != 0 ? Contact::Cross : Contact::Touch;
      return true;
    }

    /**
     *  Whether a line lies wholly to one side of a run's carrier, farther from it than a distance:
     *  both its ends on one side, the cross product of each from the run's start with the run's way
     *  past that distance times the way's length.
     *
     *  The cross products are off by thirteen steps of a double times the way's length at the most,
     *  and the run's end lies within two steps of the carrier its rounded way gives. Farther than
     *  four tolerances, of sixteen steps each at the least, no point of the line is within three
     *  tolerances of any point of a line that runs along that carrier, so the two meet in nothing,
     *  as two curves whose boxes lie as far apart do.
     *
     *  @param  line  The line.
     *  @param  run   The run whose carrier.
     *  @param  least The distance times the length of the run's way.
     *  @return       Whether it does, never where a number is not finite.
     */
    auto Beside(const Line &line, const Run &run, double least) -> bool
    {
      const double start  = Cross(line.Start - run.Start, run.Way);
      const double end    = Cross(line.End - run.Start, run.Way);
      const bool   sided  = (start > 0.0) == (end > 0.0);
      const double nearer = std::min(std::abs(start), std::abs(end));
      return sided && nearer > least && std::isfinite(start + end);
    }

  } // namespace

  auto Straightly(const Curves &curves, const Tolerance &tolerance, Intersections &into) -> bool
  {
    const std::optional<Straight> first  = Straightened(curves[0]);
    const std::optional<Straight> second = Straightened(curves[1]);
    if (!first || !second)
    {
      return false;
    }
    const Run   &one    = first->Run;
    const Run   &other  = second->Run;
    const double length = Length(one.Way);
    const double reach  = Length(other.Way);
    if (!(length > 0.0) || !(reach > 0.0))
    {
      return false;
    }
    const auto  *line   = std::get_if<Line>(&curves[0].get());
    const auto  *across = std::get_if<Line>(&curves[1].get());
    const double clear  = 4.0 * tolerance.Length;
    if (line != nullptr
        && across != nullptr
        && (Beside(*across, one, clear * length) || Beside(*line, other, clear * reach)))
    {
      return true;
    }
    /**
     *  Where the first lies along the second is asked once, for the stretch they may share here and
     *  for the crossing below, which they have only where neither lies along the other.
     */
    const std::optional<Range> span = Aligned(one, other, tolerance);
    if (span && Collinearly({.Runs = {&*first, &*second}, .Span = *span}, tolerance, into))
    {
      return true;
    }
    const auto far = [clear](const Straight &ends, const Run &run) -> bool
    {
      return std::ranges::all_of(std::span{ends.Ends}.first(ends.Count),
                                 [&run, clear](Point end) -> bool
                                 { return Away(end, run) > clear; });
    };
    const double    cross     = Cross(one.Way, other.Way);
    const Straights straights = {.Runs      = {&*first, &*second},
                                 .Curves    = {&curves[0].get(), &curves[1].get()},
                                 .Tolerance = tolerance,
                                 .Cross     = cross};
    /**
     *  An angle whose sine is below 0 takes in no pair of lines, and two parallel exactly cross
     *  nowhere: neither is crossed below.
     */
    const bool shallow =
        Zero(cross) || !(std::abs(cross) > 2.0 * Trig::Sin(tolerance.Radians) * length * reach);
    if (!far(*first, other) || !far(*second, one))
    {
      return Endly(straights, into) || (shallow && Shallowly(straights, span, into));
    }
    if (shallow)
    {
      const bool apart = line != nullptr && across != nullptr && Apart(*line, *across) > clear;
      return apart || Shallowly(straights, span, into);
    }
    if (span || Aligned(other, one, tolerance))
    {
      return false;
    }
    const Vector gap    = other.Start - one.Start;
    const double place  = Cross(gap, other.Way) / cross;
    const double beyond = Cross(gap, one.Way) / cross;
    const bool   within = one.Along.Start   <= place
                       && place             <= one.Along.End
                       && other.Along.Start <= beyond
                       && beyond            <= other.Along.End;
    if (!within)
    {
      return true;
    }
    return Exactly(
        {&curves[0].get(), &curves[1].get()},
        std::array{Approach{.Along = {place, beyond}, .At = one.Start + (one.Way * place)}}, into);
  }
} // namespace Geoxy::Met
