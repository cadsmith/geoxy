/***************************************************************************************************
 *  @file      source/met/round.cpp
 *  @brief     A straight curve and an arc, or two arcs, met directly: their crossings found
 *             roughly, and arcs of one center or one circle.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "met/round.hpp"

#include "elliptic.hpp"
#include "met/ellipse.hpp"
#include "met/linear.hpp"
#include "met/overlap.hpp"
#include "met/parts.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "straight.hpp"
#include "thorough.hpp"
#include "trig.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy::Met
{
  namespace
  {
    /**
     *  The sweep of an arc the direct ways take, in radians: an arc of a radius, sweeping at most
     *  once round; nothing for any other. No cut is worked out, so that what a circle's center and
     *  radius decide alone asks for no sine.
     *
     *  @param  arc The arc.
     *  @return     Its sweep.
     */
    auto Plain(const Arc &arc) -> std::optional<double>
    {
      const double sweep = arc.Sweep.InRadians();
      const double size  = std::abs(sweep);
      const bool   plain =
          arc.Radius > 0.0 && std::isfinite(arc.Radius) && size > 0.0 && size <= 2.0 * Pi;
      return plain ? std::optional<double>{sweep} : std::nullopt;
    }

  } // namespace

  auto Rounded(const Arc &arc, double sweep) -> Round
  {
    const auto cut = [&arc](Vector way) -> Point
    {
      return arc.Center + (way * arc.Radius);
    };
    return Round{.Center = arc.Center,
                 .Radius = arc.Radius,
                 .Sweep  = sweep,
                 .Cuts   = {cut(Trig::Way(arc.Start)), cut(Trig::Way(arc.Start, arc.Sweep))}};
  }

  auto Rounded(const Arc &arc) -> std::optional<Round>
  {
    const std::optional<double> sweep = Plain(arc);
    if (!sweep)
    {
      return std::nullopt;
    }
    return Rounded(arc, *sweep);
  }

  auto Off(Point point, const Round &round) -> double
  {
    return std::abs(Distance(point, round.Center) - round.Radius);
  }

  auto Swept(const Round &round, Point point) -> bool
  {
    const double size = std::abs(round.Sweep);
    if (size >= 2.0 * Pi)
    {
      return true;
    }
    const bool   forward = round.Sweep > 0.0;
    const Vector start   = (forward ? round.Cuts[0] : round.Cuts[1]) - round.Center;
    const Vector end     = (forward ? round.Cuts[1] : round.Cuts[0]) - round.Center;
    const Vector toward  = point - round.Center;
    if (size <= Pi)
    {
      return Cross(start, toward) >= 0.0 && Cross(toward, end) >= 0.0;
    }
    const bool rest = Cross(end, toward) > 0.0 && Cross(toward, start) > 0.0;
    return !rest;
  }

  namespace
  {
    /**
     *  The ends of two curves settled against their carriers' crossings: each end either farther
     *  than four tolerances from the other curve's carrier, or within half the tolerance of a
     *  crossing, where it is kept to make that crossing a hit at an end. A closed curve's seam must
     *  lie far, as the thorough way cuts it there. The crossings lie farther apart than the four
     *  tolerances, so that an end is near one of them at most.
     *
     *  @param  ends      Each end, with its distance from the other curve's carrier, and whether it
     *                    is a seam.
     *  @param  crossings The crossings.
     *  @param  tolerance The tolerance.
     *  @return           The ends at crossings; nothing where an end is neither far nor at one.
     */
    auto Settle(std::span<const std::pair<Point, double>> ends, const Carried &crossings,
                const Tolerance &tolerance) -> std::optional<Nearby>
    {
      const double          clear = 4.0 * tolerance.Length;
      std::optional<Nearby> near{std::in_place};
      for (const auto &[end, off] : ends)
      {
        if (off > clear)
        {
          continue;
        }
        const bool atCrossing = std::ranges::any_of(
            std::span{crossings.At}.first(crossings.Count), [end, &tolerance](Point point) -> bool
            { return Distance(end, point) <= tolerance.Length / 2; });
        if (!atCrossing)
        {
          return std::nullopt;
        }
        near->End(end);
      }
      return near;
    }

    /**
     *  Where two curves come within four tolerances at one place only, touching there if at all:
     *  each end of either farther than four tolerances from the other's carrier, so that none lies
     *  where they stay near each other; a whole circle's seam, not an end, may lie anywhere.
     *  Nothing comes of them where that place lies off either; else one touch where the carriers
     *  touch exactly, the thorough way's band about it polished to that point; and where they do
     *  not touch exactly, they are left to the thorough way.
     *
     *  @param  point  The place, roughly, its place on a straight first curve with it.
     *  @param  ends   Each end, and its distance from the other curve's carrier.
     *  @param  placed Whether the place lies on both curves, and four tolerances.
     *  @return        The touch to ask for; nothing where the thorough way must decide.
     */
    auto Singly(const Approach &point, std::span<const std::pair<Point, double>> ends,
                std::pair<bool, double> placed) -> std::optional<Nearby>
    {
      const auto &[both, clear] = placed;
      const bool far =
          std::ranges::all_of(ends, [clear](const std::pair<Point, double> &end) -> bool
                              { return end.second > clear; });
      if (!far)
      {
        return std::nullopt;
      }
      Nearby near;
      if (both)
      {
        Approach touch = point;
        touch.Touch    = true;
        near.Add(touch);
      }
      return near;
    }

    /**
     *  A straight curve near tangent to an arc's circle, met directly where it can be: at the foot
     *  of the center on its carrier, as Singly decides.
     *
     *  @param  straight  The straight curve, its run and its ends.
     *  @param  round     The arc, held as points.
     *  @param  tolerance The tolerance.
     *  @return           The touch to ask for; nothing where the thorough way must decide.
     */
    auto Touching(const Straight &straight, const Round &round, const Tolerance &tolerance)
        -> std::optional<Nearby>
    {
      const Run   &run   = straight.Run;
      const double along = Dot(round.Center - run.Start, run.Way) / Dot(run.Way, run.Way);
      const Point  foot  = run.Start + (run.Way * along);
      std::array<std::pair<Point, double>, 4> ends  = {};
      std::size_t                             count = 0;
      for (const Point end : std::span{straight.Ends}.first(straight.Count))
      {
        ends.at(count) = {end, Off(end, round)};
        ++count;
      }
      for (const Point cut : std::span{round.Cuts}.first(std::abs(round.Sweep) >= 2.0 * Pi ? 0 : 2))
      {
        ends.at(count) = {cut, Away(cut, run)};
        ++count;
      }
      const bool onBoth = run.Along.Start < along && along < run.Along.End && Swept(round, foot);
      return Singly({.Along = {along, 0.0}, .At = foot}, std::span{ends}.first(count),
                    {onBoth, 4.0 * tolerance.Length});
    }

    /**
     *  Two circles near tangent, met directly where they can be: at the point of the larger on the
     *  line of centers toward the smaller, or toward each other where neither is inside, as Singly
     *  decides.
     *
     *  @param  rounds    The two arcs, held as points, of circles some way apart.
     *  @param  tolerance The tolerance.
     *  @return           The touch to ask for; nothing where the thorough way must decide.
     */
    auto Kissing(const std::array<const Round *, 2> &rounds, const Tolerance &tolerance)
        -> std::optional<Nearby>
    {
      const Round &one   = *rounds[0];
      const Round &other = *rounds[1];
      const double apart = Distance(one.Center, other.Center);
      const double clear = 4.0 * tolerance.Length;
      if (!(apart > clear))
      {
        return std::nullopt;
      }
      const Vector toward = (other.Center - one.Center) / apart;
      const bool   inside = std::abs(apart - (one.Radius + other.Radius)) > clear;
      const Point  point  = inside && other.Radius > one.Radius
                              ? other.Center - (toward * other.Radius)
                              : one.Center + (toward * one.Radius);
      std::array<std::pair<Point, double>, 4> ends  = {};
      std::size_t                             count = 0;
      for (const auto &[round, against] : {std::pair{&one, &other}, std::pair{&other, &one}})
      {
        for (const Point cut :
             std::span{round->Cuts}.first(std::abs(round->Sweep) >= 2.0 * Pi ? 0 : 2))
        {
          ends.at(count) = {cut, Off(cut, *against)};
          ++count;
        }
      }
      return Singly({.Along = {0.0, 0.0}, .At = point}, std::span{ends}.first(count),
                    {Swept(one, point) && Swept(other, point), clear});
    }

    /**
     *  How a straight curve's carrier passes a circle, roughly, and how far from it is clear of it.
     */
    struct Passage
    {
      double  Off       = 0.0; /**< How far the carrier passes from the center. */
      Carried Crossings = {};  /**< Where it crosses the circle, if it does.    */
      double  Wide      = 0.0; /**< The distance from the circle that is clear. */
    };

    /**
     *  Whether a straight curve keeps clear of an arc's whole circle, and so meets no arc of it,
     *  whatever its sweep, with no cut worked out: its carrier passing farther outside than the
     *  clear distance; or neither crossing within its stretch and each of its ends farther from the
     *  circle than that, the nearest it comes to the circle being then at an end.
     *
     *  At eight tolerances clear, every cut of the arc is farther than four from the curve, and no
     *  end of the curve near a crossing: what the direct way would find is nothing.
     *
     *  @param  straight The straight curve, its run and its ends.
     *  @param  arc      The arc whose circle.
     *  @param  passing  How the carrier passes the circle, and the clear distance.
     *  @return          Whether it keeps clear.
     */
    auto Aloof(const Straight &straight, const Arc &arc, const Passage &passing) -> bool
    {
      const Carried &crossings = passing.Crossings;
      if (crossings.Count == 0)
      {
        return passing.Off - arc.Radius > passing.Wide;
      }
      const Run &run    = straight.Run;
      const auto within = [&run](double along) -> bool
      {
        return run.Along.Start < along && along < run.Along.End;
      };
      const auto distant = [&arc, &passing](Point end) -> bool
      {
        return std::abs(Distance(end, arc.Center) - arc.Radius) > passing.Wide;
      };
      return !within(crossings.Along[0])
          && !within(crossings.Along[1])
          && std::ranges::all_of(std::span{straight.Ends}.first(straight.Count), distant);
    }

  } // namespace

  auto Roundly(const Straight &straight, const Arc &arc, const Tolerance &tolerance)
      -> std::optional<Nearby>
  {
    const Run                  &run    = straight.Run;
    const double                clear  = 4.0 * tolerance.Length;
    const double                square = Dot(run.Way, run.Way);
    const std::optional<double> sweep  = Plain(arc);
    if (!(square > 0.0) || !sweep || !(arc.Radius > 16.0 * clear))
    {
      return std::nullopt;
    }
    const double radius = arc.Radius;
    const double middle = Dot(arc.Center - run.Start, run.Way) / square;
    const double off    = Distance(run.Start + (run.Way * middle), arc.Center);
    const double ratio  = off / radius;
    if (!(std::abs(off - radius) > clear))
    {
      return Touching(straight, Rounded(arc, *sweep), tolerance);
    }
    const bool steep = off > radius || (1.0 - (ratio * ratio)) > Steep * Steep;
    if (!steep)
    {
      return std::nullopt;
    }
    Carried crossings;
    if (off < radius)
    {
      const double half = std::sqrt((radius * radius) - (off * off)) / std::sqrt(square);
      crossings         = {
          .At = {run.Start + (run.Way * (middle - half)), run.Start + (run.Way * (middle + half))},
          .Along = {middle - half, middle + half},
          .Count = 2};
    }
    if (Aloof(straight, arc, {.Off = off, .Crossings = crossings, .Wide = 2.0 * clear}))
    {
      return std::optional<Nearby>{std::in_place};
    }
    const Round                             round = Rounded(arc, *sweep);
    const bool                              whole = std::abs(round.Sweep) >= 2.0 * Pi;
    std::array<std::pair<Point, double>, 4> ends  = {};
    std::size_t                             count = 0;
    for (const Point end : std::span{straight.Ends}.first(straight.Count))
    {
      ends.at(count) = {end, Off(end, round)};
      ++count;
    }
    for (const Point cut : round.Cuts)
    {
      const double away = Away(cut, run);
      if (whole && !(away > clear))
      {
        /**
         *  A whole circle's seam is no end: on a crossing it is that crossing, placed at the seam's
         *  first pass; near the line elsewhere, the thorough way's.
         */
        const bool seamed = std::ranges::any_of(
            std::span{crossings.At}.first(crossings.Count), [cut, &tolerance](Point point) -> bool
            { return Distance(cut, point) <= tolerance.Length / 2; });
        if (!seamed)
        {
          return std::nullopt;
        }
        continue;
      }
      ends.at(count) = {cut, away};
      ++count;
    }
    std::optional<Nearby> near = Settle(std::span{ends}.first(count), crossings, tolerance);
    if (!near)
    {
      return std::nullopt;
    }
    const auto atStraight = std::span{ends}.first(straight.Count);
    const auto atArc      = std::span{ends}.first(count).subspan(straight.Count);
    for (std::size_t index = 0; index < crossings.Count; ++index)
    {
      const Point  point = crossings.At.at(index);
      const double along = crossings.Along.at(index);
      const auto ended = [point, &tolerance](std::span<const std::pair<Point, double>> some) -> bool
      {
        return std::ranges::any_of(some, [point, &tolerance](const auto &end) -> bool
                                   { return Distance(end.first, point) <= tolerance.Length / 2; });
      };
      const bool onStraight =
          (run.Along.Start < along && along < run.Along.End) || ended(atStraight);
      if (onStraight && (Swept(round, point) || ended(atArc)))
      {
        near->Add({.Along = {along, 0.0}, .At = point});
      }
    }
    return near;
  }

  void Seamwise(const Arc &circle, const Tolerance &tolerance, Intersections &into)
  {
    const Point seam  = At(circle, 0.0);
    auto        front = into.Hits.begin();
    for (auto hit = into.Hits.begin(); hit != into.Hits.end(); ++hit)
    {
      if (Distance(hit->First.At, seam) <= tolerance.Length)
      {
        std::rotate(front, hit, std::next(hit));
        ++front;
      }
    }
  }

  auto Cuts(const Round &round) -> std::span<const Point>
  {
    return std::span{round.Cuts}.first(std::abs(round.Sweep) >= 2.0 * Pi ? 0 : 2);
  }

  auto Abutting(const std::array<const Round *, 2> &rounds, const Tolerance &tolerance)
      -> std::optional<Abutted>
  {
    Abutted             abutted;
    std::array<bool, 2> taken = {};
    const auto          ends  = Cuts(*rounds[0]);
    const auto          other = Cuts(*rounds[1]);
    for (std::size_t one = 0; one < ends.size(); ++one)
    {
      bool paired = false;
      for (std::size_t two = 0; two < other.size(); ++two)
      {
        const double apart = Distance(ends[one], other[two]);
        if (apart <= tolerance.Length / 2)
        {
          if (paired || taken.at(two))
          {
            return std::nullopt;
          }
          abutted.Pairs.at(abutted.Count) = {one, two};
          ++abutted.Count;
          paired        = true;
          taken.at(two) = true;
        }
        else if (!(apart > 4.0 * tolerance.Length))
        {
          return std::nullopt;
        }
      }
    }
    return abutted;
  }

  auto Concentrically(const std::array<const Arc *, 2> &arcs, const Tolerance &tolerance,
                      Intersections &into) -> bool
  {
    const Arc &first  = *arcs[0];
    const Arc &second = *arcs[1];
    if (!(first.Center == second.Center && Exactly(first.Radius, second.Radius)))
    {
      return false;
    }
    const std::optional<Round> one   = Rounded(first);
    const std::optional<Round> other = one ? Rounded(second) : std::nullopt;
    if (!one || !other)
    {
      return false;
    }
    const double clear = 4.0 * tolerance.Length;
    const auto   whole = [](const Round &round) -> bool
    {
      return std::abs(round.Sweep) >= 2.0 * Pi;
    };
    /**
     *  A whole circle against an arc of it shares the arc's stretch, cut nowhere where the circle's
     *  seam lies outside the arc and far from its ends, or at an end.
     */
    const auto seamless = [&whole, &tolerance, clear](const Round &round,
                                                      const Round &against) -> bool
    {
      const Point  seam = round.Cuts[0];
      const double near =
          std::min(Distance(seam, against.Cuts[0]), Distance(seam, against.Cuts[1]));
      const bool off = near > clear && !Swept(against, seam);
      return !whole(round) || (!whole(against) && (near <= tolerance.Length / 2 || off));
    };
    const std::optional<Abutted> abutted = Abutting({&*one, &*other}, tolerance);
    if (!seamless(*one, *other) || !seamless(*other, *one) || !abutted)
    {
      return false;
    }
    const Found found = Overlapped(Elliptic(first), Elliptic(second), 0.0);
    if (abutted->Count > 0)
    {
      return false;
    }
    List<Overlap> shared;
    for (const Overlap &overlap : found.Shared)
    {
      shared.push_back(
          {.First  = {.Start = 0.0 + overlap.First.Start, .End = 0.0 + overlap.First.End},
           .Second = {.Start = 0.0 + overlap.Second.Start, .End = 0.0 + overlap.Second.End},
           .Along  = overlap.Along});
    }
    const List<Overlap> joined = Joined(std::move(shared));
    into.Overlaps.assign(joined.begin(), joined.end());
    /**
     *  A point the sweeps meet at alone lies where ends meet: ends of the two arcs, met in closed
     *  form, not here; or a whole circle's seam at the other's end, where the sweeps meet once more
     *  a whole turn round, at the end of the stretch they share.
     */
    return true;
  }

  namespace
  {
    /**
     *  The ends of two arcs, each with its distance from the other's circle: the first's two, then
     *  the second's; nothing where a whole circle's seam lies within four tolerances of the other.
     *
     *  @param  rounds The two arcs, held as points.
     *  @param  clear  Four tolerances.
     *  @return        The ends.
     */
    auto Cut(const std::array<const Round *, 2> &rounds, double clear)
        -> std::optional<std::array<std::pair<Point, double>, 4>>
    {
      std::array<std::pair<Point, double>, 4> ends = {};
      for (std::size_t side = 0; side < 2; ++side)
      {
        const Round &round   = *rounds.at(side);
        const Round &against = *rounds.at(1 - side);
        for (std::size_t cut = 0; cut < 2; ++cut)
        {
          const Point  end  = round.Cuts.at(cut);
          const double away = Off(end, against);
          if (std::abs(round.Sweep) >= 2.0 * Pi && !(away > clear))
          {
            return std::nullopt;
          }
          ends.at((2 * side) + cut) = {end, away};
        }
      }
      return ends;
    }

  } // namespace

  auto Circularly(const Arc &first, const Arc &second, const Tolerance &tolerance)
      -> std::optional<Nearby>
  {
    const double                clear = 4.0 * tolerance.Length;
    const std::optional<double> sweep = Plain(first);
    const std::optional<double> swept = sweep ? Plain(second) : std::nullopt;
    if (!sweep || !swept || !(first.Radius > 16.0 * clear) || !(second.Radius > 16.0 * clear))
    {
      return std::nullopt;
    }
    const double inner = first.Radius;
    const double outer = second.Radius;
    const double apart = Distance(first.Center, second.Center);
    const double sum   = inner + outer;
    const double gap   = std::abs(inner - outer);
    /**
     *  Circles that keep eight tolerances clear of each other, apart or one inside the other, meet
     *  in nothing, whatever their arcs' sweeps, and no cut is worked out: every cut is then farther
     *  than four tolerances from the other circle.
     */
    const bool crossed = apart < sum && apart > gap;
    if (!crossed && std::min(std::abs(apart - sum), std::abs(apart - gap)) > 2.0 * clear)
    {
      return std::optional<Nearby>{std::in_place};
    }
    const Round one   = Rounded(first, *sweep);
    const Round other = Rounded(second, *swept);
    if (!(std::abs(apart - sum) > clear) || !(std::abs(apart - gap) > clear))
    {
      return Kissing({&one, &other}, tolerance);
    }
    Carried crossings;
    if (apart < sum && apart > gap)
    {
      const double cosine =
          ((inner * inner) + (outer * outer) - (apart * apart)) / (2.0 * inner * outer);
      if (!((1.0 - (cosine * cosine)) > Steep * Steep))
      {
        return std::nullopt;
      }
      /**
       *  The crossings are worked from the smaller circle. From the larger, a chord much shorter
       *  than its radius is the difference of two squares all but equal: off by more than the
       *  tolerance, so that a crossing near an end of an arc was taken for one past it.
       */
      const Round &small  = inner <= outer ? one : other;
      const Round &large  = inner <= outer ? other : one;
      const double less   = small.Radius;
      const double more   = large.Radius;
      const Vector toward = (large.Center - small.Center) / apart;
      const double along  = ((apart * apart) + (less * less) - (more * more)) / (2.0 * apart);
      const double height = std::sqrt(std::max(0.0, (less * less) - (along * along)));
      const Point  base   = small.Center + (toward * along);
      const Vector rise   = Perpendicular(toward) * height;
      crossings           = {.At = {base + rise, base - rise}, .Along = {}, .Count = 2};
    }
    const std::optional<std::array<std::pair<Point, double>, 4>> ends = Cut({&one, &other}, clear);
    if (!ends)
    {
      return std::nullopt;
    }
    std::optional<Nearby> near = Settle(*ends, crossings, tolerance);
    if (!near)
    {
      return std::nullopt;
    }
    for (const Point point : std::span{crossings.At}.first(crossings.Count))
    {
      const auto ended = [point, &tolerance](std::span<const std::pair<Point, double>> some) -> bool
      {
        return std::ranges::any_of(some, [point, &tolerance](const auto &end) -> bool
                                   { return Distance(end.first, point) <= tolerance.Length / 2; });
      };
      const bool onOne = Swept(one, point) || ended(std::span{*ends}.first(2));
      if (onOne && (Swept(other, point) || ended(std::span{*ends}.subspan(2))))
      {
        near->Add({.Along = {0.0, 0.0}, .At = point});
      }
    }
    return near;
  }
} // namespace Geoxy::Met
