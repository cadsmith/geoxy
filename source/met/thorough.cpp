/***************************************************************************************************
 *  @file      source/met/thorough.cpp
 *  @brief     Two curves met the thorough way: every pair of parts met, the stretches shared
 *             settled, the candidates judged and merged, and each hit put where the curves meet
 *             exactly.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "thorough.hpp"

#include "crossing.hpp"
#include "met/floor.hpp"
#include "met/judge.hpp"
#include "met/overlap.hpp"
#include "met/parts.hpp"
#include "met/scale.hpp"
#include "met/slight.hpp"
#include "met/solve.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "touch.hpp"
#include "trig.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy::Met
{
  namespace
  {
    /**
     *  Whether a crossing is clear: its point told by the curves' numbers alone, no coordinate but
     *  0 under 2^-40 of the size about it, where the last rounding of Newton's method hangs on
     *  where it started; and the curves' tangents there on two sides of each other, further apart
     *  than twice the tolerance's angle, and still so, by at least half as far, at each end of the
     *  band within which curves as straight as their tangents there stay within four tolerances.
     *  Any hit polished within the tolerance of a clear crossing is that crossing, the same way
     *  across; a shallow crossing, or one beside a cusp, about which the curves may cross back or
     *  turn along each other, is not clear. The step along the second is the first's over the
     *  second's length twice, not over its square, which a way too short to square loses.
     *
     *  @param  scene    The curves being met.
     *  @param  crossing The crossing, found exactly.
     *  @return          Whether it is clear.
     */
    auto Clear(const Scene &scene, const Crossings::Crossed &crossing) -> bool
    {
      const Curve &first  = *scene.First;
      const Curve &second = *scene.Second;
      const auto   sine   = [](Vector one, Vector other) -> double
      {
        const double speeds = Length(one) * Length(other);
        return Cross(one, other) / std::max(speeds, std::numeric_limits<double>::min());
      };
      const Vector way   = Derivative(first, crossing.Along[0]);
      const Vector back  = Derivative(second, crossing.Along[1]);
      const Point  point = crossing.At;
      const double size =
          std::max({std::abs(point.X), std::abs(point.Y), Length(way), Length(back)});
      const auto told = [size](double coordinate) -> bool
      {
        return Zero(coordinate) || std::abs(coordinate) >= 0x1p-40 * size;
      };
      const double least  = 2.0 * Trig::Sin(scene.Tolerance.Radians);
      const double across = sine(way, back);
      if (!told(point.X) || !told(point.Y) || !(std::abs(across) > least))
      {
        return false;
      }
      const double reach  = 4.0 * scene.Tolerance.Length / std::abs(across);
      const double ahead  = reach / Length(way);
      const double beside = ahead * (Dot(way, back) / Length(back)) / Length(back);
      const Range  one    = Domain(first);
      const Range  other  = Domain(second);
      return std::ranges::all_of(
          std::array{-1.0, 1.0},
          [&](double sign) -> bool
          {
            const double turn =
                sine(Derivative(first,
                                std::clamp(crossing.Along[0] + (sign * ahead), one.Start, one.End)),
                     Derivative(second, std::clamp(crossing.Along[1] + (sign * beside), other.Start,
                                                   other.End)));
            return std::copysign(1.0, across) * turn > std::max(0.5 * std::abs(across), least);
          });
    }

    /**
     *  Whether a touch is smooth: neither curve stalls at it, its speed there under 2^-24 of the
     *  faster of its speeds a 2^-8 part of its span to either side, as at a cusp, where where the
     *  curves come nearest hangs on where Newton's method started.
     *
     *  @param  scene The curves being met.
     *  @param  touch The touch, where the curves come nearest.
     *  @return       Whether it is smooth.
     */
    auto Smooth(const Scene &scene, const Crossings::Crossed &touch) -> bool
    {
      const std::array<const Curve *, 2> curves = {scene.First, scene.Second};
      return std::ranges::all_of(
          std::array{0UZ, 1UZ},
          [&curves, &touch](std::size_t side) -> bool
          {
            const Curve &curve  = *curves.at(side);
            const Range  domain = Domain(curve);
            const double span   = domain.End - domain.Start;
            const double step   = 0x1p-8 * (std::isfinite(span) ? span : 1.0);
            const double along  = touch.Along.at(side);
            const auto   speed  = [&curve, &domain](double place) -> double
            {
              return Length(Derivative(curve, std::clamp(place, domain.Start, domain.End)));
            };
            return speed(along) > 0x1p-24 * std::max(speed(along - step), speed(along + step));
          });
    }

    /**
     *  Whether a point lies within the tolerance of an end of either curve.
     *
     *  @param  scene The curves being met.
     *  @param  point The point.
     *  @return       Whether it does.
     */
    auto Endward(const Scene &scene, Point point) -> bool
    {
      const auto near = [&scene, point](Point end) -> bool
      {
        return Distance(point, end) <= scene.Tolerance.Length;
      };
      return std::ranges::any_of(scene.FirstEnds, near)
          || std::ranges::any_of(scene.SecondEnds, near);
    }

    /**
     *  Where two curves cross near an end hit that no exact way places, as on a polyline's bulge:
     *  polished in double-double as a crossing is, and taken only where it is clear, as a crossing
     *  is settled, and its point lies on both curves as Located places it. A crossing of the
     *  curves' carriers past an end, or a shallow one whose band runs far along both, as where an
     *  offset's pieces meet end to end near tangent, leaves the hit where its candidates found it.
     *  A place of nothing is +0, as Newton's method may leave it -0.
     *
     *  @param  scene  The curves being met.
     *  @param  places The hit's place on each.
     *  @param  near   The hit's point.
     *  @return        The crossing; nothing where it is not found so.
     */
    auto Polishing(const Scene &scene, const std::array<double, 2> &places, Point near)
        -> std::optional<Crossings::Crossed>
    {
      const std::array<const Curve *, 2> curves = {scene.First, scene.Second};
      std::optional<Crossings::Crossed>  found  = Crossings::Polished(curves, places, near);
      if (!found || !Clear(scene, *found))
      {
        return std::nullopt;
      }
      for (std::size_t side = 0; side < 2; ++side)
      {
        if (!Crossings::Located(*curves.at(side), found->Along.at(side), found->At).has_value())
        {
          return std::nullopt;
        }
      }
      found->Along = {found->Along[0] + 0.0, found->Along[1] + 0.0};
      return found;
    }

    /**
     *  A hit put where the curves meet exactly: a touch where they come nearest, polished to the
     *  bit, whatever found it; a crossing or an end, or a touch that does not settle so, where they
     *  cross exactly, a crossing or such a touch by Newton's method in double-double where no
     *  closed form has it.
     *
     *  A touch's two points may lie as far apart as twice the tolerance: two crossings are one
     *  touch where the curves stay within twice it between them, and where they cross deepest, as
     *  far apart as they come there, is where that touch is. Held to the tolerance, a crossing
     *  deeper than it was put at the crossing its candidate lay nearest, which of two hung on the
     *  order of the curves and the last bits of their numbers.
     *
     *  An end is put where the curves cross only where that lies within the tolerance of an end:
     *  the band about a shallow crossing makes an end and a crossing far along it one hit, which is
     *  an end of one of them, and stood at the crossing, 1e8 tolerances from every end; else it
     *  stays as its candidates found it. Where no exact way has that crossing, as for a polyline's
     *  bulge, it is the one Polishing finds, and the hit is not settled by it.
     *
     *  @param  scene The curves being met.
     *  @param  hit   The hit, as its candidates found it; put where the curves meet.
     *  @return       Whether no other candidate could have made it otherwise: a smooth touch put
     *                where the curves come nearest, an end where they cross exactly, a crossing
     *                found exactly and clear.
     */
    auto Exacted(const Scene &scene, Hit &hit) -> bool
    {
      const std::array<const Curve *, 2> curves = {scene.First, scene.Second};
      const std::array<double, 2>        places = {hit.First.Along, hit.Second.Along};
      std::optional<Crossings::Crossed>  crossing =
          hit.Kind == Contact::Touch
              ? Touches::Nearest(curves, places, 2.0 * scene.Tolerance.Length)
              : std::nullopt;
      bool settled = crossing.has_value() && Smooth(scene, *crossing);
      if (!crossing)
      {
        crossing = Crossings::Crossing(curves, places, hit.At);
        settled  = crossing.has_value() && hit.Kind == Contact::End;
        if (!crossing && hit.Kind == Contact::End)
        {
          crossing = Polishing(scene, places, hit.At);
        }
        if (crossing && hit.Kind == Contact::End && !Endward(scene, crossing->At))
        {
          crossing.reset();
          settled = false;
        }
      }
      if (hit.Kind != Contact::End)
      {
        crossing = crossing.or_else([&curves, &places, &hit] -> std::optional<Crossings::Crossed>
                                    { return Crossings::Polished(curves, places, hit.At); });
        settled  = crossing.has_value() && Clear(scene, *crossing);
      }
      if (crossing)
      {
        hit.At     = crossing->At;
        hit.First  = {.At = crossing->On[0], .Along = crossing->Along[0]};
        hit.Second = {.At = crossing->On[1], .Along = crossing->Along[1]};
      }
      return settled;
    }

    /**
     *  A hit as its candidates found it, and as it was put where the curves meet exactly.
     */
    struct Put
    {
      Hit  Found   = {};    /**< The hit as found.          */
      Hit  Exact   = {};    /**< The hit put exactly.       */
      bool Settled = false; /**< Whether it was settled so. */
    };

    /**
     *  What meeting two curves finds before its candidates are judged: what every pair of parts
     *  found, the stretches the curves share and the bands about touches; and which ends of parts
     *  are met with.
     */
    struct Groundwork
    {
      List<Local>   Gathered = {};    /**< The candidates every pair of parts found.      */
      List<Local>   Middles  = {};    /**< The middle of each band about a touch.         */
      List<Overlap> Shared   = {};    /**< The stretches shared, joined.                  */
      List<Overlap> Bands    = {};    /**< The bands about touches and shallow crossings. */
      bool          Spent    = false; /**< Whether a budget ran out.                      */
      bool          All      = false; /**< Whether every end of a part is met with, or only
                                           those that may lie near the other curve.         */
      List<Local> Weighed = {};       /**< The candidates the way before judged, in order. */
      List<std::optional<Seen>> Verdicts = {}; /**< How it judged each.               */
      List<Put>                 Placed   = {}; /**< The hits it put exactly.          */
    };

    /**
     *  Whether two hits are one to the bit, so that whatever is worked out from one is so from the
     *  other.
     *
     *  @param  one   A hit.
     *  @param  other Another.
     *  @return       Whether they are.
     */
    auto Twin(const Hit &one, const Hit &other) -> bool
    {
      const auto bits = [](const Hit &hit) -> std::array<std::uint64_t, 8>
      {
        return {std::bit_cast<std::uint64_t>(hit.At.X),
                std::bit_cast<std::uint64_t>(hit.At.Y),
                std::bit_cast<std::uint64_t>(hit.First.At.X),
                std::bit_cast<std::uint64_t>(hit.First.At.Y),
                std::bit_cast<std::uint64_t>(hit.First.Along),
                std::bit_cast<std::uint64_t>(hit.Second.At.X),
                std::bit_cast<std::uint64_t>(hit.Second.At.Y),
                std::bit_cast<std::uint64_t>(hit.Second.Along)};
      };
      return one.Kind == other.Kind && bits(one) == bits(other);
    }

    /**
     *  The hits put exactly, each put where another was put, its point and places to the bit,
     *  dropped: two candidates far apart along a curve, as a curve's two ends within the tolerance
     *  of one point of the other are, judged two hits and each put where the curves meet exactly,
     *  meet there in one hit, an end if either is, else a touch if either is, as Merged tells a
     *  run. Of two such, the one found nearer its place along the first curve is kept, so that the
     *  hits stay in order along it; of two found as near, the first.
     *
     *  @param placed The hits as found and as put, in order along the first curve.
     *  @param into   The meeting, its hits the ones put, in that order; each twin dropped.
     */
    void Untwinned(std::span<const Put> placed, Intersections &into)
    {
      const auto moved = [placed](std::size_t index) -> double
      {
        return std::abs(placed[index].Exact.First.Along - placed[index].Found.First.Along);
      };
      List<std::uint8_t> dropped(placed.size(), 0);
      for (std::size_t one = 0; one < placed.size(); ++one)
      {
        for (std::size_t other = one + 1; dropped[one] == 0 && other < placed.size(); ++other)
        {
          Hit twin  = placed[other].Exact;
          twin.Kind = placed[one].Exact.Kind;
          if (dropped[other] != 0 || !Twin(placed[one].Exact, twin))
          {
            continue;
          }
          const bool        later      = !(moved(other) < moved(one));
          const std::size_t kept       = later ? one : other;
          into.Hits[kept].Kind         = std::max(into.Hits[one].Kind, into.Hits[other].Kind);
          dropped[later ? other : one] = 1;
        }
      }
      std::size_t kept = 0;
      for (std::size_t index = 0; index < placed.size(); ++index)
      {
        if (dropped[index] == 0)
        {
          into.Hits[kept++] = into.Hits[index];
        }
      }
      into.Hits.resize(kept);
    }

    /**
     *  Hits put in order along an open first curve, those at one place in the order they came, in
     *  place: a run's hit is its end or the one where the curves come nearest, a touch from a band
     *  is where they come nearest, and a hit put exactly is where they meet, any of which may lie
     *  past the next hit along the first curve. A closed one's run about its seam stays first, as
     *  the direct ways give it.
     *
     *  @param scene The curves being met.
     *  @param hits  The hits, put in order.
     */
    void Placewise(const Scene &scene, std::vector<Hit> &hits)
    {
      if (scene.FirstEnds.empty() && std::isfinite(Domain(*scene.First).Start))
      {
        return;
      }
      const auto place = [](const Hit &hit) -> double
      {
        return hit.First.Along;
      };
      for (auto next = hits.begin(); next != hits.end(); ++next)
      {
        std::rotate(std::ranges::upper_bound(hits.begin(), next, place(*next), {}, place), next,
                    std::next(next));
      }
    }

    /**
     *  What two curves meet in, from what every pair of parts found and the ends of parts: the
     *  candidates judged and merged into hits, those a shared stretch runs on into dropped, and
     *  each put where the curves meet exactly. A candidate or a hit the way before worked out is
     *  taken from it, as working it out again gives the same.
     *
     *  @param  scene  The curves being met.
     *  @param  ground What was found before the candidates were judged, and what the way before
     *                 worked out, none for the first; this way's goes in its place.
     *  @param  into   Where what they meet in goes, whatever it held replaced.
     *  @return        Whether it is what every end met with would find: every end was, or every hit
     *                 is settled without the far ones and no stretch is shared.
     */
    auto Concluded(const Scene &scene, Groundwork &ground, Intersections &into) -> bool
    {
      List<Local> locals  = ground.Gathered;
      bool        skipped = false;
      for (const std::size_t side : {0UZ, 1UZ})
      {
        const Tips tips = Ended(scene, side, ground.All);
        locals.insert(locals.end(), tips.Locals.begin(), tips.Locals.end());
        skipped = skipped || tips.Skipped;
      }
      locals.insert(locals.end(), ground.Middles.begin(), ground.Middles.end());
      /**
       *  The candidates the way before judged come in the same order, the far ends' among them.
       */
      List<std::optional<Seen>> verdicts;
      List<Seen>                seen;
      std::size_t               next = 0;
      for (const Local &local : locals)
      {
        const bool                again = next < ground.Weighed.size()
                                       && Exactly(ground.Weighed[next].First, local.First)
                                       && Exactly(ground.Weighed[next].Second, local.Second);
        const std::optional<Seen> verdict =
            again ? ground.Verdicts[next] : Judged(scene, ground.Shared, local);
        next += again ? 1 : 0;
        if (verdict)
        {
          seen.push_back(*verdict);
        }
        verdicts.push_back(verdict);
      }
      List<Hit> hits = Hits(scene, ground.Shared, seen);
      Touched(scene, ground.Bands, hits);
      into.Hits.assign(hits.begin(), hits.end());
      into.Overlaps.assign(ground.Shared.begin(), ground.Shared.end());
      into.Spent = ground.Spent
                || ((!hits.empty() || !ground.Shared.empty())
                    && (Overwound(*scene.First) || Overwound(*scene.Second)));
      Absorb(scene, into);
      /**
       *  The hits the way before put come in the same order along the first curve: each is looked
       *  for from just past the last one found.
       */
      bool      settled = into.Overlaps.empty();
      List<Put> placed;
      auto      from = ground.Placed.begin();
      for (Hit &hit : into.Hits)
      {
        const auto twin = [&hit](const Put &put) -> bool
        {
          return Twin(put.Found, hit);
        };
        auto before = std::find_if(from, ground.Placed.end(), twin);
        before = before != ground.Placed.end() ? before : std::ranges::find_if(ground.Placed, twin);
        Put put = {.Found = hit, .Exact = hit, .Settled = false};
        if (before != ground.Placed.end())
        {
          put  = *before;
          from = std::next(before);
        }
        else
        {
          put.Settled = Exacted(scene, put.Exact);
        }
        hit     = put.Exact;
        settled = put.Settled && settled;
        placed.push_back(put);
      }
      Untwinned(placed, into);
      Placewise(scene, into.Hits);
      ground.Weighed             = std::move(locals);
      ground.Verdicts            = std::move(verdicts);
      ground.Placed              = std::move(placed);
      const List<Overlap> joined = Joined({into.Overlaps.begin(), into.Overlaps.end()});
      into.Overlaps.assign(joined.begin(), joined.end());
      return settled || !skipped;
    }

    /**
     *  Where two curves meet, found the thorough way in the order given. The ends of parts far from
     *  the other curve are met with only where the hits found without them are not settled: their
     *  candidates, polished from a long way off, take Newton's method to hits the others find, and
     *  change only how a hit about a touch, a cusp or a shallow crossing is told.
     *
     *  @param curves    The two curves, the first first.
     *  @param tolerance The tolerance.
     *  @param into      Where what they meet in goes, whatever it held replaced.
     */
    void Worked(const Curves &curves, const Tolerance &tolerance, Intersections &into)
    {
      const Curve &first  = curves[0];
      const Curve &second = curves[1];
      into.Hits.clear();
      into.Overlaps.clear();
      into.Spent             = false;
      const List<Part> parts = Parted(first);
      if (parts.empty())
      {
        return;
      }
      if (first == second)
      {
        const Range domain = Domain(first);
        into.Overlaps.push_back({.First = domain, .Second = domain, .Along = true});
        return;
      }
      const List<Part> others = Parted(second);
      if (others.empty())
      {
        return;
      }
      const Scene        scene    = {.First        = &first,
                                     .Second       = &second,
                                     .Tolerance    = Grained(first, second, tolerance),
                                     .FirstEnds    = Extremes(first),
                                     .SecondEnds   = Extremes(second),
                                     .Splined      = std::holds_alternative<Spline>(first)
                                                  || std::holds_alternative<Spline>(second),
                                     .FirstFloors  = Floors(parts),
                                     .SecondFloors = Floors(others),
                                     .FirstParts   = &parts,
                                     .SecondParts  = &others};
      Found              gathered = Gathered(parts, others, scene);
      const List<double> bounds   = Boundaries(scene, {&parts, &others});
      List<Overlap>      shared   = std::move(gathered.Shared);
      for (Overlap &overlap : shared)
      {
        overlap = Snapped(scene, bounds, overlap);
      }
      const auto banded =
          std::ranges::partition(shared, [&scene, &bounds](const Overlap &overlap) -> bool
                                 { return Settled(scene, bounds, overlap); });
      Groundwork ground = {.Gathered = std::move(gathered.Locals), .Spent = gathered.Spent};
      for (const Overlap &band : banded)
      {
        ground.Middles.push_back(
            {.First  = band.First.Start + ((band.First.End - band.First.Start) / 2),
             .Second = band.Second.Start + ((band.Second.End - band.Second.Start) / 2)});
      }
      ground.Bands.assign(banded.begin(), banded.end());
      shared.erase(banded.begin(), banded.end());
      ground.Shared = Joined(std::move(shared));
      if (!Concluded(scene, ground, into))
      {
        ground.All = true;
        Concluded(scene, ground, into);
      }
    }

  } // namespace
} // namespace Geoxy::Met

namespace Geoxy
{
  auto Shares(const Curves &curves, const Tolerance &tolerance) -> bool
  {
    return !Met::Matched(Met::Parted(curves[0]).front(), Met::Parted(curves[1]).front(), tolerance)
                .Shared.empty();
  }

  void Thorough(const Curves &curves, const Tolerance &tolerance, Intersections &into)
  {
    /**
     *  A slight arc within the tolerance of its chord is met as the chord, as Intersect meets it.
     */
    const auto worked = [](const void *, const Curves &stood, const Tolerance &held,
                           Intersections &found) -> void
    {
      Met::Worked(stood, held, found);
    };
    if (!Met::Slightly(curves,
                       {.Sane      = tolerance,
                        .Length    = Met::Grained(curves[0], curves[1], tolerance).Length,
                        .Intersect = worked,
                        .Deep      = Met::Worked},
                       into))
    {
      Met::Worked(curves, tolerance, into);
    }
  }
} // namespace Geoxy
