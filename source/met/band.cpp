/***************************************************************************************************
 *  @file      source/met/band.cpp
 *  @brief     Two curves near tangent met in closed form, a straight curve and an arc, two arcs or
 *             two straight curves: the touch where they come nearest, or the crossings where they
 *             cross deeper.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "elliptic.hpp"
#include "met/ellipse.hpp"
#include "met/linear.hpp"
#include "met/overlap.hpp"
#include "met/parts.hpp"
#include "met/round.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "straight.hpp"
#include "tangency.hpp"
#include "thorough.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
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
     *  An arc's place where it points a way from its center: the first of its passes there, nothing
     *  where its sweep misses it.
     *
     *  @param  arc    The arc.
     *  @param  toward The way.
     *  @return        The place.
     */
    auto Facing(const Arc &arc, Vector toward) -> std::optional<double>
    {
      return Earliest(Elliptic(arc), Heading(toward).InRadians());
    }

    /**
     *  The touch of a band in closed form: each curve's point at its place where the two come
     *  nearest, and the point halfway between where they come nearest, the double nearest it as
     *  Nearing works it out from the curves' own numbers: the same point to the bit whichever curve
     *  is first and however the plane is turned, as an arc's point at its place, worked out from
     *  its angle, is not.
     *
     *  @param  first   The first curve, and its place.
     *  @param  second  The second, and its place.
     *  @param  halfway The point halfway.
     *  @return         The touch.
     */
    auto Banding(std::pair<const Curve *, double> first, std::pair<const Curve *, double> second,
                 Point halfway) -> Hit
    {
      return {.At     = halfway,
              .First  = {.At = At(*first.first, first.second), .Along = first.second},
              .Second = {.At = At(*second.first, second.second), .Along = second.second},
              .Kind   = Contact::Touch};
    }

    /**
     *  How two curves near tangent meet where they come nearest, in closed form, given the gap
     *  between them there, less than zero where they cross: nothing farther apart than the
     *  tolerance; one touch where they stay within it, apart, or cross no deeper than twice it, as
     *  two hits the curves stay within twice the tolerance between are one; two crossings, the
     *  exact ones at the places given, where they cross deeper. Within an eighth of either bound,
     *  where the thorough way's sampled band may end either side of it, it decides.
     *
     *  @param  apart     The gap where they come nearest, in tolerances.
     *  @param  touch     The touch, made where it is wanted.
     *  @param  crossings Where they cross, roughly, for crossing deeper.
     *  @return           What they meet in; nothing where the thorough way must decide.
     */
    template <typename TouchT>
    auto Gapped(double apart, const TouchT &touch, const Carried &crossings)
        -> std::optional<Nearby>
    {
      constexpr double Below = 7.0 / 8.0;
      constexpr double Above = 9.0 / 8.0;
      const double     deep  = -apart / 2.0;
      Nearby           near;
      if (apart > Above)
      {
        return near;
      }
      if (apart >= Below || (deep >= Below && deep <= Above))
      {
        return std::nullopt;
      }
      if (deep < Below)
      {
        return touch().transform(
            [&near](const Hit &hit) -> Nearby
            {
              near.Ready.front() = hit;
              near.Made          = 1;
              return near;
            });
      }
      for (std::size_t index = 0; index < crossings.Count; ++index)
      {
        near.Add({.Along = {crossings.Along.at(index), 0.0}, .At = crossings.At.at(index)});
      }
      return near;
    }

    /**
     *  A straight curve and an arc whose carriers come within four tolerances of tangent, or cross
     *  however deep, met in closed form where every end lies farther than four tolerances from the
     *  other's carrier, and so outside where they come near: at the foot of the center on the
     *  straight carrier and the circle's point toward it, or nothing where either lies off its
     *  curve; as Gapped decides by the gap between them.
     *
     *  @param  straight  The straight curve.
     *  @param  arc       The arc.
     *  @param  tolerance The tolerance, no finer than their grain.
     *  @return           What they meet in, the straight curve first; nothing where the thorough
     *                    way must decide.
     */
    auto Bend(const Straight &straight, const Arc &arc, const Tolerance &tolerance)
        -> std::optional<Nearby>
    {
      const Run                 &run    = straight.Run;
      const double               clear  = 4.0 * tolerance.Length;
      const double               square = Dot(run.Way, run.Way);
      const std::optional<Round> round  = Rounded(arc);
      if (!(square > 0.0) || !round || !(round->Radius > 16.0 * clear))
      {
        return std::nullopt;
      }
      const double middle  = Dot(round->Center - run.Start, run.Way) / square;
      const Point  foot    = run.Start + (run.Way * middle);
      const double off     = Distance(foot, round->Center);
      const double gap     = off - round->Radius;
      const bool   shallow = gap >= -clear;
      const bool   far     = std::ranges::all_of(std::span{straight.Ends}.first(straight.Count),
                                                 [&round, clear](Point end) -> bool
                                                 { return Off(end, *round) > clear; })
                          && std::ranges::all_of(
                                 shallow ? Cuts(*round) : std::span<const Point>{round->Cuts},
                                 [&run, clear](Point cut) -> bool { return Away(cut, run) > clear; });
      if (!(gap <= clear) || !far)
      {
        return std::nullopt;
      }
      /**
       *  Near tangent, where they come near lies on both curves or off one; crossing deeper, each
       *  crossing lies on both or off one, the ends far from either.
       */
      const bool onBoth = run.Along.Start < middle && middle < run.Along.End && Swept(*round, foot);
      if (shallow && !onBoth)
      {
        return Nearby{};
      }
      Carried crossings;
      if (gap < 0.0)
      {
        const double half =
            std::sqrt((round->Radius * round->Radius) - (off * off)) / std::sqrt(square);
        for (const double along : {middle - half, middle + half})
        {
          const Point point = run.Start + (run.Way * along);
          if (shallow || (run.Along.Start < along && along < run.Along.End && Swept(*round, point)))
          {
            crossings.At.at(crossings.Count)    = point;
            crossings.Along.at(crossings.Count) = along;
            ++crossings.Count;
          }
        }
      }
      const Curve curve = arc;
      return Gapped(
          gap / tolerance.Length,
          [&] -> std::optional<Hit>
          {
            return Facing(arc, foot - round->Center)
                .transform(
                    [&](double along) -> Hit
                    {
                      return Banding({straight.Shape, middle}, {&curve, along},
                                     Crossings::Nearing(*straight.Shape, arc));
                    });
          },
          crossings);
    }

    /**
     *  Which way along the line of centers, from the first center to the second, each of two
     *  circles' points where they come nearest lies from its own center: toward each other where
     *  they kiss from outside; inside, both the far side of the larger one's.
     *
     *  @param  kissing Whether they come near from outside.
     *  @param  larger  Whether the first is the larger.
     *  @return         1 or -1 for each.
     */
    auto Facings(bool kissing, bool larger) -> std::array<double, 2>
    {
      if (kissing)
      {
        return {1.0, -1.0};
      }
      return larger ? std::array{1.0, 1.0} : std::array{-1.0, -1.0};
    }

    /**
     *  Two arcs whose circles come within four tolerances of tangent, from outside or inside, or
     *  cross however deep, met in closed form where every end lies farther than four tolerances
     *  from the other's circle: at each circle's point on the line of centers where they come
     *  nearest, or nothing where either lies off its arc; as Gapped decides by the gap between
     *  them.
     *
     *  @param  arcs      The two arcs.
     *  @param  tolerance The tolerance, no finer than their grain.
     *  @return           What they meet in; nothing where the thorough way must decide.
     */
    auto Rim(const std::array<const Arc *, 2> &arcs, const Tolerance &tolerance)
        -> std::optional<Nearby>
    {
      const double               clear = 4.0 * tolerance.Length;
      const std::optional<Round> one   = Rounded(*arcs[0]);
      const std::optional<Round> other = one ? Rounded(*arcs[1]) : std::nullopt;
      if (!one || !other || !(one->Radius > 16.0 * clear) || !(other->Radius > 16.0 * clear))
      {
        return std::nullopt;
      }
      const double apart   = Distance(one->Center, other->Center);
      const double outside = apart - (one->Radius + other->Radius);
      const double inside  = std::abs(one->Radius - other->Radius) - apart;
      const bool   kissing = std::abs(outside) < std::abs(inside);
      const double gap     = kissing ? outside : inside;
      const bool   shallow = gap >= -clear;
      const auto   ends    = [shallow](const Round &round) -> std::span<const Point>
      {
        return shallow ? Cuts(round) : std::span<const Point>{round.Cuts};
      };
      const bool far = std::ranges::all_of(ends(*one), [&other, clear](Point cut) -> bool
                                           { return Off(cut, *other) > clear; })
                    && std::ranges::all_of(ends(*other), [&one, clear](Point cut) -> bool
                                           { return Off(cut, *one) > clear; });
      if (!(apart > clear) || !(gap <= clear) || !far)
      {
        return std::nullopt;
      }
      const Vector                toward = (other->Center - one->Center) / apart;
      const std::array<double, 2> ways   = Facings(kissing, one->Radius >= other->Radius);
      const Vector                first  = toward * ways[0];
      const Vector                second = toward * ways[1];
      if (shallow
          && !(Swept(*one, one->Center + (first * one->Radius))
               && Swept(*other, other->Center + (second * other->Radius))))
      {
        return Nearby{};
      }
      Carried crossings;
      if (gap < 0.0)
      {
        const double inner  = one->Radius;
        const double outer  = other->Radius;
        const double along  = ((apart * apart) + (inner * inner) - (outer * outer)) / (2.0 * apart);
        const double height = std::sqrt(std::max(0.0, (inner * inner) - (along * along)));
        const Point  base   = one->Center + (toward * along);
        const Vector rise   = Perpendicular(toward) * height;
        for (const Point point : {base + rise, base - rise})
        {
          if (shallow || (Swept(*one, point) && Swept(*other, point)))
          {
            crossings.At.at(crossings.Count) = point;
            ++crossings.Count;
          }
        }
      }
      const Curve curve = *arcs[0];
      const Curve bent  = *arcs[1];
      return Gapped(
          gap / tolerance.Length,
          [&] -> std::optional<Hit>
          {
            return Facing(*arcs[0], first)
                .and_then(
                    [&](double place) -> std::optional<Hit>
                    {
                      return Facing(*arcs[1], second)
                          .transform(
                              [&](double across) -> Hit
                              {
                                return Banding({&curve, place}, {&bent, across},
                                               Crossings::Nearing(arcs, ways));
                              });
                    });
          },
          crossings);
    }

    /**
     *  Whether an end of a straight curve near the other's carrier stops short of it: the carriers'
     *  crossing nearest the end past it, off the straight curve, so that the curve comes nearest
     *  the circle at the end and nowhere else near. A crossing on the curve, however near the end,
     *  is met exactly before this is asked.
     *
     *  @param  straight  The straight curve, its run and its ends.
     *  @param  index     Which end: 0 its start, 1 its end.
     *  @param  crossings Where the carriers cross, at a steep angle.
     *  @return           Whether it does.
     */
    auto Stops(const Straight &straight, std::size_t index, const Carried &crossings) -> bool
    {
      const Point  end    = straight.Ends.at(index);
      const auto   nearer = Distance(end, crossings.At[0]) <= Distance(end, crossings.At[1])
                              ? std::size_t{0}
                              : std::size_t{1};
      const double along  = crossings.Along.at(nearer);
      return index == 0 ? along < straight.Run.Along.Start : along > straight.Run.Along.End;
    }

    /**
     *  A straight curve and an arc crossing at a steep angle, an end of the straight one stopping
     *  short of the arc, met in closed form: every other end, and a whole circle's seam, farther
     *  than four tolerances from the other's carrier; the end short within the tolerance of the
     *  circle a hit at an end, halfway between it, at its own place, and the arc's point nearest
     *  it, where the arc reaches there; the crossings on both, exactly. The thorough way finds that
     *  hit from its candidates, a few steps of a double either way.
     *
     *  @param  straight  The straight curve.
     *  @param  arc       The arc.
     *  @param  tolerance The tolerance, no finer than their grain.
     *  @return           What they meet in, the straight curve first; nothing where it is not so.
     */
    auto Shy(const Straight &straight, const Arc &arc, const Tolerance &tolerance)
        -> std::optional<Nearby>
    {
      const Run                 &run    = straight.Run;
      const double               clear  = 4.0 * tolerance.Length;
      const double               square = Dot(run.Way, run.Way);
      const std::optional<Round> round  = Rounded(arc);
      if (!(square > 0.0) || !round || !(round->Radius > 16.0 * clear))
      {
        return std::nullopt;
      }
      const double middle = Dot(round->Center - run.Start, run.Way) / square;
      const double off    = Distance(run.Start + (run.Way * middle), round->Center);
      const double ratio  = off / round->Radius;
      const bool   steep  = off < round->Radius - clear && (1.0 - (ratio * ratio)) > Steep * Steep;
      const bool   cut    = std::ranges::all_of(round->Cuts, [&run, clear](Point end) -> bool
                                                { return Away(end, run) > clear; });
      if (!steep || !cut)
      {
        return std::nullopt;
      }
      const double half =
          std::sqrt((round->Radius * round->Radius) - (off * off)) / std::sqrt(square);
      const Carried crossings = {
          .At = {run.Start + (run.Way * (middle - half)), run.Start + (run.Way * (middle + half))},
          .Along = {middle - half, middle + half},
          .Count = 2};
      Nearby near;
      for (std::size_t index = 0; index < straight.Count; ++index)
      {
        const Point end = straight.Ends.at(index);
        if (Off(end, *round) > clear)
        {
          continue;
        }
        if (!Stops(straight, index, crossings))
        {
          return std::nullopt;
        }
        const double place = straight.Places.at(index);
        const Curve  curve = arc;
        if (const std::optional<double> along = Facing(arc, end - round->Center))
        {
          const Point here  = At(*straight.Shape, place);
          const Point there = At(curve, *along);
          if (Distance(here, there) <= tolerance.Length)
          {
            near.Ready.at(near.Made) = {.At     = Midpoint(here, there),
                                        .First  = {.At = here, .Along = place},
                                        .Second = {.At = there, .Along = *along},
                                        .Kind   = Contact::End};
            ++near.Made;
          }
        }
      }
      for (std::size_t index = 0; index < 2; ++index)
      {
        const double along = crossings.Along.at(index);
        const Point  point = crossings.At.at(index);
        if (run.Along.Start < along && along < run.Along.End && Swept(*round, point))
        {
          near.Add({.Along = {along, 0.0}, .At = point});
        }
      }
      return near;
    }

    /**
     *  Two arcs of one circle, or nearly, held as points, neither whole, the ends they share, and
     *  the stretches their sweeps share: what Abut and Cocircle each ask, worked out once.
     */
    struct Concord
    {
      std::array<const Arc *, 2> Arcs  = {}; /**< The two arcs.                           */
      Abutted                    Pairs = {}; /**< The ends they share, as Abutting tells. */
      Found                      Lain  = {}; /**< What their sweeps share, as Overlapped
                                                  finds it.                              */
    };

    /**
     *  Two arcs as Abut and Cocircle take them: each plain and not whole, every pair of ends, one
     *  of each, within half the tolerance of each other or farther than four tolerances apart.
     *
     *  @param  arcs      The two arcs.
     *  @param  tolerance The tolerance.
     *  @return           What both ask of them; nothing where they are not so.
     */
    auto Accorded(const std::array<const Arc *, 2> &arcs, const Tolerance &tolerance)
        -> std::optional<Concord>
    {
      const std::optional<Round> one   = Rounded(*arcs[0]);
      const std::optional<Round> other = one ? Rounded(*arcs[1]) : std::nullopt;
      if (!one || !other || Cuts(*one).empty() || Cuts(*other).empty())
      {
        return std::nullopt;
      }
      const std::optional<Abutted> abutted = Abutting({&*one, &*other}, tolerance);
      if (!abutted)
      {
        return std::nullopt;
      }
      return Concord{.Arcs  = arcs,
                     .Pairs = *abutted,
                     .Lain  = Overlapped(Elliptic(*arcs[0]), Elliptic(*arcs[1]), 0.0)};
    }

    /**
     *  Two arcs of one circle, its center and radius the same to the bit, neither whole, that meet
     *  only at ends they share, every other pair of ends farther than four tolerances apart, in
     *  closed form: a hit at each shared end, each arc's end at its own place, halfway between. The
     *  thorough way finds the place a sweep's angles give, a few steps of a double short of the
     *  end.
     *
     *  @param  concord The two arcs, as Accorded takes them.
     *  @return         The hits; nothing where the arcs are not so.
     */
    auto Abut(const Concord &concord) -> std::optional<Nearby>
    {
      const Arc &first  = *concord.Arcs[0];
      const Arc &second = *concord.Arcs[1];
      if (!(first.Center == second.Center && Exactly(first.Radius, second.Radius))
          || concord.Pairs.Count == 0
          || !concord.Lain.Shared.empty())
      {
        return std::nullopt;
      }
      Nearby near;
      for (const std::array<std::size_t, 2> &pair :
           std::span{concord.Pairs.Pairs}.first(concord.Pairs.Count))
      {
        const auto place = [](std::size_t end) -> double
        {
          return static_cast<double>(end);
        };
        const Point here         = At(first, place(pair[0]));
        const Point there        = At(second, place(pair[1]));
        near.Ready.at(near.Made) = {.At     = Midpoint(here, there),
                                    .First  = {.At = here, .Along = place(pair[0])},
                                    .Second = {.At = there, .Along = place(pair[1])},
                                    .Kind   = Contact::End};
        ++near.Made;
      }
      return near;
    }

    /**
     *  Two straight curves meeting only at ends they share, near parallel, end to end along one
     *  line or turning at a corner: each pair of ends, one of each, within half the tolerance of
     *  each other, and every other end farther than four tolerances from the other curve; a hit at
     *  each shared end, each curve's end at its own place, halfway between, in closed form. Two
     *  straight curves meet at one place at most, unless they share a stretch, which would bring an
     *  end near the other.
     *
     *  @param  straights The two straight curves, each of some length.
     *  @param  tolerance The tolerance.
     *  @return           The hits; nothing where they are not so.
     */
    auto Butt(const std::array<const Straight *, 2> &straights, const Tolerance &tolerance)
        -> std::optional<Nearby>
    {
      const Straight &first  = *straights[0];
      const Straight &second = *straights[1];
      const double    clear  = 4.0 * tolerance.Length;
      if (!(Dot(first.Run.Way, first.Run.Way) > 0.0)
          || !(Dot(second.Run.Way, second.Run.Way) > 0.0))
      {
        return std::nullopt;
      }
      const auto lone = [clear](const Straight &ends, std::size_t index,
                                const Straight &against) -> bool
      {
        return Away(ends.Ends.at(index), against.Run) > clear;
      };
      Nearby              near;
      std::array<bool, 2> paired = {};
      for (std::size_t one = 0; one < first.Count; ++one)
      {
        bool shared = false;
        for (std::size_t two = 0; two < second.Count; ++two)
        {
          if (Distance(first.Ends.at(one), second.Ends.at(two)) <= tolerance.Length / 2)
          {
            /**
             *  A second pair of ends that meet: lines along each other end to end, or no longer
             *  than the tolerance, where as many as four pairs do. The thorough way's.
             */
            if (near.Made > 0)
            {
              return std::nullopt;
            }
            const Point here         = At(*first.Shape, first.Places.at(one));
            const Point there        = At(*second.Shape, second.Places.at(two));
            near.Ready.at(near.Made) = {.At     = Midpoint(here, there),
                                        .First  = {.At = here, .Along = first.Places.at(one)},
                                        .Second = {.At = there, .Along = second.Places.at(two)},
                                        .Kind   = Contact::End};
            ++near.Made;
            shared         = true;
            paired.at(two) = true;
          }
        }
        if (!shared && !lone(first, one, second))
        {
          return std::nullopt;
        }
      }
      for (std::size_t two = 0; two < second.Count; ++two)
      {
        if (!paired.at(two) && !lone(second, two, first))
        {
          return std::nullopt;
        }
      }
      if (near.Made != 1)
      {
        return std::nullopt;
      }
      return near;
    }

    /**
     *  Two arcs on one circle, to within a quarter of the tolerance, neither whole, every pair of
     *  ends, one of each, within half the tolerance of each other or farther than four tolerances
     *  apart, met in closed form: the stretches their sweeps share, joined; a hit at each pair of
     *  ends that meet, each arc's end at its own place, halfway between, but where it bounds a
     *  stretch they share. The thorough way finds the same stretches and hits, where the circles
     *  part by a hair, by its tolerance's steps.
     *
     *  @pre              The circles lie within a quarter of the tolerance of each other.
     *  @param  concord   The two arcs, as Accorded takes them.
     *  @param  tolerance The tolerance.
     *  @param  into      Where what they meet in goes, empty on entry.
     *  @return           Whether it was decided so.
     */
    auto Cocircle(const Concord &concord, const Tolerance &tolerance, Intersections &into) -> bool
    {
      const Arc   &first  = *concord.Arcs[0];
      const Arc   &second = *concord.Arcs[1];
      const Found &found  = concord.Lain;
      /**
       *  A stretch no longer than half the tolerance is where ends meet, the sweeps a rounding
       *  apart, not a stretch shared; one longer than twice it is shared; between, the thorough way
       *  decides.
       */
      const double  reach = std::abs(first.Sweep.InRadians()) * first.Radius;
      List<Overlap> shared;
      const auto    span = [reach](const Overlap &overlap) -> double
      {
        return (overlap.First.End - overlap.First.Start) * reach;
      };
      const bool   same    = first.Center == second.Center && Exactly(first.Radius, second.Radius);
      const double least   = same ? 0.0 : tolerance.Length;
      const bool   unclear = std::ranges::any_of(found.Shared,
                                                 [&span, &tolerance](const Overlap &overlap) -> bool
                                                 {
                                                 return span(overlap) > tolerance.Length / 2
                                                     && span(overlap) <= 2.0 * tolerance.Length;
                                                 });
      if (unclear && !same)
      {
        return false;
      }
      for (const Overlap &overlap : found.Shared)
      {
        if (span(overlap) > least)
        {
          shared.push_back(
              {.First  = {.Start = 0.0 + overlap.First.Start, .End = 0.0 + overlap.First.End},
               .Second = {.Start = 0.0 + overlap.Second.Start, .End = 0.0 + overlap.Second.End},
               .Along  = overlap.Along});
        }
      }
      const List<Overlap> joined = Joined(std::move(shared));
      for (const std::array<std::size_t, 2> &pair :
           std::span{concord.Pairs.Pairs}.first(concord.Pairs.Count))
      {
        const auto place = [](std::size_t end) -> double
        {
          return static_cast<double>(end);
        };
        const Point here  = At(first, place(pair[0]));
        const Point there = At(second, place(pair[1]));
        const bool  bound = std::ranges::any_of(
            joined,
            [&first, here, &tolerance](const Overlap &overlap) -> bool
            {
              return Distance(here, At(first, overlap.First.Start)) <= tolerance.Length
                  || Distance(here, At(first, overlap.First.End))   <= tolerance.Length;
            });
        if (!bound)
        {
          into.Hits.push_back({.At     = Midpoint(here, there),
                               .First  = {.At = here, .Along = place(pair[0])},
                               .Second = {.At = there, .Along = place(pair[1])},
                               .Kind   = Contact::End});
        }
      }
      std::ranges::sort(into.Hits, {}, [](const Hit &hit) -> double { return hit.First.Along; });
      into.Overlaps.assign(joined.begin(), joined.end());
      return true;
    }

    /**
     *  Two arcs of different circles meeting only at ends they share, a vertex: a pair of ends, one
     *  of each, within half the tolerance of each other, every other pair farther than four
     *  tolerances apart, every other end farther than that from the other's circle; and where the
     *  circles cross or come near elsewhere, there off one arc or the other. A hit at the shared
     *  ends, each arc's end at its own place, halfway between, in closed form, whether the circles
     *  cross there, shallow or steep, or touch. The thorough way polishes it to their crossing, a
     *  few steps of a double from the ends.
     *
     *  @param  arcs      The two arcs.
     *  @param  tolerance The tolerance.
     *  @return           The hit; nothing where they are not so.
     */
    auto Joint(const std::array<const Arc *, 2> &arcs, const Tolerance &tolerance)
        -> std::optional<Nearby>
    {
      const Arc                 &first  = *arcs[0];
      const Arc                 &second = *arcs[1];
      const double               clear  = 4.0 * tolerance.Length;
      const std::optional<Round> one    = Rounded(first);
      const std::optional<Round> other  = one ? Rounded(second) : std::nullopt;
      if (!one || !other || Cuts(*one).empty() || Cuts(*other).empty())
      {
        return std::nullopt;
      }
      const double                 apart   = Distance(one->Center, other->Center);
      const std::optional<Abutted> abutted = Abutting({&*one, &*other}, tolerance);
      if (!(apart > clear) || !abutted || abutted->Count != 1)
      {
        return std::nullopt;
      }
      const std::array<std::size_t, 2> pair   = abutted->Pairs.front();
      const Point                      vertex = one->Cuts.at(pair[0]);
      const bool                       far    = Off(one->Cuts.at(1 - pair[0]), *other) > clear
                                             && Off(other->Cuts.at(1 - pair[1]), *one) > clear;
      if (!far)
      {
        return std::nullopt;
      }
      /**
       *  Where the circles cross, or come nearest: each such place near the vertex, or off an arc.
       */
      const Vector toward = (other->Center - one->Center) / apart;
      const double inner  = one->Radius;
      const double outer  = other->Radius;
      const double along  = ((apart * apart) + (inner * inner) - (outer * outer)) / (2.0 * apart);
      const double height = std::sqrt(std::max(0.0, (inner * inner) - (along * along)));
      const Point  base   = one->Center + (toward * along);
      const Vector rise   = Perpendicular(toward) * height;
      for (const Point point : {base + rise, base - rise})
      {
        const bool there = Distance(point, vertex) <= clear;
        const bool off   = !Swept(*one, point) || !Swept(*other, point);
        if (!there && !off)
        {
          return std::nullopt;
        }
      }
      const Curve curve = first;
      const Curve bent  = second;
      const auto  place = [](std::size_t end) -> double
      {
        return static_cast<double>(end);
      };
      const Point here  = At(curve, place(pair[0]));
      const Point there = At(bent, place(pair[1]));
      Nearby      near;
      near.Ready.front() = {.At     = Midpoint(here, there),
                            .First  = {.At = here, .Along = place(pair[0])},
                            .Second = {.At = there, .Along = place(pair[1])},
                            .Kind   = Contact::End};
      near.Made          = 1;
      return near;
    }

    /**
     *  Two arcs on one circle, or nearly, as Accorded takes them, met in closed form: meeting only
     *  at ends they share, else as Cocircle meets them.
     *
     *  @param  concord   The two arcs, as Accorded takes them.
     *  @param  tolerance The tolerance.
     *  @param  into      Where what they meet in goes, empty on entry.
     *  @return           Whether it was decided so.
     */
    auto Cocircling(const Concord &concord, const Tolerance &tolerance, Intersections &into) -> bool
    {
      if (const std::optional<Nearby> near = Abut(concord))
      {
        for (const Hit &hit : std::span{near->Ready}.first(near->Made))
        {
          into.Hits.push_back(hit);
        }
        std::ranges::sort(into.Hits, {}, [](const Hit &hit) -> double { return hit.First.Along; });
        return true;
      }
      return Cocircle(concord, tolerance, into);
    }

    /**
     *  Two arcs on one circle, or nearly, met in closed form: meeting only at ends they share, else
     *  as Cocircle meets them; nothing for arcs of circles farther apart than a quarter of the
     *  tolerance, centers and radii together, which neither meets.
     *
     *  @param  arcs      The two arcs.
     *  @param  tolerance The tolerance.
     *  @param  into      Where what they meet in goes, empty on entry.
     *  @return           Whether it was decided so.
     */
    auto Cocircled(const std::array<const Arc *, 2> &arcs, const Tolerance &tolerance,
                   Intersections &into) -> bool
    {
      const Arc   &first  = *arcs[0];
      const Arc   &second = *arcs[1];
      const double drift =
          Distance(first.Center, second.Center) + std::abs(first.Radius - second.Radius);
      const bool same = first.Center == second.Center && Exactly(first.Radius, second.Radius);
      if (!same && !(drift <= tolerance.Length / 4))
      {
        return false;
      }
      const std::optional<Concord> concord = Accorded(arcs, tolerance);
      return concord && Cocircling(*concord, tolerance, into);
    }

  } // namespace

  auto Coincident(const std::array<const Arc *, 2> &arcs, const Tolerance &tolerance,
                  Intersections &into) -> bool
  {
    const Arc   &first  = *arcs[0];
    const Arc   &second = *arcs[1];
    const double drift =
        Distance(first.Center, second.Center) + std::abs(first.Radius - second.Radius);
    /**
     *  Arcs of one circle to the bit are Concentrically's first. Where it leaves them, or the
     *  circles part by no more than a quarter of the tolerance, centers and radii together, the
     *  direct way's rough crossings find nothing: Kissing asks the centers to lie more than four
     *  tolerances apart. Banded then meets them as Cocircled does.
     */
    if (!(first.Center == second.Center && Exactly(first.Radius, second.Radius)))
    {
      return drift <= tolerance.Length / 4 && Cocircled(arcs, tolerance, into);
    }
    return Concentrically(arcs, tolerance, into) || Cocircled(arcs, tolerance, into);
  }
} // namespace Geoxy::Met

namespace Geoxy
{
  auto Banded(const Curves &curves, const Tolerance &tolerance, Intersections &into) -> bool
  {
    const Curve               &first  = curves[0];
    const Curve               &second = curves[1];
    const auto                 one    = Met::Straightened(first);
    const auto                 other  = Met::Straightened(second);
    const auto                *round  = std::get_if<Arc>(&first);
    const auto                *across = std::get_if<Arc>(&second);
    std::optional<Met::Nearby> near;
    if (one && other)
    {
      near = Met::Butt({&*one, &*other}, tolerance);
    }
    else if (one && across != nullptr)
    {
      near = Met::Bend(*one, *across, tolerance);
      if (!near)
      {
        near = Met::Shy(*one, *across, tolerance);
      }
    }
    else if (round != nullptr && other)
    {
      near = Met::Bend(*other, *round, tolerance);
      if (!near)
      {
        near = Met::Shy(*other, *round, tolerance);
      }
      near = near.transform(
          [](Met::Nearby locals) -> Met::Nearby
          {
            for (Approach &local : std::span{locals.Items}.first(locals.Count))
            {
              std::swap(local.Along[0], local.Along[1]);
            }
            for (Hit &hit : std::span{locals.Ready}.first(locals.Made))
            {
              std::swap(hit.First, hit.Second);
            }
            return locals;
          });
    }
    else if (round != nullptr && across != nullptr)
    {
      if (Met::Cocircled({round, across}, tolerance, into))
      {
        return true;
      }
      near = Met::Rim({round, across}, tolerance);
      if (!near)
      {
        near = Met::Joint({round, across}, tolerance);
      }
    }
    return near
        .transform(
            [&first, &second, &into](const Met::Nearby &places) -> bool
            {
              /**
               *  A band has a touch or crossings, never both: where the crossings are not found,
               *  there is no touch to add to the meeting left empty.
               */
              const bool found = Exactly({&first, &second}, places.Held(), into);
              for (const Hit &hit : std::span{places.Ready}.first(places.Made))
              {
                into.Hits.push_back(hit);
              }
              std::ranges::sort(into.Hits, {},
                                [](const Hit &hit) -> double { return hit.First.Along; });
              return found;
            })
        .value_or(false);
  }
} // namespace Geoxy
