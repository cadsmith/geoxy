/***************************************************************************************************
 *  @file      source/swift/met.cpp
 *  @brief     Two curves met by a swift kernel where one is straight, the lane and the kernel its
 *             other curve asks for; or where both are Bézier curves, the twin kernel.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/met.hpp"

#include "crossing.hpp"
#include "straight.hpp"
#include "sure.hpp"
#include "swift/braid.hpp"
#include "swift/lane.hpp"
#include "swift/pieced.hpp"
#include "swift/ring.hpp"
#include "swift/section.hpp"
#include "swift/twin.hpp"
#include "swift/vault.hpp"

#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <type_traits>
#include <variant>

namespace Geoxy::Swift
{
  namespace
  {
    /**
     *  A lane from a point, its way exact and rounded, its stretch and its ends.
     *
     *  @param  start The point at 0.
     *  @param  way   Its way, exact.
     *  @param  along Its stretch.
     *  @return       The lane, its ends to be added; nothing where the way has no length or a
     *                number lies outside the range the kernels take.
     */
    auto MadeLane(Point start, const std::array<Sure::Number, 2> &way, Range along)
        -> std::optional<Lane>
    {
      const Vector rough = {.X = way[0].Value.High, .Y = way[1].Value.High};
      const double speed = Length(rough);
      if (!Held({start.X, start.Y, rough.X, rough.Y}) || !(speed > Narrowest))
      {
        return std::nullopt;
      }
      return Lane{.Start = start,
                  .Way   = way,
                  .Rough = rough,
                  .Speed = speed,
                  .Along = along,
                  .Ends  = {},
                  .Count = 0,
                  .Shape = nullptr,
                  .Lead  = 0.0};
    }

    /**
     *  A lane and the curve it is met with, which of the two that curve is, and the tolerance.
     */
    struct Pairing
    {
      const Lane      *Straight  = nullptr; /**< The lane.                              */
      const Curve     *Shape     = nullptr; /**< The other curve.                       */
      bool             Second    = false;   /**< Whether that is the second of the two. */
      Geoxy::Tolerance Tolerance = {};      /**< The tolerance.                         */
    };

    /**
     *  The swift kernel a lane and another curve ask for.
     *
     *  @param  pairing The lane and the curve.
     *  @return         What it found; nothing where it is not sure, or there is no kernel for the
     *                  curve.
     */
    auto Dispatched(const Pairing &pairing) -> std::optional<Spots>
    {
      const Lane      &lane      = *pairing.Straight;
      const Tolerance &tolerance = pairing.Tolerance;
      return std::visit(
          [&lane, &tolerance](const auto &shape) -> std::optional<Spots>
          {
            using Shape = std::decay_t<decltype(shape)>;
            if constexpr (std::is_same_v<Shape, Ellipse>)
            {
              return Oval(lane, shape, tolerance);
            }
            else if constexpr (std::is_same_v<Shape, Polyline>)
            {
              return Bulged(lane, shape, tolerance);
            }
            else if constexpr (std::is_same_v<Shape, Bezier::Quadratic>
                               || std::is_same_v<Shape, Bezier::Cubic>)
            {
              return Bernstein(lane, Sectioned(shape), tolerance);
            }
            else if constexpr (std::is_same_v<Shape, Spline>)
            {
              return Strung(lane, shape, tolerance);
            }
            else
            {
              return std::nullopt;
            }
          },
          *pairing.Shape);
    }

    /**
     *  A crossing a kernel found made a hit: placed as the thorough way places it, or at the places
     *  the kernel worked out, each the curve's own.
     *
     *  @param  pairing The lane and the curve.
     *  @param  spot    The crossing.
     *  @param  into    Where the hit goes.
     *  @return         Whether it was placed.
     */
    auto Placed(const Pairing &pairing, const Spot &spot, Intersections &into) -> bool
    {
      const std::array<const Curve *, 2> pair =
          pairing.Second ? std::array<const Curve *, 2>{pairing.Straight->Shape, pairing.Shape}
                         : std::array<const Curve *, 2>{pairing.Shape, pairing.Straight->Shape};
      const std::array<double, 2> places =
          pairing.Second ? spot.Along : std::array<double, 2>{spot.Along[1], spot.Along[0]};
      const std::optional<Crossings::Crossed> crossed =
          spot.Own ? std::optional<Crossings::Crossed>{Crossings::Crossed{
                         .At    = spot.At,
                         .Along = places,
                         .On    = {At(*pair[0], places[0]), At(*pair[1], places[1])}}}
                   : Crossings::Placed(pair, places, spot.At);
      if (crossed)
      {
        into.Hits.push_back({.At     = crossed->At,
                             .First  = {.At = crossed->On[0], .Along = crossed->Along[0]},
                             .Second = {.At = crossed->On[1], .Along = crossed->Along[1]},
                             .Kind   = Contact::Cross});
      }
      return crossed.has_value();
    }

    /**
     *  Whether no two hits are one, as the thorough way merges them: two are one only where they
     *  lie within the tolerance of each other or where, halfway between them along each curve, each
     *  curve is within twice the tolerance of the other. Each pair lies farther apart than the
     *  tolerance, and the other curve's point halfway between them farther than twice it from the
     *  lane, both with the margin: halfway the plain way, and the short way round where the other
     *  curve may close on itself.
     *
     *  @param  pairing The lane and the curve.
     *  @param  meeting The hits.
     *  @return         Whether none are one.
     */
    auto Apart(const Pairing &pairing, const Intersections &meeting) -> bool
    {
      const Lane  &lane   = *pairing.Straight;
      const Curve &shape  = *pairing.Shape;
      const Range  domain = Domain(shape);
      const double whole  = domain.End - domain.Start;
      const double length = pairing.Tolerance.Length;
      const double clear  = 2.0 * length * Margin;
      const auto   place  = [&pairing](const Hit &hit) -> double
      {
        return pairing.Second ? hit.Second.Along : hit.First.Along;
      };
      for (std::size_t one = 0; one < meeting.Hits.size(); ++one)
      {
        for (std::size_t two = one + 1; two < meeting.Hits.size(); ++two)
        {
          const Hit   &first  = meeting.Hits[one];
          const Hit   &second = meeting.Hits[two];
          const double from   = place(first);
          const double middle = from + ((place(second) - from) / 2.0);
          const double round  = middle + (whole / 2.0);
          const double back   = round > domain.End ? round - whole : round;
          if (!(Distance(first.At, second.At) > length * Margin)
              || !(Away(lane, At(shape, middle)) > clear)
              || !(Away(lane, At(shape, back)) > clear))
          {
            return false;
          }
        }
      }
      return true;
    }

    /**
     *  Two Bézier curves met by the twin kernel, its hits placed on each curve at the places it
     *  proved, each the double nearest the true one.
     *
     *  @param  curves    The two curves.
     *  @param  tolerance The tolerance, no finer than the curves' grain.
     *  @param  into      Where the hits go, empty on entry, and left so where it is not sure.
     *  @return           Whether it was sure; nothing where either curve is not a quadratic or
     *                    cubic Bézier curve.
     */
    auto Twinned(const Curves &curves, const Tolerance &tolerance, Intersections &into)
        -> std::optional<bool>
    {
      std::array<Section, 2> sections = {};
      for (std::size_t which = 0; which < 2; ++which)
      {
        const Curve &curve = curves.at(which);
        if (const auto *quadratic = std::get_if<Bezier::Quadratic>(&curve))
        {
          sections.at(which) = Sectioned(*quadratic);
        }
        else if (const auto *cubic = std::get_if<Bezier::Cubic>(&curve))
        {
          sections.at(which) = Sectioned(*cubic);
        }
        else
        {
          return std::nullopt;
        }
      }
      const std::optional<Twains> twains = Twined({sections.data(), &sections[1]}, tolerance);
      if (!twains)
      {
        return false;
      }
      for (const Twain &twain : *twains)
      {
        into.Hits.push_back(
            {.At     = twain.At,
             .First  = {.At = At(curves[0], twain.Along[0]), .Along = twain.Along[0]},
             .Second = {.At = At(curves[1], twain.Along[1]), .Along = twain.Along[1]},
             .Kind   = Contact::Cross});
      }
      std::ranges::sort(into.Hits, {}, [](const Hit &hit) -> double { return hit.First.Along; });
      return true;
    }
  } // namespace

  auto Laned(const Curve &curve) -> std::optional<Lane>
  {
    if (const auto *line = std::get_if<Line>(&curve))
    {
      std::optional<Lane> lane = MadeLane(
          line->Start,
          {Sure::Apart(line->End.X, line->Start.X), Sure::Apart(line->End.Y, line->Start.Y)},
          Domain(*line));
      if (lane)
      {
        lane->Shape = &curve;
        lane->Ends  = {line->Start, line->End};
        lane->Count = 2;
      }
      return lane;
    }
    if (const auto *ray = std::get_if<Ray>(&curve))
    {
      std::optional<Lane> lane =
          MadeLane(ray->Start, {Sure::Exact(ray->Way.X), Sure::Exact(ray->Way.Y)}, Domain(*ray));
      if (lane)
      {
        lane->Shape = &curve;
        lane->Ends  = {ray->Start, {}};
        lane->Count = 1;
      }
      return lane;
    }
    if (const auto *xline = std::get_if<Xline>(&curve))
    {
      std::optional<Lane> lane = MadeLane(
          xline->Through, {Sure::Exact(xline->Way.X), Sure::Exact(xline->Way.Y)}, Domain(*xline));
      if (lane)
      {
        lane->Shape = &curve;
      }
      return lane;
    }
    return std::nullopt;
  }

  auto Away(const Lane &lane, Point point) -> double
  {
    /**
     *  The gap to P and the cross product with d are each within 2^-50 of their terms, d rounded
     *  once at the most and every operation once: the carrier's distance less that, over |d| grown
     *  by its own rounding, is never above the true one.
     */
    const Vector gap   = point - lane.Start;
    const double cross = Cross(gap, lane.Rough);
    const double slack =
        0x1p-50 * (std::abs(gap.X * lane.Rough.Y) + std::abs(gap.Y * lane.Rough.X));
    const double off    = (std::abs(cross) - slack) / (lane.Speed * (1.0 + 0x1p-50));
    const double along  = Dot(gap, lane.Rough) / (lane.Speed * lane.Speed);
    const double margin = 0x1p-40 * (1.0 + std::abs(along));
    double       beyond = -std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < lane.Count; ++index)
    {
      const bool before = index == 0 && along < lane.Along.Start - margin;
      const bool after  = index == 1 && along > lane.Along.End + margin;
      if (before || after)
      {
        beyond = Distance(point, lane.Ends.at(index)) * (1.0 - 0x1p-50);
      }
    }
    return std::max(off, beyond);
  }

  auto Unaligned(const Run &first, const Run &second, const Tolerance &tolerance) -> bool
  {
    const Vector way    = first.Way;
    const double square = Dot(way, way);
    const double base   = Dot(second.Start - first.Start, way) / square;
    const double rate   = Dot(second.Way, way) / square;
    const double from   = base + (rate * second.Along.Start);
    const double until  = base + (rate * second.Along.End);
    const Range  span   = {.Start = std::max(first.Along.Start, std::min(from, until)),
                           .End   = std::min(first.Along.End, std::max(from, until))};
    const double length = (span.End - span.Start) * std::sqrt(square);
    if (length < tolerance.Length / 2.0)
    {
      return true;
    }
    const double reach = Length(second.Way);
    const auto   off   = [&second, reach](Point point) -> double
    {
      return std::abs(Cross(point - second.Start, second.Way)) / reach;
    };
    return std::max(off(first.Start + (way * span.Start)), off(first.Start + (way * span.End)))
         > 2.0 * tolerance.Length;
  }

  auto Met(const Curves &curves, const Tolerance &tolerance, Intersections &into) -> bool
  {
    if (const std::optional<bool> ringed = Ringed(curves, tolerance, into))
    {
      return *ringed;
    }
    if (const std::optional<bool> laced = Laced(curves, tolerance, into))
    {
      return *laced;
    }
    if (const std::optional<bool> vaulted = Vaulted(curves, tolerance, into))
    {
      return *vaulted;
    }
    if (const std::optional<bool> braided = Braided(curves, tolerance, into))
    {
      return *braided;
    }
    if (const std::optional<bool> pieced = Pieced(curves, tolerance, into))
    {
      return *pieced;
    }
    if (const std::optional<bool> twinned = Twinned(curves, tolerance, into))
    {
      return *twinned;
    }
    const std::optional<Lane> first  = Laned(curves[0]);
    const std::optional<Lane> second = first ? std::nullopt : Laned(curves[1]);
    if (!first && !second)
    {
      return false;
    }
    const Pairing              pairing = {.Straight  = first ? &*first : &*second,
                                          .Shape     = &curves.at(first ? 1 : 0).get(),
                                          .Second    = first.has_value(),
                                          .Tolerance = tolerance};
    const std::optional<Spots> spots   = Dispatched(pairing);
    const bool placed                  = spots
                                      && std::ranges::all_of(*spots, [&pairing, &into](const Spot &spot) -> bool
                                                             { return Placed(pairing, spot, into); });
    if (!placed || !Apart(pairing, into))
    {
      into.Hits.clear();
      return false;
    }
    std::ranges::sort(into.Hits, {}, [](const Hit &hit) -> double { return hit.First.Along; });
    return true;
  }
} // namespace Geoxy::Swift
