/***************************************************************************************************
 *  @file      source/met/overlap.cpp
 *  @brief     The stretches two curves share: snapped to the ends of parts, settled as kept or as a
 *             band about a touch, and joined where one runs on into the next.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "met/overlap.hpp"

#include "met/floor.hpp"
#include "met/judge.hpp"
#include "met/parts.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "trig.hpp"

#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <span>
#include <utility>

namespace Geoxy::Met
{
  auto Joined(List<Overlap> overlaps) -> List<Overlap>
  {
    std::ranges::sort(overlaps, {},
                      [](const Overlap &overlap) -> double { return overlap.First.Start; });
    List<Overlap> joined;
    for (const Overlap &overlap : overlaps)
    {
      if (!joined.empty() && Continues(joined.back(), overlap))
      {
        Overlap &last  = joined.back();
        last.First.End = std::max(last.First.End, overlap.First.End);
        last.Second    = {.Start = std::min(last.Second.Start, overlap.Second.Start),
                          .End   = std::max(last.Second.End, overlap.Second.End)};
      }
      else
      {
        joined.push_back(overlap);
      }
    }
    return joined;
  }

  auto Boundaries(const Scene &scene, std::array<const List<Part> *, 2> parts) -> List<double>
  {
    List<double> bounds;
    for (std::size_t side = 0; side < parts.size(); ++side)
    {
      for (const Part &part : *parts.at(side))
      {
        const Range range = Bounds(part);
        for (const double end : {range.Start, range.End})
        {
          const double along = Placed(part, end);
          if (!std::isfinite(along))
          {
            continue;
          }
          if (side == 0)
          {
            bounds.push_back(along);
            continue;
          }
          const Point point = At(*scene.Second, along);
          if (Remote(scene.FirstFloors, point, scene.Tolerance.Length))
          {
            continue;
          }
          const Place nearest = Nearest(*scene.First, point);
          if (Distance(point, nearest.At) <= scene.Tolerance.Length)
          {
            bounds.push_back(nearest.Along);
          }
        }
      }
    }
    bounds.push_back(-std::numeric_limits<double>::infinity());
    bounds.push_back(std::numeric_limits<double>::infinity());
    std::ranges::sort(bounds);
    return bounds;
  }

  namespace
  {
    /**
     *  Whether the curves stay within the tolerance between two places of the first: a quarter,
     *  half and three quarters of the way.
     *
     *  @param  scene The curves being met.
     *  @param  span  The two places, as parameters of the first.
     *  @return       Whether they do.
     */
    auto Holds(const Scene &scene, Range span) -> bool
    {
      constexpr std::array<double, 3> Parts = {0.5, 0.25, 0.75};
      return std::ranges::all_of(
          Parts,
          [&scene, span](double part) -> bool
          {
            const Point point = At(*scene.First, span.Start + ((span.End - span.Start) * part));
            return !Remote(scene.SecondFloors, point, scene.Tolerance.Length)
                && Distance(point, Nearest(*scene.Second, point).At) <= scene.Tolerance.Length;
          });
    }

    /**
     *  Whether a place of the first curve is at a boundary, within the tolerance of where a part of
     *  either curve ends, or out at infinity.
     *
     *  @param  scene  The curves being met.
     *  @param  bounds Where the curves' parts end, as parameters of the first, in order.
     *  @param  along  The place.
     *  @return        Whether it is.
     */
    auto Bounded(const Scene &scene, const List<double> &bounds, double along) -> bool
    {
      if (!std::isfinite(along))
      {
        return true;
      }
      const Point end = At(*scene.First, along);
      return std::ranges::any_of(
          bounds, [&scene, end](double bound) -> bool
          { return Distance(At(*scene.First, bound), end) <= scene.Tolerance.Length; });
    }

    /**
     *  Whether an end of a shared stretch is where it should end: at a boundary, or where either
     *  curve all but stalls, as where one turns back on itself.
     *
     *  @param  scene   The curves being met.
     *  @param  bounds  Where the curves' parts end, as parameters of the first, in order.
     *  @param  overlap The stretch, and which of its ends, 0 its start along the first.
     *  @return         Whether it is.
     */
    auto Rests(const Scene &scene, const List<double> &bounds,
               std::pair<const Overlap &, std::size_t> overlap) -> bool
    {
      constexpr double Stall    = 1e-3;
      const auto &[shared, end] = overlap;
      const auto stalls         = [](const Curve &curve, Range span, double along) -> bool
      {
        const double pace = Length(Derivative(curve, span.Start + ((span.End - span.Start) / 2)));
        return Length(Derivative(curve, along)) <= Stall * pace;
      };
      const double first  = end == 0 ? shared.First.Start : shared.First.End;
      const double second = (end == 0) == shared.Along ? shared.Second.Start : shared.Second.End;
      return Bounded(scene, bounds, first)
          || stalls(*scene.First, shared.First, first)
          || stalls(*scene.Second, shared.Second, second);
    }

  } // namespace

  auto Snapped(const Scene &scene, const List<double> &bounds, Overlap overlap) -> Overlap
  {
    const Range was    = overlap.First;
    const auto  after  = std::ranges::upper_bound(bounds, overlap.First.End);
    const auto  before = std::ranges::lower_bound(bounds, overlap.First.Start);
    if (!Rests(scene, bounds, {overlap, 1})
        && Holds(scene, {.Start = overlap.First.End, .End = *after}))
    {
      overlap.First.End = *after;
    }
    if (!Rests(scene, bounds, {overlap, 0})
        && Holds(scene, {.Start = *std::prev(before), .End = overlap.First.Start}))
    {
      overlap.First.Start = *std::prev(before);
    }
    for (const double end : {overlap.First.Start, overlap.First.End})
    {
      if (!Exactly(end, was.Start) && !Exactly(end, was.End))
      {
        const double along = Nearest(*scene.Second, At(*scene.First, end)).Along;
        overlap.Second     = {.Start = std::min(overlap.Second.Start, along),
                              .End   = std::max(overlap.Second.End, along)};
      }
    }
    return overlap;
  }

  auto Settled(const Scene &scene, const List<double> &bounds, const Overlap &overlap) -> bool
  {
    return Rests(scene, bounds, {overlap, 0}) && Rests(scene, bounds, {overlap, 1});
  }

  namespace
  {
    /**
     *  How slowly a curve moves anywhere along a stretch of its parameter: the slowest of nine
     *  places, then the least speed between its neighbours by golden-section search, so that a
     *  curve all but stalled between two places, as at a cusp rounded by its numbers, is found
     *  whichever places the stretch's ends give and whichever curve the places are taken on.
     *
     *  @param  curve   The curve.
     *  @param  stretch The stretch.
     *  @return         The least speed found.
     */
    auto Slowest(const Curve &curve, Range stretch) -> double
    {
      constexpr int Places = 8;
      const double  width  = (stretch.End - stretch.Start) / Places;
      const auto    speed  = [&curve](double along) -> double
      {
        return Length(Derivative(curve, along));
      };
      int    slowest = 0;
      double least   = speed(stretch.Start);
      for (int place = 1; place <= Places; ++place)
      {
        const double here = speed(stretch.Start + (width * place));
        if (here < least)
        {
          slowest = place;
          least   = here;
        }
      }
      constexpr double Golden = 0.6180339887498949;
      double           low    = stretch.Start + (width * std::max(slowest - 1, 0));
      double           high   = stretch.Start + (width * std::min(slowest + 1, Places));
      for (int step = 0; step < 64; ++step)
      {
        const double left  = high - (Golden * (high - low));
        const double right = low + (Golden * (high - low));
        const double one   = speed(left);
        const double other = speed(right);
        least              = std::min({least, one, other});
        if (one < other)
        {
          high = right;
        }
        else
        {
          low = left;
        }
      }
      return least;
    }

    /**
     *  Whether two curves are tangent somewhere along a band, the stretch where they stay within
     *  the tolerance about a contact: at nine places of it, their directions within the angle at
     *  one; or turning one way from each other at one place and the other way at another, so
     *  parallel between; or one of them all but stalled somewhere in it, sought between the places
     *  on each curve's own stretch, as at a cusp, whose direction turns onto the tangent as its
     *  speed falls to nothing. Where along the band a hit's candidate is polished to does not then
     *  decide its kind.
     *
     *  @param  scene The curves being met.
     *  @param  band  The band, as a stretch of each.
     *  @return       Whether they are.
     */
    auto Tangent(const Scene &scene, const Overlap &band) -> bool
    {
      constexpr int Places = 8;
      bool          left   = false;
      bool          right  = false;
      double        slow   = std::numeric_limits<double>::infinity();
      double        fast   = 0.0;
      for (int place = 0; place <= Places; ++place)
      {
        const double along =
            band.First.Start + ((band.First.End - band.First.Start) * place / Places);
        const Vector way = Derivative(*scene.First, along);
        const Vector other =
            Derivative(*scene.Second, Nearest(*scene.Second, At(*scene.First, along)).Along);
        const double turn = Cross(way, other);
        if (std::abs(turn) <= Trig::Sin(scene.Tolerance.Radians) * Length(way) * Length(other))
        {
          return true;
        }
        left  = left || turn > 0.0;
        right = right || turn < 0.0;
        for (const double speed : {Length(way), Length(other)})
        {
          slow = std::min(slow, speed);
          fast = std::max(fast, speed);
        }
      }
      constexpr double Stall = 1e-3;
      slow =
          std::min({slow, Slowest(*scene.First, band.First), Slowest(*scene.Second, band.Second)});
      return (left && right) || slow <= Stall * fast;
    }

  } // namespace

  void Touched(const Scene &scene, std::span<const Overlap> bands, std::span<Hit> hits)
  {
    for (const Overlap &band : bands)
    {
      if (!Tangent(scene, band))
      {
        continue;
      }
      for (Hit &hit : hits)
      {
        const bool within =
            (band.First.Start <= hit.First.Along && hit.First.Along <= band.First.End)
            || Distance(hit.First.At, At(*scene.First, band.First.Start)) <= scene.Tolerance.Length
            || Distance(hit.First.At, At(*scene.First, band.First.End)) <= scene.Tolerance.Length;
        if (within && hit.Kind == Contact::Cross)
        {
          hit.Kind = Contact::Touch;
        }
      }
    }
  }
} // namespace Geoxy::Met
