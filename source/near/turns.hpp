/***************************************************************************************************
 *  @file      source/near/turns.hpp
 *  @brief     The place of an arc nearest a point, told by its angles where they tell it, and by
 *             its feet and ends where not.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_NEAR_TURNS_HPP
#define GEOXY_NEAR_TURNS_HPP

#include "near/feet.hpp"
#include "scratch.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace Geoxy::Near
{
  /**
   *  A point asked about a polyline and the least distance found so far.
   */
  struct Asking
  {
    Point  From  = {};  /**< The point.                                           */
    double Least = 0.0; /**< The least distance so far.                           */
    double Size  = 0.0; /**< The point's largest coordinate and twice the poles'. */
  };

  /**
   *  Whether a polyline's segment can be passed by as its nearest place is sought: each foot it has
   *  lying clearly off it, or clearly on it, away from its ends and clearly farther from the point
   *  than the nearest so far. Solved, it would give no nearer place, nor a foot by a vertex.
   *
   *  A straight segment's foot is at (p - a)·(b - a) / |b - a|² of it, |(b - a) × (p - a)| / |b -
   *  a| from the point. A bulge's feet are where the line through its circle's center and the point
   *  meets the circle, |d - r| and d + r from the point, each on the arc where it lies on the
   *  bulge's side of the chord, (b - a) × (f - a) of the bulge's sign the other way. The margins
   *  are 2^-40 of the sizes at hand and more, beside roundings of some 2^-50 of them, in these
   *  doubles, in Bent's circle and angles and in At's points.
   *
   *  @param  from   The vertex the segment starts at.
   *  @param  end    Where it ends.
   *  @param  asking The point and the least distance so far.
   *  @return        Whether it can be passed by; not where anything is not finite, nor at sizes
   *                 past 2^400 or below 2^-400, where squares leave the doubles.
   */
  auto Bypassed(const Vertex &from, Point end, const Asking &asking) -> bool;

  /**
   *  A polyline segment's feet as Around finds them, with no list: at most two, a straight
   *  segment's one and a bulge's near and far.
   */
  struct Brace
  {
    std::array<double, 2> Alongs = {}; /**< The feet's parameters, the first Count. */
    std::array<double, 2> Sides  = {}; /**< A bulge's: 1 toward the point, -1 away. */
    std::size_t           Count  = 0;  /**< How many there are.                     */
  };

  /**
   *  A polyline segment's feet as Around finds them, to the bit, with no list: of an angle a
   *  bulge's circle passes twice, at its two ends, only the first, which is by a vertex as the
   *  second is, and has the polyline weighed as Everywhere weighs it all the same.
   *
   *  @param  segment The segment.
   *  @param  point   The point.
   *  @return         The feet, unordered.
   */
  auto Braced(const Segment &segment, Point point) -> Brace;

  /**
   *  The place of a curve with two ends nearest a point: the nearest of its feet and its ends.
   *
   *  @tparam CurveT The curve's kind.
   *  @param  curve  The curve.
   *  @param  point  The point.
   *  @return        The place.
   */
  template <typename CurveT>
  auto Bounded(const CurveT &curve, Point point) -> Place
  {
    List<Place> places = Listed(curve, point);
    const Range domain = Domain(curve);
    places.push_back({.At = At(curve, domain.Start), .Along = domain.Start});
    places.push_back({.At = At(curve, domain.End), .Along = domain.End});
    return Least(point, places);
  }

  /**
   *  The nearest place on an arc to a point: by the angles where they tell it; else its start, its
   *  feet in order and its end, each weighed in turn, the first of the nearest along the arc, as
   *  Least takes it from them sorted, with nothing gathered.
   *
   *  @param  arc   The arc.
   *  @param  point The point.
   *  @param  ends  Its start and end, as At gives them, or none.
   *  @return       The place.
   */
  auto Endpoint(const Arc &arc, Point point, const std::array<Point, 2> *ends) -> Place;
} // namespace Geoxy::Near

#endif /* GEOXY_NEAR_TURNS_HPP */
