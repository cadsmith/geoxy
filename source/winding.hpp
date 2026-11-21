/***************************************************************************************************
 *  @file      source/winding.hpp
 *  @brief     How many times a curve run between two points crosses the ray from a point, each
 *             crossing decided exactly on the numbers stored.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_WINDING_HPP
#define GEOXY_WINDING_HPP

#include "geoxy/curve.hpp"
#include "geoxy/point.hpp"

#include <cstdint>

namespace Geoxy
{
  /**
   *  A curve run from one point to another, the two points standing for its ends: a loop's edge
   *  from its start to its end, or an arrangement's piece between the vertices its ends merged
   *  into, so that each part of a chain is counted from the very point the next starts from.
   */
  struct Traverse
  {
    const Curve *Shape = nullptr; /**< The curve, run forward from its start. */
    Point        From  = {};      /**< The point that stands for its start.   */
    Point        To    = {};      /**< The point that stands for its end.     */
  };

  /**
   *  How many times a traverse crosses the ray from a point in the direction of x, upward +1 and
   *  downward -1, a crossing at the height of an end counted on the side above it alone, as
   *  Sunday's rule counts a polygon's sides.
   *
   *  A line is counted by the exact side of it the point lies on. An arc is cut at its quarter
   *  turns into parts that rise or fall all along, each within one half of its circle; where a part
   *  spans the point's height, its crossing is right of the point by the side of the center the
   *  point lies on and, past the center, by whether the point is inside the circle, the square of
   *  its distance against the radius's worked out in sums of doubles that never round. Every whole
   *  turn of an arc counts once where the point is inside its circle. Any other curve is halved
   *  where a part is not shorter than half the distance from the point to where it starts, so that
   *  no part can wind about the point, and each part counted as the line between its ends.
   *
   *  The point must be farther from the traverse than its ends are from the curve's, so that the
   *  count is that of a closed chain of courses about it: the chain's winding number.
   *
   *  @param  traverse The traverse: a line, an arc, or any curve with an end.
   *  @param  point    The point.
   *  @return          The crossings, signed; for an arc past 2^24 turns, as of 2^24.
   */
  auto Ascents(const Traverse &traverse, Point point) -> std::int64_t;
} // namespace Geoxy

#endif /* GEOXY_WINDING_HPP */
