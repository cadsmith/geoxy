/***************************************************************************************************
 *  @file      source/axis.hpp
 *  @brief     Whether a point lies exactly on a line's carrier or an arc's circle: what tells that
 *             two carriers cross on an axis, a coordinate of zero, which no bound on a rounded
 *             value can tell.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_AXIS_HPP
#define GEOXY_AXIS_HPP

#include "geoxy/arc.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"

/**
 *  Where a crossing's coordinate is not certain, its bound reaching zero, and the other is certain,
 *  the crossing lies on the axis exactly where the point of zero and that other coordinate lies on
 *  both carriers exactly: carriers that meet at one point there, and cross at none other near it,
 *  cross there. Each test is the sign of a sum of exact products of the numbers' differences, of
 *  degree two, held in a few doubles with no rounding: zero exactly where the point is on the
 *  carrier. A number past 2^100 in size, or under 2^-100 but not zero, is not asked: there every
 *  product keeps clear of underflow and overflow, as exactness asks.
 */
namespace Geoxy::Axis
{
  /**
   *  Whether a number may be asked of: zero, or between 2^-100 and 2^100 in size.
   *
   *  @param  value The number.
   *  @return       Whether it may.
   */
  auto Held(double value) -> bool;

  /**
   *  Whether the numerator of a coordinate that is zero, of this size, takes Intersect's short sums
   *  and exact way to their zero, +0: past 2^-280 and below 2^280, well inside the range where they
   *  work it out, 2^-300 to 2^300.
   *
   *  @param  size The numerator's size: M's coordinate times the quotient's denominator.
   *  @return      Whether it does.
   */
  auto Fits(double size) -> bool;

  /**
   *  Whether a point lies on a line's carrier exactly: the cross product of the line's way and the
   *  point's gap from its start zero.
   *
   *  @param  line  The line, of some length.
   *  @param  point The point.
   *  @return       Whether it does; never where a number is not held.
   */
  auto On(const Line &line, Point point) -> bool;

  /**
   *  Whether a point lies on an arc's circle exactly: its gap from the center squared, less the
   *  radius squared, zero.
   *
   *  @param  arc   The arc whose circle.
   *  @param  point The point.
   *  @return       Whether it does; never where a number is not held.
   */
  auto On(const Arc &arc, Point point) -> bool;
} // namespace Geoxy::Axis

#endif /* GEOXY_AXIS_HPP */
