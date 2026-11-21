/***************************************************************************************************
 *  @file      source/touch.hpp
 *  @brief     Where two curves come nearest about a touch, polished in double-double and rounded
 *             once: one point, whichever curve is first and however the plane is turned.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TOUCH_HPP
#define GEOXY_TOUCH_HPP

#include "crossing.hpp"

#include "geoxy/curve.hpp"

#include <array>
#include <optional>

namespace Geoxy::Touches
{
  /**
   *  Where two curves come nearest about a touch: the places, one on each, where their tangents are
   *  parallel and the line between their points meets both square, found by Newton's method from
   *  the touch's places on the part of each the touch is on, a line's, a polyline's segment, a
   *  spline's or a Bézier curve's points and slopes worked out in double-double by de Casteljau's
   *  or de Boor's algorithm, an arc's, an ellipse's or a bulge's from the point of its unit circle
   *  turned; and the point halfway between them rounded once. It hangs on the curves' own numbers
   *  only, not on where the touch was found: the same point to the bit whichever curve is first,
   *  and turned, mirrored or scaled by a power of two with the plane, as the double nearest a point
   *  is; a coordinate under 2^-80 of the size about the touch, lost in double-double's width, is 0.
   *
   *  @param  curves The two curves.
   *  @param  places The touch's place on each.
   *  @param  length How near the two places must come: twice the tolerance, as deep as two
   *                 crossings that are one touch cross.
   *  @return        The point halfway, each place and each curve's point there; nothing where
   *                 Newton's method does not settle, or settles off the part of either curve the
   *                 touch is on, or with the places farther apart than the length; nor where the
   *                 curves bend alike, two lines or two arcs of one circle, which come nearest all
   *                 along, nor on an arc of no sweep.
   */
  auto Nearest(const std::array<const Curve *, 2> &curves, const std::array<double, 2> &places,
               double length) -> std::optional<Crossings::Crossed>;
} // namespace Geoxy::Touches

#endif /* GEOXY_TOUCH_HPP */
