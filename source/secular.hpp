/***************************************************************************************************
 *  @file      source/secular.hpp
 *  @brief     The place of an ellipse nearest a point, by Eberly's secular equation in the
 *             ellipse's principal frame, polished on its own map.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SECULAR_HPP
#define GEOXY_SECULAR_HPP

#include "geoxy/ellipse.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"

#include <array>
#include <cstddef>
#include <optional>

namespace Geoxy
{
  /**
   *  The place of an ellipse swept round no more than once nearest a point, with no foot listed.
   *
   *  In its principal frame, where its unit circle turned by some phi goes to its longest semi-axis
   *  a and turned a quarter more to its shortest b, a point (x, y) has its nearest place at (a² x /
   *  (s + a² − b²), b² y / s), s the one root past zero of the secular equation G(s) = (ax / (s +
   *  a² − b²))² + (by / s)² − 1 (Eberly, "Distance from a Point to an Ellipse, an Ellipsoid, or a
   *  Hyperellipsoid", 2013). G falls and bends up there, so Newton's method run from a place G is
   *  not below never passes the root. Inside the evolute, (ax)^(2/3) + (by)^(2/3) < (a² −
   *  b²)^(2/3), the distance has a second least, across the major axis, at a root of the same
   *  equation with s turned below zero: it is sought only where an arc leaves out the first, or
   *  where the point is so near the axis that which side is nearer is rounding's to say. Each place
   *  is polished by a step of Newton's method on the map itself, as the frame is the ellipse's to a
   *  rounding only, and the nearest of them and of the arc's ends is the one.
   *
   *  @param  ellipse The ellipse, swept round no more than once.
   *  @param  point   The point.
   *  @return         The nearest place, the first along it of those as near as rounding tells;
   *                  nothing for a circle, which has no principal frame, every turn of it one, and
   *                  for a sweep of nothing.
   */
  auto Closest(const Ellipse &ellipse, Point point) -> std::optional<Place>;

  /**
   *  The angles of an ellipse's unit circle at which the perpendiculars from a point meet it, every
   *  one of them, at most four.
   */
  struct Footholds
  {
    std::array<double, 4> Angles = {}; /**< The angles, in radians, the first Count of them. */
    std::size_t           Count  = 0;  /**< How many there are, 2 or 4.                      */
  };

  /**
   *  Every foot of the perpendiculars from a point to a whole ellipse, by the roots of Eberly's
   *  secular equation, where they are sure.
   *
   *  In the principal frame the feet are (ax / (s + a² − b²), by / s) on the unit circle, for every
   *  real root s of G(s) = (ax / (s + a² − b²))² + (by / s)² − 1 but its poles, the point's terms
   *  taken not below nothing and the signs given back: G falls and bends up between its poles, so
   *  there is one root past nothing, the point's own quarter, and one below b² − a², the quarter
   *  across from it; between them G bends up from both poles to its least, (cbrt((ax)²) +
   *  cbrt((by)²))³ / (a² − b²)² − 1, and has two roots where that is below nothing, the point
   *  within the evolute, in the quarter across the major axis. Each is found by Newton's method
   *  from a place where G is not below nothing, as Closest finds the first, the last two each from
   *  its pole's side, and each foot is turned to the map's unit circle and polished by a step of
   *  Newton's method on the map itself.
   *
   *  @param  ellipse The ellipse; its sweep is not looked at.
   *  @param  point   The point.
   *  @return         The angles, in no order; nothing where they are not sure: for a circle, an
   *                  ellipse so thin, below 2^-20 of its length across, that its frame is lost in
   *                  rounding, a point within 2^-40 of the sizes of an axis, where a root comes to
   *                  a pole, or within 2^-30 of the evolute, where two roots come together.
   */
  auto Plumbs(const Ellipse &ellipse, Point point) -> std::optional<Footholds>;
} // namespace Geoxy

#endif /* GEOXY_SECULAR_HPP */
