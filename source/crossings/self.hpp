/***************************************************************************************************
 *  @file      source/crossings/self.hpp
 *  @brief     Where a curve crosses or touches itself, and the stretches it runs along itself.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CROSSINGS_SELF_HPP
#define GEOXY_CROSSINGS_SELF_HPP

#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/tolerance.hpp"

namespace Geoxy::Crossings
{
  /**
   *  Where a curve crosses or touches itself, as PEDIT and a region's check ask.
   *
   *  The curve is cut into pieces that cannot cross themselves: a polyline into its segments; a
   *  Bézier curve, and a spline span by span, halved until each piece's control points run one way
   *  in x or in y; an arc or an ellipse of more than a turn into turns. Every two pieces are met
   *  through Intersect's exact kernels, the pairs found by a Grid's cells, and where two pieces
   *  follow one another the end they share is no hit. Hits of several pairs at one place, as where
   *  the curve passes through one of its vertices again, are one; where Intersect tells only that
   *  pieces end there, the hit is a crossing where one pass's ways in and out lie either side of
   *  the other's, decided by Orient, a touch where not, and an end at an open curve's own ends.
   *
   *  @param  curve     The curve.
   *  @param  tolerance Within which two places of it are one.
   *  @return           The hits, the first place of each before its second, both in the curve's own
   *                    parameter, in order along it; the stretches it runs along itself, as pairs
   *                    of stretches; and whether a part of the answer may be missing.
   */
  [[nodiscard]] auto Self(const Curve &curve, const Tolerance &tolerance) -> Intersections;
} // namespace Geoxy::Crossings

#endif /* GEOXY_CROSSINGS_SELF_HPP */
