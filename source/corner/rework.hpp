/***************************************************************************************************
 *  @file      source/corner/rework.hpp
 *  @brief     A polyline's corners each replaced, in order, where what replaces it fits within what
 *             is left of its two segments.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CORNER_REWORK_HPP
#define GEOXY_CORNER_REWORK_HPP

#include "corner/choose.hpp"

#include "geoxy/fillet.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/tolerance.hpp"

#include <array>
#include <expected>
#include <functional>

namespace Geoxy::Corner
{
  /**
   *  What replaces a corner of a polyline: where the segment before it now ends and the one after
   *  now starts, each with its parameter on its segment, and the bulge of what joins them.
   */
  struct Joint
  {
    std::array<Point, 2>  Ends   = {};  /**< Where the two segments are cut.             */
    std::array<double, 2> Alongs = {};  /**< Each cut's parameter on its segment.        */
    double                Bulge  = 0.0; /**< The bulge from the first cut to the second. */
  };

  /**
   *  What replaces a corner, given the guides of its two segments, the first picked at its start
   *  and the second at its end, and the corner.
   */
  using Joiner = std::function<std::expected<Joint, Misfit>(const Guides &, Point)>;

  /**
   *  A polyline with its corners replaced, in order, where each fits.
   *
   *  @param  polyline  The polyline, of finite numbers.
   *  @param  joiner    What replaces each corner.
   *  @param  tolerance Within which a cut still fits, two segments run on smoothly, and a segment
   *                    left is of no length.
   *  @return           The polyline and the corners left, with why.
   */
  auto Reworked(const Polyline &polyline, const Joiner &joiner, const Tolerance &tolerance)
      -> Reshaped;
} // namespace Geoxy::Corner

#endif /* GEOXY_CORNER_REWORK_HPP */
