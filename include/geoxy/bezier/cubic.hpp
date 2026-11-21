/***************************************************************************************************
 *  @file      geoxy/bezier/cubic.hpp
 *  @brief     Cubic Bézier curves, as PostScript fonts and SVG paths draw them.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BEZIER_CUBIC_HPP
#define GEOXY_BEZIER_CUBIC_HPP

#include "geoxy/point.hpp"

namespace Geoxy::Bezier
{
  /**
   *  A cubic Bézier curve: from a start, leaving toward one control point, arriving from another,
   *  to an end.
   *
   *  It stays a cubic under any affine map, its points mapped.
   */
  struct Cubic
  {
    Point Start = {}; /**< Where it starts.                              */
    Point Out   = {}; /**< The control point it leaves its start toward. */
    Point In    = {}; /**< The control point it arrives at its end from. */
    Point End   = {}; /**< Where it ends.                                */

    /**
     *  Whether two curves are equal exactly.
     *
     *  @param  lhs One curve.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Cubic &lhs, const Cubic &rhs) -> bool
    {
      return lhs.Start == rhs.Start && lhs.Out == rhs.Out && lhs.In == rhs.In && lhs.End == rhs.End;
    }
  };
} // namespace Geoxy::Bezier

#endif /* GEOXY_BEZIER_CUBIC_HPP */
