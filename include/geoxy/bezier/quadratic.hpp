/***************************************************************************************************
 *  @file      geoxy/bezier/quadratic.hpp
 *  @brief     Quadratic Bézier curves, as TrueType draws its glyphs.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BEZIER_QUADRATIC_HPP
#define GEOXY_BEZIER_QUADRATIC_HPP

#include "geoxy/point.hpp"

namespace Geoxy::Bezier
{
  /**
   *  A quadratic Bézier curve: from a start toward a control point to an end.
   *
   *  It stays a quadratic under any affine map, its points mapped.
   */
  struct Quadratic
  {
    Point Start   = {}; /**< Where it starts.                            */
    Point Control = {}; /**< The point it heads toward, and leaves from. */
    Point End     = {}; /**< Where it ends.                              */

    /**
     *  Whether two curves are equal exactly.
     *
     *  @param  lhs One curve.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Quadratic &lhs, const Quadratic &rhs)
        -> bool
    {
      return lhs.Start == rhs.Start && lhs.Control == rhs.Control && lhs.End == rhs.End;
    }
  };
} // namespace Geoxy::Bezier

#endif /* GEOXY_BEZIER_QUADRATIC_HPP */
