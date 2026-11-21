/***************************************************************************************************
 *  @file      geoxy/ellipse.hpp
 *  @brief     Elliptical arcs, as the unit circle under an affine map.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_ELLIPSE_HPP
#define GEOXY_ELLIPSE_HPP

#include "geoxy/angle.hpp"
#include "geoxy/matrix.hpp"

namespace Geoxy
{
  /**
   *  An arc of an ellipse: an arc of the unit circle under an affine map.
   *
   *  The map sends the unit circle's center to the ellipse's and its x and y axes to two conjugate
   *  semi-diameters, so any affine map of an ellipse, or of a circle, is an ellipse: a block scaled
   *  unevenly, a plane facing down, a circle seen edge on, which is flat and still drawable. Its
   *  start and sweep are parameters on the unit circle, DXF's ELLIPSE's own, and the sweep is
   *  signed as an arc's is.
   */
  struct Ellipse
  {
    Matrix Map   = {}; /**< The map of the unit circle.                        */
    Angle  Start = {}; /**< Where it starts, a parameter on the unit circle.   */
    Angle  Sweep = {}; /**< How far it runs, a whole turn for a whole ellipse. */

    /**
     *  Whether two ellipses are equal exactly.
     *
     *  @param  lhs One ellipse.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Ellipse &lhs, const Ellipse &rhs) -> bool
    {
      return lhs.Map == rhs.Map && lhs.Start == rhs.Start && lhs.Sweep == rhs.Sweep;
    }
  };
} // namespace Geoxy

#endif /* GEOXY_ELLIPSE_HPP */
