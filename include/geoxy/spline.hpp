/***************************************************************************************************
 *  @file      geoxy/spline.hpp
 *  @brief     Splines: non-uniform rational B-splines.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SPLINE_HPP
#define GEOXY_SPLINE_HPP

#include "geoxy/api.hpp"
#include "geoxy/point.hpp"

#include <vector>

namespace Geoxy
{
  /**
   *  The highest degree of spline geoxy holds: 31, where AutoCAD writes 11 at most.
   *
   *  A point of a spline is worked out in a fixed array of a degree and one points, so that no
   *  question of a spline allocates.
   */
  inline constexpr int MaxDegree = 31;

  /**
   *  A spline: a non-uniform rational B-spline, as DXF's SPLINE holds one.
   *
   *  Its parameter runs over its knots, from the one at its degree to the one past its last pole.
   *  With no weights it is polynomial, every weight 1.
   */
  struct Spline
  {
    int                 Degree  = 3;  /**< The degree.                 */
    std::vector<double> Knots   = {}; /**< The knots, none falling.    */
    std::vector<Point>  Poles   = {}; /**< The control points.         */
    std::vector<double> Weights = {}; /**< A weight per pole, or none. */

    /**
     *  Whether two splines are equal exactly, knot for knot and pole for pole.
     *
     *  @param  lhs One spline.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend auto operator==(const Spline &lhs, const Spline &rhs) -> bool
    {
      return lhs.Degree  == rhs.Degree
          && lhs.Knots   == rhs.Knots
          && lhs.Poles   == rhs.Poles
          && lhs.Weights == rhs.Weights;
    }
  };

  /**
   *  Whether a spline is one. A spline that is not is a curve of nothing to every other question,
   *  as a polyline with no vertices is: no length, no area, no box, and it meets nothing.
   *
   *  A degree of at least 1 and at most MaxDegree, more poles than its degree, a knot more than its
   *  poles and its degree, its knots finite and none falling, a span of its parameter not empty,
   *  and a weight per pole, finite and positive, or none.
   *
   *  It is also of one piece: no knot within its domain repeated more than its degree, where the
   *  curve would jump from one pole to another, and none at its ends repeated more than its degree
   *  and one, which would leave a pole nothing weighs. And its weights are held by doubles: the
   *  weights about each span, brought by one power of two as near 1 as keeps the lightest a double
   *  of all its digits, leave every pole times its weight finite, so that every point of it is
   *  worked out from numbers. Weights some 2^2045 apart in one span are past that, and weights less
   *  far apart beside poles far from the origin.
   *
   *  @param  spline The spline.
   *  @return        Whether it is valid.
   */
  [[nodiscard]] GEOXY_API auto Valid(const Spline &spline) -> bool;
} // namespace Geoxy

#endif /* GEOXY_SPLINE_HPP */
