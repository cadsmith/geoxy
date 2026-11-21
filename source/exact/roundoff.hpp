/***************************************************************************************************
 *  @file      source/exact/roundoff.hpp
 *  @brief     What rounding takes from the difference of two points.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_EXACT_ROUNDOFF_HPP
#define GEOXY_EXACT_ROUNDOFF_HPP

#include "exact/two.hpp"

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <cmath>

namespace Geoxy::Exact
{
  /**
   *  What rounding took from the difference of two points, component by component.
   *
   *  @param  end   The point arrived at.
   *  @param  start The point left.
   *  @return       The roundoff of each component of `end - start`; zero where it is exact.
   */
  inline auto Roundoff(Point end, Point start) -> Vector
  {
    const Vector rounded = end - start;
    return {.X = TwoDiffTail(end.X, start.X, rounded.X),
            .Y = TwoDiffTail(end.Y, start.Y, rounded.Y)};
  }

  /**
   *  Whether a difference was exact: both components of its roundoff zero, of either sign, and
   *  neither NaN, as an overflow on the way to a roundoff leaves it near the top of the doubles.
   *
   *  @param  roundoff The difference's roundoff.
   *  @return          Whether it is the zero vector, NaN excluded.
   */
  inline auto Unrounded(Vector roundoff) -> bool
  {
    return std::fpclassify(roundoff.X) == FP_ZERO && std::fpclassify(roundoff.Y) == FP_ZERO;
  }

  /**
   *  Whether both components of a displacement are zero, of either sign.
   *
   *  @param  vector The displacement.
   *  @return        Whether it is the zero vector.
   */
  constexpr auto IsZero(Vector vector) -> bool
  {
    return IsZero(vector.X) && IsZero(vector.Y);
  }
} // namespace Geoxy::Exact

#endif /* GEOXY_EXACT_ROUNDOFF_HPP */
