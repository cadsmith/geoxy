/***************************************************************************************************
 *  @file      source/finite.hpp
 *  @brief     Whether a point's coordinates or a displacement's components are finite.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_FINITE_HPP
#define GEOXY_FINITE_HPP

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <cmath>

namespace Geoxy
{
  /**
   *  Whether both of a point's coordinates are finite.
   *
   *  @param  point The point.
   *  @return       Whether they are.
   */
  inline auto Finite(Point point) -> bool
  {
    return std::isfinite(point.X) && std::isfinite(point.Y);
  }

  /**
   *  Whether both of a displacement's components are finite.
   *
   *  @param  vector The displacement.
   *  @return        Whether they are.
   */
  inline auto Finite(Vector vector) -> bool
  {
    return std::isfinite(vector.X) && std::isfinite(vector.Y);
  }
} // namespace Geoxy

#endif /* GEOXY_FINITE_HPP */
