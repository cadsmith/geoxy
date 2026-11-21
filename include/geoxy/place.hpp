/***************************************************************************************************
 *  @file      geoxy/place.hpp
 *  @brief     Places on a curve: a point and the curve's own parameter there, to hand back.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_PLACE_HPP
#define GEOXY_PLACE_HPP

#include "geoxy/point.hpp"

#include <functional>

namespace Geoxy
{
  /**
   *  A place on a curve: the point, and the curve's parameter there, as At takes it back.
   */
  struct Place
  {
    Point  At    = {};  /**< The point.                   */
    double Along = 0.0; /**< The curve's parameter there. */

    /**
     *  Whether both the points and the parameters are equal exactly.
     *
     *  @param  lhs One place.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Place &lhs, const Place &rhs) -> bool
    {
      return lhs.At == rhs.At && std::equal_to<double>{}(lhs.Along, rhs.Along);
    }
  };
} // namespace Geoxy

#endif /* GEOXY_PLACE_HPP */
