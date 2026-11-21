/***************************************************************************************************
 *  @file      geoxy/ray.hpp
 *  @brief     Rays: a line from a point that never ends.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_RAY_HPP
#define GEOXY_RAY_HPP

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

namespace Geoxy
{
  /**
   *  A ray: a line from a point that runs on one way without end, as DXF's RAY.
   *
   *  Its way gives the direction and the unit of its parameter, whatever its length; a way of no
   *  length leaves the ray its start alone.
   */
  struct Ray
  {
    Point  Start = {}; /**< Where it starts. */
    Vector Way   = {}; /**< The way it runs. */

    /**
     *  Whether two rays are equal exactly.
     *
     *  @param  lhs One ray.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Ray &lhs, const Ray &rhs) -> bool
    {
      return lhs.Start == rhs.Start && lhs.Way == rhs.Way;
    }
  };
} // namespace Geoxy

#endif /* GEOXY_RAY_HPP */
