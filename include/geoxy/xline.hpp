/***************************************************************************************************
 *  @file      geoxy/xline.hpp
 *  @brief     Xlines: a line through a point that never ends either way.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_XLINE_HPP
#define GEOXY_XLINE_HPP

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

namespace Geoxy
{
  /**
   *  An xline: a line through a point that runs both ways without end, as DXF's XLINE.
   *
   *  Its way gives the direction and the unit of its parameter, whatever its length.
   */
  struct Xline
  {
    Point  Through = {}; /**< A point it runs through, where its parameter is 0. */
    Vector Way     = {}; /**< The way it runs.                                   */

    /**
     *  Whether two xlines are equal exactly.
     *
     *  @param  lhs One xline.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Xline &lhs, const Xline &rhs) -> bool
    {
      return lhs.Through == rhs.Through && lhs.Way == rhs.Way;
    }
  };
} // namespace Geoxy

#endif /* GEOXY_XLINE_HPP */
