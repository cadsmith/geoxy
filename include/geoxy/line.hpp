/***************************************************************************************************
 *  @file      geoxy/line.hpp
 *  @brief     Lines: the stretch between two points.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_LINE_HPP
#define GEOXY_LINE_HPP

#include "geoxy/point.hpp"

namespace Geoxy
{
  /**
   *  A line: the stretch from one point to another.
   *
   *  Its two ends may be one point; a line of no length stays a line.
   */
  struct Line
  {
    Point Start = {}; /**< Where it starts. */
    Point End   = {}; /**< Where it ends.   */

    /**
     *  Whether two lines are equal exactly, end for end.
     *
     *  @param  lhs One line.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Line &lhs, const Line &rhs) -> bool
    {
      return lhs.Start == rhs.Start && lhs.End == rhs.End;
    }
  };
} // namespace Geoxy

#endif /* GEOXY_LINE_HPP */
