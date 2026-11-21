/***************************************************************************************************
 *  @file      geoxy/orient.hpp
 *  @brief     Which way three points turn, decided exactly.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_ORIENT_HPP
#define GEOXY_ORIENT_HPP

#include "geoxy/api.hpp"
#include "geoxy/point.hpp"

#include <cstdint>

namespace Geoxy
{
  /**
   *  Which way three points turn.
   */
  enum class Orientation : std::int8_t
  {
    Clockwise        = -1, /**< The third is right of the line from the first to the second. */
    Collinear        =  0, /**< The three are on one line, or two or more are one point.     */
    Counterclockwise =  1, /**< The third is left of the line from the first to the second.  */
  };

  /**
   *  Which way three points turn, decided exactly, not within a tolerance.
   *
   *  Shewchuk's adaptive orientation test, exact for every finite double. The same three points in
   *  any order give answers that agree, which no test with a tolerance can promise.
   *
   *  @param  first  The first point.
   *  @param  second The second point.
   *  @param  third  The point whose side of the line from the first to the second is asked.
   *  @return        The orientation; Collinear where a coordinate is NaN or infinite.
   */
  [[nodiscard]] GEOXY_API auto Orient(Point first, Point second, Point third) -> Orientation;
} // namespace Geoxy

#endif /* GEOXY_ORIENT_HPP */
