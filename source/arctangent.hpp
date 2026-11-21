/***************************************************************************************************
 *  @file      source/arctangent.hpp
 *  @brief     The arc tangent of a double as a pair, as headings work it out before rounding.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_ARCTANGENT_HPP
#define GEOXY_ARCTANGENT_HPP

#include "twofold.hpp"

namespace Geoxy
{
  /**
   *  The arc tangent of a double as a pair of doubles, within 2^-68 of itself: the first octant's
   *  pair Heading settles to the nearest double, for t up to 1, and a quarter turn less that of 1 /
   *  t past it.
   *
   *  @param  value The tangent, any double.
   *  @return       Its arc tangent, in [-pi / 2, pi / 2], the low part within a step of a double of
   *                the high one; not a number for none.
   */
  auto ArcTangent(double value) -> Twofold::Two;
} // namespace Geoxy

#endif /* GEOXY_ARCTANGENT_HPP */
