/***************************************************************************************************
 *  @file      source/exact/triangle.hpp
 *  @brief     The three corners a predicate reads.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_EXACT_TRIANGLE_HPP
#define GEOXY_EXACT_TRIANGLE_HPP

#include "geoxy/point.hpp"

namespace Geoxy::Exact
{
  /**
   *  Three corners, in the order a predicate reads them: the triangle whose orientation is asked,
   *  or the circle through them.
   */
  struct Triangle
  {
    Point First;  /**< The first corner.  */
    Point Second; /**< The second corner. */
    Point Third;  /**< The third corner.  */
  };
} // namespace Geoxy::Exact

#endif /* GEOXY_EXACT_TRIANGLE_HPP */
