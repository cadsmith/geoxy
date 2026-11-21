/***************************************************************************************************
 *  @file      source/exact/batch.hpp
 *  @brief     The orientation test of many points against one line.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_EXACT_BATCH_HPP
#define GEOXY_EXACT_BATCH_HPP

#include "geoxy/point.hpp"

#include <cstdint>
#include <span>

namespace Geoxy::Exact
{
  /**
   *  A line through two points, directed from the first to the second.
   */
  struct Line
  {
    Point First;  /**< The point it leaves.   */
    Point Second; /**< The point it heads to. */
  };

  /**
   *  The exact side of a line each point lies on: Orient(line.First, line.Second, point).
   *
   *  The shape of a hull, a clip or a point location. The first stage runs over blocks of points in
   *  one branch-free loop, which the compiler gives to the vector unit; only the points it leaves
   *  open take the later stages.
   *
   *  @pre          signs has room for a sign per point.
   *  @param line   The line.
   *  @param points The points.
   *  @param signs  Where the signs go, one per point: 1 left, -1 right, 0 on the line or where a
   *                coordinate is NaN or infinite.
   */
  void OrientSigns(const Line &line, std::span<const Point> points, std::span<std::int8_t> signs);
} // namespace Geoxy::Exact

#endif /* GEOXY_EXACT_BATCH_HPP */
