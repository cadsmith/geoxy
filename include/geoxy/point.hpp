/***************************************************************************************************
 *  @file      geoxy/point.hpp
 *  @brief     Places in the plane, which move by displacements and never add.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_POINT_HPP
#define GEOXY_POINT_HPP

#include "geoxy/vector.hpp"

#include <functional>
#include <numeric>

namespace Geoxy
{
  /**
   *  A place in the plane, in drawing units.
   *
   *  Two places do not add: `point + point` does not compile. Their difference is a Vector, and a
   *  Point plus a Vector is a Point.
   */
  struct Point
  {
    double X = 0.0; /**< Along the x axis. */
    double Y = 0.0; /**< Along the y axis. */

    /**
     *  Whether both coordinates are equal exactly, as IEEE 754 compares.
     *
     *  @param  lhs One point.
     *  @param  rhs The other.
     *  @return     Whether they are equal; Same compares within a tolerance.
     */
    [[nodiscard]] friend constexpr auto operator==(const Point &lhs, const Point &rhs) -> bool
    {
      return std::equal_to<double>{}(lhs.X, rhs.X) && std::equal_to<double>{}(lhs.Y, rhs.Y);
    }
  };

  /**
   *  The displacement between two points.
   *
   *  @param  lhs The point arrived at.
   *  @param  rhs The point left.
   *  @return     The displacement that takes `rhs` to `lhs`.
   */
  [[nodiscard]] constexpr auto operator-(Point lhs, Point rhs) -> Vector
  {
    return {.X = lhs.X - rhs.X, .Y = lhs.Y - rhs.Y};
  }

  /**
   *  The point a displacement away.
   *
   *  @param  lhs The point.
   *  @param  rhs The displacement.
   *  @return     The point moved by it.
   */
  [[nodiscard]] constexpr auto operator+(Point lhs, Vector rhs) -> Point
  {
    return {.X = lhs.X + rhs.X, .Y = lhs.Y + rhs.Y};
  }

  /**
   *  The point a displacement back.
   *
   *  @param  lhs The point.
   *  @param  rhs The displacement.
   *  @return     The point moved back by it.
   */
  [[nodiscard]] constexpr auto operator-(Point lhs, Vector rhs) -> Point
  {
    return {.X = lhs.X - rhs.X, .Y = lhs.Y - rhs.Y};
  }

  /**
   *  The point halfway between two.
   *
   *  Each coordinate the double nearest it, as std::midpoint rounds it: no finite input overflows,
   *  and below the normal doubles, where a half rounds, a point is its own midpoint.
   *
   *  @param  first  One point.
   *  @param  second The other.
   *  @return        The point halfway.
   */
  [[nodiscard]] constexpr auto Midpoint(Point first, Point second) -> Point
  {
    return {.X = std::midpoint(first.X, second.X), .Y = std::midpoint(first.Y, second.Y)};
  }

  /**
   *  The point a fraction of the way from one point to another.
   *
   *  @param  start    The point at 0.
   *  @param  end      The point at 1.
   *  @param  fraction How far along, past 0 and 1 too.
   *  @return          The point: exactly `start` at 0, and exactly `end` at 1.
   */
  [[nodiscard]] constexpr auto Lerp(Point start, Point end, double fraction) -> Point
  {
    const Vector across = end - start;
    if (fraction <= 0.5)
    {
      return start + (across * fraction);
    }
    return end - (across * (1.0 - fraction));
  }

  /**
   *  The distance between two points.
   *
   *  @param  first  One point.
   *  @param  second The other.
   *  @return        The length of the displacement between them.
   */
  [[nodiscard]] inline auto Distance(Point first, Point second) -> double
  {
    return Length(second - first);
  }
} // namespace Geoxy

#endif /* GEOXY_POINT_HPP */
