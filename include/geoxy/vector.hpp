/***************************************************************************************************
 *  @file      geoxy/vector.hpp
 *  @brief     Displacements in the plane: a direction and a length, with no place of their own.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_VECTOR_HPP
#define GEOXY_VECTOR_HPP

#include "geoxy/api.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <optional>

namespace Geoxy
{
  /**
   *  A displacement in the plane, in drawing units.
   *
   *  A direction and a length, with no place of its own: a Point plus a Vector is a Point, and two
   *  Points subtract to a Vector.
   */
  struct Vector
  {
    double X = 0.0; /**< Along the x axis. */
    double Y = 0.0; /**< Along the y axis. */

    /**
     *  Whether both components are equal exactly, as IEEE 754 compares.
     *
     *  @param  lhs One displacement.
     *  @param  rhs The other.
     *  @return     Whether they are equal: 0 equals -0, and NaN equals nothing.
     */
    [[nodiscard]] friend constexpr auto operator==(const Vector &lhs, const Vector &rhs) -> bool
    {
      return std::equal_to<double>{}(lhs.X, rhs.X) && std::equal_to<double>{}(lhs.Y, rhs.Y);
    }
  };

  /**
   *  The sum of two displacements.
   *
   *  @param  lhs The first displacement.
   *  @param  rhs The second displacement.
   *  @return     Their sum, component by component.
   */
  [[nodiscard]] constexpr auto operator+(Vector lhs, Vector rhs) -> Vector
  {
    return {.X = lhs.X + rhs.X, .Y = lhs.Y + rhs.Y};
  }

  /**
   *  The difference of two displacements.
   *
   *  @param  lhs The displacement subtracted from.
   *  @param  rhs The displacement subtracted.
   *  @return     Their difference, component by component.
   */
  [[nodiscard]] constexpr auto operator-(Vector lhs, Vector rhs) -> Vector
  {
    return {.X = lhs.X - rhs.X, .Y = lhs.Y - rhs.Y};
  }

  /**
   *  The displacement the other way.
   *
   *  @param  vector The displacement.
   *  @return        The same length, the opposite way.
   */
  [[nodiscard]] constexpr auto operator-(Vector vector) -> Vector
  {
    return {.X = -vector.X, .Y = -vector.Y};
  }

  /**
   *  A displacement scaled.
   *
   *  @param  lhs The displacement.
   *  @param  rhs The scale, negative to turn it about.
   *  @return     Each component times the factor.
   */
  [[nodiscard]] constexpr auto operator*(Vector lhs, double rhs) -> Vector
  {
    return {.X = lhs.X * rhs, .Y = lhs.Y * rhs};
  }

  /**
   *  A displacement scaled.
   *
   *  @param  lhs The scale, negative to turn it about.
   *  @param  rhs The displacement.
   *  @return     Each component times the factor.
   */
  [[nodiscard]] constexpr auto operator*(double lhs, Vector rhs) -> Vector
  {
    return rhs * lhs;
  }

  /**
   *  A displacement divided.
   *
   *  @param  lhs The displacement.
   *  @param  rhs The divisor.
   *  @return     Each component over the divisor.
   */
  [[nodiscard]] constexpr auto operator/(Vector lhs, double rhs) -> Vector
  {
    return {.X = lhs.X / rhs, .Y = lhs.Y / rhs};
  }

  /**
   *  The displacement turned a quarter turn counterclockwise, exactly.
   *
   *  @param  vector The displacement.
   *  @return        The same length, a quarter turn to the left.
   */
  [[nodiscard]] constexpr auto Perpendicular(Vector vector) -> Vector
  {
    return {.X = -vector.Y, .Y = vector.X};
  }

  /**
   *  The larger magnitude of the two components.
   *
   *  A length within a factor of the square root of 2, worked out with no arithmetic that can
   *  round.
   *
   *  @param  vector The displacement.
   *  @return        The larger of |X| and |Y|.
   */
  [[nodiscard]] constexpr auto Largest(Vector vector) -> double
  {
    const double x = vector.X < 0.0 ? -vector.X : vector.X;
    const double y = vector.Y < 0.0 ? -vector.Y : vector.Y;
    return std::max(x, y);
  }

  /**
   *  The dot product.
   *
   *  @param  first  The first displacement.
   *  @param  second The second displacement.
   *  @return        The product of their lengths and the cosine of the angle between them.
   */
  [[nodiscard]] constexpr auto Dot(Vector first, Vector second) -> double
  {
    return (first.X * second.X) + (first.Y * second.Y);
  }

  /**
   *  The z component of the cross product.
   *
   *  @param  first  The displacement turned from.
   *  @param  second The displacement turned to.
   *  @return        The product of their lengths and the sine of the angle from the first to the
   *                 second, positive counterclockwise.
   */
  [[nodiscard]] constexpr auto Cross(Vector first, Vector second) -> double
  {
    return (first.X * second.Y) - (first.Y * second.X);
  }

  /**
   *  The length of a displacement.
   *
   *  Components past 1e154 or under 1e-154 are scaled first, so that the square neither overflows
   *  nor underflows.
   *
   *  @param  vector The displacement.
   *  @return        Its length, finite for every finite displacement.
   */
  [[nodiscard]] inline auto Length(Vector vector) -> double
  {
    constexpr double Low     = 0x1.0p-500;
    constexpr double High    = 0x1.0p500;
    const double     largest = Largest(vector);
    if (largest >= Low && largest <= High)
    {
      return std::sqrt(Dot(vector, vector));
    }
    if (!(largest > 0.0) || !std::isfinite(largest))
    {
      return largest;
    }
    const Vector scaled = vector / largest;
    return largest * std::sqrt(Dot(scaled, scaled));
  }

  /**
   *  The displacement of length 1 the same way.
   *
   *  Scaled by its largest component first, so that no finite input overflows or underflows.
   *
   *  @param  vector The displacement.
   *  @return        A vector of length 1, or none for a zero or non-finite displacement.
   */
  [[nodiscard]] GEOXY_API auto Unit(Vector vector) -> std::optional<Vector>;
} // namespace Geoxy

#endif /* GEOXY_VECTOR_HPP */
