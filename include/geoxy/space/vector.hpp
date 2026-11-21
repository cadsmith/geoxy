/***************************************************************************************************
 *  @file      geoxy/space/vector.hpp
 *  @brief     Vectors in space, for what leaves the plane: an extrusion, a normal.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SPACE_VECTOR_HPP
#define GEOXY_SPACE_VECTOR_HPP

#include "geoxy/api.hpp"

#include <functional>
#include <optional>

namespace Geoxy::Space
{
  /**
   *  A displacement or a direction in space.
   */
  struct Vector
  {
    double X = 0.0; /**< Along the x axis. */
    double Y = 0.0; /**< Along the y axis. */
    double Z = 0.0; /**< Along the z axis. */

    /**
     *  Whether all three components are equal exactly.
     *
     *  @param  lhs One vector.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Vector &lhs, const Vector &rhs) -> bool
    {
      return std::equal_to<double>{}(lhs.X, rhs.X)
          && std::equal_to<double>{}(lhs.Y, rhs.Y)
          && std::equal_to<double>{}(lhs.Z, rhs.Z);
    }
  };

  /**
   *  The unit x axis.
   */
  inline constexpr Vector East{.X = 1.0, .Y = 0.0, .Z = 0.0};

  /**
   *  The unit y axis.
   */
  inline constexpr Vector North{.X = 0.0, .Y = 1.0, .Z = 0.0};

  /**
   *  The unit z axis: the normal of the plan, and of an entity whose file gives none.
   */
  inline constexpr Vector Up{.X = 0.0, .Y = 0.0, .Z = 1.0};

  /**
   *  The dot product.
   *
   *  @param  first  One vector.
   *  @param  second The other.
   *  @return        The product of their lengths and the cosine of the angle between them.
   */
  [[nodiscard]] GEOXY_API auto Dot(const Vector &first, const Vector &second) -> double;

  /**
   *  The cross product, right-handed.
   *
   *  @param  first  The vector turned from.
   *  @param  second The vector turned to.
   *  @return        The vector at right angles to both, its length the area they span.
   */
  [[nodiscard]] GEOXY_API auto Cross(const Vector &first, const Vector &second) -> Vector;

  /**
   *  The length of a vector.
   *
   *  @param  vector The vector.
   *  @return        Its length, finite for every finite vector.
   */
  [[nodiscard]] GEOXY_API auto Length(const Vector &vector) -> double;

  /**
   *  The vector of length 1 the same way.
   *
   *  Scaled by its largest component first, so that an extrusion of 1e300, which AutoCAD accepts,
   *  has one.
   *
   *  @param  vector The vector.
   *  @return        A vector of length 1, or none for a zero or non-finite vector.
   */
  [[nodiscard]] GEOXY_API auto Unit(const Vector &vector) -> std::optional<Vector>;
} // namespace Geoxy::Space

#endif /* GEOXY_SPACE_VECTOR_HPP */
