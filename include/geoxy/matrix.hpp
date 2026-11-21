/***************************************************************************************************
 *  @file      geoxy/matrix.hpp
 *  @brief     Affine maps of the plane: moves, turns, scales, mirrors and what composes from them.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_MATRIX_HPP
#define GEOXY_MATRIX_HPP

#include "geoxy/angle.hpp"
#include "geoxy/api.hpp"
#include "geoxy/box.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <optional>

namespace Geoxy
{
  /**
   *  An affine map of the plane, the identity until set.
   *
   *  Held as where it sends the two unit axes and the origin: a point (x, y) goes to `Origin + x X
   *  + y Y`.
   */
  struct Matrix
  {
    Vector X      = {.X = 1.0, .Y = 0.0}; /**< Where the unit x axis goes. */
    Vector Y      = {.X = 0.0, .Y = 1.0}; /**< Where the unit y axis goes. */
    Point  Origin = {};                   /**< Where the origin goes.      */

    /**
     *  Whether all six numbers are equal exactly.
     *
     *  @param  lhs One map.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Matrix &lhs, const Matrix &rhs) -> bool
    {
      return lhs.X == rhs.X && lhs.Y == rhs.Y && lhs.Origin == rhs.Origin;
    }
  };

  /**
   *  The map that moves every point by a displacement.
   *
   *  @param  displacement The displacement.
   *  @return              The move.
   */
  [[nodiscard]] constexpr auto Translation(Vector displacement) -> Matrix
  {
    return {.Origin = Point{} + displacement};
  }

  /**
   *  The map that scales about the origin, by one factor along each axis.
   *
   *  @param  alongX The scale along x, negative to mirror.
   *  @param  alongY The scale along y, negative to mirror.
   *  @return        The scale.
   */
  [[nodiscard]] constexpr auto Scaling(double alongX, double alongY) -> Matrix
  {
    return {.X = {.X = alongX, .Y = 0.0}, .Y = {.X = 0.0, .Y = alongY}};
  }

  /**
   *  The map that scales about the origin by one factor.
   *
   *  @param  factor The scale.
   *  @return        The scale.
   */
  [[nodiscard]] constexpr auto Scaling(double factor) -> Matrix
  {
    return Scaling(factor, factor);
  }

  /**
   *  The map that turns about the origin, counterclockwise.
   *
   *  @param  angle The angle turned.
   *  @return       The turn, exact at every quarter turn.
   */
  [[nodiscard]] GEOXY_API auto Rotation(Angle angle) -> Matrix;

  /**
   *  The map that turns about a point, counterclockwise.
   *
   *  @param  angle The angle turned.
   *  @param  about The point that stays.
   *  @return       The turn.
   */
  [[nodiscard]] GEOXY_API auto Rotation(Angle angle, Point about) -> Matrix;

  /**
   *  The map that mirrors across the line through two points.
   *
   *  @param  first  One point on the line.
   *  @param  second Another.
   *  @return        The mirror, or none when the points are the same point and make no line.
   */
  [[nodiscard]] GEOXY_API auto Mirror(Point first, Point second) -> std::optional<Matrix>;

  /**
   *  Two maps composed.
   *
   *  @param  lhs The map applied second.
   *  @param  rhs The map applied first.
   *  @return     The map that applies the right one, then the left one.
   */
  [[nodiscard]] GEOXY_API auto operator*(const Matrix &lhs, const Matrix &rhs) -> Matrix;

  /**
   *  A point mapped.
   *
   *  @param  lhs The map.
   *  @param  rhs The point.
   *  @return     Where the map sends it.
   */
  [[nodiscard]] GEOXY_API auto operator*(const Matrix &lhs, Point rhs) -> Point;

  /**
   *  A displacement mapped: turned and scaled, never moved.
   *
   *  @param  lhs The map.
   *  @param  rhs The displacement.
   *  @return     Where the map sends it.
   */
  [[nodiscard]] GEOXY_API auto operator*(const Matrix &lhs, Vector rhs) -> Vector;

  /**
   *  A box mapped: the box of its four corners mapped.
   *
   *  A curve's own box under a map is the curve's business: a mapped box holds it but may be
   *  larger, as AutoCAD's extents of an INSERT are.
   *
   *  @param  lhs The map.
   *  @param  rhs The box.
   *  @return     The box of the mapped corners; empty stays empty.
   */
  [[nodiscard]] GEOXY_API auto operator*(const Matrix &lhs, const Box &rhs) -> Box;

  /**
   *  The factor by which a map scales areas.
   *
   *  @param  map The map.
   *  @return     The determinant, negative when the map mirrors.
   */
  [[nodiscard]] GEOXY_API auto Determinant(const Matrix &map) -> double;

  /**
   *  The map that undoes another.
   *
   *  @param  map The map.
   *  @return     The inverse, or none when the map's axes are parallel, its determinant 0 as worked
   *              out exactly, or the inverse does not fit in a double.
   */
  [[nodiscard]] GEOXY_API auto Inverse(const Matrix &map) -> std::optional<Matrix>;

  /**
   *  The most a map lengthens a displacement.
   *
   *  Its largest singular value: the radius a circle of radius 1 takes at its widest. Not its
   *  longest axis, which falls short on a shear and on a circle seen edge on.
   *
   *  @param  map The map.
   *  @return     The stretch, finite for every finite map.
   */
  [[nodiscard]] GEOXY_API auto Stretch(const Matrix &map) -> double;
} // namespace Geoxy

#endif /* GEOXY_MATRIX_HPP */
