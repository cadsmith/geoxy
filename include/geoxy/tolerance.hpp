/***************************************************************************************************
 *  @file      geoxy/tolerance.hpp
 *  @brief     Tolerances, passed to every call that decides whether two things are one, and the
 *             step of a double.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TOLERANCE_HPP
#define GEOXY_TOLERANCE_HPP

#include "geoxy/api.hpp"
#include "geoxy/box.hpp"
#include "geoxy/point.hpp"

namespace Geoxy
{
  /**
   *  How far apart two things may be and still count as one.
   *
   *  Passed to every call that decides it; geoxy holds none of its own, since a model's, a pick
   *  aperture's and a command's fuzz are three scales in one program.
   */
  struct Tolerance
  {
    double Length  = 0.0; /**< The distance within which two points are one, in drawing units. */
    double Radians = 0.0; /**< The angle within which two directions are one.                  */

    /**
     *  A tolerance for a model whose extents are a box.
     *
     *  1e-10 of its larger side, the ratio of AutoCAD's own default to a drawing of unit size, and
     *  no less than 16 steps of a double at its farthest corner, so that it means something at
     *  survey coordinates. The angle is that length over that side.
     *
     *  @param  extents The model's extents.
     *  @return         The tolerance, finite always: a box with no extent gives the steps alone,
     *                  and one with a side or a corner past every double the largest double, at the
     *                  ratio's angle.
     */
    [[nodiscard]] GEOXY_API static auto Of(const Box &extents) -> Tolerance;
  };

  /**
   *  Whether two points are one within a tolerance.
   *
   *  @param  first     One point.
   *  @param  second    The other.
   *  @param  tolerance The tolerance; a length of 0 means exactly.
   *  @return           Whether they are within its length of each other, its edge included.
   */
  [[nodiscard]] GEOXY_API auto Same(Point first, Point second, const Tolerance &tolerance) -> bool;

  /**
   *  The step of a double: the smallest difference that means anything there.
   *
   *  @param  value The double.
   *  @return       The step from it to the next double away from zero.
   */
  [[nodiscard]] GEOXY_API auto Grain(double value) -> double;

  /**
   *  The step of a point's coordinates.
   *
   *  @param  point The point.
   *  @return       The larger of the steps of its two coordinates.
   */
  [[nodiscard]] GEOXY_API auto Grain(Point point) -> double;
} // namespace Geoxy

#endif /* GEOXY_TOLERANCE_HPP */
