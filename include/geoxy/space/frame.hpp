/***************************************************************************************************
 *  @file      geoxy/space/frame.hpp
 *  @brief     Frames of planes in space, by AutoCAD's arbitrary axis, and their maps to the plan.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SPACE_FRAME_HPP
#define GEOXY_SPACE_FRAME_HPP

#include "geoxy/api.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/space/vector.hpp"

#include <optional>

namespace Geoxy::Space
{
  /**
   *  Three unit axes at right angles, right-handed: a plane's own x and y, and its normal.
   */
  struct Frame
  {
    Vector X = East;  /**< The plane's x axis. */
    Vector Y = North; /**< The plane's y axis. */
    Vector Z = Up;    /**< The plane's normal. */

    /**
     *  Whether the three axes are equal exactly.
     *
     *  @param  lhs One frame.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Frame &lhs, const Frame &rhs) -> bool
    {
      return lhs.X == rhs.X && lhs.Y == rhs.Y && lhs.Z == rhs.Z;
    }
  };

  /**
   *  The frame of a plane by its normal, by AutoCAD's arbitrary axis algorithm.
   *
   *  The normal is made unit. Its x axis is the world y axis crossed with it when both its x and y
   *  are under 1/64, strictly, and the world z axis crossed with it otherwise; its y axis is the
   *  normal crossed with that.
   *
   *  @param  normal The plane's normal, of any length.
   *  @return        The frame, or none for a zero or non-finite normal, which has no plane.
   */
  [[nodiscard]] GEOXY_API auto Facing(const Vector &normal) -> std::optional<Frame>;

  /**
   *  The map from a plane to the plan: where its points fall seen from above, their world z
   *  dropped.
   *
   *  @param  frame     The plane's frame.
   *  @param  elevation How far along its normal the plane lies.
   *  @return           The map from the plane's own coordinates; mirrored for a plane facing down,
   *                    flat for one seen edge on.
   */
  [[nodiscard]] GEOXY_API auto Plan(const Frame &frame, double elevation) -> Geoxy::Matrix;
} // namespace Geoxy::Space

#endif /* GEOXY_SPACE_FRAME_HPP */
