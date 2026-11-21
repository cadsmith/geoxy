/***************************************************************************************************
 *  @file      source/space/frame.cpp
 *  @brief     Frames by AutoCAD's arbitrary axis, and their plans.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/space/frame.hpp"

#include "geoxy/matrix.hpp"
#include "geoxy/space/vector.hpp"

#include <cmath>
#include <optional>

namespace Geoxy::Space
{
  namespace
  {
    constexpr double Arbitrary = 1.0 / 64.0;

    /**
     *  A vector made unit that is known not to be zero.
     */
    auto Normalized(const Vector &vector) -> Vector
    {
      const double length = Length(vector);
      return {.X = vector.X / length, .Y = vector.Y / length, .Z = vector.Z / length};
    }
  } // namespace

  auto Facing(const Vector &normal) -> std::optional<Frame>
  {
    const std::optional<Vector> z = Unit(normal);
    if (!z)
    {
      return std::nullopt;
    }
    const bool   nearPole = std::abs(z->X) < Arbitrary && std::abs(z->Y) < Arbitrary;
    const Vector x        = Normalized(Cross(nearPole ? North : Up, *z));
    const Vector y        = Normalized(Cross(*z, x));
    return Frame{.X = x, .Y = y, .Z = *z};
  }

  auto Plan(const Frame &frame, double elevation) -> Geoxy::Matrix
  {
    return {.X      = {.X = frame.X.X, .Y = frame.X.Y},
            .Y      = {.X = frame.Y.X, .Y = frame.Y.Y},
            .Origin = {.X = frame.Z.X * elevation, .Y = frame.Z.Y * elevation}};
  }
} // namespace Geoxy::Space
