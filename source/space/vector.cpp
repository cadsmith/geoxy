/***************************************************************************************************
 *  @file      source/space/vector.cpp
 *  @brief     Vectors in space: products, lengths and units.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/space/vector.hpp"

#include <algorithm>
#include <cmath>
#include <optional>

namespace Geoxy::Space
{
  namespace
  {
    constexpr double Low  = 0x1.0p-500;
    constexpr double High = 0x1.0p500;

    auto Largest(const Vector &vector) -> double
    {
      return std::max({std::abs(vector.X), std::abs(vector.Y), std::abs(vector.Z)});
    }

    auto Scaled(const Vector &vector, double divisor) -> Vector
    {
      return {.X = vector.X / divisor, .Y = vector.Y / divisor, .Z = vector.Z / divisor};
    }
  } // namespace

  auto Dot(const Vector &first, const Vector &second) -> double
  {
    return (first.X * second.X) + (first.Y * second.Y) + (first.Z * second.Z);
  }

  auto Cross(const Vector &first, const Vector &second) -> Vector
  {
    return {.X = (first.Y * second.Z) - (first.Z * second.Y),
            .Y = (first.Z * second.X) - (first.X * second.Z),
            .Z = (first.X * second.Y) - (first.Y * second.X)};
  }

  auto Length(const Vector &vector) -> double
  {
    const double largest = Largest(vector);
    if (largest >= Low && largest <= High)
    {
      return std::sqrt(Dot(vector, vector));
    }
    if (!(largest > 0.0) || !std::isfinite(largest))
    {
      return largest;
    }
    const Vector scaled = Scaled(vector, largest);
    return largest * std::sqrt(Dot(scaled, scaled));
  }

  auto Unit(const Vector &vector) -> std::optional<Vector>
  {
    const double largest = Largest(vector);
    if (!std::isfinite(vector.X)
        || !std::isfinite(vector.Y)
        || !std::isfinite(vector.Z)
        || !(largest > 0.0))
    {
      return std::nullopt;
    }
    const Vector scaled = Scaled(vector, largest);
    return Scaled(scaled, std::sqrt(Dot(scaled, scaled)));
  }
} // namespace Geoxy::Space
