/***************************************************************************************************
 *  @file      source/vector.cpp
 *  @brief     Vectors: products, lengths and units, finite at every scale.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/vector.hpp"

#include <cmath>
#include <optional>

namespace Geoxy
{
  auto Unit(Vector vector) -> std::optional<Vector>
  {
    const double largest = Largest(vector);
    if (!std::isfinite(vector.X) || !std::isfinite(vector.Y) || !(largest > 0.0))
    {
      return std::nullopt;
    }
    const Vector scaled = vector / largest;
    return scaled / std::sqrt(Dot(scaled, scaled));
  }
} // namespace Geoxy
