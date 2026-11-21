/***************************************************************************************************
 *  @file      source/box.cpp
 *  @brief     Boxes: their middles and their margins.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/box.hpp"

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

namespace Geoxy
{
  auto Center(const Box &box) -> Point
  {
    return Midpoint(box.Min, box.Max);
  }

  auto Padded(const Box &box, double margin) -> Box
  {
    const Vector pad = {.X = margin, .Y = margin};
    return {.Min = box.Min - pad, .Max = box.Max + pad};
  }
} // namespace Geoxy
