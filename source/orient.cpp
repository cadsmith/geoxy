/***************************************************************************************************
 *  @file      source/orient.cpp
 *  @brief     Orientation, decided by the exact predicate.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/orient.hpp"

#include "exact/orient.hpp"

#include "geoxy/point.hpp"

namespace Geoxy
{
  /**
   *  The enumerators are the signs, so the exact sign is the orientation.
   */
  auto Orient(Point first, Point second, Point third) -> Orientation
  {
    return static_cast<Orientation>(Exact::Orient(first, second, third));
  }
} // namespace Geoxy
