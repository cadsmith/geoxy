/***************************************************************************************************
 *  @file      source/tolerance.cpp
 *  @brief     Tolerances of a model, and the step of a double.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/tolerance.hpp"

#include "real.hpp"

#include "geoxy/box.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Geoxy
{
  namespace
  {
    constexpr double Relative = 1e-10;
    constexpr double Steps    = 16.0;
  } // namespace

  auto Tolerance::Of(const Box &extents) -> Tolerance
  {
    if (IsEmpty(extents))
    {
      return {.Length = Steps * Grain(0.0), .Radians = Relative};
    }
    const Vector size = Size(extents);
    /**
     *  A side or a corner past every double makes a length past them too, or no number at all: the
     *  length is the largest double, and the angle the ratio alone.
     */
    if (!std::isfinite(size.X) || !std::isfinite(size.Y))
    {
      return {.Length = std::numeric_limits<double>::max(), .Radians = Relative};
    }
    const double span   = std::max(size.X, size.Y);
    const double far    = std::max(Grain(extents.Min), Grain(extents.Max));
    const double length = std::max(Relative * span, Steps * far);
    return {.Length = length, .Radians = span > 0.0 ? length / span : Relative};
  }

  auto Same(Point first, Point second, const Tolerance &tolerance) -> bool
  {
    return Distance(first, second) <= tolerance.Length;
  }

  auto Grain(double value) -> double
  {
    return Above(std::abs(value));
  }

  auto Grain(Point point) -> double
  {
    return std::max(Grain(point.X), Grain(point.Y));
  }
} // namespace Geoxy
