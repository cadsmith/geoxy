/***************************************************************************************************
 *  @file      source/domain.cpp
 *  @brief     The range of a polyline's, a spline's and any curve's parameter.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/domain.hpp"

#include "basis.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"

#include <cstddef>
#include <variant>

namespace Geoxy
{
  auto Domain(const Polyline &polyline) -> Range
  {
    return {.Start = 0.0, .End = static_cast<double>(Count(polyline))};
  }

  auto Domain(const Spline &spline) -> Range
  {
    if (!Counted(spline))
    {
      return {};
    }
    const auto degree = static_cast<std::size_t>(spline.Degree);
    return {.Start = spline.Knots[degree], .End = spline.Knots[spline.Poles.size()]};
  }

  auto Domain(const Curve &curve) -> Range
  {
    return std::visit([](const auto &shape) -> Range { return Domain(shape); }, curve);
  }
} // namespace Geoxy
