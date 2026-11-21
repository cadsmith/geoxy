/***************************************************************************************************
 *  @file      source/ends.cpp
 *  @brief     Where a curve starts and where it ends.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/ends.hpp"

#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/place.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/spline.hpp"

#include <cmath>
#include <optional>
#include <variant>

namespace Geoxy
{
  auto Ends(const Curve &curve) -> Endpoints
  {
    const auto *polyline = std::get_if<Polyline>(&curve);
    const auto *spline   = std::get_if<Spline>(&curve);
    if ((polyline != nullptr && polyline->Vertices.empty())
        || (spline != nullptr && !Valid(*spline)))
    {
      return {};
    }
    const Range domain = Domain(curve);
    const auto  end    = [&curve](double along) -> std::optional<Place>
    {
      if (std::isinf(along))
      {
        return std::nullopt;
      }
      return Place{.At = At(curve, along), .Along = along};
    };
    return {.Start = end(domain.Start), .End = end(domain.End)};
  }
} // namespace Geoxy
