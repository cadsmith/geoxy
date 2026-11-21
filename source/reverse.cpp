/***************************************************************************************************
 *  @file      source/reverse.cpp
 *  @brief     Each curve run the other way.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/reverse.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <type_traits>
#include <variant>

namespace Geoxy
{
  auto Reversed(const Line &line) -> Line
  {
    return {.Start = line.End, .End = line.Start};
  }

  auto Reversed(const Xline &xline) -> Xline
  {
    return {.Through = xline.Through, .Way = -xline.Way};
  }

  auto Reversed(const Arc &arc) -> Arc
  {
    return {.Center = arc.Center,
            .Radius = arc.Radius,
            .Start  = arc.Start + (arc.Sweep * 1.0),
            .Sweep  = -arc.Sweep};
  }

  auto Reversed(const Ellipse &ellipse) -> Ellipse
  {
    return {.Map   = ellipse.Map,
            .Start = ellipse.Start + (ellipse.Sweep * 1.0),
            .Sweep = -ellipse.Sweep};
  }

  auto Reversed(const Polyline &polyline) -> Polyline
  {
    /**
     *  The vertex the reversed one at an index is: from the last back to the first, a closed
     *  polyline's first kept first. Its bulge is that of the segment it now leaves, the one that
     *  arrived at it, held by the vertex before; an open polyline's last, which joins nothing, is
     *  the first's own, kept so that reversing twice gives back every bulge.
     */
    const std::size_t count = polyline.Vertices.size();
    const std::size_t shift = polyline.Closed ? 0 : 1;
    Polyline          reversed{.Vertices = {}, .Closed = polyline.Closed};
    reversed.Vertices.reserve(count);
    for (std::size_t index = 0; index < count; ++index)
    {
      const Vertex &placed = polyline.Vertices[((2 * count) - shift - index) % count];
      const Vertex &bent   = polyline.Vertices[((2 * count) - 1 - shift - index) % count];
      reversed.Vertices.push_back({.Position = placed.Position, .Bulge = -bent.Bulge});
    }
    return reversed;
  }

  auto Reversed(const Spline &spline) -> Spline
  {
    Spline reversed = spline;
    std::ranges::reverse(reversed.Poles);
    std::ranges::reverse(reversed.Weights);
    std::ranges::reverse(reversed.Knots);
    std::ranges::transform(reversed.Knots, reversed.Knots.begin(),
                           [](double knot) -> double { return -knot; });
    return reversed;
  }

  auto Reversed(const Bezier::Quadratic &quadratic) -> Bezier::Quadratic
  {
    return {.Start = quadratic.End, .Control = quadratic.Control, .End = quadratic.Start};
  }

  auto Reversed(const Bezier::Cubic &cubic) -> Bezier::Cubic
  {
    return {.Start = cubic.End, .Out = cubic.In, .In = cubic.Out, .End = cubic.Start};
  }

  auto Reversed(const Curve &curve) -> std::optional<Curve>
  {
    return std::visit(
        [](const auto &shape) -> std::optional<Curve>
        {
          if constexpr (std::is_same_v<std::decay_t<decltype(shape)>, Ray>)
          {
            return std::nullopt;
          }
          else
          {
            return Curve{Reversed(shape)};
          }
        },
        curve);
  }
} // namespace Geoxy
