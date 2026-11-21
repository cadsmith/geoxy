/***************************************************************************************************
 *  @file      source/segment.cpp
 *  @brief     A polyline's segments.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/segment.hpp"

#include "bow.hpp"
#include "trig.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"

#include <cmath>
#include <cstddef>

namespace Geoxy
{
  auto Bent(Point start, Point end, double bulge) -> Segment
  {
    if (!(std::abs(bulge) >= Straight) || start == end)
    {
      return Line{.Start = start, .End = end};
    }
    const Hoop hoop = Hooped(start, end, bulge);
    return Arc{
        .Center = hoop.Center,
        .Radius = hoop.Radius,
        .Start  = Heading(start - hoop.Center),
        .Sweep  = Radians(4.0 * Trig::Atan(bulge)),
    };
  }

  auto Count(const Polyline &polyline) -> std::size_t
  {
    const std::size_t vertices = polyline.Vertices.size();
    if (vertices < 2)
    {
      return 0;
    }
    return polyline.Closed ? vertices : vertices - 1;
  }

  auto Nth(const Polyline &polyline, std::size_t index) -> Segment
  {
    const Vertex &start = polyline.Vertices[index];
    const Vertex &end   = polyline.Vertices[(index + 1) % polyline.Vertices.size()];
    return Bent(start.Position, end.Position, start.Bulge);
  }
} // namespace Geoxy
