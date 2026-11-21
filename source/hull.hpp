/***************************************************************************************************
 *  @file      source/hull.hpp
 *  @brief     A box a polyline's segment cannot leave, worked out from its vertices and its bulge
 *             alone, and how far a point lies from a box.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_HULL_HPP
#define GEOXY_HULL_HPP

#include "geoxy/box.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>

namespace Geoxy
{
  /**
   *  A box no point worked out along a polyline's segment lies outside: its chord's, grown by how
   *  far its arc bows out, the bulge times half the chord up to a half turn and the diameter past
   *  it, the chord taken no shorter than it is, as the sum of its two extents; and grown again by
   *  far more than a point's rounding, the steps of a double at the size of its ends and its reach,
   *  and past a half turn the sine of half the sweep, which is lost in its own rounding as the arc
   *  closes.
   *
   *  @param  polyline The polyline.
   *  @param  index    The segment.
   *  @return          The box.
   */
  inline auto Hull(const Polyline &polyline, std::size_t index) -> Box
  {
    const std::span<const Vertex> vertices = polyline.Vertices;
    const Point                   start    = vertices[index].Position;
    const Point                   end      = vertices[(index + 1) % vertices.size()].Position;
    const double                  bulge    = std::abs(vertices[index].Bulge);
    const double                  reach    = std::abs(end.X - start.X) + std::abs(end.Y - start.Y);
    double                        height   = 0.0;
    double                        spread   = 0.0;
    if (bulge >= Straight && bulge <= 1.0)
    {
      height = bulge * reach / 2.0;
    }
    else if (bulge > 1.0)
    {
      height = reach * (1.0 + (bulge * bulge)) / (2.0 * bulge);
      spread = height * bulge * 0x1p-44;
    }
    const double size =
        std::max({std::abs(start.X), std::abs(start.Y), std::abs(end.X), std::abs(end.Y)})
        + reach
        + (2.0 * height);
    const double grow = (height * (1.0 + 0x1p-20)) + spread + (size * 0x1p-36);
    return {.Min = {.X = std::min(start.X, end.X) - grow, .Y = std::min(start.Y, end.Y) - grow},
            .Max = {.X = std::max(start.X, end.X) + grow, .Y = std::max(start.Y, end.Y) + grow}};
  }

  /**
   *  How far a point lies past a box along x and along y: nothing within it.
   *
   *  @param  box   The box.
   *  @param  point The point.
   *  @return       The gap along each axis.
   */
  inline auto Overhang(const Box &box, Point point) -> Vector
  {
    return {.X = std::max({box.Min.X - point.X, point.X - box.Max.X, 0.0}),
            .Y = std::max({box.Min.Y - point.Y, point.Y - box.Max.Y, 0.0})};
  }

  /**
   *  How far a point lies from a box: nothing within it.
   *
   *  @param  box   The box.
   *  @param  point The point.
   *  @return       The distance, no more than from the point to any point of the box, as Distance
   *                works it out.
   */
  inline auto Outlying(const Box &box, Point point) -> double
  {
    return Length(Overhang(box, point));
  }
} // namespace Geoxy

#endif /* GEOXY_HULL_HPP */
