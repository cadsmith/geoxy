/***************************************************************************************************
 *  @file      geoxy/polyline.hpp
 *  @brief     Polylines: vertices joined by lines and arcs, an arc by its bulge.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_POLYLINE_HPP
#define GEOXY_POLYLINE_HPP

#include "geoxy/point.hpp"

#include <functional>
#include <vector>

namespace Geoxy
{
  /**
   *  A polyline's vertex, and how the segment from it to the next one bends.
   */
  struct Vertex
  {
    Point  Position = {}; /**< Where it is. */
    double Bulge =
        0.0; /**< The tangent of a quarter of the sweep to the next vertex; 0 is a line. */

    /**
     *  Whether two vertices are equal exactly.
     *
     *  @param  lhs One vertex.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Vertex &lhs, const Vertex &rhs) -> bool
    {
      return lhs.Position == rhs.Position && std::equal_to<double>{}(lhs.Bulge, rhs.Bulge);
    }
  };

  /**
   *  A polyline: vertices joined in order, each segment a line or an arc by its bulge, as DXF's
   *  LWPOLYLINE.
   *
   *  A positive bulge bends counterclockwise. An open polyline's last bulge joins nothing; a closed
   *  one's joins the last vertex to the first.
   */
  struct Polyline
  {
    std::vector<Vertex> Vertices = {};    /**< The vertices, in order.                  */
    bool                Closed   = false; /**< Whether the last vertex joins the first. */

    /**
     *  Whether two polylines are equal exactly, vertex for vertex.
     *
     *  @param  lhs One polyline.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend auto operator==(const Polyline &lhs, const Polyline &rhs) -> bool
    {
      return lhs.Closed == rhs.Closed && lhs.Vertices == rhs.Vertices;
    }
  };
} // namespace Geoxy

#endif /* GEOXY_POLYLINE_HPP */
