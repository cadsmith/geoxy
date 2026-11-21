/***************************************************************************************************
 *  @file      geoxy/triangulate.hpp
 *  @brief     Fills cut into triangles: a region given as an outline and its holes, by ears clipped
 *             with exact turns, and made Delaunay if asked.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TRIANGULATE_HPP
#define GEOXY_TRIANGULATE_HPP

#include "geoxy/api.hpp"
#include "geoxy/point.hpp"

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace Geoxy
{
  /**
   *  How the triangles of a fill are shaped.
   */
  enum class Fill : std::uint8_t
  {
    Ears,     /**< As the ears are clipped: fastest, slivers and all. */
    Delaunay, /**< Then flipped to the constrained Delaunay triangulation: no slivers but
                   what the outline forces.                                              */
  };

  /**
   *  A fill's triangles, and whether they cover all of it.
   */
  struct Mesh
  {
    std::vector<std::array<std::uint32_t, 3>> Triangles = {}; /**< Each triangle's corners, as
                                                                   indices into the loops' points
                                                                   taken in order, turning
                                                                   counterclockwise.            */
    bool Whole = true; /**< Whether they cover the region: false where loops cross, which leaves
                            the part no ear could be cut from without a triangle.                */
  };

  /**
   *  The triangles of a region: its outline and its holes, as earcut takes them.
   *
   *  Each loop is a polygon, either way round, its last point joined to its first. Points one after
   *  the other the same, or on a line with their neighbours, are no corners. Each hole, its
   *  leftmost point first, is bridged to the outline by the vertex nearest the ray to its left that
   *  no other corner hides; then ears are cut where a corner turns counterclockwise and no corner
   *  that turns the other way lies within or on the triangle, every turn decided by Orient,
   *  exactly. Large regions look for such corners by their order along a Z curve, as earcut does.
   *
   *  The triangles' areas sum to the region's exactly, as only ears and zero areas are cut away,
   *  and for loops that neither cross nor touch, every triangle is within the region.
   *
   *  @param  loops The outline first, then the holes; holes within the outline.
   *  @param  fill  How the triangles are shaped.
   *  @return       The triangles; none for an outline of no area.
   */
  [[nodiscard]] GEOXY_API auto Triangulate(std::span<const std::vector<Point>> loops, Fill fill)
      -> Mesh;
} // namespace Geoxy

#endif /* GEOXY_TRIANGULATE_HPP */
