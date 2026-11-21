/***************************************************************************************************
 *  @file      geoxy/region.hpp
 *  @brief     Regions: what loops fill, where a point lies against one, and the union, difference
 *             and intersection of two made of lines and arcs.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_REGION_HPP
#define GEOXY_REGION_HPP

#include "geoxy/accuracy.hpp"
#include "geoxy/api.hpp"
#include "geoxy/error.hpp"
#include "geoxy/loop.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"

#include <cstdint>
#include <expected>
#include <vector>

namespace Geoxy
{
  /**
   *  A region: what its loops fill, every point about which they wind a number of times other than
   *  0, as AutoCAD's REGION and a HATCH fill.
   *
   *  An outer loop runs counterclockwise and a hole clockwise, as the booleans give them; loops
   *  that cross themselves or each other fill all they wind about, whichever way.
   */
  struct Region
  {
    std::vector<Loop> Loops = {}; /**< The loops, outer ones and holes. */

    /**
     *  Whether two regions are equal exactly, loop for loop.
     *
     *  @param  lhs One region.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend auto operator==(const Region &lhs, const Region &rhs) -> bool
    {
      return lhs.Loops == rhs.Loops;
    }
  };

  /**
   *  Where a point lies against a region.
   */
  enum class Location : std::uint8_t
  {
    Outside, /**< Outside it: its loops wind about the point 0 times in all. */
    On,      /**< On it: within the tolerance of a loop.                     */
    Inside,  /**< Inside it.                                                 */
  };

  /**
   *  Where a point lies against a region: on it within the tolerance of any loop, else inside where
   *  its loops wind about it in all a number of times other than 0, as Winding counts them.
   *
   *  @param  region    The region.
   *  @param  point     The point.
   *  @param  tolerance Within which a point is on a loop.
   *  @return           Where; On also for a point or a loop of numbers that are not finite, or a
   *                    loop with an edge without end, of which no side can be told.
   */
  [[nodiscard]] GEOXY_API auto Locate(const Region &region, Point point, const Tolerance &tolerance)
      -> Location;

  /**
   *  The area of a region: the sum of its loops' signed areas, each part counted as many times as
   *  the loops wind about it. For loops that do not cross, an outer one counterclockwise and its
   *  holes clockwise, its area; for any, the area Filled gives is the area counted once.
   *
   *  @param  region   The region.
   *  @param  accuracy How far each edge's area may be off, for a kind that is approximated.
   *  @return          The signed area.
   */
  [[nodiscard]] GEOXY_API auto Area(const Region &region, Accuracy accuracy) -> double;

  /**
   *  What a region fills, as loops that never cross: outer ones counterclockwise, holes clockwise,
   *  each part filled once. A loop that crosses itself, a bow-tie or a loop that winds twice, comes
   *  out as the loops of what it encloses.
   *
   *  The loops' lines and arcs are cut where any two meet, Intersect deciding every crossing and
   *  every stretch two share; vertices within the tolerance merged, in the order of the curves'
   *  ends and then of their crossings, each into the first before it; the faces traced, and the
   *  edges kept that have the region on one side alone. Pieces of one curve that follow each other
   *  are joined back, so an arc kept whole is the same arc.
   *
   *  @param  region    The region, of lines, arcs and polylines.
   *  @param  tolerance Within which two points are one, as Intersect takes it.
   *  @return           The region; Degenerate for a loop that does not close within the tolerance
   *                    or with an edge that is not a line, an arc of up to a turn or a polyline;
   *                    NotFinite for a number that is not finite; Spent where an intersection's
   *                    budget ran out.
   */
  [[nodiscard]] GEOXY_API auto Filled(const Region &region, const Tolerance &tolerance)
      -> std::expected<Region, Error>;

  /**
   *  What either of two regions fills, as UNION gives it.
   *
   *  As Filled works: both regions' lines and arcs cut where any two meet, the faces traced, and
   *  the edges kept with a face of either on one side and of neither on the other.
   *
   *  @param  first     One region, of lines, arcs and polylines.
   *  @param  second    The other.
   *  @param  tolerance Within which two points are one, as Intersect takes it.
   *  @return           The union, as Filled returns.
   */
  [[nodiscard]] GEOXY_API auto Union(const Region &first, const Region &second,
                                     const Tolerance &tolerance) -> std::expected<Region, Error>;

  /**
   *  What one region fills and another does not, as SUBTRACT gives it.
   *
   *  @param  first     The region taken from, of lines, arcs and polylines.
   *  @param  second    The region taken away.
   *  @param  tolerance Within which two points are one, as Intersect takes it.
   *  @return           The difference, as Filled returns.
   */
  [[nodiscard]] GEOXY_API auto Difference(const Region &first, const Region &second,
                                          const Tolerance &tolerance)
      -> std::expected<Region, Error>;

  /**
   *  What both of two regions fill, as INTERSECT gives it.
   *
   *  @param  first     One region, of lines, arcs and polylines.
   *  @param  second    The other.
   *  @param  tolerance Within which two points are one, as Intersect takes it.
   *  @return           The intersection, as Filled returns.
   */
  [[nodiscard]] GEOXY_API auto Intersection(const Region &first, const Region &second,
                                            const Tolerance &tolerance)
      -> std::expected<Region, Error>;
} // namespace Geoxy

#endif /* GEOXY_REGION_HPP */
