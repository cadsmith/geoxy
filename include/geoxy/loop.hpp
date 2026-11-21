/***************************************************************************************************
 *  @file      geoxy/loop.hpp
 *  @brief     Loops: closed chains of curves, how many times one winds about a point, the area it
 *             encloses and its perimeter.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_LOOP_HPP
#define GEOXY_LOOP_HPP

#include "geoxy/accuracy.hpp"
#include "geoxy/api.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"

#include <optional>
#include <vector>

namespace Geoxy
{
  /**
   *  A loop: curves in a chain, each ending where the next starts and the last where the first
   *  starts, as a HATCH boundary path or a region's edge.
   *
   *  Any finite curve may be an edge: a line, an arc, a polyline, an ellipse, a Bézier curve, a
   *  spline; a closed polyline or a circle alone is a loop. Where an edge's end is not the next
   *  one's start, the gap between them is closed by the line from the one to the other, as AREA
   *  closes an open curve. A loop may cross itself and run about a point more than once.
   */
  struct Loop
  {
    std::vector<Curve> Edges = {}; /**< The edges, in order along it. */

    /**
     *  Whether two loops are equal exactly, edge for edge.
     *
     *  @param  lhs One loop.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend auto operator==(const Loop &lhs, const Loop &rhs) -> bool
    {
      return lhs.Edges == rhs.Edges;
    }
  };

  /**
   *  How many times a loop runs about a point, counterclockwise.
   *
   *  A ray is cast from the point and every edge counted where it crosses it, upward +1 and
   *  downward -1, by tests decided exactly on the numbers stored: a line by the side of it the
   *  point lies on, an arc cut into parts that rise or fall all along, each by the side of its
   *  center the point lies on and by whether the point is inside its circle, worked out in sums of
   *  doubles that never round; any other curve cut where a part is longer than half its distance
   *  from the point, so that no part can wind about it, and each part counted as the line between
   *  its ends. Each part is counted from the very ends the next one starts from, so a ray through a
   *  vertex counts once, and the count is the winding number of a closed curve within the tolerance
   *  of the loop.
   *
   *  @param  loop      The loop.
   *  @param  point     The point.
   *  @param  tolerance Within which a point is on the loop.
   *  @return           The winding number; none where the point is within the tolerance of an edge
   *                    or of a gap's line, where the loop has an edge without end, or where a
   *                    number of either is not finite.
   */
  [[nodiscard]] GEOXY_API auto Winding(const Loop &loop, Point point, const Tolerance &tolerance)
      -> std::optional<int>;

  /**
   *  The area a loop encloses, each part as many times as the loop winds about it.
   *
   *  The area of each edge with its chord, as Area gives it, and the area of the polygon of the
   *  edges' ends and starts, summed by the shoelace formula carried with its rounding error: a gap
   *  is closed by its line, and a loop that crosses itself counts each lobe by its winding, one
   *  that runs the other way negative.
   *
   *  @param  loop     The loop.
   *  @param  accuracy How far each edge's area may be off, for a kind that is approximated.
   *  @return          The signed area, positive for a loop that runs counterclockwise; 0 for one
   *                   with an edge without end.
   */
  [[nodiscard]] GEOXY_API auto Area(const Loop &loop, Accuracy accuracy) -> double;

  /**
   *  The perimeter of a loop: the length of every edge, and of every gap's line.
   *
   *  @param  loop     The loop.
   *  @param  accuracy How far each edge's length may be off, for a kind that is approximated.
   *  @return          The length; infinite for a loop with an edge without end.
   */
  [[nodiscard]] GEOXY_API auto Length(const Loop &loop, Accuracy accuracy) -> double;

  /**
   *  Whether a loop closes: each edge's end within a tolerance of the next one's start, the last's
   *  of the first's, and every edge of finite numbers and with an end.
   *
   *  @param  loop      The loop.
   *  @param  tolerance Within which an end and a start are one.
   *  @return           Whether it does; false for a loop of no edge.
   */
  [[nodiscard]] GEOXY_API auto Closed(const Loop &loop, const Tolerance &tolerance) -> bool;
} // namespace Geoxy

#endif /* GEOXY_LOOP_HPP */
