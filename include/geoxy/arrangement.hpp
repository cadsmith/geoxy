/***************************************************************************************************
 *  @file      geoxy/arrangement.hpp
 *  @brief     Curves cut where they meet into the faces they bound, once, and the loop around a
 *             picked point among them, with its islands, as HATCH and BOUNDARY find it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_ARRANGEMENT_HPP
#define GEOXY_ARRANGEMENT_HPP

#include "geoxy/api.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/error.hpp"
#include "geoxy/point.hpp"
#include "geoxy/region.hpp"
#include "geoxy/tolerance.hpp"

#include <expected>
#include <memory>
#include <span>

namespace Geoxy
{
  struct Arranged;

  /**
   *  How curves are arranged: within which two points are one, and how wide a gap is closed,
   *  between two loose ends or a loose end and a curve.
   */
  struct Closing
  {
    Tolerance Within = {};  /**< Within which two points are one, as Intersect takes it. */
    double    Gap    = 0.0; /**< The widest gap closed, as HPGAPTOL; 0 closes none.      */
  };

  /**
   *  Curves cut where they meet, and the faces they bound: what every pick of a HATCH or a BOUNDARY
   *  reads, worked out once.
   *
   *  Every crossing comes from a Grid of the curves, intersected all at once, each placed on both
   *  curves where it is nearest them, as Intersect places a hit; where two curves share a stretch,
   *  every cut of one in it is made in the other, and the two pieces are one edge. Vertices within
   *  the tolerance are merged in a stated order, the curves' ends in the order given and then the
   *  crossings as they are met, each into the first before it. Edges leave each vertex in the exact
   *  order of their tangents, those within the tolerance's angle by their curvature; the faces are
   *  traced from that order, with every edge with a loose end taken away. A gap up to the closing's
   *  between two loose ends, nearest first, or between a loose end and the nearest curve, is
   *  bridged by a line first.
   *
   *  The faces are kept, not changed: a copy shares them.
   */
  class Arrangement
  {
  public:
    /**
     *  The arrangement of curves.
     *
     *  @param curves  The curves: lines, arcs of up to a turn, and polylines.
     *  @param closing The tolerance, and the widest gap to close.
     */
    GEOXY_API Arrangement(std::span<const Curve> curves, const Closing &closing);

  private:
    friend struct Arranging;

    std::shared_ptr<const Arranged> m_Faces = {}; /**< The faces, and what they were made from. */
  };

  /**
   *  The loop around a point among the curves of an arrangement, and its islands: the face the
   *  point lies in, as HATCH and BOUNDARY pick it.
   *
   *  The face is the one bounded by the cycle of least area, among those that run counterclockwise
   *  and wind about the point, as Winding counts; its islands are the outer cycles of the curves
   *  that lie in it untouched by it, each nearest enclosed by that same cycle. The loop runs
   *  counterclockwise and the islands clockwise; pieces of one curve that follow each other are
   *  joined back.
   *
   *  @param  arrangement The arrangement.
   *  @param  point       The point picked.
   *  @return             The face, as a region of its outer loop and its islands; a region of no
   *                      loop where no loop is around the point, or the point is on a curve within
   *                      the tolerance. Degenerate where a curve is a ray, an endless line, an arc
   *                      past a turn or another kind not arranged yet; NotFinite where a number of
   *                      the curves or the point is not finite; Spent where an intersection's
   *                      budget ran out.
   */
  [[nodiscard]] GEOXY_API auto Boundary(const Arrangement &arrangement, Point point)
      -> std::expected<Region, Error>;
} // namespace Geoxy

#endif /* GEOXY_ARRANGEMENT_HPP */
