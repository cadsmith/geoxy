/***************************************************************************************************
 *  @file      geoxy/contain.hpp
 *  @brief     What a window picks, and what a crossing window picks: whether a curve lies wholly
 *             within a box or a polygon, and whether it meets one.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CONTAIN_HPP
#define GEOXY_CONTAIN_HPP

#include "geoxy/api.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/region.hpp"
#include "geoxy/tolerance.hpp"

namespace Geoxy
{
  /**
   *  Whether a curve lies wholly within a box, its edges included: what a window picks.
   *
   *  Told from the curve's extents, which hold every point At gives on it: a curve they put inside
   *  is inside. They are exact for a line, a ray, an endless line and a polyline of no bulge; for
   *  any other curve they may lie a few steps of a double past where it reaches, so a curve whose
   *  farthest point is within those few steps of an edge, inside, may be told it is not.
   *
   *  @param  curve The curve.
   *  @param  box   The box.
   *  @return       Whether it lies within; never for a curve of no point, nor one of a number that
   *                is not a number.
   */
  [[nodiscard]] GEOXY_API auto Inside(const Curve &curve, const Box &box) -> bool;

  /**
   *  Whether a curve meets a box, its edges included: has a point within it or on an edge, as a
   *  crossing window picks.
   *
   *  A line is told exactly, by the box of its ends and the side of it each corner of the box lies
   *  on, as Orient tells them: a line that only touches a corner meets the box. Any other curve
   *  meets it where it lies within, as Inside tells it, or where Clip cuts a stretch of it longer
   *  than nothing out of the box; one that only touches the box at one point may be told it does
   *  not.
   *
   *  @param  curve The curve.
   *  @param  box   The box.
   *  @return       Whether it meets it; never an empty box.
   */
  [[nodiscard]] GEOXY_API auto Intersects(const Curve &curve, const Box &box) -> bool;

  /**
   *  Whether a curve lies wholly within a region, its edges included: what WPOLYGON picks, the
   *  region a polygon's.
   *
   *  The curve is cut wherever it meets an edge of the region, as Intersect finds it within the
   *  tolerance, crossings, touches and both ends of every stretch the two share; it lies within
   *  where the middle of every piece does, as Locate tells it, on an edge counting as within.
   *
   *  @param  curve     The curve.
   *  @param  region    The region.
   *  @param  tolerance Within which the curve meets an edge, and a point is on one.
   *  @return           Whether it lies within; never for a ray, an endless line, a curve of no
   *                    point or one of a number that is not finite.
   */
  [[nodiscard]] GEOXY_API auto Inside(const Curve &curve, const Region &region,
                                      const Tolerance &tolerance) -> bool;

  /**
   *  Whether a curve meets a region, its edges included: what CPOLYGON picks, the region a
   *  polygon's.
   *
   *  It does where it meets an edge, as Intersect finds it within the tolerance, or else where its
   *  start lies within, as Locate tells it: a curve that meets no edge lies wholly within or wholly
   *  without.
   *
   *  @param  curve     The curve.
   *  @param  region    The region.
   *  @param  tolerance Within which the curve meets an edge, and a point is on one.
   *  @return           Whether it meets it; never a curve of a number that is not finite.
   */
  [[nodiscard]] GEOXY_API auto Intersects(const Curve &curve, const Region &region,
                                          const Tolerance &tolerance) -> bool;

  /**
   *  Whether a curve lies wholly within a polygon, its edges included: what WPOLYGON picks.
   *
   *  @param  curve     The curve.
   *  @param  polygon   The polygon, its segments lines or arcs; open, it is closed by the line from
   *                    its last vertex to its first.
   *  @param  tolerance Within which the curve meets an edge, and a point is on one.
   *  @return           As Inside of the region of the polygon's one loop.
   */
  [[nodiscard]] GEOXY_API auto Inside(const Curve &curve, const Polyline &polygon,
                                      const Tolerance &tolerance) -> bool;

  /**
   *  Whether a curve meets a polygon, its edges included: what CPOLYGON picks.
   *
   *  @param  curve     The curve.
   *  @param  polygon   The polygon, closed as Inside closes it.
   *  @param  tolerance Within which the curve meets an edge, and a point is on one.
   *  @return           As Intersects of the region of the polygon's one loop.
   */
  [[nodiscard]] GEOXY_API auto Intersects(const Curve &curve, const Polyline &polygon,
                                          const Tolerance &tolerance) -> bool;
} // namespace Geoxy

#endif /* GEOXY_CONTAIN_HPP */
