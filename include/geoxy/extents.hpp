/***************************************************************************************************
 *  @file      geoxy/extents.hpp
 *  @brief     The box of a curve: the smallest axis-aligned box that holds it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_EXTENTS_HPP
#define GEOXY_EXTENTS_HPP

#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

namespace Geoxy
{
  /**
   *  The box of a line: the smallest that holds it, its corners found where it reaches them, not
   *  sampled.
   *
   *  @param  line The line.
   *  @return      The box of its ends.
   */
  [[nodiscard]] GEOXY_API auto Extents(const Line &line) -> Box;

  /**
   *  The box of a ray: the smallest that holds it, its corners found where it reaches them, not
   *  sampled.
   *
   *  @param  ray The ray.
   *  @return     Its start's box, run without end along each axis its way leaves: infinite there.
   */
  [[nodiscard]] GEOXY_API auto Extents(const Ray &ray) -> Box;

  /**
   *  The box of an xline: the smallest that holds it, its corners found where it reaches them, not
   *  sampled.
   *
   *  @param  xline The xline.
   *  @return       Its point's box, run without end both ways along each axis its way leaves.
   */
  [[nodiscard]] GEOXY_API auto Extents(const Xline &xline) -> Box;

  /**
   *  The box of an arc: the smallest that holds it and every point At gives on it, its corners
   *  found where it reaches them, not sampled.
   *
   *  @param  arc The arc.
   *  @return     The box of every point a quarter turn round its circle that it passes, exact there
   *              or the next double out, and of its ends, each padded by its rounding but at a
   *              quarter turn.
   */
  [[nodiscard]] GEOXY_API auto Extents(const Arc &arc) -> Box;

  /**
   *  The box of an ellipse: the smallest that holds it and every point At gives on it, within a few
   *  steps of a double, its corners found where it reaches them, not sampled.
   *
   *  @param  ellipse The ellipse.
   *  @return         The box of each extreme of a coordinate it passes, the center's plus or less
   *                  the length of the map's row, and of its ends, each padded by its rounding.
   */
  [[nodiscard]] GEOXY_API auto Extents(const Ellipse &ellipse) -> Box;

  /**
   *  The box of a polyline: the smallest that holds it and every point At gives on it, within a few
   *  steps of a double, its corners found where it reaches them, not sampled.
   *
   *  @param  polyline The polyline.
   *  @return          The box of its vertices and of each bulged segment's extremes, worked out
   *                   from its chord and its bulge with no trigonometry, a bulged segment's padded
   *                   by four steps of a double of its reach; empty for no vertices.
   */
  [[nodiscard]] GEOXY_API auto Extents(const Polyline &polyline) -> Box;

  /**
   *  The box of a quadratic Bézier curve: one that holds it and every point At gives on it, its
   *  corners found where it reaches them, not sampled; exact at its ends, and a few steps of a
   *  double beyond where it turns, as At's own rounding may reach.
   *
   *  @param  quadratic The curve.
   *  @return           The box of its ends and of where its derivative is zero along an axis.
   */
  [[nodiscard]] GEOXY_API auto Extents(const Bezier::Quadratic &quadratic) -> Box;

  /**
   *  The box of a cubic Bézier curve: one that holds it and every point At gives on it, its corners
   *  found where it reaches them, not sampled; exact at its ends, and a few steps of a double
   *  beyond where it turns, as At's own rounding may reach.
   *
   *  @param  cubic The curve.
   *  @return       The box of its ends and of where its derivative is zero along an axis, by the
   *                stable quadratic formula.
   */
  [[nodiscard]] GEOXY_API auto Extents(const Bezier::Cubic &cubic) -> Box;

  /**
   *  The box of a spline: the smallest that holds it and every point At gives on it, within a few
   *  steps of a double of its poles.
   *
   *  Its spans' ends, and where a piece's coordinate turns, its derivative's numerator zero, are
   *  worked out as At works them out, each with a bound on its rounding; each is then padded by
   *  that bound and by how far At strays from the spline on its span, (4 d + 1) steps of a double
   *  of its poles' largest coordinate. A piece whose control points, so padded, lie within the box
   *  of the ends adds nothing.
   *
   *  @param  spline The spline.
   *  @return        The box; the empty one for a spline that is not valid.
   */
  [[nodiscard]] GEOXY_API auto Extents(const Spline &spline) -> Box;

  /**
   *  The box of any curve.
   *
   *  @param  curve The curve.
   *  @return       The box.
   */
  [[nodiscard]] GEOXY_API auto Extents(const Curve &curve) -> Box;
} // namespace Geoxy

#endif /* GEOXY_EXTENTS_HPP */
