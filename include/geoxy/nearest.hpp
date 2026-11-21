/***************************************************************************************************
 *  @file      geoxy/nearest.hpp
 *  @brief     A curve's places nearest a point: every foot of a perpendicular from it, and the
 *             nearest of those and the curve's ends.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_NEAREST_HPP
#define GEOXY_NEAREST_HPP

#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <vector>

namespace Geoxy
{
  /**
   *  The foot of the perpendicular from a point to a line.
   *
   *  @param  line  The line.
   *  @param  point The point.
   *  @return       The foot, if it lies between the line's ends; none for a line of no length.
   */
  [[nodiscard]] GEOXY_API auto Perpendiculars(const Line &line, Point point) -> std::vector<Place>;

  /**
   *  The foot of the perpendicular from a point to a ray.
   *
   *  @param  ray   The ray.
   *  @param  point The point.
   *  @return       The foot, if it lies past the ray's start; none for a ray of no way.
   */
  [[nodiscard]] GEOXY_API auto Perpendiculars(const Ray &ray, Point point) -> std::vector<Place>;

  /**
   *  The foot of the perpendicular from a point to an xline.
   *
   *  @param  xline The xline.
   *  @param  point The point.
   *  @return       The foot; none for an xline of no way.
   */
  [[nodiscard]] GEOXY_API auto Perpendiculars(const Xline &xline, Point point)
      -> std::vector<Place>;

  /**
   *  The feet of the perpendiculars from a point to an arc, in closed form: where the arc meets the
   *  line through its center and the point.
   *
   *  An arc swept round more than 65,536 times has its feet of its first 65,536 turns only, as Clip
   *  seeks its crossings: past that, too many to list. Swept so far that one turn is less than a
   *  few steps of a double of its parameter, past some 2^46 turns, the feet of turns so near are
   *  given once.
   *
   *  @param  arc   The arc.
   *  @param  point The point.
   *  @return       The feet, in order, every turn of the arc's to its 65,536th, each at the double
   *                nearest the circle's exact point on that line, within a few steps of a double of
   *                At at its parameter; none when the point is the center, where every point of the
   *                arc is as far.
   */
  [[nodiscard]] GEOXY_API auto Perpendiculars(const Arc &arc, Point point) -> std::vector<Place>;

  /**
   *  The feet of the perpendiculars from a point to an ellipse: up to four in a turn, found on its
   *  quarters as exact conics.
   *
   *  An ellipse thin past what doubles hold leaves rounding noise among them, hundreds in a turn:
   *  the four nearest square to the way to the point are kept. One swept round more than 65,536
   *  times has its feet of its first 65,536 turns only, as an arc has.
   *
   *  @param  ellipse The ellipse.
   *  @param  point   The point.
   *  @return         The feet, in order, every turn of the ellipse's to its 65,536th; none when it
   *                  is a circle and the point its center, where every point of it is as far.
   */
  [[nodiscard]] GEOXY_API auto Perpendiculars(const Ellipse &ellipse, Point point)
      -> std::vector<Place>;

  /**
   *  The feet of the perpendiculars from a point to a polyline's segments, and the vertices where
   *  the distance from it turns.
   *
   *  @param  polyline The polyline.
   *  @param  point    The point.
   *  @return          The places where the distance stops falling or rising, in order; a foot at a
   *                   vertex shared by two segments once; a foot within a segment at the double
   *                   nearest the exact one, within a few steps of a double of At at its parameter.
   */
  [[nodiscard]] GEOXY_API auto Perpendiculars(const Polyline &polyline, Point point)
      -> std::vector<Place>;

  /**
   *  The feet of the perpendiculars from a point to a spline, found on its pieces, and the knots
   *  where the distance from it turns.
   *
   *  Beside a weight far heavier than its neighbours a spline runs along its poles' legs within a
   *  sliver of its parameter, thinner than doubles are apart for a weight 2^52 times the rest. A
   *  foot there is told as it is found: its point to a few steps of a double, and its parameter the
   *  double nearest, which may be the end of the span for feet apart on the curve. So for a spline
   *  whose weights are more than twice each other, a foot's point is not the spline's at the foot's
   *  parameter worked out again, but the foot. Weights more than some 2^1070 times apart in one
   *  span are more than doubles hold beside each other: such a spline is answered as its own
   *  parameter tells it, which there is not to be relied on.
   *
   *  @param  spline The spline.
   *  @param  point  The point.
   *  @return        The places where the distance stops falling or rising, in order; none for a
   *                 spline that is not valid.
   */
  [[nodiscard]] GEOXY_API auto Perpendiculars(const Spline &spline, Point point)
      -> std::vector<Place>;

  /**
   *  The feet of the perpendiculars from a point to a quadratic Bézier curve.
   *
   *  @param  quadratic The curve.
   *  @param  point     The point.
   *  @return           The feet, in order.
   */
  [[nodiscard]] GEOXY_API auto Perpendiculars(const Bezier::Quadratic &quadratic, Point point)
      -> std::vector<Place>;

  /**
   *  The feet of the perpendiculars from a point to a cubic Bézier curve.
   *
   *  @param  cubic The curve.
   *  @param  point The point.
   *  @return       The feet, in order.
   */
  [[nodiscard]] GEOXY_API auto Perpendiculars(const Bezier::Cubic &cubic, Point point)
      -> std::vector<Place>;

  /**
   *  The feet of the perpendiculars from a point to any curve.
   *
   *  @param  curve The curve.
   *  @param  point The point.
   *  @return       The feet, as its own kind gives them.
   */
  [[nodiscard]] GEOXY_API auto Perpendiculars(const Curve &curve, Point point)
      -> std::vector<Place>;

  /**
   *  The place of a line nearest a point.
   *
   *  @param  line  The line.
   *  @param  point The point.
   *  @return       Its foot, or the nearer end; the start of a line of no length.
   */
  [[nodiscard]] GEOXY_API auto Nearest(const Line &line, Point point) -> Place;

  /**
   *  The place of a ray nearest a point.
   *
   *  @param  ray   The ray.
   *  @param  point The point.
   *  @return       Its foot, or its start.
   */
  [[nodiscard]] GEOXY_API auto Nearest(const Ray &ray, Point point) -> Place;

  /**
   *  The place of an xline nearest a point.
   *
   *  @param  xline The xline.
   *  @param  point The point.
   *  @return       Its foot; the point it passes through, for an xline of no way.
   */
  [[nodiscard]] GEOXY_API auto Nearest(const Xline &xline, Point point) -> Place;

  /**
   *  The place of an arc nearest a point.
   *
   *  @param  arc   The arc.
   *  @param  point The point.
   *  @return       The nearest of its feet and ends, the first along it of those as near; a foot at
   *                its point as Perpendiculars gives it.
   */
  [[nodiscard]] GEOXY_API auto Nearest(const Arc &arc, Point point) -> Place;

  /**
   *  The place of an ellipse nearest a point.
   *
   *  @param  ellipse The ellipse.
   *  @param  point   The point.
   *  @return         The nearest of its feet and ends, the first along it of those as near.
   */
  [[nodiscard]] GEOXY_API auto Nearest(const Ellipse &ellipse, Point point) -> Place;

  /**
   *  The place of a polyline nearest a point.
   *
   *  @param  polyline The polyline.
   *  @param  point    The point.
   *  @return          The nearest of its feet and vertices, the first along it of those as near, a
   *                   foot at its point as Perpendiculars gives it; the origin at 0 for no
   *                   vertices.
   */
  [[nodiscard]] GEOXY_API auto Nearest(const Polyline &polyline, Point point) -> Place;

  /**
   *  The place of a spline nearest a point.
   *
   *  @param  spline The spline.
   *  @param  point  The point.
   *  @return        The nearest of its feet and knots, the first along it of those as near, told as
   *                 its feet are; the origin at 0 for a spline that is not valid.
   */
  [[nodiscard]] GEOXY_API auto Nearest(const Spline &spline, Point point) -> Place;

  /**
   *  The place of a quadratic Bézier curve nearest a point.
   *
   *  @param  quadratic The curve.
   *  @param  point     The point.
   *  @return           The nearest of its feet and ends, the first along it of those as near.
   */
  [[nodiscard]] GEOXY_API auto Nearest(const Bezier::Quadratic &quadratic, Point point) -> Place;

  /**
   *  The place of a cubic Bézier curve nearest a point.
   *
   *  @param  cubic The curve.
   *  @param  point The point.
   *  @return       The nearest of its feet and ends, the first along it of those as near.
   */
  [[nodiscard]] GEOXY_API auto Nearest(const Bezier::Cubic &cubic, Point point) -> Place;

  /**
   *  The place of any curve nearest a point.
   *
   *  @param  curve The curve.
   *  @param  point The point.
   *  @return       The place, as its own kind gives it.
   */
  [[nodiscard]] GEOXY_API auto Nearest(const Curve &curve, Point point) -> Place;
} // namespace Geoxy

#endif /* GEOXY_NEAREST_HPP */
