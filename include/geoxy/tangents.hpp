/***************************************************************************************************
 *  @file      geoxy/tangents.hpp
 *  @brief     The places of a curve whose tangent passes through a point: where the lines from the
 *             point touch it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TANGENTS_HPP
#define GEOXY_TANGENTS_HPP

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
   *  The places of a line whose tangent passes through a point.
   *
   *  @param  line  The line.
   *  @param  point The point.
   *  @return       None: a line's tangent is the line, which touches nothing.
   */
  [[nodiscard]] GEOXY_API auto Tangents(const Line &line, Point point) -> std::vector<Place>;

  /**
   *  The places of a ray whose tangent passes through a point.
   *
   *  @param  ray   The ray.
   *  @param  point The point.
   *  @return       None: a ray's tangent is its line, which touches nothing.
   */
  [[nodiscard]] GEOXY_API auto Tangents(const Ray &ray, Point point) -> std::vector<Place>;

  /**
   *  The places of an xline whose tangent passes through a point.
   *
   *  @param  xline The xline.
   *  @param  point The point.
   *  @return       None: an xline's tangent is itself, which touches nothing.
   */
  [[nodiscard]] GEOXY_API auto Tangents(const Xline &xline, Point point) -> std::vector<Place>;

  /**
   *  The places of an arc whose tangent passes through a point, in closed form.
   *
   *  Seen from its center, scaled by its radius, the point is at a distance d and a heading h, and
   *  the arc touches the lines from it at the angles h plus or less acos(1 / d). An arc swept round
   *  more than 65,536 times has the places of its first 65,536 turns only, as its feet are.
   *
   *  @param  arc   The arc.
   *  @param  point The point.
   *  @return       The places, in order, every turn of the arc's to its 65,536th: two for a point
   *                outside its circle, the point itself on it, none within or for an arc of no
   *                radius.
   */
  [[nodiscard]] GEOXY_API auto Tangents(const Arc &arc, Point point) -> std::vector<Place>;

  /**
   *  The places of an ellipse whose tangent passes through a point, found on its quarters as exact
   *  conics.
   *
   *  Of an ellipse thin past what doubles hold, four in a turn at most are kept, those whose
   *  tangent passes nearest the point; one swept round more than 65,536 times has the places of its
   *  first 65,536 turns only, as its feet are.
   *
   *  @param  ellipse The ellipse.
   *  @param  point   The point.
   *  @return         The places, in order, every turn of the ellipse's to its 65,536th: two for a
   *                  point outside it, none within; the point itself on it may or may not be among
   *                  them.
   */
  [[nodiscard]] GEOXY_API auto Tangents(const Ellipse &ellipse, Point point) -> std::vector<Place>;

  /**
   *  The places of a polyline whose tangent passes through a point: on its bulged segments in
   *  closed form, and at the vertices where the polyline turns about the point the other way.
   *
   *  @param  polyline The polyline.
   *  @param  point    The point.
   *  @return          The places, in order.
   */
  [[nodiscard]] GEOXY_API auto Tangents(const Polyline &polyline, Point point)
      -> std::vector<Place>;

  /**
   *  The places of a spline whose tangent passes through a point, found on its pieces, and the
   *  knots where it turns about the point the other way.
   *
   *  A place beside a weight far heavier than its neighbours is told as it is found, as a foot of
   *  the spline is: its point to a few steps of a double, its parameter the double nearest.
   *
   *  @param  spline The spline.
   *  @param  point  The point.
   *  @return        The places, in order; the point itself, where the spline passes through it, may
   *                 or may not be among them. None for a spline that is not valid.
   */
  [[nodiscard]] GEOXY_API auto Tangents(const Spline &spline, Point point) -> std::vector<Place>;

  /**
   *  The places of a quadratic Bézier curve whose tangent passes through a point.
   *
   *  @param  quadratic The curve.
   *  @param  point     The point.
   *  @return           The places, in order; the point itself, where the curve passes through it,
   *                    may or may not be among them.
   */
  [[nodiscard]] GEOXY_API auto Tangents(const Bezier::Quadratic &quadratic, Point point)
      -> std::vector<Place>;

  /**
   *  The places of a cubic Bézier curve whose tangent passes through a point.
   *
   *  @param  cubic The curve.
   *  @param  point The point.
   *  @return       The places, in order; the point itself, where the curve passes through it, may
   *                or may not be among them.
   */
  [[nodiscard]] GEOXY_API auto Tangents(const Bezier::Cubic &cubic, Point point)
      -> std::vector<Place>;

  /**
   *  The places of any curve whose tangent passes through a point.
   *
   *  @param  curve The curve.
   *  @param  point The point.
   *  @return       The places, as its own kind gives them.
   */
  [[nodiscard]] GEOXY_API auto Tangents(const Curve &curve, Point point) -> std::vector<Place>;
} // namespace Geoxy

#endif /* GEOXY_TANGENTS_HPP */
