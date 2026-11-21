/***************************************************************************************************
 *  @file      geoxy/at.hpp
 *  @brief     The point of a curve at a value of its parameter.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_AT_HPP
#define GEOXY_AT_HPP

#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

namespace Geoxy
{
  /**
   *  The point of a line at a value of its parameter.
   *
   *  @param  line  The line.
   *  @param  along 0 at its start, 1 at its end, and past them the line run on.
   *  @return       The point: exactly the start at 0 and exactly the end at 1.
   */
  [[nodiscard]] GEOXY_API auto At(const Line &line, double along) -> Point;

  /**
   *  The point of a ray at a value of its parameter.
   *
   *  @param  ray   The ray.
   *  @param  along How many of its ways from its start.
   *  @return       The point.
   */
  [[nodiscard]] GEOXY_API auto At(const Ray &ray, double along) -> Point;

  /**
   *  The point of an xline at a value of its parameter.
   *
   *  @param  xline The xline.
   *  @param  along How many of its ways from its point.
   *  @return       The point.
   */
  [[nodiscard]] GEOXY_API auto At(const Xline &xline, double along) -> Point;

  /**
   *  The point of an arc at a value of its parameter.
   *
   *  @param  arc   The arc.
   *  @param  along The fraction of its sweep, and past 0 and 1 its circle run on.
   *  @return       The point, exact where its angle is a whole number of quarter turns.
   */
  [[nodiscard]] GEOXY_API auto At(const Arc &arc, double along) -> Point;

  /**
   *  The point of an ellipse at a value of its parameter.
   *
   *  @param  ellipse The ellipse.
   *  @param  along   The fraction of its sweep, and past 0 and 1 the ellipse run on.
   *  @return         The point.
   */
  [[nodiscard]] GEOXY_API auto At(const Ellipse &ellipse, double along) -> Point;

  /**
   *  The point of a polyline at a value of its parameter.
   *
   *  A bulged segment's points are worked out from its start along the chord, not from its center,
   *  which lies far off for a small bulge: at survey coordinates the center's rounding would move
   *  them.
   *
   *  @param  polyline The polyline.
   *  @param  along    A segment's index and the fraction along it; before the first segment or past
   *                   the last, that segment run on.
   *  @return          The point, exactly a vertex at a whole number; the origin for no vertices,
   *                   and the vertex for one.
   */
  [[nodiscard]] GEOXY_API auto At(const Polyline &polyline, double along) -> Point;

  /**
   *  The point of a spline at a value of its parameter, by de Boor's algorithm.
   *
   *  One point must not cost a pass over the whole spline, so only what the point is read from is
   *  checked: the counts Valid asks for, and around the parameter's span knots that are finite and
   *  do not fall and weights that are finite and positive.
   *
   *  @param  spline The spline.
   *  @param  along  A knot value; before its first span or past its last, that span's polynomial
   *                 run on.
   *  @return        The point; the origin where the spline cannot be read.
   */
  [[nodiscard]] GEOXY_API auto At(const Spline &spline, double along) -> Point;

  /**
   *  The point of a quadratic Bézier curve at a value of its parameter, by de Casteljau's
   *  algorithm.
   *
   *  @param  quadratic The curve.
   *  @param  along     0 at its start, 1 at its end.
   *  @return           The point: exactly the start at 0 and exactly the end at 1.
   */
  [[nodiscard]] GEOXY_API auto At(const Bezier::Quadratic &quadratic, double along) -> Point;

  /**
   *  The point of a cubic Bézier curve at a value of its parameter, by de Casteljau's algorithm.
   *
   *  @param  cubic The curve.
   *  @param  along 0 at its start, 1 at its end.
   *  @return       The point: exactly the start at 0 and exactly the end at 1.
   */
  [[nodiscard]] GEOXY_API auto At(const Bezier::Cubic &cubic, double along) -> Point;

  /**
   *  The point of any curve at a value of its parameter.
   *
   *  @param  curve The curve.
   *  @param  along The value, as its own kind reads one.
   *  @return       The point.
   */
  [[nodiscard]] GEOXY_API auto At(const Curve &curve, double along) -> Point;
} // namespace Geoxy

#endif /* GEOXY_AT_HPP */
