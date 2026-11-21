/***************************************************************************************************
 *  @file      geoxy/derivative.hpp
 *  @brief     The derivative of a curve with respect to its parameter: its tangent, and its speed.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_DERIVATIVE_HPP
#define GEOXY_DERIVATIVE_HPP

#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

namespace Geoxy
{
  /**
   *  The derivative of a line with respect to its parameter.
   *
   *  @param  line  The line.
   *  @param  along 0 at its start, 1 at its end.
   *  @return       Its end less its start, the same everywhere.
   */
  [[nodiscard]] GEOXY_API auto Derivative(const Line &line, double along) -> Vector;

  /**
   *  The derivative of a ray with respect to its parameter.
   *
   *  @param  ray   The ray.
   *  @param  along How many of its ways from its start.
   *  @return       Its way.
   */
  [[nodiscard]] GEOXY_API auto Derivative(const Ray &ray, double along) -> Vector;

  /**
   *  The derivative of an xline with respect to its parameter.
   *
   *  @param  xline The xline.
   *  @param  along How many of its ways from its point.
   *  @return       Its way.
   */
  [[nodiscard]] GEOXY_API auto Derivative(const Xline &xline, double along) -> Vector;

  /**
   *  The derivative of an arc with respect to its parameter.
   *
   *  @param  arc   The arc.
   *  @param  along The fraction of its sweep.
   *  @return       The tangent, its length the radius times the sweep in radians.
   */
  [[nodiscard]] GEOXY_API auto Derivative(const Arc &arc, double along) -> Vector;

  /**
   *  The derivative of an ellipse with respect to its parameter.
   *
   *  @param  ellipse The ellipse.
   *  @param  along   The fraction of its sweep.
   *  @return         The tangent, its length that of the mapped unit tangent times the sweep in
   *                  radians.
   */
  [[nodiscard]] GEOXY_API auto Derivative(const Ellipse &ellipse, double along) -> Vector;

  /**
   *  The derivative of a polyline with respect to its parameter.
   *
   *  @param  polyline The polyline.
   *  @param  along    A segment's index and the fraction along it.
   *  @return          The tangent of the segment the value falls in, a segment's end taking the
   *                   segment it ends; its length the segment's own speed.
   */
  [[nodiscard]] GEOXY_API auto Derivative(const Polyline &polyline, double along) -> Vector;

  /**
   *  The derivative of a spline with respect to its parameter.
   *
   *  @param  spline The spline.
   *  @param  along  A knot value.
   *  @return        The tangent; a rational spline's by the quotient rule; none where the spline
   *                 cannot be read, as At tells it.
   */
  [[nodiscard]] GEOXY_API auto Derivative(const Spline &spline, double along) -> Vector;

  /**
   *  The derivative of a quadratic Bézier curve with respect to its parameter.
   *
   *  @param  quadratic The curve.
   *  @param  along     0 at its start, 1 at its end.
   *  @return           The tangent.
   */
  [[nodiscard]] GEOXY_API auto Derivative(const Bezier::Quadratic &quadratic, double along)
      -> Vector;

  /**
   *  The derivative of a cubic Bézier curve with respect to its parameter.
   *
   *  @param  cubic The curve.
   *  @param  along 0 at its start, 1 at its end.
   *  @return       The tangent.
   */
  [[nodiscard]] GEOXY_API auto Derivative(const Bezier::Cubic &cubic, double along) -> Vector;

  /**
   *  The derivative of any curve with respect to its parameter.
   *
   *  @param  curve The curve.
   *  @param  along The value, as its own kind reads one.
   *  @return       The tangent.
   */
  [[nodiscard]] GEOXY_API auto Derivative(const Curve &curve, double along) -> Vector;

} // namespace Geoxy

#endif /* GEOXY_DERIVATIVE_HPP */
