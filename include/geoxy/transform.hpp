/***************************************************************************************************
 *  @file      geoxy/transform.hpp
 *  @brief     Each curve under an affine map, its kind kept where the map allows.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TRANSFORM_HPP
#define GEOXY_TRANSFORM_HPP

#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

namespace Geoxy
{
  /**
   *  A line mapped.
   *
   *  @param  line The line.
   *  @param  map  The map.
   *  @return      The line between its ends mapped, its parameter kept.
   */
  [[nodiscard]] GEOXY_API auto Transformed(const Line &line, const Matrix &map) -> Line;

  /**
   *  A ray mapped.
   *
   *  @param  ray The ray.
   *  @param  map The map.
   *  @return     Its start and its way mapped, its parameter kept.
   */
  [[nodiscard]] GEOXY_API auto Transformed(const Ray &ray, const Matrix &map) -> Ray;

  /**
   *  An xline mapped.
   *
   *  @param  xline The xline.
   *  @param  map   The map.
   *  @return       Its point and its way mapped, its parameter kept.
   */
  [[nodiscard]] GEOXY_API auto Transformed(const Xline &xline, const Matrix &map) -> Xline;

  /**
   *  An arc mapped: an arc under a similarity, an ellipse under any other map.
   *
   *  Under a similarity, its center is mapped, its radius times the length of the map's x axis, and
   *  its start turned by the heading of that axis, or, where the map mirrors, mirrored about it
   *  with its sweep negated: an angle held in degrees stays in degrees, exactly so for maps by
   *  quarter turns. Under any other map it is the ellipse of the map composed with the circle's
   *  own, its angles kept as its parameters.
   *
   *  @param  arc The arc.
   *  @param  map The map.
   *  @return     The arc or the ellipse, its parameter kept.
   */
  [[nodiscard]] GEOXY_API auto Transformed(const Arc &arc, const Matrix &map) -> Curve;

  /**
   *  An ellipse mapped.
   *
   *  @param  ellipse The ellipse.
   *  @param  map     The map.
   *  @return         The map composed with its own, its start and sweep kept, and so its parameter.
   */
  [[nodiscard]] GEOXY_API auto Transformed(const Ellipse &ellipse, const Matrix &map) -> Ellipse;

  /**
   *  A polyline mapped: a polyline where the map keeps its segments lines and arcs, else the spline
   *  that is its image exactly.
   *
   *  Under a similarity its vertices are mapped and its bulges kept, negated where the map mirrors;
   *  under any other map so are those of a polyline whose segments are all lines, as Bent tells
   *  them. Any other is a rational quadratic spline: its lines raised to quadratics, its arcs cut
   *  into pieces of at most a quarter turn, each piece over a unit of its parameter, so that it
   *  passes through every vertex mapped, to the bit, at the knots between.
   *
   *  @param  polyline The polyline.
   *  @param  map      The map.
   *  @return          The polyline, its parameter kept, or the spline.
   */
  [[nodiscard]] GEOXY_API auto Transformed(const Polyline &polyline, const Matrix &map) -> Curve;

  /**
   *  A spline mapped.
   *
   *  @param  spline The spline.
   *  @param  map    The map.
   *  @return        Its poles mapped, its knots and weights kept, and so its parameter.
   */
  [[nodiscard]] GEOXY_API auto Transformed(const Spline &spline, const Matrix &map) -> Spline;

  /**
   *  A quadratic Bézier curve mapped.
   *
   *  @param  quadratic The curve.
   *  @param  map       The map.
   *  @return           Its control points mapped, its parameter kept.
   */
  [[nodiscard]] GEOXY_API auto Transformed(const Bezier::Quadratic &quadratic, const Matrix &map)
      -> Bezier::Quadratic;

  /**
   *  A cubic Bézier curve mapped.
   *
   *  @param  cubic The curve.
   *  @param  map   The map.
   *  @return       Its control points mapped, its parameter kept.
   */
  [[nodiscard]] GEOXY_API auto Transformed(const Bezier::Cubic &cubic, const Matrix &map)
      -> Bezier::Cubic;

  /**
   *  Any curve mapped.
   *
   *  @param  curve The curve.
   *  @param  map   The map.
   *  @return       It mapped as its own kind is.
   */
  [[nodiscard]] GEOXY_API auto Transformed(const Curve &curve, const Matrix &map) -> Curve;
} // namespace Geoxy

#endif /* GEOXY_TRANSFORM_HPP */
