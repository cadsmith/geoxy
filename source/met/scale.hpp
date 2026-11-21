/***************************************************************************************************
 *  @file      source/met/scale.hpp
 *  @brief     The size of a curve's numbers and the grain a meeting is held to; curves too large or
 *             too small, or of ways too long or too short, met at a size where they are not.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_MET_SCALE_HPP
#define GEOXY_MET_SCALE_HPP

#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/xline.hpp"

#include <array>
#include <cstddef>
#include <span>

namespace Geoxy::Met
{
  /**
   *  A line's points: its ends; the points a curve is given by, where they are all it is.
   *
   *  @param  shape The curve.
   *  @return       Its points.
   */
  auto Given(const Line &shape) -> std::array<Point, 2>;

  /**
   *  A ray's points: its start.
   *
   *  @param  shape The curve.
   *  @return       Its points.
   */
  auto Given(const Ray &shape) -> std::array<Point, 1>;

  /**
   *  An endless line's points: the point it runs through.
   *
   *  @param  shape The curve.
   *  @return       Its points.
   */
  auto Given(const Xline &shape) -> std::array<Point, 1>;

  /**
   *  A spline's points: its control points.
   *
   *  @param  shape The curve.
   *  @return       Its points.
   */
  auto Given(const Spline &shape) -> std::span<const Point>;

  /**
   *  A quadratic Bézier curve's points: its control points.
   *
   *  @param  shape The curve.
   *  @return       Its points.
   */
  auto Given(const Bezier::Quadratic &shape) -> std::array<Point, 3>;

  /**
   *  A cubic Bézier curve's points: its control points.
   *
   *  @param  shape The curve.
   *  @return       Its points.
   */
  auto Given(const Bezier::Cubic &shape) -> std::array<Point, 4>;

  /**
   *  The size of a curve's numbers: the largest magnitude of the points it is given by, its center
   *  and reach for an arc or an ellipse, and for a polyline each bulge's circle's too, finite ones
   *  only.
   *
   *  @param  curve The curve.
   *  @return       The size, 0 for none.
   */
  auto Size(const Curve &curve) -> double;

  /**
   *  A tolerance no finer than the curves' numbers can hold: its length at least 16 steps of a
   *  double at their size, below which where two curves lie, and so whether they stay together or
   *  part, is not known.
   *
   *  @param  first     One curve.
   *  @param  second    The other.
   *  @param  tolerance The tolerance asked for.
   *  @return           The tolerance used.
   */
  auto Grained(const Curve &first, const Curve &second, Tolerance tolerance) -> Tolerance;
} // namespace Geoxy::Met

#endif /* GEOXY_MET_SCALE_HPP */
