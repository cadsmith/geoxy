/***************************************************************************************************
 *  @file      source/swift/met.hpp
 *  @brief     Two curves met by a swift kernel where one is straight, or both are Bézier curves,
 *             and the thorough way would leave nothing to decide: each crossing the double nearest
 *             the true point.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_MET_HPP
#define GEOXY_SWIFT_MET_HPP

#include "geoxy/intersect.hpp"
#include "geoxy/tolerance.hpp"

namespace Geoxy::Swift
{
  /**
   *  Two curves met by a swift kernel: a line, a ray or an endless line, either way round, with an
   *  ellipse, a polyline or a quadratic or cubic Bézier curve; two quadratic or cubic Bézier
   *  curves; an arc or an ellipse with a Bézier curve, a polynomial spline or another conic, two
   *  arcs aside; a quadratic Bézier curve with a rational spline or one past a cubic, on its
   *  parabola; a cubic one with such a spline to degree five, on its implicit cubic; or a polyline
   *  with a spline, an ellipse or a Bézier curve, leg by leg. Each kernel answers that they do not
   *  meet, that they cross at these points, each the double nearest the true one, or that it is not
   *  sure; it is sure only where its error bounds prove its points and the thorough way's tolerance
   *  would find these crossings and nothing else: every crossing steep, far from either curve's
   *  ends and from any other, and wherever the curves come within eight tolerances of each other, a
   *  crossing there. The hits are placed as the thorough way places them.
   *
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @param  into      Where the hits go, empty on entry, and left so where it is not sure.
   *  @return           Whether it was sure.
   */
  auto Met(const Curves &curves, const Tolerance &tolerance, Intersections &into) -> bool;
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_MET_HPP */
