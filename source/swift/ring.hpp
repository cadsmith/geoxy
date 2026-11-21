/***************************************************************************************************
 *  @file      source/swift/ring.hpp
 *  @brief     An arc or an ellipse met with a quadratic or cubic Bézier curve, a spline or another
 *             arc or ellipse, or a quadratic or cubic Bézier curve with a spline, by a swift
 *             kernel: the conic's or the curve's implicit equation on the other curve's piece, its
 *             roots isolated and each crossing the double nearest the true point.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_RING_HPP
#define GEOXY_SWIFT_RING_HPP

#include "swift/lane.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"

#include <cstddef>
#include <optional>

namespace Geoxy::Swift
{
  /**
   *  A bulged segment of a polyline, as the ring kernel meets a spline with its circle.
   */
  struct Hump
  {
    Point       Start = {};  /**< Its first vertex, P0.          */
    Point       End   = {};  /**< Its second vertex, P1.         */
    double      Bend  = 0.0; /**< Its bulge, b, not nothing.     */
    std::size_t Index = 0;   /**< Which segment of the polyline. */
  };

  /**
   *  Two curves met by the ring kernel, where one is an arc or an ellipse and the other a quadratic
   *  or cubic Bézier curve, a spline of degree two to eight, polynomial or rational, an ellipse, or
   *  an arc met with an ellipse: the conic's equation |A (p - o)|^2 = rho^2 on the other curve, of
   *  twice its piece's degree in its parameter, its roots isolated in doubles and each polished
   *  once in double-double; a spline taken span by span, an ellipse or an arc by its quarters, each
   *  a rational quadratic piece of exact weights. It is sure only where the thorough way would find
   *  these crossings and nothing else: wherever the curves come within eight tolerances of each
   *  other, a steep crossing there, far from either curve's ends and from any other. Each crossing
   *  is the double nearest the true point, a Bézier curve's place the double nearest the true one,
   *  and the hits are placed as the thorough way places them.
   *
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @param  into      Where the hits go, empty on entry, and left so where it is not sure.
   *  @return           Whether it was sure; nothing where the curves are not such a pair, or are
   *                    two arcs, which the chord form meets.
   */
  auto Ringed(const Curves &curves, const Tolerance &tolerance, Intersections &into)
      -> std::optional<bool>;

  /**
   *  A bulge of a polyline met with a spline, an ellipse or a quadratic or cubic Bézier curve by
   *  the ring kernel: its circle |X - C|^2 = r^2, the center and the radius worked out in
   *  double-double from its ends and bulge, C = M + n (1 - b^2) / (4 b) for the chord's middle M
   *  and the chord turned a quarter left n, r^2 = |c|^2 (1 + b^2)^2 / (16 b^2), the other curve's
   *  pieces taken about the true center in double-double and the rounded center charged its drift
   *  wherever else it is weighed; the arc the circle's points on the bulge's side of the chord, and
   *  only pieces near its box met. Sure as the kernel is with an arc, the arc's sweep told by the
   *  chord, each crossing's side of it clear by the margin about it, so that the bulge's ends lie
   *  farther than the reach from the other curve; a spline's stride the polyline's at the segment.
   *
   *  @param  hump      The bulge, its numbers within the range the kernels take.
   *  @param  other     The other curve.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @return           The crossings, the other curve's place second: a spline's or a Bézier
   *                    curve's own, the double nearest the true one, an ellipse's not, its point
   *                    placed by its angle; nothing where the kernel is not sure or the other curve
   *                    is none of these.
   */
  auto Humped(const Hump &hump, const Curve &other, const Tolerance &tolerance)
      -> std::optional<Spots>;
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_RING_HPP */
