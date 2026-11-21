/***************************************************************************************************
 *  @file      source/swift/pieced.hpp
 *  @brief     Polynomial splines met piece by piece by the swift kernels: each span a quadratic or
 *             cubic piece of its own, met with each piece of the other curve whose box comes near
 *             it, where the thorough way, every pair of parts weighed at once, would give the same.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_PIECED_HPP
#define GEOXY_SWIFT_PIECED_HPP

#include "swift/lane.hpp"
#include "swift/section.hpp"

#include "geoxy/intersect.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"

#include <optional>
#include <span>

namespace Geoxy::Swift
{
  /**
   *  Two curves met piece by piece, where one is a polynomial spline of degree two or three and the
   *  other is one too, or a quadratic or cubic Bézier curve: each pair of pieces whose boxes come
   *  near met by the twin kernel, the rest far apart; then every two crossings told apart as the
   *  thorough way tells them, by a place halfway between them along either curve that the kernels'
   *  bands show is far from the other curve, and every crossing's places such that the thorough way
   *  holds its candidates to the tolerance alone.
   *
   *  The thorough way weighs each pair of parts, a spline's parts being its spans, on its own, and
   *  then the hits of all of them together; the twin kernel answers for a pair of pieces as the
   *  thorough way would for the two pieces alone, so for the curves where nothing ties one pair's
   *  crossings to another's.
   *
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @param  into      Where the hits go, empty on entry, and left so where it is not sure.
   *  @return           Whether it was sure; nothing where neither curve is a spline, or either is
   *                    straight or of a kind the kernels do not meet with a spline.
   */
  auto Pieced(const Curves &curves, const Tolerance &tolerance, Intersections &into)
      -> std::optional<bool>;

  /**
   *  A spline, an ellipse or a quadratic or cubic Bézier curve met with a polyline, either first:
   *  each straight segment near the curve met by the kernel of a line and that curve, its place on
   *  the polyline its index ahead of the line's; each bulge by the ring kernel on its circle;
   *  segments far from the curve passed; every two crossings told apart as the thorough way tells
   *  them, by the curve halfway between them far from the polyline. Every vertex then lies farther
   *  than the reach from the curve, each being a segment's end or a far one's, and the hits are
   *  placed as the thorough way places them: on a straight segment at the place proved the double
   *  nearest the true one, or an ellipse's crossing by the point's foot; on a bulge and an ellipse
   *  by the point's angle; on a spline or a Bézier curve at the place proved.
   *
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @param  into      Where the hits go, empty on entry, and left so where it is not sure.
   *  @return           Whether it was sure; nothing where the two are not such a curve and a
   *                    polyline.
   */
  auto Laced(const Curves &curves, const Tolerance &tolerance, Intersections &into)
      -> std::optional<bool>;

  /**
   *  A line, a ray or an endless line met with a polynomial spline of degree two or three: each
   *  span whose box comes near the straight curve met by the kernel of a line and a Bézier curve,
   *  the rest far from it.
   *
   *  @param  lane      The straight curve.
   *  @param  spline    The spline.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @return           The crossings, each placed on the spline by its own parameter; nothing where
   *                    the thorough way must decide.
   */
  auto Strung(const Lane &lane, const Spline &spline, const Tolerance &tolerance)
      -> std::optional<Spots>;

  /**
   *  A line, a ray or an endless line met with a polynomial spline of degree two or three already
   *  cut into its sections, as Strung meets it.
   *
   *  @param  lane      The straight curve.
   *  @param  sections  The spline's sections, in order.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @return           The crossings, each placed on the spline by its own parameter; nothing where
   *                    the thorough way must decide.
   */
  auto Strung(const Lane &lane, std::span<const Section> sections, const Tolerance &tolerance)
      -> std::optional<Spots>;

  /**
   *  A line, a ray or an endless line met with a spline of any degree from two to fifteen, rational
   *  or not: the line's equation times the spline's weight is a spline function of the spline's own
   *  knots and degree, a coefficient per pole; each span whose poles come near the straight curve
   *  has its Bernstein coefficients worked out as blossoms of those, its roots isolated and each
   *  proved and placed as the kernel of a line and a Bézier curve proves and places them, the point
   *  worked out from the spline's own poles and knots in double-double.
   *
   *  @param  lane      The straight curve.
   *  @param  spline    The spline.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @return           The crossings, each placed on the spline by its own parameter; nothing where
   *                    the thorough way must decide.
   */
  auto Skeined(const Lane &lane, const Spline &spline, const Tolerance &tolerance)
      -> std::optional<Spots>;

  /**
   *  Whether the thorough way holds a crossing's candidates on a spline's span to the tolerance
   *  alone: where either curve is a spline it holds them to twice what a step of each place moves
   *  its curve where that is more, and here it is less, the span's step no more than at twice its
   *  farthest end and the other curve's as given.
   *
   *  @param  section   The span the crossing is on.
   *  @param  stride    The most a step of the other curve's place may move it there.
   *  @param  tolerance The tolerance.
   *  @return           Whether it does.
   */
  auto Fine(const Broad &section, double stride, const Tolerance &tolerance) -> bool;
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_PIECED_HPP */
