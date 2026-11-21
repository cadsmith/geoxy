/***************************************************************************************************
 *  @file      source/swift/rim.hpp
 *  @brief     The circle an arc or an ellipse lies on, as its kernel sees it: its frame, its
 *             quarters and the side of it a point is on.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_RIM_HPP
#define GEOXY_SWIFT_RIM_HPP

#include "swift/lane.hpp"
#include "swift/ring.hpp"
#include "swift/terms.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <utility>

namespace Geoxy::Swift
{
  /**
   *  An arc or an ellipse as the kernel takes it: the points p where |A (p - o)|^2 = rho^2, A the
   *  adjugate of its map M and rho the map's determinant D, so that A (p - o) is D times p's point
   *  on the unit circle's plane; an arc is its circle's map, its radius times the identity. A point
   *  at distance s from the carrier has |A (p - o)| within s / kappa of rho, kappa = 1 /
   *  sigma_max(M) at the least, which is no less than 1 / |M|_F; and the carrier's points that far
   *  from one of it lie within s / sigma_min of it on the unit circle, sigma_min = |D| kappa.
   */
  struct Ring
  {
    Matrix                Map    = {};      /**< The map of the unit circle.               */
    std::array<Vector, 2> Rows   = {};      /**< A, the map's adjugate, by rows.           */
    Number                Level  = {};      /**< rho^2 = D^2.                              */
    double                High   = 0.0;     /**< rho, at the most.                         */
    double                Low    = 0.0;     /**< rho, at the least.                        */
    double                Least  = 0.0;     /**< kappa, at the least.                      */
    double                Most   = 0.0;     /**< |M|_F, no less than sigma_max(M) and |A|. */
    double                Facing = 1.0;     /**< D's sign.                                 */
    Vector                From   = {};      /**< Its start's direction on the unit circle. */
    Vector                To     = {};      /**< Its end's.                                */
    double                Start  = 0.0;     /**< Its start, in radians.                    */
    double                Sweep  = 0.0;     /**< Its sweep, in radians.                    */
    bool                  Whole  = false;   /**< Whether it is whole.                      */
    const Curve          *Shape  = nullptr; /**< The curve it is.                          */
    std::array<Number, 2> Hub    = {};      /**< Its center, with its bound: the map's origin
                                                 exactly but for a bulge's, in double-double. */
    double Drift = 0.0;          /**< How far its true center may be from the map's origin, a
                                      coordinate: nothing but for a bulge's.                    */
    double Stride = 0.0;         /**< The most a step of its curve's place moves its point, for a
                                      bulge's: a polyline's.                                    */
    const Hump *Chord = nullptr; /**< The bulge it is the circle of, whose arc is the
                                      circle's side of its chord; none for an arc or an
                                      ellipse.                                              */
  };

  /**
   *  An arc or an ellipse as the kernel takes it. D is worked out in double-double, its bound a
   *  step of the doubles below the normal ones at the most, so that rho's two ends are the
   *  double-double's less and more that; the Frobenius norm rounded and grown is no less than the
   *  true one.
   *
   *  @param  curve The curve, an arc or an ellipse.
   *  @return       It as a ring; nothing for an arc of a radius not above zero, one swept no way or
   *                round more than once, one too flat to invert, or one whose numbers lie outside
   *                the range the kernels take.
   */
  auto Rimmed(const Curve &curve) -> std::optional<Ring>;

  /**
   *  A bulge's circle as a ring: A the identity, so that kappa is 1 and rho the radius, the center
   *  C = M + n (1 - b^2) / (4 b) and rho^2 = |c|^2 (1 + b^2)^2 / (16 b^2) worked out in
   *  double-double from the exact differences of the ends, the center rounded and its bound the
   *  drift. Its map, the radius times the identity about the center, only sizes the box about it. A
   *  step of the polyline's place on the segment moves its point by no more than r times the arc's
   *  sweep 4 atan |b|, no more than |c| (1 + b^2), times a step at the segment's end.
   *
   *  @param  hump The bulge, its numbers within the range the kernels take, its bulge not below
   *               Straight.
   *  @return      Its circle; nothing where its center or radius lies outside that range, or its
   *               radius may be nothing.
   */
  auto Encircled(const Hump &hump) -> std::optional<Ring>;

  /**
   *  A point on the unit circle's plane of a ring, A (p - o) and its slack, in doubles: the gap
   *  rounded once, within u of itself a coordinate, and each product and the sum once more, all
   *  within 2^-50 of the terms' sizes.
   *
   *  @param  ring  The ring.
   *  @param  point The point.
   *  @return       A (p - o), and how far each coordinate may be from the exact one.
   */
  auto Seen(const Ring &ring, Point point) -> std::pair<Vector, double>;

  /**
   *  A coordinate of a ring's center, with how far the true one may be from it: exact for an arc or
   *  an ellipse, a bulge's circle's in double-double within its bound.
   *
   *  @param  ring The ring.
   *  @param  axis Which coordinate: 0 for x.
   *  @return      The coordinate, with its bound.
   */
  auto Hubbed(const Ring &ring, std::size_t axis) -> Number;

  /**
   *  Which side of a bulge's chord a crossing lies on, as its arc is the circle's points on the
   *  bulge's side, right of the chord for a bulge above zero: by c x (X - P0), within 2^-48 of its
   *  terms rounded, and |c| times how far X may be. Every point of the circle within the margin of
   *  the crossing lies on its side where c x (X - P0) is clear of |c| times the margin: the circle
   *  meets the chord's line only at the ends, which then lie farther than the margin from it.
   *
   *  @param  hump     The bulge.
   *  @param  crossing The crossing's point, each coordinate with its bound.
   *  @param  margin   The margin, and how far the rounded point may be from the true one.
   *  @return          In where on the arc, out where off it; near where the side is not clear.
   */
  auto Sidled(const Hump &hump, const std::array<Number, 2> &crossing, double margin) -> Side;

  /**
   *  How far a point is from a ring's carrier, at the least: kappa times how far |A (p - o)| is
   *  from rho, less what rounding may take off, its root within a step of itself.
   *
   *  @param  ring  The ring.
   *  @param  point The point.
   *  @return       The distance, at the least.
   */
  auto RingGap(const Ring &ring, Point point) -> double;

  /**
   *  Where a direction on a ring's unit circle lies in its sweep: by the signs of its cross
   *  products with the sweep's ends, each decided where it is clear of a margin, an angle on the
   *  circle as a share of the direction's length, and of rounding. A whole ring's seam is an end.
   *
   *  @param  ring  The ring.
   *  @param  seen  The direction, of any length, and how far it may be from the exact one.
   *  @param  angle The margin, the most the angle may be off.
   *  @return       Inside or out; near where a sign is not clear.
   */
  auto Swept(const Ring &ring, const std::pair<Vector, double> &seen, double angle) -> Side;

  /**
   *  An arc or an ellipse cut into quarters, each met with the ring as a rational quadratic piece,
   *  and what is worked out once about the two: its map taken to the ring's unit circle's plane, B
   *  = A M', and its center's gap there, k = A (o' - o), in double-double, so that a point of it at
   *  a direction c of its unit circle is there at k + B c. A bulge's circle's center o is its true
   *  one, in double-double, and the gap rounded is off by the double-double's low part and bound
   *  past the rounding QuarterPiece charges.
   */
  struct Bend
  {
    const Ring           *Bent   = nullptr; /**< The arc or the ellipse cut. */
    std::array<Number, 4> Turned = {};      /**< B, by rows.                 */
    std::array<Number, 2> Offset = {};      /**< k.                          */
    Vector                Gap    = {};      /**< o' - o, rounded.            */
    double                Drift  = 0.0;     /**< How far the rounded gap may be off,
                                                 a coordinate, past its rounding: a
                                                 bulge's circle's, nothing else.    */
  };

  /**
   *  An arc or an ellipse cut into quarters, about a ring.
   *
   *  @param  bent The one cut.
   *  @param  ring The ring.
   *  @return      What is worked out once about the two.
   */
  auto Bendings(const Ring &bent, const Ring &ring) -> Bend;
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_RIM_HPP */
