/***************************************************************************************************
 *  @file      source/swift/gather.hpp
 *  @brief     The hits of a ring and a curve gathered piece by piece, and told apart.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_GATHER_HPP
#define GEOXY_SWIFT_GATHER_HPP

#include "scratch.hpp"
#include "swift/circuit.hpp"
#include "swift/rim.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/tolerance.hpp"

#include <cstddef>
#include <optional>

namespace Geoxy::Swift
{
  /**
   *  Whether a curve is an arc or an ellipse.
   *
   *  @param  curve The curve.
   *  @return       Whether it is.
   */
  auto Conic(const Curve &curve) -> bool;

  /**
   *  Whether a curve is a Bézier curve or a spline, whose pieces are polynomial.
   *
   *  @param  curve The curve.
   *  @return       Whether it is.
   */
  auto Polynomial(const Curve &curve) -> bool;

  /**
   *  The crossings of a ring with the other curve: a Bézier curve as one polynomial piece, a spline
   *  span by span, an arc or an ellipse by its quarters that meet its sweep. A bulge's circle is
   *  met only with pieces near its arc, as Bulging bounds it.
   *
   *  @param  carrier   The ring.
   *  @param  other     The other curve.
   *  @param  tolerance The tolerance.
   *  @return           The crossings; nothing where it is not sure.
   */
  auto Gathered(const Ring &carrier, const Curve &other, const Tolerance &tolerance)
      -> std::optional<List<Found>>;

  /**
   *  Two curves as the ring kernel meets them: which is the piece's curve, and the tolerance.
   */
  struct RingPair
  {
    const Curves    *Pair      = nullptr; /**< The two curves.                          */
    std::size_t      Side      = 0;       /**< Which is the piece's curve: 0 the first. */
    Geoxy::Tolerance Tolerance = {};      /**< The tolerance.                           */
  };

  /**
   *  A ring met with the other curve and its crossings made hits. A crossing the thorough way's
   *  placing would not place leaves fewer hits than crossings, and the kernel not sure, as where it
   *  found none.
   *
   *  @param  carrier The ring.
   *  @param  pairing The two curves and the tolerance.
   *  @param  into    Where the hits go.
   *  @return         Whether it was sure.
   */
  auto Rung(const Ring &carrier, const RingPair &pairing, Intersections &into) -> bool;
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_GATHER_HPP */
