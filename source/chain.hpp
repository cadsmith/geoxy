/***************************************************************************************************
 *  @file      source/chain.hpp
 *  @brief     Curves as runs of rational Bézier pieces, and runs of pieces as one spline.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CHAIN_HPP
#define GEOXY_CHAIN_HPP

#include "kernel/piece.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/spline.hpp"

#include <span>
#include <vector>

namespace Geoxy
{
  /**
   *  A polyline as rational Bézier pieces, end to end: a line for each straight segment, and each
   *  bulged one as rational quadratic pieces of at most a quarter turn each.
   *
   *  Each piece starts and ends at the polyline's own points, At's at its vertices and at the
   *  fractions of a segment the pieces cut it at, and their weights there are 1, so that a joint
   *  holds the point to the bit.
   *
   *  @param  polyline The polyline.
   *  @return          The pieces, in order; none for a polyline of no segment.
   */
  auto Sliced(const Polyline &polyline) -> std::vector<Kernel::Piece>;

  /**
   *  Any curve as rational Bézier pieces, end to end: a line, a Bézier curve as its one piece; an
   *  arc or an ellipse as quadratic pieces of at most a quarter turn; a polyline as Sliced gives
   *  it; a spline as its spans' pieces.
   *
   *  @param  curve The curve.
   *  @return       The pieces, in order; none for a curve with no end, a ray or an xline, for an
   *                arc or ellipse swept more than a turn or not by a number, a polyline of no
   *                segment and a spline that is not valid.
   */
  auto Sliced(const Curve &curve) -> std::vector<Kernel::Piece>;

  /**
   *  A piece run the other way.
   *
   *  @param  piece The piece.
   *  @return       Its control points in the other order.
   */
  auto Flipped(Kernel::Piece piece) -> Kernel::Piece;

  /**
   *  Pieces end to end as one spline, each lifted to the highest degree among them and over a unit
   *  of its parameter, a knot of full multiplicity at each joint.
   *
   *  Each piece's weights are scaled, which leaves it the same curve, so that its first is the last
   *  of the piece before: the two then share that control point, the earlier one's end. Where every
   *  weight comes to 1 the spline has none.
   *
   *  @pre           There is a piece, and each has a degree of 1 at least.
   *  @param  pieces The pieces, in order.
   *  @return        The spline.
   */
  auto Joined(std::span<const Kernel::Piece> pieces) -> Spline;
} // namespace Geoxy

#endif /* GEOXY_CHAIN_HPP */
