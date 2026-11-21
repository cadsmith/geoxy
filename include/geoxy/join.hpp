/***************************************************************************************************
 *  @file      geoxy/join.hpp
 *  @brief     Curves end to end, joined into one: JOIN.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_JOIN_HPP
#define GEOXY_JOIN_HPP

#include "geoxy/api.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/tolerance.hpp"

#include <optional>
#include <span>

namespace Geoxy
{
  /**
   *  Curves end to end, joined into one curve, as JOIN joins them.
   *
   *  The curves come in order, each with an end within the tolerance of an end of the next, and
   *  each is turned as it must be: the first so that its end, not its start, touches the second,
   *  where only its start does, and each after so that it starts where the one before ends. Each
   *  joint takes the earlier curve's end.
   *
   *  Where every curve is a line, an arc or a polyline, they join into a polyline, an arc as its
   *  bulge or, swept a whole turn, as two halves; else into a spline of the highest degree among
   *  them, each a run of rational Bézier pieces, an arc or an ellipse cut into pieces of at most a
   *  quarter turn, each piece over a unit of the spline's parameter and joined to the next at a
   *  knot of full multiplicity. Where the chain ends within the tolerance of its start it is
   *  closed: a polyline by its flag, a spline by its last pole set to its first.
   *
   *  @param  curves    The curves, in order.
   *  @param  tolerance How far apart two ends may be and still join.
   *  @return           The joined curve; none for no curves, for a ray or an xline, an arc or an
   *                    ellipse swept more than a turn, a polyline of no segment or a spline that is
   *                    not valid, or where two curves in a row do not touch.
   */
  [[nodiscard]] GEOXY_API auto Join(std::span<const Curve> curves, const Tolerance &tolerance)
      -> std::optional<Curve>;
} // namespace Geoxy

#endif /* GEOXY_JOIN_HPP */
