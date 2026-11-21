/***************************************************************************************************
 *  @file      source/lap.hpp
 *  @brief     A curve walked by length, piece by piece: a polyline's segments, a spline's spans, so
 *             that places at many lengths cost one pass, not a search each.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_LAP_HPP
#define GEOXY_LAP_HPP

#include "geoxy/accuracy.hpp"
#include "geoxy/curve.hpp"

#include <cstddef>
#include <vector>

namespace Geoxy
{
  /**
   *  A piece of a curve on a walk along it: the piece as a curve of its own, whose parameter plus a
   *  base is the curve's, and how far along the walk it starts and how long it is.
   */
  struct Lap
  {
    Curve  Shape  = {};  /**< The piece.                                  */
    double Base   = 0.0; /**< What its parameter is short of the curve's. */
    double From   = 0.0; /**< How far along the walk it starts.           */
    double Length = 0.0; /**< Its length.                                 */
  };

  /**
   *  A walk along a curve: its pieces in order.
   */
  using Walk = std::vector<Lap>;

  /**
   *  The walks along a curve: a polyline's segments, in one walk or each a walk of its own; a
   *  spline's spans in one walk; any other curve one lap; none for a curve without end, a polyline
   *  with no segment or a spline that is not valid.
   *
   *  @param  curve    The curve.
   *  @param  across   Whether a polyline's segments are one walk.
   *  @param  accuracy How far each lap's length may be off.
   *  @return          The walks.
   */
  auto Walks(const Curve &curve, bool across, Accuracy accuracy) -> std::vector<Walk>;

  /**
   *  How far a walk has come, to find places at lengths that never fall.
   */
  struct Walker
  {
    const Walk     *Laps     = nullptr; /**< The walk.                       */
    Geoxy::Accuracy Accuracy = {};      /**< How far a place may be off.     */
    std::size_t     Index    = 0;       /**< The lap the last place fell in. */
  };

  /**
   *  The curve's parameter at a length along a walk, no less than the last one asked.
   *
   *  @param  walker   The walker, moved on to the lap the length falls in.
   *  @param  distance The length, within the walk.
   *  @return          The parameter.
   */
  auto Trodden(Walker &walker, double distance) -> double;

  /**
   *  How long a walk is.
   *
   *  @param  walk The walk, not empty.
   *  @return      The length.
   */
  auto Long(const Walk &walk) -> double;
} // namespace Geoxy

#endif /* GEOXY_LAP_HPP */
