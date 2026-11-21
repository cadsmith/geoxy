/***************************************************************************************************
 *  @file      source/parallel/smooth.hpp
 *  @brief     A curve that is neither line nor arc, an ellipse or a spline, laid with pairs of arcs
 *             tangent to its offset, within an accuracy: the raw offset's pieces, and the curve's
 *             own as what the offset is measured against.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_PARALLEL_SMOOTH_HPP
#define GEOXY_PARALLEL_SMOOTH_HPP

#include "kernel/piece.hpp"
#include "scratch.hpp"

#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"

#include <cstddef>
#include <span>
#include <variant>

namespace Geoxy::Parallel
{
  /**
   *  An arc laid along a curve's offset: its ends, center, radius and sweep.
   */
  struct Moved
  {
    Point  Start  = {};  /**< Where it starts.                            */
    Point  End    = {};  /**< Where it ends.                              */
    Point  Center = {};  /**< Its center.                                 */
    double Radius = 0.0; /**< Its radius.                                 */
    double Sweep  = 0.0; /**< Its sweep in radians, signed; 0 for a line. */
  };

  /**
   *  A stretch of an ellipse: its map, and where it starts and how far it runs on the unit circle.
   */
  struct Oval
  {
    Matrix Map   = {};  /**< The map of the unit circle.          */
    double Start = 0.0; /**< Where it starts, in radians.         */
    double Sweep = 0.0; /**< How far it runs, in radians, signed. */
  };

  /**
   *  A smooth stretch of a curve, its parameter from 0 to 1: an ellipse's, or a rational Bézier
   *  piece of a spline.
   */
  using Run = std::variant<Oval, Kernel::Piece>;

  /**
   *  A run's point and derivative.
   *
   *  @param  run   The run.
   *  @param  along Where, from 0 to 1.
   *  @return       The point and its derivative by the parameter.
   */
  auto Sampled(const Run &run, double along) -> Kernel::Sample;

  /**
   *  How a smooth curve is laid with arcs.
   */
  struct Laying
  {
    double      Distance = 0.0; /**< The distance, to the right when positive; 0 for the curve. */
    double      Accuracy = 0.0; /**< How near the arcs must come to the offset.                 */
    std::size_t Budget   = 0;   /**< The most arcs that may be laid.                            */
  };

  /**
   *  A smooth curve laid with arcs, run by run.
   */
  struct Layout
  {
    List<Moved> Pieces = {};    /**< The arcs, as moved pieces, in order.               */
    bool        Spent  = false; /**< Whether the budget ran out, the arcs not all laid. */
  };

  /**
   *  The offset of a smooth curve's runs laid with arcs: each run's offset point c(t) + d n(t), n
   *  the unit normal to the right, laid with pairs of arcs tangent to it at both ends, a pair
   *  through each stretch halved until, at a quarter, a half and three quarters of its way, the
   *  offset lies within the accuracy of the pair, each arc of which turns less than a quarter turn.
   *  A stretch the offset runs back along, past a cusp where the radius of curvature is the
   *  distance, is halved down to a chord.
   *
   *  @param  runs   The runs, in order.
   *  @param  laying The distance, the accuracy and the budget.
   *  @return        The arcs, and the curve's points they start from.
   */
  auto Laid(std::span<const Run> runs, const Laying &laying) -> Layout;
} // namespace Geoxy::Parallel

#endif /* GEOXY_PARALLEL_SMOOTH_HPP */
