/***************************************************************************************************
 *  @file      source/met/parts.hpp
 *  @brief     A curve cut into parts as the thorough way meets it, runs, elliptic arcs and rational
 *             Bézier pieces, each with how its parameter maps to the curve's; what meeting two
 *             parts finds; how far round a curve winds; its parts' bounds and its ends.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_MET_PARTS_HPP
#define GEOXY_MET_PARTS_HPP

#include "kernel/piece.hpp"
#include "scratch.hpp"
#include "straight.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <variant>

namespace Geoxy::Met
{
  using Kernel::Piece;

  /**
   *  A part of a curve: a straight run, an elliptic arc or a rational Bézier piece, and how its own
   *  parameter maps to the curve's.
   */
  struct Part
  {
    std::variant<Run, Ellipse, Piece> Shape  = Run{}; /**< What the part is.                      */
    double                            Offset = 0.0;   /**< The curve's parameter at the part's 0. */
    double                            Scale  = 1.0;   /**< The curve's for each 1 of the part's.  */
    double                            Bend   = 0.0;   /**< The bulge of a polyline's slight arc
                                                           held as a piece; else 0.               */
  };

  /**
   *  The curve's parameter at a part's own.
   *
   *  A slight arc held as a piece runs by the piece's parameter t, and its curve by the share of
   *  the sweep, which is 1 / 2 + atan(b (2 t - 1)) / (2 atan b).
   *
   *  @param  part  The part.
   *  @param  along The part's parameter.
   *  @return       The curve's.
   */
  inline auto Placed(const Part &part, double along) -> double
  {
    if (std::equal_to<double>{}(part.Bend, 0.0))
    {
      return part.Offset + (part.Scale * along);
    }
    const double share =
        0.5 + (std::atan(part.Bend * ((2.0 * along) - 1.0)) / (2.0 * std::atan(part.Bend)));
    return part.Offset + (part.Scale * share);
  }

  /**
   *  A candidate meeting of two parts, in their own parameters.
   */
  struct Local
  {
    double First  = 0.0; /**< The first part's parameter. */
    double Second = 0.0; /**< The second's.               */
  };

  /**
   *  What meeting two parts found, in their own parameters.
   */
  struct Found
  {
    List<Local>   Locals = {};    /**< The candidate meetings.   */
    List<Overlap> Shared = {};    /**< The stretches they share. */
    bool          Spent  = false; /**< Whether a budget ran out. */
  };

  /**
   *  What two parts found, the parts' roles swapped.
   *
   *  @param  found What was found with the parts one way round.
   *  @return       The same, the other way round.
   */
  auto Swapped(Found found) -> Found;

  /**
   *  A curve cut into parts.
   *
   *  @param  curve The curve.
   *  @return       Its parts, in order.
   */
  auto Parted(const Curve &curve) -> List<Part>;

  /**
   *  The most turns of an arc or an ellipse swept round many times that a meeting looks at, the
   *  first along it; and the most pairs of turns of two such curves at one point. Every turn that
   *  passes a crossing is a hit, two such curves meet once for every pair of their turns, and
   *  telling which hits are one costs the square of their count: past these the meeting is cut
   *  short, and says so by Spent.
   */
  inline constexpr std::size_t MostLaps  = 1024;
  inline constexpr std::size_t MostPairs = 2048;

  /**
   *  How far round an arc or an ellipse sweeps, in radians: nothing for any other curve.
   *
   *  @param  curve The curve.
   *  @return       Its sweep.
   */
  auto Wound(const Curve &curve) -> double;

  /**
   *  Whether a curve is swept round more turns than a meeting looks at.
   *
   *  @param  curve The curve.
   *  @return       Whether it is.
   */
  auto Overwound(const Curve &curve) -> bool;

  /**
   *  Where an ellipse passes an angle of its unit circle, on the first MostLaps turns that do.
   *
   *  @param  ellipse The ellipse.
   *  @param  angle   The angle, in radians, any turn of it.
   *  @return         Its parameters there, in order for one swept round more than MostLaps turns.
   */
  auto Laps(const Ellipse &ellipse, double angle) -> List<double>;

  /**
   *  The stretch of a part's own parameter.
   *
   *  @param  part The part.
   *  @return      A run's stretch; 0 to 1 for an arc or a piece.
   */
  auto Bounds(const Part &part) -> Range;

  /**
   *  A curve's ends, unless it closes on itself.
   *
   *  @param  curve The curve.
   *  @return       Its ends that are finite; none for a closed polyline, or for any other curve
   *                whose ends are one point but that is not a point itself, as a whole arc, a whole
   *                ellipse or a loop of a spline is.
   */
  auto Extremes(const Curve &curve) -> List<Point>;
} // namespace Geoxy::Met

#endif /* GEOXY_MET_PARTS_HPP */
