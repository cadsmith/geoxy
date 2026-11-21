/***************************************************************************************************
 *  @file      source/swift/lane.hpp
 *  @brief     A straight curve as the swift kernels take it, what they find on it, and the kernels
 *             that meet it with an ellipse, a polyline or a Bézier curve.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_LANE_HPP
#define GEOXY_SWIFT_LANE_HPP

#include "scratch.hpp"
#include "straight.hpp"
#include "sure.hpp"
#include "swift/section.hpp"
#include "trig.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <optional>

namespace Geoxy::Swift
{
  /**
   *  A line, a ray or an endless line as the swift kernels take it: the point P where its parameter
   *  is 0 and its way d, exact in two doubles, the stretch of its parameter and its ends.
   */
  struct Lane
  {
    Point                       Start = {};      /**< P, where its parameter is 0.             */
    std::array<Sure::Number, 2> Way   = {};      /**< d, exactly, a line's end less its start. */
    Vector                      Rough = {};      /**< d rounded, the double nearest it.        */
    double                      Speed = 0.0;     /**< |d| as doubles work it out.              */
    Range                       Along = {};      /**< Its parameter's stretch.                 */
    std::array<Point, 2>        Ends  = {};      /**< Its ends.                                */
    std::size_t                 Count = 0;       /**< How many ends it has.                    */
    const Curve                *Shape = nullptr; /**< The curve it is.                         */
    double Lead = 0.0; /**< How far the place of the curve it stands for runs ahead of its own:
                            a polyline's segment's index, added where a crossing's place is
                            rounded, so that it is the double nearest that curve's place. */
  };

  /**
   *  A lane's place at a crossing as the curve it stands for places it: its own, its lead added
   *  exactly where it has one.
   *
   *  @param  lane  The lane.
   *  @param  along Its own place, with its bound.
   *  @return       The place, with its bound.
   */
  inline auto Led(const Lane &lane, const Sure::Number &along) -> Sure::Number
  {
    return lane.Lead > 0.0 ? along + Sure::Exact(lane.Lead) : along;
  }

  /**
   *  A crossing a swift kernel found: its point, the double nearest the true one, and its place on
   *  each curve, the straight one's first.
   */
  struct Spot
  {
    Point                 At    = {};    /**< The point.                              */
    std::array<double, 2> Along = {};    /**< Its places, the straight curve's first. */
    bool                  Own   = false; /**< Whether the places are the curves' own, each
                                              the double nearest the true one; else they are
                                              rough, to choose a segment by, and the point
                                              is placed as the exact crossing places it.     */
  };

  /**
   *  The crossings a swift kernel found.
   */
  using Spots = List<Spot>;

  /**
   *  The farthest from zero any number a kernel takes may lie, so that nothing it works out, of
   *  four factors at the most, runs past the doubles: 2^200.
   */
  inline constexpr double Widest = 0x1p200;

  /**
   *  The nearest to zero a length a kernel takes may lie, so that nothing it works out falls among
   *  the doubles below the normal ones where a bound must hold: 2^-200.
   */
  inline constexpr double Narrowest = 0x1p-200;

  /**
   *  How far past the reach of eight tolerances a kernel keeps every decision about the tolerance:
   *  a quarter more, many steps of a double past where the thorough way's own roundings could put
   *  the curves, its tolerance being sixteen such steps at the least.
   */
  inline constexpr double Margin = 1.25;

  /**
   *  How far apart the curves must be wherever they do not cross, and every end from the other
   *  curve: eight tolerances, and the margin.
   *
   *  @param  tolerance The tolerance.
   *  @return           The length.
   */
  inline auto Reach(const Tolerance &tolerance) -> double
  {
    return 8.0 * tolerance.Length * Margin;
  }

  /**
   *  Whether a crossing is steep: the sine of the angle between the curves there above twice the
   *  tolerance's, with room for rounding, so that the thorough way finds the curves pass to either
   *  side of each other there, and finds no band about it.
   *
   *  @param  sine      The sine, as worked out.
   *  @param  tolerance The tolerance.
   *  @return           Whether it is.
   */
  inline auto Steep(double sine, const Tolerance &tolerance) -> bool
  {
    return sine > (2.0 * Trig::Sin(tolerance.Radians)) + 0x1p-30;
  }

  /**
   *  Where a value lies in a stretch.
   */
  enum class Side : std::uint8_t
  {
    Out,  /**< Outside, by more than the margin. */
    In,   /**< Inside, by more than the margin.  */
    Near, /**< Within the margin of an end.      */
  };

  /**
   *  Where a value lies in a stretch, to within a margin.
   *
   *  @param  range  The stretch, either end infinite or not.
   *  @param  value  The value.
   *  @param  margin The margin.
   *  @return        Where it lies; near where the margin or the value is not a finite size.
   */
  inline auto Within(Range range, double value, double margin) -> Side
  {
    if (!(margin >= 0.0) || !std::isfinite(margin) || !std::isfinite(value))
    {
      return Side::Near;
    }
    if (value < range.Start - margin || value > range.End + margin)
    {
      return Side::Out;
    }
    if (value - range.Start <= margin || range.End - value <= margin)
    {
      return Side::Near;
    }
    return Side::In;
  }

  /**
   *  Whether numbers lie within the range the kernels take.
   *
   *  @param  values The numbers.
   *  @return        Whether every one is no farther from zero than Widest; never where one is not a
   *                 number.
   */
  inline auto Held(std::initializer_list<double> values) -> bool
  {
    return std::ranges::all_of(values,
                               [](double value) -> bool { return std::abs(value) <= Widest; });
  }

  /**
   *  A line, a ray or an endless line as the kernels take it.
   *
   *  @param  curve The curve.
   *  @return       It as a lane; nothing for any other curve, one of no length, or one whose
   *                numbers lie outside the range the kernels take.
   */
  auto Laned(const Curve &curve) -> std::optional<Lane>;

  /**
   *  How far a point is from a lane, at the least: from its carrier, or, where its foot lies past
   *  an end, from that end, each less what rounding may take off.
   *
   *  @param  lane  The lane.
   *  @param  point The point.
   *  @return       The distance, at the least.
   */
  auto Away(const Lane &lane, Point point) -> double;

  /**
   *  Whether two runs certainly do not lie along each other as the thorough way's Aligned finds it:
   *  the stretch of the first the second's projection covers no longer than half the tolerance, or
   *  an end of it more than twice the tolerance from the second's line.
   *
   *  @param  first     One run.
   *  @param  second    The other, not parallel to it.
   *  @param  tolerance The tolerance.
   *  @return           Whether they certainly do not.
   */
  auto Unaligned(const Run &first, const Run &second, const Tolerance &tolerance) -> bool;

  /**
   *  A line, a ray or an endless line met with an ellipse, as the thorough way meets them, where
   *  that leaves nothing to decide: the carrier misses the ellipse by more than eight tolerances,
   *  or cuts into it deeper than that, every crossing far from either curve's ends and steeper than
   *  twice the tolerance's angle; each crossing the double nearest the true point.
   *
   *  @param  lane      The straight curve.
   *  @param  ellipse   The ellipse.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @return           The crossings, roughly placed; nothing where the thorough way must decide.
   */
  auto Oval(const Lane &lane, const Ellipse &ellipse, const Tolerance &tolerance)
      -> std::optional<Spots>;

  /**
   *  How far a point is from a polyline, at the least: from each segment, straight or bulged, as
   *  the kernel of a lane and a polyline weighs its ends; endless where it has no segment of any
   *  length.
   *
   *  @param  polyline The polyline, its numbers within the range the kernels take.
   *  @param  point    The point.
   *  @return          The distance, at the least.
   */
  auto Shunned(const Polyline &polyline, Point point) -> double;

  /**
   *  A line, a ray or an endless line met with a polyline, as the thorough way meets them, where
   *  that leaves nothing to decide: every vertex and every end far from the other curve, each
   *  segment missing the carrier by more than eight tolerances or crossing it steeply, a bulge's
   *  circle cut deeper than that; each crossing the double nearest the true point.
   *
   *  @param  lane      The straight curve.
   *  @param  polyline  The polyline.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @return           The crossings, roughly placed; nothing where the thorough way must decide.
   */
  auto Bulged(const Lane &lane, const Polyline &polyline, const Tolerance &tolerance)
      -> std::optional<Spots>;

  /**
   *  A line, a ray or an endless line met with a quadratic or cubic piece, a Bézier curve or a
   *  polynomial spline's span, as the thorough way meets them, where that leaves nothing to decide:
   *  the piece far from the carrier but about crossings, each a simple root of the carrier's
   *  equation on the piece, steep and far from either curve's ends and the piece's; each crossing
   *  the double nearest the true point, its places the doubles nearest the true ones. The piece's
   *  control points are charged their slack wherever they are weighed, and a crossing is proved on
   *  the true piece.
   *
   *  @param  lane      The straight curve.
   *  @param  section   The piece.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @return           The crossings, placed, the curve's place its own parameter; nothing where
   *                    the thorough way must decide.
   */
  auto Bernstein(const Lane &lane, const Section &section, const Tolerance &tolerance)
      -> std::optional<Spots>;
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_LANE_HPP */
