/***************************************************************************************************
 *  @file      source/seek.hpp
 *  @brief     What a point asks of a curve, found piece by piece: an ellipse's quarters and a
 *             polyline's segments solved and chained, and the places at the parameters found.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SEEK_HPP
#define GEOXY_SEEK_HPP

#include "eased.hpp"
#include "elliptic.hpp"
#include "kernel/piece.hpp"
#include "real.hpp"
#include "scratch.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  /**
   *  The nearest place on an arc to a point, as Nearest gives it, the arc's start and end known as
   *  At gives them: Nearest's steps on the same numbers, those points taken as given where it would
   *  work them out again.
   *
   *  @param  arc   The arc.
   *  @param  point The point.
   *  @param  ends  The arc's start and end, as At gives them at 0 and 1.
   *  @return       The place.
   */
  auto Nearest(const Arc &arc, Point point, const std::array<Point, 2> &ends) -> Place;

  /**
   *  A curve's places at parameters.
   *
   *  @tparam CurveT The curve's kind.
   *  @param  curve  The curve.
   *  @param  alongs The parameters, unordered, any repeated or a few steps of a double apart.
   *  @return        The places, in order, each once.
   */
  template <typename CurveT>
  auto Placed(const CurveT &curve, List<double> alongs) -> List<Place>
  {
    Gathered(alongs);
    List<Place> places;
    places.reserve(alongs.size());
    for (const double along : alongs)
    {
      places.push_back({.At = At(curve, along), .Along = along});
    }
    return places;
  }

  /**
   *  The most places a point asks of an ellipse in a turn: a conic has four feet of the
   *  perpendiculars from a point at most, and two places whose tangent passes through it.
   */
  inline constexpr std::size_t MostSpots = 4;

  /**
   *  How many spots in a turn mean noise: four times what a conic has. A solver on an ellipse thin
   *  or far beside its size finds a real place more than once, a step of a double apart, and which
   *  of those is truest no residual tells; past this many, the spots are a thin ellipse's rounding,
   *  and only MostSpots of them are kept.
   */
  inline constexpr std::size_t NoisySpots = 4 * MostSpots;

  /**
   *  No residual: every spot a solver finds is kept, as for the places whose tangent passes through
   *  a point, whose solver leaves no such noise.
   */
  struct Uncapped
  {
  };

  /**
   *  How far from answering a place of a curve is, as a measure of the way to it from a point and
   *  its tangent, both of length 1, tells it: 0 where it answers.
   *
   *  @tparam MeasureT What measures the two, 0 where the place answers.
   *  @param  point    The point.
   *  @param  measure  The measure.
   *  @return          What gives the residual of a place's point and tangent: the measure's size,
   *                   taken of no way where either is none, as where the point is on the curve.
   */
  template <typename MeasureT>
  auto Residual(Point point, const MeasureT &measure)
  {
    return [point, measure](Point place, Vector tangent) -> double
    {
      return std::abs(
          measure(Unit(place - point).value_or(Vector{}), Unit(tangent).value_or(Vector{})));
    };
  }

  /**
   *  What a solver finds on an ellipse's quarters, chained round, as the ellipse's parameters in
   *  every turn it sweeps, to its first 65,536.
   *
   *  A turn has MostSpots places at most. Beside an ellipse thin past what doubles hold, rounding
   *  leaves the solver many more, over a thousand in a turn, each then passed on every turn: past
   *  NoisySpots in a turn, the MostSpots whose point and tangent are nearest answering are kept, in
   *  the order they came, and the rest are noise. Fewer are all kept, as a real place found twice
   *  is one of them.
   *
   *  @tparam SolveT    What solves on a piece, giving a Kernel::Solution.
   *  @tparam ResidualT What tells how far a place's point and tangent are from answering.
   *  @param  ellipse   The ellipse.
   *  @param  solve     The solver.
   *  @param  residual  The residual, as Residual gives it; Uncapped keeps every spot.
   *  @return           The parameters, unordered.
   */
  template <typename SolveT, typename ResidualT>
  auto Rounded(const Ellipse &ellipse, const SolveT &solve, const ResidualT &residual)
      -> List<double>
  {
    std::array<Kernel::Piece, 4>    pieces;
    std::array<Kernel::Solution, 4> quarters;
    for (std::size_t index = 0; index < quarters.size(); ++index)
    {
      pieces.at(index)   = Kernel::Quarter(ellipse.Map, index);
      quarters.at(index) = solve(pieces.at(index));
    }

    /**
     *  A spot found and how far from answering it is, its angle on the unit circle, and where it
     *  came in the order the quarters gave.
     */
    struct Found
    {
      Kernel::Spot Spot     = {};  /**< Where on the quarters.     */
      double       Angle    = 0.0; /**< Its angle.                 */
      double       Residual = 0.0; /**< How far from answering.    */
      std::size_t  Order    = 0;   /**< Its place among the spots. */
    };

    /**
     *  The spots in the order the quarters give them, which is the order the places come in and so
     *  which of two as near Nearest keeps: a spot on a seam, found from both its sides, is one
     *  angle, kept where it first comes.
     */
    List<Found> found;
    for (const Kernel::Spot spot : Kernel::Chained(quarters, true))
    {
      const double angle = Kernel::Quartered(spot.Index, spot.Along);
      if (std::ranges::none_of(found, [angle](const Found &one) -> bool
                               { return Exactly(one.Angle, angle); }))
      {
        found.push_back({.Spot = spot, .Angle = angle, .Order = found.size()});
      }
    }
    if constexpr (!std::is_same_v<ResidualT, Uncapped>)
    {
      if (found.size() > NoisySpots)
      {
        for (Found &one : found)
        {
          const Kernel::Sample sample = Kernel::Sampled(pieces.at(one.Spot.Index), one.Spot.Along);
          one.Residual                = residual(sample.At, sample.Way);
        }
        /**
         *  Ties kept in the order they came, as each spot's order breaks them: the order a stable
         *  sort gives, with no room taken.
         */
        std::ranges::sort(
            found, [](const Found &lhs, const Found &rhs) -> bool
            { return std::pair{lhs.Residual, lhs.Order} < std::pair{rhs.Residual, rhs.Order}; });
        found.resize(MostSpots);
        std::ranges::sort(found, {}, &Found::Order);
      }
    }
    List<double> alongs;
    for (const Found &one : found)
    {
      const List<double> passes = Passes(ellipse, one.Angle);
      alongs.insert(alongs.end(), passes.begin(), passes.end());
    }
    return alongs;
  }

  /**
   *  What a solver finds on a polyline's segments, chained, as the polyline's parameters.
   *
   *  @tparam SolveT   What solves on a segment, giving a Kernel::Solution.
   *  @param  polyline The polyline.
   *  @param  solve    The solver.
   *  @return          The parameters, unordered.
   */
  template <typename SolveT>
  auto Segmented(const Polyline &polyline, const SolveT &solve) -> List<double>
  {
    List<Kernel::Solution> segments;
    const std::size_t      count = Count(polyline);
    segments.reserve(count);
    for (std::size_t index = 0; index < count; ++index)
    {
      segments.push_back(solve(Nth(polyline, index)));
    }
    List<double> alongs;
    for (const Kernel::Spot spot : Kernel::Chained(segments, polyline.Closed))
    {
      alongs.push_back(static_cast<double>(spot.Index) + spot.Along);
    }
    return alongs;
  }

  /**
   *  What a solver finds on the pieces of a spline's paving, chained: each place's point worked out
   *  on the piece it was found on, whose parameter tells it to a few steps of a double, and its
   *  parameter the spline's own, the double nearest.
   *
   *  Beside a heavy weight the spline's parameter tells a place coarsely or not at all, so the
   *  point is not the spline's at the parameter given, but the place found; places apart on the
   *  curve may so share a parameter.
   *
   *  @tparam SolveT What solves on a Kernel::Tread, giving a Kernel::Solution.
   *  @param  paving The paving.
   *  @param  solve  The solver.
   *  @return        The places, in order, one found from either side of a seam once.
   */
  template <typename SolveT>
  auto Stood(const Paving &paving, const SolveT &solve) -> List<Place>
  {
    List<Kernel::Solution> solutions;
    solutions.reserve(paving.Treads.size());
    for (const Kernel::Tread &tread : paving.Treads)
    {
      solutions.push_back(solve(tread));
    }
    List<Place>  places;
    Kernel::Spot last  = {};
    const double close = 64.0 * std::numeric_limits<double>::epsilon();
    for (const Kernel::Spot spot : Kernel::Chained(solutions, false))
    {
      const double apart = static_cast<double>(spot.Index - last.Index) + (spot.Along - last.Along);
      if (!places.empty() && apart <= close)
      {
        continue;
      }
      places.push_back({.At    = Kernel::Reached(paving.Treads[spot.Index].Shape, spot.Along),
                        .Along = Homed(paving.Pieces[spot.Index], spot.Along)});
      last = spot;
    }
    return places;
  }

  /**
   *  What is found on a rational spline: on its paving where its weights make it stiff, and on its
   *  spans otherwise.
   *
   *  @tparam TreadT What solves on a Kernel::Tread, a piece of the paving, giving a
   *                 Kernel::Solution.
   *  @tparam SpanT  What solves on a span of the spline, told by the index of the knot that starts
   *                 it, giving a Kernel::Solution.
   *  @param  spline The spline.
   *  @param  tread  The paving's solver.
   *  @param  span   The spans' solver.
   *  @return        The places, in order, each once.
   */
  template <typename TreadT, typename SpanT>
  auto Sought(const Spline &spline, const TreadT &tread, const SpanT &span) -> List<Place>
  {
    if (const std::optional<Paving> paving = Paved(spline); paving.has_value())
    {
      return Stood(*paving, tread);
    }
    return Placed(spline, Kernel::Stepped(spline, span));
  }

  /**
   *  A segment's solution for a point: roots found in closed form, and a measure of the way to the
   *  point and the segment's tangent at its ends, for the turns at its vertices.
   *
   *  @tparam MeasureT What measures the way from the point and the tangent.
   *  @param  segment  The segment.
   *  @param  point    The point.
   *  @param  measure  The measure.
   *  @return          No roots, and the measure at the segment's ends.
   */
  template <typename MeasureT>
  auto Leaned(const Segment &segment, Point point, const MeasureT &measure) -> Kernel::Solution
  {
    const auto measured = [&segment, point, &measure](double along) -> double
    {
      return std::visit([point, along, &measure](const auto &shape) -> double
                        { return measure(At(shape, along) - point, Derivative(shape, along)); },
                        segment);
    };
    return {.Roots = {}, .Start = measured(0.0), .End = measured(1.0)};
  }

  /**
   *  A list drawn from the thread's pool, copied out as a plain vector for the public interface.
   *
   *  @tparam ValueT What it holds.
   *  @param  list   The list.
   *  @return        Its copy.
   */
  template <typename ValueT>
  auto Vectored(const List<ValueT> &list) -> std::vector<ValueT>
  {
    return {list.begin(), list.end()};
  }
} // namespace Geoxy

#endif /* GEOXY_SEEK_HPP */
