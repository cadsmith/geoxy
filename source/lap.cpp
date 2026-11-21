/***************************************************************************************************
 *  @file      source/lap.cpp
 *  @brief     A curve walked by length, piece by piece.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "lap.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/length.hpp"
#include "geoxy/part.hpp"
#include "geoxy/place.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A lap added to the end of a walk, its start and length worked out.
     *
     *  @param walk     The walk.
     *  @param lap      The piece, and what its parameter is short of the curve's.
     *  @param accuracy How far its length may be off.
     */
    void Appended(Walk &walk, Lap lap, Accuracy accuracy)
    {
      lap.From   = walk.empty() ? 0.0 : Long(walk);
      lap.Length = Length(lap.Shape, accuracy);
      walk.push_back(std::move(lap));
    }

    /**
     *  The walks along a polyline: its segments.
     *
     *  @param  polyline The polyline.
     *  @param  across   Whether they are one walk.
     *  @param  accuracy Unused: a segment's length is exact.
     *  @return          The walks.
     */
    auto Walking(const Polyline &polyline, bool across, Accuracy accuracy) -> std::vector<Walk>
    {
      std::vector<Walk> walks;
      const std::size_t count = Count(polyline);
      for (std::size_t index = 0; index < count; ++index)
      {
        if (walks.empty() || !across)
        {
          walks.emplace_back();
        }
        const Curve shape =
            std::visit([](const auto &segment) -> Curve { return segment; }, Nth(polyline, index));
        Appended(walks.back(), {.Shape = shape, .Base = static_cast<double>(index)}, accuracy);
      }
      return walks;
    }

    /**
     *  The walk along a spline: its spans, each clamped.
     *
     *  @param  spline   The spline.
     *  @param  across   Unused: a spline is one walk.
     *  @param  accuracy How far each span's length may be off.
     *  @return          The walk; none for a spline that is not valid.
     */
    auto Walking(const Spline &spline, bool /*across*/, Accuracy accuracy) -> std::vector<Walk>
    {
      if (!Valid(spline))
      {
        return {};
      }
      const auto        degree = static_cast<std::size_t>(spline.Degree);
      const std::size_t poles  = spline.Poles.size();
      Walk              walk;
      for (std::size_t span = degree; span < poles; ++span)
      {
        const double start = spline.Knots[span];
        const double end   = spline.Knots[span + 1];
        if (start < end)
        {
          Appended(walk, {.Shape = Part(spline, Range{.Start = start, .End = end})}, accuracy);
        }
      }
      return {walk};
    }

    /**
     *  The walk along a ray or an xline: none, as it has no end.
     *
     *  @tparam ShapeT A ray or an xline.
     *  @return        No walk.
     */
    template <typename ShapeT>
    auto Nowhere(const ShapeT & /*shape*/) -> std::vector<Walk>
    {
      return {};
    }

    /**
     *  The walk along any other curve: the curve, one lap.
     *
     *  @tparam ShapeT   A line, an arc, an ellipse or a Bézier curve.
     *  @param  shape    The curve.
     *  @param  accuracy How far its length may be off.
     *  @return          The walk.
     */
    template <typename ShapeT>
    auto OneLap(const ShapeT &shape, Accuracy accuracy) -> std::vector<Walk>
    {
      Walk walk;
      Appended(walk, {.Shape = shape}, accuracy);
      return {walk};
    }
  } // namespace

  auto Walks(const Curve &curve, bool across, Accuracy accuracy) -> std::vector<Walk>
  {
    return std::visit(
        [across, accuracy](const auto &shape) -> std::vector<Walk>
        {
          using ShapeT = std::decay_t<decltype(shape)>;
          if constexpr (std::is_same_v<ShapeT, Ray> || std::is_same_v<ShapeT, Xline>)
          {
            return Nowhere(shape);
          }
          else if constexpr (std::is_same_v<ShapeT, Polyline> || std::is_same_v<ShapeT, Spline>)
          {
            return Walking(shape, across, accuracy);
          }
          else
          {
            return OneLap(shape, accuracy);
          }
        },
        curve);
  }

  auto Trodden(Walker &walker, double distance) -> double
  {
    const Walk &walk = *walker.Laps;
    while (walker.Index + 1 < walk.size()
           && distance      > walk[walker.Index].From + walk[walker.Index].Length)
    {
      ++walker.Index;
    }
    const Lap   &lap   = walk[walker.Index];
    const double local = std::clamp(distance - lap.From, 0.0, lap.Length);
    return lap.Base + PlaceAt(lap.Shape, local, walker.Accuracy).Along;
  }

  auto Long(const Walk &walk) -> double
  {
    return walk.back().From + walk.back().Length;
  }
} // namespace Geoxy
