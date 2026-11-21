/***************************************************************************************************
 *  @file      source/parallel/bend.hpp
 *  @brief     A path's segment held for the questions an offset asks of it many times over: its
 *             box, its ways at its ends, its points along it and its length, each in a few
 *             operations of doubles.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_PARALLEL_BEND_HPP
#define GEOXY_PARALLEL_BEND_HPP

#include "trig.hpp"

#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>

namespace Geoxy::Parallel
{
  /**
   *  A path's segment: a line, or an arc by its bulge, read as Bent reads it.
   *
   *  What is worked out here only decides which pairs are met, and where a slice is sampled: every
   *  crossing is Intersect's, on the segment as Bent makes it.
   */
  struct Bend
  {
    Point  Start  = {};  /**< Where it starts.                        */
    Point  End    = {};  /**< Where it ends.                          */
    double Bulge  = 0.0; /**< Its bulge, 0 for a line.                */
    double Sweep  = 0.0; /**< Its sweep in radians, four atan(bulge). */
    Point  Center = {};  /**< An arc's center.                        */
    double Radius = 0.0; /**< An arc's radius.                        */
  };

  /**
   *  A path's segment.
   *
   *  @param  vertices The path's vertices.
   *  @param  index    The segment, from its vertex to the next, the last to the first.
   *  @return          The segment: a line where the bulge is under Straight or its ends are one
   *                   point, as Bent has it.
   */
  inline auto Bent(std::span<const Vertex> vertices, std::size_t index) -> Bend
  {
    const Point  start = vertices[index].Position;
    const Point  end   = vertices[(index + 1) % vertices.size()].Position;
    const double bulge = vertices[index].Bulge;
    if (!(std::abs(bulge) >= Straight) || start == end)
    {
      return {.Start = start, .End = end};
    }
    const Vector chord  = end - start;
    const double square = bulge * bulge;
    return {
        .Start  = start,
        .End    = end,
        .Bulge  = bulge,
        .Sweep  = 4.0 * Trig::Atan(bulge),
        .Center = Midpoint(start, end) + (Perpendicular(chord) * ((1.0 - square) / (4.0 * bulge))),
        .Radius = Length(chord) * ((1.0 + square) / (4.0 * std::abs(bulge))),
    };
  }

  /**
   *  Whether a segment is an arc.
   *
   *  @param  bend The segment.
   *  @return      Whether it is.
   */
  inline auto Round(const Bend &bend) -> bool
  {
    return bend.Radius > 0.0;
  }

  /**
   *  How far a point lies from a segment.
   *
   *  @param  bend  The segment.
   *  @param  point The point.
   *  @return       The distance to the nearest point of the segment.
   */
  inline auto Apart(const Bend &bend, Point point) -> double
  {
    if (!Round(bend))
    {
      const Vector way    = bend.End - bend.Start;
      const double square = std::max(Dot(way, way), std::numeric_limits<double>::min());
      const double along  = std::clamp(Dot(point - bend.Start, way) / square, 0.0, 1.0);
      return Distance(point, bend.Start + (way * along));
    }
    const Vector off   = point - bend.Center;
    const double side  = bend.Sweep < 0.0 ? -1.0 : 1.0;
    const double after = Cross(bend.Start - bend.Center, off) * side;
    const double until = Cross(off, bend.End - bend.Center) * side;
    const bool   upon  = std::abs(bend.Sweep) <= std::numbers::pi ? std::min(after, until) >= 0.0
                                                                  : std::max(after, until) >= 0.0;
    if (upon)
    {
      return std::abs(Length(off) - bend.Radius);
    }
    return std::min(Distance(point, bend.Start), Distance(point, bend.End));
  }

} // namespace Geoxy::Parallel

#endif /* GEOXY_PARALLEL_BEND_HPP */
