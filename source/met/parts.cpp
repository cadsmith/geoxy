/***************************************************************************************************
 *  @file      source/met/parts.cpp
 *  @brief     A curve cut into parts as the thorough way meets it, runs, elliptic arcs and rational
 *             Bézier pieces, each with how its parameter maps to the curve's; what meeting two
 *             parts finds; how far round a curve winds; its parts' bounds and its ends.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "met/parts.hpp"

#include "bow.hpp"
#include "elliptic.hpp"
#include "kernel/piece.hpp"
#include "scratch.hpp"
#include "straight.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy::Met
{
  auto Swapped(Found found) -> Found
  {
    for (Local &local : found.Locals)
    {
      std::swap(local.First, local.Second);
    }
    for (Overlap &overlap : found.Shared)
    {
      std::swap(overlap.First, overlap.Second);
    }
    return found;
  }

  auto Parted(const Curve &curve) -> List<Part>
  {
    List<Part> parts;
    if (const auto *line = std::get_if<Line>(&curve))
    {
      parts.push_back({.Shape = Running(*line)});
    }
    else if (const auto *ray = std::get_if<Ray>(&curve))
    {
      parts.push_back({.Shape = Running(*ray)});
    }
    else if (const auto *xline = std::get_if<Xline>(&curve))
    {
      parts.push_back({.Shape = Running(*xline)});
    }
    else if (const auto *arc = std::get_if<Arc>(&curve))
    {
      parts.push_back({.Shape = Elliptic(*arc)});
    }
    else if (const auto *ellipse = std::get_if<Ellipse>(&curve))
    {
      parts.push_back({.Shape = *ellipse});
    }
    else if (const auto *polyline = std::get_if<Polyline>(&curve))
    {
      const std::size_t count = Count(*polyline);
      for (std::size_t index = 0; index < count; ++index)
      {
        const auto    offset = static_cast<double>(index);
        const Vertex &from   = polyline->Vertices[index];
        const Point   until  = polyline->Vertices[(index + 1) % polyline->Vertices.size()].Position;
        const double  bulge  = from.Bulge;
        if (Remote(from.Position, until, bulge))
        {
          /**
           *  The arc as a rational piece: its middle pole where its ends' tangents meet, off the
           *  chord's middle by b / (1 - b^2) of the chord turned a quarter, its weight the cosine
           *  of half the sweep, (1 - b^2) / (1 + b^2).
           */
          const double square = bulge * bulge;
          const Vector chord  = until - from.Position;
          const Point  bend   = Midpoint(from.Position, until)
                              + (Vector{.X = chord.Y, .Y = -chord.X} * (bulge / (1.0 - square)));
          const double weight = (1.0 - square) / (1.0 + square);
          Piece        piece  = Kernel::Polynomial(std::array{from.Position, bend, until});
          piece.Poles[1]      = {.X = bend.X * weight, .Y = bend.Y * weight, .W = weight};
          parts.push_back({.Shape = piece, .Offset = offset, .Bend = bulge});
          continue;
        }
        const Segment segment = Nth(*polyline, index);
        if (const Line *straight = std::get_if<Line>(&segment))
        {
          parts.push_back({.Shape = Running(*straight), .Offset = offset});
        }
        else
        {
          parts.push_back({.Shape = Elliptic(std::get<Arc>(segment)), .Offset = offset});
        }
      }
      if (count == 0 && !polyline->Vertices.empty())
      {
        const Point lone = polyline->Vertices.front().Position;
        parts.push_back({.Shape = Run{.Start = lone, .Way = {}, .Along = {}}});
      }
    }
    else if (const auto *spline = std::get_if<Spline>(&curve))
    {
      Kernel::Sections(
          *spline,
          [&parts](const Piece &piece, Range span) -> void
          {
            parts.push_back({.Shape = piece, .Offset = span.Start, .Scale = span.End - span.Start});
          });
    }
    else if (const auto *quadratic = std::get_if<Bezier::Quadratic>(&curve))
    {
      parts.push_back({.Shape = Kernel::Polynomial(*quadratic)});
    }
    else
    {
      parts.push_back({.Shape = Kernel::Polynomial(std::get<Bezier::Cubic>(curve))});
    }
    return parts;
  }

  auto Wound(const Curve &curve) -> double
  {
    if (const auto *arc = std::get_if<Arc>(&curve))
    {
      return arc->Sweep.InRadians();
    }
    if (const auto *ellipse = std::get_if<Ellipse>(&curve))
    {
      return ellipse->Sweep.InRadians();
    }
    return 0.0;
  }

  auto Overwound(const Curve &curve) -> bool
  {
    return std::abs(Wound(curve)) > static_cast<double>(MostLaps) * Turn.InRadians();
  }

  auto Laps(const Ellipse &ellipse, double angle) -> List<double>
  {
    List<double> passes = Passes(ellipse, angle);
    if (passes.size() > MostLaps)
    {
      std::ranges::sort(passes);
      passes.resize(MostLaps);
    }
    return passes;
  }

  auto Bounds(const Part &part) -> Range
  {
    if (const Run *run = std::get_if<Run>(&part.Shape))
    {
      return run->Along;
    }
    return {.Start = 0.0, .End = 1.0};
  }

  auto Extremes(const Curve &curve) -> List<Point>
  {
    const Range domain   = Domain(curve);
    const auto *polyline = std::get_if<Polyline>(&curve);
    List<Point> ends;
    for (const double end : {domain.Start, domain.End})
    {
      if (std::isfinite(end))
      {
        ends.push_back(At(curve, end));
      }
    }
    if (polyline != nullptr)
    {
      return polyline->Closed ? List<Point>{} : ends;
    }
    /**
     *  An arc or an ellipse swept a whole number of turns is closed, as the direct ways take it,
     *  though its ends, worked out from a start in radians and a sweep in degrees, may lie a step
     *  of a double apart.
     */
    const auto *arc     = std::get_if<Arc>(&curve);
    const auto *ellipse = std::get_if<Ellipse>(&curve);
    if (arc != nullptr || ellipse != nullptr)
    {
      const double degrees = (arc != nullptr ? arc->Sweep : ellipse->Sweep).InDegrees();
      if (!(std::fpclassify(degrees)                               == FP_ZERO)
          && std::fpclassify(std::fmod(degrees, Turn.InDegrees())) == FP_ZERO)
      {
        return {};
      }
    }
    /**
     *  Only a curve whose ends are one point may close on itself: its box, which costs a pass over
     *  a spline's spans, is asked of no other.
     */
    if (ends.size() != 2 || !(ends[0] == ends[1]))
    {
      return ends;
    }
    const Vector size = Size(Extents(curve));
    return size.X > 0.0 || size.Y > 0.0 ? List<Point>{} : ends;
  }
} // namespace Geoxy::Met
