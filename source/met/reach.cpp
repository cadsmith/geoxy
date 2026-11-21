/***************************************************************************************************
 *  @file      source/met/reach.cpp
 *  @brief     Whether two curves may meet at all, told by their boxes before anything is worked
 *             out, and the length they are met within.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "met/scale.hpp"
#include "real.hpp"
#include "thorough.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <variant>

namespace Geoxy::Met
{
  namespace
  {
    /**
     *  The box some points are hemmed in by, and so all that lies within their hull: a Bézier curve
     *  by its control points, a spline of weights above nothing by its poles.
     *
     *  @param  points The points.
     *  @return        The box, exact, empty for no point; nothing where a coordinate is not finite,
     *                 a sum of them past every double taken as that too.
     */
    auto Hemmed(std::span<const Point> points) -> std::optional<Box>
    {
      Box    box;
      double sum = 0.0;
      for (const Point point : points)
      {
        box  = Grown(box, point);
        sum += point.X + point.Y;
      }
      return std::isfinite(sum) ? std::optional<Box>{box} : std::nullopt;
    }

    /**
     *  A box a polyline cannot leave: about its vertices, grown by how far its arcs bow out.
     *
     *  An arc of a bulge lies within its chord's box grown by its height, the bulge times half the
     *  chord: up to a half turn it stays that near its chord, and past it that near the chord's
     *  middle. The chord is taken no shorter than it is, as the sum of its two extents. The arc
     *  held, of a center and a radius rounded and angles rounded, strays from that by a few steps
     *  of a double at the size of its center and radius, sixteen here, which is much where a slight
     *  bulge makes them large.
     *
     *  @param  polyline The polyline.
     *  @return          The box, of edges not finite for no vertex; nothing where a number is not
     *                   finite.
     */
    auto Strung(const Polyline &polyline) -> std::optional<Box>
    {
      const std::span<const Vertex> vertices = polyline.Vertices;
      Box                           box;
      double                        sum    = 0.0;
      double                        height = 0.0;
      double                        radius = 0.0;
      for (std::size_t index = 0; index < vertices.size(); ++index)
      {
        const Vertex &vertex  = vertices[index];
        const Point   next    = vertices[(index + 1) % vertices.size()].Position;
        const Vector  chord   = next - vertex.Position;
        const double  bulge   = std::abs(vertex.Bulge);
        box                   = Grown(box, vertex.Position);
        sum                  += vertex.Position.X + vertex.Position.Y + bulge;
        if (bulge >= Geoxy::Straight)
        {
          const double reach = std::abs(chord.X) + std::abs(chord.Y);
          height             = std::max(height, bulge * reach / 2.0);
          radius             = std::max(radius, reach * (1.0 + (bulge * bulge)) / (4.0 * bulge));
        }
      }
      if (!std::isfinite(sum))
      {
        return std::nullopt;
      }
      const double size = std::max(
          {std::abs(box.Min.X), std::abs(box.Max.X), std::abs(box.Min.Y), std::abs(box.Max.Y)});
      const double grow = height + (16.0 * Above(2.0 * (size + radius)));
      return Box{.Min = {.X = box.Min.X - grow, .Y = box.Min.Y - grow},
                 .Max = {.X = box.Max.X + grow, .Y = box.Max.Y + grow}};
    }

    /**
     *  A line's box, between its ends, or an arc's, about its whole circle whatever its sweep, read
     *  off its numbers with nothing worked out: the commonest curves, kept apart to be told far
     *  from each other in a few steps.
     *
     *  @param  curve The curve.
     *  @return       The box, a line's exact, each edge of an arc's within half a step of a double
     *                of its circle's; nothing for any other curve, nor where a number of it is not
     *                finite, a sum of them past every double taken as that too.
     */
    auto Read(const Curve &curve) -> std::optional<Box>
    {
      if (const auto *line = std::get_if<Line>(&curve))
      {
        const Point start = line->Start;
        const Point end   = line->End;
        if (!std::isfinite(start.X + start.Y + end.X + end.Y))
        {
          return std::nullopt;
        }
        return Box{.Min = {.X = std::min(start.X, end.X), .Y = std::min(start.Y, end.Y)},
                   .Max = {.X = std::max(start.X, end.X), .Y = std::max(start.Y, end.Y)}};
      }
      if (const auto *arc = std::get_if<Arc>(&curve))
      {
        const Point  center = arc->Center;
        const double radius = std::abs(arc->Radius);
        if (!std::isfinite(center.X + center.Y + radius))
        {
          return std::nullopt;
        }
        return Box{.Min = {.X = center.X - radius, .Y = center.Y - radius},
                   .Max = {.X = center.X + radius, .Y = center.Y + radius}};
      }
      return std::nullopt;
    }

    /**
     *  A box a curve cannot leave, read off its numbers with little worked out: a line's or an
     *  arc's as Read gives it; about an ellipse's whole oval, whatever its sweep; about a Bézier
     *  curve's control points, a sound spline's poles, a polyline's vertices and arcs.
     *
     *  Nothing for a curve that runs on without end, nor where a number of it is not finite, which
     *  the size of the curves, and so their grain, leaves out. A sum past every double is taken as
     *  not finite too, and costs only the box.
     *
     *  @param  curve The curve.
     *  @return       The box: exact for a line, a Bézier curve and a spline, which their numbers'
     *                rounding carries a few steps of a double past it; within a step of a double
     *                for an arc and an ellipse.
     */
    auto Boxed(const Curve &curve) -> std::optional<Box>
    {
      if (const std::optional<Box> read = Read(curve))
      {
        return read;
      }
      if (const auto *polyline = std::get_if<Polyline>(&curve))
      {
        return Strung(*polyline);
      }
      if (const auto *cubic = std::get_if<Bezier::Cubic>(&curve))
      {
        return Hemmed(Given(*cubic));
      }
      if (const auto *quadratic = std::get_if<Bezier::Quadratic>(&curve))
      {
        return Hemmed(Given(*quadratic));
      }
      if (const auto *spline = std::get_if<Spline>(&curve))
      {
        return Valid(*spline) ? Hemmed(spline->Poles) : std::nullopt;
      }
      if (const auto *ellipse = std::get_if<Ellipse>(&curve))
      {
        const Matrix &map    = ellipse->Map;
        const Vector  reach  = {.X = std::abs(map.X.X) + std::abs(map.Y.X),
                                .Y = std::abs(map.X.Y) + std::abs(map.Y.Y)};
        const Point   center = map.Origin;
        if (!std::isfinite(center.X + center.Y + reach.X + reach.Y))
        {
          return std::nullopt;
        }
        return Box{.Min = center - reach, .Max = center + reach};
      }
      return std::nullopt;
    }

  } // namespace
} // namespace Geoxy::Met

namespace Geoxy
{
  auto Confines(const Curve &curve) -> std::optional<Box>
  {
    return Met::Boxed(curve);
  }

  /**
   *  The length two curves are met within, where they may meet at all: asked before anything is
   *  worked out, as most pairs of a drawing lie far apart.
   *
   *  Every hit and every shared stretch has a point of each curve within the tolerance of a point
   *  of the other, and a band is weighed within four. The boxes' edges of a line or an arc are off
   *  by half a step of a double at the most, and their gap by as much again, where the tolerance is
   *  sixteen steps at the least: lines and arcs whose boxes lie more than four tolerances apart
   *  have no two points within three, and so meet in nothing, whatever the thorough way would
   *  weigh. Any other curve's points are worked out by sums of up to thirty-two terms, each
   *  rounded, a few steps past its box: it is kept sixteen tolerances clear.
   *
   *  The boxes give the grain of lines and arcs too: the largest magnitude of their edges is the
   *  largest of a line's coordinates, and of an arc's center's with its radius added, to the bit,
   *  so the length is the one Grained gives, the curves' numbers read once.
   *
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance asked for.
   *  @return           The length used, no finer than the curves' grain; nothing where neither runs
   *                    on without end and their boxes lie more than four times that length apart,
   *                    sixteen where either is not a line or an arc, along x or along y.
   */
  auto Reachable(const Curves &curves, const Tolerance &tolerance) -> std::optional<double>
  {
    std::optional<Box> one   = Met::Read(curves[0]);
    std::optional<Box> other = one ? Met::Read(curves[1]) : std::nullopt;
    const bool         plain = one && other;
    if (!plain)
    {
      one   = Met::Boxed(curves[0]);
      other = one ? Met::Boxed(curves[1]) : std::nullopt;
    }
    if (one && other)
    {
      const double size =
          std::max({std::abs(one->Min.X), std::abs(one->Max.X), std::abs(one->Min.Y),
                    std::abs(one->Max.Y), std::abs(other->Min.X), std::abs(other->Max.X),
                    std::abs(other->Min.Y), std::abs(other->Max.Y)});
      /**
       *  An arc's box may run past every double where its center and radius do not: its size is
       *  then left to Grained, which leaves such sums out.
       */
      if (std::isfinite(size))
      {
        /**
         *  Two lines or arcs are read once, their size their boxes'; any other curve's is its own,
         *  and it is kept sixteen tolerances clear, not four, as its points are worked out by sums
         *  of up to thirty-two terms, each rounded.
         */
        const double length = plain ? std::max(tolerance.Length, 16.0 * Above(size))
                                    : Met::Grained(curves[0], curves[1], tolerance).Length;
        const double clear  = (plain ? 4.0 : 16.0) * length;
        const bool   remote = one->Min.X - other->Max.X > clear
                           || other->Min.X - one->Max.X > clear
                           || one->Min.Y - other->Max.Y > clear
                           || other->Min.Y - one->Max.Y > clear;
        return remote ? std::nullopt : std::optional<double>{length};
      }
    }
    return Met::Grained(curves[0], curves[1], tolerance).Length;
  }
} // namespace Geoxy
