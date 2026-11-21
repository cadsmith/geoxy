/***************************************************************************************************
 *  @file      source/nearest.cpp
 *  @brief     A curve's places nearest a point: straight runs and arcs in closed form, every other
 *             curve by the kernel's feet on its pieces.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/nearest.hpp"

#include "basis.hpp"
#include "elliptic.hpp"
#include "hull.hpp"
#include "kernel/piece.hpp"
#include "near/feet.hpp"
#include "near/plumb.hpp"
#include "near/turns.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "secular.hpp"
#include "seek.hpp"
#include "straight.hpp"
#include "walk.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/place.hpp"
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
#include <limits>
#include <optional>
#include <span>
#include <variant>
#include <vector>

namespace Geoxy
{
  auto Perpendiculars(const Line &line, Point point) -> std::vector<Place>
  {
    return Vectored(Near::Listed(line, point));
  }

  auto Perpendiculars(const Ray &ray, Point point) -> std::vector<Place>
  {
    return Vectored(Near::Listed(ray, point));
  }

  auto Perpendiculars(const Xline &xline, Point point) -> std::vector<Place>
  {
    return Vectored(Near::Listed(xline, point));
  }

  auto Perpendiculars(const Arc &arc, Point point) -> std::vector<Place>
  {
    return Vectored(Near::Listed(arc, point));
  }

  auto Perpendiculars(const Ellipse &ellipse, Point point) -> std::vector<Place>
  {
    return Vectored(Near::Listed(ellipse, point));
  }

  auto Perpendiculars(const Polyline &polyline, Point point) -> std::vector<Place>
  {
    return Vectored(Near::Listed(polyline, point));
  }

  auto Perpendiculars(const Spline &spline, Point point) -> std::vector<Place>
  {
    return Vectored(Near::Listed(spline, point));
  }

  auto Perpendiculars(const Bezier::Quadratic &quadratic, Point point) -> std::vector<Place>
  {
    return Vectored(Near::Listed(quadratic, point));
  }

  auto Perpendiculars(const Bezier::Cubic &cubic, Point point) -> std::vector<Place>
  {
    return Vectored(Near::Listed(cubic, point));
  }

  auto Perpendiculars(const Curve &curve, Point point) -> std::vector<Place>
  {
    return std::visit([point](const auto &shape) -> std::vector<Place>
                      { return Perpendiculars(shape, point); }, curve);
  }

  auto Nearest(const Line &line, Point point) -> Place
  {
    /**
     *  The start, the foot and the end, in order along the line: the first of the nearest, as Least
     *  would take it, with nothing gathered.
     */
    Place  best  = {.At = Lerp(line.Start, line.End, 0.0), .Along = 0.0};
    double least = Distance(point, best.At);
    if (const std::optional<double> along = Footing(Running(line), point))
    {
      const Place foot = {.At = Lerp(line.Start, line.End, *along), .Along = *along};
      if (Distance(point, foot.At) < least)
      {
        best  = foot;
        least = Distance(point, foot.At);
      }
    }
    const Place end = {.At = Lerp(line.Start, line.End, 1.0), .Along = 1.0};
    return Distance(point, end.At) < least ? end : best;
  }

  auto Nearest(const Ray &ray, Point point) -> Place
  {
    const Place start = {.At = ray.Start, .Along = 0.0};
    if (const std::optional<double> along = Footing(Running(ray), point))
    {
      const Place foot = {.At = At(ray, *along), .Along = *along};
      return Distance(point, foot.At) < Distance(point, start.At) ? foot : start;
    }
    return start;
  }

  auto Nearest(const Xline &xline, Point point) -> Place
  {
    if (const std::optional<double> along = Footing(Running(xline), point))
    {
      return {.At = At(xline, *along), .Along = *along};
    }
    return {.At = xline.Through, .Along = 0.0};
  }

  auto Nearest(const Arc &arc, Point point) -> Place
  {
    return Near::Endpoint(arc, point, nullptr);
  }

  auto Nearest(const Arc &arc, Point point, const std::array<Point, 2> &ends) -> Place
  {
    return Near::Endpoint(arc, point, &ends);
  }

  auto Nearest(const Ellipse &ellipse, Point point) -> Place
  {
    if (Coiled(ellipse))
    {
      /**
       *  Swept round more than once it passes every point of itself within its first turn, and the
       *  first of the nearest along it is there: that turn's nearest place, at as much of the whole
       *  sweep as the turn is, with none of the later passes listed or weighed.
       */
      const double  sweep = ellipse.Sweep.InRadians();
      const Ellipse once  = {
          .Map = ellipse.Map, .Start = ellipse.Start, .Sweep = sweep > 0.0 ? Turn : -Turn};
      const Place  first = Nearest(once, point);
      const double along = first.Along * (Turn.InRadians() / std::abs(sweep));
      return {.At = At(ellipse, along), .Along = along};
    }
    if (const std::optional<Place> place = Closest(ellipse, point))
    {
      return *place;
    }
    return Near::Bounded(ellipse, point);
  }

  auto Nearest(const Polyline &polyline, Point point) -> Place
  {
    const std::size_t             count    = Count(polyline);
    const std::span<const Vertex> vertices = polyline.Vertices;
    if (count == 0)
    {
      return {.At = At(polyline, 0.0), .Along = 0.0};
    }
    /**
     *  The vertices first, at their own positions, which At gives back where every number is
     *  finite: the nearest of them is as far as the answer can be.
     */
    Place  best  = {.At = vertices.front().Position, .Along = 0.0};
    double least = Distance(point, best.At);
    double sum   = 0.0;
    double most  = 0.0;
    for (std::size_t index = 0; index <= count; ++index)
    {
      const Point  corner  = vertices[index % vertices.size()].Position;
      const double away    = Distance(point, corner);
      sum                 += corner.X + corner.Y;
      most                 = std::max(most, Largest(corner - Point{}));
      if (away < least)
      {
        best  = {.At = corner, .Along = static_cast<double>(index)};
        least = away;
      }
    }
    if (!std::isfinite(sum))
    {
      return Near::Everywhere(polyline, point);
    }
    /**
     *  Then the feet of each segment that may come as near as the nearest place so far, the rest
     *  farther than it or off it (Bypassed, then the box). A foot within a few steps of a double of
     *  a vertex may be the one of two found from either side of it that Everywhere keeps once, or
     *  not, as their neighbours fall: such a point is weighed as Everywhere weighs it, every
     *  segment solved.
     */
    const double close  = 256.0
                        * std::numeric_limits<double>::epsilon()
                        * std::max(1.0, static_cast<double>(count) + 1.0);
    const double size   = Largest(point - Point{}) + (2.0 * most);
    const auto   inside = [close](std::span<double> alongs, double offset) -> bool
    {
      for (double &along : alongs)
      {
        along            = offset + along;
        const bool clear = offset + close < along && along < offset + 1.0 - close;
        if (!clear)
        {
          return false;
        }
      }
      return true;
    };
    const auto weigh = [&best, &least, point](const Place &place) -> void
    {
      const double away = Distance(point, place.At);
      if (away < least || (Exactly(away, least) && place.Along < best.Along))
      {
        best  = place;
        least = away;
      }
    };
    for (std::size_t index = 0; index < count; ++index)
    {
      const std::size_t next = index + 1 < vertices.size() ? index + 1 : 0;
      if (Near::Bypassed(vertices[index], vertices[next].Position,
                         {.From = point, .Least = least, .Size = size})
          || Outlying(Hull(polyline, index), point) > least * (1.0 + 0x1p-40))
      {
        continue;
      }
      Near::Brace             brace = Near::Braced(Nth(polyline, index), point);
      const std::span<double> feet  = std::span{brace.Alongs}.first(brace.Count);
      const auto              start = static_cast<double>(index);
      if (!inside(feet, start))
      {
        return Near::Everywhere(polyline, point);
      }
      std::array<Place, 2> places = {};
      for (std::size_t foot = 0; foot < brace.Count; ++foot)
      {
        const Step step = {.Start    = vertices[index].Position,
                           .End      = vertices[next].Position,
                           .Bulge    = vertices[index].Bulge,
                           .Fraction = feet[foot] - start};
        places.at(foot) = {.At    = Near::Trod(step, point, brace.Sides.at(foot)),
                           .Along = feet[foot]};
      }
      const std::span<Place> placed = std::span{places}.first(brace.Count);
      Near::Lineup(placed);
      std::ranges::for_each(placed, weigh);
    }
    return best;
  }

  auto Nearest(const Spline &spline, Point point) -> Place
  {
    if (!Valid(spline))
    {
      return {};
    }
    if (const std::optional<Place> place =
            Rational(spline) ? Near::Rationally(spline, point) : Near::Pruned(spline, point))
    {
      return *place;
    }
    List<Place> places = Near::Listed(spline, point);
    Kernel::Sections(spline,
                     [&places, &spline](const Kernel::Piece & /*piece*/, Range span) -> void
                     {
                       places.push_back({.At = At(spline, span.Start), .Along = span.Start});
                       places.push_back({.At = At(spline, span.End), .Along = span.End});
                     });
    return Near::Least(point, places);
  }

  auto Nearest(const Bezier::Quadratic &quadratic, Point point) -> Place
  {
    return Near::Bounded(quadratic, point);
  }

  auto Nearest(const Bezier::Cubic &cubic, Point point) -> Place
  {
    return Near::Bounded(cubic, point);
  }

  auto Nearest(const Curve &curve, Point point) -> Place
  {
    return std::visit([point](const auto &shape) -> Place { return Nearest(shape, point); }, curve);
  }
} // namespace Geoxy
