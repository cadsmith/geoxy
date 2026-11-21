/***************************************************************************************************
 *  @file      source/met/scale.cpp
 *  @brief     The size of a curve's numbers and the grain a meeting is held to; curves too large or
 *             too small, or of ways too long or too short, met at a size where they are not.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "met/scale.hpp"

#include "bow.hpp"
#include "crossing.hpp"
#include "met/linear.hpp"
#include "real.hpp"
#include "thorough.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <type_traits>
#include <variant>

namespace Geoxy::Met
{
  auto Given(const Line &shape) -> std::array<Point, 2>
  {
    return {shape.Start, shape.End};
  }

  auto Given(const Ray &shape) -> std::array<Point, 1>
  {
    return {shape.Start};
  }

  auto Given(const Xline &shape) -> std::array<Point, 1>
  {
    return {shape.Through};
  }

  auto Given(const Spline &shape) -> std::span<const Point>
  {
    return shape.Poles;
  }

  auto Given(const Bezier::Quadratic &shape) -> std::array<Point, 3>
  {
    return {shape.Start, shape.Control, shape.End};
  }

  auto Given(const Bezier::Cubic &shape) -> std::array<Point, 4>
  {
    return {shape.Start, shape.Out, shape.In, shape.End};
  }

  namespace
  {
    /**
     *  The largest magnitude of some points' coordinates, each grown by a reach, finite ones only.
     *
     *  @param  points The points.
     *  @param  reach  How far the curve may run from each.
     *  @return        The magnitude, 0 for none.
     */
    auto Largest(std::span<const Point> points, double reach) -> double
    {
      double size = 0.0;
      for (const Point point : points)
      {
        for (const double value : {std::abs(point.X) + reach, std::abs(point.Y) + reach})
        {
          size = std::isfinite(value) ? std::max(size, value) : size;
        }
      }
      return size;
    }

    /**
     *  The size of a polyline's numbers: the largest magnitude of its vertices and of each bulged
     *  segment's circle, its center's grown by its radius, as each segment met as an arc of its own
     *  is sized. A bulge near Straight puts the center far off, its grain far coarser than the
     *  vertices'.
     *
     *  @param  polyline The polyline.
     *  @return          The size, 0 for none.
     */
    auto Spanned(const Polyline &polyline) -> double
    {
      const std::span<const Vertex> vertices = polyline.Vertices;
      const std::size_t             segments = Count(polyline);
      double                        size     = 0.0;
      for (std::size_t index = 0; index < vertices.size(); ++index)
      {
        const Point  start = vertices[index].Position;
        const Point  end   = vertices[(index + 1) % vertices.size()].Position;
        const double bulge = vertices[index].Bulge;
        size               = std::max(size, Largest(std::array{start}, 0.0));
        if (index < segments
            && std::abs(bulge) >= Geoxy::Straight
            && start != end
            && !Remote(start, end, bulge))
        {
          const Hoop hoop = Hooped(start, end, bulge);
          size            = std::max(size, Largest(std::array{hoop.Center}, std::abs(hoop.Radius)));
        }
      }
      return size;
    }

  } // namespace

  auto Size(const Curve &curve) -> double
  {
    return std::visit(
        [](const auto &shape) -> double
        {
          using Shape = std::decay_t<decltype(shape)>;
          if constexpr (std::is_same_v<Shape, Arc>)
          {
            return Largest(std::array{shape.Center}, std::abs(shape.Radius));
          }
          else if constexpr (std::is_same_v<Shape, Ellipse>)
          {
            return Largest(std::array{shape.Map.Origin}, Length(shape.Map.X) + Length(shape.Map.Y));
          }
          else if constexpr (std::is_same_v<Shape, Polyline>)
          {
            return Spanned(shape);
          }
          else
          {
            return Largest(Given(shape), 0.0);
          }
        },
        curve);
  }

  namespace
  {
    /**
     *  A point with both its coordinates times a power of two: exactly, but for one that leaves the
     *  doubles at either end.
     *
     *  @param  point The point.
     *  @param  power The power.
     *  @return       The point scaled.
     */
    auto Scaled(Point point, int power) -> Point
    {
      return {.X = std::ldexp(point.X, power), .Y = std::ldexp(point.Y, power)};
    }

    /**
     *  A vector with both its coordinates times a power of two.
     *
     *  @param  vector The vector.
     *  @param  power  The power.
     *  @return        The vector scaled.
     */
    auto Scaled(Vector vector, int power) -> Vector
    {
      return {.X = std::ldexp(vector.X, power), .Y = std::ldexp(vector.Y, power)};
    }

    /**
     *  A curve with every length of it times a power of two, about the origin: its points, its ways
     *  and its radius, and no angle, bulge, knot or weight. A place along it stays the place it
     *  was.
     *
     *  @param  curve The curve.
     *  @param  power The power.
     *  @return       The curve scaled.
     */
    auto Scaled(const Curve &curve, int power) -> Curve
    {
      const auto placed = [power](Point point) -> Point
      {
        return Scaled(point, power);
      };
      const auto headed = [power](Vector vector) -> Vector
      {
        return Scaled(vector, power);
      };
      if (const auto *line = std::get_if<Line>(&curve))
      {
        return Line{.Start = placed(line->Start), .End = placed(line->End)};
      }
      if (const auto *ray = std::get_if<Ray>(&curve))
      {
        return Ray{.Start = placed(ray->Start), .Way = headed(ray->Way)};
      }
      if (const auto *xline = std::get_if<Xline>(&curve))
      {
        return Xline{.Through = placed(xline->Through), .Way = headed(xline->Way)};
      }
      if (const auto *arc = std::get_if<Arc>(&curve))
      {
        return Arc{.Center = placed(arc->Center),
                   .Radius = std::ldexp(arc->Radius, power),
                   .Start  = arc->Start,
                   .Sweep  = arc->Sweep};
      }
      if (const auto *ellipse = std::get_if<Ellipse>(&curve))
      {
        return Ellipse{.Map   = {.X      = headed(ellipse->Map.X),
                                 .Y      = headed(ellipse->Map.Y),
                                 .Origin = placed(ellipse->Map.Origin)},
                       .Start = ellipse->Start,
                       .Sweep = ellipse->Sweep};
      }
      if (const auto *polyline = std::get_if<Polyline>(&curve))
      {
        Polyline scaled = *polyline;
        for (Vertex &vertex : scaled.Vertices)
        {
          vertex.Position = placed(vertex.Position);
        }
        return scaled;
      }
      if (const auto *spline = std::get_if<Spline>(&curve))
      {
        Spline scaled = *spline;
        for (Point &pole : scaled.Poles)
        {
          pole = placed(pole);
        }
        return scaled;
      }
      if (const auto *quadratic = std::get_if<Bezier::Quadratic>(&curve))
      {
        return Bezier::Quadratic{.Start   = placed(quadratic->Start),
                                 .Control = placed(quadratic->Control),
                                 .End     = placed(quadratic->End)};
      }
      const auto &cubic = std::get<Bezier::Cubic>(curve);
      return Bezier::Cubic{.Start = placed(cubic.Start),
                           .Out   = placed(cubic.Out),
                           .In    = placed(cubic.In),
                           .End   = placed(cubic.End)};
    }

    /**
     *  Past this, or under one over it, curves are met at a size near 1: 2^96, well within the
     *  2^100 either way up to which every pair of kinds was measured to meet right as it is.
     */
    constexpr double Vast = 0x1p96;

    /**
     *  The way an endless curve runs: a ray's or an endless line's, which is theirs to choose a
     *  length for; nothing for a curve with an end at each side.
     *
     *  @param  curve The curve.
     *  @return       Its way.
     */
    auto Heads(Curve &curve) -> Vector *
    {
      if (auto *ray = std::get_if<Ray>(&curve))
      {
        return &ray->Way;
      }
      auto *xline = std::get_if<Xline>(&curve);
      return xline != nullptr ? &xline->Way : nullptr;
    }

    /**
     *  The way an endless curve runs, read: a ray's or an endless line's; nothing for a curve with
     *  an end at each side.
     *
     *  @param  curve The curve.
     *  @return       Its way.
     */
    auto Headed(const Curve &curve) -> const Vector *
    {
      if (const auto *ray = std::get_if<Ray>(&curve))
      {
        return &ray->Way;
      }
      const auto *xline = std::get_if<Xline>(&curve);
      return xline != nullptr ? &xline->Way : nullptr;
    }

    /**
     *  The power of two that brings an endless curve's way to a length near 1, where the way is so
     *  long or so short that its length squared is past every double or nothing: 0 for a way that
     *  is not, for one that is nothing or not finite, and for any other curve.
     *
     *  @param  curve The curve.
     *  @return       The power.
     */
    auto Pace(const Curve &curve) -> int
    {
      const Vector *way = Headed(curve);
      if (way == nullptr)
      {
        return 0;
      }
      const double reach = std::max(std::abs(way->X), std::abs(way->Y));
      return Outscaled(reach) ? -std::ilogb(reach) : 0;
    }

    /**
     *  An endless curve's way times a power of two, the curve itself the same: a place along it is
     *  then that many times less.
     *
     *  @param curve The curve.
     *  @param pace  The power.
     */
    void Repace(Curve &curve, int pace)
    {
      if (Vector *way = Heads(curve); way != nullptr)
      {
        *way = Scaled(*way, pace);
      }
    }

    /**
     *  A curve met scaled, scaled back: its endless way's pace undone, then the power.
     *
     *  @param  curve The curve as met.
     *  @param  power The power it was scaled by.
     *  @param  pace  The power its way was scaled by besides.
     *  @return       The curve, as given where no number was lost below the doubles.
     */
    auto Unscaled(Curve curve, int power, int pace) -> Curve
    {
      Repace(curve, -pace);
      return Scaled(curve, -power);
    }

    /**
     *  Two curves met scaled: as given and as met, the power they were scaled by, whether every
     *  number of them was kept, and whether both are straight.
     */
    struct Scaling
    {
      const Curves *Given    = nullptr; /**< The curves as given.           */
      const Curves *Met      = nullptr; /**< The curves as met, scaled.     */
      int           Power    = 0;       /**< The power they were scaled by. */
      bool          Kept     = false;   /**< Whether no number was lost.    */
      bool          Straight = false;   /**< Whether both are straight.     */
    };

    /**
     *  A hit of two curves met scaled that lies at their scaled carriers' crossing, that crossing
     *  worked out again where scaling lost a bit, of a number or of the point: from the curves as
     *  met, scaled up with every number kept, and rounded once at the size they were given, as
     *  where a crossing lies below the normal doubles; else from two straight curves as given, a
     *  number of them lost, as the exact crossing finds them at any size.
     *
     *  @param  scaling How the curves were met.
     *  @param  hit     The hit, as met.
     *  @pre            Kept, the curves were scaled up.
     *  @return         The crossing at the size given; nothing for a hit elsewhere, or for curves
     *                  not straight that lost a number, which no crossing tells.
     */
    auto Recrossed(const Scaling &scaling, const Hit &hit) -> std::optional<Point>
    {
      const Curves             &given = *scaling.Given;
      const Curves             &met   = *scaling.Met;
      const Crossings::Approach near = {.Along = {hit.First.Along, hit.Second.Along}, .At = hit.At};
      const bool crossing = Crossings::Crossing({&met[0].get(), &met[1].get()}, near.Along, near.At)
                                .transform([&hit](const Crossings::Crossed &found) -> bool
                                           { return found.At == hit.At; })
                                .value_or(false);
      if (!crossing)
      {
        return std::nullopt;
      }
      if (scaling.Kept)
      {
        return Crossings::Spotted({&met[0].get(), &met[1].get()}, near, scaling.Power);
      }
      return scaling.Straight ? Crossings::Spotted({&given[0].get(), &given[1].get()}, near, 0)
                              : std::nullopt;
    }
  } // namespace

  auto Grained(const Curve &first, const Curve &second, Tolerance tolerance) -> Tolerance
  {
    const double size = std::max(Size(first), Size(second));
    const double step = Above(size);
    tolerance.Length  = std::max(tolerance.Length, 16.0 * step);
    return tolerance;
  }
} // namespace Geoxy::Met

namespace Geoxy
{
  auto Usable(const Curves &curves, const Tolerance &tolerance) -> Tolerance
  {
    return Met::Grained(curves[0], curves[1], tolerance);
  }

  auto Sized(const Curve &curve) -> double
  {
    return Met::Size(curve);
  }

  auto Outscaled(double size) -> bool
  {
    return std::isfinite(size) && size > 0.0 && (size > Met::Vast || size < 1.0 / Met::Vast);
  }

  auto Shrunk(const Curves &curves, const Tolerance &tolerance, Intersections &into) -> bool
  {
    const double size  = std::max(Met::Size(curves[0]), Met::Size(curves[1]));
    const int    power = Outscaled(size) ? -std::ilogb(size) : 0;
    /**
     *  At a power of nothing a copy is the curve itself, its way's pace its own: both are read off
     *  the curves as they are, and nothing is copied where nothing is scaled.
     */
    if (power == 0 && Met::Pace(curves[0]) == 0 && Met::Pace(curves[1]) == 0)
    {
      return false;
    }
    Curve                    first  = Met::Scaled(curves[0], power);
    Curve                    second = Met::Scaled(curves[1], power);
    const std::array<int, 2> paces  = {Met::Pace(first), Met::Pace(second)};
    Met::Repace(first, paces[0]);
    Met::Repace(second, paces[1]);
    Intersect(Curves{first, second},
              {.Length = std::ldexp(tolerance.Length, power), .Radians = tolerance.Radians}, into);
    const auto paced = [](Range range, int pace) -> Range
    {
      return {.Start = std::ldexp(range.Start, pace), .End = std::ldexp(range.End, pace)};
    };
    /**
     *  A number of the curves lost below the doubles, as a coordinate a sliver of the others' size
     *  scaled down, or a point that loses a bit scaled back, as a crossing below the normal
     *  doubles: the meeting is of curves a hair from these, its points a step off theirs. A
     *  crossing of two straight curves is worked out again on the curves as given, exactly, and
     *  where a number was lost each curve's point is taken again at its place.
     */
    const bool         kept    = Met::Unscaled(first, power, paces[0])  == curves[0].get()
                              && Met::Unscaled(second, power, paces[1]) == curves[1].get();
    const Curves       scaled  = {first, second};
    const Met::Scaling scaling = {.Given    = &curves,
                                  .Met      = &scaled,
                                  .Power    = power,
                                  .Kept     = kept,
                                  .Straight = Met::Straightened(curves[0]).has_value()
                                           && Met::Straightened(curves[1]).has_value()};
    for (Hit &hit : into.Hits)
    {
      const Point back = Met::Scaled(hit.At, -power);
      const bool  lost = !kept || (power > 0 && !(Met::Scaled(back, power) == hit.At));
      const std::optional<Point> exact = lost ? Met::Recrossed(scaling, hit) : std::nullopt;
      hit.At                           = exact.value_or(back);
      hit.First.At                     = Met::Scaled(hit.First.At, -power);
      hit.Second.At                    = Met::Scaled(hit.Second.At, -power);
      hit.First.Along                  = std::ldexp(hit.First.Along, paces[0]);
      hit.Second.Along                 = std::ldexp(hit.Second.Along, paces[1]);
      if (exact && !kept)
      {
        hit.First.At  = At(curves[0], hit.First.Along);
        hit.Second.At = At(curves[1], hit.Second.Along);
      }
    }
    for (Overlap &overlap : into.Overlaps)
    {
      overlap.First  = paced(overlap.First, paces[0]);
      overlap.Second = paced(overlap.Second, paces[1]);
    }
    return true;
  }
} // namespace Geoxy
