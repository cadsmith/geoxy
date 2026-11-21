/***************************************************************************************************
 *  @file      tests/bench/shapes.hpp
 *  @brief     Curves of every kind drawn from fixed seeds, in a drawing ten units square.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TESTS_BENCH_SHAPES_HPP
#define GEOXY_TESTS_BENCH_SHAPES_HPP

#include "../geoxy/draws.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace Geoxy::Bench
{
  /**
   *  How many curves a benchmark goes round: too many for a branch predictor to learn by heart.
   */
  inline constexpr std::size_t Pool = 256;

  /**
   *  The drawing every curve lies in.
   */
  inline constexpr Box Sheet = {.Min = {.X = 0, .Y = 0}, .Max = {.X = 10, .Y = 10}};

  /**
   *  The tolerance of the drawing, as an application takes it.
   *
   *  @return The tolerance.
   */
  inline auto Within() -> Tolerance
  {
    return Tolerance::Of(Sheet);
  }

  /**
   *  A point anywhere in the drawing.
   *
   *  @param  draws The draws.
   *  @return       The point.
   */
  inline auto Anywhere(Draws &draws) -> Point
  {
    return {.X = draws.Near(5, 5), .Y = draws.Near(5, 5)};
  }

  /**
   *  A way of any heading, of a length near one.
   *
   *  @param  draws The draws.
   *  @return       The way.
   */
  inline auto Anyway(Draws &draws) -> Vector
  {
    const Vector way    = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)};
    const double length = std::sqrt((way.X * way.X) + (way.Y * way.Y));
    return way / length;
  }

  /**
   *  A line between two points anywhere.
   *
   *  @param  draws The draws.
   *  @return       The line.
   */
  inline auto AnyLine(Draws &draws) -> Curve
  {
    const Point start = Anywhere(draws);
    return Line{.Start = start, .End = Anywhere(draws)};
  }

  /**
   *  A ray from anywhere, any way.
   *
   *  @param  draws The draws.
   *  @return       The ray.
   */
  inline auto AnyRay(Draws &draws) -> Curve
  {
    const Point start = Anywhere(draws);
    return Ray{.Start = start, .Way = Anyway(draws)};
  }

  /**
   *  An endless line through anywhere, any way.
   *
   *  @param  draws The draws.
   *  @return       The endless line.
   */
  inline auto AnyXline(Draws &draws) -> Curve
  {
    const Point through = Anywhere(draws);
    return Xline{.Through = through, .Way = Anyway(draws)};
  }

  /**
   *  An arc in degrees, as drawings hold them: a radius from 1 to 4, a sweep from 30 to 300 degrees
   *  counterclockwise.
   *
   *  @param  draws The draws.
   *  @return       The arc.
   */
  inline auto AnyArc(Draws &draws) -> Curve
  {
    const Point  center = {.X = draws.Near(5, 2), .Y = draws.Near(5, 2)};
    const double radius = draws.Near(2.5, 1.5);
    const double start  = draws.Near(180, 180);
    return Arc{.Center = center,
               .Radius = radius,
               .Start  = Degrees(start),
               .Sweep  = Degrees(draws.Near(165, 135))};
  }

  /**
   *  An elliptic arc turned any way: a major radius from 1 to 4, a ratio of its radii from 0.1 to
   *  1, a sweep from a quarter of a turn to all but a quarter.
   *
   *  @param  draws The draws.
   *  @return       The elliptic arc.
   */
  inline auto AnyEllipse(Draws &draws) -> Curve
  {
    const Point  center = {.X = draws.Near(5, 2), .Y = draws.Near(5, 2)};
    const Vector major  = Anyway(draws) * draws.Near(2.5, 1.5);
    const double ratio  = draws.Near(0.55, 0.45);
    const double start  = draws.Near(Pi, Pi);
    return Ellipse{
        .Map   = {.X = major, .Y = {.X = -major.Y * ratio, .Y = major.X * ratio}, .Origin = center},
        .Start = Radians(start),
        .Sweep = Radians(draws.Near(Pi, Pi / 2))};
  }

  /**
   *  A quadratic Bézier curve of three points anywhere.
   *
   *  @param  draws The draws.
   *  @return       The curve.
   */
  inline auto AnyQuadratic(Draws &draws) -> Curve
  {
    const Point start   = Anywhere(draws);
    const Point control = Anywhere(draws);
    return Bezier::Quadratic{.Start = start, .Control = control, .End = Anywhere(draws)};
  }

  /**
   *  A cubic Bézier curve of four points anywhere, a loop or a cusp now and then.
   *
   *  @param  draws The draws.
   *  @return       The curve.
   */
  inline auto AnyCubic(Draws &draws) -> Curve
  {
    const Point start = Anywhere(draws);
    const Point leave = Anywhere(draws);
    const Point reach = Anywhere(draws);
    return Bezier::Cubic{.Start = start, .Out = leave, .In = reach, .End = Anywhere(draws)};
  }

  /**
   *  A cubic spline of eight poles from left to right, up and down at random, its knots even and
   *  clamped.
   *
   *  @param  draws The draws.
   *  @return       The spline.
   */
  inline auto AnySpline(Draws &draws) -> Curve
  {
    constexpr std::size_t Poles = 8;

    Spline spline = {.Degree = 3, .Knots = {0, 0, 0, 0}, .Poles = {}, .Weights = {}};
    for (std::size_t index = 0; index < Poles; ++index)
    {
      const double along = 10.0 * static_cast<double>(index) / static_cast<double>(Poles - 1);
      spline.Poles.push_back({.X = along + draws.Near(0, 0.5), .Y = draws.Near(5, 4)});
    }
    for (std::size_t index = 1; index < Poles - 3; ++index)
    {
      spline.Knots.push_back(static_cast<double>(index));
    }
    spline.Knots.insert(spline.Knots.end(), 4, static_cast<double>(Poles - 3));
    return spline;
  }

  /**
   *  A polyline zigzagging from left to right between heights 3 and 7, every third vertex a bulge
   *  of 0.2: a line across it at height 5 crosses every segment.
   *
   *  @param  vertices How many vertices, two at least.
   *  @return          The polyline.
   */
  inline auto Zigzag(std::size_t vertices) -> Polyline
  {
    Polyline polyline;
    for (std::size_t index = 0; index < vertices; ++index)
    {
      const double along = 10.0 * static_cast<double>(index) / static_cast<double>(vertices - 1);
      polyline.Vertices.push_back({.Position = {.X = along, .Y = (index % 2) == 0 ? 3.0 : 7.0},
                                   .Bulge    = (index % 3) == 0 ? 0.2 : 0.0});
    }
    return polyline;
  }

  /**
   *  A polyline of sixteen vertices zigzagging from left to right, its heights and bulges at
   *  random, a third of its segments arcs.
   *
   *  @param  draws The draws.
   *  @return       The polyline.
   */
  inline auto AnyPolyline(Draws &draws) -> Curve
  {
    Polyline polyline = Zigzag(16);
    for (Vertex &vertex : polyline.Vertices)
    {
      vertex.Position.Y += draws.Near(0, 1);
      vertex.Bulge      *= draws.Near(0, 2);
    }
    return polyline;
  }

  /**
   *  A line across the drawing from left to right, near its middle height.
   *
   *  @param  draws The draws.
   *  @return       The line.
   */
  inline auto Across(Draws &draws) -> Curve
  {
    const double start = draws.Near(5, 1);
    return Line{.Start = {.X = 0, .Y = start}, .End = {.X = 10, .Y = draws.Near(5, 1)}};
  }

  /**
   *  A pool of curves drawn from a seed.
   *
   *  @param  seed The seed.
   *  @param  make How to draw one.
   *  @return      The curves.
   */
  template <typename MakeT>
  auto Drawn(std::uint64_t seed, MakeT make) -> std::vector<Curve>
  {
    Draws              draws{seed};
    std::vector<Curve> curves;
    curves.reserve(Pool);
    for (std::size_t index = 0; index < Pool; ++index)
    {
      curves.push_back(make(draws));
    }
    return curves;
  }
} // namespace Geoxy::Bench

#endif /* GEOXY_TESTS_BENCH_SHAPES_HPP */
