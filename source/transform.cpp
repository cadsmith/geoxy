/***************************************************************************************************
 *  @file      source/transform.cpp
 *  @brief     Each curve under an affine map, its kind kept where the map allows.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/transform.hpp"

#include "chain.hpp"
#include "walk.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <cstddef>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether every segment of a polyline is a line, as Bent tells it.
     *
     *  @param  polyline The polyline.
     *  @return          Whether it is.
     */
    auto Unbent(const Polyline &polyline) -> bool
    {
      const std::size_t segments = Count(polyline);
      for (std::size_t index = 0; index < segments; ++index)
      {
        if (!Straightness(Locate(polyline, static_cast<double>(index))))
        {
          return false;
        }
      }
      return true;
    }

    /**
     *  Whether a map is a similarity geoxy keeps an arc an arc under: its y axis its x axis turned
     *  a quarter turn either way, exactly. Every map built of turns, mirrors, uniform scales and
     *  moves, and every product of them, is one to the bit, as composing keeps that shape.
     *
     *  @param  map The map.
     *  @return     1 for a turn and a scale, -1 for one that mirrors too, 0 for any other map.
     */
    auto Similar(const Matrix &map) -> int
    {
      if (map.Y == Perpendicular(map.X))
      {
        return 1;
      }
      return map.Y == -Perpendicular(map.X) ? -1 : 0;
    }
  } // namespace

  auto Transformed(const Line &line, const Matrix &map) -> Line
  {
    return {.Start = map * line.Start, .End = map * line.End};
  }

  auto Transformed(const Ray &ray, const Matrix &map) -> Ray
  {
    return {.Start = map * ray.Start, .Way = map * ray.Way};
  }

  auto Transformed(const Xline &xline, const Matrix &map) -> Xline
  {
    return {.Through = map * xline.Through, .Way = map * xline.Way};
  }

  auto Transformed(const Arc &arc, const Matrix &map) -> Curve
  {
    const int similar = Similar(map);
    if (similar == 0)
    {
      const Matrix circle = {
          .X = {.X = arc.Radius, .Y = 0.0}, .Y = {.X = 0.0, .Y = arc.Radius}, .Origin = arc.Center};
      return Ellipse{.Map = map * circle, .Start = arc.Start, .Sweep = arc.Sweep};
    }
    const Angle turn   = Heading(map.X);
    const Arc   mapped = {.Center = map * arc.Center,
                          .Radius = arc.Radius * Length(map.X),
                          .Start  = similar > 0 ? arc.Start + turn : turn - arc.Start,
                          .Sweep  = similar > 0 ? arc.Sweep : -arc.Sweep};
    return mapped;
  }

  auto Transformed(const Ellipse &ellipse, const Matrix &map) -> Ellipse
  {
    return {.Map = map * ellipse.Map, .Start = ellipse.Start, .Sweep = ellipse.Sweep};
  }

  auto Transformed(const Polyline &polyline, const Matrix &map) -> Curve
  {
    const int similar = Similar(map);
    if (similar == 0 && !Unbent(polyline))
    {
      return Transformed(Joined(Sliced(polyline)), map);
    }
    Polyline mapped = polyline;
    for (Vertex &vertex : mapped.Vertices)
    {
      vertex.Position = map * vertex.Position;
      vertex.Bulge    = similar < 0 ? -vertex.Bulge : vertex.Bulge;
    }
    return mapped;
  }

  auto Transformed(const Spline &spline, const Matrix &map) -> Spline
  {
    Spline mapped = spline;
    std::ranges::transform(mapped.Poles, mapped.Poles.begin(),
                           [&map](Point pole) -> Point { return map * pole; });
    return mapped;
  }

  auto Transformed(const Bezier::Quadratic &quadratic, const Matrix &map) -> Bezier::Quadratic
  {
    return {.Start   = map * quadratic.Start,
            .Control = map * quadratic.Control,
            .End     = map * quadratic.End};
  }

  auto Transformed(const Bezier::Cubic &cubic, const Matrix &map) -> Bezier::Cubic
  {
    return {.Start = map * cubic.Start,
            .Out   = map * cubic.Out,
            .In    = map * cubic.In,
            .End   = map * cubic.End};
  }

  auto Transformed(const Curve &curve, const Matrix &map) -> Curve
  {
    return std::visit([&map](const auto &shape) -> Curve { return Transformed(shape, map); },
                      curve);
  }
} // namespace Geoxy
