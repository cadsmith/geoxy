/***************************************************************************************************
 *  @file      source/at.cpp
 *  @brief     The point of each curve at a value of its parameter.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/at.hpp"

#include "basis.hpp"
#include "bow.hpp"
#include "real.hpp"
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
#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <variant>

namespace Geoxy
{
  auto At(const Line &line, double along) -> Point
  {
    return Lerp(line.Start, line.End, along);
  }

  auto At(const Ray &ray, double along) -> Point
  {
    return ray.Start + (ray.Way * along);
  }

  auto At(const Xline &xline, double along) -> Point
  {
    return xline.Through + (xline.Way * along);
  }

  auto At(const Arc &arc, double along) -> Point
  {
    return arc.Center + (Direction(arc.Start + (arc.Sweep * along)) * arc.Radius);
  }

  auto At(const Ellipse &ellipse, double along) -> Point
  {
    return ellipse.Map * (Point{} + Direction(ellipse.Start + (ellipse.Sweep * along)));
  }

  auto At(const Polyline &polyline, double along) -> Point
  {
    if (polyline.Vertices.empty())
    {
      return {};
    }
    if (Count(polyline) == 0)
    {
      return polyline.Vertices.front().Position;
    }
    const Step step = Locate(polyline, along);
    if (Straightness(step))
    {
      return Lerp(step.Start, step.End, step.Fraction);
    }
    if (Exactly(step.Fraction, 0.0))
    {
      return step.Start;
    }
    if (Exactly(step.Fraction, 1.0))
    {
      return step.End;
    }
    const Bight bow = Bowed(step.Start, step.End, step.Bulge);
    return Along(bow, bow.Sweep * step.Fraction);
  }

  auto At(const Spline &spline, double along) -> Point
  {
    const std::optional<std::size_t> span = Legible(spline, along);
    if (!span)
    {
      return {};
    }
    const int  power = Leveling(spline, *span);
    const auto poles = [&spline, power](std::size_t index) -> Homogeneous
    {
      return Pole(spline, index, power);
    };
    const Basis                   basis = Spanned(spline);
    std::array<double, MaxDegree> values{};
    const std::span<double>       taken = std::span{values}.first(basis.Degree);
    std::ranges::fill(taken, along);
    return Projected(Weighted(spline, Blossom(basis, poles, {.Span = *span, .Values = taken})));
  }

  auto At(const Bezier::Quadratic &quadratic, double along) -> Point
  {
    const Point start = Lerp(quadratic.Start, quadratic.Control, along);
    const Point end   = Lerp(quadratic.Control, quadratic.End, along);
    return Lerp(start, end, along);
  }

  auto At(const Bezier::Cubic &cubic, double along) -> Point
  {
    const Point first  = Lerp(cubic.Start, cubic.Out, along);
    const Point second = Lerp(cubic.Out, cubic.In, along);
    const Point third  = Lerp(cubic.In, cubic.End, along);
    return Lerp(Lerp(first, second, along), Lerp(second, third, along), along);
  }

  auto At(const Curve &curve, double along) -> Point
  {
    return std::visit([along](const auto &shape) -> Point { return At(shape, along); }, curve);
  }
} // namespace Geoxy
