/***************************************************************************************************
 *  @file      source/derivative.cpp
 *  @brief     The derivative of each curve with respect to its parameter.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/derivative.hpp"

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

#include <cstddef>
#include <optional>
#include <span>
#include <variant>

namespace Geoxy
{
  auto Derivative(const Line &line, double /*along*/) -> Vector
  {
    return line.End - line.Start;
  }

  auto Derivative(const Ray &ray, double /*along*/) -> Vector
  {
    return ray.Way;
  }

  auto Derivative(const Xline &xline, double /*along*/) -> Vector
  {
    return xline.Way;
  }

  auto Derivative(const Arc &arc, double along) -> Vector
  {
    const Vector way = Direction(arc.Start + (arc.Sweep * along));
    return Perpendicular(way) * (arc.Radius * arc.Sweep.InRadians());
  }

  auto Derivative(const Ellipse &ellipse, double along) -> Vector
  {
    const Vector way = Direction(ellipse.Start + (ellipse.Sweep * along));
    return ellipse.Map * (Perpendicular(way) * ellipse.Sweep.InRadians());
  }

  auto Derivative(const Polyline &polyline, double along) -> Vector
  {
    if (Count(polyline) == 0)
    {
      return {};
    }
    const Step step = Locate(polyline, along);
    if (Straightness(step))
    {
      return step.End - step.Start;
    }
    const Bight bow = Bowed(step.Start, step.End, step.Bulge);
    return Tangent(bow, bow.Sweep * step.Fraction);
  }

  auto Derivative(const Spline &spline, double along) -> Vector
  {
    const std::optional<std::size_t> read = Legible(spline, along);
    if (!read)
    {
      return {};
    }
    const int                     power  = Leveling(spline, *read);
    const auto                    degree = static_cast<std::size_t>(spline.Degree);
    const std::size_t             count  = spline.Poles.size();
    const std::span<const double> knots  = spline.Knots;

    const Basis basis = Spanned(spline);
    const auto  poles = [&spline, power](std::size_t index) -> Homogeneous
    {
      return Pole(spline, index, power);
    };

    /**
     *  The derivative's control points are its slopes, each the difference of two poles over the
     *  span of knots they share. De Boor reads only those whose span holds the parameter's, and
     *  that span is never empty, so no division here is by zero.
     */
    const Basis slopes = {
        .Degree = degree - 1, .Knots = knots.subspan(1, knots.size() - 2), .Count = count - 1};
    const auto steps = [&spline, &knots, degree, power](std::size_t index) -> Homogeneous
    {
      const double      span  = knots[index + degree + 1] - knots[index + 1];
      const Homogeneous start = Pole(spline, index, power);
      const Homogeneous end   = Pole(spline, index + 1, power);
      const double      scale = static_cast<double>(degree) / span;
      return {.X = (end.X - start.X) * scale,
              .Y = (end.Y - start.Y) * scale,
              .W = (end.W - start.W) * scale};
    };

    const Homogeneous point = Weighted(spline, Boor(basis, poles, along));
    /**
     *  Where the weight sums to nothing, past the domain or beside weights too far apart for a
     *  double to hold their sum, the spline is at no point and has no tangent.
     */
    if (Zero(point.W))
    {
      return {};
    }
    const Homogeneous slope = Boor(slopes, steps, along);
    const Point       place = Projected(point);
    return {.X = (slope.X - (slope.W * place.X)) / point.W,
            .Y = (slope.Y - (slope.W * place.Y)) / point.W};
  }

  auto Derivative(const Bezier::Quadratic &quadratic, double along) -> Vector
  {
    const Vector start = quadratic.Control - quadratic.Start;
    const Vector end   = quadratic.End - quadratic.Control;
    return ((start * (1.0 - along)) + (end * along)) * 2.0;
  }

  auto Derivative(const Bezier::Cubic &cubic, double along) -> Vector
  {
    const double rest   = 1.0 - along;
    const Vector first  = cubic.Out - cubic.Start;
    const Vector second = cubic.In - cubic.Out;
    const Vector third  = cubic.End - cubic.In;
    const Vector start  = (first * rest) + (second * along);
    const Vector end    = (second * rest) + (third * along);
    return ((start * rest) + (end * along)) * 3.0;
  }

  auto Derivative(const Curve &curve, double along) -> Vector
  {
    return std::visit([along](const auto &shape) -> Vector { return Derivative(shape, along); },
                      curve);
  }
} // namespace Geoxy
