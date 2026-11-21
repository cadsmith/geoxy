/***************************************************************************************************
 *  @file      source/curvature.cpp
 *  @brief     The unit tangent, the normal and the signed curvature of each curve at a place: in
 *             closed form for lines, arcs, ellipses and bulges, from the first two derivatives for
 *             Bézier curves and splines.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/curvature.hpp"

#include "basis.hpp"
#include "kernel/piece.hpp"
#include "real.hpp"
#include "walk.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  A curve's first two derivatives by its parameter at a place.
     */
    struct Moving
    {
      Vector First  = {}; /**< The first derivative.  */
      Vector Second = {}; /**< The second derivative. */
    };

    /**
     *  The signed curvature from a curve's first two derivatives, x' × x'' / |x'|^3: each divided
     *  by |x'| before they are crossed, and the cross once more, so that no cube overflows.
     *
     *  @param  moving The derivatives.
     *  @return        The curvature; none where the first is nothing.
     */
    auto Curving(const Moving &moving) -> std::optional<double>
    {
      const double speed = Length(moving.First);
      if (!(speed > 0.0))
      {
        return std::nullopt;
      }
      return Cross(moving.First / speed, moving.Second / speed) / speed;
    }

    /**
     *  A displacement turned a quarter turn counterclockwise, to its left.
     *
     *  @param  vector The displacement.
     *  @return        It turned.
     */
    auto Left(Vector vector) -> Vector
    {
      return Perpendicular(vector);
    }

    /**
     *  The curvature of a straight curve, where it has a tangent.
     *
     *  @return 0.
     */
    auto Flat(Vector /*tangent*/) -> double
    {
      return 0.0;
    }

    /**
     *  A homogeneous point's difference from another, times a factor.
     *
     *  @param  lhs    The point.
     *  @param  rhs    The point taken from it.
     *  @param  factor The factor.
     *  @return        (lhs - rhs) times the factor, coordinate by coordinate and the weight.
     */
    auto Differ(const Homogeneous &lhs, const Homogeneous &rhs, double factor) -> Homogeneous
    {
      return {.X = (lhs.X - rhs.X) * factor,
              .Y = (lhs.Y - rhs.Y) * factor,
              .W = (lhs.W - rhs.W) * factor};
    }

    /**
     *  The first two derivatives of a rational Bézier piece by its parameter, by de Casteljau's
     *  algorithm run to its last three points, then the quotient rule.
     *
     *  Of degree n, its homogeneous point is the last mix, its first derivative n times the
     *  difference of the last two points, its second n (n - 1) times the second difference of the
     *  last three. Then x = X / W, x' = (X' - W' x) / W and x'' = (X'' - 2 W' x' - W'' x) / W.
     *
     *  @param  piece The piece.
     *  @param  along Where, 0 at its start and 1 at its end.
     *  @return       The derivatives.
     */
    auto Moved(const Kernel::Piece &piece, double along) -> Moving
    {
      std::array<Homogeneous, MaxDegree + 1> work   = piece.Poles;
      const std::span<Homogeneous>           cells  = work;
      const std::size_t                      degree = piece.Degree;
      const auto                             count  = static_cast<double>(degree);
      for (std::size_t round = 1; round + 2 <= degree; ++round)
      {
        for (std::size_t index = 0; index + round <= degree; ++index)
        {
          cells[index] = Mix(cells[index], cells[index + 1], along);
        }
      }
      Homogeneous point  = Mix(cells[0], cells[1], along);
      Homogeneous first  = Differ(cells[1], cells[0], 1.0);
      Homogeneous second = {};
      if (degree > 1)
      {
        const Homogeneous start = Mix(cells[0], cells[1], along);
        const Homogeneous end   = Mix(cells[1], cells[2], along);
        const Homogeneous ahead = Differ(cells[2], cells[1], 1.0);
        second = Differ(ahead, Differ(cells[1], cells[0], 1.0), count * (count - 1));
        first  = Differ(end, start, count);
        point  = Mix(start, end, along);
      }
      const double weight = point.W;
      const Vector place  = {.X = point.X / weight, .Y = point.Y / weight};
      const Vector speed  = {.X = (first.X - (first.W * place.X)) / weight,
                             .Y = (first.Y - (first.W * place.Y)) / weight};
      const Vector bend   = {
          .X = (second.X - (2.0 * first.W * speed.X) - (second.W * place.X)) / weight,
          .Y = (second.Y - (2.0 * first.W * speed.Y) - (second.W * place.Y)) / weight};
      return {.First = speed, .Second = bend};
    }

    /**
     *  The first two derivatives of a spline by its parameter, from its span's piece.
     *
     *  @param  spline The spline.
     *  @param  along  A knot value.
     *  @return        The derivatives; none where the spline cannot be read, as At tells it.
     */
    auto Moved(const Spline &spline, double along) -> std::optional<Moving>
    {
      const std::optional<std::size_t> span = Legible(spline, along);
      if (!span)
      {
        return std::nullopt;
      }
      const double start  = spline.Knots[*span];
      const double width  = spline.Knots[*span + 1] - start;
      const Moving moving = Moved(Kernel::Section(spline, *span), (along - start) / width);
      return Moving{.First = moving.First / width, .Second = moving.Second / (width * width)};
    }

    /**
     *  The first two derivatives of a quadratic Bézier curve.
     *
     *  @param  quadratic The curve.
     *  @param  along     Where.
     *  @return           The derivatives.
     */
    auto Moved(const Bezier::Quadratic &quadratic, double along) -> Moving
    {
      const Vector bend =
          (quadratic.Start - quadratic.Control) + (quadratic.End - quadratic.Control);
      return {.First = Derivative(quadratic, along), .Second = bend * 2.0};
    }

    /**
     *  The first two derivatives of a cubic Bézier curve.
     *
     *  @param  cubic The curve.
     *  @param  along Where.
     *  @return       The derivatives.
     */
    auto Moved(const Bezier::Cubic &cubic, double along) -> Moving
    {
      const Vector start = (cubic.Start - cubic.Out) + (cubic.In - cubic.Out);
      const Vector end   = (cubic.Out - cubic.In) + (cubic.End - cubic.In);
      return {.First  = Derivative(cubic, along),
              .Second = ((start * (1.0 - along)) + (end * along)) * 6.0};
    }

    /**
     *  The way an ellipse runs, unscaled by its sweep: its map's image of the unit circle's
     *  tangent.
     *
     *  @param  ellipse The ellipse.
     *  @param  along   The fraction of its sweep.
     *  @return         The way.
     */
    auto Way(const Ellipse &ellipse, double along) -> Vector
    {
      return ellipse.Map * Perpendicular(Direction(ellipse.Start + (ellipse.Sweep * along)));
    }
  } // namespace

  auto Tangent(const Line &line, double along) -> std::optional<Vector>
  {
    return Unit(Derivative(line, along));
  }

  auto Tangent(const Ray &ray, double along) -> std::optional<Vector>
  {
    return Unit(Derivative(ray, along));
  }

  auto Tangent(const Xline &xline, double along) -> std::optional<Vector>
  {
    return Unit(Derivative(xline, along));
  }

  auto Tangent(const Arc &arc, double along) -> std::optional<Vector>
  {
    const double sweep = arc.Sweep.InRadians();
    if (Zero(arc.Radius) || Zero(sweep))
    {
      return std::nullopt;
    }
    const Vector way = Perpendicular(Direction(arc.Start + (arc.Sweep * along)));
    return (arc.Radius > 0.0) == (sweep > 0.0) ? way : -way;
  }

  auto Tangent(const Ellipse &ellipse, double along) -> std::optional<Vector>
  {
    const double sweep = ellipse.Sweep.InRadians();
    if (Zero(sweep))
    {
      return std::nullopt;
    }
    const Vector way = Way(ellipse, along);
    return Unit(sweep > 0.0 ? way : -way);
  }

  auto Tangent(const Polyline &polyline, double along) -> std::optional<Vector>
  {
    return Unit(Derivative(polyline, along));
  }

  auto Tangent(const Spline &spline, double along) -> std::optional<Vector>
  {
    return Moved(spline, along)
        .and_then([](const Moving &moving) -> std::optional<Vector> { return Unit(moving.First); });
  }

  auto Tangent(const Bezier::Quadratic &quadratic, double along) -> std::optional<Vector>
  {
    return Unit(Derivative(quadratic, along));
  }

  auto Tangent(const Bezier::Cubic &cubic, double along) -> std::optional<Vector>
  {
    return Unit(Derivative(cubic, along));
  }

  auto Tangent(const Curve &curve, double along) -> std::optional<Vector>
  {
    return std::visit([along](const auto &shape) -> std::optional<Vector>
                      { return Tangent(shape, along); }, curve);
  }

  auto Normal(const Line &line, double along) -> std::optional<Vector>
  {
    return Tangent(line, along).transform(Left);
  }

  auto Normal(const Ray &ray, double along) -> std::optional<Vector>
  {
    return Tangent(ray, along).transform(Left);
  }

  auto Normal(const Xline &xline, double along) -> std::optional<Vector>
  {
    return Tangent(xline, along).transform(Left);
  }

  auto Normal(const Arc &arc, double along) -> std::optional<Vector>
  {
    return Tangent(arc, along).transform(Left);
  }

  auto Normal(const Ellipse &ellipse, double along) -> std::optional<Vector>
  {
    return Tangent(ellipse, along).transform(Left);
  }

  auto Normal(const Polyline &polyline, double along) -> std::optional<Vector>
  {
    return Tangent(polyline, along).transform(Left);
  }

  auto Normal(const Spline &spline, double along) -> std::optional<Vector>
  {
    return Tangent(spline, along).transform(Left);
  }

  auto Normal(const Bezier::Quadratic &quadratic, double along) -> std::optional<Vector>
  {
    return Tangent(quadratic, along).transform(Left);
  }

  auto Normal(const Bezier::Cubic &cubic, double along) -> std::optional<Vector>
  {
    return Tangent(cubic, along).transform(Left);
  }

  auto Normal(const Curve &curve, double along) -> std::optional<Vector>
  {
    return Tangent(curve, along).transform(Left);
  }

  auto Curvature(const Line &line, double along) -> std::optional<double>
  {
    return Tangent(line, along).transform(Flat);
  }

  auto Curvature(const Ray &ray, double along) -> std::optional<double>
  {
    return Tangent(ray, along).transform(Flat);
  }

  auto Curvature(const Xline &xline, double along) -> std::optional<double>
  {
    return Tangent(xline, along).transform(Flat);
  }

  auto Curvature(const Arc &arc, double along) -> std::optional<double>
  {
    if (!Tangent(arc, along))
    {
      return std::nullopt;
    }
    return std::copysign(1.0 / std::abs(arc.Radius), arc.Sweep.InRadians());
  }

  auto Curvature(const Ellipse &ellipse, double along) -> std::optional<double>
  {
    const double sweep = ellipse.Sweep.InRadians();
    const double speed = Length(Way(ellipse, along));
    if (Zero(sweep) || !(speed > 0.0))
    {
      return std::nullopt;
    }
    const double turn = Determinant(ellipse.Map) / speed / speed / speed;
    return sweep > 0.0 ? turn : -turn;
  }

  auto Curvature(const Polyline &polyline, double along) -> std::optional<double>
  {
    if (!Tangent(polyline, along))
    {
      return std::nullopt;
    }
    const Step step = Locate(polyline, along);
    if (Straightness(step))
    {
      return 0.0;
    }
    const double bulge = step.Bulge;
    return 4.0 * bulge / (Length(step.End - step.Start) * (1.0 + (bulge * bulge)));
  }

  auto Curvature(const Spline &spline, double along) -> std::optional<double>
  {
    return Moved(spline, along).and_then(Curving);
  }

  auto Curvature(const Bezier::Quadratic &quadratic, double along) -> std::optional<double>
  {
    return Curving(Moved(quadratic, along));
  }

  auto Curvature(const Bezier::Cubic &cubic, double along) -> std::optional<double>
  {
    return Curving(Moved(cubic, along));
  }

  auto Curvature(const Curve &curve, double along) -> std::optional<double>
  {
    return std::visit([along](const auto &shape) -> std::optional<double>
                      { return Curvature(shape, along); }, curve);
  }
} // namespace Geoxy
