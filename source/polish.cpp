/***************************************************************************************************
 *  @file      source/polish.cpp
 *  @brief     Where two curves cross, polished in double-double: each curve worked out to about 106
 *             bits, Newton's method run on the gap between them, the point rounded once.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "basis.hpp"
#include "crossing.hpp"
#include "eased.hpp"
#include "twofold.hpp"
#include "walk.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
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
#include <utility>
#include <variant>

namespace Geoxy::Crossings
{
  namespace
  {
    using Twofold::Add;
    using Twofold::Difference;
    using Twofold::Divide;
    using Twofold::Double;
    using Twofold::Lifted;
    using Twofold::Multiply;
    using Twofold::Two;
    using Twofold::TwoProduct;
    using Twofold::TwoSum;

    /**
     *  A point in double-double.
     */
    struct Spot
    {
      Double X = {}; /**< Along x. */
      Double Y = {}; /**< Along y. */
    };

    /**
     *  A homogeneous point in double-double: a rational spline's control point times its weight.
     */
    struct Weighted
    {
      Double X = {}; /**< Along x, times the weight. */
      Double Y = {}; /**< Along y, times the weight. */
      Double W = {}; /**< The weight.                */
    };

    /**
     *  A difference of two doubles, exactly.
     *
     *  @param  first  The double.
     *  @param  second The double taken from it.
     *  @return        The difference.
     */
    auto Subtracted(double first, double second) -> Double
    {
      const Two difference = TwoSum(first, -second);
      return Add(Lifted(difference.High), Lifted(difference.Low));
    }

    /**
     *  A point of the segment between two points.
     *
     *  @param  start Where it starts.
     *  @param  end   Where it ends.
     *  @param  along How far along it.
     *  @return       The point.
     */
    auto Lerped(const Spot &start, const Spot &end, Double along) -> Spot
    {
      return {.X = Add(start.X, Multiply(Difference(end.X, start.X), along)),
              .Y = Add(start.Y, Multiply(Difference(end.Y, start.Y), along))};
    }

    /**
     *  A point, lifted.
     *
     *  @param  point The point.
     *  @return       It in double-double.
     */
    auto Spotted(Point point) -> Spot
    {
      return {.X = Lifted(point.X), .Y = Lifted(point.Y)};
    }

    /**
     *  A point of a line through a point along a way.
     *
     *  @param  through The point at 0.
     *  @param  way     The way, per unit.
     *  @param  along   How far.
     *  @return         The point.
     */
    auto Walked(Point through, Vector way, Double along) -> Spot
    {
      return {.X = Add(Lifted(through.X), Multiply(Lifted(way.X), along)),
              .Y = Add(Lifted(through.Y), Multiply(Lifted(way.Y), along))};
    }

    /**
     *  A polyline's point, where its segment there is straight.
     *
     *  @param  polyline The polyline.
     *  @param  along    The place.
     *  @return          The point; nothing on a bulge or with no segment.
     */
    auto Evaluated(const Polyline &polyline, Double along) -> std::optional<Spot>
    {
      if (Count(polyline) == 0)
      {
        return std::nullopt;
      }
      const Step step = Locate(polyline, along.High);
      if (!Straightness(step))
      {
        return std::nullopt;
      }
      const double index = along.High - step.Fraction;
      return Lerped(Spotted(step.Start), Spotted(step.End),
                    Add(Subtracted(along.High, index), Lifted(along.Low)));
    }

    /**
     *  A spline's control point, homogeneous, its products exact.
     *
     *  @param  spline The spline.
     *  @param  index  Which control point.
     *  @param  power  A power of two to weigh by as well, as Leveling gives it.
     *  @return        The point times its weight, and the weight.
     */
    auto Lifted(const Spline &spline, std::size_t index, int power) -> Weighted
    {
      const Point  pole   = spline.Poles[index];
      const double weight = spline.Weights.empty() ? 1.0 : std::ldexp(spline.Weights[index], power);
      const Two    x      = TwoProduct(pole.X, weight);
      const Two    y      = TwoProduct(pole.Y, weight);
      return {.X = {.High = x.High, .Low = x.Low},
              .Y = {.High = y.High, .Low = y.Low},
              .W = Lifted(weight)};
    }

    /**
     *  A spline's point at a place told from both ends of its span, by de Boor's algorithm in
     *  double-double.
     *
     *  Each mix is the one pole times how far the place is from the other's knot and the other
     *  times how far it is from the one's, each distance worked out on its own from the end of the
     *  span it is told from: beside a heavy weight the place is a sliver from a knot, and one
     *  distance taken as what the other leaves of 1 would lose the sliver.
     *
     *  @param  spline The spline, valid.
     *  @param  site   The place.
     *  @return        The point.
     */
    auto Evaluated(const Spline &spline, const Sited &site) -> Spot
    {
      const Basis       basis  = Spanned(spline);
      const std::size_t degree = basis.Degree;
      const std::size_t span   = site.Span;
      const int         power  = Leveling(spline, span);

      std::array<Weighted, MaxDegree + 1> work{};
      const std::span<Weighted>           cells{work};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        cells[index] = Lifted(spline, span - degree + index, power);
      }
      for (std::size_t round = 1; round <= degree; ++round)
      {
        for (std::size_t index = degree; index >= round; --index)
        {
          const double start = basis.Knots[span - degree + index];
          const double end   = basis.Knots[span + 1 + index - round];
          const Double width = Subtracted(end, start);
          const Double fraction =
              Divide(Add(site.Ahead, Subtracted(basis.Knots[span], start)), width);
          const Double rest =
              Divide(Add(site.Short, Subtracted(end, basis.Knots[span + 1])), width);
          const Weighted lhs = cells[index - 1];
          const Weighted rhs = cells[index];
          cells[index]       = {.X = Add(Multiply(lhs.X, rest), Multiply(rhs.X, fraction)),
                                .Y = Add(Multiply(lhs.Y, rest), Multiply(rhs.Y, fraction)),
                                .W = Add(Multiply(lhs.W, rest), Multiply(rhs.W, fraction))};
        }
      }
      const Weighted point = cells[degree];
      return {.X = Divide(point.X, point.W), .Y = Divide(point.Y, point.W)};
    }

    /**
     *  A spline's point at a place of its parameter, in double-double.
     *
     *  @param  spline The spline, valid.
     *  @param  along  The place.
     *  @return        The point.
     */
    auto Evaluated(const Spline &spline, Double along) -> Spot
    {
      /**
       *  A place told in two doubles may lie just short of a knot its leading double is on: it is
       *  read on the span it lies in, not on the one that starts there.
       */
      const Basis  basis     = Spanned(spline);
      const double below     = std::nextafter(along.High, -std::numeric_limits<double>::infinity());
      const std::size_t span = Span(basis, along.Low < 0.0 ? below : along.High);
      return Evaluated(spline, Sited{.Span  = span,
                                     .Ahead = Difference(along, Lifted(basis.Knots[span])),
                                     .Short = Difference(Lifted(basis.Knots[span + 1]), along)});
    }

    /**
     *  A curve's point in double-double, where it is a polynomial or a ratio of two: a line, a ray,
     *  an endless line, a polyline's straight segment, a Bézier curve, a spline.
     *
     *  @param  curve The curve.
     *  @param  along The place.
     *  @return       The point; nothing on an arc, an ellipse or a bulge.
     */
    auto Evaluated(const Curve &curve, Double along) -> std::optional<Spot>
    {
      if (const auto *line = std::get_if<Line>(&curve))
      {
        return Lerped(Spotted(line->Start), Spotted(line->End), along);
      }
      if (const auto *ray = std::get_if<Ray>(&curve))
      {
        return Walked(ray->Start, ray->Way, along);
      }
      if (const auto *xline = std::get_if<Xline>(&curve))
      {
        return Walked(xline->Through, xline->Way, along);
      }
      if (const auto *polyline = std::get_if<Polyline>(&curve))
      {
        return Evaluated(*polyline, along);
      }
      if (const auto *quadratic = std::get_if<Bezier::Quadratic>(&curve))
      {
        const Spot control = Spotted(quadratic->Control);
        return Lerped(Lerped(Spotted(quadratic->Start), control, along),
                      Lerped(control, Spotted(quadratic->End), along), along);
      }
      if (const auto *cubic = std::get_if<Bezier::Cubic>(&curve))
      {
        const Spot out   = Spotted(cubic->Out);
        const Spot inner = Spotted(cubic->In);
        const Spot one   = Lerped(Spotted(cubic->Start), out, along);
        const Spot two   = Lerped(out, inner, along);
        return Lerped(Lerped(one, two, along),
                      Lerped(two, Lerped(inner, Spotted(cubic->End), along), along), along);
      }
      if (const auto *spline = std::get_if<Spline>(&curve))
      {
        return Evaluated(*spline, along);
      }
      return std::nullopt;
    }

    /**
     *  A curve's point at a place, as the curve is read: at the place itself, or at the place of
     *  the curve that a stand-in's stands for.
     *
     *  @param  read  The curve, as it is read.
     *  @param  along The place given.
     *  @return       The point; nothing for a bulge.
     */
    auto Evaluated(const Read &read, Double along) -> std::optional<Spot>
    {
      if (read.Way != nullptr)
      {
        return Evaluated(std::get<Spline>(*read.Shape), Back(*read.Way, along));
      }
      return Evaluated(*read.Shape, along);
    }

    /**
     *  Whether a place lies in a curve's domain, to within a few steps of a double.
     *
     *  @param  curve The curve.
     *  @param  along The place.
     *  @return       Whether it does.
     */
    auto Interior(const Curve &curve, double along) -> bool
    {
      const Range  range = Domain(curve);
      const double slack = 16.0 * std::numeric_limits<double>::epsilon();
      return along >= range.Start - (slack * std::max(1.0, std::abs(range.Start)))
          && along <= range.End + (slack * std::max(1.0, std::abs(range.End)));
    }

    /**
     *  A conic as the points where a quadratic is zero: |A (p - o)|^2 - level, A a matrix of the
     *  curve's own numbers, the level worked out in double-double. A circle has A the identity and
     *  its radius squared for level; an ellipse the adjugate of its map and its determinant
     *  squared. A polyline's bulge from P0 to P1 of bulge b is its circle as 4 b (p - P0).(p - P1)
     *  - 2 (1 - b^2) (p - P0).n, n the chord turned a quarter to the left: its exact rational
     *    circle, whose center no double holds.
     */
    struct Ovoid
    {
      Point                 Origin = {};    /**< The center, o; a bulge's first vertex. */
      std::array<Vector, 2> Rows   = {};    /**< The matrix A, by rows.                 */
      Double                Level  = {};    /**< What |A (p - o)|^2 is on it.           */
      double                Size   = 0.0;   /**< The root of its level, or a bulge's
                                                 chord, against which a step is lost.     */
      double Span = 0.0;                    /**< How long it is: its radius, its longer
                                                 axis, or a bulge's chord.                */
      bool                  Bowed  = false; /**< Whether a bulge's circle.    */
      Point                 End    = {};    /**< A bulge's second vertex, P1. */
      double                Bend   = 0.0;   /**< A bulge's bulge, b.          */
      std::array<Double, 2> Across = {};    /**< A bulge's chord turned a quarter to the
                                                 left, n, exactly.                        */
      Double Twice = {};                    /**< 2 (1 - b^2), as Valued works it out. */
    };

    /**
     *  A curve as a conic at a place, where it lies on one there: an arc, an ellipse not flat, a
     *  polyline's bulge.
     *
     *  @param  curve The curve.
     *  @param  along The place, to choose a polyline's segment.
     *  @return       The conic; nothing for other curves, or a polyline's straight segment.
     */
    auto Implicit(const Curve &curve, double along) -> std::optional<Ovoid>
    {
      if (const auto *arc = std::get_if<Arc>(&curve))
      {
        const Two squared = TwoProduct(arc->Radius, arc->Radius);
        return Ovoid{.Origin = arc->Center,
                     .Rows   = {Vector{.X = 1.0, .Y = 0.0}, Vector{.X = 0.0, .Y = 1.0}},
                     .Level  = {.High = squared.High, .Low = squared.Low},
                     .Size   = std::abs(arc->Radius),
                     .Span   = std::abs(arc->Radius)};
      }
      if (const auto *ellipse = std::get_if<Ellipse>(&curve))
      {
        const Vector x   = ellipse->Map.X;
        const Vector y   = ellipse->Map.Y;
        const Two    one = TwoProduct(x.X, y.Y);
        const Two    two = TwoProduct(x.Y, y.X);
        const Double det =
            Difference({.High = one.High, .Low = one.Low}, {.High = two.High, .Low = two.Low});
        if (std::fpclassify(det.High) == FP_ZERO)
        {
          return std::nullopt;
        }
        const Double level = Multiply(det, det);
        return Ovoid{.Origin = ellipse->Map.Origin,
                     .Rows   = {Vector{.X = y.Y, .Y = -y.X}, Vector{.X = -x.Y, .Y = x.X}},
                     .Level  = level,
                     .Size   = std::sqrt(level.High),
                     .Span   = std::max(Length(x), Length(y))};
      }
      const auto *polyline = std::get_if<Polyline>(&curve);
      if (polyline == nullptr || Count(*polyline) == 0)
      {
        return std::nullopt;
      }
      const Step step = Locate(*polyline, along);
      if (Straightness(step))
      {
        return std::nullopt;
      }
      /**
       *  What Valued takes of the bulge at every step of the polish, worked out once.
       */
      const Two square = TwoProduct(step.Bulge, step.Bulge);
      return Ovoid{
          .Origin = step.Start,
          .Rows   = {},
          .Level  = {},
          .Size   = Length(step.End - step.Start),
          .Span   = Length(step.End - step.Start),
          .Bowed  = true,
          .End    = step.End,
          .Bend   = step.Bulge,
          .Across = {Subtracted(step.Start.Y, step.End.Y), Subtracted(step.End.X, step.Start.X)},
          .Twice  = Multiply(Lifted(2.0),
                             Difference(Lifted(1.0), {.High = square.High, .Low = square.Low}))};
    }

    /**
     *  A conic's quadratic at a point, in double-double: zero on the curve.
     *
     *  @param  conic The conic.
     *  @param  point The point.
     *  @return       The value.
     */
    auto Valued(const Ovoid &conic, const Spot &point) -> Double
    {
      const Double x = Difference(point.X, Lifted(conic.Origin.X));
      const Double y = Difference(point.Y, Lifted(conic.Origin.Y));
      if (conic.Bowed)
      {
        const Double ahead = Difference(point.X, Lifted(conic.End.X));
        const Double above = Difference(point.Y, Lifted(conic.End.Y));
        const Double inner = Add(Multiply(x, ahead), Multiply(y, above));
        const Double side  = Add(Multiply(x, conic.Across[0]), Multiply(y, conic.Across[1]));
        return Difference(Multiply(Lifted(4.0 * conic.Bend), inner), Multiply(conic.Twice, side));
      }
      const Vector first = conic.Rows[0];
      const Vector other = conic.Rows[1];
      const Double u     = Add(Multiply(Lifted(first.X), x), Multiply(Lifted(first.Y), y));
      const Double v     = Add(Multiply(Lifted(other.X), x), Multiply(Lifted(other.Y), y));
      return Difference(Add(Multiply(u, u), Multiply(v, v)), conic.Level);
    }

    /**
     *  The gradient of a conic's quadratic at a point: 2 A^T A (p - o); a bulge's, 4 b (2 p - P0 -
     *  P1) - 2 (1 - b^2) n.
     *
     *  @param  conic The conic.
     *  @param  point The point.
     *  @return       The gradient.
     */
    auto Gradient(const Ovoid &conic, Point point) -> Vector
    {
      const Vector apart = point - conic.Origin;
      if (conic.Bowed)
      {
        const Vector ahead  = point - conic.End;
        const Vector across = Perpendicular(conic.End - conic.Origin);
        const double rest   = 1.0 - (conic.Bend * conic.Bend);
        return ((apart + ahead) * (4.0 * conic.Bend)) - (across * (2.0 * rest));
      }
      const Vector first = conic.Rows[0];
      const Vector other = conic.Rows[1];
      const double u     = Dot(first, apart);
      const double v     = Dot(other, apart);
      return {.X = 2.0 * ((first.X * u) + (other.X * v)),
              .Y = 2.0 * ((first.Y * u) + (other.Y * v))};
    }

    /**
     *  A point rounded to doubles.
     *
     *  @param  spot The point in double-double.
     *  @return      Each coordinate's leading part, the double nearest it.
     */
    auto Pinned(const Spot &spot) -> Point
    {
      return {.X = spot.X.High, .Y = spot.Y.High};
    }

    /**
     *  A point moved a last, tiny step along a curve's tangent: within the curve's bend times the
     *  step squared, some 2^-128 of its size, of the curve's point there.
     *
     *  @param  spot  The point.
     *  @param  way   The tangent.
     *  @param  along The step.
     *  @return       The point moved.
     */
    auto Ahead(const Spot &spot, Vector way, double along) -> Spot
    {
      return {.X = Add(spot.X, Multiply(Lifted(way.X), Lifted(along))),
              .Y = Add(spot.Y, Multiply(Lifted(way.Y), Lifted(along)))};
    }

    /**
     *  Whether a step moves a point by under 2^-64 of the size about it: a 2^-11 part of a step of
     *  a double there, lost in double-double's width.
     *
     *  @param  reach How far the step moves the point.
     *  @param  size  The size: the point's largest coordinate, or a curve's speed.
     *  @return       Whether it is lost so.
     */
    auto Lost(double reach, double size) -> bool
    {
      return reach <= 0x1p-64 * size;
    }

    /**
     *  Whether a place moved by no more than 2^-20 of its size: polished, not jumped to another
     *  crossing.
     *
     *  @param  from Where it was.
     *  @param  now  Where it is.
     *  @return      Whether it stayed.
     */
    auto Stayed(double from, double now) -> bool
    {
      return std::abs(now - from) <= 0x1p-20 * std::max(1.0, std::abs(from));
    }

    /**
     *  Whether two directions are further than 2^-48 radians from parallel: nearer than that, the
     *  doubles that hold them tell no angle between them at all.
     *
     *  @param  one   One.
     *  @param  other The other.
     *  @return       Whether they are.
     */
    auto Across(Vector one, Vector other) -> bool
    {
      return std::abs(Cross(one, other)) > 0x1p-48 * Length(one) * Length(other);
    }

    /**
     *  Two curves each a polynomial or a ratio of two: Newton's method on both places.
     *
     *  The crossing is the hit's where the places stayed, or where the hit lay on both curves and
     *  every step was no more than half the one before: then it is the one crossing within twice
     *  the first step, however far that is, as it is far along two curves that cross at a slight
     *  angle from a hit a step of a double off either.
     *
     *  @param  curves The curves.
     *  @param  places Where the hit is on each.
     *  @return        The crossing.
     */
    auto Parametric(const std::array<Read, 2> &curves, const std::array<double, 2> &places)
        -> std::optional<Crossed>
    {
      std::array<Double, 2> held   = {Lifted(places[0]), Lifted(places[1])};
      double                last   = std::numeric_limits<double>::infinity();
      bool                  halved = true;
      for (int round = 0; round < 8; ++round)
      {
        const std::optional<Spot> here  = Evaluated(curves[0], held[0]);
        const std::optional<Spot> there = Evaluated(curves[1], held[1]);
        if (!here || !there)
        {
          return std::nullopt;
        }
        const Vector apart = {.X = Difference(here->X, there->X).High,
                              .Y = Difference(here->Y, there->Y).High};
        const Vector way   = Derivative(*curves[0].Stand, held[0].High);
        const Vector back  = Derivative(*curves[1].Stand, held[1].High);
        if (!Across(back, way))
        {
          return std::nullopt;
        }
        const double turn   = Cross(back, way);
        const double first  = Cross(apart, back) / turn;
        const double second = Cross(apart, way) / turn;
        held                = {Add(held[0], Lifted(first)), Add(held[1], Lifted(second))};
        const double size =
            std::max({std::abs(here->X.High), std::abs(here->Y.High), Length(way), Length(back)});
        const double reach =
            std::max(std::abs(first) * Length(way), std::abs(second) * Length(back));
        halved = halved && reach <= 0.5 * last && (round > 0 || Length(apart) <= 0x1p-20 * size);
        last   = reach;
        if (Lost(reach, size))
        {
          if (!halved && (!Stayed(places[0], held[0].High) || !Stayed(places[1], held[1].High)))
          {
            return std::nullopt;
          }
          return Crossed{.At    = Pinned(Ahead(*here, way, first)),
                         .Along = {held[0].High, held[1].High}};
        }
      }
      return std::nullopt;
    }

    /**
     *  A polynomial curve and a conic: Newton's method on the first's place alone, the conic's
     *  quadratic at its point driven to zero; the crossing is the hit's as Parametric tells.
     *
     *  @param  curve The curve.
     *  @param  conic The conic.
     *  @param  along Where the hit is on the curve.
     *  @return       The crossing's place on the curve and its point.
     */
    auto Mixed(const Read &curve, const Ovoid &conic, double along)
        -> std::optional<std::pair<double, Point>>
    {
      Double held   = Lifted(along);
      double last   = std::numeric_limits<double>::infinity();
      bool   halved = true;
      for (int round = 0; round < 8; ++round)
      {
        const std::optional<Spot> here = Evaluated(curve, held);
        if (!here)
        {
          return std::nullopt;
        }
        const Vector way   = Derivative(*curve.Stand, held.High);
        const Vector slope = Gradient(conic, Pinned(*here));
        if (!(std::abs(Dot(slope, way)) > 0x1p-48 * Length(slope) * Length(way)))
        {
          return std::nullopt;
        }
        const double value = Valued(conic, *here).High;
        const double step  = -value / Dot(slope, way);
        held               = Add(held, Lifted(step));
        const double size = std::max({std::abs(here->X.High), std::abs(here->Y.High), Length(way)});
        const double reach = std::abs(step) * Length(way);
        const bool   near  = round > 0 || std::abs(value) <= 0x1p-20 * size * Length(slope);
        halved             = halved && near && reach <= 0.5 * last;
        last               = reach;
        if (Lost(reach, size))
        {
          if (!halved && !Stayed(along, held.High))
          {
            return std::nullopt;
          }
          return std::pair{held.High, Pinned(Ahead(*here, way, step))};
        }
      }
      return std::nullopt;
    }

    /**
     *  Two conics: Newton's method on the point, both quadratics driven to zero.
     *
     *  @param  conics The conics.
     *  @param  near   Where the hit is.
     *  @return        The point.
     */
    auto Conics(const std::array<Ovoid, 2> &conics, Point near) -> std::optional<Point>
    {
      Spot held = Spotted(near);
      for (int round = 0; round < 8; ++round)
      {
        const Point  point = Pinned(held);
        const Vector one   = Gradient(conics[0], point);
        const Vector other = Gradient(conics[1], point);
        if (!Across(one, other))
        {
          return std::nullopt;
        }
        const double first  = Valued(conics[0], held).High;
        const double second = Valued(conics[1], held).High;
        const double det    = Cross(one, other);
        const Vector step   = {.X = ((second * one.Y) - (first * other.Y)) / det,
                               .Y = ((first * other.X) - (second * one.X)) / det};
        held                = {.X = Add(held.X, Lifted(step.X)), .Y = Add(held.Y, Lifted(step.Y))};
        if (Lost(Length(step),
                 std::max({std::abs(point.X), std::abs(point.Y), conics[0].Size, conics[1].Size})))
        {
          /**
           *  Polished, not jumped to another crossing: moved by no more than 2^-20 of the size
           *  about it, the point's or the conics', so that the same crossing scaled by a power of
           *  two is polished the same.
           */
          const Point  out = Pinned(held);
          const double about =
              std::max({std::abs(near.X), std::abs(near.Y), conics[0].Span, conics[1].Span});
          if (!(std::abs(out.X - near.X) <= 0x1p-20 * about)
              || !(std::abs(out.Y - near.Y) <= 0x1p-20 * about))
          {
            return std::nullopt;
          }
          return out;
        }
      }
      return std::nullopt;
    }

    /**
     *  Whether Newton's method may leave a crossing's coordinate short of its last bit: its steps
     *  stop under 2^-64 of the size about the point, the point's or the conics', which settles a
     *  coordinate to the bit only where it is not under 2^-8 of that size, as at a crossing near an
     *  axis far out, or near the origin on conics far larger.
     *
     *  @param  conics The conics.
     *  @param  near   Where the crossing is, roughly.
     *  @return        Whether it may.
     */
    auto Lopsided(const std::array<Ovoid, 2> &conics, Point near) -> bool
    {
      const double size =
          std::max({std::abs(near.X), std::abs(near.Y), conics[0].Span, conics[1].Span});
      return std::min(std::abs(near.X), std::abs(near.Y)) < 0x1p-8 * size;
    }

    /**
     *  Where two conics cross near a hit: by Newton's method on the point; but two bulges' circles
     *  crossed exactly, the double nearest the crossing nearer the hit, first where Newton's method
     *  may leave a coordinate short of its last bit, and last where it finds no crossing, as from a
     *  hit found so far off that it would jump past what it may move.
     *
     *  @param  curves The two curves, each a conic at its place.
     *  @param  conics Their conics.
     *  @param  hit    Where the hit is, and its place on each.
     *  @return        The point; nothing where neither finds it.
     */
    auto Conjoined(const std::array<Read, 2> &curves, const std::array<Ovoid, 2> &conics,
                   const Crossed &hit) -> std::optional<Point>
    {
      const bool bowed   = conics[0].Bowed && conics[1].Bowed;
      const auto exactly = [&curves, &hit] -> std::optional<Point>
      {
        return Bowed({&std::get<Polyline>(*curves[0].Shape), &std::get<Polyline>(*curves[1].Shape)},
                     hit.Along, hit.At);
      };
      const bool lopsided = bowed && Lopsided(conics, hit.At);
      if (const std::optional<Point> point = lopsided ? exactly() : std::nullopt)
      {
        return point;
      }
      if (const std::optional<Point> point = Conics(conics, hit.At))
      {
        return point;
      }
      return bowed && !lopsided ? exactly() : std::nullopt;
    }
  } // namespace

  auto Retold(const std::array<Read, 2> &curves, const std::array<double, 2> &places, Point near)
      -> std::optional<Crossed>
  {
    const std::optional<Ovoid> one   = Implicit(*curves[0].Shape, places[0]);
    const std::optional<Ovoid> other = Implicit(*curves[1].Shape, places[1]);
    std::optional<Crossed>     out;
    if (one && other)
    {
      if (const std::optional<Point> point =
              Conjoined(curves, {*one, *other}, {.At = near, .Along = places, .On = {}}))
      {
        out = Crossed{.At = *point, .Along = places};
      }
    }
    else if (one)
    {
      if (const auto found = Mixed(curves[1], *one, places[1]))
      {
        out = Crossed{.At = found->second, .Along = {places[0], found->first}};
      }
    }
    else if (other)
    {
      if (const auto found = Mixed(curves[0], *other, places[0]))
      {
        out = Crossed{.At = found->second, .Along = {found->first, places[1]}};
      }
    }
    else
    {
      out = Parametric(curves, places);
    }
    if (!out)
    {
      return std::nullopt;
    }
    /**
     *  A crossing is placed on each conic by its point, as Crossing places one, so that its place
     *  there is the point's whichever candidate it was polished from: on an arc or an ellipse by
     *  the point's angle, on a polyline's arc by its angle from the chord. A conic the point is not
     *  placed on, as just past its end, keeps the place it was polished from; a polynomial keeps
     *  the place Newton's method settled on, which no start changes.
     */
    /**
     *  Each curve's point at its place, as At works it out: on a polyline, read as itself as only a
     *  spline has a stand-in, the one Located worked out by At already, at a place on the segment
     *  it names and so in the domain, which clamping leaves as it is.
     */
    const std::array<bool, 2>           conic  = {one.has_value(), other.has_value()};
    std::array<std::optional<Point>, 2> points = {};
    for (std::size_t side = 0; side < conic.size(); ++side)
    {
      const Curve &shape = *curves.at(side).Shape;
      const double along = out->Along.at(side);
      if (const std::optional<Place> place =
              conic.at(side) ? Located(shape, along, out->At) : std::nullopt)
      {
        out->Along.at(side) = place->Along;
        points.at(side) = std::holds_alternative<Polyline>(shape) ? std::optional<Point>{place->At}
                                                                  : std::nullopt;
      }
    }
    if (!Interior(*curves[0].Stand, out->Along[0]) || !Interior(*curves[1].Stand, out->Along[1]))
    {
      return std::nullopt;
    }
    for (std::size_t side = 0; side < conic.size(); ++side)
    {
      const Curve &stand               = *curves.at(side).Stand;
      const Range  domain              = Domain(stand);
      out->Along.at(side)              = std::clamp(out->Along.at(side), domain.Start, domain.End);
      const std::optional<Point> &kept = points.at(side);
      out->On.at(side)                 = kept ? *kept : At(stand, out->Along.at(side));
    }
    return out;
  }

  auto Polished(const std::array<const Curve *, 2> &curves, const std::array<double, 2> &places,
                Point near) -> std::optional<Crossed>
  {
    return Retold(std::array<Read, 2>{Read{.Shape = curves[0], .Stand = curves[0], .Way = nullptr},
                                      Read{.Shape = curves[1], .Stand = curves[1], .Way = nullptr}},
                  places, near);
  }
} // namespace Geoxy::Crossings
