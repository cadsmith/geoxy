/***************************************************************************************************
 *  @file      source/extents.cpp
 *  @brief     The box of each curve, from where it reaches its extremes.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/extents.hpp"

#include "basis.hpp"
#include "exact/two.hpp"
#include "kernel/piece.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "trig.hpp"
#include "walk.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
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
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <optional>
#include <span>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  The unit roundoff: a rounding is off by at most this times its result.
     */
    constexpr double Roundoff = std::numeric_limits<double>::epsilon() / 2.0;

    /**
     *  Room for a row of binomial coefficients up to twice the highest degree.
     */
    constexpr std::size_t Twice = 2 * static_cast<std::size_t>(MaxDegree);

    /**
     *  How many cells a point of a spline of a low degree takes room for: a degree below it is
     *  worked out in that many, so that what is cleared is no more than it uses.
     */
    constexpr int Scant = 8;

    /**
     *  A double at least a value and a pad above it: the sum rounded, then a step of a double up,
     *  read from its bits rather than asked of nextafter; nothing added where there is no pad.
     *
     *  @param  value The value.
     *  @param  pad   The pad, not negative.
     *  @return       The bound.
     */
    auto Raised(double value, double pad) -> double
    {
      const double sum = value + pad;
      return pad > 0.0 ? std::fmax(sum, sum + Above(std::abs(sum))) : value;
    }

    /**
     *  A double at most a value less a pad below it.
     *
     *  @param  value The value.
     *  @param  pad   The pad, not negative.
     *  @return       The bound.
     */
    auto Lowered(double value, double pad) -> double
    {
      const double difference = value - pad;
      return pad > 0.0 ? std::fmin(difference, difference - Above(std::abs(difference))) : value;
    }

    /**
     *  A box grown to hold a point, each coordinate a pad either way, rounded outward.
     *
     *  @param  box   The box.
     *  @param  point The point.
     *  @param  pad   The pad along each axis.
     *  @return       The box grown.
     */
    auto Engulfed(const Box &box, Point point, Vector pad) -> Box
    {
      return Grown(Grown(box, {.X = Lowered(point.X, pad.X), .Y = Lowered(point.Y, pad.Y)}),
                   {.X = Raised(point.X, pad.X), .Y = Raised(point.Y, pad.Y)});
    }

    /**
     *  A coordinate plus an offset, exactly, rounded away from the coordinate: the double at least
     *  as far out as the exact sum.
     *
     *  @param  base   The coordinate.
     *  @param  offset The offset.
     *  @return        The sum, rounded outward.
     */
    auto Beyond(double base, double offset) -> double
    {
      const Exact::Two sum = Exact::TwoSum(base, offset);
      if (sum.Small * offset > 0.0)
      {
        const double step = Above(std::abs(sum.Big));
        return offset > 0.0 ? sum.Big + step : sum.Big - step;
      }
      return sum.Big;
    }

    /**
     *  A box with the side an offset from a center reaches set to the center plus that offset,
     *  rounded outward: the side of an extreme the curve passes, which no point of it, nor any At
     *  gives, lies past.
     *
     *  @param  box    The box.
     *  @param  center The center.
     *  @param  offset The offset, along one axis.
     *  @return        The box.
     */
    auto Reaching(Box box, Point center, Vector offset) -> Box
    {
      if (offset.X > 0.0)
      {
        box.Max.X = Beyond(center.X, offset.X);
      }
      if (offset.X < 0.0)
      {
        box.Min.X = Beyond(center.X, offset.X);
      }
      if (offset.Y > 0.0)
      {
        box.Max.Y = Beyond(center.Y, offset.Y);
      }
      if (offset.Y < 0.0)
      {
        box.Min.Y = Beyond(center.Y, offset.Y);
      }
      return box;
    }

    /**
     *  Where a sweep from a start first reaches a direction, as a part of the sweep.
     *
     *  @param  start  Where the sweep starts.
     *  @param  sweep  The sweep, its sign its way round.
     *  @param  target The direction.
     *  @return        The part of the sweep, in radians, its sign the sweep's; none when the sweep
     *                 stops short of it.
     */
    auto Arrival(Angle start, Angle sweep, Angle target) -> std::optional<double>
    {
      if (!Within(target, start, sweep))
      {
        return std::nullopt;
      }
      double forward = Wrapped(target - start).InDegrees();
      if (sweep.InDegrees() >= 0.0 && forward < 0.0)
      {
        forward += 360.0;
      }
      else if (sweep.InDegrees() < 0.0 && forward > 0.0)
      {
        forward -= 360.0;
      }
      return Degrees(forward).InRadians();
    }

    /**
     *  The range along one axis that a start reaches running without end one way.
     *
     *  @param  start Where it starts along the axis.
     *  @param  way   The way's component along the axis.
     *  @return       The least and the greatest it reaches.
     */
    auto Ranged(double start, double way) -> std::array<double, 2>
    {
      if (way > 0.0)
      {
        return {start, Infinite};
      }
      if (way < 0.0)
      {
        return {-Infinite, start};
      }
      return {start, start};
    }

    /**
     *  The range along one axis that a point reaches running without end both ways.
     *
     *  @param  through Where it is along the axis.
     *  @param  way     The way's component along the axis.
     *  @return         The least and the greatest it reaches: all of the axis when the way leaves
     *                  it.
     */
    auto Endwise(double through, double way) -> std::array<double, 2>
    {
      if (!Zero(way))
      {
        return {-Infinite, Infinite};
      }
      return {through, through};
    }

    /**
     *  A point of a curve as At works it out, and a bound on how far each of its coordinates is off
     *  from the curve's point there.
     */
    struct Reckoning
    {
      Point  At    = {}; /**< The point.                      */
      Vector Error = {}; /**< How far each coordinate is off. */
    };

    /**
     *  How far a point At gives on a Bézier curve may stray from the curve, along each axis,
     *  whatever the parameter: de Casteljau's algorithm by Lerp over a degree d is off by at most
     *  three steps of a double of the largest control coordinate a level, the difference, its
     *  scaling and the sum each rounded once, and the mixes carry their inputs' errors as convex
     *  sums do; (4 d + 1) steps bound it, and 2^-40 of that covers the second-order terms. A
     *  product that falls among the doubles below the normal ones is off by half a step of the
     *  least more, d of them at the most: a step of the least a level (tests/proof/pads.py).
     *
     *  @param  controls The curve's control points.
     *  @return          The bound along each axis.
     */
    auto Wandering(std::initializer_list<Point> controls) -> Vector
    {
      Vector largest = {};
      for (const Point control : controls)
      {
        largest = {.X = std::max(largest.X, std::abs(control.X)),
                   .Y = std::max(largest.Y, std::abs(control.Y))};
      }
      const auto   degree = static_cast<double>(controls.size() - 1);
      const double under  = degree * 0x1p-1074;
      return (largest * (((4.0 * degree) + 1.0) * Roundoff * (1.0 + 0x1p-40)))
           + Vector{.X = under, .Y = under};
    }

    /**
     *  The box a Bézier curve reaches, grown by its point at a parameter where an axis turns, if
     *  that parameter lies strictly between 0 and 1. The point At gives there is within its bound
     *  of the curve, and so of the curve's true extreme, the parameter's own rounding moving the
     *  point along the curve to second order only; any other point At gives is within the bound of
     *  the curve, inside that extreme. Twice the bound beyond the point holds them all.
     *
     *  @param  box   The box so far.
     *  @param  point The point there, and how far At may stray.
     *  @param  along The parameter.
     *  @return       The box.
     */
    auto Turnpoint(const Box &box, const Reckoning &point, double along) -> Box
    {
      return along > 0.0 && along < 1.0 ? Engulfed(box, point.At, point.Error * 2.0) : box;
    }

    /**
     *  The roots of a t^2 + b t + c, by the stable quadratic formula, that lie strictly between 0
     *  and 1.
     *
     *  The roots are those of the coefficients times any power of two, which is exact: past 2^500
     *  or under 2^-500, where the discriminant's squares would overflow or come to nothing, a
     *  quadratic's are scaled so that the largest is near 1.
     *
     *  @param  a The square's coefficient.
     *  @param  b The linear coefficient.
     *  @param  c The constant.
     *  @return   Up to two roots; a slot that holds none holds -1.
     */
    auto Roots(double a, double b, double c) -> std::array<double, 2>
    {
      if (Zero(a))
      {
        return {Zero(b) ? -1.0 : -c / b, -1.0};
      }
      const double largest      = std::max({std::abs(a), std::abs(b), std::abs(c)});
      const bool   far          = largest > 0x1p500 || largest < 0x1p-500;
      const int    power        = far ? -std::ilogb(largest) : 0;
      const double square       = std::ldexp(a, power);
      const double linear       = std::ldexp(b, power);
      const double constant     = std::ldexp(c, power);
      const double discriminant = (linear * linear) - (4.0 * square * constant);
      if (discriminant < 0.0)
      {
        return {-1.0, -1.0};
      }
      const double q = -0.5 * (linear + std::copysign(std::sqrt(discriminant), linear));
      return {q / square, Zero(q) ? -1.0 : constant / q};
    }

    /**
     *  A bound on an angle's magnitude in radians, read without converting it: an angle held in
     *  degrees is taken at its degrees times 0.0175, a little more than pi / 180.
     *
     *  @param  angle The angle.
     *  @return       The bound.
     */
    auto Spanning(Angle angle) -> double
    {
      return angle.HeldInDegrees() ? std::abs(angle.InDegrees()) * 0.0175
                                   : std::abs(angle.InRadians());
    }

    /**
     *  How far an arc's end may lie from where At puts it, along one axis, and how far At's points
     *  near it may stray past it: its direction's coordinate off by a step of a double, the product
     *  and the sum by half a step each, both ways; and the end's angle off by half a step of the
     *  start's and the sweep's magnitudes, which moves it along the circle. Nothing where the
     *  direction's coordinate is -1, 0 or 1, as at a quarter turn: the coordinate is then the
     *  center's or an extreme's, and the arc stays on its side of it.
     *
     *  @param  center    The center's coordinate.
     *  @param  component The end direction's coordinate.
     *  @param  swing     What the radius adds to the pad: the radius times 8 and the start's and
     *                    the sweep's magnitudes in radians, times a step of a double.
     *  @return           The pad.
     */
    auto Wobbled(double center, double component, double swing) -> double
    {
      if (Zero(component) || Exactly(std::abs(component), 1.0))
      {
        return 0.0;
      }
      return (2.0 * Roundoff * std::abs(center)) + swing;
    }

    /**
     *  A box with the sides an ellipse reaches along one axis set where it passes its extremes
     *  there: the origin's coordinate plus or less sqrt(x^2 + y^2), x and y its map's columns'
     *  coordinates, rounded outward, then padded by the square root's rounding and by how far At's
     *  points, the sum of the two columns times a cosine and a sine, may pass it: three steps of a
     *  double of it. With a column's coordinate 0 the extreme is the other's, exactly, and no point
     *  At gives passes it.
     *
     *  @param  ellipse The ellipse.
     *  @param  axis    The axis.
     *  @param  box     The box so far, from its ends.
     *  @return         The box.
     */
    auto Pivoted(const Ellipse &ellipse, Kernel::Axis axis, Box box) -> Box
    {
      const double across = Kernel::Coordinate(ellipse.Map.X, axis);
      const double rising = Kernel::Coordinate(ellipse.Map.Y, axis);
      const double origin = Kernel::Coordinate(ellipse.Map.Origin, axis);
      const double reach  = Trig::Hypot(across, rising);
      const double pad    = Zero(across) || Zero(rising) ? 0.0 : 8.0 * Roundoff * reach;
      const Angle  toward = Radians(Trig::Atan2(rising, across));
      const bool   x      = axis == Kernel::Axis::X;
      if (Arrival(ellipse.Start, ellipse.Sweep, toward))
      {
        (x ? box.Max.X : box.Max.Y) = Raised(Beyond(origin, reach), pad);
      }
      if (Arrival(ellipse.Start, ellipse.Sweep, toward + (Turn / 2.0)))
      {
        (x ? box.Min.X : box.Min.Y) = Lowered(Beyond(origin, -reach), pad);
      }
      return box;
    }

    /**
     *  How far past its start a bulged segment reaches along an axis's way, where its arc passes
     *  its extreme that way, worked out from the chord and the bulge, with no trigonometry.
     *
     *  With c the chord's length, b the bulge's magnitude, s its sign, p the chord's and q its
     *  quarter turn's coordinates along the way, the arc's center lies c (1 - b^2) / (4 b) along
     *  that quarter turn from the chord's middle and its radius is c (1 + b^2) / (4 b). The extreme
     *  center + radius lies on the arc, on the side of the chord the bulge takes, exactly when A <=
     *  b^2 B, A = c + s q and B = c - s q; it lies p / 2 + (A + b^2 B) / (4 b) past the start. A
     *  and B are not negative, and their product is p^2: the smaller is taken as p^2 over the
     *  larger, so that neither cancels.
     *
     *  @param  chord The chord, not nothing.
     *  @param  bulge The bulge, not 0.
     *  @param  way   The way, a unit vector along an axis.
     *  @return       How far, or nothing where the arc does not pass the extreme.
     */
    auto Bulging(Vector chord, double bulge, Vector way) -> std::optional<double>
    {
      const double length = Length(chord);
      const double turn   = bulge > 0.0 ? Cross(chord, way) : -Cross(chord, way);
      const double along  = Dot(chord, way);
      const double square = bulge * bulge;
      const double outer  = turn >= 0.0 ? length + turn : along * along / (length - turn);
      const double inner  = turn >= 0.0 ? along * along / (length + turn) : length - turn;
      if (!(outer <= square * inner))
      {
        return std::nullopt;
      }
      return (along / 2.0) + ((outer + (square * inner)) / (4.0 * std::abs(bulge)));
    }

    /**
     *  The box of a polyline's bulged segment: its ends and the extremes its arc passes, each
     *  padded by four steps of a double of its reach, its start's largest coordinate plus its
     *  chord's length times the bulge's magnitude, or 1 where that is more: no point of the arc
     *  lies further than that from the origin. That holds every point At gives on it, by measure,
     *  not by proof: on 60,000 random bulges from 1e-6 to 1e4, chords from 1e-3 to 1e3 and starts
     *  up to 1e8 away, At's points passed the extremes by 2.1 such steps at most.
     *
     *  @param  step The segment, bulged.
     *  @return      The box.
     */
    auto Bowing(const Step &step) -> Box
    {
      const Vector chord = step.End - step.Start;
      const double size =
          4.0
          * Roundoff
          * (Largest(step.Start - Point{}) + (Length(chord) * std::max(1.0, std::abs(step.Bulge))));
      Box box = Engulfed(Engulfed(Box{}, step.Start, {.X = size, .Y = size}), step.End,
                         {.X = size, .Y = size});
      for (const Vector way : {Vector{.X = 1.0, .Y = 0.0}, Vector{.X = -1.0, .Y = 0.0},
                               Vector{.X = 0.0, .Y = 1.0}, Vector{.X = 0.0, .Y = -1.0}})
      {
        const std::optional<double> reach = Bulging(chord, step.Bulge, way);
        if (reach)
        {
          box = Engulfed(box, step.Start + (way * *reach), {.X = size, .Y = size});
        }
      }
      return box;
    }

    /**
     *  The binomial coefficients of a count, as doubles: exact to 2^53, rounded past it.
     *
     *  @param  count The count, below 2 MaxDegree.
     *  @return       C(count, k) for k from 0 to the count, then nothing.
     */
    auto Pascal(std::size_t count) -> std::array<double, Twice>
    {
      std::array<double, Twice> row{};
      row.at(0) = 1.0;
      for (std::size_t index = 0; index < count; ++index)
      {
        row.at(index + 1) =
            row.at(index) * static_cast<double>(count - index) / static_cast<double>(index + 1);
      }
      return row;
    }

    /**
     *  Where a piece's coordinate along an axis turns: the roots of the numerator of its
     *  derivative. With one weight throughout it is the coordinate's differences, of degree d - 1;
     *  else X' W - X W' for X and W its homogeneous coordinate and weight, of degree 2 d - 1 in
     *  Bernstein form, each product of Bernstein polynomials taken coefficient by coefficient: (f
     *  g)_k = sum over i + j = k of C(d - 1, i) C(d, j) / C(2 d - 1, k) f_i g_j.
     *
     *  @param  piece The piece, of degree 1 at least.
     *  @param  axis  The axis.
     *  @return       The roots, from 0 to 1.
     */
    auto Turns(const Kernel::Piece &piece, Kernel::Axis axis) -> List<double>
    {
      const std::size_t                  degree = piece.Degree;
      const std::span<const Homogeneous> poles  = Kernel::Poles(piece);
      std::array<double, Kernel::Terms>  cells{};
      const bool even = std::ranges::all_of(poles, [&poles](const Homogeneous &pole) -> bool
                                            { return Exactly(pole.W, poles[0].W); });
      if (even)
      {
        for (std::size_t index = 0; index < degree; ++index)
        {
          cells.at(index) =
              Kernel::Coordinate(poles[index + 1], axis) - Kernel::Coordinate(poles[index], axis);
        }
        return Kernel::Roots(std::span{cells}.first(degree), Kernel::Budget);
      }
      const std::array<double, Twice> lower = Pascal(degree - 1);
      const std::array<double, Twice> upper = Pascal(degree);
      const std::array<double, Twice> whole = Pascal((2 * degree) - 1);
      for (std::size_t first = 0; first < degree; ++first)
      {
        const double slope =
            Kernel::Coordinate(poles[first + 1], axis) - Kernel::Coordinate(poles[first], axis);
        const double weight = poles[first + 1].W - poles[first].W;
        for (std::size_t second = 0; second <= degree; ++second)
        {
          cells.at(first + second) +=
              lower.at(first)
              * upper.at(second)
              / whole.at(first + second)
              * ((slope * poles[second].W) - (weight * Kernel::Coordinate(poles[second], axis)));
        }
      }
      return Kernel::Roots(std::span{cells}.first(2 * degree), Kernel::Budget);
    }

    /**
     *  A mix of two homogeneous points, as de Boor's algorithm takes it, and a bound on its error:
     *  the inputs' errors carried by the weights, and its own roundings, each found exactly as the
     *  roundoff of its product or its sum, so that a mix that rounds nothing, at a knot, adds
     *  nothing.
     *
     *  @param  values The two points.
     *  @param  errors Their errors.
     *  @param  ratio  How far toward the second.
     *  @return        The mix, as Mix gives it, and its error.
     */
    auto Mixed(const std::array<Homogeneous, 2> &values, const std::array<Homogeneous, 2> &errors,
               double ratio) -> std::array<Homogeneous, 2>
    {
      const double rest  = 1.0 - ratio;
      Homogeneous  mixed = {};
      Homogeneous  error = {};
      for (double Homogeneous::*const part : {&Homogeneous::X, &Homogeneous::Y, &Homogeneous::W})
      {
        const Exact::Two lower = Exact::TwoProduct(values[0].*part, rest);
        const Exact::Two upper = Exact::TwoProduct(values[1].*part, ratio);
        const Exact::Two sum   = Exact::TwoSum(lower.Big, upper.Big);
        mixed.*part            = sum.Big;
        error.*part = (rest * errors[0].*part)
                    + (ratio * errors[1].*part)
                    + (std::abs(lower.Small) + std::abs(upper.Small) + std::abs(sum.Small));
      }
      return {mixed, error};
    }

    /**
     *  A spline's point at a parameter over a span, as At works it out, with a bound on its
     *  rounding. The fractions' own rounding is left out: it moves the point off the spline, by
     *  what Straying allows for it twice, at this point and at any other At gives.
     *
     *  @tparam CountT How many cells to take room for: more than the degree.
     *  @param  spline The spline, valid.
     *  @param  span   The index of the knot that starts the span.
     *  @param  along  The parameter, within the span.
     *  @return        The point and its bound.
     */
    template <std::size_t CountT>
    auto Reckon(const Spline &spline, std::size_t span, double along) -> Reckoning
    {
      const auto                      degree = static_cast<std::size_t>(spline.Degree);
      const int                       power  = Leveling(spline, span);
      std::array<Homogeneous, CountT> values{};
      std::array<Homogeneous, CountT> errors{};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        values.at(index) = Pole(spline, span - degree + index, power);
      }
      for (std::size_t round = 1; round <= degree; ++round)
      {
        for (std::size_t index = degree; index >= round; --index)
        {
          const double                     start = spline.Knots[span - degree + index];
          const double                     end   = spline.Knots[span + 1 + index - round];
          const std::array<Homogeneous, 2> mixed =
              Mixed({values.at(index - 1), values.at(index)},
                    {errors.at(index - 1), errors.at(index)}, (along - start) / (end - start));
          values.at(index) = mixed[0];
          errors.at(index) = mixed[1];
        }
      }
      const Homogeneous last  = values.at(degree);
      const Homogeneous wrong = errors.at(degree);
      const Point       point = Projected(last);
      const double      under = (last.W - wrong.W) * (1.0 - 0x1p-40);
      return {.At    = point,
              .Error = {.X = ((wrong.X + (std::abs(point.X) * wrong.W)) / under)
                           + (Roundoff * std::abs(point.X)),
                        .Y = ((wrong.Y + (std::abs(point.Y) * wrong.W)) / under)
                           + (Roundoff * std::abs(point.Y))}};
    }

    /**
     *  A spline's point at a parameter over a span, as At works it out, with a bound on its
     *  rounding: Reckon in as many cells as the degree asks.
     *
     *  @param  spline The spline, valid.
     *  @param  span   The index of the knot that starts the span.
     *  @param  along  The parameter, within the span.
     *  @return        The point and its bound.
     */
    auto Reckoned(const Spline &spline, std::size_t span, double along) -> Reckoning
    {
      return spline.Degree < Scant ? Reckon<Scant>(spline, span, along)
                                   : Reckon<MaxDegree + 1>(spline, span, along);
    }

    /**
     *  How far a point At gives on a span of a spline may stray from the spline, along each axis,
     *  whatever the parameter, and how far the point Reckoned works out at a turn may besides its
     *  own bound, which takes only the mixes' roundings.
     *
     *  De Boor's mixes over a degree d are off by at most gamma(2 d) times the sum of their inputs'
     *  magnitudes weighted by the coefficients the fractions give, as de Casteljau's is (Farouki
     *  and Rajan): so the homogeneous coordinate X by gamma(2 d) sum |w_i x_i| c_i, and the weight
     *  W by gamma(2 d) W, the weights being positive. Their quotient is then off by gamma(2 d) (sum
     *  |w_i x_i| c_i / W + |x|), its own rounding and the poles' products' own: at most (4 d + 2)
     *  steps of a double of the largest |x_i| of the span's poles, whatever the weights.
     *
     *  The fractions round too, each mix's by its own few steps, which moves the point off the
     *  curve, not along it: each coefficient c_i is off by D_i, the |D_i| summing to 6.0002 d u at
     *  most, and the point by sum w_i D_i (x_i - x) / sum w_i c_i, at most (w_max / w_min) times
     *  that times the poles' spread (tests/proof/pads.py). At's point and Reckoned's are each off
     *  so: 12.5 d (w_max / w_min) steps of the spread cover both. 2^-40 of all that covers the
     *  second-order excess.
     *
     *  @param  spline The spline, valid.
     *  @param  span   The index of the knot that starts the span.
     *  @return        The bound along each axis.
     */
    auto Straying(const Spline &spline, std::size_t span) -> Vector
    {
      const auto degree  = static_cast<std::size_t>(spline.Degree);
      Vector     largest = {};
      Point      least   = spline.Poles[span];
      Point      most    = spline.Poles[span];
      double     ratio   = 1.0;
      if (!spline.Weights.empty())
      {
        const auto [lightest, heaviest] =
            std::ranges::minmax(std::span{spline.Weights}.subspan(span - degree, degree + 1));
        ratio = heaviest / lightest;
      }
      for (std::size_t index = span - degree; index <= span; ++index)
      {
        const Point pole = spline.Poles[index];
        largest          = {.X = std::max(largest.X, std::abs(pole.X)),
                            .Y = std::max(largest.Y, std::abs(pole.Y))};
        least            = {.X = std::min(least.X, pole.X), .Y = std::min(least.Y, pole.Y)};
        most             = {.X = std::max(most.X, pole.X), .Y = std::max(most.Y, pole.Y)};
      }
      const auto   power  = static_cast<double>(degree);
      const double mixes  = ((4.0 * power) + 2.0) * Roundoff * (1.0 + 0x1p-40);
      const double moving = 12.5 * power * ratio * Roundoff * (1.0 + 0x1p-40);
      return (largest * mixes) + ((most - least) * moving);
    }

    /**
     *  Whether a box holds another grown by a pad either way, rounded outward.
     *
     *  @param  box   The box.
     *  @param  inner The other.
     *  @param  pad   The pad along each axis.
     *  @return       Whether it does.
     */
    auto Embraced(const Box &box, const Box &inner, Vector pad) -> bool
    {
      return Contains(box, {.X = Lowered(inner.Min.X, pad.X), .Y = Lowered(inner.Min.Y, pad.Y)})
          && Contains(box, {.X = Raised(inner.Max.X, pad.X), .Y = Raised(inner.Max.Y, pad.Y)});
    }

    /**
     *  The box of a span's poles: the span's Bézier control points, each a convex mix of them by
     *  fractions worked out from 0 to 1, lie in it but for what their mixes round, which Straying's
     *  pad bounds; a span whose poles' box, grown by eight pads, the box holds has its control
     *  points' box, grown by two, held as well, a pad and more to spare for the roundings of the
     *  growing, and is passed over as it would be once worked out.
     *
     *  @param  spline The spline.
     *  @param  span   The index of the knot that starts the span.
     *  @return        The box.
     */
    auto Poled(const Spline &spline, std::size_t span) -> Box
    {
      const auto degree = static_cast<std::size_t>(spline.Degree);
      Box        box    = {.Min = spline.Poles[span], .Max = spline.Poles[span]};
      for (std::size_t index = span - degree; index < span; ++index)
      {
        box = Grown(box, spline.Poles[index]);
      }
      return box;
    }
  } // namespace

  auto Extents(const Line &line) -> Box
  {
    return Grown(Grown(Box{}, line.Start), line.End);
  }

  auto Extents(const Ray &ray) -> Box
  {
    const std::array<double, 2> x = Ranged(ray.Start.X, ray.Way.X);
    const std::array<double, 2> y = Ranged(ray.Start.Y, ray.Way.Y);
    return {.Min = {.X = x[0], .Y = y[0]}, .Max = {.X = x[1], .Y = y[1]}};
  }

  auto Extents(const Xline &xline) -> Box
  {
    const std::array<double, 2> x = Endwise(xline.Through.X, xline.Way.X);
    const std::array<double, 2> y = Endwise(xline.Through.Y, xline.Way.Y);
    return {.Min = {.X = x[0], .Y = y[0]}, .Max = {.X = x[1], .Y = y[1]}};
  }

  auto Extents(const Arc &arc) -> Box
  {
    const double swing =
        Roundoff * std::abs(arc.Radius) * (8.0 + Spanning(arc.Start) + Spanning(arc.Sweep));
    Box box;
    for (const double along : {0.0, 1.0})
    {
      const Vector way = Direction(arc.Start + (arc.Sweep * along));
      box              = Engulfed(
          box, arc.Center + (way * arc.Radius),
          {.X = Wobbled(arc.Center.X, way.X, swing), .Y = Wobbled(arc.Center.Y, way.Y, swing)});
    }
    for (const double degrees : {0.0, 90.0, 180.0, 270.0})
    {
      if (Arrival(arc.Start, arc.Sweep, Degrees(degrees)))
      {
        box = Reaching(box, arc.Center, Direction(Degrees(degrees)) * arc.Radius);
      }
    }
    return box;
  }

  auto Extents(const Ellipse &ellipse) -> Box
  {
    const double size   = Largest(ellipse.Map.X) + Largest(ellipse.Map.Y);
    const double angles = Spanning(ellipse.Start) + Spanning(ellipse.Sweep);
    const auto   pad    = [size, angles](double origin) -> double
    {
      return Roundoff * ((2.0 * (std::abs(origin) + (5.0 * size))) + (size * angles));
    };
    Box box;
    for (const double along : {0.0, 1.0})
    {
      box = Engulfed(box, At(ellipse, along),
                     {.X = pad(ellipse.Map.Origin.X), .Y = pad(ellipse.Map.Origin.Y)});
    }
    return Pivoted(ellipse, Kernel::Axis::Y, Pivoted(ellipse, Kernel::Axis::X, box));
  }

  auto Extents(const Polyline &polyline) -> Box
  {
    Box box;
    for (const Vertex &vertex : polyline.Vertices)
    {
      box = Grown(box, vertex.Position);
    }
    const std::size_t count = Count(polyline);
    for (std::size_t index = 0; index < count; ++index)
    {
      const Step step = Locate(polyline, static_cast<double>(index));
      if (!Straightness(step))
      {
        box = Joined(box, Bowing(step));
      }
    }
    return box;
  }

  auto Extents(const Bezier::Quadratic &quadratic) -> Box
  {
    Box          box  = Grown(Grown(Box{}, quadratic.Start), quadratic.End);
    const Vector bend = (quadratic.Start - quadratic.Control) + (quadratic.End - quadratic.Control);
    const Vector lean = quadratic.Start - quadratic.Control;
    const Vector pad  = Wandering({quadratic.Start, quadratic.Control, quadratic.End});
    for (const double along : {Roots(0.0, bend.X, -lean.X)[0], Roots(0.0, bend.Y, -lean.Y)[0]})
    {
      box = Turnpoint(box, {.At = At(quadratic, along), .Error = pad}, along);
    }
    return box;
  }

  auto Extents(const Bezier::Cubic &cubic) -> Box
  {
    Box          box    = Grown(Grown(Box{}, cubic.Start), cubic.End);
    const Vector first  = cubic.Out - cubic.Start;
    const Vector second = cubic.In - cubic.Out;
    const Vector third  = cubic.End - cubic.In;
    const Vector square = (first - second) + (third - second);
    const Vector linear = (second - first) * 2.0;
    const Vector pad    = Wandering({cubic.Start, cubic.Out, cubic.In, cubic.End});
    for (const std::array<double, 2> &roots :
         {Roots(square.X, linear.X, first.X), Roots(square.Y, linear.Y, first.Y)})
    {
      for (const double along : roots)
      {
        box = Turnpoint(box, {.At = At(cubic, along), .Error = pad}, along);
      }
    }
    return box;
  }

  auto Extents(const Spline &spline) -> Box
  {
    if (!Valid(spline))
    {
      return {};
    }
    Box         box;
    const Basis basis  = Spanned(spline);
    const auto  degree = static_cast<std::size_t>(spline.Degree);
    Vector      before = {};
    for (std::size_t span = degree; span < spline.Poles.size(); ++span)
    {
      if (spline.Knots[span] < spline.Knots[span + 1])
      {
        const Vector    pad       = Straying(spline, span);
        const Reckoning reckoning = Reckoned(spline, span, spline.Knots[span]);
        box =
            Engulfed(box, reckoning.At,
                     reckoning.Error
                         + Vector{.X = std::max(pad.X, before.X), .Y = std::max(pad.Y, before.Y)});
        before = pad;
      }
    }
    const std::size_t last      = Span(basis, spline.Knots[spline.Poles.size()]);
    const Reckoning   reckoning = Reckoned(spline, last, spline.Knots[spline.Poles.size()]);
    box                         = Engulfed(box, reckoning.At, reckoning.Error + before);
    for (std::size_t span = degree; span < spline.Poles.size(); ++span)
    {
      if (!(spline.Knots[span] < spline.Knots[span + 1]))
      {
        continue;
      }
      const Vector pad = Straying(spline, span);
      if (Embraced(box, Poled(spline, span), pad * 8.0))
      {
        continue;
      }
      const Kernel::Piece piece = Kernel::Section(spline, span);
      if (Embraced(box, Kernel::Hull(piece), pad * 2.0))
      {
        continue;
      }
      const Range stretch = {.Start = spline.Knots[span], .End = spline.Knots[span + 1]};
      for (const Kernel::Axis axis : {Kernel::Axis::X, Kernel::Axis::Y})
      {
        for (const double root : Turns(piece, axis))
        {
          const Reckoning turn = Reckoned(spline, span, Kernel::Across(stretch, root));
          box                  = Engulfed(box, turn.At, turn.Error + pad);
        }
      }
    }
    return box;
  }

  auto Extents(const Curve &curve) -> Box
  {
    return std::visit([](const auto &shape) -> Box { return Extents(shape); }, curve);
  }
} // namespace Geoxy
