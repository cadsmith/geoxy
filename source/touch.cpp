/***************************************************************************************************
 *  @file      source/touch.cpp
 *  @brief     Where two curves come nearest about a touch, polished in double-double and rounded
 *             once: one point, whichever curve is first and however the plane is turned.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "touch.hpp"

#include "basis.hpp"
#include "bow.hpp"
#include "crossing.hpp"
#include "trig.hpp"
#include "twofold.hpp"
#include "walk.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
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
#include <optional>
#include <variant>

namespace Geoxy::Touches
{
  namespace
  {
    using Twofold::Add;
    using Twofold::Difference;
    using Twofold::Divide;
    using Twofold::Double;
    using Twofold::Lifted;
    using Twofold::Multiply;
    using Twofold::Root;
    using Twofold::Two;
    using Twofold::TwoProduct;
    using Twofold::TwoSum;

    /**
     *  A point or a vector in double-double.
     */
    struct Spot
    {
      Double X = {}; /**< Along x. */
      Double Y = {}; /**< Along y. */
    };

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
     *  A difference of two doubles, exactly.
     *
     *  @param  first  The double.
     *  @param  second The double taken from it.
     *  @return        The difference.
     */
    auto Apart(double first, double second) -> Double
    {
      const Two difference = TwoSum(first, -second);
      return Add(Lifted(difference.High), Lifted(difference.Low));
    }

    /**
     *  The way from one point to another, exactly.
     *
     *  @param  from   Where it starts.
     *  @param  target Where it ends.
     *  @return        The way.
     */
    auto Between(Point from, Point target) -> Spot
    {
      return {.X = Apart(target.X, from.X), .Y = Apart(target.Y, from.Y)};
    }

    /**
     *  A sum of two points or vectors.
     *
     *  @param  first  One.
     *  @param  second The other.
     *  @return        The sum.
     */
    auto Plus(const Spot &first, const Spot &second) -> Spot
    {
      return {.X = Add(first.X, second.X), .Y = Add(first.Y, second.Y)};
    }

    /**
     *  A difference of two points or vectors.
     *
     *  @param  first  The one.
     *  @param  second The one taken from it.
     *  @return        The difference.
     */
    auto Minus(const Spot &first, const Spot &second) -> Spot
    {
      return {.X = Difference(first.X, second.X), .Y = Difference(first.Y, second.Y)};
    }

    /**
     *  A vector times a number.
     *
     *  @param  spot   The vector.
     *  @param  factor The number.
     *  @return        The product.
     */
    auto Times(const Spot &spot, Double factor) -> Spot
    {
      return {.X = Multiply(spot.X, factor), .Y = Multiply(spot.Y, factor)};
    }

    /**
     *  The dot product of two vectors.
     *
     *  @param  first  One.
     *  @param  second The other.
     *  @return        The product.
     */
    auto Dotted(const Spot &first, const Spot &second) -> Double
    {
      return Add(Multiply(first.X, second.X), Multiply(first.Y, second.Y));
    }

    /**
     *  The cross product of two vectors.
     *
     *  @param  first  One.
     *  @param  second The other.
     *  @return        The product.
     */
    auto Across(const Spot &first, const Spot &second) -> Double
    {
      return Difference(Multiply(first.X, second.Y), Multiply(first.Y, second.X));
    }

    /**
     *  A vector's leading parts.
     *
     *  @param  spot The vector.
     *  @return      It in doubles.
     */
    auto Leading(const Spot &spot) -> Vector
    {
      return {.X = spot.X.High, .Y = spot.Y.High};
    }

    /**
     *  A point between two, part of the way.
     *
     *  @param  start Where it starts.
     *  @param  end   Where it ends.
     *  @param  along How far.
     *  @return       The point.
     */
    auto Lerped(const Spot &start, const Spot &end, Double along) -> Spot
    {
      return Plus(start, Times(Minus(end, start), along));
    }

    /**
     *  A vector made one long.
     *
     *  @param  spot The vector, not nothing.
     *  @return      It, one long.
     */
    auto Unit(const Spot &spot) -> Spot
    {
      const Double length = Root(Dotted(spot, spot));
      return {.X = Divide(spot.X, length), .Y = Divide(spot.Y, length)};
    }

    /**
     *  A curve's point and slope at a place in double-double, and how its slope turns there in
     *  doubles: what Newton's method on where two curves come nearest asks of each.
     */
    struct Jet
    {
      Spot   At   = {}; /**< The point.                */
      Spot   Way  = {}; /**< The slope.                */
      Vector Bend = {}; /**< How fast the slope turns. */
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
     *  A mix of two homogeneous points.
     *
     *  @param  first  The one weighed by the rest.
     *  @param  second The one weighed by the fraction.
     *  @param  weighs The rest and the fraction.
     *  @return        The mix.
     */
    auto Mixed(const Weighted &first, const Weighted &second, const std::array<Double, 2> &weighs)
        -> Weighted
    {
      const auto mix = [&weighs](Double one, Double other) -> Double
      {
        return Add(Multiply(one, weighs[0]), Multiply(other, weighs[1]));
      };
      return {
          .X = mix(first.X, second.X), .Y = mix(first.Y, second.Y), .W = mix(first.W, second.W)};
    }

    /**
     *  A spline's point and slope at a place, by de Boor's algorithm in double-double on its
     *  homogeneous control points: the slope from the two points the last round mixes, as their
     *  difference times the degree over the span, by the quotient rule where the weights differ.
     *
     *  @param  spline The spline, valid.
     *  @param  along  The place.
     *  @return        The point and the slope.
     */
    auto Jetted(const Spline &spline, Double along) -> std::array<Spot, 2>
    {
      const Basis                         basis  = Spanned(spline);
      const std::size_t                   degree = basis.Degree;
      const std::size_t                   span   = Span(basis, along.High);
      const int                           power  = Leveling(spline, span);
      std::array<Weighted, MaxDegree + 1> cells{};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        const Point  pole = spline.Poles[span - degree + index];
        const double weight =
            spline.Weights.empty() ? 1.0 : std::ldexp(spline.Weights[span - degree + index], power);
        const Two x     = TwoProduct(pole.X, weight);
        const Two y     = TwoProduct(pole.Y, weight);
        cells.at(index) = {.X = {.High = x.High, .Low = x.Low},
                           .Y = {.High = y.High, .Low = y.Low},
                           .W = Lifted(weight)};
      }
      std::array<Weighted, 2> last = {cells.at(degree - 1), cells.at(degree)};
      for (std::size_t round = 1; round <= degree; ++round)
      {
        last = {cells.at(degree - 1), cells.at(degree)};
        for (std::size_t index = degree; index >= round; --index)
        {
          const double start = basis.Knots[span - degree + index];
          const double end   = basis.Knots[span + 1 + index - round];
          const Double width = Apart(end, start);
          cells.at(index)    = Mixed(cells.at(index - 1), cells.at(index),
                                     {Divide(Difference(Lifted(end), along), width),
                                      Divide(Difference(along, Lifted(start)), width)});
        }
      }
      const Weighted point = cells.at(degree);
      const Double   scale = Divide(Lifted(static_cast<double>(degree)),
                                    Apart(basis.Knots[span + 1], basis.Knots[span]));
      const Weighted slope = {.X = Multiply(Difference(last[1].X, last[0].X), scale),
                              .Y = Multiply(Difference(last[1].Y, last[0].Y), scale),
                              .W = Multiply(Difference(last[1].W, last[0].W), scale)};
      const Spot     place = {.X = Divide(point.X, point.W), .Y = Divide(point.Y, point.W)};
      return {place, Spot{.X = Divide(Difference(slope.X, Multiply(place.X, slope.W)), point.W),
                          .Y = Divide(Difference(slope.Y, Multiply(place.Y, slope.W)), point.W)}};
    }

    /**
     *  The sine and the cosine of a small angle in double-double, by their series: six terms, which
     *  under an angle of a thousandth leave less than double-double holds.
     *
     *  @param  angle The angle, in radians, under a thousandth.
     *  @return       Its sine and its cosine.
     */
    auto Waved(Double angle) -> std::array<Double, 2>
    {
      const Double square = Multiply(angle, angle);
      Double       sine   = Lifted(1.0);
      Double       cosine = Lifted(1.0);
      for (int term = 5; term >= 1; --term)
      {
        const double even = 2.0 * term;
        sine = Difference(Lifted(1.0), Multiply(Divide(square, Lifted(even * (even + 1.0))), sine));
        cosine =
            Difference(Lifted(1.0), Multiply(Divide(square, Lifted((even - 1.0) * even)), cosine));
      }
      return {Multiply(angle, sine), cosine};
    }

    /**
     *  A displacement turned by an angle told by its sine and cosine.
     *
     *  @param  spot The displacement.
     *  @param  wave The angle's sine and cosine.
     *  @return      The displacement turned, counterclockwise.
     */
    auto Swung(const Spot &spot, const std::array<Double, 2> &wave) -> Spot
    {
      return {.X = Difference(Multiply(wave[1], spot.X), Multiply(wave[0], spot.Y)),
              .Y = Add(Multiply(wave[0], spot.X), Multiply(wave[1], spot.Y))};
    }

    /**
     *  A polyline's slight arc's point and slope at a part of its sweep, in double-double, from its
     *  ends and its bulge as Along and Tangent work them out, with no center: the point is the
     *  start plus the chord turned by half the sweep left and scaled by the ratio of the half
     *  sines, sin(s / 2) / sin(S / 2); the slope is the chord turned by s - S / 2 and scaled by S /
     *  (2 sin(S / 2)). The quarter sweep is atan b by its series, and sin(S / 2) is 2 b / (1 +
     *  b^2).
     *
     *  @param  step The segment, a remote arc, its bulge under Slight.
     *  @param  part The part of the sweep.
     *  @return      The point and the slope.
     */
    auto Bowed(const Step &step, Double part) -> std::array<Spot, 2>
    {
      const Double bulge  = Lifted(step.Bulge);
      const Two    twice  = TwoProduct(step.Bulge, step.Bulge);
      const Double square = {.High = twice.High, .Low = twice.Low};
      Double       series = Lifted(1.0 / 9.0);
      for (const double odd : {7.0, 5.0, 3.0, 1.0})
      {
        series = Difference(Divide(Lifted(1.0), Lifted(odd)), Multiply(square, series));
      }
      const Double sweep  = Multiply(Lifted(4.0), Multiply(bulge, series));
      const Double half   = Divide(Multiply(Lifted(2.0), bulge), Add(Lifted(1.0), square));
      const Double swept  = Multiply(sweep, part);
      const Double middle = Multiply(sweep, Lifted(0.5));
      const Spot   chord  = Between(step.Start, step.End);
      const Double scale  = Divide(Waved(Multiply(swept, Lifted(0.5)))[0], half);
      const Spot   reached =
          Swung(Times(chord, scale), Waved(Multiply(Difference(swept, sweep), Lifted(0.5))));
      const Double speed = Divide(middle, half);
      return {Plus(Spotted(step.Start), reached),
              Swung(Times(chord, speed), Waved(Difference(swept, middle)))};
    }

    /**
     *  A curve's point and slope at a place in double-double, where the curve is held by its place:
     *  a line, a ray, an endless line, a polyline's straight segment or slight arc, a spline, a
     *  Bézier curve, by de Casteljau's algorithm for the last.
     *
     *  @param  curve The curve, not an arc nor an ellipse, nor a polyline at a bulge not slight.
     *  @param  along The place.
     *  @return       The point and the slope.
     */
    auto Jetted(const Curve &curve, Double along) -> std::array<Spot, 2>
    {
      if (const auto *line = std::get_if<Line>(&curve))
      {
        return {Lerped(Spotted(line->Start), Spotted(line->End), along),
                Between(line->Start, line->End)};
      }
      if (const auto *ray = std::get_if<Ray>(&curve))
      {
        const Spot way = Spotted(Point{} + ray->Way);
        return {Plus(Spotted(ray->Start), Times(way, along)), way};
      }
      if (const auto *xline = std::get_if<Xline>(&curve))
      {
        const Spot way = Spotted(Point{} + xline->Way);
        return {Plus(Spotted(xline->Through), Times(way, along)), way};
      }
      if (const auto *polyline = std::get_if<Polyline>(&curve))
      {
        const Step   step  = Locate(*polyline, along.High);
        const double index = along.High - step.Fraction;
        if (!Straightness(step))
        {
          return Bowed(step, Difference(along, Lifted(index)));
        }
        return {Lerped(Spotted(step.Start), Spotted(step.End), Difference(along, Lifted(index))),
                Between(step.Start, step.End)};
      }
      if (const auto *spline = std::get_if<Spline>(&curve))
      {
        return Jetted(*spline, along);
      }
      if (const auto *quadratic = std::get_if<Bezier::Quadratic>(&curve))
      {
        const Spot control = Spotted(quadratic->Control);
        const Spot one     = Lerped(Spotted(quadratic->Start), control, along);
        const Spot two     = Lerped(control, Spotted(quadratic->End), along);
        return {Lerped(one, two, along), Times(Minus(two, one), Lifted(2.0))};
      }
      const auto &cubic = std::get<Bezier::Cubic>(curve);
      const Spot  out   = Spotted(cubic.Out);
      const Spot  inner = Spotted(cubic.In);
      const Spot  one =
          Lerped(Lerped(Spotted(cubic.Start), out, along), Lerped(out, inner, along), along);
      const Spot two =
          Lerped(Lerped(out, inner, along), Lerped(inner, Spotted(cubic.End), along), along);
      return {Lerped(one, two, along), Times(Minus(two, one), Lifted(3.0))};
    }

    /**
     *  A conic as Newton's method moves along it: its frame, the point of its unit circle it is at,
     *  turned as the method goes, and how a turn tells the conic's place.
     */
    struct Frame
    {
      Spot   Origin = {};  /**< The center.                                       */
      Spot   X      = {};  /**< Where the unit x axis goes.                       */
      Spot   Y      = {};  /**< Where the unit y axis goes.                       */
      Spot   Unit   = {};  /**< The point of the unit circle it is at.            */
      Spot   Start  = {};  /**< The point it started at.                          */
      double From   = 0.0; /**< The conic's place there.                          */
      double Rate   = 0.0; /**< How far round its unit circle a unit of place is. */
    };

    /**
     *  A curve as Newton's method moves along it near a touch: by its place, or by the point of its
     *  unit circle for a conic; and the stretch of its place the touch must stay in.
     */
    struct Chart
    {
      const Curve *Shape = nullptr;    /**< The curve.                 */
      Range        Part  = {};         /**< Where its place must stay. */
      bool         Wraps = false;      /**< Whether it closes on itself, its place
                                            taken round past its ends.              */
      Double               Along = {}; /**< Its place, held by its place. */
      std::optional<Frame> Round = {}; /**< Its frame, a conic.           */
    };

    /**
     *  A conic's frame at a point of it.
     *
     *  @param  axes   Its center and the two axes of its map, in double-double.
     *  @param  unit   The point of its unit circle there, as its angle gives it.
     *  @param  places Its place there and how far round its unit circle a unit of place is.
     *  @return        The frame.
     */
    auto Framed(const std::array<Spot, 3> &axes, const Spot &unit,
                const std::array<double, 2> &places) -> Frame
    {
      const Spot start = Unit(unit);
      return {.Origin = axes[0],
              .X      = axes[1],
              .Y      = axes[2],
              .Unit   = start,
              .Start  = start,
              .From   = places[0],
              .Rate   = places[1]};
    }

    /**
     *  A curve as Newton's method moves along it from a place: an arc, an ellipse and a polyline's
     *  bulge by the point of their unit circle, a bulge's circle worked out in double-double from
     *  its vertices and its bulge, as its own numbers give it; any other curve by its place. A
     *  polyline's place stays on its segment; an arc's or an ellipse's that closes on itself goes
     *  round past its seam.
     *
     *  @param  curve The curve.
     *  @param  along The place.
     *  @return       The chart; nothing for a conic of no sweep.
     */
    auto Charted(const Curve &curve, double along) -> std::optional<Chart>
    {
      const Range domain = Domain(curve);
      const bool  closed = At(curve, domain.Start) == At(curve, domain.End);
      Chart       chart = {.Shape = &curve, .Part = domain, .Wraps = false, .Along = Lifted(along)};
      if (const auto *arc = std::get_if<Arc>(&curve))
      {
        const Spot radius = Spotted(Point{.X = arc->Radius, .Y = 0.0});
        chart.Round =
            Framed({Spotted(arc->Center), radius, Spot{.X = Lifted(0.0), .Y = Lifted(arc->Radius)}},
                   Spotted(Point{} + Direction(arc->Start + (arc->Sweep * along))),
                   {along, arc->Sweep.InRadians()});
        chart.Wraps = closed;
      }
      else if (const auto *ellipse = std::get_if<Ellipse>(&curve))
      {
        const Matrix &map = ellipse->Map;
        chart.Round =
            Framed({Spotted(map.Origin), Spotted(Point{} + map.X), Spotted(Point{} + map.Y)},
                   Spotted(Point{} + Direction(ellipse->Start + (ellipse->Sweep * along))),
                   {along, ellipse->Sweep.InRadians()});
        chart.Wraps = closed;
      }
      else if (const auto *polyline = std::get_if<Polyline>(&curve))
      {
        const Step    step    = Locate(*polyline, along);
        const double  index   = along - step.Fraction;
        const Segment segment = Nth(*polyline, static_cast<std::size_t>(index));
        chart.Part            = {.Start = index, .End = index + 1.0};
        /**
         *  A remote arc's center lies too far to frame it by: it is held by its place.
         */
        const auto *bow =
            Remote(step.Start, step.End, step.Bulge) ? nullptr : std::get_if<Arc>(&segment);
        if (bow != nullptr)
        {
          /**
           *  A bulge b from P0 to P1 bows about the center (P0 + P1) / 2 + n (1 - b^2) / (4 b), n
           *  the chord turned a quarter to the left, at the radius |P1 - P0| (1 + b^2) / (4 |b|).
           */
          const Spot   chord  = Between(step.Start, step.End);
          const Two    square = TwoProduct(step.Bulge, step.Bulge);
          const Double bend   = {.High = square.High, .Low = square.Low};
          const Spot   middle = Plus(Spotted(step.Start), Times(chord, Lifted(0.5)));
          const Spot   normal = {.X = Twofold::Negated(chord.Y), .Y = chord.X};
          const Double reach  = Divide(Difference(Lifted(1.0), bend), Lifted(4.0 * step.Bulge));
          const Spot   center = Plus(middle, Times(normal, reach));
          const Double radius = Divide(Multiply(Root(Dotted(chord, chord)), Add(Lifted(1.0), bend)),
                                       Lifted(4.0 * std::abs(step.Bulge)));
          chart.Round         = Framed(
              {center, Spot{.X = radius, .Y = Lifted(0.0)}, Spot{.X = Lifted(0.0), .Y = radius}},
              Spotted(Point{} + Direction(bow->Start + (bow->Sweep * step.Fraction))),
              {along, bow->Sweep.InRadians()});
        }
      }
      /**
       *  A conic of no sweep is a point: no turn of it is a place.
       */
      if (chart.Round && !(std::abs(chart.Round->Rate) > 0.0))
      {
        return std::nullopt;
      }
      return chart;
    }

    /**
     *  A chart's place.
     *
     *  @param  chart The chart.
     *  @return       The curve's place: a conic's from how far its unit circle's point has turned.
     */
    auto Placed(const Chart &chart) -> double
    {
      if (!chart.Round)
      {
        return chart.Along.High;
      }
      const Frame &frame = *chart.Round;
      const double turned =
          Trig::Atan2(Across(frame.Start, frame.Unit).High, Dotted(frame.Start, frame.Unit).High);
      return frame.From + (turned / frame.Rate);
    }

    /**
     *  A chart's point, slope and bend.
     *
     *  @param  chart The chart.
     *  @return       Them: a conic's by the turn of its unit circle's point, in radians.
     */
    auto Jetted(const Chart &chart) -> Jet
    {
      if (const std::optional<Frame> &frame = chart.Round)
      {
        const Spot   along  = Times(frame->X, frame->Unit.X);
        const Spot   across = Times(frame->Y, frame->Unit.Y);
        const Spot   toward = Plus(along, across);
        const Vector bend   = Leading(toward);
        return {.At   = Plus(frame->Origin, toward),
                .Way  = Minus(Times(frame->Y, frame->Unit.X), Times(frame->X, frame->Unit.Y)),
                .Bend = {.X = -bend.X, .Y = -bend.Y}};
      }
      const auto [point, way] = Jetted(*chart.Shape, chart.Along);
      /**
       *  The bend is only for Newton's method's steps, not for where it settles: the slope a step
       *  of 2^-20 either way changes by.
       */
      const double place = chart.Along.High;
      const double step  = 0x1p-20 * std::max(1.0, std::abs(place));
      const Vector turn =
          Derivative(*chart.Shape, place + step) - Derivative(*chart.Shape, place - step);
      return {.At = point, .Way = way, .Bend = turn / (2.0 * step)};
    }

    /**
     *  A chart moved by a step of Newton's method: its place, or its unit circle's point turned.
     *
     *  @param chart The chart.
     *  @param step  The step: of the place, or in radians.
     */
    void Move(Chart &chart, double step)
    {
      if (std::optional<Frame> &frame = chart.Round)
      {
        const Spot normal = {.X = Twofold::Negated(frame->Unit.Y), .Y = frame->Unit.X};
        frame->Unit       = Unit(Plus(frame->Unit, Times(normal, Lifted(step))));
        return;
      }
      chart.Along = Add(chart.Along, Lifted(step));
    }

    /**
     *  Whether a chart's place lies in the stretch the touch must stay in: anywhere on a conic that
     *  closes on itself, as long as it is a number.
     *
     *  @param  chart The chart.
     *  @return       Whether it does.
     */
    auto Inside(const Chart &chart) -> bool
    {
      const double place = Placed(chart);
      return chart.Wraps ? !std::isnan(place)
                         : chart.Part.Start <= place && place <= chart.Part.End;
    }

    /**
     *  A chart's place, taken round into its domain on a conic that closes on itself.
     *
     *  @param  chart The chart.
     *  @return       The place.
     */
    auto Kept(const Chart &chart) -> double
    {
      const double place = Placed(chart);
      return chart.Wraps ? place - std::floor(place) : place;
    }

    /**
     *  The most rounds of Newton's method, and the part of the size about the touch below which a
     *  step is lost in double-double's width: 2^16 of its steps, far below a step of the double the
     *  point halfway is rounded to.
     */
    constexpr int    Rounds = 32;
    constexpr double Settle = 0x1p-88;

    /**
     *  A coordinate of the point halfway, 0 where double-double's width about the touch leaves it
     *  no number: under 2^-80 of the size about it, where what is left is the rounding of the way
     *  it was worked out, and its bits would hang on which curve was first.
     *
     *  @param  value The coordinate.
     *  @param  size  The size about the touch.
     *  @return       The coordinate, or 0.
     */
    auto Cleared(double value, double size) -> double
    {
      return std::abs(value) <= 0x1p-80 * size ? 0.0 : value;
    }

    /**
     *  How much two curves must bend apart for where they come nearest to be one place: their
     *  curvatures 2^-26 of themselves apart at the least, as Newton's system weighs them. Two
     *  lines, two arcs of one circle or two curves along each other come nearest all along.
     */
    constexpr double Bending = 0x1p-26;
  } // namespace

  auto Nearest(const std::array<const Curve *, 2> &curves, const std::array<double, 2> &places,
               double length) -> std::optional<Crossings::Crossed>
  {
    std::optional<Chart> one   = Charted(*curves[0], places[0]);
    std::optional<Chart> other = Charted(*curves[1], places[1]);
    if (!one || !other)
    {
      return std::nullopt;
    }
    for (int round = 0; round < Rounds; ++round)
    {
      const Jet  first  = Jetted(*one);
      const Jet  second = Jetted(*other);
      const Spot apart  = Minus(first.At, second.At);
      /**
       *  The tangents parallel and the way between the points square to the first's: their cross
       *  product and the dot product worked out in double-double, Newton's step in doubles.
       */
      const double parallel = Across(first.Way, second.Way).High;
      const double square   = Dotted(apart, first.Way).High;
      const Vector way      = Leading(first.Way);
      const Vector back     = Leading(second.Way);
      const Vector gap      = Leading(apart);
      const double turn     = Cross(first.Bend, back);
      const double turning  = Cross(way, second.Bend);
      const double slide    = Dot(way, way) + Dot(gap, first.Bend);
      const double sliding  = -Dot(back, way);
      const double det      = (turn * sliding) - (turning * slide);
      if (!(std::abs(det) > Bending * (std::abs(turn * sliding) + std::abs(turning * slide))))
      {
        return std::nullopt;
      }
      const double ahead  = ((turning * square) - (sliding * parallel)) / det;
      const double behind = ((slide * parallel) - (turn * square)) / det;
      Move(*one, ahead);
      Move(*other, behind);
      if (!Inside(*one) || !Inside(*other))
      {
        return std::nullopt;
      }
      const Vector here = Leading(first.At);
      const double size = std::max({std::abs(here.X), std::abs(here.Y), Length(way), Length(back)});
      if (std::abs(ahead) * Length(way)      <= Settle * size
          && std::abs(behind) * Length(back) <= Settle * size)
      {
        const Jet    near  = Jetted(*one);
        const Jet    far   = Jetted(*other);
        const Vector split = Leading(Minus(near.At, far.At));
        if (!(Length(split) <= length))
        {
          return std::nullopt;
        }
        const std::array<double, 2> kept   = {Kept(*one), Kept(*other)};
        const Vector                middle = Leading(Times(Plus(near.At, far.At), Lifted(0.5)));
        return Crossings::Crossed{
            .At    = {.X = Cleared(middle.X, size), .Y = Cleared(middle.Y, size)},
            .Along = kept,
            .On    = {At(*curves[0], kept[0]), At(*curves[1], kept[1])}};
      }
    }
    return std::nullopt;
  }
} // namespace Geoxy::Touches
