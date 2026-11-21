/***************************************************************************************************
 *  @file      source/swift/bernstein.cpp
 *  @brief     A line, a ray or an endless line met with a quadratic or cubic Bézier curve: the
 *             line's equation on the curve held by its Bernstein coefficients, its roots isolated
 *             in doubles under one running bound, each polished by one double-double Newton step
 *             whose error the mean value theorem bounds, and the point and places proved the
 *             doubles nearest the true ones.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "straight.hpp"
#include "sure.hpp"
#include "swift/casteljau.hpp"
#include "swift/lane.hpp"
#include "swift/roots.hpp"
#include "swift/section.hpp"
#include "twofold.hpp"

#include "geoxy/domain.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <span>

namespace Geoxy::Swift
{
  namespace
  {
    using Sure::Number;
  } // namespace

  auto Bernstein(const Lane &lane, const Section &section, const Tolerance &tolerance)
      -> std::optional<Spots>
  {
    const std::span<const Point> poles  = std::span{section.Poles}.first(section.Degree + 1);
    const std::size_t            degree = section.Degree;
    const double                 reach  = Reach(tolerance);
    const double                 speedD = lane.Speed * (1.0 + 0x1p-50);
    if (!std::ranges::all_of(poles, [](Point pole) -> bool { return Held({pole.X, pole.Y}); }))
    {
      return std::nullopt;
    }
    /**
     *  f(t) = d x (B(t) - P) has Bernstein coefficients d x (B_i - P), exact differences and
     *  double-double products; |f| / |d| is B(t)'s distance from the carrier. The reach T is eight
     *  tolerances and the margin, times |d|. Control points within a slack s of the true ones move
     *  each coefficient by (|dx| + |dy|) s at most.
     */
    const Number                    dx       = lane.Way[0];
    const Number                    dy       = lane.Way[1];
    const double                    drift    = (Sure::Most(dx) + Sure::Most(dy)) * section.Slack;
    std::array<Number, Highest + 1> equation = {};
    Stretch<Highest + 1>            whole;
    double                          first   = 0.0;
    double                          largest = 0.0;
    for (std::size_t index = 0; index <= degree; ++index)
    {
      const Point pole = poles[index];
      equation.at(index) =
          (dx * Sure::Apart(pole.Y, lane.Start.Y)) - (dy * Sure::Apart(pole.X, lane.Start.X));
      equation.at(index).Bound += drift * Sure::Grow;
      whole.C.at(index)         = equation.at(index).Value.High;
      first   = std::max(first, (equation.at(index).Bound + std::abs(equation.at(index).Value.Low))
                                    * Sure::Grow);
      largest = std::max(largest, std::abs(equation.at(index).Value.High));
    }
    const double                            cut   = reach * speedD;
    const std::optional<Roots<Highest + 1>> roots = Isolated(whole, {first, largest, cut}, degree);
    if (!roots)
    {
      return std::nullopt;
    }
    /**
     *  The curve's speed, at most n times its longest control step, and its bend, n (n
     *  - 1) times its largest second difference; f'' likewise, from f's own, with what their
     *    rounding adds; the slack moving a step's coordinates by twice itself, a second
     *    difference's by four times.
     */
    double speed = 0.0;
    double bend  = 0.0;
    double turn  = 0.0;
    for (std::size_t index = 0; index < degree; ++index)
    {
      speed = std::max(speed, Length(poles[index + 1] - poles[index]));
    }
    for (std::size_t index = 0; index + 2 <= degree; ++index)
    {
      const Vector second =
          (poles[index + 2] - poles[index + 1]) - (poles[index + 1] - poles[index]);
      bend = std::max({bend, std::abs(second.X), std::abs(second.Y)});
      turn = std::max(
          turn, std::abs(whole.C.at(index + 2) - (2.0 * whole.C.at(index + 1)) + whole.C.at(index))
                    + (8.0 * first)
                    + (8.0 * Unit * largest));
    }
    const auto order  = static_cast<double>(degree);
    speed            += 2.0 * std::numbers::sqrt2 * section.Slack * Sure::Grow;
    bend             += 4.0 * section.Slack;
    speed            *= order * (1.0 + 0x1p-48);
    bend             *= order * (order - 1.0) * 2.0 * (1.0 + 0x1p-48);
    turn             *= order * (order - 1.0) * (1.0 + 0x1p-48);
    const Run run     = {.Start = lane.Start, .Way = lane.Rough, .Along = lane.Along};
    const Run chord   = {.Start = poles.front(),
                         .Way   = poles.back() - poles.front(),
                         .Along = {.Start = 0.0, .End = 1.0}};
    if (roots->Count > 0 && Dot(chord.Way, chord.Way) > 0.0 && !Unaligned(run, chord, tolerance))
    {
      return std::nullopt;
    }
    Spots spots;
    for (const Root &root : std::span{roots->At}.first(roots->Count))
    {
      if (!Steep(root.Slope / (speedD * speed), tolerance))
      {
        return std::nullopt;
      }
      const double along = Guessed(root, whole.C, degree);
      /**
       *  f at the grid place: from its coefficients where the control points are exact, else from
       *  the true point, which the slack does not blur.
       */
      const std::array<Number, 2> point = Exactly(section, along);
      const Number                value = section.Source == nullptr
                                            ? Valued(equation, degree, along)
                                            : (dx * (point[1] - Sure::Exact(lane.Start.Y)))
                                                  - (dy * (point[0] - Sure::Exact(lane.Start.X)));
      /**
       *  f' at the grid place in doubles, within twice n times each coefficient's bound and the
       *  roundings of the differences and of de Casteljau's rounds.
       */
      Coefficients slopes = {};
      for (std::size_t index = 0; index < degree; ++index)
      {
        slopes.at(index) = order * (whole.C.at(index + 1) - whole.C.at(index));
      }
      const double slope = Steady(Valued(slopes, degree - 1, along));
      const double slack =
          (2.0 * order * (first + (16.0 * Unit * largest))) + std::numeric_limits<double>::min();
      /**
       *  The mean value theorem: the root t* lies within R = (|f(t0)| + E) / g of t0, g the least
       *  |f'| on the stretch, and the step t1 = t0 - f0 / f'(t0) is within [E + 2u |f0| + F (slack
       *  + |f''| R) / g] / |f'(t0)| of it, F the most |f(t0)| may be.
       */
      const double most  = (Sure::Size(value) + value.Bound) * Sure::Grow;
      const double far   = most / root.Slope * Sure::Grow;
      const double step  = (value.Value.High + value.Value.Low) / slope;
      const double off   = ((((value.Bound + (2.0 * Unit * Sure::Size(value)))
                              + (most * (slack + (turn * far)) / root.Slope))
                             / std::abs(slope))
                            + (2.0 * Unit * std::abs(step)))
                         * Sure::Grow;
      const Number place = Number{
          .Value = Twofold::Add(Twofold::Lifted(along), Twofold::Lifted(-step)), .Bound = off};
      /**
       *  The point by Taylor from the grid place: B(t*) = B(t0) + B'(t0) D + B''(xi) D^2 / 2 for D
       *  = t* - t0, which -step is within off of; B'(t0) in doubles off by n 8u times the longest
       *  control step's size per coordinate, and by 2 n times the slack.
       */
      Coefficients wayX = {};
      Coefficients wayY = {};
      for (std::size_t index = 0; index < degree; ++index)
      {
        wayX.at(index) = order * (poles[index + 1].X - poles[index].X);
        wayY.at(index) = order * (poles[index + 1].Y - poles[index].Y);
      }
      const double reachD = std::abs(step) + off;
      const double tail   = ((speed * off)
                             + (speed * 0x1p-48 * std::abs(step))
                             + (2.0 * order * section.Slack * std::abs(step))
                             + (bend * reachD * reachD * 0.5))
                          * Sure::Grow;
      const auto   moved  = [&](const Number &axis, const Coefficients &way) -> Number
      {
        const Twofold::Two lean = Twofold::TwoProduct(Valued(way, degree - 1, along), -step);
        return axis + Number{.Value = {.High = lean.High, .Low = lean.Low}, .Bound = tail};
      };
      const Number x = moved(point[0], wayX);
      const Number y = moved(point[1], wayY);
      /**
       *  The line's place: (X - P).d / |d|^2, its numbers within their bounds.
       */
      const Number line = Sure::Quotient(((x - Sure::Exact(lane.Start.X)) * dx)
                                             + ((y - Sure::Exact(lane.Start.Y)) * dy),
                                         (dx * dx) + (dy * dy));
      /**
       *  A line's end within eight tolerances of the curve lies within R' = T / g times the speed,
       *  and eight tolerances, of a crossing, as |f| <= T only near one: the crossing then is that
       *  far from the line's end along it, or it is not sure.
       */
      const double near = (reach + (speed * cut / root.Slope)) / (lane.Speed * Sure::Shrink);
      const double margin =
          near + (line.Bound * Sure::Grow) + (0x1p-50 * std::abs(line.Value.High));
      const Side where = Within(lane.Along, line.Value.High, margin);
      if (where == Side::Near)
      {
        return std::nullopt;
      }
      if (where == Side::Out)
      {
        continue;
      }
      /**
       *  The root's stretch has its ends clear of the line, so the root lies inside the curve, and
       *  its place rounds inside it: the curve's own parameter there, a spline's span's start and
       *  width times the piece's.
       */
      const auto rounded = Sure::Rounded(std::array{x, y, Led(lane, line), Spread(section, place)});
      if (!rounded)
      {
        return std::nullopt;
      }
      const auto [pointX, pointY, lineT, curve] = *rounded;
      spots.push_back({.At = {.X = pointX, .Y = pointY}, .Along = {lineT, curve}, .Own = true});
    }
    return spots;
  }
} // namespace Geoxy::Swift
