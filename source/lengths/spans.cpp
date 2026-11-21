/***************************************************************************************************
 *  @file      source/lengths/spans.cpp
 *  @brief     The length along a spline or a Bézier curve, span by span, and along a polyline, leg
 *             by leg.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "lengths/spans.hpp"

#include "basis.hpp"
#include "eased.hpp"
#include "exact/two.hpp"
#include "hodograph.hpp"
#include "rational.hpp"
#include "real.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <span>
#include <vector>

namespace Geoxy::Lengths
{
  namespace
  {
    /**
     *  The length of a stretch of a strip, each strip's share of the accuracy its share of the
     *  stretch.
     */
    struct Stage
    {
      Range  Along = {};  /**< The stretch of the strip's own parameter. */
      double Error = 0.0; /**< How far off its length may be.            */
    };

    /**
     *  The length a strip runs over a stretch of its parameter.
     *
     *  @param  strip The strip.
     *  @param  stage The stretch, and how far off its length may be.
     *  @return       The length.
     */
    auto Strode(const Strip &strip, Stage stage) -> double
    {
      return Measured(Paced(strip), stage.Along, {.Within = {.Error = stage.Error}});
    }
  } // namespace

  auto Forever(Vector way) -> double
  {
    return Zero(Dot(way, way)) ? 0.0 : Infinite;
  }

  auto Evenly(Range stretch, double speed) -> double
  {
    return speed * (stretch.End - stretch.Start);
  }

  auto Evenly(double start, double speed, double distance) -> double
  {
    if (Zero(speed))
    {
      return Zero(distance) ? start : std::numeric_limits<double>::quiet_NaN();
    }
    return start + (distance / speed);
  }

  auto Forward(Range stretch) -> Range
  {
    return {.Start = std::min(stretch.Start, stretch.End),
            .End   = std::max(stretch.Start, stretch.End)};
  }

  auto AllFinite(std::initializer_list<Point> points) -> bool
  {
    double sum = 0.0;
    for (const Point point : points)
    {
      sum += (point.X * 0.0) + (point.Y * 0.0);
    }
    return !std::isnan(sum);
  }

  auto Unweighted(const Spline &spline) -> bool
  {
    return spline.Degree <= 3
        && std::ranges::all_of(spline.Weights, [&spline](double weight) -> bool
                               { return Exactly(weight, spline.Weights.front()); });
  }

  auto Hauls(const Spline &spline) -> std::vector<Haul>
  {
    const auto                    degree = static_cast<std::size_t>(spline.Degree);
    const std::span<const double> knots  = spline.Knots;
    const Basis                   slopes = {.Degree = degree - 1,
                                            .Knots  = knots.subspan(1, knots.size() - 2),
                                            .Count  = spline.Poles.size() - 1};
    std::vector<Haul>             legs;
    for (std::size_t span = degree; span < spline.Poles.size(); ++span)
    {
      const Range  stretch = {.Start = knots[span], .End = knots[span + 1]};
      const double width   = stretch.End - stretch.Start;
      if (!(width > 0.0))
      {
        continue;
      }
      const auto steps = [&spline, knots, degree, width](std::size_t index) -> Homogeneous
      {
        const Vector step = (spline.Poles[index + 1] - spline.Poles[index])
                          * (width / (knots[index + degree + 1] - knots[index + 1]));
        return {.X = step.X, .Y = step.Y, .W = 1.0};
      };
      Haul leg = {
          .Shape = {.Degree = degree - 1, .Points = {}, .Factor = static_cast<double>(degree)},
          .Span  = stretch};
      std::array<double, MaxDegree> values{};
      for (std::size_t point = 0; point < degree; ++point)
      {
        for (std::size_t index = 0; index + 1 < degree; ++index)
        {
          values.at(index) = index + 1 + point < degree ? stretch.Start : stretch.End;
        }
        const Homogeneous blossom = Blossom(
            slopes, steps, {.Span = span - 1, .Values = std::span{values}.first(degree - 1)});
        leg.Shape.Points.at(point) = {.X = blossom.X, .Y = blossom.Y};
      }
      legs.push_back(leg);
    }
    return legs;
  }

  auto Traced(const Haul &leg, Range stretch, double error) -> double
  {
    const double width = leg.Span.End - leg.Span.Start;
    const Range  along = {.Start = (stretch.Start - leg.Span.Start) / width,
                          .End   = (stretch.End - leg.Span.Start) / width};
    return Measured(leg.Shape, along, {.Within = {.Error = error}});
  }

  auto Legwise(const std::vector<Haul> &legs, Range stretch, Accuracy accuracy) -> double
  {
    const double width = stretch.End - stretch.Start;
    Exact::Two   total = {};
    for (std::size_t index = 0; index < legs.size(); ++index)
    {
      const Haul  &leg  = legs[index];
      const double from = index == 0 ? stretch.Start : std::max(stretch.Start, leg.Span.Start);
      const double until =
          index + 1 == legs.size() ? stretch.End : std::min(stretch.End, leg.Span.End);
      if (!(from < until))
      {
        continue;
      }
      const double     share = accuracy.Error * ((until - from) / width);
      const Exact::Two sum =
          Exact::TwoSum(total.Big, Traced(leg, {.Start = from, .End = until}, share));
      total = {.Big = sum.Big, .Small = total.Small + sum.Small};
    }
    return total.Big + total.Small;
  }

  auto Spanwise(const Spline &spline, Range stretch, double error) -> double
  {
    const std::vector<Strip> strips = Stripped(spline);
    const Range              domain = Domain(spline);
    const double             width  = stretch.End - stretch.Start;
    const auto               share  = [error, width](Range along, const Slice &where) -> Stage
    {
      return {.Along = along,
              .Error = error
                     * ((along.End - along.Start) * where.Width * (where.To - where.From) / width)};
    };
    std::vector<double> lengths;
    for (const bool after : {false, true})
    {
      const Strip &end  = after ? strips.back() : strips.front();
      const Strip  span = Whole(spline, end.Where.Span);
      const Range  run  = {.Start = after ? std::max(stretch.Start, domain.End) : stretch.Start,
                           .End   = after ? stretch.End : std::min(stretch.End, domain.Start)};
      if (run.Start < run.End)
      {
        lengths.push_back(Strode(
            span, share({.Start = Local(span.Where, run.Start), .End = Local(span.Where, run.End)},
                        span.Where)));
      }
    }
    for (const Strip &strip : strips)
    {
      const Slice &where = strip.Where;
      if (Ahead(where, stretch.End) || Past(where, stretch.Start))
      {
        continue;
      }
      const Range along = {.Start = Ahead(where, stretch.Start) ? 0.0 : Local(where, stretch.Start),
                           .End   = Past(where, stretch.End) ? 1.0 : Local(where, stretch.End)};
      lengths.push_back(Strode(strip, share(along, where)));
    }
    Exact::Two total = {};
    for (const double length : lengths)
    {
      const Exact::Two sum = Exact::TwoSum(total.Big, length);
      total                = {.Big = sum.Big, .Small = total.Small + sum.Small};
    }
    return std::isfinite(total.Big) ? total.Big + total.Small : total.Big;
  }

  auto Signed(Range stretch, double length) -> double
  {
    return stretch.End < stretch.Start ? -length : length;
  }
} // namespace Geoxy::Lengths
