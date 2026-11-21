/***************************************************************************************************
 *  @file      source/lengths/reach.cpp
 *  @brief     The place a distance along a curve: the stretch it falls in bracketed, and the place
 *             in it found.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "lengths/reach.hpp"

#include "arctangent.hpp"
#include "basis.hpp"
#include "eased.hpp"
#include "exact/two.hpp"
#include "kernel/piece.hpp"
#include "lean.hpp"
#include "lengths/pairs.hpp"
#include "lengths/spans.hpp"
#include "rational.hpp"
#include "real.hpp"
#include "twofold.hpp"
#include "walk.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/length.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy::Lengths
{
  namespace
  {
    constexpr double Epsilon = std::numeric_limits<double>::epsilon();

    /**
     *  Where a weight in Bernstein form, followed on past the end its last coefficient stands at,
     *  first comes to nothing.
     *
     *  Past that end the piece's parameter u is 1 / (1 - v) for a v from 0 to 1, and the weight
     *  times (1 - v) to the degree is a polynomial in v whose Bernstein coefficients are the
     *  differences of the weight's own, taken from its last, their signs turn about: its first
     *  root, found as any piece's roots are, is where the weight first is nothing.
     *
     *  @param  weights The piece's weights, in order, the last where the piece is left from.
     *  @return         How far that is from the piece's start, in its own lengths, its end at 1;
     *                  nothing where the weight stays above nothing for ever.
     */
    auto Vanishing(std::span<const double> weights) -> std::optional<double>
    {
      std::array<double, MaxDegree + 1> cells{};
      const std::span<double>           turned = std::span{cells}.first(weights.size());
      std::ranges::reverse_copy(weights, turned.begin());
      for (std::size_t round = 1; round < turned.size(); ++round)
      {
        for (std::size_t index = turned.size() - 1; index >= round; --index)
        {
          turned[index] -= turned[index - 1];
        }
      }
      for (std::size_t index = 1; index < turned.size(); index += 2)
      {
        turned[index] = -turned[index];
      }
      for (const double root : Kernel::Roots(turned, Kernel::Budget))
      {
        if (root < 1.0)
        {
          return 1.0 / (1.0 - root);
        }
      }
      return std::nullopt;
    }

    /**
     *  A stretch of a curve's parameter and the lengths from the curve's start to its ends.
     */
    struct Bracket
    {
      double Low   = 0.0; /**< Where it starts.         */
      double High  = 0.0; /**< Where it ends.           */
      double Below = 0.0; /**< The length to its start. */
      double Above = 0.0; /**< The length to its end.   */
    };

    /**
     *  A stretch of a curve's parameter that a distance from its start falls in, where there is
     *  one: the domain, or the stretch as long before it, doubled past it while the distance lies
     *  beyond, as the curve runs on, for as long as doubles last.
     *
     *  A curve whose length stops growing as it runs on never gets to a distance beyond it: one
     *  that is a point, or a spline with weights run on toward the one point it nears for ever. A
     *  length past where such a spline stops running is without end, as Length tells it, and so
     *  brackets every distance.
     *
     *  @param  gauge    What the curve goes by.
     *  @param  distance The distance, negative before the start.
     *  @return          The stretch; one the distance still lies beyond where the curve, run on,
     *                   never gets that far.
     */
    auto Bracketed(const Gauge &gauge, double distance) -> Bracket
    {
      const auto  &measure = gauge.Measure;
      const Range  domain  = gauge.Domain;
      const double width   = domain.End - domain.Start;
      Bracket      bracket = {
          .Low = domain.Start, .High = domain.End, .Below = 0.0, .Above = measure(domain.End)};
      if (distance < 0.0)
      {
        const double low = domain.Start - width;
        bracket          = {.Low = low, .High = domain.Start, .Below = measure(low), .Above = 0.0};
      }
      for (int doubling = 0; bracket.Above < distance; ++doubling)
      {
        const double high  = bracket.High + (width * std::ldexp(1.0, doubling + 1));
        const double above = std::isfinite(high) ? measure(high) : bracket.Above;
        if (!(above > bracket.Above))
        {
          break;
        }
        bracket = {.Low = bracket.High, .High = high, .Below = bracket.Above, .Above = above};
      }
      for (int doubling = 0; distance < bracket.Below; ++doubling)
      {
        const double low   = bracket.Low - (width * std::ldexp(1.0, doubling + 1));
        const double below = std::isfinite(low) ? measure(low) : bracket.Below;
        if (!(below < bracket.Below))
        {
          break;
        }
        bracket = {.Low = low, .High = bracket.Low, .Below = below, .Above = bracket.Below};
      }
      return bracket;
    }
  } // namespace

  auto Stopped(const Spline &spline, Range stretch) -> bool
  {
    const Range  domain = Domain(spline);
    const double low    = std::min(stretch.Start, stretch.End);
    const double high   = std::max(stretch.Start, stretch.End);
    const bool   before = low < domain.Start;
    const bool   after  = high > domain.End;
    if (spline.Weights.empty() || !(before || after))
    {
      return false;
    }
    std::optional<Kernel::Piece> first;
    Kernel::Piece                last = {};
    Range                        head = {};
    Range                        tail = {};
    Kernel::Sections(spline,
                     [&first, &last, &head, &tail](const Kernel::Piece &piece, Range span) -> void
                     {
                       if (!first)
                       {
                         first = piece;
                         head  = span;
                       }
                       last = piece;
                       tail = span;
                     });
    const auto vanishing = [](const Kernel::Piece &piece, bool backward) -> std::optional<double>
    {
      std::array<double, MaxDegree + 1> weights{};
      for (std::size_t index = 0; index <= piece.Degree; ++index)
      {
        weights.at(index) = piece.Poles.at(backward ? piece.Degree - index : index).W;
      }
      return Vanishing(std::span{weights}.first(piece.Degree + 1));
    };
    bool stopped = false;
    if (after)
    {
      const std::optional<double> stop = vanishing(last, false);
      stopped                          = !(Weighed(spline, high) > 0.0)
                                      || (stop && !(high < tail.Start + (*stop * (tail.End - tail.Start))));
    }
    if (before && !stopped)
    {
      const std::optional<double> stop = vanishing(first.value_or(last), true);
      stopped                          = !(Weighed(spline, low) > 0.0)
                                      || (stop && !(low > head.End - (*stop * (head.End - head.Start))));
    }
    return stopped;
  }

  auto Situated(Point point, double along) -> Place
  {
    if (!std::isfinite(along))
    {
      constexpr double None = std::numeric_limits<double>::quiet_NaN();
      return {.At = {.X = None, .Y = None}, .Along = None};
    }
    return {.At = point, .Along = along};
  }

  auto Coursed(const Spline &spline, Accuracy accuracy) -> Gauge
  {
    std::vector<Haul>   legs   = Hauls(spline);
    const Range         domain = Domain(spline);
    const double        half   = accuracy.Error / 2.0;
    std::vector<double> before;
    Exact::Two          total = {};
    for (const Haul &leg : legs)
    {
      before.push_back(total.Big + total.Small);
      const double share   = half * ((leg.Span.End - leg.Span.Start) / (domain.End - domain.Start));
      const Exact::Two sum = Exact::TwoSum(total.Big, Traced(leg, leg.Span, share));
      total                = {.Big = sum.Big, .Small = total.Small + sum.Small};
    }
    const auto measure = [legs = std::move(legs), before = std::move(before),
                          half](double along) -> double
    {
      const auto after = std::ranges::upper_bound(legs, along, {}, [](const Haul &leg) -> double
                                                  { return leg.Span.Start; });
      const auto index =
          static_cast<std::size_t>(std::max(after - legs.begin(), std::ptrdiff_t{1}) - 1);
      const Haul &leg     = legs[index];
      const Range stretch = {.Start = leg.Span.Start, .End = along};
      return before[index] + Signed(stretch, Traced(leg, Forward(stretch), half));
    };
    return {.Measure = measure,
            .Speed   = [&spline](double along) -> double
            { return Length(Derivative(spline, along)); },
            .Domain = domain};
  }

  auto Charted(const Spline &spline, double error) -> Course
  {
    Course      course = {.Strips = Stripped(spline), .Starts = {}, .Before = {}};
    const Range domain = Domain(spline);
    Exact::Two  total  = {};
    for (const Strip &strip : course.Strips)
    {
      course.Starts.push_back(Homed(strip.Where, 0.0));
      course.Before.push_back(total.Big + total.Small);
      const double     share = error
                             * strip.Where.Width
                             * ((strip.Where.To - strip.Where.From) / (domain.End - domain.Start));
      const Exact::Two sum =
          Exact::TwoSum(total.Big, Measured(Paced(strip), {.Start = 0.0, .End = 1.0},
                                            {.Within = {.Error = share}}));
      total = {.Big = sum.Big, .Small = total.Small + sum.Small};
    }
    course.Before.push_back(total.Big + total.Small);
    return course;
  }

  auto Strewn(const Spline &spline, Course course, double error) -> Gauge
  {
    const Range domain  = Domain(spline);
    const auto  measure = [&spline, course = std::move(course), domain,
                           error](double along) -> double
    {
      const Range whole = {.Start = domain.Start, .End = along};
      if (Stopped(spline, whole))
      {
        return Signed(whole, Infinite);
      }
      const auto after = std::ranges::upper_bound(course.Starts, along);
      const auto index =
          static_cast<std::size_t>(std::max(after - course.Starts.begin(), std::ptrdiff_t{1}) - 1);
      const Strip &strip   = course.Strips[index];
      const Range  stretch = {.Start = 0.0, .End = Local(strip.Where, along)};
      return course.Before[index]
           + Signed(stretch,
                    Measured(Paced(strip), Forward(stretch), {.Within = {.Error = error}}));
    };
    return {.Measure = measure,
            .Speed   = [&spline](double along) -> double
            { return Length(Derivative(spline, along)); },
            .Domain = domain};
  }

  auto Reached(const Gauge &gauge, double distance, double tolerance) -> double
  {
    constexpr int    Steps = 100;
    constexpr double None  = std::numeric_limits<double>::quiet_NaN();
    if (!std::isfinite(distance))
    {
      return None;
    }
    const Bracket bracket = Bracketed(gauge, distance);
    double        low     = bracket.Low;
    double        high    = bracket.High;
    const bool    told    = !std::isnan(bracket.Below) && !std::isnan(bracket.Above);
    if (!told || bracket.Above < distance || distance < bracket.Below)
    {
      return None;
    }
    double along = low + ((high - low) / 2.0);
    if (bracket.Below < bracket.Above && std::isfinite(bracket.Above - bracket.Below))
    {
      along = low + ((high - low) * ((distance - bracket.Below) / (bracket.Above - bracket.Below)));
    }
    for (int step = 0; step < Steps; ++step)
    {
      const double off = gauge.Measure(along) - distance;
      if (!(std::abs(off) > tolerance))
      {
        break;
      }
      if (off < 0.0)
      {
        low = along;
      }
      else
      {
        high = along;
      }
      const double pace = gauge.Speed(along);
      const double next = Zero(pace) ? low : along - (off / pace);
      along             = low < next && next < high ? next : low + ((high - low) / 2.0);
    }
    return along;
  }

  auto Allowed(double distance, double accuracy) -> double
  {
    return std::max(accuracy, 8.0 * Epsilon * std::abs(distance));
  }

  auto Lengthed(const Step &step) -> Double
  {
    const Lean::Way   way    = Lean::Apart(step.End, step.Start);
    const Lean::Sized square = Lean::Dot(way, way);
    if (!(square.Value.High > 0x1p-900) || !(square.Value.High < 0x1p900))
    {
      const Segment segment = Bent(step.Start, step.End, step.Bulge);
      return {.High =
                  std::visit([](const auto &shape) -> double { return Length(shape); }, segment),
              .Low = 0.0};
    }
    const Double chord = Lean::Root(square.Value).Value;
    if (Straightness(step))
    {
      return chord;
    }
    const double       bulge  = step.Bulge;
    const Twofold::Two shrunk = Twofold::TwoProduct(bulge, bulge);
    const Double       spread =
        Lean::Plus({.High = 1.0, .Low = 0.0}, {.High = shrunk.High, .Low = shrunk.Low});
    const Twofold::Two tangent = ArcTangent(bulge);
    const Double       turn    = Lean::Over({.High = tangent.High, .Low = tangent.Low},
                                            {.High = bulge, .Low = 0.0}, 1.0 / bulge);
    return Lean::Times(chord, Lean::Times(spread, turn));
  }
} // namespace Geoxy::Lengths
