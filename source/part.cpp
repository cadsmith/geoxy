/***************************************************************************************************
 *  @file      source/part.cpp
 *  @brief     A curve's stretch between two places, cut at a place, run on to one, split at many.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/part.hpp"

#include "basis.hpp"
#include "real.hpp"
#include "trig.hpp"
#include "walk.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/reverse.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The stretch of an xline between two places, or of the line a ray runs along.
     *
     *  @param  line    The xline, or a ray's line, the ray's start its point.
     *  @param  stretch Where the part starts and ends.
     *  @return         A line, a ray or an xline, as Part tells it.
     */
    auto Stretched(const Xline &line, Range stretch) -> Curve
    {
      const bool startless = std::isinf(stretch.Start);
      const bool endless   = std::isinf(stretch.End);
      if (!startless && !endless)
      {
        return Line{.Start = At(line, stretch.Start), .End = At(line, stretch.End)};
      }
      const Vector way = stretch.End < stretch.Start ? -line.Way : line.Way;
      if (startless && endless)
      {
        return Xline{.Through = line.Through, .Way = way};
      }
      if (startless)
      {
        return Ray{.Start = At(line, stretch.End), .Way = -way};
      }
      return Ray{.Start = At(line, stretch.Start), .Way = way};
    }

    /**
     *  The blossom of a cubic Bézier curve at three values, by de Casteljau's algorithm, each round
     *  at its own value.
     *
     *  @param  cubic  The curve.
     *  @param  values The values.
     *  @return        The blossom: the curve's point where the three are one value.
     */
    auto Blossomed(const Bezier::Cubic &cubic, std::array<double, 3> values) -> Point
    {
      const Point first  = Lerp(cubic.Start, cubic.Out, values[0]);
      const Point second = Lerp(cubic.Out, cubic.In, values[0]);
      const Point third  = Lerp(cubic.In, cubic.End, values[0]);
      return Lerp(Lerp(first, second, values[1]), Lerp(second, third, values[1]), values[2]);
    }

    /**
     *  The bulge of a part of a polyline's segment.
     *
     *  @param  step The segment, as Locate gives it.
     *  @param  kept The fractions of it the part runs from and to.
     *  @return      The segment's own where all of it is kept; 0 for a line, as Bent tells it; else
     *               tan(atan(b) f) for f the fraction kept, none where that sweeps a whole turn or
     *               more.
     */
    auto Shortened(const Step &step, Range kept) -> std::optional<double>
    {
      if (Exactly(kept.Start, 0.0) && Exactly(kept.End, 1.0))
      {
        return step.Bulge;
      }
      if (Straightness(step))
      {
        return 0.0;
      }
      const double quarter = Trig::Atan(step.Bulge) * (kept.End - kept.Start);
      if (!(std::abs(quarter) < Pi / 2.0))
      {
        return std::nullopt;
      }
      return Trig::Tan(quarter);
    }

    /**
     *  The index of a polyline's segment a place falls in, kept to its segments.
     *
     *  @param  along The place.
     *  @param  last  The index of its last segment.
     *  @param  ended Whether the place ends a part, so that a whole number takes the segment
     *                before.
     *  @return       The index.
     */
    auto Indexed(double along, double last, bool ended) -> std::size_t
    {
      const double index = ended ? std::ceil(along) - 1.0 : std::floor(along);
      return static_cast<std::size_t>(index > 0.0 ? std::min(index, last) : 0.0);
    }
  } // namespace

  auto Part(const Line &line, Range stretch) -> Line
  {
    return {.Start = At(line, stretch.Start), .End = At(line, stretch.End)};
  }

  auto Part(const Ray &ray, Range stretch) -> Curve
  {
    return Stretched({.Through = ray.Start, .Way = ray.Way}, stretch);
  }

  auto Part(const Xline &xline, Range stretch) -> Curve
  {
    return Stretched(xline, stretch);
  }

  auto Part(const Arc &arc, Range stretch) -> Arc
  {
    const Angle start = arc.Start + (arc.Sweep * stretch.Start);
    const Angle end   = arc.Start + (arc.Sweep * stretch.End);
    return {.Center = arc.Center, .Radius = arc.Radius, .Start = start, .Sweep = end - start};
  }

  auto Part(const Ellipse &ellipse, Range stretch) -> Ellipse
  {
    const Angle start = ellipse.Start + (ellipse.Sweep * stretch.Start);
    const Angle end   = ellipse.Start + (ellipse.Sweep * stretch.End);
    return {.Map = ellipse.Map, .Start = start, .Sweep = end - start};
  }

  auto Part(const Polyline &polyline, Range stretch) -> std::optional<Polyline>
  {
    if (stretch.Start > stretch.End)
    {
      return Part(polyline, {.Start = stretch.End, .End = stretch.Start})
          .transform([](const Polyline &part) -> Polyline { return Reversed(part); });
    }
    const std::size_t segments = Count(polyline);
    if (segments == 0)
    {
      return polyline;
    }
    const auto        last  = static_cast<double>(segments - 1);
    const std::size_t first = Indexed(stretch.Start, last, false);
    const std::size_t stop  = std::max(first, Indexed(stretch.End, last, true));
    Polyline          part  = {.Vertices = {}, .Closed = false};
    for (std::size_t index = first; index <= stop; ++index)
    {
      const auto                  whole = static_cast<double>(index);
      const Step                  step  = Locate(polyline, whole);
      const double                from  = index == first ? stretch.Start - whole : 0.0;
      const double                until = index == stop ? stretch.End - whole : 1.0;
      const std::optional<double> bulge = Shortened(step, {.Start = from, .End = until});
      if (!bulge)
      {
        return std::nullopt;
      }
      part.Vertices.push_back(
          {.Position = index == first ? At(polyline, stretch.Start) : step.Start, .Bulge = *bulge});
    }
    part.Vertices.push_back({.Position = At(polyline, stretch.End), .Bulge = 0.0});
    return part;
  }

  auto Part(const Spline &spline, Range stretch) -> Spline
  {
    if (stretch.Start > stretch.End)
    {
      return Reversed(Part(spline, {.Start = stretch.End, .End = stretch.Start}));
    }
    if (!Valid(spline)
        || !std::isfinite(stretch.Start)
        || !std::isfinite(stretch.End)
        || !(stretch.Start < stretch.End))
    {
      return {};
    }
    const auto degree = static_cast<std::size_t>(spline.Degree);
    Spline     part   = {.Degree = spline.Degree, .Knots = {}, .Poles = {}, .Weights = {}};
    part.Knots.assign(degree + 1, stretch.Start);
    for (std::size_t index = degree + 1; index < spline.Poles.size(); ++index)
    {
      const double knot = spline.Knots[index];
      if (stretch.Start < knot && knot < stretch.End)
      {
        part.Knots.push_back(knot);
      }
    }
    part.Knots.insert(part.Knots.end(), degree + 1, stretch.End);

    /**
     *  Each pole is the blossom, at the knots it weighs, of a span of the part it weighs that is
     *  not empty: one of the spline's spans, or a run of one past the domain. The last is the point
     *  at the stretch's end, and is worked on the span At reads there, which is the next where the
     *  end is one of the spline's knots: the same point and weight but for rounding, and At's.
     */
    const Basis                   basis = Spanned(spline);
    const std::size_t             count = part.Knots.size() - degree - 1;
    const std::span<const double> knots = part.Knots;
    for (std::size_t index = 0; index < count; ++index)
    {
      std::size_t span = index;
      while (!(knots[span] < knots[span + 1]))
      {
        ++span;
      }
      const std::size_t over  = Span(basis, index + 1 == count ? stretch.End : knots[span]);
      const int         power = Leveling(spline, over);
      const auto        poles = [&spline, power](std::size_t pole) -> Homogeneous
      {
        return Pole(spline, pole, power);
      };
      const Homogeneous pole =
          Blossom(basis, poles, {.Span = over, .Values = knots.subspan(index + 1, degree)});
      part.Poles.push_back(Projected(pole));
      part.Weights.push_back(std::ldexp(pole.W, -power));
    }
    if (spline.Weights.empty())
    {
      part.Weights.clear();
    }

    /**
     *  Its end poles are At's own points there; a rational spline's are worked on the spans At
     *  reads, as the last pole above is, so that its end weights are the weights At divides by.
     */
    part.Poles.front() = At(spline, stretch.Start);
    part.Poles.back()  = At(spline, stretch.End);
    return part;
  }

  auto Part(const Bezier::Quadratic &quadratic, Range stretch) -> Bezier::Quadratic
  {
    const Point start = Lerp(quadratic.Start, quadratic.Control, stretch.Start);
    const Point end   = Lerp(quadratic.Control, quadratic.End, stretch.Start);
    return {.Start   = At(quadratic, stretch.Start),
            .Control = Lerp(start, end, stretch.End),
            .End     = At(quadratic, stretch.End)};
  }

  auto Part(const Bezier::Cubic &cubic, Range stretch) -> Bezier::Cubic
  {
    const double start = stretch.Start;
    const double end   = stretch.End;
    return {.Start = At(cubic, start),
            .Out   = Blossomed(cubic, {start, start, end}),
            .In    = Blossomed(cubic, {start, end, end}),
            .End   = At(cubic, end)};
  }

  auto Part(const Curve &curve, Range stretch) -> std::optional<Curve>
  {
    return std::visit(
        [stretch](const auto &shape) -> std::optional<Curve>
        {
          if constexpr (std::is_same_v<std::decay_t<decltype(shape)>, Polyline>)
          {
            return Part(shape, stretch)
                .transform([](Polyline part) -> Curve { return Curve{std::move(part)}; });
          }
          else
          {
            return Curve{Part(shape, stretch)};
          }
        },
        curve);
  }

  auto Trim(const Curve &curve, double along, Keep keep) -> std::optional<Curve>
  {
    const Range domain = Domain(curve);
    return Part(curve, keep == Keep::Start ? Range{.Start = domain.Start, .End = along}
                                           : Range{.Start = along, .End = domain.End});
  }

  auto Extend(const Curve &curve, double along) -> std::optional<Curve>
  {
    const Range domain = Domain(curve);
    if (!(along < domain.Start) && !(along > domain.End))
    {
      return curve;
    }
    return Part(curve,
                {.Start = std::min(domain.Start, along), .End = std::max(domain.End, along)});
  }

  auto Split(const Curve &curve, std::span<const double> alongs) -> std::vector<Curve>
  {
    const Range         domain = Domain(curve);
    std::vector<double> cuts   = {domain.Start};
    std::ranges::copy_if(alongs, std::back_inserter(cuts), [domain](double along) -> bool
                         { return domain.Start < along && along < domain.End; });
    Ordered(cuts);
    cuts.push_back(domain.End);
    std::vector<Curve> parts;
    for (std::size_t index = 0; index + 1 < cuts.size(); ++index)
    {
      parts.push_back(Part(curve, {.Start = cuts[index], .End = cuts[index + 1]}).value_or(curve));
    }
    return parts;
  }
} // namespace Geoxy
