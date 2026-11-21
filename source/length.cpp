/***************************************************************************************************
 *  @file      source/length.cpp
 *  @brief     The length of a curve: lines, arcs and polylines exactly, ellipses by Carlson's
 *             elliptic integrals, splines and Bézier curves by Gauss-Legendre rules whose error is
 *             bounded from the roots of their speeds.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/length.hpp"

#include "eased.hpp"
#include "hodograph.hpp"
#include "lean.hpp"
#include "lengths/ellipse.hpp"
#include "lengths/pairs.hpp"
#include "lengths/reach.hpp"
#include "lengths/spans.hpp"
#include "rational.hpp"
#include "real.hpp"
#include "twofold.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
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
#include <cmath>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <limits>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  auto Length(const Line &line) -> double
  {
    return Length(line, Domain(line));
  }

  auto Length(const Ray &ray) -> double
  {
    return Lengths::Forever(ray.Way);
  }

  auto Length(const Xline &xline) -> double
  {
    return Lengths::Forever(xline.Way);
  }

  auto Length(const Arc &arc) -> double
  {
    return Length(arc, Domain(arc));
  }

  auto Length(const Ellipse &ellipse) -> double
  {
    return Length(ellipse, Domain(ellipse));
  }

  auto Length(const Polyline &polyline) -> double
  {
    return Length(polyline, Domain(polyline));
  }

  auto Length(const Spline &spline, Accuracy accuracy) -> double
  {
    return Length(spline, Domain(spline), accuracy);
  }

  auto Length(const Bezier::Quadratic &quadratic, Accuracy accuracy) -> double
  {
    return Length(quadratic, Domain(quadratic), accuracy);
  }

  auto Length(const Bezier::Cubic &cubic, Accuracy accuracy) -> double
  {
    return Length(cubic, Domain(cubic), accuracy);
  }

  auto Length(const Curve &curve, Accuracy accuracy) -> double
  {
    return std::visit(
        [accuracy](const auto &shape) -> double
        {
          if constexpr (requires { Length(shape); })
          {
            return Length(shape);
          }
          else
          {
            return Length(shape, accuracy);
          }
        },
        curve);
  }

  auto Length(const Line &line, Range stretch) -> double
  {
    return Lengths::Evenly(stretch, Length(line.End - line.Start));
  }

  auto Length(const Ray &ray, Range stretch) -> double
  {
    return Lengths::Evenly(stretch, Length(ray.Way));
  }

  auto Length(const Xline &xline, Range stretch) -> double
  {
    return Lengths::Evenly(stretch, Length(xline.Way));
  }

  auto Length(const Arc &arc, Range stretch) -> double
  {
    return Lengths::Evenly(stretch, std::abs(arc.Radius * arc.Sweep.InRadians()));
  }

  auto Length(const Ellipse &ellipse, Range stretch) -> double
  {
    return Lengths::Signed(stretch, Lengths::Arclength(ellipse, stretch));
  }

  auto Length(const Polyline &polyline, Range stretch) -> double
  {
    const std::size_t count = Count(polyline);
    if (count == 0)
    {
      return 0.0;
    }
    if (std::isnan(stretch.Start) || std::isnan(stretch.End))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double low    = std::min(stretch.Start, stretch.End);
    const double high   = std::max(stretch.Start, stretch.End);
    const auto   last   = static_cast<double>(count - 1);
    const auto   first  = static_cast<std::size_t>(std::clamp(std::floor(low), 0.0, last));
    const auto   ending = static_cast<std::size_t>(std::clamp(std::floor(high), 0.0, last));
    /**
     *  The high parts summed by TwoSum, each roundoff and low part summed apart: only the running
     *  sum waits on the one before.
     */
    const std::span<const Vertex> vertices = polyline.Vertices;
    double                        sum      = 0.0;
    double                        rest     = 0.0;
    for (std::size_t index = first; index <= ending; ++index)
    {
      const auto      offset = static_cast<double>(index);
      const double    from   = index == first ? low - offset : 0.0;
      const double    until  = index == ending ? high - offset : 1.0;
      const double    share  = until - from;
      Lengths::Double piece =
          Lengths::Lengthed({.Start = vertices[index].Position,
                             .End = vertices[index + 1 < vertices.size() ? index + 1 : 0].Position,
                             .Bulge    = vertices[index].Bulge,
                             .Fraction = 0.0});
      if (!Exactly(share, 1.0))
      {
        piece = Lean::Times(piece, {.High = share, .Low = 0.0});
      }
      const Twofold::Two added  = Twofold::TwoSum(sum, piece.High);
      sum                       = added.High;
      rest                     += added.Low + piece.Low;
    }
    const double total = sum + rest;
    return Lengths::Signed(stretch, total);
  }

  auto Length(const Spline &spline, Range stretch, Accuracy accuracy) -> double
  {
    if (!Valid(spline))
    {
      return 0.0;
    }
    if (std::isnan(stretch.Start) || std::isnan(stretch.End))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    if (Lengths::Stopped(spline, stretch))
    {
      return Lengths::Signed(stretch, std::numeric_limits<double>::infinity());
    }
    if (Lengths::Unweighted(spline))
    {
      return Lengths::Signed(
          stretch, Lengths::Legwise(Lengths::Hauls(spline), Lengths::Forward(stretch), accuracy));
    }
    return Lengths::Signed(stretch,
                           Lengths::Spanwise(spline, Lengths::Forward(stretch), accuracy.Error));
  }

  auto Length(const Bezier::Quadratic &quadratic, Range stretch, Accuracy accuracy) -> double
  {
    if (std::isnan(stretch.Start)
        || std::isnan(stretch.End)
        || !Lengths::AllFinite({quadratic.Start, quadratic.Control, quadratic.End}))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return Lengths::Signed(
        stretch, Measured(Hodographed(quadratic), Lengths::Forward(stretch), {.Within = accuracy}));
  }

  auto Length(const Bezier::Cubic &cubic, Range stretch, Accuracy accuracy) -> double
  {
    if (std::isnan(stretch.Start)
        || std::isnan(stretch.End)
        || !Lengths::AllFinite({cubic.Start, cubic.Out, cubic.In, cubic.End}))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return Lengths::Signed(
        stretch, Measured(Hodographed(cubic), Lengths::Forward(stretch), {.Within = accuracy}));
  }

  auto Length(const Curve &curve, Range stretch, Accuracy accuracy) -> double
  {
    return std::visit(
        [stretch, accuracy](const auto &shape) -> double
        {
          if constexpr (requires { Length(shape, stretch); })
          {
            return Length(shape, stretch);
          }
          else
          {
            return Length(shape, stretch, accuracy);
          }
        },
        curve);
  }

  auto PlaceAt(const Line &line, double distance) -> Place
  {
    const double along = Lengths::Evenly(0.0, Length(line.End - line.Start), distance);
    return Lengths::Situated(At(line, along), along);
  }

  auto PlaceAt(const Ray &ray, double distance) -> Place
  {
    const double along = Lengths::Evenly(0.0, Length(ray.Way), distance);
    return Lengths::Situated(At(ray, along), along);
  }

  auto PlaceAt(const Xline &xline, double distance) -> Place
  {
    const double along = Lengths::Evenly(0.0, Length(xline.Way), distance);
    return Lengths::Situated(At(xline, along), along);
  }

  auto PlaceAt(const Arc &arc, double distance) -> Place
  {
    const double along = Lengths::Evenly(0.0, Length(arc), distance);
    return Lengths::Situated(At(arc, along), along);
  }

  auto PlaceAt(const Ellipse &ellipse, double distance) -> Place
  {
    const Lengths::Gauge gauge =
        Lengths::Gauged(ellipse, [&ellipse](double along) -> double
                        { return Length(ellipse, {.Start = 0.0, .End = along}); });
    const double along = Lengths::Reached(gauge, distance, Lengths::Allowed(distance, 0.0));
    return Lengths::Situated(At(ellipse, along), along);
  }

  auto PlaceAt(const Polyline &polyline, double distance) -> Place
  {
    const std::size_t count = Count(polyline);
    double            left  = distance;
    double            along =
        count == 0 && !polyline.Vertices.empty() ? Lengths::Evenly(0.0, 0.0, distance) : 0.0;
    for (std::size_t index = 0; index < count; ++index)
    {
      const double segment = std::visit([](const auto &shape) -> double { return Length(shape); },
                                        Nth(polyline, index));
      if (left <= segment || index + 1 == count)
      {
        along = Lengths::Evenly(static_cast<double>(index), segment, left);
        break;
      }
      left -= segment;
    }
    return Lengths::Situated(At(polyline, along), along);
  }

  auto PlaceAt(const Spline &spline, double distance, Accuracy accuracy) -> Place
  {
    if (!Valid(spline))
    {
      return {};
    }
    const Accuracy half = {.Error = accuracy.Error / 2.0};
    if (Lengths::Unweighted(spline))
    {
      const double along = Lengths::Reached(Lengths::Coursed(spline, half), distance,
                                            Lengths::Allowed(distance, half.Error));
      return Lengths::Situated(At(spline, along), along);
    }
    Lengths::Course course = Lengths::Charted(spline, half.Error);
    if (!(distance < 0.0) && distance <= course.Before.back())
    {
      const auto after = std::ranges::upper_bound(course.Before, distance);
      const auto index = static_cast<std::size_t>(
          std::min(after - course.Before.begin(), static_cast<std::ptrdiff_t>(course.Strips.size()))
          - 1);
      const Strip &strip = course.Strips[index];
      const double along = Homed(strip.Where, Walked(Paced(strip), distance - course.Before[index],
                                                     Lengths::Allowed(distance, half.Error)));
      return Lengths::Situated(At(spline, along), along);
    }
    const double along = Lengths::Reached(Lengths::Strewn(spline, std::move(course), half.Error),
                                          distance, Lengths::Allowed(distance, half.Error));
    return Lengths::Situated(At(spline, along), along);
  }

  auto PlaceAt(const Bezier::Quadratic &quadratic, double distance, Accuracy accuracy) -> Place
  {
    const Accuracy       half = {.Error = accuracy.Error / 2.0};
    const Lengths::Gauge gauge =
        Lengths::Gauged(quadratic, [&quadratic, half](double along) -> double
                        { return Length(quadratic, {.Start = 0.0, .End = along}, half); });
    const double along = Lengths::Reached(gauge, distance, Lengths::Allowed(distance, half.Error));
    return Lengths::Situated(At(quadratic, along), along);
  }

  auto PlaceAt(const Bezier::Cubic &cubic, double distance, Accuracy accuracy) -> Place
  {
    const Accuracy       half = {.Error = accuracy.Error / 2.0};
    const Lengths::Gauge gauge =
        Lengths::Gauged(cubic, [&cubic, half](double along) -> double
                        { return Length(cubic, {.Start = 0.0, .End = along}, half); });
    const double along = Lengths::Reached(gauge, distance, Lengths::Allowed(distance, half.Error));
    return Lengths::Situated(At(cubic, along), along);
  }

  auto PlaceAt(const Curve &curve, double distance, Accuracy accuracy) -> Place
  {
    return std::visit(
        [distance, accuracy](const auto &shape) -> Place
        {
          if constexpr (requires { PlaceAt(shape, distance); })
          {
            return PlaceAt(shape, distance);
          }
          else
          {
            return PlaceAt(shape, distance, accuracy);
          }
        },
        curve);
  }
} // namespace Geoxy
