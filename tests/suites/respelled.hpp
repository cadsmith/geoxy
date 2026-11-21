/***************************************************************************************************
 *  @file      tests/suites/respelled.hpp
 *  @brief     A curve spelled another way, the same geometry: a line as a quadratic Bézier curve,
 *             an arc, an ellipse or a polyline as a spline, a Bézier curve as a spline.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TESTS_SUITES_RESPELLED_HPP
#define GEOXY_TESTS_SUITES_RESPELLED_HPP

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
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy::Suites
{
  /**
   *  A spline of degree 2 being built piece by piece, each piece three poles, the first shared with
   *  the piece before, its inner knots twofold.
   */
  struct Quadratics
  {
    Spline Built = {.Degree = 2, .Knots = {0, 0, 0}, .Poles = {}, .Weights = {}}; /**< So far. */
  };

  /**
   *  One piece added: its middle pole and weight, and its end.
   *
   *  @param spline The spline so far, its start pole set.
   *  @param middle The middle pole and its weight.
   *  @param end    The piece's end.
   */
  inline void Add(Quadratics &spline, std::pair<Point, double> middle, Point end)
  {
    Spline &built = spline.Built;
    if (built.Poles.size() > 1)
    {
      built.Knots.push_back(built.Knots.back());
    }
    built.Poles.push_back(middle.first);
    built.Weights.push_back(middle.second);
    built.Poles.push_back(end);
    built.Weights.push_back(1.0);
    built.Knots.push_back(built.Knots.back() + 1.0);
  }

  /**
   *  The spline closed off: its last knot made threefold, and its last pole made its first where
   *  they are one but for rounding, as a whole turn's sine and cosine leave them.
   *
   *  @param  spline The spline.
   *  @return        It, valid.
   */
  inline auto Closed(Quadratics spline) -> Spline
  {
    std::vector<Point> &poles = spline.Built.Poles;
    if (poles.size()                             > 2
        && Distance(poles.front(), poles.back()) <= 1e-12 * Distance(poles.front(), poles[1]))
    {
      poles.back() = poles.front();
    }
    spline.Built.Knots.push_back(spline.Built.Knots.back());
    spline.Built.Knots.push_back(spline.Built.Knots.back());
    return spline.Built;
  }

  /**
   *  An elliptic arc added as rational quadratic pieces of at most a quarter turn: each from the
   *  map of the unit circle at its ends, its middle pole where their tangents meet, weighted by the
   *  cosine of half its turn.
   *
   *  @param spline The spline so far, its last pole the arc's start.
   *  @param map    The map of the unit circle.
   *  @param angles The arc's start and sweep, in radians.
   */
  inline void Add(Quadratics &spline, const Matrix &map, std::pair<double, double> angles)
  {
    const auto [start, sweep] = angles;
    const auto count =
        static_cast<std::size_t>(std::ceil(std::abs(sweep) / (std::numbers::pi / 2)));
    const double step = sweep / static_cast<double>(std::max<std::size_t>(count, 1));
    for (std::size_t piece = 0; piece < std::max<std::size_t>(count, 1); ++piece)
    {
      const double from   = start + (step * static_cast<double>(piece));
      const double half   = step / 2;
      const double middle = from + half;
      const Point  apex   = {.X = std::cos(middle) / std::cos(half),
                             .Y = std::sin(middle) / std::cos(half)};
      const Point  end    = {.X = std::cos(from + step), .Y = std::sin(from + step)};
      Add(spline, {map * apex, std::cos(half)}, map * end);
    }
  }

  /**
   *  A curve spelled another way, the same geometry within rounding, or nothing where it has no
   *  other spelling here: a ray, an endless line, a spline.
   *
   *  @param  curve The curve.
   *  @return       It spelled another way.
   */
  inline auto Respelled(const Curve &curve) -> std::optional<Curve>
  {
    if (const auto *line = std::get_if<Line>(&curve))
    {
      return Bezier::Quadratic{
          .Start = line->Start, .Control = Midpoint(line->Start, line->End), .End = line->End};
    }
    if (const auto *quadratic = std::get_if<Bezier::Quadratic>(&curve))
    {
      return Spline{.Degree  = 2,
                    .Knots   = {0, 0, 0, 1, 1, 1},
                    .Poles   = {quadratic->Start, quadratic->Control, quadratic->End},
                    .Weights = {}};
    }
    if (const auto *cubic = std::get_if<Bezier::Cubic>(&curve))
    {
      return Spline{.Degree  = 3,
                    .Knots   = {0, 0, 0, 0, 1, 1, 1, 1},
                    .Poles   = {cubic->Start, cubic->Out, cubic->In, cubic->End},
                    .Weights = {}};
    }
    if (const auto *arc = std::get_if<Arc>(&curve))
    {
      const Matrix map = {
          .X = {.X = arc->Radius, .Y = 0}, .Y = {.X = 0, .Y = arc->Radius}, .Origin = arc->Center};
      Quadratics spline;
      spline.Built.Poles = {
          map
          * Point{.X = std::cos(arc->Start.InRadians()), .Y = std::sin(arc->Start.InRadians())}};
      spline.Built.Weights = {1.0};
      Add(spline, map, {arc->Start.InRadians(), arc->Sweep.InRadians()});
      return Closed(spline);
    }
    if (const auto *ellipse = std::get_if<Ellipse>(&curve))
    {
      Quadratics spline;
      spline.Built.Poles   = {ellipse->Map
                              * Point{.X = std::cos(ellipse->Start.InRadians()),
                                      .Y = std::sin(ellipse->Start.InRadians())}};
      spline.Built.Weights = {1.0};
      Add(spline, ellipse->Map, {ellipse->Start.InRadians(), ellipse->Sweep.InRadians()});
      return Closed(spline);
    }
    if (const auto *polyline = std::get_if<Polyline>(&curve))
    {
      const std::size_t count = Count(*polyline);
      if (count == 0)
      {
        return std::nullopt;
      }
      Quadratics spline;
      for (std::size_t index = 0; index < count; ++index)
      {
        const Segment segment = Nth(*polyline, index);
        if (const auto *straight = std::get_if<Line>(&segment))
        {
          if (index == 0)
          {
            spline.Built.Poles   = {straight->Start};
            spline.Built.Weights = {1.0};
          }
          Add(spline, {Midpoint(straight->Start, straight->End), 1.0}, straight->End);
          continue;
        }
        const Arc &bow = std::get<Arc>(segment);
        if (index == 0)
        {
          spline.Built.Poles   = {bow.Center + (Direction(bow.Start) * bow.Radius)};
          spline.Built.Weights = {1.0};
        }
        Add(spline,
            Matrix{.X      = {.X = bow.Radius, .Y = 0},
                   .Y      = {.X = 0, .Y = bow.Radius},
                   .Origin = bow.Center},
            {bow.Start.InRadians(), bow.Sweep.InRadians()});
      }
      return Closed(spline);
    }
    return std::nullopt;
  }
} // namespace Geoxy::Suites

#endif /* GEOXY_TESTS_SUITES_RESPELLED_HPP */
