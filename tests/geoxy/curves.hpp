/***************************************************************************************************
 *  @file      tests/geoxy/curves.hpp
 *  @brief     Curves the tests share: one of each kind, and splines with an independent evaluator.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TESTS_CURVES_HPP
#define GEOXY_TESTS_CURVES_HPP

#include "geoxy/bezier/cubic.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"

#include <cmath>
#include <cstddef>

namespace Geoxy::Tests
{
  /**
   *  A cubic Bézier curve, the same curve as Clamped.
   */
  inline constexpr Bezier::Cubic Wave = {.Start = {.X = 0, .Y = 0},
                                         .Out   = {.X = 1, .Y = 3},
                                         .In    = {.X = 3, .Y = -2},
                                         .End   = {.X = 4, .Y = 1}};

  /**
   *  The spline of Wave: its poles, knots clamped at 0 and 1.
   *
   *  @return The spline.
   */
  inline auto Clamped() -> Spline
  {
    return {.Degree  = 3,
            .Knots   = {0, 0, 0, 0, 1, 1, 1, 1},
            .Poles   = {Wave.Start, Wave.Out, Wave.In, Wave.End},
            .Weights = {}};
  }

  /**
   *  A rational quadratic spline that is the quarter of the unit circle from 1,0 to 0,1.
   *
   *  @return The spline.
   */
  inline auto Quarter() -> Spline
  {
    return {.Degree  = 2,
            .Knots   = {0, 0, 0, 1, 1, 1},
            .Poles   = {{.X = 1, .Y = 0}, {.X = 1, .Y = 1}, {.X = 0, .Y = 1}},
            .Weights = {1, std::sqrt(0.5), 1}};
  }

  /**
   *  A cubic spline of several spans, its knots uneven, one doubled, and its weights uneven.
   *
   *  @return The spline.
   */
  inline auto Uneven() -> Spline
  {
    return {
        .Degree  = 3,
        .Knots   = {0, 0, 0, 0, 0.3, 0.3, 1.1, 2, 2, 2, 2},
        .Poles   = {{.X = 0, .Y = 0},
                    {.X = 1, .Y =  2},
                    {.X = 2, .Y = -1},
                    {.X = 4, .Y =  3},
                    {.X = 5, .Y =  0},
                    {.X = 6, .Y =  2},
                    {.X = 8, .Y = 1}},
        .Weights = {1, 2, 0.5, 1, 3, 1, 1},
    };
  }

  /**
   *  Which basis function: its index and its degree.
   */
  struct Function
  {
    std::size_t Index  = 0; /**< Its index.  */
    int         Degree = 0; /**< Its degree. */
  };

  /**
   *  The Cox-de Boor basis function, by its recursion, 0/0 taken as 0: an evaluator that shares
   *  nothing with de Boor's algorithm.
   *
   *  @param  spline   The spline.
   *  @param  function Which basis function.
   *  @param  along    The parameter, inside a span.
   *  @return          Its value.
   */
  inline auto Basis(const Spline &spline, Function function, double along) -> double
  {
    const auto       &knots = spline.Knots;
    const std::size_t index = function.Index;
    if (function.Degree == 0)
    {
      return knots[index] <= along && along < knots[index + 1] ? 1.0 : 0.0;
    }
    const auto     order = static_cast<std::size_t>(function.Degree);
    const Function lower = {.Index = index, .Degree = function.Degree - 1};
    const Function next  = {.Index = index + 1, .Degree = function.Degree - 1};
    const double   left  = knots[index + order] - knots[index];
    const double   right = knots[index + order + 1] - knots[index + 1];
    const double   rise =
        left > 0.0 ? (along - knots[index]) / left * Basis(spline, lower, along) : 0.0;
    const double fall =
        right > 0.0 ? (knots[index + order + 1] - along) / right * Basis(spline, next, along) : 0.0;
    return rise + fall;
  }

  /**
   *  A spline's point by the sum of its weighted basis functions.
   *
   *  @param  spline The spline.
   *  @param  along  The parameter, inside a span.
   *  @return        The point.
   */
  inline auto Summed(const Spline &spline, double along) -> Point
  {
    double x     = 0.0;
    double y     = 0.0;
    double total = 0.0;
    for (std::size_t index = 0; index < spline.Poles.size(); ++index)
    {
      const double weight = spline.Weights.empty() ? 1.0 : spline.Weights[index];
      const double basis = Basis(spline, {.Index = index, .Degree = spline.Degree}, along) * weight;
      x     += basis * spline.Poles[index].X;
      y     += basis * spline.Poles[index].Y;
      total += basis;
    }
    return {.X = x / total, .Y = y / total};
  }
} // namespace Geoxy::Tests

#endif /* GEOXY_TESTS_CURVES_HPP */
