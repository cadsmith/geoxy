/***************************************************************************************************
 *  @file      source/fit.cpp
 *  @brief     Cubic splines through points, by a tridiagonal system solved and its residual solved
 *             again; splines framed by a polyline's vertices, and made polylines again.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/fit.hpp"

#include "finite.hpp"
#include "real.hpp"
#include "twofold.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/error.hpp"
#include "geoxy/flatten.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A row of a tridiagonal system: the poles before, on and after the diagonal, times their
     *  coefficients, summing to the right side. In a cyclic system the first row's pole before is
     *  the last one, and the last row's pole after is the first.
     */
    struct Row
    {
      double Before = 0.0; /**< The coefficient of the pole before. */
      double On     = 0.0; /**< The coefficient on the diagonal.    */
      double After  = 0.0; /**< The coefficient of the pole after.  */
      Vector Right  = {};  /**< The right side.                     */
    };

    /**
     *  A tridiagonal system, cyclic or not.
     */
    struct System
    {
      std::vector<Row> Rows   = {};    /**< The rows.                        */
      bool             Cyclic = false; /**< Whether its corners are its own. */
    };

    /**
     *  A tridiagonal system solved by Thomas's algorithm, its corners left out.
     *
     *  @param  rows The rows, at least one.
     *  @return      The unknowns.
     */
    auto Thomas(std::span<const Row> rows) -> std::vector<Vector>
    {
      const std::size_t   count = rows.size();
      std::vector<double> after(count);
      std::vector<Vector> right(count);
      after[0] = rows[0].After / rows[0].On;
      right[0] = rows[0].Right / rows[0].On;
      for (std::size_t index = 1; index < count; ++index)
      {
        const Row   &row   = rows[index];
        const double pivot = row.On - (row.Before * after[index - 1]);
        after[index]       = row.After / pivot;
        right[index]       = (row.Right - (right[index - 1] * row.Before)) / pivot;
      }
      for (std::size_t index = count - 1; index > 0; --index)
      {
        right[index - 1] = right[index - 1] - (right[index] * after[index - 1]);
      }
      return right;
    }

    /**
     *  A cyclic tridiagonal system solved by Sherman and Morrison's formula: the system without its
     *  corners, its diagonal's ends changed to make up for them, solved for the right side and for
     *  the corners' column, and the second taken from the first.
     *
     *  @param  rows The rows, at least three.
     *  @return      The unknowns.
     */
    auto Wraparound(std::span<const Row> rows) -> std::vector<Vector>
    {
      const std::size_t last   = rows.size() - 1;
      const double      corner = rows[last].After;
      const double      top    = rows[0].Before;
      const double      shift  = -rows[0].On;
      std::vector<Row>  plain(rows.begin(), rows.end());
      plain[0].On                      -= shift;
      plain[last].On                   -= corner * top / shift;
      const std::vector<Vector> solved  = Thomas(plain);
      for (Row &row : plain)
      {
        row.Right = {};
      }
      plain[0].Right                   = {.X = shift, .Y = 0.0};
      plain[last].Right                = {.X = corner, .Y = 0.0};
      const std::vector<Vector> column = Thomas(plain);
      const double              below  = 1.0 + column[0].X + (top * column[last].X / shift);
      const Vector              share  = (solved[0] + (solved[last] * (top / shift))) / below;
      std::vector<Vector>       out(rows.size());
      for (std::size_t index = 0; index <= last; ++index)
      {
        out[index] =
            solved[index] - Vector{.X = share.X * column[index].X, .Y = share.Y * column[index].X};
      }
      return out;
    }

    /**
     *  A system solved, cyclic or not.
     *
     *  @param  system The system.
     *  @return        The unknowns.
     */
    auto Resolved(const System &system) -> std::vector<Vector>
    {
      return system.Cyclic ? Wraparound(system.Rows) : Thomas(system.Rows);
    }

    /**
     *  What a row's right side leaves of three products, worked out in twice the precision and
     *  rounded once.
     *
     *  @param  right    The right side.
     *  @param  products The coefficients and the unknowns they multiply.
     *  @return          The residual, rounded.
     */
    auto Leftover(double right, const std::array<std::array<double, 2>, 3> &products) -> double
    {
      double sum   = right;
      double error = 0.0;
      for (const auto &[coefficient, unknown] : products)
      {
        const Twofold::Two product  = Twofold::TwoProduct(coefficient, unknown);
        const Twofold::Two total    = Twofold::TwoSum(sum, -product.High);
        sum                         = total.High;
        error                      += total.Low - product.Low;
      }
      return sum + error;
    }

    /**
     *  A system solved, then solved again for what its answer leaves, every residual in twice the
     *  precision: the answer as near the system's as its coefficients allow.
     *
     *  @param  system The system.
     *  @return        The unknowns.
     */
    auto Refined(System system) -> std::vector<Vector>
    {
      std::vector<Vector> solved = Resolved(system);
      const std::size_t   count  = solved.size();
      for (std::size_t index = 0; index < count; ++index)
      {
        Row         &row    = system.Rows[index];
        const Vector before = solved[(index + count - 1) % count];
        const Vector middle = solved[index];
        const Vector after  = solved[(index + 1) % count];
        row.Right           = {
            .X = Leftover(row.Right.X,
                          {{{row.Before, before.X}, {row.On, middle.X}, {row.After, after.X}}}),
            .Y = Leftover(row.Right.Y,
                          {{{row.Before, before.Y}, {row.On, middle.Y}, {row.After, after.Y}}}),
        };
      }
      const std::vector<Vector> fix = Resolved(system);
      for (std::size_t index = 0; index < count; ++index)
      {
        solved[index] = solved[index] + fix[index];
      }
      return solved;
    }

    /**
     *  The three cubic B-splines not nothing at a knot that starts a span, by Cox and de Boor's
     *  recurrence: those of the span's poles less three, two and one.
     *
     *  @param  knots The knots.
     *  @param  span  The index of the knot, which starts a span that is not empty.
     *  @return       Their values, from the first.
     */
    auto CoxDeBoor(std::span<const double> knots, std::size_t span) -> std::array<double, 3>
    {
      const double            along  = knots[span];
      std::array<double, 4>   held   = {1.0, 0.0, 0.0, 0.0};
      std::array<double, 4>   lefts  = {};
      std::array<double, 4>   rights = {};
      const std::span<double> values{held};
      const std::span<double> left{lefts};
      const std::span<double> right{rights};
      for (std::size_t round = 1; round <= 3; ++round)
      {
        left[round]  = along - knots[span + 1 - round];
        right[round] = knots[span + round] - along;
        double saved = 0.0;
        for (std::size_t index = 0; index < round; ++index)
        {
          const double share = values[index] / (right[index + 1] + left[round - index]);
          values[index]      = saved + (right[index + 1] * share);
          saved              = left[round - index] * share;
        }
        values[round] = saved;
      }
      return {values[0], values[1], values[2]};
    }

    /**
     *  The step from one point's knot to the next's, by the knot rule.
     *
     *  @param  chord The chord between the points.
     *  @param  rule  The rule.
     *  @return       The step.
     */
    auto KnotStep(double chord, Knotting rule) -> double
    {
      if (rule == Knotting::Root)
      {
        return std::sqrt(chord);
      }
      if (rule == Knotting::Uniform)
      {
        return 1.0;
      }
      return chord;
    }

    /**
     *  The points a fit runs through, each but the first less the first, and the chords between
     *  them, the closing one too when it closes.
     */
    struct Fitted
    {
      std::vector<Vector> Points = {}; /**< Each point less the first.     */
      std::vector<double> Chords = {}; /**< Each chord's length, in order. */
    };

    /**
     *  The way a fit leaves or arrives at an end, at the speed of the chord there over its knot
     *  step, a third of a knot step along: how far the pole next to the end is from it.
     *
     *  @param  tangent The tangent asked for.
     *  @param  chord   The chord at that end.
     *  @return         The step from the end to its neighbouring pole; none for a tangent of no
     *                  length.
     */
    auto Leaving(Vector tangent, double chord) -> std::optional<Vector>
    {
      const std::optional<Vector> way = Unit(tangent);
      return way ? std::optional{*way * (chord / 3.0)} : std::nullopt;
    }

    /**
     *  The row that sets an end of an open fit: its tangent's pole, or its second derivative
     *  nothing.
     *
     *  @param  tangent The pole next to the end, less the end, if a tangent was asked for.
     *  @param  steps   The knot step at the end, then the one from that end's knot to the knot two
     *                  points in; and the end point, less the first.
     *  @return         The row, with the end's pole before it and the one two in after it.
     */
    auto Ending(const std::optional<Vector> &tangent, std::array<double, 2> steps) -> Row
    {
      if (tangent)
      {
        return {.Before = 0.0, .On = 1.0, .After = 0.0, .Right = *tangent};
      }
      return {.Before = steps[1], .On = -(steps[0] + steps[1]), .After = steps[0], .Right = {}};
    }

    /**
     *  The open fit's knots and poles from its fitted.
     *
     *  @param  fitted  The points and chords.
     *  @param  knots   The knots, clamped.
     *  @param  fitting The fitting.
     *  @return         The poles less the first point; none for a tangent of no length.
     */
    auto Unclosed(const Fitted &fitted, std::span<const double> knots, const Fitting &fitting)
        -> std::optional<std::vector<Vector>>
    {
      const std::size_t           last = fitted.Points.size() - 1;
      const std::optional<Vector> start =
          fitting.Start ? Leaving(*fitting.Start, fitted.Chords.front()) : std::nullopt;
      const std::optional<Vector> end =
          fitting.End ? Leaving(-*fitting.End, fitted.Chords.back()) : std::nullopt;
      if (fitting.Start.has_value()  != start.has_value()
          || fitting.End.has_value() != end.has_value())
      {
        return std::nullopt;
      }
      System system  = {.Rows = std::vector<Row>(last + 3), .Cyclic = false};
      system.Rows[0] = {.Before = 0.0, .On = 1.0, .After = 0.0, .Right = {}};
      system.Rows[1] = Ending(start, {knots[4] - knots[3], knots[5] - knots[2]});
      for (std::size_t index = 1; index < last; ++index)
      {
        const std::array<double, 3> weights = CoxDeBoor(knots, index + 3);
        system.Rows[index + 1]              = {.Before = weights[0],
                                               .On     = weights[1],
                                               .After  = weights[2],
                                               .Right  = fitted.Points[index]};
      }
      Row ending =
          Ending(end, {knots[last + 3] - knots[last + 2], knots[last + 4] - knots[last + 1]});
      std::swap(ending.Before, ending.After);
      ending.Right          = ending.Right + (end ? fitted.Points[last] : Vector{});
      system.Rows[last + 1] = ending;
      system.Rows[last + 2] = {
          .Before = 0.0, .On = 1.0, .After = 0.0, .Right = fitted.Points[last]};
      return Refined(std::move(system));
    }

    /**
     *  The closed fit's poles from its fitted, the first three repeated after the rest.
     *
     *  @param  fitted The points and chords.
     *  @param  knots  The knots, periodic.
     *  @return        The poles less the first point.
     */
    auto Periodic(const Fitted &fitted, std::span<const double> knots) -> std::vector<Vector>
    {
      const std::size_t count  = fitted.Points.size();
      System            system = {.Rows = std::vector<Row>(count), .Cyclic = true};
      for (std::size_t index = 0; index < count; ++index)
      {
        const std::array<double, 3> weights = CoxDeBoor(knots, index + 3);
        system.Rows[index]                  = {.Before = weights[0],
                                               .On     = weights[1],
                                               .After  = weights[2],
                                               .Right  = fitted.Points[index]};
      }
      const std::vector<Vector> solved = Refined(std::move(system));
      std::vector<Vector>       poles  = {solved.back()};
      poles.insert(poles.end(), solved.begin(), solved.end());
      poles.insert(poles.end(), solved.begin(), solved.begin() + 2);
      return poles;
    }

    /**
     *  The knots of a fit: the running sum of its steps, clamped when open, run on by the loop's
     *  length either side when closed.
     *
     *  @param  steps  The knot steps, one per chord.
     *  @param  closed Whether it closes.
     *  @return        The knots.
     */
    auto Knotted(std::span<const double> steps, bool closed) -> std::vector<double>
    {
      std::vector<double> sums = {0.0};
      for (const double step : steps)
      {
        sums.push_back(sums.back() + step);
      }
      const std::size_t count = steps.size();
      const double      loop  = sums.back();
      if (closed)
      {
        std::vector<double> knots = {sums[count - 3] - loop, sums[count - 2] - loop,
                                     sums[count - 1] - loop};
        knots.insert(knots.end(), sums.begin(), sums.end());
        knots.insert(knots.end(), {sums[1] + loop, sums[2] + loop, sums[3] + loop});
        return knots;
      }
      std::vector<double> knots(3, 0.0);
      knots.insert(knots.end(), sums.begin(), sums.end());
      knots.insert(knots.end(), 3, loop);
      return knots;
    }
  } // namespace

  auto Fit(std::span<const Point> points, const Fitting &fitting) -> std::expected<Spline, Error>
  {
    const bool finite =
        std::ranges::all_of(points, [](Point point) -> bool { return Finite(point); })
        && (!fitting.Start || (std::isfinite(fitting.Start->X) && std::isfinite(fitting.Start->Y)))
        && (!fitting.End || (std::isfinite(fitting.End->X) && std::isfinite(fitting.End->Y)));
    if (!finite)
    {
      return std::unexpected{Error::NotFinite};
    }
    std::span<const Point> through = points;
    if (fitting.Closed && through.size() > 1 && through.front() == through.back())
    {
      through = through.first(through.size() - 1);
    }
    const std::size_t fewest = fitting.Closed ? 3 : 2;
    if (through.size() < fewest)
    {
      return std::unexpected{Error::Degenerate};
    }
    const Point origin = through.front();
    Fitted      fitted = {};
    for (const Point point : through)
    {
      fitted.Points.push_back(point - origin);
    }
    const std::size_t   count = fitting.Closed ? through.size() : through.size() - 1;
    std::vector<double> steps;
    for (std::size_t index = 0; index < count; ++index)
    {
      const double chord = Length(through[(index + 1) % through.size()] - through[index]);
      if (Zero(chord))
      {
        return std::unexpected{Error::Degenerate};
      }
      fitted.Chords.push_back(chord);
      steps.push_back(KnotStep(chord, fitting.Knots));
    }
    std::vector<double> knots = Knotted(steps, fitting.Closed);
    const auto          first = knots.begin() + 3;
    const auto          last  = first + static_cast<std::ptrdiff_t>(count) + 1;
    if (!std::isfinite(*(last - 1)))
    {
      return std::unexpected{Error::NotFinite};
    }
    if (std::adjacent_find(first, last, [](double lhs, double rhs) -> bool { return !(lhs < rhs); })
        != last)
    {
      return std::unexpected{Error::Degenerate};
    }
    std::vector<Vector> poles;
    if (fitting.Closed)
    {
      poles = Periodic(fitted, knots);
    }
    else
    {
      std::optional<std::vector<Vector>> opened = Unclosed(fitted, knots, fitting);
      if (!opened)
      {
        return std::unexpected{Error::Degenerate};
      }
      poles = std::move(*opened);
    }
    Spline spline = {.Degree = 3, .Knots = std::move(knots), .Poles = {}, .Weights = {}};
    for (const Vector pole : poles)
    {
      spline.Poles.push_back(origin + pole);
    }
    if (!fitting.Closed)
    {
      spline.Poles.front() = through.front();
      spline.Poles.back()  = through.back();
    }
    if (!std::ranges::all_of(spline.Poles, [](Point point) -> bool { return Finite(point); }))
    {
      return std::unexpected{Error::NotFinite};
    }
    return spline;
  }

  auto Framed(const Polyline &polyline, int degree) -> std::optional<Spline>
  {
    const std::size_t count = polyline.Vertices.size();
    if (degree < 1 || degree > MaxDegree || std::cmp_less_equal(count, degree))
    {
      return std::nullopt;
    }
    const auto size   = static_cast<std::size_t>(degree);
    Spline     spline = {.Degree = degree, .Knots = {}, .Poles = {}, .Weights = {}};
    for (const Vertex &vertex : polyline.Vertices)
    {
      spline.Poles.push_back(vertex.Position);
    }
    if (polyline.Closed)
    {
      for (std::size_t index = 0; index < size; ++index)
      {
        spline.Poles.push_back(polyline.Vertices[index].Position);
      }
      for (std::size_t index = 0; index < count + (2 * size) + 1; ++index)
      {
        spline.Knots.push_back(static_cast<double>(index) - static_cast<double>(size));
      }
      return spline;
    }
    spline.Knots.assign(size, 0.0);
    for (std::size_t index = 0; index + size < count + 1; ++index)
    {
      spline.Knots.push_back(static_cast<double>(index));
    }
    spline.Knots.insert(spline.Knots.end(), size, static_cast<double>(count - size));
    return spline;
  }

  auto Explode(const Spline &spline, Accuracy accuracy) -> Polyline
  {
    if (!Valid(spline))
    {
      return {};
    }
    std::vector<Point> points   = Flatten(spline, accuracy);
    Polyline           polyline = {};
    polyline.Closed             = points.size() > 2 && points.front() == points.back();
    if (polyline.Closed)
    {
      points.pop_back();
    }
    for (const Point point : points)
    {
      polyline.Vertices.push_back({.Position = point, .Bulge = 0.0});
    }
    return polyline;
  }
} // namespace Geoxy
