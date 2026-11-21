/***************************************************************************************************
 *  @file      source/knots.cpp
 *  @brief     Knots inserted into a spline and its degree raised, by blossoming its spans.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "knots.hpp"

#include "basis.hpp"
#include "real.hpp"

#include "geoxy/domain.hpp"
#include "geoxy/fit.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A new pole: its index, and the new degree.
     */
    struct Seat
    {
      std::size_t Index  = 0; /**< Which pole.     */
      std::size_t Degree = 0; /**< The new degree. */
    };

    /**
     *  The span of a spline whose polynomial a new pole is a blossom of: the first span of the new
     *  knots that is not empty, in the pole's support and the new domain, as the old spline's span
     *  that holds it.
     *
     *  @param  basis The old spline's basis.
     *  @param  knots The new knots.
     *  @param  seat  The new pole's index, and the new degree.
     *  @return       The old span.
     */
    auto Holding(const Basis &basis, std::span<const double> knots, Seat seat) -> std::size_t
    {
      const std::size_t last = knots.size() - seat.Degree - 2;
      std::size_t       span = std::max(seat.Index, seat.Degree);
      const std::size_t end  = std::min(seat.Index + seat.Degree, last);
      while (span < end && !(knots[span] < knots[span + 1]))
      {
        ++span;
      }
      return Span(basis, knots[span]);
    }

    /**
     *  The knots of a spline raised a degree: its domain's ends repeated a degree and two times,
     *  every knot within repeated once more.
     *
     *  @param  spline The spline, valid.
     *  @return        The knots.
     */
    auto Reknotted(const Spline &spline) -> std::vector<double>
    {
      const Range         domain = Domain(spline);
      const auto          ends   = static_cast<std::size_t>(spline.Degree) + 2;
      std::vector<double> knots(ends, domain.Start);
      for (std::size_t index = 0; index < spline.Knots.size(); ++index)
      {
        const double knot = spline.Knots[index];
        const bool   last =
            index + 1 == spline.Knots.size() || !Exactly(spline.Knots[index + 1], knot);
        if (domain.Start < knot && knot < domain.End)
        {
          knots.push_back(knot);
          if (last)
          {
            knots.push_back(knot);
          }
        }
      }
      knots.insert(knots.end(), ends, domain.End);
      return knots;
    }
  } // namespace

  auto Respanned(const Spline &spline, std::span<const double> knots, int raise) -> Spline
  {
    const Basis       basis  = Spanned(spline);
    const auto        raised = static_cast<std::size_t>(raise);
    const std::size_t degree = basis.Degree + raised;
    const std::size_t count  = knots.size() - degree - 1;
    const auto        poles  = [&spline](std::size_t index) -> Homogeneous
    {
      return Pole(spline, index);
    };

    Spline              made = {.Degree  = static_cast<int>(degree),
                                .Knots   = {knots.begin(), knots.end()},
                                .Poles   = std::vector<Point>(count),
                                .Weights = {}};
    std::vector<double> weights(count);
    for (std::size_t index = 0; index < count; ++index)
    {
      const std::size_t span = Holding(basis, knots, {.Index = index, .Degree = degree});
      const std::span<const double> values = knots.subspan(index + 1, degree);
      Homogeneous                   sum    = {};
      for (std::size_t left = 0; left <= raised * basis.Degree; ++left)
      {
        std::array<double, MaxDegree> held{};
        const std::span<double>       kept{held};
        std::size_t                   size = 0;
        for (std::size_t value = 0; value < degree; ++value)
        {
          if (raise == 0 || value != left)
          {
            kept[size] = values[value];
            ++size;
          }
        }
        const Homogeneous blossom =
            Blossom(basis, poles, {.Span = span, .Values = kept.first(basis.Degree)});
        sum = {.X = sum.X + blossom.X, .Y = sum.Y + blossom.Y, .W = sum.W + blossom.W};
      }
      const double share = 1.0 / static_cast<double>((raised * basis.Degree) + 1);
      sum                = {.X = sum.X * share, .Y = sum.Y * share, .W = sum.W * share};
      made.Poles[index]  = spline.Weights.empty() ? Point{.X = sum.X, .Y = sum.Y} : Projected(sum);
      weights[index]     = sum.W;
    }
    if (!spline.Weights.empty())
    {
      made.Weights = std::move(weights);
    }
    return made;
  }

  auto Inserted(const Spline &spline, double knot, int times) -> std::optional<Spline>
  {
    if (!Valid(spline) || times < 0)
    {
      return std::nullopt;
    }
    const Range domain  = Domain(spline);
    const auto  already = static_cast<int>(std::ranges::count_if(
        spline.Knots, [knot](double each) -> bool { return Exactly(each, knot); }));
    if (std::isnan(knot)
        || knot < domain.Start
        || domain.End < knot
        || already + times > spline.Degree)
    {
      return std::nullopt;
    }
    if (times == 0)
    {
      return spline;
    }
    std::vector<double> knots = spline.Knots;
    knots.insert(std::ranges::upper_bound(knots, knot), static_cast<std::size_t>(times), knot);
    return Respanned(spline, knots, 0);
  }

  auto Elevated(const Spline &spline, int times) -> std::optional<Spline>
  {
    if (!Valid(spline) || times < 0 || spline.Degree + times > MaxDegree)
    {
      return std::nullopt;
    }
    Spline raised = spline;
    for (int round = 0; round < times; ++round)
    {
      raised = Respanned(raised, Reknotted(raised), 1);
    }
    return raised;
  }
} // namespace Geoxy
