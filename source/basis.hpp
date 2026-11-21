/***************************************************************************************************
 *  @file      source/basis.hpp
 *  @brief     De Boor's algorithm: a point of a B-spline from its knots and its control points.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BASIS_HPP
#define GEOXY_BASIS_HPP

#include "real.hpp"

#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <optional>
#include <span>

namespace Geoxy
{
  /**
   *  A point in homogeneous coordinates, the weight held apart: a rational spline's control point
   *  times its weight.
   */
  struct Homogeneous
  {
    double X = 0.0; /**< Along x, times the weight. */
    double Y = 0.0; /**< Along y, times the weight. */
    double W = 0.0; /**< The weight.                */
  };

  /**
   *  A B-spline's basis: its degree, its knots and how many control points they weigh.
   */
  struct Basis
  {
    std::size_t             Degree = 0;  /**< The degree.                                   */
    std::span<const double> Knots  = {}; /**< The knots, a degree and one more than points. */
    std::size_t             Count  = 0;  /**< How many control points.                      */
  };

  /**
   *  A coordinate over a weight of nothing, as a division gives it, with no division, which the
   *  sanitizers stop on: infinite, signed as the coordinate and the weight are, and no number for a
   *  coordinate of nothing or of no number.
   *
   *  @param  value  The coordinate.
   *  @param  weight The weight, zero of either sign.
   *  @return        The quotient.
   */
  inline auto Unbounded(double value, double weight) -> double
  {
    if (Zero(value))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return std::copysign(std::numeric_limits<double>::infinity(), value)
         * std::copysign(1.0, weight);
  }

  /**
   *  The point a homogeneous point stands for: at a weight of nothing, as a run-on spline's may sum
   *  to, a point past every number, as dividing gives it.
   *
   *  @param  point The homogeneous point.
   *  @return       Its coordinates over its weight.
   */
  inline auto Projected(const Homogeneous &point) -> Point
  {
    if (Exactly(point.W, 0.0))
    {
      return {.X = Unbounded(point.X, point.W), .Y = Unbounded(point.Y, point.W)};
    }
    return {.X = point.X / point.W, .Y = point.Y / point.W};
  }

  /**
   *  A spline's control point, homogeneous.
   *
   *  @param  spline The spline.
   *  @param  index  Which control point.
   *  @param  power  A power of two to weigh by as well, as Leveling gives it: the weights about a
   *                 span all times one power of two are the same curve, to the bit.
   *  @return        The point times its weight, and the weight: 1 when the spline has none.
   */
  inline auto Pole(const Spline &spline, std::size_t index, int power = 0) -> Homogeneous
  {
    const Point  pole   = spline.Poles[index];
    const double weight = spline.Weights.empty() ? 1.0 : std::ldexp(spline.Weights[index], power);
    return {.X = pole.X * weight, .Y = pole.Y * weight, .W = weight};
  }

  /**
   *  A spline's basis.
   *
   *  @param  spline The spline.
   *  @return        Its degree, its knots and the count of its control points.
   */
  inline auto Spanned(const Spline &spline) -> Basis
  {
    return {.Degree = static_cast<std::size_t>(spline.Degree),
            .Knots  = spline.Knots,
            .Count  = spline.Poles.size()};
  }

  /**
   *  Whether a spline's counts agree, so that reading it stays within it: a degree from 1 to
   *  MaxDegree, more poles than its degree, a knot more than its poles and its degree, and a weight
   *  per pole or none. What Valid asks of its counts, and nothing of what it asks of its values.
   *
   *  @param  spline The spline.
   *  @return        Whether its counts agree.
   */
  inline auto Counted(const Spline &spline) -> bool
  {
    if (spline.Degree < 1 || spline.Degree > MaxDegree)
    {
      return false;
    }
    const auto        degree = static_cast<std::size_t>(spline.Degree);
    const std::size_t poles  = spline.Poles.size();
    return poles > degree
        && spline.Knots.size() == poles + degree + 1
        && (spline.Weights.empty() || spline.Weights.size() == poles);
  }

  /**
   *  The weighted mix of two homogeneous points.
   *
   *  @param  lhs      The point at 0.
   *  @param  rhs      The point at 1.
   *  @param  fraction How far toward the second.
   *  @return          The mix.
   */
  inline auto Mix(const Homogeneous &lhs, const Homogeneous &rhs, double fraction) -> Homogeneous
  {
    const double rest = 1.0 - fraction;
    return {
        .X = (lhs.X * rest) + (rhs.X * fraction),
        .Y = (lhs.Y * rest) + (rhs.Y * fraction),
        .W = (lhs.W * rest) + (rhs.W * fraction),
    };
  }

  /**
   *  The span a parameter falls in: the last whose start is at most the parameter, within the
   *  domain's spans and not empty. Before a domain whose first spans are empty, a start knot
   *  repeated more than the degree asks, it is the first that is not.
   *
   *  @param  basis The basis.
   *  @param  along The parameter.
   *  @return       The index of the knot that starts the span.
   */
  inline auto Span(const Basis &basis, double along) -> std::size_t
  {
    const std::span<const double> knots = basis.Knots;
    const auto  first = knots.begin() + static_cast<std::ptrdiff_t>(basis.Degree + 1);
    const auto  last  = knots.begin() + static_cast<std::ptrdiff_t>(basis.Count);
    std::size_t span =
        static_cast<std::size_t>(std::distance(knots.begin(), std::upper_bound(first, last, along)))
        - 1;
    while (span > basis.Degree && !(knots[span] < knots[span + 1]))
    {
      --span;
    }
    while (span + 1 < basis.Count && !(knots[span] < knots[span + 1]))
    {
      ++span;
    }
    return span;
  }

  /**
   *  The span a spline is read at for a parameter, where it can be read there: its counts agree,
   *  the span is not empty, the knots de Boor reads around it are finite and none falls, and the
   *  weights it reads are finite and positive.
   *
   *  A valid spline can be read at every parameter. One point of a spline must not cost a pass over
   *  all of it, so nothing farther off is asked: a spline that is not valid elsewhere is still read
   *  where it is sound.
   *
   *  @param  spline The spline, valid or not.
   *  @param  along  The parameter.
   *  @return        The index of the knot that starts the span, or nothing.
   */
  inline auto Legible(const Spline &spline, double along) -> std::optional<std::size_t>
  {
    if (!Counted(spline))
    {
      return std::nullopt;
    }
    const Basis                   basis  = Spanned(spline);
    const std::size_t             degree = basis.Degree;
    const std::size_t             span   = Span(basis, along);
    const std::span<const double> knots  = basis.Knots.subspan(span + 1 - degree, 2 * degree);
    const auto                    sound  = [](double value) -> bool
    {
      return std::isfinite(value) && value > 0.0;
    };
    const bool rising =
        std::ranges::all_of(knots, [](double knot) -> bool { return std::isfinite(knot); })
        && std::ranges::is_sorted(knots)
        && knots[degree - 1] < knots[degree];
    const bool weighted =
        spline.Weights.empty()
        || std::ranges::all_of(std::span{spline.Weights}.subspan(span - degree, degree + 1), sound);
    return rising && weighted ? std::optional{span} : std::nullopt;
  }

  /**
   *  Where a blossom is taken: the span whose control points it reads, and the values it is taken
   *  at, one per round.
   */
  struct Blossoming
  {
    std::size_t             Span   = 0;  /**< The span, as Span gives it.        */
    std::span<const double> Values = {}; /**< The values, as many as the degree. */
  };

  /**
   *  The blossom of a B-spline over a span: de Boor's triangle, each round taken at its own value.
   *
   *  At one value repeated it is the point there; at the span's ends, repeated in turn, it is the
   *  span's Bézier control points.
   *
   *  @pre            The basis is a valid spline's, or its derivative's, and its degree at most
   *                  MaxDegree.
   *  @tparam PointsT What gives the control point at an index.
   *  @param  basis   The basis.
   *  @param  points  The control points, one per index under the basis's count.
   *  @param  where   The span and the values.
   *  @return         The blossom, homogeneous.
   */
  template <typename PointsT>
  auto Blossom(const Basis &basis, const PointsT &points, const Blossoming &where) -> Homogeneous
  {
    const std::size_t degree = basis.Degree;
    const std::size_t span   = where.Span;

    std::array<Homogeneous, MaxDegree + 1> work{};
    const std::span<Homogeneous>           cells{work};
    for (std::size_t index = 0; index <= degree; ++index)
    {
      cells[index] = points(span - degree + index);
    }
    for (std::size_t round = 1; round <= degree; ++round)
    {
      const double along = where.Values[round - 1];
      for (std::size_t index = degree; index >= round; --index)
      {
        const double start = basis.Knots[span - degree + index];
        const double end   = basis.Knots[span + 1 + index - round];
        cells[index]       = Mix(cells[index - 1], cells[index], (along - start) / (end - start));
      }
    }
    return cells[degree];
  }

  /**
   *  The point of a B-spline at a parameter, by de Boor's algorithm: its blossom at the parameter
   *  repeated.
   *
   *  @pre            The basis is a valid spline's, or its derivative's, and its degree at most
   *                  MaxDegree.
   *  @tparam PointsT What gives the control point at an index.
   *  @param  basis   The basis.
   *  @param  points  The control points, one per index under the basis's count.
   *  @param  along   The parameter; before the first span or past the last, that span's polynomial
   *                  run on.
   *  @return         The point, homogeneous.
   */
  template <typename PointsT>
  auto Boor(const Basis &basis, const PointsT &points, double along) -> Homogeneous
  {
    std::array<double, MaxDegree> values{};
    values.fill(along);
    return Blossom(basis, points,
                   {.Span = Span(basis, along), .Values = std::span{values}.first(basis.Degree)});
  }

  /**
   *  Whether a spline's weights differ: a rational spline, where one whose weights are all one
   *  number, or that has none, is a polynomial one.
   *
   *  @param  spline The spline.
   *  @return        Whether it has two weights that are not the same.
   */
  inline auto Rational(const Spline &spline) -> bool
  {
    return std::ranges::adjacent_find(spline.Weights, [](double lhs, double rhs) -> bool
                                      { return lhs < rhs || lhs > rhs; })
        != spline.Weights.end();
  }

  /**
   *  The power of two that brings the weights about a span of a spline near 1, where they are so
   *  heavy or so light, past 2^64 either way, that a pole times its weight, or a weight times
   *  itself, would run past every double or come to nothing, or where the lightest is so light,
   *  under 2^-1021, that it has lost digits and its products come to nothing: 0 for any other
   *  weights, and for none. Weights are told only against each other, so all of them times one
   *  number are the same curve. The heaviest is brought to 1 unless that would leave the lightest
   *  less than a number: beside a weight 2^1100 times it, a weight is still kept, the heaviest left
   *  above 1 by as much.
   *
   *  @param  spline The spline, its weights about the span finite and above nothing, as a valid
   *                 spline's are and as Legible tells them.
   *  @param  span   The index of the knot that starts the span.
   *  @return        The power.
   */
  inline auto Leveling(const Spline &spline, std::size_t span) -> int
  {
    constexpr double Heavy = 0x1p64;
    if (spline.Weights.empty())
    {
      return 0;
    }
    /**
     *  The least power of two a weight is left at: a number still, with all its digits.
     */
    constexpr int Least  = -1021;
    const auto    degree = static_cast<std::size_t>(spline.Degree);
    const auto [lightest, heaviest] =
        std::ranges::minmax(std::span{spline.Weights}.subspan(span - degree, degree + 1));
    if (heaviest <= Heavy && heaviest >= 1.0 / Heavy && std::ilogb(lightest) >= Least)
    {
      return 0;
    }
    return std::max(-std::ilogb(heaviest), Least - std::ilogb(lightest));
  }

  /**
   *  A spline's point as de Boor sums it, its weight set right for a spline that has none.
   *
   *  Such a spline's weight is 1, and what is summed for it is 1 but for rounding: in and about the
   *  domain it is kept as summed, the coordinates being off by as much and the quotient free of it.
   *  Far past the domain the sum is of two numbers far larger than 1 and comes to nothing, or to
   *  twice as much, and the quotient to no number at all: there the weight is 1, as it is.
   *
   *  @param  spline The spline.
   *  @param  point  Its point at a parameter, as summed.
   *  @return        The point, its weight 1 where the spline has none and the sum lost it.
   */
  inline auto Weighted(const Spline &spline, Homogeneous point) -> Homogeneous
  {
    constexpr double Rounding = 0x1p-40;
    if (spline.Weights.empty() && !(std::abs(point.W - 1.0) <= Rounding))
    {
      point.W = 1.0;
    }
    return point;
  }

  /**
   *  The weight a spline's point is divided by at a parameter, its end spans run on past its
   *  domain: where it is nothing, the spline is at no point, and past that its run-on is no curve
   *  to follow.
   *
   *  @param  spline The spline.
   *  @param  along  The parameter.
   *  @return        The weight: 1 for a spline of no weights, 0 where it cannot be read.
   */
  inline auto Weighed(const Spline &spline, double along) -> double
  {
    const std::optional<std::size_t> span = Legible(spline, along);
    if (!span)
    {
      return 0.0;
    }
    const int  power = Leveling(spline, *span);
    const auto poles = [&spline, power](std::size_t index) -> Homogeneous
    {
      return Pole(spline, index, power);
    };
    const Basis                   basis = Spanned(spline);
    std::array<double, MaxDegree> values{};
    const std::span<double>       taken = std::span{values}.first(basis.Degree);
    std::ranges::fill(taken, along);
    return Blossom(basis, poles, {.Span = *span, .Values = taken}).W;
  }
} // namespace Geoxy

#endif /* GEOXY_BASIS_HPP */
