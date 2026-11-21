/***************************************************************************************************
 *  @file      source/aberth.cpp
 *  @brief     The complex roots of a polynomial in Bernstein form, by the Aberth-Ehrlich iteration.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "aberth.hpp"

#include "real.hpp"
#include "trig.hpp"

#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>

namespace Geoxy::Aberth
{
  namespace
  {
    constexpr double Epsilon = std::numeric_limits<double>::epsilon();

    /**
     *  The binomial coefficients of every degree up to Most, each row as the multiplicative formula
     *  gives it.
     */
    constexpr std::array<std::array<double, Most + 1>, Most + 1> Binomials =
        [] -> std::array<std::array<double, Most + 1>, Most + 1>
    {
      std::array<std::array<double, Most + 1>, Most + 1> table = {};
      for (std::size_t degree = 0; degree <= Most; ++degree)
      {
        std::array<double, Most + 1> &row = table.at(degree);
        row.at(0)                         = 1.0;
        for (std::size_t index = 1; index <= degree; ++index)
        {
          row.at(index) = row.at(index - 1)
                        * static_cast<double>(degree - index + 1)
                        / static_cast<double>(index);
        }
      }
      return table;
    }();

    /**
     *  A complex number's product with another.
     *
     *  @param  lhs One number.
     *  @param  rhs The other.
     *  @return     Their product.
     */
    auto Times(Vector lhs, Vector rhs) -> Vector
    {
      return {.X = (lhs.X * rhs.X) - (lhs.Y * rhs.Y), .Y = (lhs.X * rhs.Y) + (lhs.Y * rhs.X)};
    }

    /**
     *  A complex number's reciprocal: its conjugate over its length squared where that is well
     *  within the doubles, and otherwise worked out on the number scaled by the power of two that
     *  brings it near 1, so that its length squared neither overflows nor underflows; 0 for 0, a
     *  step that moves nothing where a step would be without end.
     *
     *  @param  value The number.
     *  @return       Its reciprocal.
     */
    auto Inverse(Vector value) -> Vector
    {
      constexpr double Least  = std::numeric_limits<double>::min();
      constexpr double Safe   = 0x1p-500;
      const double     square = Dot(value, value);
      if (square > Safe && square < 1.0 / Safe)
      {
        return Vector{.X = value.X, .Y = -value.Y} / square;
      }
      const int    power  = std::ilogb(std::max(Largest(value), Least));
      const Vector scaled = {.X = std::ldexp(value.X, -power), .Y = std::ldexp(value.Y, -power)};
      const double size   = std::max(Dot(scaled, scaled), Least);
      return {.X = std::ldexp(scaled.X / size, -power), .Y = std::ldexp(-scaled.Y / size, -power)};
    }

    /**
     *  Newton's step for a polynomial at a point, and whether the polynomial is there within its
     *  rounding.
     */
    struct Step
    {
      Vector Ratio   = {};    /**< p(z) / p'(z).                        */
      bool   Settled = false; /**< Whether p(z) is within its rounding. */
    };

    /**
     *  Newton's step at a point, p and p' in Bernstein form: with b the ratio of the nearer of z
     *  and 1 - z to the other and B the farther, p = B^m S_m and p' = m B^(m-1) S_(m-1), S the sums
     *  of the coefficients, and of their differences, times binomials and powers of b, taken from
     *  the end z is nearer; so p / p' = B S_m / (m S_(m-1)), with the sign of the differences
     *  turned where the sums are taken from 1.
     *
     *  @param  coefficients The Bernstein coefficients.
     *  @param  point        The point.
     *  @return              The step.
     */
    auto Stepped(std::span<const Vector> coefficients, Vector point) -> Step
    {
      const std::size_t degree  = coefficients.size() - 1;
      const Vector      rest    = Vector{.X = 1.0, .Y = 0.0} - point;
      const bool        forward = Length(point) <= Length(rest);
      const Vector      base    = forward ? rest : point;
      const Vector      ratio   = Times(forward ? point : rest, Inverse(base));
      const double      size    = Length(ratio);
      const auto        nth     = [&coefficients, forward, degree](std::size_t index) -> Vector
      {
        return coefficients[forward ? index : degree - index];
      };
      const std::array<double, Most + 1> &row   = Binomials.at(degree);
      const std::array<double, Most + 1> &lower = Binomials.at(degree - 1);
      Vector                              value = nth(degree) * row.at(degree);
      Vector slope = (nth(degree) - nth(degree - 1)) * lower.at(degree - 1);
      double bound = Length(value);
      for (std::size_t index = degree; index-- > 0;)
      {
        value = Times(value, ratio) + (nth(index) * row.at(index));
        bound = (bound * size) + (Length(nth(index)) * row.at(index));
        if (index > 0)
        {
          slope = Times(slope, ratio) + ((nth(index) - nth(index - 1)) * lower.at(index - 1));
        }
      }
      const double sign    = forward ? 1.0 : -1.0;
      const bool   settled = Length(value) <= 4.0 * static_cast<double>(degree) * Epsilon * bound;
      return {.Ratio =
                  Times(Times(base, value), Inverse(slope * (sign * static_cast<double>(degree)))),
              .Settled = settled};
    }

    /**
     *  Where the approximations start: the coefficients in powers of v = 2 t - 1, by the forward
     *  differences' Bernstein forms at the middle, and for each edge of the upper convex hull of
     *  the points (k, log2 |c_k|), as many approximations as the edge is wide, evenly about a
     *  circle whose radius in v is the edge's slope as a power of two, each circle turned on by its
     *  own angle; as many as the first point's k at the middle itself, where the coefficients below
     *  it are 0.
     *
     *  @param coefficients The Bernstein coefficients.
     *  @param roots        Where the approximations go, in t.
     */
    void Started(std::span<const Vector> coefficients, Roots &roots)
    {
      const std::size_t                 degree = coefficients.size() - 1;
      std::array<Vector, Most + 1>      steps{};
      std::array<double, Most + 1>      sizes{};
      std::array<std::size_t, Most + 1> hull{};
      std::size_t                       count = 0;
      std::ranges::copy(coefficients, steps.begin());
      for (std::size_t order = 0; order <= degree; ++order)
      {
        const std::array<double, Most + 1> &row = Binomials.at(degree - order);
        Vector                              sum = {};
        for (std::size_t index = 0; index + order <= degree; ++index)
        {
          sum = sum + (steps.at(index) * row.at(index));
        }
        for (std::size_t index = 0; index + order < degree; ++index)
        {
          steps.at(index) = steps.at(index + 1) - steps.at(index);
        }
        const double size = Length(sum) * Binomials.at(degree).at(order);
        if (!(size > 0.0))
        {
          continue;
        }
        sizes.at(order) = Trig::Log2(size);
        while (count >= 2)
        {
          const std::size_t before = hull.at(count - 2);
          const std::size_t last   = hull.at(count - 1);
          const double      turn =
              (static_cast<double>(last - before) * (sizes.at(order) - sizes.at(before)))
              - ((sizes.at(last) - sizes.at(before)) * static_cast<double>(order - before));
          if (turn < 0.0)
          {
            break;
          }
          --count;
        }
        hull.at(count++) = order;
      }
      for (std::size_t order = 0; order < hull.at(0); ++order)
      {
        roots.Values.at(roots.Count++) = {.X = 0.5, .Y = 0.0};
      }
      for (std::size_t edge = 0; edge + 1 < count; ++edge)
      {
        const std::size_t low    = hull.at(edge);
        const std::size_t high   = hull.at(edge + 1);
        const auto        wide   = static_cast<double>(high - low);
        const double      radius = Trig::Exp2((sizes.at(low) - sizes.at(high)) / wide) / 2.0;
        const double      turned =
            (2.0 * std::numbers::pi * static_cast<double>(low) / static_cast<double>(degree)) + 0.4;
        for (std::size_t step = 0; step < high - low; ++step)
        {
          const double angle = turned + (2.0 * std::numbers::pi * static_cast<double>(step) / wide);
          roots.Values.at(roots.Count++) = {.X = 0.5 + (radius * Trig::Cos(angle)),
                                            .Y = radius * Trig::Sin(angle)};
        }
      }
    }

    /**
     *  A complex number's principal square root.
     *
     *  @param  value The number.
     *  @return       Its root, the real part not negative.
     */
    auto Root(Vector value) -> Vector
    {
      const double size = Length(value);
      return {.X = std::sqrt(std::max(0.0, (size + value.X) / 2.0)),
              .Y = std::copysign(std::sqrt(std::max(0.0, (size - value.X) / 2.0)), value.Y)};
    }

    /**
     *  A root in w = t / (1 - t) as a root in t, w / (1 + w): far off, past every root that
     *  matters, where w is -1 and t without end, the degree in t less than the form's.
     *
     *  @param  ratio The root in w.
     *  @return       The root in t.
     */
    auto Unleaned(Vector ratio) -> Vector
    {
      const Vector sum = Vector{.X = 1.0, .Y = 0.0} + ratio;
      return Zero(Largest(sum)) ? Vector{.X = 0x1p+500, .Y = 0.0} : Times(ratio, Inverse(sum));
    }

    /**
     *  The roots of a polynomial in Bernstein form of degree 1 or 2, in closed form: (1 - t)^m
     *  times b_0 + m b_1 w + b_2 w^2 in w = t / (1 - t), the coefficients as they stand, no
     *  difference of them taken, solved by the stable quadratic formula, the root of the larger
     *  magnitude from -(l + sqrt(l^2 - 4 s c)) / 2 over s, the sign of the root of the discriminant
     *  that of l, so that nothing cancels, and the other from c over the same. A root in w without
     *  end is 1 in t: where b_m is 0, and twice where b_1 is too.
     *
     *  @param coefficients The Bernstein coefficients, two or three, not all 0.
     *  @param roots        Where the roots go.
     */
    void Closed(std::span<const Vector> coefficients, Roots &roots)
    {
      const Vector one      = {.X = 1.0, .Y = 0.0};
      const Vector constant = coefficients.front();
      const Vector square   = coefficients.size() == 3 ? coefficients[2] : Vector{};
      const Vector linear   = coefficients[1] * static_cast<double>(coefficients.size() - 1);
      if (!Zero(Largest(square)))
      {
        Vector root       = Root(Times(linear, linear) - (Times(square, constant) * 4.0));
        root              = Dot(linear, root) < 0.0 ? -root : root;
        const Vector half = (linear + root) * -0.5;
        const Vector far  = Times(half, Inverse(square));
        roots.Values.at(roots.Count++) = Unleaned(far);
        roots.Values.at(roots.Count++) =
            Unleaned(Zero(Largest(half)) ? far : Times(constant, Inverse(half)));
        return;
      }
      if (coefficients.size() == 3)
      {
        roots.Values.at(roots.Count++) = one;
      }
      roots.Values.at(roots.Count++) =
          Zero(Largest(linear)) ? one : Unleaned(-Times(constant, Inverse(linear)));
    }
  } // namespace

  auto Solved(std::span<const Vector> coefficients, int rounds) -> Roots
  {
    Roots roots = {};
    if (coefficients.size() <= 3)
    {
      Closed(coefficients, roots);
      return roots;
    }
    Started(coefficients, roots);
    std::array<bool, Most> settled{};
    std::size_t            left = roots.Count;
    for (int round = 0; round < rounds && left > 0; ++round)
    {
      for (std::size_t index = 0; index < roots.Count; ++index)
      {
        if (settled.at(index))
        {
          continue;
        }
        const Vector point = roots.Values.at(index);
        const Step   step  = Stepped(coefficients, point);
        Vector       sum   = {};
        for (std::size_t other = 0; other < roots.Count; ++other)
        {
          sum = sum + Inverse(point - roots.Values.at(other));
        }
        const Vector move =
            Times(step.Ratio, Inverse(Vector{.X = 1.0, .Y = 0.0} - Times(step.Ratio, sum)));
        if (!step.Settled)
        {
          roots.Values.at(index) = point - move;
        }
        if (step.Settled || !(Length(move) > Epsilon * Length(point)))
        {
          settled.at(index) = true;
          --left;
        }
      }
    }
    return roots;
  }
} // namespace Geoxy::Aberth
