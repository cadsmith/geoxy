/***************************************************************************************************
 *  @file      source/corner/surd.cpp
 *  @brief     The signs of numbers with nested square roots, by squaring away one root at a time,
 *             and the double nearest one, by halving the doubles about a guess.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "corner/surd.hpp"

#include "expansion.hpp"

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>

namespace Geoxy::Corner
{
  namespace
  {
    /**
     *  A double's place among the doubles, so that neighbours are neighbouring whole numbers: its
     *  bits as a magnitude, negated below zero; both zeros at 0.
     *
     *  @param  value The double, finite.
     *  @return       Its place.
     */
    auto Rung(double value) -> std::int64_t
    {
      const auto bits = std::bit_cast<std::uint64_t>(value);
      const auto size = static_cast<std::int64_t>(bits & ~(std::uint64_t{1} << 63U));
      return (bits >> 63U) == 0 ? size : -size;
    }

    /**
     *  The double at a place among the doubles.
     *
     *  @param  place The place.
     *  @return       The double; +0 at 0.
     */
    auto Stepped(std::int64_t place) -> double
    {
      if (place >= 0)
      {
        return std::bit_cast<double>(place);
      }
      return -std::bit_cast<double>(-place);
    }

    /**
     *  The point halfway between the double at a place and the next one up, exactly.
     *
     *  @param  place The place, below the largest double's.
     *  @return       The halfway point.
     */
    auto Midway(std::int64_t place) -> Expansion::Number
    {
      const double low  = Stepped(place);
      const double high = Stepped(place + 1);
      return Expansion::Sum(Expansion::Of(low), Expansion::Of((high - low) / 2.0));
    }

    /**
     *  The sign of a + b sqrt(p) + d sqrt(g + h sqrt(p)), exactly: of the first two by squaring
     *  away the root where they differ, of the last by its factor where its root is above nothing,
     *  and where the two parts differ in sign, of the difference of their squares, which holds the
     *  one root alone.
     *
     *  @param  whole a + b sqrt(p) as its three numbers.
     *  @param  outer d.
     *  @param  base  g + h sqrt(p) as its three numbers, p shared.
     *  @return       1, -1 or 0.
     */
    auto Nested(const Form<Big> &whole, const Big &outer, const Form<Big> &base) -> int
    {
      const int first  = Signed(whole.Whole, whole.Inner, whole.Root);
      const int under  = Signed(base.Whole, base.Inner, whole.Root);
      const int second = under > 0 ? Expansion::Sign(outer.Value) : 0;
      if (second == 0)
      {
        return first;
      }
      if (first == 0 || first == second)
      {
        return second;
      }
      /**
       *  (a + b sqrt(p))^2 - d^2 (g + h sqrt(p)) = (a^2 + b^2 p - d^2 g) + (2 a b - d^2 h) sqrt(p).
       */
      const Big squared = outer * outer;
      const Big two     = {.Value = Expansion::Of(2.0)};
      const Big rest    = (whole.Whole * whole.Whole)
                        + (whole.Inner * whole.Inner * whole.Root)
                        - (squared * base.Whole);
      const Big share   = (two * whole.Whole * whole.Inner) - (squared * base.Inner);
      const int larger  = Signed(rest, share, whole.Root);
      return larger * first;
    }
  } // namespace

  auto Signed(const Big &whole, const Big &share, const Big &root) -> int
  {
    const int first  = Expansion::Sign(whole.Value);
    const int second = Expansion::Sign(root.Value) > 0 ? Expansion::Sign(share.Value) : 0;
    if (second == 0)
    {
      return first;
    }
    if (first == 0 || first == second)
    {
      return second;
    }
    const int larger = Expansion::Sign(Expansion::Difference(
        Expansion::Product(whole.Value, whole.Value),
        Expansion::Product(Expansion::Product(share.Value, share.Value), root.Value)));
    return larger * first;
  }

  auto Signed(const Form<Big> &form, const Expansion::Number &mark) -> int
  {
    const Big       less  = {.Value = Expansion::Product(mark, form.Below.Value)};
    const Form<Big> whole = {.Whole = form.Whole - less, .Inner = form.Inner, .Root = form.Root};
    const Form<Big> base  = {.Whole = form.Base, .Inner = form.Bend};
    return Nested(whole, form.Outer, base) * Expansion::Sign(form.Below.Value);
  }

  auto Pinpointed(const Form<Big> &form, double guess, double bound) -> double
  {
    /**
     *  A value of nothing exactly, as a coordinate on an axis is, is told first: the halfway points
     *  about it are below the normal doubles, where a product with them is lost.
     */
    if (Signed(form, Expansion::Number{}) == 0)
    {
      return 0.0;
    }
    const std::int64_t top   = Rung(std::numeric_limits<double>::max());
    const auto         above = [&form](std::int64_t place) -> bool
    {
      return Signed(form, Midway(place)) <= 0;
    };
    std::int64_t low  = -top - 1;
    std::int64_t high = top;
    if (std::abs(guess) + (2.0 * bound) < 0x1p1000)
    {
      const std::int64_t near = Rung(guess - (2.0 * bound)) - 2;
      const std::int64_t far  = Rung(guess + (2.0 * bound)) + 2;
      if (!above(near) && above(far))
      {
        low  = near;
        high = far;
      }
    }
    while (low + 1 < high)
    {
      const std::int64_t middle    = std::midpoint(low, high);
      (above(middle) ? high : low) = middle;
    }
    if (high < top && Signed(form, Midway(high)) == 0)
    {
      const double here = Stepped(high);
      return (std::bit_cast<std::uint64_t>(here) & 1U) == 0 ? here : Stepped(high + 1);
    }
    return Stepped(high);
  }
} // namespace Geoxy::Corner
