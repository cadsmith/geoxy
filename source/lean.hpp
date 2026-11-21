/***************************************************************************************************
 *  @file      source/lean.hpp
 *  @brief     Lean double-double for exact kernels: exact differences, cross products that keep
 *             only what can matter, a quotient from one fused multiply-add, and the certificate
 *             that a double is the one nearest a number known within a bound.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_LEAN_HPP
#define GEOXY_LEAN_HPP

#include "real.hpp"
#include "twofold.hpp"

#include "geoxy/point.hpp"

#include <array>
#include <cmath>

/**
 *  Each step's error is bounded here, by hand, for its callers to add up: u is 2^-53, a rounding to
 *  nearest is off by at most u times its result, an error-free transformation is exact, and a
 *  number in two doubles holds its low part within u of its high part, as TwoSum leaves it. A
 *  product too small for the normal doubles is off by at most 2^-1075 more; the callers keep every
 *  term they bound far above that.
 */
namespace Geoxy::Lean
{
  using Twofold::Double;

  /**
   *  The finest grain any lean step is bounded by: four times the worst of the bounds below.
   */
  inline constexpr double Grain = 0x1p-98;

  /**
   *  A factor above 1 by more than a few roundings can take a number down: a product past a bound
   *  times it lies past the bound however its own few roundings went.
   */
  inline constexpr double Loose = 1.0 + 0x1p-48;

  /**
   *  A factor below 1 by more than a few roundings can take a number up: a product short of a bound
   *  times it lies short of the bound however its own few roundings went.
   */
  inline constexpr double Tight = 1.0 - 0x1p-48;

  /**
   *  The least magnitude a coordinate is certified at: below it products may lose bits to
   *  underflow, and an exact zero can never be told by a bound.
   */
  inline constexpr double Tiny = 0x1p-960;

  /**
   *  A displacement exactly, each component a double and what rounding it left.
   */
  using Way = std::array<Double, 2>;

  /**
   *  The displacement from one point to another, exactly: TwoSum on each axis, whose high parts are
   *  the rounded differences plain arithmetic gives.
   *
   *  @param  head The point reached.
   *  @param  tail The point left.
   *  @return      head - tail.
   */
  inline auto Apart(Point head, Point tail) -> Way
  {
    const Twofold::Two x = Twofold::TwoSum(head.X, -tail.X);
    const Twofold::Two y = Twofold::TwoSum(head.Y, -tail.Y);
    return {Double{.High = x.High, .Low = x.Low}, Double{.High = y.High, .Low = y.Low}};
  }

  /**
   *  A number in double-double, and the size of the terms it was summed from, which bounds its
   *  error.
   */
  struct Sized
  {
    Double Value = {};  /**< The number.                       */
    double Size  = 0.0; /**< Its terms' size, from high parts. */
  };

  /**
   *  A cross product u x v of exact displacements in double-double, within 2^-100 P of it, P its
   *  terms' size |u.x v.y| + |u.y v.x| from the high parts.
   *
   *  The high products are worked exactly (TwoProduct) and their difference exactly (TwoSum); the
   *  six terms of a low part with a high one, each at most u(1 + u) times a term of P, sum to at
   *  most 3.01 u P, and are added in nine roundings, within 27.1 u^2 P; the two products of low
   *  parts, at most u^2 P, are dropped; adding the low part to the difference's roundoff rounds
   *  once more, within 4.03 u^2 P. In all 32.2 u^2 P, under 2^-100 P. The last TwoSum leaves the
   *  low part within u of the high one however much the terms cancel. The difference's high part is
   *  the cross product of the high parts as plain arithmetic rounds it.
   *
   *  @param  first  u.
   *  @param  second v.
   *  @return        u x v, and P.
   */
  inline auto Cross(const Way &first, const Way &second) -> Sized
  {
    const Twofold::Two lead  = Twofold::TwoProduct(first[0].High, second[1].High);
    const Twofold::Two trail = Twofold::TwoProduct(first[1].High, second[0].High);
    const Twofold::Two head  = Twofold::TwoSum(lead.High, -trail.High);
    const double       low = ((lead.Low - trail.Low)
                              + ((first[0].High * second[1].Low) + (first[0].Low * second[1].High)))
                           - ((first[1].High * second[0].Low) + (first[1].Low * second[0].High));
    const Twofold::Two out = Twofold::TwoSum(head.High, head.Low + low);
    return {.Value = {.High = out.High, .Low = out.Low},
            .Size  = std::abs(lead.High) + std::abs(trail.High)};
  }

  /**
   *  A dot product u . v of exact displacements in double-double, within 2^-100 P of it, P its
   *  terms' size |u.x v.x| + |u.y v.y| from the high parts: worked as Cross works its cross
   *  product, the same steps with a sum where it takes a difference, within the same bound.
   *
   *  @param  first  u.
   *  @param  second v.
   *  @return        u . v, and P.
   */
  inline auto Dot(const Way &first, const Way &second) -> Sized
  {
    const Twofold::Two lead  = Twofold::TwoProduct(first[0].High, second[0].High);
    const Twofold::Two trail = Twofold::TwoProduct(first[1].High, second[1].High);
    const Twofold::Two head  = Twofold::TwoSum(lead.High, trail.High);
    const double       low = ((lead.Low + trail.Low)
                              + ((first[0].High * second[0].Low) + (first[0].Low * second[0].High)))
                           + ((first[1].High * second[1].Low) + (first[1].Low * second[1].High));
    const Twofold::Two out = Twofold::TwoSum(head.High, head.Low + low);
    return {.Value = {.High = out.High, .Low = out.Low},
            .Size  = std::abs(lead.High) + std::abs(trail.High)};
  }

  /**
   *  A sum in double-double, within 2^-104 (|a| + |b|) of it for a and b the terms' high parts: the
   *  high parts summed exactly, the low parts, at most u (|a| + |b|), in two roundings, and a
   *  TwoSum that leaves the low part within u of the high one however much the terms cancel.
   *
   *  @param  first  One term, its low part within u of its high one.
   *  @param  second The other, likewise.
   *  @return        The sum.
   */
  inline auto Plus(Double first, Double second) -> Double
  {
    const Twofold::Two head = Twofold::TwoSum(first.High, second.High);
    const Twofold::Two out  = Twofold::TwoSum(head.High, head.Low + (first.Low + second.Low));
    return {.High = out.High, .Low = out.Low};
  }

  /**
   *  A product in double-double, within 2^-102 |a b| of it for a and b the factors' high parts: the
   *  high parts' product exact (TwoProduct); the three terms with a low part, at most 3.01 u |a b|,
   *  added in four roundings; the low parts' product, at most u^2 |a b|, dropped; and the high
   *  product, which they cannot reach, summed with them exactly (FastTwoSum).
   *
   *  @param  first  One factor, its low part within u of its high one.
   *  @param  second The other, likewise.
   *  @return        The product.
   */
  inline auto Times(Double first, Double second) -> Double
  {
    return Twofold::Multiply(first, second);
  }

  /**
   *  A double-double with its sign changed, exactly.
   *
   *  @param  value The number.
   *  @return       Its negative.
   */
  inline auto Minus(Double value) -> Double
  {
    return {.High = -value.High, .Low = -value.Low};
  }

  /**
   *  A square root in double-double, and the reciprocal of its high part, which its bound is worked
   *  out with.
   */
  struct Rooted
  {
    Double Value   = {};  /**< The root.                   */
    double Inverse = 0.0; /**< 1 / its high part, rounded. */
  };

  /**
   *  A square root in double-double and the reciprocal of its high part.
   *
   *  The high part s is the root of the number's high part rounded, so the remainder v - s^2 of the
   *  high parts is exact by one fused multiply-add and at most 2.01 u s^2; with the low part, at
   *  most 3.02 u s^2. Half of it over s, by the rounded reciprocal, is within 3.02 u^2 s of the
   *  root's first-order step, and the step after it at most 1.14 u^2 s: in all within 2^-103 s, the
   *  low part at most 1.52 u s and then within u of s, as a last FastTwoSum leaves it.
   *
   *  @param  value The number, its high part above zero, its low part within u of it.
   *  @return       The root, and 1 / s rounded.
   */
  inline auto Root(Double value) -> Rooted
  {
    const double       high    = std::sqrt(value.High);
    const double       inverse = 1.0 / high;
    const double       rest    = std::fma(-high, high, value.High);
    const Twofold::Two out     = Twofold::FastTwoSum(high, (rest + value.Low) * (0.5 * inverse));
    return {.Value = {.High = out.High, .Low = out.Low}, .Inverse = inverse};
  }

  /**
   *  A quotient in double-double from a reciprocal of the divisor's high part: the quotient's high
   *  part, then its remainder by one fused multiply-add, then the remainder's quotient.
   *
   *  The high part q is within 2.01 u of the high parts' quotient, so the remainder of the high
   *  parts is at most 2.02 u |q d|, and rounds within u of that. With the low parts it is at most
   *  4.03 u |q d|, added in three roundings within 10.1 u^2 |q d|, and the reciprocal, its product
   *  and the divisor's low part are each within u: the low part is within 22.2 u^2 |q| of what
   *  makes the quotient exact, under 2^-101 |q|, and at most 4.1 u |q| itself.
   *
   *  @param  top     The dividend, its low part within u of its high one.
   *  @param  bottom  The divisor, likewise, not zero.
   *  @param  inverse 1 / bottom.High, rounded.
   *  @return         top / bottom, within 2^-101 of its high part; not normalized.
   */
  inline auto Over(Double top, Double bottom, double inverse) -> Double
  {
    const double high = top.High * inverse;
    const double rest = std::fma(-high, bottom.High, top.High);
    return {.High = high, .Low = ((rest + top.Low) - (high * bottom.Low)) * inverse};
  }

  /**
   *  A number rounded to a double, with what rounding left, and the size of the terms it is made
   *  of, which its error is bounded by.
   */
  struct Rounded
  {
    double Value = 0.0; /**< The double nearest the sum worked out.     */
    double Rest  = 0.0; /**< What the double leaves of that sum.        */
    double Size  = 0.0; /**< The size of its terms, for bounding error. */
  };

  /**
   *  A start moved along an exact step by a double-double amount: start + step along, as a double
   *  and what it leaves.
   *
   *  The high product (TwoProduct) and its sum with the start (TwoSum) are exact; the product of
   *  the low parts, at most 4.1 u^2 |p| for p the high product, is dropped; the four other low
   *  terms, at most u |s| + 6.2 u |p| for s the exact sum's high part, are added in five roundings;
   *  the double and the rest it leaves round twice more. In all within 50 u^2 (|start| + |p|),
   *  under 2^-100 times Size, besides what the amount's own error moves.
   *
   *  @param  start The start.
   *  @param  step  The step, exact, its low part within u of its high one.
   *  @param  along The amount, its low part within 4.1 u of its high one.
   *  @return       The point's coordinate, its rest, and |start| + |p|.
   */
  inline auto Moved(double start, Double step, Double along) -> Rounded
  {
    const Twofold::Two product = Twofold::TwoProduct(step.High, along.High);
    const Twofold::Two sum     = Twofold::TwoSum(start, product.High);
    const double       low =
        (sum.Low + product.Low) + ((step.High * along.Low) + (step.Low * along.High));
    const double value = sum.High + low;
    return {.Value = value,
            .Rest  = (sum.High - value) + low,
            .Size  = std::abs(start) + std::abs(product.High)};
  }

  /**
   *  Whether a double is the one nearest a number, the hardware's own rounding telling it: with the
   *  number known as value + rest within a bound, the double is the nearest exactly when value +
   *  (rest + bound) and value + (rest - bound) both round back to it.
   *
   *  Let mu be half the gap from value to its neighbour above, a double. If value + (rest + bound)
   *  rounds to value, then rest + bound rounds to at most mu, so is at most mu (1 + u); the number
   *  is at most value + rest + 0.6 bound, below value + mu by 0.4 bound - u mu, which is more than
   *  nothing where the bound is at least 2^-99 |value|, as mu is at most u |value|. So the number
   *  lies below the midpoint above, and likewise above the one below: value is the nearest, and a
   *  number exactly halfway is never certified. The hardware rounds ties to even, so the test holds
   *  at the boundaries as well.
   *
   *  @pre          The number is within 0.6 bound of value + rest, and the bound is at least 2^-99
   *                |value|.
   *  @param  value The double.
   *  @param  rest  What the number leaves of it, roughly.
   *  @param  bound How far value + rest may be from the number.
   *  @return       Whether value is the double nearest the number; never for a value under Tiny, an
   *                exact zero included, nor where the rest or the bound is not a number.
   */
  inline auto Certain(double value, double rest, double bound) -> bool
  {
    return std::abs(value) >= Tiny
        && Exactly(value + (rest + bound), value)
        && Exactly(value + (rest - bound), value);
  }
} // namespace Geoxy::Lean

#endif /* GEOXY_LEAN_HPP */
