/***************************************************************************************************
 *  @file      source/twofold.hpp
 *  @brief     Sums and products of doubles with their rounding errors, and numbers held in two
 *             doubles: double-double arithmetic, about 106 bits.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TWOFOLD_HPP
#define GEOXY_TWOFOLD_HPP

#include "exact/two.hpp"

#include <cmath>

/**
 *  Double-double arithmetic, built on the error-free transformations of exact/two.hpp, the one home
 *  of Knuth's, Dekker's and Shewchuk's sums and products of doubles: here as results of the high
 *  part and the low.
 */
namespace Geoxy::Twofold
{
  /**
   *  A sum or a product, rounded, and what rounding left out: their sum is exact.
   */
  struct Two
  {
    double High = 0.0; /**< The result rounded. */
    double Low  = 0.0; /**< What it left out.   */
  };

  /**
   *  Knuth's two-sum: a + b exactly, as its rounding and the error.
   *
   *  @param  first  a.
   *  @param  second b.
   *  @return        The sum and its error.
   */
  inline auto TwoSum(double first, double second) -> Two
  {
    const Exact::Two sum = Exact::TwoSum(first, second);
    return {.High = sum.Big, .Low = sum.Small};
  }

  /**
   *  Dekker's fast two-sum, where |a| is at least |b|.
   *
   *  @param  big   a.
   *  @param  small b.
   *  @return       The sum and its error.
   */
  inline auto FastTwoSum(double big, double small) -> Two
  {
    const Exact::Two sum = Exact::FastTwoSum(big, small);
    return {.High = sum.Big, .Low = sum.Small};
  }

  /**
   *  a b exactly, its rounding and the error: by a fused multiply-add where the hardware has one,
   *  else by Dekker's split, as exact/two.hpp works it.
   *
   *  @param  first  a.
   *  @param  second b.
   *  @return        The product and its error.
   */
  inline auto TwoProduct(double first, double second) -> Two
  {
    const Exact::Two product = Exact::TwoProduct(first, second);
    return {.High = product.Big, .Low = product.Small};
  }

  /**
   *  A number in two doubles, their sum within about 2^-104 of it: double-double arithmetic. The
   *  leading part is the double nearest the number.
   */
  struct Double
  {
    double High = 0.0; /**< The leading part.             */
    double Low  = 0.0; /**< What the leading part leaves. */
  };

  /**
   *  A double as a double-double.
   *
   *  @param  value The double.
   *  @return       It, nothing left.
   */
  inline auto Lifted(double value) -> Double
  {
    return {.High = value, .Low = 0.0};
  }

  /**
   *  A sum in two doubles, its error within about 2^-104 of the larger term.
   *
   *  @param  first  One.
   *  @param  second The other.
   *  @return        The sum.
   */
  inline auto Add(Double first, Double second) -> Double
  {
    const Two sum = TwoSum(first.High, second.High);
    const Two out = FastTwoSum(sum.High, sum.Low + first.Low + second.Low);
    return {.High = out.High, .Low = out.Low};
  }

  /**
   *  A double-double with its sign changed.
   *
   *  @param  value The number.
   *  @return       Its negative.
   */
  inline auto Negated(Double value) -> Double
  {
    return {.High = -value.High, .Low = -value.Low};
  }

  /**
   *  A difference in two doubles.
   *
   *  @param  first  The number.
   *  @param  second The number taken from it.
   *  @return        The difference.
   */
  inline auto Difference(Double first, Double second) -> Double
  {
    return Add(first, Negated(second));
  }

  /**
   *  A product in two doubles.
   *
   *  @param  first  One.
   *  @param  second The other.
   *  @return        The product.
   */
  inline auto Multiply(Double first, Double second) -> Double
  {
    const Two product = TwoProduct(first.High, second.High);
    const Two out = FastTwoSum(product.High,
                               product.Low + (first.High * second.Low) + (first.Low * second.High));
    return {.High = out.High, .Low = out.Low};
  }

  /**
   *  A quotient in two doubles: one step of long division on the remainder.
   *
   *  @param  first  The dividend.
   *  @param  second The divisor, not zero.
   *  @return        The quotient.
   */
  inline auto Divide(Double first, Double second) -> Double
  {
    const double guess = first.High / second.High;
    const Double rest  = Add(first, Multiply(Lifted(-guess), second));
    const Two    out   = FastTwoSum(guess, rest.High / second.High);
    return {.High = out.High, .Low = out.Low};
  }

  /**
   *  A square root in two doubles: one step of Newton's method on the remainder.
   *
   *  @param  value The number, not negative.
   *  @return       Its square root.
   */
  inline auto Root(Double value) -> Double
  {
    const double guess = std::sqrt(value.High);
    if (std::fpclassify(guess) == FP_ZERO)
    {
      return {};
    }
    const Double rest = Add(value, Multiply(Lifted(-guess), Lifted(guess)));
    const Two    out  = FastTwoSum(guess, rest.High / (2.0 * guess));
    return {.High = out.High, .Low = out.Low};
  }
} // namespace Geoxy::Twofold

#endif /* GEOXY_TWOFOLD_HPP */
