/***************************************************************************************************
 *  @file      source/exact/two.hpp
 *  @brief     Exact sums and products of two doubles, each as a rounded result and its roundoff.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_EXACT_TWO_HPP
#define GEOXY_EXACT_TWO_HPP

#include "exact/strict.hpp" // IWYU pragma: keep

#include <cmath>

/**
 *  Shewchuk's macros of predicates.c as functions, step for step. They are exact only while the
 *  compiler rounds each product and sum on its own: the build turns contraction off, and
 *  tests/geoxy/bits.cpp checks it.
 *
 *  Products are the one change. The roundoff of a product is itself a double, so a fused
 *  multiply-add finds it in one rounding, fma(a, b, -p), where Shewchuk splits both factors into
 *  halves in seventeen operations. The value is the same; only a zero roundoff may come out +0
 *  where his split gives -0, which changes no sign.
 */

/**
 *  Whether the hardware fuses a multiply and an add in one rounding: every arm64, and x86-64 built
 *  for FMA. Elsewhere std::fma is a slow call into the C library, so products split their factors
 *  as Shewchuk does.
 */

#ifndef GEOXY_EXACT_FMA
#  if defined(__aarch64__) || defined(_M_ARM64) || defined(__FMA__)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage): #if reads it, which no constant can serve
#    define GEOXY_EXACT_FMA 1
#  else
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage): #if reads it, which no constant can serve
#    define GEOXY_EXACT_FMA 0
#  endif
#endif

#if !GEOXY_EXACT_FMA
#  include "exact/bounds.hpp"
#endif

namespace Geoxy::Exact
{
  /**
   *  A value as a rounded result and its exact roundoff, which sum to it.
   */
  struct Two
  {
    double Big   = 0.0; /**< The rounded result. */
    double Small = 0.0; /**< What rounding took. */
  };

  /**
   *  Whether a double is zero, of either sign, without comparing doubles for equality.
   *
   *  @param  value The double.
   *  @return       Whether it is +0 or -0; true for NaN too, which no exact test is given.
   */
  constexpr auto IsZero(double value) -> bool
  {
    return !(value < 0.0) && !(value > 0.0);
  }

  /**
   *  The sign of a double.
   *
   *  @param  value The double, not NaN.
   *  @return       1, 0 or -1.
   */
  constexpr auto Sign(double value) -> int
  {
    return static_cast<int>(value > 0.0) - static_cast<int>(value < 0.0);
  }

  /**
   *  The sum of two doubles, exactly.
   *
   *  @param  first  One double.
   *  @param  second The other.
   *  @return        Their sum and its roundoff.
   */
  constexpr auto TwoSum(double first, double second) -> Two
  {
    const double sum           = first + second;
    const double virtualSecond = sum - first;
    const double virtualFirst  = sum - virtualSecond;
    const double secondRound   = second - virtualSecond;
    const double firstRound    = first - virtualFirst;
    return {.Big = sum, .Small = firstRound + secondRound};
  }

  /**
   *  Dekker's fast sum of two doubles, exactly, where the first is at least the second in
   *  magnitude.
   *
   *  @param  big   The larger double.
   *  @param  small The smaller.
   *  @return       Their sum and its roundoff.
   */
  constexpr auto FastTwoSum(double big, double small) -> Two
  {
    const double sum = big + small;
    return {.Big = sum, .Small = small - (sum - big)};
  }

  /**
   *  The roundoff of a difference, given the difference rounded.
   *
   *  @param  first      The double subtracted from.
   *  @param  second     The double subtracted.
   *  @param  difference `first - second`, rounded.
   *  @return            What rounding took.
   */
  constexpr auto TwoDiffTail(double first, double second, double difference) -> double
  {
    const double virtualSecond = first - difference;
    const double virtualFirst  = difference + virtualSecond;
    const double secondRound   = virtualSecond - second;
    const double firstRound    = first - virtualFirst;
    return firstRound + secondRound;
  }

  /**
   *  The difference of two doubles, exactly.
   *
   *  @param  first  The double subtracted from.
   *  @param  second The double subtracted.
   *  @return        Their difference and its roundoff.
   */
  constexpr auto TwoDiff(double first, double second) -> Two
  {
    const double difference = first - second;
    return {.Big = difference, .Small = TwoDiffTail(first, second, difference)};
  }

#if !GEOXY_EXACT_FMA
  /**
   *  A double as two halves of 26 bits each, which sum to it exactly.
   *
   *  @param  value The double.
   *  @return       Its high half and its low half.
   */
  constexpr auto Split(double value) -> Two
  {
    const double scaled = Splitter * value;
    const double big    = scaled - value;
    const double high   = scaled - big;
    return {.Big = high, .Small = value - high};
  }
#endif

  /**
   *  The product of two doubles, exactly.
   *
   *  @param  first  One double.
   *  @param  second The other.
   *  @return        Their product and its roundoff.
   */
  inline auto TwoProduct(double first, double second) -> Two
  {
    const double product = first * second;
#if GEOXY_EXACT_FMA
    return {.Big = product, .Small = std::fma(first, second, -product)};
#else
    const Two    firstHalves  = Split(first);
    const Two    secondHalves = Split(second);
    const double error1       = product - (firstHalves.Big * secondHalves.Big);
    const double error2       = error1 - (firstHalves.Small * secondHalves.Big);
    const double error3       = error2 - (firstHalves.Big * secondHalves.Small);
    return {.Big = product, .Small = (firstHalves.Small * secondHalves.Small) - error3};
#endif
  }

  /**
   *  The square of a double, exactly: the product of it with itself.
   *
   *  Shewchuk's Two_Square merges two of Dekker's terms where products are split, three operations
   *  fewer, and no proof covers the merge; tests/proof/rounding.v proves Dekker's product as
   *  TwoProduct computes it, so the square is that product. Its value is his; only a zero roundoff
   *  may come out +0 where his gives -0.
   *
   *  @param  value The double.
   *  @return       Its square and the roundoff.
   */
  inline auto Square(double value) -> Two
  {
    return TwoProduct(value, value);
  }

} // namespace Geoxy::Exact

#endif /* GEOXY_EXACT_TWO_HPP */
