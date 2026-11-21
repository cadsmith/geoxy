/***************************************************************************************************
 *  @file      source/exact/bounds.hpp
 *  @brief     Shewchuk's error bounds, and the test that a sign is beyond them.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_EXACT_BOUNDS_HPP
#define GEOXY_EXACT_BOUNDS_HPP

#include "exact/strict.hpp" // IWYU pragma: keep

/**
 *  The bounds are Shewchuk's, from "Adaptive Precision Floating-Point Arithmetic and Fast Robust
 *  Geometric Predicates" (1997) and its public domain predicates.c, where exactinit computes them
 *  at run time; here they are the same doubles, computed by the compiler. The circle test's second
 *  stage is the one exception: it is compensated, not exact, and its bounds are proved where it is
 *  (exact/circle.cpp). Two more guard the first stages against underflow, Ozaki's u_N and its like
 *  for the circle test.
 */

#include <bit>
#include <cstdint>

namespace Geoxy::Exact
{
  /**
   *  2^-53, half the distance from 1 to the next double.
   */
  inline constexpr double Epsilon = 0x1.0p-53;

  /**
   *  2^27 + 1, which splits a double into two halves of 26 bits.
   */
  inline constexpr double Splitter = 0x1.0p27 + 1.0;

  /**
   *  The relative error of the last approximation of an adaptive predicate.
   */
  inline constexpr double ResultBound = (3.0 + (8.0 * Epsilon)) * Epsilon;

  /**
   *  The relative error of the orientation test's first stage.
   */
  inline constexpr double OrientBoundA = (3.0 + (16.0 * Epsilon)) * Epsilon;

  /**
   *  The smallest normal double, u_N: the orientation test's first stage adds it to its magnitude,
   *  so that no underflow fools it (exact/orient.hpp).
   */
  inline constexpr double Normal = 0x1.0p-1022;

  /**
   *  What underflow can take from the circle test's first stage, per the sum of its lifts plus one:
   *  32 times the smallest double (exact/circle.hpp).
   */
  inline constexpr double CircleBoundUnderflow = 0x1.0p-1069;

  /**
   *  The relative error of the circle test's first stage.
   */
  inline constexpr double CircleBoundA = (10.0 + (96.0 * Epsilon)) * Epsilon;

  /**
   *  The relative error of the circle test's second stage, compensated where Shewchuk's is exact:
   *  its estimate's own error added to his (exact/circle.cpp).
   */
  inline constexpr double CircleBoundB = (4.0 + (128.0 * Epsilon)) * Epsilon;

  /**
   *  The most the circle test's compensated estimate of the determinant of the rounded differences
   *  is off, relative to the permanent and beyond its own roundoff (exact/circle.cpp).
   */
  inline constexpr double CircleBoundCompensated = 64.0 * Epsilon * Epsilon;

  /**
   *  The relative error of the circle test's third stage, Shewchuk's; the compensated estimate it
   *  starts from adds CircleBoundCompensated to it.
   */
  inline constexpr double CircleBoundC = (44.0 + (576.0 * Epsilon)) * Epsilon * Epsilon;

  /**
   *  Where a double's exponent lies in a range of magnitudes from 2^-exponent up to, not including,
   *  2^exponent: read from its bits, a few integer operations and no branch.
   *
   *  @param  value    The double.
   *  @param  exponent The range, below 1023.
   *  @return          0 for zero; below 2 exponent within the range; at least 2 exponent outside
   *                   it, and for NaN, infinity and subnormals.
   */
  constexpr auto Rung(double value, unsigned exponent) -> std::uint64_t
  {
    const std::uint64_t bits = std::bit_cast<std::uint64_t>(value) << 1U;
    const std::uint64_t rung = (bits >> 53U) - (1023U - exponent);
    return bits == 0 ? 0 : rung;
  }

  /**
   *  Whether an approximation of a determinant has its sign beyond doubt.
   *
   *  @param  determinant The approximation.
   *  @param  bound       The most its error can be.
   *  @return             Whether it is at least the bound away from zero; false for NaN.
   */
  constexpr auto Decided(double determinant, double bound) -> bool
  {
    return determinant >= bound || -determinant >= bound;
  }
} // namespace Geoxy::Exact

#endif /* GEOXY_EXACT_BOUNDS_HPP */
