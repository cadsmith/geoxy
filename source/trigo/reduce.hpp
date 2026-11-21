/***************************************************************************************************
 *  @file      source/trigo/reduce.hpp
 *  @brief     An angle less its whole quarter turns: radians as a pair of doubles, degrees exactly.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TRIGO_REDUCE_HPP
#define GEOXY_TRIGO_REDUCE_HPP

#include <bit>
#include <cmath>
#include <cstdint>

namespace Geoxy::Trigo
{
  /**
   *  An angle as whole quarter turns and what is left: the angle is Quarter pi / 2 + High + Low.
   */
  struct Reduced
  {
    double        High    = 0.0; /**< What is left, rounded, at most pi / 4 in size.     */
    double        Low     = 0.0; /**< What rounding left out of it.                      */
    std::uint64_t Quarter = 0;   /**< The quarter turns; only their last two bits count. */
  };

  /**
   *  1.5 2^52: a value of less than 2^51 added to it rounds to a whole number, held in the low
   *  bits.
   */
  inline constexpr double Shifter = 0x1.8p52;

  /**
   *  2 / pi, nearest.
   */
  inline constexpr double TwoOverPi = 0x1.45f306dc9c883p-1;

  /**
   *  pi / 2, nearest: the first of three doubles whose sum is within 2^-160 of it.
   */
  inline constexpr double HalfPiHigh = 0x1.921fb54442d18p+0;

  /**
   *  The second of pi / 2's three doubles: what the first leaves out, nearest.
   */
  inline constexpr double HalfPiMiddle = 0x1.1a62633145c07p-54;

  /**
   *  The third of pi / 2's three doubles: what the first two leave out, nearest.
   */
  inline constexpr double HalfPiLow = -0x1.f1976b7ed8fbcp-110;

  /**
   *  pi / 180, nearest.
   */
  inline constexpr double DegreeHigh = 0x1.1df46a2529d39p-6;

  /**
   *  What DegreeHigh leaves out of pi / 180, nearest.
   */
  inline constexpr double DegreeLow = 0x1.5c1d8becdd291p-62;

  /**
   *  A double's bits.
   *
   *  @param  value The double.
   *  @return       Its bits.
   */
  inline auto BitsOf(double value) -> std::uint64_t
  {
    return std::bit_cast<std::uint64_t>(value);
  }

  /**
   *  The double of some bits.
   *
   *  @param  bits The bits.
   *  @return      Their double.
   */
  inline auto OfBits(std::uint64_t bits) -> double
  {
    return std::bit_cast<double>(bits);
  }

  /**
   *  An angle in radians less its whole quarter turns, in eight operations, its pair within 2^-140
   *  of the exact remainder. k, the nearest whole number to the angle times 2 / pi, comes from the
   *  shifter's low bits. x - k p0 is exact in one fused multiply-add: it lies on the grid of 2^-53
   *  and below 1. Its difference with k p1 is exact too, for the remainder's grid is coarser than
   *  that product's, whichever is larger; and one more fused multiply-add gives what rounding the
   *  product and the difference left out, together.
   *
   *  @param  radians The angle, less than 2^30 in size.
   *  @return         The quarter turns and what is left.
   */
  inline auto Reduce(double radians) -> Reduced
  {
    const double shifted  = std::fma(radians, TwoOverPi, Shifter);
    const double quarters = shifted - Shifter;
    const double first    = std::fma(-quarters, HalfPiHigh, radians);
    const double product  = quarters * HalfPiMiddle;
    const double high     = first - product;
    const double rest     = std::fma(-quarters, HalfPiMiddle, first - high);
    return {.High = high, .Low = std::fma(-quarters, HalfPiLow, rest), .Quarter = BitsOf(shifted)};
  }

  /**
   *  An angle in radians of any size less its whole quarter turns, by Payne and Hanek's reduction:
   *  the angle's 53 bits times the 192 bits of 2 / pi that matter to it, in whole numbers.
   *
   *  @param  magnitude The angle's magnitude, at least 2^30.
   *  @return           The quarter turns and what is left; not a number left for an infinity or a
   *                    not a number.
   */
  auto ReduceFar(double magnitude) -> Reduced;

  /**
   *  An angle in degrees modulo a full turn, exactly: its significand and 2^e each taken modulo
   *  360, the second by 8 times 2^(e - 3) modulo 45, a cycle of twelve.
   *
   *  @param  magnitude The angle's magnitude, finite and at least 2^40.
   *  @return           It modulo 360, in [0, 360), a whole number.
   */
  auto TurnFar(double magnitude) -> double;
} // namespace Geoxy::Trigo

#endif /* GEOXY_TRIGO_REDUCE_HPP */
