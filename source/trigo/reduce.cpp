/***************************************************************************************************
 *  @file      source/trigo/reduce.cpp
 *  @brief     Angles too large for the quick reductions: radians by Payne and Hanek, degrees by
 *             cycles of powers of two.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "trigo/reduce.hpp"

#include "exact/grid.hpp"
#include "radians.hpp"
#include "real.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace Geoxy::Trigo
{
  namespace
  {
    /**
     *  2 / pi's first 1280 bits after the point, 64 to a word, after a word of zeros: the words a
     *  significand's product reads from start where its bits do.
     */
    constexpr std::array<std::uint64_t, 21> FarWords = [] -> std::array<std::uint64_t, 21>
    {
      std::array<std::uint64_t, 21> words = {};
      for (std::size_t word = 1; word < words.size(); ++word)
      {
        words.at(word) = (std::uint64_t{Table::TwoOverPi.at((2 * word) - 2)} << 32U)
                       | Table::TwoOverPi.at((2 * word) - 1);
      }
      return words;
    }();

    /**
     *  2^n modulo 45, n from 0 to 11: the cycle 2 makes modulo 45.
     */
    constexpr std::array<std::uint64_t, 12> FortyFifths = {1,  2,  4,  8,  16, 32,
                                                           19, 38, 31, 17, 34, 23};

    /**
     *  The 64 bits of 2 / pi that start a given number of bits after the point.
     *
     *  @param  start Where they start, from 0 to 1215.
     *  @return       Them.
     */
    auto FarWindow(int start) -> std::uint64_t
    {
      const auto word  = static_cast<std::size_t>(start) >> 6U;
      const auto shift = static_cast<unsigned>(start) & 63U;
      if (shift == 0)
      {
        return FarWords.at(word);
      }
      return (FarWords.at(word) << shift) | (FarWords.at(word + 1) >> (64U - shift));
    }
  } // namespace

  auto ReduceFar(double magnitude) -> Reduced
  {
    if (!(magnitude < Infinite))
    {
      return {.High = std::numeric_limits<double>::quiet_NaN(), .Low = 0.0, .Quarter = 0};
    }
    constexpr std::uint64_t Fraction = (std::uint64_t{1} << 52U) - 1U;

    const std::uint64_t bits        = BitsOf(magnitude);
    const int           exponent    = static_cast<int>(bits >> 52U) - 1075;
    const std::uint64_t significand = (bits & Fraction) | (std::uint64_t{1} << 52U);

    /**
     *  The significand times 192 bits of 2 / pi, from 2 bits above the point of the product on: two
     *  bits of quarter turns, then 128 of the fraction of one, then guard bits.
     */
    const int            start = 62 + exponent;
    const Exact::Product low   = Exact::Multiply(significand, FarWindow(start + 128));
    const Exact::Product mid   = Exact::Multiply(significand, FarWindow(start + 64));
    const Exact::Product high  = Exact::Multiply(significand, FarWindow(start));

    const std::uint64_t carried = low.High + mid.Low;
    const std::uint64_t third   = low.Low;
    const std::uint64_t second  = carried;
    const std::uint64_t first   = (carried < low.High ? 1U : 0U) + mid.High + high.Low;

    /**
     *  The fraction rounded to the nearest quarter turn: past a half, the next turn, less.
     */
    std::uint64_t quarter   = first >> 62U;
    std::uint64_t upper     = (first << 2U) | (second >> 62U);
    std::uint64_t lower     = (second << 2U) | (third >> 62U);
    const bool    negative  = (upper >> 63U) != 0;
    quarter                += upper >> 63U;
    if (negative)
    {
      lower = ~lower + 1U;
      upper = ~upper + static_cast<std::uint64_t>(lower == 0);
    }
    /**
     *  No double lies within 2^-62 of a quarter turn's multiple, so the fraction, at most a half,
     *  has its leading one from the second bit to the 62nd: shifted up to the top.
     */
    const int lead = std::clamp(std::countl_zero(upper), 1, 63);
    upper   = (upper << static_cast<unsigned>(lead)) | (lower >> static_cast<unsigned>(64 - lead));
    lower <<= static_cast<unsigned>(lead);

    /**
     *  The fraction as a pair, times pi / 2 as a pair.
     */
    const double scale = OfBits(static_cast<std::uint64_t>(1023 - 64 - lead) << 52U);
    const double head  = static_cast<double>(upper & ~std::uint64_t{0x7FF}) * scale;
    const double tail =
        (static_cast<double>(upper & 0x7FFU) + (static_cast<double>(lower) * 0x1p-64)) * scale;
    const double product = head * HalfPiHigh;
    double       error =
        std::fma(head, HalfPiHigh, -product) + std::fma(head, HalfPiMiddle, tail * HalfPiHigh);
    const double sum = product + error;
    error            = error - (sum - product);
    return {.High = negative ? -sum : sum, .Low = negative ? -error : error, .Quarter = quarter};
  }

  auto TurnFar(double magnitude) -> double
  {
    constexpr std::uint64_t Fraction = (std::uint64_t{1} << 52U) - 1U;
    constexpr std::uint64_t Turn     = 360;

    const std::uint64_t bits        = BitsOf(magnitude);
    const int           exponent    = static_cast<int>(bits >> 52U) - 1075;
    const std::uint64_t significand = (bits & Fraction) | (std::uint64_t{1} << 52U);
    if (exponent < 0)
    {
      const auto drop = static_cast<unsigned>(-exponent);
      return static_cast<double>(significand % (Turn << drop))
           * OfBits(static_cast<std::uint64_t>(1023 + exponent) << 52U);
    }
    const auto          climb = static_cast<unsigned>(exponent);
    const std::uint64_t power =
        climb < 3U ? std::uint64_t{1} << climb : 8U * FortyFifths.at((climb - 3U) % 12U);
    return static_cast<double>((significand % Turn) * power % Turn);
  }
} // namespace Geoxy::Trigo
