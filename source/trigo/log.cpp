/***************************************************************************************************
 *  @file      source/trigo/log.cpp
 *  @brief     Logarithms, and what is built on them: the inverse hyperbolic sine and powers.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "real.hpp"
#include "trigo.hpp"
#include "trigo/reduce.hpp"
#include "twofold.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>

namespace Geoxy::Trigo
{
  namespace
  {
    /**
     *  ln 2 as a pair, its first part of 32 bits, so that a whole number of up to 21 bits times it
     *  is exact.
     */
    constexpr double LogTwoHigh = std::bit_cast<double>(std::uint64_t{0x3FE6'2E42'FEFA'3800});

    /**
     *  What LogTwoHigh leaves out of ln 2, nearest.
     */
    constexpr double LogTwoLow = 0x1.ef35793c7673p-45;

    /**
     *  1 / ln 2 as a pair.
     */
    constexpr double InverseLogTwoHigh = std::numbers::log2e;

    /**
     *  What InverseLogTwoHigh leaves out of 1 / ln 2, nearest.
     */
    constexpr double InverseLogTwoLow = 0x1.777d0ffda0d24p-56;

    /**
     *  2 / (2k + 1), k from 1 to 11: the series of 2 atanh(u) past its first term, in u^2.
     */
    constexpr std::array<double, 11> AtanhTerms = {
        2.0 / 3.0,  2.0 / 5.0,  2.0 / 7.0,  2.0 / 9.0,  2.0 / 11.0, 2.0 / 13.0,
        2.0 / 15.0, 2.0 / 17.0, 2.0 / 19.0, 2.0 / 21.0, 2.0 / 23.0,
    };

    /**
     *  ln(h + l) as a pair, h positive, finite and normal, l at most a step of a double of it: h =
     *  2^e m with m in [sqrt(1/2), sqrt(2)), f = m - 1 plus l's share, and ln m = 2 atanh(u), u = f
     *  / (2 + f) as a pair, |u| below 0.172, its series in u^2 summed in doubles past the first
     *  term; then e ln 2 added as a pair.
     *
     *  @param  high h.
     *  @param  low  l.
     *  @return      The logarithm, its pair within about 2^-60 of it past e ln 2's share.
     */
    auto Logarithm(double high, double low) -> Twofold::Two
    {
      constexpr std::uint64_t Fraction = (std::uint64_t{1} << 52U) - 1U;

      const std::uint64_t bits = BitsOf(high);
      /**
       *  m in [sqrt(1/2), sqrt(2)): the exponent rounded by the significand's first bits.
       */
      const std::uint64_t above = (bits & Fraction) > 0x6A09'E667'F3BCDU ? 1U : 0U;
      const int           power = static_cast<int>(bits >> 52U) - 1023 + static_cast<int>(above);
      const double        significand = OfBits((bits & Fraction) | ((1023U - above) << 52U));
      /**
       *  l is 0 but from Asinh, whose h is below 2^30: 2^-e is then a normal double.
       */
      const double scale =
          OfBits(static_cast<std::uint64_t>(std::clamp(1023 - power, 1, 2046)) << 52U);
      const double offset    = significand - 1.0;
      const double offsetLow = low * scale;
      const double below     = 2.0 + offset;
      const double under     = ((2.0 - below) + offset) + offsetLow;
      const double u         = offset / below;
      const double uLow      = (std::fma(-u, below, offset) + offsetLow - (u * under)) / below;
      const double v         = u * u;
      double       sum       = AtanhTerms.back();
      for (auto term = AtanhTerms.rbegin() + 1; term != AtanhTerms.rend(); ++term)
      {
        sum = std::fma(v, sum, *term);
      }
      const double       tail  = u * v * sum;
      const auto         whole = static_cast<double>(power);
      const Twofold::Two lead  = Twofold::FastTwoSum(whole * LogTwoHigh, 2.0 * u);
      return Twofold::FastTwoSum(lead.High, lead.Low + ((2.0 * uLow) + tail + (whole * LogTwoLow)));
    }

    /**
     *  2^r for r in [-1/2, 1/2]: e^t, t = r ln 2 as a pair, by its Taylor series to t^14 / 14!, 1
     *  kept apart from the rest and the pair's low part times e^t to first order.
     *
     *  @param  rest r.
     *  @return      2^r.
     */
    auto Exponential(double rest) -> double
    {
      constexpr std::array<double, 13> Terms = {
          1.0 / 2.0,          1.0 / 6.0,        1.0 / 24.0,        1.0 / 120.0,
          1.0 / 720.0,        1.0 / 5040.0,     1.0 / 40320.0,     1.0 / 362880.0,
          1.0 / 3628800.0,    1.0 / 39916800.0, 1.0 / 479001600.0, 1.0 / 6227020800.0,
          1.0 / 87178291200.0};

      /**
       *  ln 2 as a pair of two doubles of 53 bits.
       */
      constexpr double LogTwo     = std::numbers::ln2;
      constexpr double LogTwoRest = 0x1.abc9e3b39803fp-56;

      const double t    = rest * LogTwo;
      const double tLow = std::fma(rest, LogTwo, -t) + (rest * LogTwoRest);
      double       sum  = Terms.back();
      for (auto term = Terms.rbegin() + 1; term != Terms.rend(); ++term)
      {
        sum = std::fma(t, sum, *term);
      }
      const double square = t * t;
      /**
       *  e^(t + l) - 1 = (e^t - 1) + l e^t, l's share taken to its first order in t.
       */
      const double lead = 1.0 + t;
      return lead + (((1.0 - lead) + t) + ((square * sum) + (tLow + (tLow * t))));
    }

    /**
     *  log2 of a positive, finite and normal value as a pair: its natural logarithm times 1 / ln 2,
     *  both as pairs.
     *
     *  @param  value The value.
     *  @return       The logarithm.
     */
    auto BinaryLogarithm(double value) -> Twofold::Two
    {
      const Twofold::Two natural = Logarithm(value, 0.0);
      const double       high    = natural.High * InverseLogTwoHigh;
      const double       low     = std::fma(natural.High, InverseLogTwoHigh, -high)
                                 + (natural.Low * InverseLogTwoHigh)
                                 + (natural.High * InverseLogTwoLow);
      return Twofold::FastTwoSum(high, low);
    }

    /**
     *  2^(h + l), |h| below 1100 and l at most a step of a double of it: 2 to the nearest whole
     *  number n to h, exactly, in two halves, so that neither overflows nor leaves the normals
     *  before the product, times 2^((h - n) + l).
     *
     *  @param  high h.
     *  @param  low  l.
     *  @return      The power.
     */
    auto Power(double high, double low) -> double
    {
      const double whole  = (high + Shifter) - Shifter;
      const double scaled = Exponential((high - whole) + low);
      const auto   count  = static_cast<int>(whole);
      const int    half   = count / 2;
      return scaled
           * OfBits(static_cast<std::uint64_t>(1023 + half) << 52U)
           * OfBits(static_cast<std::uint64_t>(1023 + (count - half)) << 52U);
    }
  } // namespace

  auto Log2(double value) -> double
  {
    if (!(value > 0.0) || !(value < Infinite))
    {
      if (value < 0.0 || std::isnan(value))
      {
        return std::numeric_limits<double>::quiet_NaN();
      }
      return value > 0.0 ? value : -Infinite;
    }
    const bool         sub = value < 0x1p-1022;
    const Twofold::Two log = BinaryLogarithm(sub ? value * 0x1p64 : value);
    return (log.High + log.Low) - (sub ? 64.0 : 0.0);
  }

  auto Exp2(double power) -> double
  {
    if (!(std::abs(power) < 1100.0))
    {
      if (std::isnan(power))
      {
        return power;
      }
      return power > 0.0 ? Infinite : 0.0;
    }
    return Power(power, 0.0);
  }

  auto Pow(double base, double power) -> double
  {
    if (!(base >= 0x1p-1022) || !(base < Infinite))
    {
      return Exp2(power * Log2(base));
    }
    const Twofold::Two log     = BinaryLogarithm(base);
    const double       product = power * log.High;
    if (!(std::abs(product) < 1100.0))
    {
      return Exp2(product);
    }
    return Power(product, std::fma(power, log.High, -product) + (power * log.Low));
  }

  auto Asinh(double value) -> double
  {
    const double size = std::abs(value);
    if (!(size >= 0x1p-26) || !(size < Infinite))
    {
      return value;
    }
    /**
     *  ln(|x| + sqrt(x^2 + 1)), the sum as a pair; past 2^28, ln(2 |x|).
     */
    Twofold::Two sum = {.High = 2.0 * size, .Low = 0.0};
    if (size < 0x1p28)
    {
      const double       square = size * size;
      const Twofold::Two inside = Twofold::TwoSum(1.0, square);
      const double       root   = std::sqrt(inside.High);
      const double       ahead =
          (std::fma(-root, root, inside.High) + (inside.Low + std::fma(size, size, -square)))
          / (2.0 * root);
      const Twofold::Two total = Twofold::FastTwoSum(root, size);
      sum                      = Twofold::FastTwoSum(total.High, total.Low + ahead);
    }
    const Twofold::Two natural = Logarithm(sum.High, sum.Low);
    return std::copysign(natural.High + natural.Low, value);
  }

  auto Cbrt(double value) -> double
  {
    const double size = std::abs(value);
    if (!(size > 0.0) || !(size < Infinite))
    {
      return value + value;
    }
    /**
     *  Subnormals scaled up by 2^54 first, whose cube root is 2^18.
     */
    const bool   sub    = size < 0x1p-1022;
    const double scaled = sub ? size * 0x1p54 : size;
    /**
     *  A first guess from the bits, a third of the exponent, within 4%; Newton's steps square the
     *  error each; the last works from the cube's exact remainder.
     */
    double root = OfBits((BitsOf(scaled) / 3U) + 0x2A9F'7893'782D'A1CEU);
    for (int step = 0; step < 4; ++step)
    {
      root -= ((root * root * root) - scaled) / (3.0 * root * root);
    }
    /**
     *  Fused multiply-adds, not Twofold's products, which split their factors on a machine without
     *  them and round otherwise near underflow.
     */
    const double square = root * root;
    const double cube   = square * root;
    const double excess =
        ((cube - scaled) + std::fma(square, root, -cube)) + (std::fma(root, root, -square) * root);
    root -= excess / (3.0 * square);
    return std::copysign(sub ? root * 0x1p-18 : root, value);
  }
} // namespace Geoxy::Trigo
