/***************************************************************************************************
 *  @file      source/turning.hpp
 *  @brief     What angles held in degrees are worked with: whole turns taken off exactly, and the
 *             exact sign of a short sum of doubles.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TURNING_HPP
#define GEOXY_TURNING_HPP

#include "exact/strict.hpp" // IWYU pragma: keep
#include "twofold.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <optional>

/**
 *  The steps behind angles held in degrees. Each exact step is proved in tests/proof/angle.v, or
 *  rests on Knuth's two-sum, which tests/proof/rounding.v proves.
 */
namespace Geoxy::Turning
{
  /**
   *  Degrees less whole turns, exactly, into (-180, 180].
   *
   *  @param  degrees The degrees, finite.
   *  @return         The same direction; -0 for a whole number of turns from -0.
   */
  inline auto Wrap(double degrees) -> double
  {
    /**
     *  In (-540, 540] one turn at most, exact by tests/proof/angle.v, Wrap; past it, fmod is exact.
     */
    if (degrees <= -540.0 || degrees > 540.0)
    {
      degrees = std::fmod(degrees, 360.0);
    }
    /**
     *  A select, not a branch: the step is 0, exact, when none is needed.
     */
    const double rise = degrees <= -180.0 ? 360.0 : 0.0;
    const double step = degrees > 180.0 ? -360.0 : rise;
    return degrees + step;
  }

  /**
   *  The double nearest a value known as a pair, where the pair's error cannot carry it across a
   *  halfway point between doubles: Ziv's test.
   *
   *  @param  high  The pair's leading part, not smaller than the rest.
   *  @param  low   The rest.
   *  @param  error How far the pair may be from the value, relatively; 2^-84 or more.
   *  @return       The double nearest the value; none where the pair is too near a halfway point or
   *                the error bound, times the sum, is below 2^-1000. A pair of zeros is an exact
   *                zero.
   */
  inline auto Settled(double high, double low, double error) -> std::optional<double>
  {
    const double rounded  = high + low;
    const double residual = (high - rounded) + low;
    const double size     = std::abs(rounded);
    /**
     *  A sum that rounds to zero is exact, the residual zero too.
     */
    if (!(size > 0.0))
    {
      return rounded;
    }
    /**
     *  The residual pushed on by the error bound and a hair, 2^-20 of it: the value is nearer the
     *  rounded sum than that, on the residual's side, and within the bound on the other. Where the
     *  push still rounds back onto the rounded sum, no halfway point lies between, rounding itself
     *  minding that the doubles are twice as close below a power of two. The hair, at least 2^-104
     *  of the sum for a bound of 2^-84, covers this test's own roundings, below 2^-104 of it. All
     *  of it must be a normal double, or the hair would round away and a tie would pass: so the
     *  pushed bound must be 2^-1000 or more.
     */
    const double bound = (error * (1.0 + 0x1p-20)) * size;
    const double reach = residual + std::copysign(bound, residual);
    if (bound >= 0x1p-1000 && std::equal_to<double>{}(rounded + reach, rounded))
    {
      return rounded;
    }
    return std::nullopt;
  }

  /**
   *  Radians as degrees less whole turns, in [-180, 180], as a pair.
   *
   *  Payne and Hanek's reduction, as Direction makes it: within 2^-104 of the exact degrees, which
   *  no double of radians can give exactly, a turn being 2 pi.
   *
   *  @param  radians The radians, finite.
   *  @return         The degrees as a pair, -180 only from just above it.
   */
  auto Reduced(double radians) -> Twofold::Two;

  /**
   *  The cosine and sine of degrees held as a pair, each a pair not yet rounded.
   *
   *  Reduced as Direction reduces degrees, exactly to an eighth of a degree, and the table's pairs
   *  turned by what is left, each summed into a pair whose small part is at most half an ulp of its
   *  large: within 2^-68 of the exact value relatively. Below 2^-500 degrees the cosine is 1 and
   *  the sine the degrees' radians, to 2^-990 of themselves.
   *
   *  @param  degrees The degrees, finite, the small part at most half an ulp of the large.
   *  @return         The cosine, then the sine.
   */
  auto CosinePair(Twofold::Two degrees) -> std::array<Twofold::Two, 2>;

  /**
   *  The sign of a sum of doubles, exactly.
   *
   *  Shewchuk's Grow-Expansion adds each term into a nonoverlapping expansion, whose sign is that
   *  of its largest part not zero (Shewchuk 1997, Theorems 10 and 13). Sums are exact under
   *  underflow too; the terms must not overflow.
   *
   *  @param  terms The terms.
   *  @return       1, 0 or -1.
   */
  template <std::size_t CountT>
  [[nodiscard]] auto SumSign(const std::array<double, CountT> &terms) -> int
  {
    std::array<double, CountT> parts = {};
    std::size_t                count = 0;
    for (const double term : terms)
    {
      double      running = term;
      std::size_t kept    = 0;
      for (std::size_t index = 0; index < count; ++index)
      {
        const Twofold::Two sum = Twofold::TwoSum(running, parts.at(index));
        if (sum.Low < 0.0 || sum.Low > 0.0)
        {
          parts.at(kept++) = sum.Low;
        }
        running = sum.High;
      }
      parts.at(kept++) = running;
      count            = kept;
    }
    for (std::size_t index = count; index > 0; --index)
    {
      if (parts.at(index - 1) > 0.0)
      {
        return 1;
      }
      if (parts.at(index - 1) < 0.0)
      {
        return -1;
      }
    }
    return 0;
  }
} // namespace Geoxy::Turning

#endif /* GEOXY_TURNING_HPP */
