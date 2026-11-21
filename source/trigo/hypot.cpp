/***************************************************************************************************
 *  @file      source/trigo/hypot.cpp
 *  @brief     The length of a displacement, the same bits on every machine.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "real.hpp"
#include "trigo.hpp"

#include <cmath>

namespace Geoxy::Trigo
{
  namespace
  {
    /**
     *  sqrt(most^2 + least^2), most at least least, both finite and the squares neither overflowing
     *  nor, where it matters, underflowing: the sum of squares as a pair, its root, then one Newton
     *  step from the pair.
     *
     *  @param  most  The larger magnitude.
     *  @param  least The smaller.
     *  @return       The length.
     */
    auto Reach(double most, double least) -> double
    {
      const double square = most * most;
      const double over   = std::fma(most, most, -square);
      const double other  = least * least;
      const double under  = std::fma(least, least, -other);
      const double sum    = square + other;
      const double rest   = ((square - sum) + other) + (over + under);
      const double root   = std::sqrt(sum);
      return root + ((std::fma(-root, root, sum) + rest) / (2.0 * root));
    }
  } // namespace

  auto Hypot(double run, double rise) -> double
  {
    const double across = std::abs(run);
    const double rising = std::abs(rise);
    const bool   steep  = rising > across;
    const double most   = steep ? rising : across;
    const double least  = steep ? across : rising;
    /**
     *  A not a number fails the comparison, so it is the smaller only where the other is too and
     *  never reaches the root. From 2^-500 to 2^500 no square overflows, and one that underflows is
     *  less than a step of a double of the other: no scaling, which costs as much as the rest.
     */
    if (most >= 0x1p-500 && most < 0x1p500) [[likely]]
    {
      return Reach(most, least);
    }
    if (std::isinf(across) || std::isinf(rising))
    {
      return Infinite;
    }
    if (std::isnan(across) || std::isnan(rising))
    {
      return across + rising;
    }
    if (Zero(most))
    {
      return 0.0;
    }
    /**
     *  Scaled by a power of two, which is exact, into the range above, and back.
     */
    const double scale = most < 1.0 ? 0x1p600 : 0x1p-600;
    return Reach(most * scale, least * scale) / scale;
  }
} // namespace Geoxy::Trigo
