/***************************************************************************************************
 *  @file      source/trigo/arc.cpp
 *  @brief     Arc tangents, sines and cosines: one division, a region found from the bits, and a
 *             short kernel.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "real.hpp"
#include "trigo.hpp"
#include "trigo/reduce.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>

namespace Geoxy::Trigo
{
  namespace
  {
    /**
     *  n pi / 4 as a pair, n from 0 to 8.
     */
    constexpr std::array<std::array<double, 2>, 9> Octants = {{
        {0.0, 0.0},
        {0x1.921fb54442d18p-1, 0x1.1a62633145c07p-55},
        {0x1.921fb54442d18p+0, 0x1.1a62633145c07p-54},
        {0x1.2d97c7f3321d2p+1, 0x1.a79394c9e8a0ap-54},
        {std::numbers::pi, 0x1.1a62633145c07p-53},
        {0x1.f6a7a2955385ep+1, 0x1.60fafbfd97309p-53},
        {0x1.2d97c7f3321d2p+2, 0x1.a79394c9e8a0ap-53},
        {0x1.5fdbbe9bba775p+2, 0x1.ee2c2d963a10cp-53},
        {0x1.921fb54442d18p+2, 0x1.1a62633145c07p-52},
    }};

    /**
     *  Per region j, c_j and atan(c_j) as a pair: c of two significant bits, so that the numerator
     *  min - c max is exact in one fused multiply-add, and the reduced tangent is within 0.1365 of
     *  0. The last region, beyond 2^-3.5, has c = 0.
     */
    constexpr std::array<std::array<double, 3>, 8> ArcRegions = {{
        {0x1.0000000000000p+0, 0x1.921fb54442d18p-1, 0x1.1a62633145c07p-55},
        {0x1.8000000000000p-1, 0x1.4978fa3269ee1p-1, 0x1.2419a87f2a458p-56},
        {0x1.0000000000000p-1, 0x1.dac670561bb4fp-2, 0x1.a2b7f222f65e2p-56},
        {0x1.8000000000000p-2, 0x1.6f61941e4def1p-2, -0x1.c63aae6f6e918p-56},
        {0x1.0000000000000p-2, 0x1.f5b75f92c80ddp-3, 0x1.8ab6e3cf7afbdp-57},
        {0x1.8000000000000p-3, 0x1.7b97b4bce5b02p-3, 0x1.347b0b4f881cap-58},
        {0x1.0000000000000p-3, 0x1.fd5ba9aac2f6ep-4, -0x1.cd37686760c17p-59},
        {0.0, 0.0, 0.0},
    }};

    /**
     *  1 as bits.
     */
    constexpr std::uint64_t OneBits = 0x3FF0'0000'0000'0000U;

    /**
     *  The sign bit.
     */
    constexpr std::uint64_t SignBit = 0x8000'0000'0000'0000U;

    /**
     *  n pi / 4 as a pair, n from -4 to 4.
     *
     *  @param  octant n.
     *  @return        It.
     */
    auto Octant(int octant) -> std::array<double, 2>
    {
      const std::array<double, 2> &pair = Octants.at(static_cast<std::size_t>(std::abs(octant)));
      return octant < 0 ? std::array<double, 2>{-pair[0], -pair[1]} : pair;
    }

    /**
     *  An offset plus a signed pair, rounded once: o + s (h + l), where |o| is at least |s h| or o
     *  is 0.
     *
     *  @param  pair   h and l.
     *  @param  sign   s, a power of two of either sign.
     *  @param  offset o as a pair.
     *  @return        The sum.
     */
    auto Offset(std::array<double, 2> pair, double sign, std::array<double, 2> offset) -> double
    {
      const double sum   = std::fma(sign, pair[0], offset[0]);
      const double error = std::fma(sign, pair[0], offset[0] - sum);
      return sum + std::fma(sign, pair[1], error + offset[1]);
    }

    /**
     *  n pi / 4 + s atan(least / most), least <= most both positive: region j = round(-2 log2(least
     *  / most)) read from the difference of their bits, a log2 in pieces; then atan(c) +
     *  atan((least
     *  - c most) / (most + c least)), the second a short odd polynomial of the reduced tangent.
     *
     *  @param  ratio  least and most.
     *  @param  octant n.
     *  @param  sign   s.
     *  @return        The angle.
     */
    auto ArcCore(std::array<double, 2> ratio, int octant, double sign) -> double
    {
      const double        least  = ratio[0];
      const double        most   = ratio[1];
      const std::uint64_t spread = BitsOf(most) - BitsOf(least);
      const std::uint64_t region =
          std::min<std::uint64_t>((spread + (std::uint64_t{1} << 50U)) >> 51U, 7U);
      const std::array<double, 3> &row = ArcRegions.at(region);

      const double c           = row[0];
      const double numerator   = std::fma(-c, most, least);
      const double denominator = std::fma(c, least, most);
      const double below       = std::fma(c, least, most - denominator);
      const double inverse     = 1.0 / denominator;
      const double high        = numerator * inverse;
      const double remainder   = std::fma(-high, below, std::fma(-high, denominator, numerator));
      const double low         = remainder * inverse;
      const double square      = high * high;
      const double fourth      = square * square;
      const double eighth      = fourth * fourth;
      const double odd =
          std::fma(eighth,
                   std::fma(fourth, -0x1.ff59a5d1eb36ap-5,
                            std::fma(square, 0x1.3a8e0c07e6b72p-4, -0x1.745aff71c6627p-4)),
                   std::fma(fourth, std::fma(square, 0x1.c71c6d2750047p-4, -0x1.24924921d5fecp-3),
                            std::fma(square, 0x1.9999999998108p-3, -0x1.5555555555553p-2)));
      const double rest = std::fma(high * square, odd, std::fma(-low, square, low));
      const double sum  = row[1] + high;
      return Offset({sum, ((row[1] - sum) + high) + (row[2] + rest)}, sign, Octant(octant));
    }

    /**
     *  The multiple of pi / 4 a displacement with a zero or an infinite side, and no not a number,
     *  lies at, before its rise's sign.
     *
     *  @param  rise The rise.
     *  @param  run  The run.
     *  @return      The multiple, from 0 to 4.
     */
    auto EdgeOctant(double rise, double run) -> int
    {
      const double across = std::abs(run);
      const double rising = std::abs(rise);
      if (Zero(rising) || (std::isinf(across) && !std::isinf(rising)))
      {
        return std::signbit(run) ? 4 : 0;
      }
      if (std::isinf(rising) && std::isinf(across))
      {
        return std::signbit(run) ? 3 : 1;
      }
      return 2;
    }

    /**
     *  The angle of a displacement whose sides lie 2^60 or more apart in size, where atan(v) is v
     *  to the last bit: n pi / 4 + s v, v the smaller side over the larger.
     *
     *  @param  rise The rise.
     *  @param  run  The run.
     *  @return      Its angle.
     */
    auto ArcFar(double rise, double run) -> double
    {
      const double across = std::abs(run);
      const double rising = std::abs(rise);
      const bool   steep  = rising > across;
      const double ratio  = steep ? across / rising : rising / across;
      int          far    = steep ? 2 : 0;
      double       sign   = steep ? -1.0 : 1.0;
      if (std::signbit(run))
      {
        far  = 4 - far;
        sign = -sign;
      }
      if (std::signbit(rise))
      {
        far  = -far;
        sign = -sign;
      }
      if (far == 0)
      {
        return sign * ratio;
      }
      return Offset({ratio, 0.0}, sign, Octant(far));
    }

    /**
     *  The angle of a displacement outside the quick path: a zero, an infinity or a not a number in
     *  it, sides far apart in size, or the larger beyond [2^-960, 2^1000), brought near 1 by a
     *  power of two, which changes no angle.
     *
     *  @param  rise The rise.
     *  @param  run  The run.
     *  @return      Its angle.
     */
    auto ArcOutside(double rise, double run) -> double
    {
      if (std::isnan(run) || std::isnan(rise))
      {
        return run + rise;
      }
      const double across = std::abs(run);
      const double rising = std::abs(rise);
      if (Zero(rising) || Zero(across) || std::isinf(rising) || std::isinf(across))
      {
        return std::copysign(Octants.at(static_cast<std::size_t>(EdgeOctant(rise, run)))[0], rise);
      }
      int runPower  = 0;
      int risePower = 0;
      static_cast<void>(std::frexp(across, &runPower));
      static_cast<void>(std::frexp(rising, &risePower));
      if (std::abs(runPower - risePower) >= 60)
      {
        return ArcFar(rise, run);
      }
      const int power = std::max(runPower, risePower) - 1;
      return Atan2(std::ldexp(rise, -power), std::ldexp(run, -power));
    }

    /**
     *  asin(uh + ul) as a pair, u in [0, 0.5001]: u + u^3 / 6 as a pair, every product's roundoff
     *  kept, plus u^3 times a polynomial of u^2 and what ul adds, ul / sqrt(1 - u^2).
     *
     *  @param  high uh.
     *  @param  low  ul.
     *  @return      The angle as a pair.
     */
    auto SineArc(double high, double low) -> std::array<double, 2>
    {
      constexpr double SixthHigh = 0x1.5555555555555p-3;
      constexpr double SixthLow  = 0x1.5555555555555p-57;

      const double square    = high * high;
      const double squareLow = std::fma(high, high, -square);
      const double cube      = high * square;
      const double cubeLow   = std::fma(high, squareLow, std::fma(high, square, -cube));
      const double first     = cube * SixthHigh;
      const double firstLow  = std::fma(cubeLow, SixthHigh, std::fma(cube, SixthHigh, -first));
      const double sum       = high + first;
      const double sumLow    = (high - sum) + first;
      const double fourth    = square * square;
      const double eighth    = fourth * fourth;
      const double sixteenth = eighth * eighth;
      const double zero      = std::fma(square, 0x1.6db6db6d7ff95p-5, 0x1.333333333356ep-4);
      const double two       = std::fma(square, 0x1.6e8b9d58e7aabp-6, 0x1.f1c71c934e7d1p-6);
      const double four      = std::fma(square, 0x1.c985a7b17e0b9p-7, 0x1.1c4f572b00884p-6);
      const double six       = std::fma(square, 0x1.37e617718c38fp-7, 0x1.7b77e9fb11f22p-7);
      const double eight     = std::fma(square, 0x1.38afd90d03315p-10, 0x1.42c3f4c8f4ad9p-7);
      const double ten       = std::fma(square, -0x1.96af4a3e8fba0p-6, 0x1.8534451376101p-6);
      const double lower     = std::fma(fourth, two, zero);
      const double middle    = std::fma(fourth, six, four);
      const double upper     = std::fma(fourth, ten, eight);
      const double series    = std::fma(sixteenth, std::fma(eighth, 0x1.0d93fd748e1a6p-5, upper),
                                        std::fma(eighth, middle, lower));
      const double tail      = cube * std::fma(square, series, SixthLow);
      const double slope     = std::fma(square, std::fma(square, 0.375, 0.5), 1.0);
      return {sum, sumLow + (firstLow + std::fma(low, slope, tail))};
    }

    /**
     *  The arc sine or arc cosine of a value in [2^-27, 1] in size: of |x| itself up to 1/2, and
     *  past it as pi / 2 - 2 asin(u), u = sqrt((1 - |x|) / 2) as a pair, exact under the root.
     *
     *  @param  value  x.
     *  @param  cosine Whether the arc cosine.
     *  @return        The angle.
     */
    auto SineArcOf(double value, bool cosine) -> double
    {
      const std::uint64_t bits      = BitsOf(value);
      const std::uint64_t magnitude = bits & ~SignBit;
      const bool          big       = magnitude > 0x3FE0'0000'0000'0000U;
      const bool          negative  = (bits >> 63U) != 0;
      const double        size      = OfBits(magnitude);
      const double        half      = std::fma(size, -0.5, 0.5);
      const double        root      = std::sqrt(half);
      const double        high      = big ? root : size;
      const double        low = big ? std::fma(-root, root, half) / (root + root + 0x1p-1000) : 0.0;
      const std::array<double, 2> arc = SineArc(high, low);

      int  octant = 0;
      bool flip   = false;
      if (!cosine)
      {
        octant = big ? 2 : 0;
        flip   = big != negative;
        octant = negative ? -octant : octant;
      }
      else if (!big)
      {
        octant = 2;
        flip   = !negative;
      }
      else
      {
        octant = negative ? 4 : 0;
        flip   = negative;
      }
      const double sign = (flip ? -1.0 : 1.0) * (big ? 2.0 : 1.0);
      return Offset(arc, sign, Octant(octant));
    }

    /**
     *  The arc sine or arc cosine outside [2^-27, 1] in size: not a number past 1 or for one; for a
     *  tiny value, the value itself, or pi / 2 less it.
     *
     *  @param  value  x.
     *  @param  cosine Whether the arc cosine.
     *  @return        The angle.
     */
    auto SineArcOutside(double value, bool cosine) -> double
    {
      if (!(std::abs(value) <= 1.0))
      {
        return std::numeric_limits<double>::quiet_NaN();
      }
      return cosine ? Offset({value, 0.0}, -1.0, Octants[2]) : value;
    }
  } // namespace

  auto Atan(double tangent) -> double
  {
    const std::uint64_t bits      = BitsOf(tangent);
    const std::uint64_t magnitude = bits & ~SignBit;
    if (magnitude - 0x3E40'0000'0000'0000U >= 0x41A0'0000'0000'0000U - 0x3E40'0000'0000'0000U)
        [[unlikely]]
    {
      if (std::abs(tangent) < 0x1p-27)
      {
        return tangent;
      }
      return ArcOutside(tangent, 1.0);
    }
    const bool steep    = magnitude > OneBits;
    const bool negative = (bits >> 63U) != 0;
    int        octant   = steep ? 2 : 0;
    if (negative)
    {
      octant = -octant;
    }
    return ArcCore({OfBits(steep ? OneBits : magnitude), OfBits(steep ? magnitude : OneBits)},
                   octant, steep != negative ? -1.0 : 1.0);
  }

  auto Atan2(double rise, double run) -> double
  {
    const std::uint64_t runBits  = BitsOf(run);
    const std::uint64_t riseBits = BitsOf(rise);
    const std::uint64_t across   = runBits & ~SignBit;
    const std::uint64_t rising   = riseBits & ~SignBit;
    const bool          steep    = rising > across;
    const std::uint64_t most     = steep ? rising : across;
    const std::uint64_t least    = steep ? across : rising;
    if (least           == 0
        || most         >= 0x7E70'0000'0000'0000U
        || most         < 0x03F0'0000'0000'0000U
        || most - least >= (std::uint64_t{60} << 52U)) [[unlikely]]
    {
      return ArcOutside(rise, run);
    }
    const bool back    = (runBits >> 63U) != 0;
    const bool down    = (riseBits >> 63U) != 0;
    int        octant  = steep ? 2 : 0;
    octant             = back ? 4 - octant : octant;
    octant             = down ? -octant : octant;
    const bool flipped = (steep != back) != down;
    return ArcCore({OfBits(least), OfBits(most)}, octant, flipped ? -1.0 : 1.0);
  }

  auto Asin(double sine) -> double
  {
    const std::uint64_t magnitude = BitsOf(sine) & ~SignBit;
    if (magnitude < 0x3E40'0000'0000'0000U || magnitude > OneBits) [[unlikely]]
    {
      return SineArcOutside(sine, false);
    }
    return SineArcOf(sine, false);
  }

  auto Acos(double cosine) -> double
  {
    const std::uint64_t magnitude = BitsOf(cosine) & ~SignBit;
    if (magnitude < 0x3E40'0000'0000'0000U || magnitude > OneBits) [[unlikely]]
    {
      return SineArcOutside(cosine, true);
    }
    return SineArcOf(cosine, true);
  }
} // namespace Geoxy::Trigo
