/***************************************************************************************************
 *  @file      source/trigo/sine.cpp
 *  @brief     Sines, cosines and tangents, of radians and of degrees, from one kernel.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "real.hpp"
#include "trigo.hpp"
#include "trigo/reduce.hpp"
#include "twofold.hpp"

#include "geoxy/vector.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

#if defined(__aarch64__) && defined(__ARM_NEON)
#  include <arm_neon.h>
#endif

namespace Geoxy::Trigo
{
  namespace
  {
    /**
     *  sin r and cos r alike as m + m z (c1 + z P(z)), z = r^2, m = r for the sine and 1 for the
     *  cosine, plus what the remainder's low part adds, rl G(z) times the other row's m: two rows
     *  of coefficients, weighted minimax fits rounded one coefficient at a time. Each row holds c1
     *  as a pair, P's six coefficients, then G's three.
     */
    constexpr std::array<std::array<double, 11>, 2> SineRows = {{
        {-0x1.5555555555555p-3, -0x1.5555555555555p-57, 0x1.111111111110fp-7,
         -0x1.a01a01a01944dp-13, 0x1.71de3a53e0911p-19, -0x1.ae6453570b563p-26,
         0x1.6120ff2796473p-33, -0x1.aad5ed7220f30p-41, 1.0, -0.5, 0x1.5555555555555p-5},
        {-0.5, 0.0, 0x1.555555555554bp-5, -0x1.6c16c16c15015p-10, 0x1.a01a019c8f254p-16,
         -0x1.27e4f7f19148bp-22, 0x1.1ee9dbcefbddep-29, -0x1.8fa684873ff4bp-37, -1.0,
         0x1.5555555555555p-3, -0x1.1111111111111p-7},
    }};

    /**
     *  Per quadrant, m = a r + b and the other row's m, a' r + b', the quadrant's sign folded in:
     *  the sine of quadrant q is the sine row's for even q and the cosine row's for odd, so one
     *  kernel serves all four with no branch.
     */
    constexpr std::array<std::array<double, 4>, 4> SineQuadrants = {{
        {1.0, -0.0, 0.0, 1.0},
        {0.0, 1.0, 1.0, 0.0},
        {-1.0, 0.0, 0.0, -1.0},
        {0.0, -1.0, -1.0, 0.0},
    }};

    /**
     *  2^-27, below which the sine of an angle is the angle and its cosine 1, rounded.
     */
    constexpr std::uint64_t SineFloor = 0x3E40'0000'0000'0000U;

    /**
     *  2^30, from which Payne and Hanek reduce an angle.
     */
    constexpr std::uint64_t SineCeiling = 0x41D0'0000'0000'0000U;

    /**
     *  2^-900, below which an angle in degrees is scaled on its own, its radians near underflow.
     */
    constexpr std::uint64_t DegreeFloor = 0x07B0'0000'0000'0000U;

    /**
     *  2^40, from which an angle in degrees is taken modulo a turn in whole numbers.
     */
    constexpr std::uint64_t DegreeCeiling = 0x4270'0000'0000'0000U;

    /**
     *  Whether a double's magnitude is in [floor, ceiling): one comparison of bits, which a zero,
     *  an infinity and a not a number fail.
     *
     *  @param  bits    The double's bits.
     *  @param  floor   The floor's bits.
     *  @param  ceiling The ceiling's bits.
     *  @return         Whether it is.
     */
    auto Within(std::uint64_t bits, std::uint64_t floor, std::uint64_t ceiling) -> bool
    {
      return (bits << 1U) - (floor << 1U) < ((ceiling - floor) << 1U);
    }

    /**
     *  The sine of q pi / 2 + rh + rl as a pair, |rh| at most pi / 4: the leading m + m z c1 worked
     *  as a pair, every product's roundoff kept, which is what takes the error from libm's 0.87
     *  steps of a double to 0.56.
     *
     *  @param  high    rh.
     *  @param  low     rl.
     *  @param  quarter q.
     *  @return         The sine as a pair, its parts' sum within 2^-56.5 of it in relative terms.
     */
    auto SineKernel(double high, double low, std::uint64_t quarter) -> Twofold::Two
    {
      const std::array<double, 11> &row  = SineRows.at(quarter & 1U);
      const std::array<double, 4>  &sign = SineQuadrants.at(quarter & 3U);

      const double lead      = std::fma(sign[0], high, sign[1]);
      const double other     = std::fma(sign[2], high, sign[3]);
      const double square    = high * high;
      const double squareLow = std::fma(high, high, -square);
      const double cube      = lead * square;
      const double cubeLow   = std::fma(lead, squareLow, std::fma(lead, square, -cube));
      const double first     = cube * row[0];
      const double firstLow  = std::fma(cubeLow, row[0], std::fma(cube, row[0], -first));
      const double sum       = lead + first;
      const double sumLow    = (lead - sum) + first;
      const double fourth    = square * square;
      const double rest      = std::fma(
          fourth,
          std::fma(fourth, std::fma(square, row[7], row[6]), std::fma(square, row[5], row[4])),
          std::fma(square, row[3], row[2]));
      const double tail  = cube * std::fma(square, rest, row[1]);
      const double slope = other * std::fma(square, std::fma(square, row[10], row[9]), row[8]);
      return {.High = sum, .Low = sumLow + (firstLow + std::fma(low, slope, tail))};
    }

    /**
     *  The sine of q pi / 2 + rh + rl, rounded.
     *
     *  @param  reduced rh, rl and q.
     *  @param  shift   1 for the cosine, 0 for the sine.
     *  @return         It.
     */
    auto SineOf(const Reduced &reduced, std::uint64_t shift) -> double
    {
      const Twofold::Two sine = SineKernel(reduced.High, reduced.Low, reduced.Quarter + shift);
      return sine.High + sine.Low;
    }

#if defined(__aarch64__) && defined(__ARM_NEON)
    /**
     *  The two rows side by side, coefficient by coefficient: lane 0 the sine's, lane 1 the
     *  cosine's.
     */
    constexpr std::array<std::array<double, 2>, 11> SineLanes =
        [] -> std::array<std::array<double, 2>, 11>
    {
      std::array<std::array<double, 2>, 11> lanes = {};
      for (std::size_t term = 0; term < lanes.size(); ++term)
      {
        lanes.at(term) = {SineRows[0].at(term), SineRows[1].at(term)};
      }
      return lanes;
    }();

    /**
     *  Per quadrant, the signs of the sine and the cosine of the remainder, as bits.
     */
    constexpr std::array<std::array<std::uint64_t, 2>, 4> SineSigns = {{
        {0, 0},
        {0, 0x8000'0000'0000'0000U},
        {0x8000'0000'0000'0000U, 0x8000'0000'0000'0000U},
        {0x8000'0000'0000'0000U, 0},
    }};

    /**
     *  A coefficient of both rows.
     *
     *  @param  term Which.
     *  @return      Its two lanes.
     */
    auto SineLane(std::size_t term) -> float64x2_t
    {
      return vld1q_f64(SineLanes.at(term).data());
    }

    /**
     *  A pair in two lanes.
     */
    struct SinePair
    {
      float64x2_t High; /**< The leading parts. */
      float64x2_t Low;  /**< What they leave.   */
    };

    /**
     *  The sine and the cosine of rh + rl, both rows of SineKernel in one pass, lane by lane, the
     *  same operations: the same bits, for about the price of one.
     *
     *  @param  high rh.
     *  @param  low  rl.
     *  @return      The sine in lane 0 and the cosine in lane 1, each as a pair.
     */
    auto SineLanesOf(double high, double low) -> SinePair
    {
      const float64x2_t lead    = vsetq_lane_f64(high, vdupq_n_f64(1.0), 0);
      const float64x2_t other   = vextq_f64(lead, lead, 1);
      const double      square  = high * high;
      const double      below   = std::fma(high, high, -square);
      const double      fourth  = square * square;
      const float64x2_t first   = SineLane(0);
      const float64x2_t cube    = vmulq_n_f64(lead, square);
      const float64x2_t cubeOff = vfmsq_n_f64(vfmsq_n_f64(cube, lead, square), lead, below);
      const float64x2_t term    = vmulq_f64(cube, first);
      const float64x2_t termOff = vfmaq_f64(vfmsq_f64(term, cube, first), cubeOff, first);
      const float64x2_t sum     = vaddq_f64(lead, term);
      const float64x2_t sumLow  = vaddq_f64(vsubq_f64(lead, sum), term);
      const float64x2_t rest =
          vfmaq_n_f64(vfmaq_n_f64(SineLane(2), SineLane(3), square),
                      vfmaq_n_f64(vfmaq_n_f64(SineLane(4), SineLane(5), square),
                                  vfmaq_n_f64(SineLane(6), SineLane(7), square), fourth),
                      fourth);
      const float64x2_t tail  = vmulq_f64(cube, vfmaq_n_f64(SineLane(1), rest, square));
      const float64x2_t slope = vmulq_f64(
          other, vfmaq_n_f64(SineLane(8), vfmaq_n_f64(SineLane(9), SineLane(10), square), square));
      return {.High = sum,
              .Low  = vaddq_f64(sumLow, vsubq_f64(vfmaq_n_f64(tail, slope, low), termOff))};
    }

    /**
     *  The cosine and the sine of q pi / 2 + rh + rl: both lanes summed, swapped in an odd
     *  quadrant, then signed.
     *
     *  @param  reduced rh, rl and q.
     *  @return         The cosine and the sine.
     */
    auto WayOf(const Reduced &reduced) -> Vector
    {
      const SinePair    pair   = SineLanesOf(reduced.High, reduced.Low);
      const float64x2_t sum    = vaddq_f64(pair.High, pair.Low);
      const float64x2_t odd    = vreinterpretq_f64_u64(vdupq_n_u64(0U - (reduced.Quarter & 1U)));
      const float64x2_t turned = vbslq_f64(vreinterpretq_u64_f64(odd), vextq_f64(sum, sum, 1), sum);
      const float64x2_t way    = vreinterpretq_f64_u64(veorq_u64(
          vreinterpretq_u64_f64(turned), vld1q_u64(SineSigns.at(reduced.Quarter & 3U).data())));
      return {.X = vgetq_lane_f64(way, 1), .Y = vgetq_lane_f64(way, 0)};
    }

    /**
     *  The tangent of q pi / 2 + rh + rl: lane 0 the sine over the cosine, lane 1 the cosine over
     *  the sine, each corrected by the pairs' remainder: the scalar path's operations, lane by
     *  lane.
     *
     *  @param  reduced rh, rl and q.
     *  @return         The tangent.
     */
    auto TangentOf(const Reduced &reduced) -> double
    {
      const SinePair    pair     = SineLanesOf(reduced.High, reduced.Low);
      const float64x2_t sum      = vaddq_f64(pair.High, pair.Low);
      const float64x2_t below    = vextq_f64(sum, sum, 1);
      const float64x2_t quotient = vdivq_f64(sum, below);
      const float64x2_t inverse  = vdivq_f64(vdupq_n_f64(1.0), below);
      const float64x2_t remainder =
          vaddq_f64(vfmsq_f64(pair.High, quotient, vextq_f64(pair.High, pair.High, 1)),
                    vfmsq_f64(pair.Low, quotient, vextq_f64(pair.Low, pair.Low, 1)));
      const float64x2_t tangent = vfmaq_f64(quotient, remainder, inverse);
      return (reduced.Quarter & 1U) != 0 ? -vgetq_lane_f64(tangent, 1) : vgetq_lane_f64(tangent, 0);
    }
#else
    /**
     *  The cosine and the sine of q pi / 2 + rh + rl.
     *
     *  @param  reduced rh, rl and q.
     *  @return         The cosine and the sine.
     */
    auto WayOf(const Reduced &reduced) -> Vector
    {
      return {.X = SineOf(reduced, 1), .Y = SineOf(reduced, 0)};
    }

    /**
     *  The tangent of q pi / 2 + rh + rl: the sine over the cosine of the remainder, or less the
     *  cosine over the sine in an odd quadrant, the quotient corrected by the pairs' remainder.
     *
     *  @param  reduced rh, rl and q.
     *  @return         It.
     */
    auto TangentOf(const Reduced &reduced) -> double
    {
      const Twofold::Two sine   = SineKernel(reduced.High, reduced.Low, 0);
      const Twofold::Two cosine = SineKernel(reduced.High, reduced.Low, 1);
      const double       sum    = sine.High + sine.Low;
      const double       cosum  = cosine.High + cosine.Low;
      if ((reduced.Quarter & 1U) != 0)
      {
        const double quotient = cosum / sum;
        const double inverse  = 1.0 / sum;
        const double remainder =
            std::fma(-quotient, sine.High, cosine.High) + std::fma(-quotient, sine.Low, cosine.Low);
        return -std::fma(remainder, inverse, quotient);
      }
      const double quotient = sum / cosum;
      const double inverse  = 1.0 / cosum;
      const double remainder =
          std::fma(-quotient, cosine.High, sine.High) + std::fma(-quotient, cosine.Low, sine.Low);
      return std::fma(remainder, inverse, quotient);
    }
#endif

    /**
     *  An angle in radians outside [2^-27, 2^30) in size reduced: a tiny one is its own sine and
     *  tangent, with a cosine of 1; a larger one goes through Payne and Hanek.
     *
     *  @param  radians The angle.
     *  @return         Its reduction, or nothing for a tiny one.
     */
    auto Outside(double radians) -> Reduced
    {
      Reduced reduced = ReduceFar(std::abs(radians));
      if (radians < 0.0)
      {
        reduced.High    = -reduced.High;
        reduced.Low     = -reduced.Low;
        reduced.Quarter = 0U - reduced.Quarter;
      }
      return reduced;
    }

    /**
     *  A tiny angle in degrees in radians, correctly rounded below the normal doubles too: the
     *  product scaled up by 2^200, then down, and where it lands below the normals rounded by hand
     *  to the nearest, by every step of the subnormals the pair's low part is worth; pi / 180 being
     *  irrational, no product lies halfway.
     *
     *  @param  degrees The angle, less than 2^-900 in size.
     *  @return         Its radians.
     */
    auto Tiny(double degrees) -> double
    {
      if (Zero(degrees))
      {
        return degrees;
      }
      const double scaled  = degrees * 0x1p200;
      const double product = scaled * DegreeHigh;
      const double error   = std::fma(scaled, DegreeLow, std::fma(scaled, DegreeHigh, -product));
      const double down    = product * 0x1p-200;
      if (std::abs(down) >= 0x1p-1022)
      {
        return (product + error) * 0x1p-200;
      }
      /**
       *  What rounding the product into the subnormals left, with the pair's low part: the result
       *  moved by the step, if any, it rounds to.
       */
      const double left = (product - (down * 0x1p200)) + error;
      return down + (std::round(left * 0x1p874) * 0x1p-1074);
    }
  } // namespace

  auto Sin(double radians) -> double
  {
    if (!Within(BitsOf(radians), SineFloor, SineCeiling)) [[unlikely]]
    {
      if (std::abs(radians) < 0x1p-27)
      {
        return radians;
      }
      return SineOf(Outside(radians), 0);
    }
    return SineOf(Reduce(radians), 0);
  }

  auto Cos(double radians) -> double
  {
    if (!Within(BitsOf(radians), SineFloor, SineCeiling)) [[unlikely]]
    {
      if (std::abs(radians) < 0x1p-27)
      {
        return 1.0;
      }
      return SineOf(Outside(radians), 1);
    }
    return SineOf(Reduce(radians), 1);
  }

  auto Tan(double radians) -> double
  {
    if (!Within(BitsOf(radians), SineFloor, SineCeiling)) [[unlikely]]
    {
      if (std::abs(radians) < 0x1p-27)
      {
        return radians;
      }
      return TangentOf(Outside(radians));
    }
    return TangentOf(Reduce(radians));
  }

  auto Way(double radians) -> Vector
  {
    if (!Within(BitsOf(radians), SineFloor, SineCeiling)) [[unlikely]]
    {
      if (std::abs(radians) < 0x1p-27)
      {
        return {.X = 1.0, .Y = radians};
      }
      return WayOf(Outside(radians));
    }
    return WayOf(Reduce(radians));
  }

  auto WayDegrees(double degrees) -> Vector
  {
    if (!Within(BitsOf(degrees), DegreeFloor, DegreeCeiling)) [[unlikely]]
    {
      const double magnitude = std::abs(degrees);
      if (!(magnitude < Infinite))
      {
        return {.X = std::numeric_limits<double>::quiet_NaN(),
                .Y = std::numeric_limits<double>::quiet_NaN()};
      }
      if (magnitude < 0x1p-900)
      {
        return {.X = 1.0, .Y = Tiny(degrees)};
      }
      return WayDegrees(std::copysign(TurnFar(magnitude), degrees));
    }
    /**
     *  k quarter turns off, exactly: the remainder is a whole number of the angle's last bit,
     *  smaller than the angle, so one fused multiply-add gives it with no rounding.
     */
    const double  shifted   = std::fma(degrees, 0x1.6c16c16c16c17p-7, Shifter);
    const double  quarters  = shifted - Shifter;
    const double  remainder = std::fma(-quarters, 90.0, degrees);
    const double  high      = remainder * DegreeHigh;
    const Reduced reduced   = {
        .High    = high,
        .Low     = std::fma(remainder, DegreeLow, std::fma(remainder, DegreeHigh, -high)),
        .Quarter = BitsOf(shifted)};
    const Vector way = WayOf(reduced);
    return {.X = way.X + 0.0, .Y = Zero(way.Y) ? std::copysign(0.0, degrees) : way.Y};
  }
} // namespace Geoxy::Trigo
