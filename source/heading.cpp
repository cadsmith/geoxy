/***************************************************************************************************
 *  @file      source/heading.cpp
 *  @brief     Headings in radians: an exact fold to the first octant, a table of arctangents at
 *             every 1/128, and a Taylor polynomial whose remainder is below 2^-83; and the exact
 *             order of two headings, from signs alone.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "arctangent.hpp"
#include "finite.hpp"
#include "precise.hpp"
#include "radians.hpp"
#include "real.hpp"
#include "turning.hpp"
#include "twofold.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>

namespace Geoxy
{
  namespace
  {
    using Twofold::FastTwoSum;
    using Twofold::Two;
    using Twofold::TwoProduct;
    using Twofold::TwoSum;

    constexpr double TableSteps = 128.0;    /**< Table steps in [0, 1].          */
    constexpr double AtanStep   = 0x1p-7;   /**< A table step.                   */
    constexpr double Faint      = 0x1p-400; /**< Below it atan t is t to 2^-800. */
    constexpr double Up         = 0x1p600;  /**< Lifts a tiny ratio.             */
    constexpr double Down       = 0x1p-600; /**< Lowers it back.                 */

    /**
     *  atan u = u + u^3 (ArcThird + u^2 (ArcFifth + u^2 (ArcSeventh + u^2 ArcNinth))) + O(u^11):
     *  for |u| <= 2^-8 the remainder is below 2^-83 |u|.
     */
    constexpr double ArcThird   = -1.0 / 3.0;
    constexpr double ArcFifth   = 1.0 / 5.0;
    constexpr double ArcSeventh = -1.0 / 7.0;
    constexpr double ArcNinth   = 1.0 / 9.0;

    /**
     *  A pair negated.
     *
     *  @param  pair The pair.
     *  @return      Its negative.
     */
    auto Countered(Two pair) -> Two
    {
      return {.High = -pair.High, .Low = -pair.Low};
    }

    /**
     *  A power of two as a double.
     *
     *  @param  exponent The exponent, from -1022 to 1023.
     *  @return          2^exponent.
     */
    auto PowerOfTwo(int exponent) -> double
    {
      return std::bit_cast<double>(static_cast<std::uint64_t>(exponent + 1023) << 52U);
    }

    /**
     *  A displacement scaled by a power of two, exactly, so its larger component is in [1, 2).
     *
     *  @pre           Its smaller nonzero component does not fall below 2^-1022 doing so.
     *  @param  vector The displacement, finite and not zero.
     *  @return        The same heading, scaled.
     */
    auto Unitized(Vector vector) -> Vector
    {
      double larger = std::max(std::abs(vector.X), std::abs(vector.Y));
      if (larger < 0x1p-1000)
      {
        vector  = {.X = vector.X * Up, .Y = vector.Y * Up};
        larger *= Up;
      }
      else if (larger > 0x1p1000)
      {
        vector  = {.X = vector.X * Down, .Y = vector.Y * Down};
        larger *= Down;
      }
      const auto   bits     = std::bit_cast<std::uint64_t>(larger);
      const int    exponent = static_cast<int>((bits >> 52U) & 0x7ffU) - 1023;
      const double scale    = PowerOfTwo(-exponent);
      return {.X = vector.X * scale, .Y = vector.Y * scale};
    }

    /**
     *  A tiny angle: atan(over / under) where the ratio is below 2^-400, so atan t is t to 2^-800.
     *  Worked on the two scaled into [1, 2) and scaled back by the difference of their exponents,
     *  exact wherever the angle is a normal double.
     *
     *  @param  over  The smaller, a pair.
     *  @param  under The larger, a pair, at least 2^-500.
     *  @return       The angle in radians, a pair within 2^-100 of it where it is normal.
     */
    auto Slight(Two over, Two under) -> Two
    {
      if (Zero(over.High))
      {
        return {};
      }
      int lift = 0;
      if (over.High < 0x1p-1000)
      {
        over = {.High = over.High * Up, .Low = over.Low * Up};
        lift = 600;
      }
      const int    overExponent  = std::ilogb(over.High);
      const int    underExponent = std::ilogb(under.High);
      const double overBig       = std::scalbn(over.High, -overExponent);
      const double overSmall     = std::scalbn(over.Low, -overExponent);
      const double underBig      = std::scalbn(under.High, -underExponent);
      const double underSmall    = std::scalbn(under.Low, -underExponent);
      const double ratio         = overBig / underBig;
      const double remainder     = std::fma(-ratio, underBig, overBig);
      const double ratioLow      = (remainder + overSmall - (ratio * underSmall)) / underBig;
      const int    shift         = overExponent - underExponent - lift;
      return {.High = std::scalbn(ratio, shift), .Low = std::scalbn(ratioLow, shift)};
    }

    /**
     *  What a quotient of pairs leaves of the dividend, over / under less the ratio times under, to
     *  first order: exact for the big parts by a fused multiply-add.
     *
     *  @param  over  The dividend, a pair; its small part 0 unless PairT.
     *  @param  under The divisor, likewise.
     *  @param  ratio The quotient of the big parts, rounded.
     *  @return       The remainder.
     */
    template <bool PairT>
    auto Remainder(Two over, Two under, double ratio) -> double
    {
      const double remainder = std::fma(-ratio, under.High, over.High);
      if constexpr (PairT)
      {
        return remainder + over.Low - (ratio * under.Low);
      }
      else
      {
        return remainder;
      }
    }

    template <bool PairT>
    /**
     *  The first octant's heading, atan(over / under), as a pair.
     *
     *  @param  over  The smaller, a pair; its small part 0 unless PairT.
     *  @param  under The larger, a pair, at least 2^-500 and 2^400 times the smaller.
     *  @return       The angle in radians, as a pair.
     */
    auto Octant(Two over, Two under) -> Two
    {
      /**
       *  t = over / under as a pair, from one reciprocal: the pair need only be within 2^-100 of t,
       *  so the quotient need not be the nearest, only at most 1 (tests/proof/angle.v, Arc).
       */
      const double inverse  = 1.0 / under.High;
      const double ratio    = std::min(over.High * inverse, 1.0);
      const double ratioLow = Remainder<PairT>(over, under, ratio) * inverse;

      /**
       *  Exact: tests/proof/angle.v, Arc. u = (t - c) / (1 + t c), |u| <= 2^-8.
       */
      const double index       = std::rint(ratio * TableSteps);
      const double near        = index * AtanStep;
      const Two    top         = TwoSum(ratio - near, ratioLow);
      const Two    cross       = TwoProduct(ratio, near);
      const Two    denominator = FastTwoSum(1.0, cross.High);
      const double baseLow     = denominator.Low + cross.Low + (ratioLow * near);
      const double reverse     = 1.0 / denominator.High;
      const double uHigh       = top.High * reverse;
      const double uRem        = std::fma(-uHigh, denominator.High, top.High);
      const double uLow        = (uRem + top.Low - (uHigh * baseLow)) * reverse;

      const double square = uHigh * uHigh;
      const double tail =
          std::fma(-uLow, square, uLow)
          + ((uHigh * square)
             * (ArcThird + (square * (ArcFifth + (square * (ArcSeventh + (square * ArcNinth)))))));

      const Two &arc = Table::ArcRadians.at(static_cast<std::size_t>(index));
      const Two  sum = TwoSum(arc.High, uHigh);
      return {.High = sum.High, .Low = sum.Low + (arc.Low + tail)};
    }

    /**
     *  A displacement's components as sums of doubles, exactly, for the slow way.
     */
    struct Terms
    {
      std::span<const double> Rise = {}; /**< Doubles that sum to the y component. */
      std::span<const double> Run  = {}; /**< Doubles that sum to the x component. */
    };

    /**
     *  The heading of a displacement held in pairs, in radians.
     *
     *  @pre          Both finite, not both zero.
     *  @param  y     Along the y axis, a pair whose small part is at most half its big part's ulp;
     *                0 unless PairT.
     *  @param  x     Along the x axis, likewise.
     *  @param  terms The components as sums, for the slow way.
     *  @return       The double nearest the heading, in [-pi, pi]: -pi only where a heading just
     *                above it rounds there, y < 0 and x < 0 with y below an ulp of x.
     */
    template <bool PairT>
    auto Bearing(Two y, Two x, Terms terms) -> double
    {
      /**
       *  The pair settled by Ziv's test, or, near a halfway point, the heading to any precision, by
       *  a call to another file, so the fast path's frame stays in registers.
       */
      const auto nearest = [terms](double high, double low) -> double
      {
        if (const std::optional<double> settled = Turning::Settled(high, low, 0x1p-68)) [[likely]]
        {
          return *settled + 0.0;
        }
        return Precise::NearestArc({.High = high, .Low = low}, terms.Rise, terms.Run) + 0.0;
      };
      const bool below = y.High < 0.0;
      const bool back  = x.High < 0.0;
      const Two  high  = below ? Countered(y) : y;
      const Two  wide  = back ? Countered(x) : x;
      /**
       *  Whether over the diagonal: the sign of high - wide, exact, the big parts' difference
       *  leading wherever they differ.
       */
      const bool steep = (high.High - wide.High) + (high.Low - wide.Low) > 0.0;
      Two        over  = steep ? wide : high;
      Two        under = steep ? high : wide;
      if (under.High < 0x1p-500)
      {
        /**
         *  Both lifted, exactly, so the division's remainder below is a double.
         */
        over  = {.High = over.High * Up, .Low = over.Low * Up};
        under = {.High = under.High * Up, .Low = under.Low * Up};
      }

      /**
       *  The first octant: atan(over / under), as a pair.
       */
      double angle = 0.0;
      double rest  = 0.0;
      if (over.High < under.High * Faint)
      {
        /**
         *  Off the axis by less than 2^-390 degrees: the heading is that angle on the positive x
         *  axis, and the axis's own angle otherwise, being within an ulp's tiniest part of it.
         */
        if (!steep && !back)
        {
          const Two slight = Slight(over, under);
          return below ? nearest(-slight.High, -slight.Low) : nearest(slight.High, slight.Low);
        }
      }
      else
      {
        const Two octant = Octant<PairT>(over, under);
        angle            = octant.High;
        rest             = octant.Low;
      }

      /**
       *  Unfold the octant: base + side * angle, rounded once; selects, not branches, as the octant
       *  of a heading is as good as random.
       */
      const Two    sideways = back ? Table::Pi : Two{};
      const Two    base     = steep ? Table::HalfPi : sideways;
      const double side     = steep == back ? 1.0 : -1.0;
      const Two    sum      = TwoSum(base.High, side * angle);
      const double tail     = sum.Low + (base.Low + (side * rest));
      return below ? nearest(-sum.High, -tail) : nearest(sum.High, tail);
    }

    /**
     *  Which half-open stretch of headings a displacement is in, in order: (-180, 0), 0, (0, 180)
     *  and 180.
     *
     *  @param  vector The displacement.
     *  @return        0 to 3; no displacement is 1, at 0.
     */
    auto Half(Vector vector) -> int
    {
      if (vector.Y < 0.0)
      {
        return 0;
      }
      if (vector.Y > 0.0)
      {
        return 2;
      }
      return vector.X < 0.0 ? 3 : 1;
    }

    /**
     *  The sign of the cross product of two displacements, exactly.
     *
     *  @param  lhs One displacement, scaled.
     *  @param  rhs The other, scaled.
     *  @return     1 when the second is counterclockwise of the first, 0 when along it, else -1.
     */
    auto Sidedness(Vector lhs, Vector rhs) -> int
    {
      /**
       *  Rounding is monotone, so unequal rounded products order as the exact ones do; equal ones
       *  differ by their roundoffs, whose difference keeps its sign through rounding.
       */
      const Two first  = TwoProduct(lhs.X, rhs.Y);
      const Two second = TwoProduct(lhs.Y, rhs.X);
      if (first.High > second.High)
      {
        return 1;
      }
      if (first.High < second.High)
      {
        return -1;
      }
      const double apart = first.Low - second.Low;
      return static_cast<int>(apart > 0.0) - static_cast<int>(apart < 0.0);
    }

    /**
     *  The difference of two pairs, Joldes, Muller and Popescu's AccurateDWPlusDW: within 3 u^2 of
     *  the exact difference relatively, so exactly zero when it is.
     *
     *  @param  lhs One pair.
     *  @param  rhs The pair taken from it.
     *  @return     The difference, a pair.
     */
    auto Minus(Two lhs, Two rhs) -> Two
    {
      const Two big   = TwoSum(lhs.High, -rhs.High);
      const Two small = TwoSum(lhs.Low, -rhs.Low);
      const Two first = FastTwoSum(big.High, big.Low + small.High);
      return FastTwoSum(first.High, small.Low + first.Low);
    }

    /**
     *  Whether a displacement's products with another's keep exact roundoffs and do not overflow:
     *  each component zero or between 2^-480 and 2^480 in size.
     *
     *  @param  vector The displacement.
     *  @return        Whether it needs no scaling.
     */
    auto Moderate(Vector vector) -> bool
    {
      const auto fits = [](double value) -> bool
      {
        const double size = std::abs(value);
        return Zero(size) || (size >= 0x1p-480 && size <= 0x1p480);
      };
      return fits(vector.X) && fits(vector.Y);
    }

    /**
     *  Whether a displacement is none.
     *
     *  @param  vector The displacement.
     *  @return        Whether both components are zero.
     */
    auto Null(Vector vector) -> bool
    {
      return Zero(vector.X) && Zero(vector.Y);
    }

    constexpr double NotANumber = std::numeric_limits<double>::quiet_NaN();
  } // namespace

  auto Heading(Vector vector) -> Angle
  {
    if (!Finite(vector))
    {
      if (std::isnan(vector.X) || std::isnan(vector.Y))
      {
        return Radians(NotANumber);
      }
      const auto unit = [](double value) -> double
      {
        return std::isinf(value) ? std::copysign(1.0, value) : 0.0;
      };
      vector = {.X = unit(vector.X), .Y = unit(vector.Y)};
    }
    if (Null(vector))
    {
      return Radians(0.0);
    }
    const std::array<double, 1> rise = {vector.Y};
    const std::array<double, 1> run  = {vector.X};
    return Radians(Bearing<false>({.High = vector.Y + 0.0, .Low = 0.0},
                                  {.High = vector.X + 0.0, .Low = 0.0},
                                  {.Rise = rise, .Run = run}));
  }

  auto Between(Vector start, Vector end) -> Angle
  {
    if (!Finite(start) || !Finite(end))
    {
      return Radians(NotANumber);
    }
    if (Null(start) || Null(end))
    {
      return Radians(0.0);
    }
    const bool   moderate   = Moderate(start) && Moderate(end);
    const Vector from       = moderate ? start : Unitized(start);
    const Vector onto       = moderate ? end : Unitized(end);
    const Two    crossLeft  = TwoProduct(from.X, onto.Y);
    const Two    crossRight = TwoProduct(from.Y, onto.X);
    const Two    dotLeft    = TwoProduct(from.X, onto.X);
    const Two    dotRight   = TwoProduct(from.Y, onto.Y);
    const Two    cross      = Minus(crossLeft, crossRight);
    const Two    dot        = Minus(dotLeft, {.High = -dotRight.High, .Low = -dotRight.Low});
    const std::array<double, 4> rise = {crossLeft.High, crossLeft.Low, -crossRight.High,
                                        -crossRight.Low};
    const std::array<double, 4> run  = {dotLeft.High, dotLeft.Low, dotRight.High, dotRight.Low};
    return Radians(Bearing<true>({.High = cross.High + 0.0, .Low = cross.Low},
                                 {.High = dot.High + 0.0, .Low = dot.Low},
                                 {.Rise = rise, .Run = run}));
  }

  auto ArcTangent(double value) -> Twofold::Two
  {
    const Two    one   = {.High = 1.0, .Low = 0.0};
    const Two    size  = {.High = std::abs(value), .Low = 0.0};
    const bool   steep = size.High > 1.0;
    const Two    angle = steep ? Octant<false>(one, size) : Octant<false>(size, one);
    const Two    turn  = steep ? Table::HalfPi : Two{};
    const double side  = steep ? -1.0 : 1.0;
    const Two    sum   = TwoSum(turn.High, side * angle.High);
    const Two    out   = FastTwoSum(sum.High, sum.Low + (turn.Low + (side * angle.Low)));
    return value < 0.0 ? Two{.High = -out.High, .Low = -out.Low} : out;
  }

  auto Order(Vector lhs, Vector rhs) -> std::strong_ordering
  {
    const int lhsHalf = Half(lhs);
    const int rhsHalf = Half(rhs);
    if (lhsHalf != rhsHalf || lhsHalf == 1 || lhsHalf == 3)
    {
      return lhsHalf <=> rhsHalf;
    }
    if (Moderate(lhs) && Moderate(rhs))
    {
      return 0 <=> Sidedness(lhs, rhs);
    }
    return 0 <=> Sidedness(Unitized(lhs), Unitized(rhs));
  }
} // namespace Geoxy
