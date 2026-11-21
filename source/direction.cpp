/***************************************************************************************************
 *  @file      source/direction.cpp
 *  @brief     The cosine and sine: degrees reduced exactly to an eighth of a degree, radians by
 *             Payne and Hanek's reduction to 2^-201, then a table of pairs and Taylor polynomials
 *             whose remainder is below 2^-94.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "degrees.hpp"
#include "precise.hpp"
#include "radians.hpp"
#include "real.hpp"
#include "turning.hpp"
#include "twofold.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>

namespace Geoxy
{
  namespace
  {
    using Twofold::TwoProduct;
    using Twofold::TwoSum;

    constexpr double Ninetieth  = 1.0 / 90.0; /**< RN(1/90): picks the quarter turn.         */
    constexpr double Reducible  = 0x1p45;     /**< Below it the quarter turn is picked well. */
    constexpr double RightAngle = 90.0;       /**< A quarter turn.                           */
    constexpr double Eighths    = 8.0;        /**< Table steps in a degree.                  */
    constexpr double TableStep  = 0.125;      /**< A table step, in degrees.                 */
    constexpr double Minute     = 0x1p-500;   /**< Below it cos is 1 and sin x is x.         */
    constexpr double Lift       = 0x1p600;    /**< Keeps a tiny product's roundoff a double. */
    constexpr double Drop       = 0x1p-600;

    /**
     *  sin x = x + x^3 (SineThird + x^2 (SineFifth + x^2 SineSeventh)) + O(x^9): for |x| <= 2^-9.8
     *  the remainder is below 2^-97 |x|. cos x = 1 + x^2 (CosineSecond + x^2 (CosineFourth + x^2
     *  CosineSixth)) + O(x^8): below 2^-94.
     */
    constexpr double SineThird    = -1.0 / 6.0;
    constexpr double SineFifth    = 1.0 / 120.0;
    constexpr double SineSeventh  = -1.0 / 5040.0;
    constexpr double CosineSecond = -0.5;
    constexpr double CosineFourth = 1.0 / 24.0;
    constexpr double CosineSixth  = -1.0 / 720.0;

    constexpr double PairError =
        0x1p-68; /**< How far the pairs are from the exact values, at most. */

    /**
     *  A cosine and a sine, each a pair of doubles not yet rounded.
     */
    struct Pairs
    {
      Twofold::Two X = {}; /**< The cosine. */
      Twofold::Two Y = {}; /**< The sine.   */
    };

    auto Core(std::int64_t index, double xHigh, double xLow) -> Pairs;
    auto Quarterly(Pairs pairs, int quarter) -> Pairs;

    /**
     *  The nearest doubles to a cosine and a sine: by Ziv's test on the pairs, else slowly.
     *
     *  @param  pairs The pairs, within 2^-68 of the exact values relatively.
     *  @param  slow  Calls Precise for the nearest doubles, given the pairs, only when the test
     *                fails: a call to another file, so the fast path's frame stays in registers.
     *  @return       The nearest doubles; a zero component +0.
     */
    template <typename SlowT>
    auto Settle(const Pairs &pairs, const SlowT &slow) -> Vector
    {
      const std::optional<double> cosine = Turning::Settled(pairs.X.High, pairs.X.Low, PairError);
      const std::optional<double> sine   = Turning::Settled(pairs.Y.High, pairs.Y.Low, PairError);
      if (cosine && sine) [[likely]]
      {
        return {.X = *cosine + 0.0, .Y = *sine + 0.0};
      }
      const std::array<double, 2> nearest = slow(Precise::Way{.Cosine = pairs.X, .Sine = pairs.Y});
      return {.X = nearest[0] + 0.0, .Y = nearest[1] + 0.0};
    }

    /**
     *  The part of a degree left and the pair's small part, summed exactly: the part alone when the
     *  degrees are one double.
     *
     *  @param  part The part left.
     *  @param  low  The pair's small part; 0 unless PairT.
     *  @return      Their sum as a pair.
     */
    template <bool PairT>
    auto Totaled(double part, double low) -> Twofold::Two
    {
      if constexpr (PairT)
      {
        return TwoSum(part, low);
      }
      else
      {
        return {.High = part, .Low = low};
      }
    }

    /**
     *  Degrees reduced: what is left of a quarter turn's multiple, exactly, and its cosine and sine
     *  as pairs.
     */
    struct Reduction
    {
      Pairs  Turned  = {};  /**< The cosine and sine, quarter turns added. */
      double Rest    = 0.0; /**< The degrees left, in [-45, 45].           */
      double Low     = 0.0; /**< What the pair's leading part leaves out.  */
      int    Quarter = 0;   /**< Quarter turns, 0 to 3.                    */
    };

    /**
     *  Finite degrees of 2^-500 or more reduced, and their cosine and sine as pairs.
     *
     *  @param  high The leading degrees.
     *  @param  low  The rest, at most half an ulp of the leading part; 0 unless PairT.
     *  @return      The reduction.
     */
    template <bool PairT>
    auto Reduced(double high, double low) -> Reduction
    {
      if (!(std::abs(high) < Reducible))
      {
        const Twofold::Two whole = TwoSum(std::fmod(high, 360.0), std::fmod(low, 360.0));
        high                     = whole.High;
        low                      = whole.Low;
      }

      /**
       *  Exact: tests/proof/angle.v, Quadrant and Split. rest is in [-45, 45], part in [-1/16,
       *  1/16].
       */
      const double quarter = std::rint(high * Ninetieth);
      const double rest    = high - (RightAngle * quarter);
      const double eighth  = std::rint(rest * Eighths);
      const double split   = rest - (eighth * TableStep);

      /**
       *  The part left, in radians, as a pair x = xHigh + xLow, |xLow| <= 2^-52 |xHigh|.
       */
      const Twofold::Two part  = Totaled<PairT>(split, low);
      const Twofold::Two lead  = TwoProduct(part.High, Table::ToRadians.High);
      const double       xHigh = lead.High;
      const double       xLow =
          lead.Low + ((part.High * Table::ToRadians.Low) + (part.Low * Table::ToRadians.High));

      const int turns =
          static_cast<int>(static_cast<std::uint64_t>(static_cast<std::int64_t>(quarter)) & 3U);
      return {.Turned  = Quarterly(Core(static_cast<std::int64_t>(eighth), xHigh, xLow), turns),
              .Rest    = rest,
              .Low     = low,
              .Quarter = turns};
    }

    /**
     *  The nearest doubles to the cosine and sine of degrees held in two doubles.
     *
     *  @param  high The leading degrees.
     *  @param  low  The rest, at most half an ulp of the leading part; 0 unless PairT.
     *  @return      The cosine and sine; a zero component +0.
     */
    template <bool PairT>
    auto Cosined(double high, double low) -> Vector
    {
      /**
       *  Infinite or NaN: the pair's sum is, the small part being below an ulp of the big.
       */
      if (!std::isfinite(high + low))
      {
        constexpr double None = std::numeric_limits<double>::quiet_NaN();
        return {.X = None, .Y = None};
      }
      if (std::abs(high) < Minute)
      {
        /**
         *  cos x = 1 - x^2 / 2 and sin x = x (1 - x^2 / 6) with x^2 below 2^-990: 1 and, for one
         *  double, its radians, the nearest double (tests/proof/convert.py) being that of sin x
         *  too; for a pair, its radians to 2^-50 and the slow way where that falls short.
         */
        if (!PairT || Zero(low))
        {
          return {.X = 1.0, .Y = Degrees(high).InRadians() + 0.0};
        }
        const double lifted = high * Lift;
        const double sine =
            std::fma(lifted, Table::ToRadians.High,
                     (lifted * Table::ToRadians.Low) + (low * Lift * Table::ToRadians.High));
        const double nearest = Precise::NearestSine(sine * Drop, {.High = high, .Low = low});
        return {.X = 1.0, .Y = nearest + 0.0};
      }
      const Reduction reduced = Reduced<PairT>(high, low);
      return Settle(reduced.Turned,
                    [&reduced](Precise::Way way) -> std::array<double, 2>
                    {
                      return Precise::NearestDegree(way, {.High = reduced.Rest, .Low = reduced.Low},
                                                    reduced.Quarter);
                    });
    }

    /**
     *  The cosine and sine of a + x: a an eighth of a degree's multiple, from the table, and x a
     *  pair of radians within 1/16 of a degree and a hair.
     *
     *  @param  index The eighths of a degree in a, from -360 to 360.
     *  @param  xHigh x's leading part.
     *  @param  xLow  The rest, at most 2^-52 of the leading part.
     *  @return       The cosine and sine as pairs, within 2^-68 relatively.
     */
    auto Core(std::int64_t index, double xHigh, double xLow) -> Pairs
    {
      const Table::Eighth &entry =
          Table::Eighths.at(static_cast<std::size_t>(index < 0 ? -index : index));
      const double sign     = index < 0 ? -1.0 : 1.0;
      const double sinBig   = sign * entry.Sine.High;
      const double sinSmall = sign * entry.Sine.Low;
      const double cosBig   = entry.Cosine.High;
      const double cosSmall = entry.Cosine.Low;

      const double square = xHigh * xHigh;
      const double sinTail =
          (xHigh * square) * (SineThird + (square * (SineFifth + (square * SineSeventh))));
      const double cosTail =
          square * (CosineSecond + (square * (CosineFourth + (square * CosineSixth))));
      const double sinX = xLow + sinTail;           // sin x - xHigh
      const double cosX = cosTail - (xHigh * xLow); // cos x - 1

      /**
       *  sin(a + x) = sin a + cos a xHigh + (cos a (sin x - xHigh) + sin a (cos x - 1)), the two
       *  leading terms summed exactly, all else rounded once into the final sum.
       */
      const Twofold::Two sinLead = TwoProduct(cosBig, xHigh);
      const double       sinRest =
          sinLead.Low + sinSmall + (cosSmall * xHigh) + (cosBig * sinX) + (sinBig * cosX);
      const Twofold::Two sinSum = TwoSum(sinBig, sinLead.High);

      /**
       *  cos(a + x) = cos a - sin a xHigh + (cos a (cos x - 1) - sin a (sin x - xHigh)).
       */
      const Twofold::Two cosLead = TwoProduct(sinBig, xHigh);
      const double       cosRest =
          cosSmall - cosLead.Low - (sinSmall * xHigh) + (cosBig * cosX) - (sinBig * sinX);
      const Twofold::Two cosSum = TwoSum(cosBig, -cosLead.High);
      return {.X = {.High = cosSum.High, .Low = cosSum.Low + cosRest},
              .Y = {.High = sinSum.High, .Low = sinSum.Low + sinRest}};
    }
  } // namespace

  namespace
  {
    constexpr std::size_t Words    = 9;    /**< Words of 2 / pi a reduction multiplies by. */
    constexpr double      Snapping = 16.0; /**< Below it, a quarter turn's nearest double is it. */
    constexpr double      EighthsPerRadian = 0x1.ca5dc1a63c1f8p+8; /**< RN(1440 / pi).          */
    constexpr double      EighthFirst      = 0x1.1df46a2529d39p-9; /**< pi / 1440's first double. */
    constexpr double      EighthSecond     = 0x1.5c1d8becdd291p-65; /**< Its second.              */
    constexpr double      EighthThird    = -0x1.1d937fa428858p-119; /**< Its third, to 2^-172.    */
    constexpr double      InverseQuarter = 0x1.45f306dc9c883p-1;  /**< RN(2 / pi).                */
    constexpr double      QuarterFirst   = 0x1.921fb54442d18p+0;  /**< pi / 2's first double.     */
    constexpr double      QuarterSecond  = 0x1.1a62633145c07p-54; /**< Its second.                */
    constexpr double      QuarterThird = -0x1.f1976b7ed8fbcp-110; /**< Its third, to 2^-163.     */
    constexpr std::uint64_t Mask = 0xffffffffU; /**< A word's 32 bits.                         */

    /**
     *  Radians reduced: whole quarter turns, mod 4, and what is left in degrees, as a pair.
     */
    struct Folded
    {
      int    Quarter = 0;   /**< Quarter turns, 0 to 3.             */
      double High    = 0.0; /**< The rest in degrees, in [-45, 45]. */
      double Low     = 0.0; /**< What the leading part leaves out.  */
    };

    /**
     *  A bit of a whole number held in 32-bit words, least significant first.
     *
     *  @param  words The words.
     *  @param  bit   The bit's place.
     *  @return       The bit.
     */
    auto BitOf(const std::array<std::uint64_t, 12> &words, int bit) -> std::uint64_t
    {
      return (words.at(static_cast<std::size_t>(bit / 32)) >> static_cast<unsigned>(bit % 32)) & 1U;
    }

    /**
     *  So many bits of a whole number held in 32-bit words, from a place up.
     *
     *  @param  words The words, least significant first.
     *  @param  from  The lowest bit's place, 0 or more.
     *  @param  count How many, at most 53.
     *  @return       The bits as a whole number.
     */
    auto BitsOf(const std::array<std::uint64_t, 12> &words, int from, int count) -> std::uint64_t
    {
      std::uint64_t bits = 0;
      for (int place = from + count - 1; place >= from; --place)
      {
        bits = (bits << 1U) | BitOf(words, place);
      }
      return bits;
    }

    /**
     *  Positive radians of at least 3/4 reduced by quarter turns, Payne and Hanek's way.
     *
     *  radians = m 2^e with m a whole number of 53 bits, and radians 2 / pi = m sum_j w_j 2^(e - 32
     *  (j + 1)), w_j the words of 2 / pi. Words whose least weight is 4 or more add whole turns and
     *  are skipped; nine words after them are multiplied in whole numbers, exactly, and the words
     *  left out weigh below 2^-201 of a quarter turn. No double is nearer a multiple of a quarter
     *  turn than 2^-62 of one, so the rest is known to 2^-139 of itself, and its degrees as a pair
     *  to 2^-104.
     *
     *  @param  radians The radians, finite, at least 3/4.
     *  @return         The quarter turns and the degrees left.
     */
    auto Hanek(double radians) -> Folded
    {
      int          exponent = 0;
      const double fraction = std::frexp(radians, &exponent);
      const auto   whole    = static_cast<std::uint64_t>(std::ldexp(fraction, 53));
      const int    power    = exponent - 53;
      const int    first    = power >= 2 ? (power - 2) / 32 : 0;

      std::array<std::uint64_t, 12> sum   = {};
      const auto                    start = static_cast<std::size_t>(first);
      for (std::size_t index = 0; index < Words; ++index)
      {
        const std::uint64_t word   = Table::TwoOverPi.at(start + index);
        const std::uint64_t lower  = (whole & Mask) * word;
        const std::uint64_t upper  = (whole >> 32U) * word;
        const std::size_t   place  = Words - 1 - index;
        sum.at(place)             += lower & Mask;
        sum.at(place + 1)         += (lower >> 32U) + (upper & Mask);
        sum.at(place + 2)         += upper >> 32U;
      }
      for (std::size_t place = 0; place + 1 < sum.size(); ++place)
      {
        sum.at(place + 1) += sum.at(place) >> 32U;
        sum.at(place)     &= Mask;
      }

      /**
       *  The fraction's bits are the lowest `shift`; the two above them are the quarter turns.
       */
      const int  shift   = (32 * (first + static_cast<int>(Words))) - power;
      int        quarter = static_cast<int>(BitsOf(sum, shift, 2));
      const bool over    = BitOf(sum, shift - 1) != 0U;
      if (over)
      {
        /**
         *  The rest as 1 - fraction, a quarter turn on: the low `shift` bits complemented, plus 1.
         */
        quarter             = static_cast<int>((static_cast<unsigned>(quarter) + 1U) & 3U);
        std::uint64_t carry = 1;
        for (std::uint64_t &word : sum)
        {
          word   = ((~word) & Mask) + carry;
          carry  = word >> 32U;
          word  &= Mask;
        }
      }
      /**
       *  The fraction's highest bit set, a word at a time.
       */
      /**
       *  It is never 0, nor below 2^-63, no double being nearer a multiple of a quarter turn than
       *  2^-62 of one (Muller, Elementary Functions); so it lies in the top three words, the bits
       *  read below it are at 2^-170 and more, and the loop ends there.
       */
      int           word = (shift - 1) / 32;
      std::uint64_t bits = sum.at(static_cast<std::size_t>(word))
                         & ((std::uint64_t{2} << static_cast<unsigned>((shift - 1) % 32)) - 1U);
      while (bits == 0U)
      {
        --word;
        bits = sum.at(static_cast<std::size_t>(word));
      }
      const int    top = (32 * word) + 63 - std::countl_zero(bits);
      const double rest =
          std::ldexp(static_cast<double>(BitsOf(sum, top - 52, 53)), top - 52 - shift);
      /**
       *  Below 16 radians the double nearest a quarter turn was taken for it before: by Waite, and
       *  by Middling, which hand on only rests of 2^-30 and more.
       */
      const double high = rest;
      const double low =
          std::ldexp(static_cast<double>(BitsOf(sum, top - 105, 53)), top - 105 - shift);
      const Twofold::Two degrees = TwoProduct(high, RightAngle);
      const Twofold::Two pair    = TwoSum(degrees.High, degrees.Low + (low * RightAngle));
      const double       sign    = over ? -1.0 : 1.0;
      return {.Quarter = quarter, .High = sign * pair.High, .Low = sign * pair.Low};
    }

    /**
     *  Positive radians of at least 3/4 and below 2^20 reduced by quarter turns, Cody and Waite's
     *  way, when the rest is not tiny.
     *
     *  k = rint(radians 2 / pi) is below 2^20; pi / 2 is held in three doubles to 2^-163. The first
     *  product is taken off by one fused multiply-add, exactly: radians and k p1 are multiples of
     *  2^-53 and their difference is below 1. The second product is exact as a pair, and the third
     *  and pi / 2's truncation leave below 2^-142. Where the rest is at least 2^-30, that is 2^-112
     *  of it; nearer a quarter turn, Hanek decides.
     *
     *  @param  radians The radians, at least 3/4 and below 2^20.
     *  @return         The quarter turns and the radians left, as a pair, or none where the rest is
     *                  tiny.
     */
    auto Waite(double radians) -> std::optional<Folded>
    {
      const double       count  = std::rint(radians * InverseQuarter);
      const double       first  = std::fma(-count, QuarterFirst, radians);
      const Twofold::Two second = TwoProduct(count, QuarterSecond);
      const Twofold::Two head   = TwoSum(first, -second.High);
      const double       tail   = (head.Low - second.Low) - (count * QuarterThird);
      const Twofold::Two rest   = TwoSum(head.High, tail);
      const double       size   = std::abs(rest.High);
      const int          turn =
          static_cast<int>(static_cast<std::uint64_t>(static_cast<std::int64_t>(count)) & 3U);

      /**
       *  Half an ulp of the radians: their power of two, 2^-53 of it.
       */
      const double halfUlp =
          std::bit_cast<double>(std::bit_cast<std::uint64_t>(radians) & 0x7ff0000000000000U)
          * 0x1p-53;
      if (radians < Snapping && size < halfUlp)
      {
        return Folded{.Quarter = turn, .High = 0.0, .Low = 0.0};
      }
      if (size < 0x1p-30)
      {
        return std::nullopt;
      }
      return Folded{.Quarter = turn, .High = rest.High, .Low = rest.Low};
    }

    /**
     *  Radians as degrees, a pair to a pair, within 2^-104.
     *
     *  @param  radians The radians, as a pair.
     *  @return         The degrees, as a pair.
     */
    auto InDegrees(Twofold::Two radians) -> Twofold::Two
    {
      const Twofold::Two lead = TwoProduct(radians.High, Table::ToDegrees.High);
      return TwoSum(lead.High, lead.Low
                                   + (radians.High * Table::ToDegrees.Low)
                                   + (radians.Low * Table::ToDegrees.High));
    }

    /**
     *  The cosine and sine of radians within a quarter turn and a hair of 0, as a pair: the nearest
     *  eighth of a degree, pi / 1440, taken off and the table's core given the rest.
     *
     *  The step's first product is taken off exactly, by Sterbenz's lemma, as the radians are
     *  within half a step of it; the step as a pair is within 2^-110 of pi / 1440, the rest within
     *  2^-100 of the angle the table's entries are added to.
     *
     *  @param  radians The radians, a pair, within pi / 4 and 2^-31 of 0.
     *  @return         The cosine and sine as pairs.
     */
    auto Tabled(Twofold::Two radians) -> Pairs
    {
      constexpr double   StepHigh  = Table::ToRadians.High / Eighths;
      constexpr double   StepLow   = Table::ToRadians.Low / Eighths;
      constexpr double   PerRadian = Table::ToDegrees.High * Eighths;
      const double       eighth    = std::rint(radians.High * PerRadian);
      const Twofold::Two step      = TwoProduct(eighth, StepHigh);
      const Twofold::Two part =
          TwoSum(radians.High - step.High, (radians.Low - step.Low) - (eighth * StepLow));
      return Core(static_cast<std::int64_t>(eighth), part.High, part.Low);
    }

    /**
     *  Positive radians of at least 3/4 reduced by quarter turns: Cody and Waite's way below 2^20,
     *  Payne and Hanek's past it and where the rest is tiny.
     *
     *  @param  radians The radians, finite, at least 3/4.
     *  @return         The quarter turns and the degrees left.
     */
    auto Fold(double radians) -> Folded
    {
      if (radians < 0x1p20)
      {
        if (const std::optional<Folded> near = Waite(radians))
        {
          const Twofold::Two degrees = InDegrees({.High = near->High, .Low = near->Low});
          return {.Quarter = near->Quarter, .High = degrees.High, .Low = degrees.Low};
        }
      }
      return Hanek(radians);
    }

    /**
     *  A pair the other way.
     *
     *  @param  pair The pair.
     *  @return      Its negative, exactly.
     */
    auto Negated(Twofold::Two pair) -> Twofold::Two
    {
      return {.High = -pair.High, .Low = -pair.Low};
    }

    /**
     *  A cosine and a sine turned by whole quarter turns, exactly.
     *
     *  @param  pairs   The cosine and sine.
     *  @param  quarter Counterclockwise quarter turns, 0 to 3.
     *  @return         The pairs turned.
     */
    auto Quarterly(Pairs pairs, int quarter) -> Pairs
    {
      /**
       *  Selects, not branches: the quarter of a random angle is as good as random.
       */
      const auto         turn   = static_cast<unsigned>(quarter);
      const bool         odd    = (turn & 1U) != 0U;
      const bool         left   = ((turn + 1U) & 2U) != 0U; /* quarters 1 and 2 */
      const bool         below  = (turn & 2U) != 0U;        /* quarters 2 and 3 */
      const Twofold::Two along  = odd ? pairs.Y : pairs.X;
      const Twofold::Two upward = odd ? pairs.X : pairs.Y;
      return {.X = left ? Negated(along) : along, .Y = below ? Negated(upward) : upward};
    }

    /**
     *  The cosine and sine of positive radians of at least 3/4 and below 2^20, reduced in one step
     *  to the nearest eighth of a degree, Cody and Waite's way, when they are not too near a
     *  quarter turn.
     *
     *  n = rint(radians 1440 / pi) is below 2^29, and pi / 1440 is held in three doubles to 2^-172
     *  of itself. The first product is taken off by one fused multiply-add, exactly: radians and n
     *  s1 are multiples of 2^-61, s1's ulp, and their difference is below 2^-9. The second product
     *  is exact as a pair; the third and the truncation leave below 2^-142. Off a quarter turn the
     *  table's entry carries the result, and that is far below its 2^-104; at one, where the rest
     *  is the result, it is 2^-112 of a rest of at least 2^-30, and Hanek decides below that.
     *
     *  @param  radians The radians, at least 3/4 and below 2^20.
     *  @return         The cosine and sine as pairs, or none within 2^-30 of a quarter turn.
     */
    auto Middling(double radians) -> std::optional<Pairs>
    {
      const double       count  = std::rint(radians * EighthsPerRadian);
      const double       first  = std::fma(-count, EighthFirst, radians);
      const Twofold::Two second = TwoProduct(count, EighthSecond);
      const Twofold::Two head   = TwoSum(first, -second.High);
      const double       tail   = (head.Low - second.Low) - (count * EighthThird);
      const Twofold::Two rest   = TwoSum(head.High, tail);
      const auto         whole  = static_cast<std::int64_t>(count);
      const std::int64_t turns  = (whole + 360) / 720;
      const std::int64_t index  = whole - (720 * turns);
      const int          turn   = static_cast<int>(static_cast<std::uint64_t>(turns) & 3U);
      if (index == 0)
      {
        /**
         *  Half an ulp of the radians: their power of two, 2^-53 of it.
         */
        const double halfUlp =
            std::bit_cast<double>(std::bit_cast<std::uint64_t>(radians) & 0x7ff0000000000000U)
            * 0x1p-53;
        const double size = std::abs(rest.High);
        if (radians < Snapping && size < halfUlp)
        {
          return Quarterly({.X = {.High = 1.0}, .Y = {}}, turn);
        }
        if (size < 0x1p-30)
        {
          return std::nullopt;
        }
      }
      return Quarterly(Core(index, rest.High, rest.Low), turn);
    }

    /**
     *  The nearest doubles to the cosine and sine of radians; at the double nearest a multiple of a
     *  quarter turn, below 16 radians, that quarter turn's.
     *
     *  @param  radians The radians.
     *  @return         cos and sin; a zero component +0.
     */
    auto Radial(double radians) -> Vector
    {
      if (!std::isfinite(radians))
      {
        constexpr double None = std::numeric_limits<double>::quiet_NaN();
        return {.X = None, .Y = None};
      }
      const double size = std::abs(radians);
      if (size < 0x1p-500)
      {
        /**
         *  cos r = 1 - r^2 / 2 and sin r = r - r^3 / 6, r^2 below 2^-1000: 1 and r, nearest.
         */
        return {.X = 1.0, .Y = radians + 0.0};
      }
      Pairs way;
      if (size < 0.75)
      {
        way = Tabled({.High = size, .Low = 0.0});
      }
      else if (const std::optional<Pairs> near =
                   size < 0x1p20 ? Middling(size) : std::optional<Pairs>{})
      {
        way = *near;
      }
      else
      {
        const Folded reduced = Hanek(size);
        way = Quarterly(Reduced<true>(reduced.High, reduced.Low).Turned, reduced.Quarter);
      }
      way.Y = radians < 0.0 ? Negated(way.Y) : way.Y;
      return Settle(way, [radians](Precise::Way pairs) -> std::array<double, 2>
                    { return Precise::NearestRadian(pairs, radians); });
    }
  } // namespace

  auto Direction(Angle angle) -> Vector
  {
    if (angle.HeldInDegrees())
    {
      return Cosined<false>(angle.InDegrees(), 0.0);
    }
    return Radial(angle.InRadians());
  }

  auto Direction(Angle start, Angle sweep, double along) -> Vector
  {
    if (!start.HeldInDegrees() || !sweep.HeldInDegrees())
    {
      return Direction(start + (sweep * along));
    }
    const Twofold::Two product = TwoProduct(sweep.InDegrees(), along);
    const Twofold::Two sum     = TwoSum(start.InDegrees(), product.High);
    const Twofold::Two angle   = TwoSum(sum.High, sum.Low + product.Low);
    return Cosined<true>(angle.High, angle.Low);
  }

  auto Turning::CosinePair(Twofold::Two degrees) -> std::array<Twofold::Two, 2>
  {
    if (std::abs(degrees.High) < Minute)
    {
      const Twofold::Two lead = TwoProduct(degrees.High, Table::ToRadians.High);
      return {Twofold::Two{.High = 1.0, .Low = 0.0},
              TwoSum(lead.High, lead.Low
                                    + (degrees.High * Table::ToRadians.Low)
                                    + (degrees.Low * Table::ToRadians.High))};
    }
    const Reduction reduced = ::Geoxy::Reduced<true>(degrees.High, degrees.Low);
    return {TwoSum(reduced.Turned.X.High, reduced.Turned.X.Low),
            TwoSum(reduced.Turned.Y.High, reduced.Turned.Y.Low)};
  }

  auto Turning::Reduced(double radians) -> Twofold::Two
  {
    if (std::abs(radians) < 0.75)
    {
      const double high = radians * Table::ToDegrees.High;
      return TwoSum(high, std::fma(radians, Table::ToDegrees.High, -high)
                              + (radians * Table::ToDegrees.Low));
    }
    const Folded       reduced = Fold(std::abs(radians));
    const double       turn    = radians < 0.0 ? -1.0 : 1.0;
    const double       base    = reduced.Quarter == 2 && reduced.High > 0.0
                                   ? -180.0
                                   : RightAngle * (reduced.Quarter == 3 ? -1.0 : reduced.Quarter);
    const Twofold::Two sum     = TwoSum(base, reduced.High);
    const Twofold::Two pair    = TwoSum(sum.High, sum.Low + reduced.Low);
    return {.High = turn * pair.High, .Low = turn * pair.Low};
  }
} // namespace Geoxy
