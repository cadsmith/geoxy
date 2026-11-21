/***************************************************************************************************
 *  @file      source/precise.cpp
 *  @brief     Angles to any precision: whole numbers in 32-bit words, pi by Machin's formula, sines
 *             and cosines by Taylor's series, arctangents by Euler's.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "precise.hpp"

#include "radians.hpp"
#include "scratch.hpp"
#include "turning.hpp"
#include "twofold.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>

namespace Geoxy::Precise
{
  namespace
  {
    using Words = List<std::uint32_t>;

    constexpr std::uint64_t Mask  = 0xffffffffU; /**< A word's 32 bits.                */
    constexpr int           Exact =        1074; /**< Bits that hold any double whole. */
    constexpr double Settling = 0x1p-68; /**< How far the fast pairs are, at most, relatively. */

    /**
     *  A whole number without its leading zero words.
     *
     *  @param  words The words, least significant first.
     *  @return       The same number, trimmed.
     */
    auto Trimmed(Words words) -> Words
    {
      while (!words.empty() && words.back() == 0U)
      {
        words.pop_back();
      }
      return words;
    }

    /**
     *  A word of a whole number, 0 past its last.
     *
     *  @param  words The words.
     *  @param  index The word's place.
     *  @return       The word.
     */
    auto At(const Words &words, std::size_t index) -> std::uint64_t
    {
      return index < words.size() ? words[index] : 0U;
    }

    /**
     *  The order of two trimmed whole numbers.
     *
     *  @param  lhs One.
     *  @param  rhs The other.
     *  @return     -1, 0 or 1.
     */
    auto Compare(const Words &lhs, const Words &rhs) -> int
    {
      if (lhs.size() != rhs.size())
      {
        return lhs.size() < rhs.size() ? -1 : 1;
      }
      for (std::size_t index = lhs.size(); index > 0; --index)
      {
        if (lhs[index - 1] != rhs[index - 1])
        {
          return lhs[index - 1] < rhs[index - 1] ? -1 : 1;
        }
      }
      return 0;
    }

    /**
     *  The sum of two whole numbers.
     *
     *  @param  lhs One.
     *  @param  rhs The other.
     *  @return     The sum, trimmed.
     */
    auto Add(const Words &lhs, const Words &rhs) -> Words
    {
      Words         sum(std::max(lhs.size(), rhs.size()) + 1, 0U);
      std::uint64_t carry = 0;
      for (std::size_t index = 0; index < sum.size(); ++index)
      {
        carry       += At(lhs, index) + At(rhs, index);
        sum[index]   = static_cast<std::uint32_t>(carry & Mask);
        carry      >>= 32U;
      }
      return Trimmed(sum);
    }

    /**
     *  The difference of two whole numbers, the first the larger.
     *
     *  @param  lhs The larger.
     *  @param  rhs The smaller.
     *  @return     The difference, trimmed.
     */
    auto Difference(const Words &lhs, const Words &rhs) -> Words
    {
      Words         difference(lhs.size(), 0U);
      std::uint64_t borrow = 0;
      for (std::size_t index = 0; index < lhs.size(); ++index)
      {
        const std::uint64_t value = (lhs[index] + (Mask + 1U)) - At(rhs, index) - borrow;
        difference[index]         = static_cast<std::uint32_t>(value & Mask);
        borrow                    = 1U - (value >> 32U);
      }
      return Trimmed(difference);
    }

    /**
     *  The product of two whole numbers, exactly.
     *
     *  @param  lhs One.
     *  @param  rhs The other.
     *  @return     The product, trimmed.
     */
    auto Multiply(const Words &lhs, const Words &rhs) -> Words
    {
      Words product(lhs.size() + rhs.size() + 1, 0U);
      for (std::size_t outer = 0; outer < lhs.size(); ++outer)
      {
        std::uint64_t carry = 0;
        for (std::size_t inner = 0; inner < rhs.size(); ++inner)
        {
          carry += (static_cast<std::uint64_t>(lhs[outer]) * rhs[inner]) + product[outer + inner];
          product[outer + inner]   = static_cast<std::uint32_t>(carry & Mask);
          carry                  >>= 32U;
        }
        product[outer + rhs.size()] = static_cast<std::uint32_t>(carry);
      }
      return Trimmed(product);
    }

    /**
     *  A whole number times 2^places.
     *
     *  @param  words  The number.
     *  @param  places The places, 0 or more.
     *  @return        The number shifted, trimmed.
     */
    auto Raised(const Words &words, int places) -> Words
    {
      const auto whole = static_cast<std::size_t>(places / 32);
      const auto part  = static_cast<unsigned>(places % 32);
      Words      raised(words.size() + whole + 1, 0U);
      for (std::size_t index = 0; index < words.size(); ++index)
      {
        const std::uint64_t value  = static_cast<std::uint64_t>(words[index]) << part;
        raised[index + whole]     |= static_cast<std::uint32_t>(value & Mask);
        raised[index + whole + 1] |= static_cast<std::uint32_t>(value >> 32U);
      }
      return Trimmed(raised);
    }

    /**
     *  A whole number over 2^places, truncated.
     *
     *  @param  words  The number.
     *  @param  places The places, 0 or more.
     *  @return        The number shifted, trimmed.
     */
    auto Lowered(const Words &words, int places) -> Words
    {
      const auto whole = static_cast<std::size_t>(places / 32);
      const auto part  = static_cast<unsigned>(places % 32);
      Words      lowered(words.size() - std::min(whole, words.size()), 0U);
      for (std::size_t index = 0; index < lowered.size(); ++index)
      {
        const std::uint64_t value =
            (At(words, index + whole) | (At(words, index + whole + 1) << 32U)) >> part;
        lowered[index] = static_cast<std::uint32_t>(value & Mask);
      }
      return Trimmed(lowered);
    }

    /**
     *  A whole number times a word.
     *
     *  @param  words  The number.
     *  @param  factor The word.
     *  @return        The product, trimmed.
     */
    auto Times(const Words &words, std::uint32_t factor) -> Words
    {
      Words         product(words.size() + 1, 0U);
      std::uint64_t carry = 0;
      for (std::size_t index = 0; index < words.size(); ++index)
      {
        carry           += static_cast<std::uint64_t>(words[index]) * factor;
        product[index]   = static_cast<std::uint32_t>(carry & Mask);
        carry          >>= 32U;
      }
      product[words.size()] = static_cast<std::uint32_t>(carry);
      return Trimmed(product);
    }

    /**
     *  A whole number over a word, truncated.
     *
     *  @param  words   The number.
     *  @param  divisor The word, above 0.
     *  @return         The quotient, trimmed.
     */
    auto Over(const Words &words, std::uint32_t divisor) -> Words
    {
      Words         quotient(words.size(), 0U);
      std::uint64_t remainder = 0;
      for (std::size_t index = words.size(); index > 0; --index)
      {
        remainder            = (remainder << 32U) | words[index - 1];
        quotient[index - 1]  = static_cast<std::uint32_t>(remainder / divisor);
        remainder           %= divisor;
      }
      return Trimmed(quotient);
    }

    /**
     *  A whole number over another, truncated, a bit at a time.
     *
     *  @param  dividend The number divided.
     *  @param  divisor  The number it is divided by, above 0.
     *  @return          The quotient, trimmed.
     */
    auto Divided(const Words &dividend, const Words &divisor) -> Words
    {
      Words quotient(dividend.size(), 0U);
      Words remainder(divisor.size() + 1, 0U);
      for (std::size_t bit = dividend.size() * 32; bit > 0; --bit)
      {
        /**
         *  remainder = 2 remainder + the next bit, in place; it stays below 2 divisor.
         */
        const std::size_t place = bit - 1;
        std::uint32_t     carry = (dividend[place / 32] >> (place % 32)) & 1U;
        for (std::uint32_t &word : remainder)
        {
          const std::uint32_t next = word >> 31U;
          word                     = (word << 1U) | carry;
          carry                    = next;
        }
        /**
         *  Take the divisor off where it fits, in place, the borrow running through every word.
         */
        bool fits = true;
        for (std::size_t index = remainder.size(); index > 0; --index)
        {
          const std::uint64_t has  = remainder[index - 1];
          const std::uint64_t less = At(divisor, index - 1);
          if (has != less)
          {
            fits = has > less;
            break;
          }
        }
        if (fits)
        {
          std::uint64_t borrow = 0;
          for (std::size_t index = 0; index < remainder.size(); ++index)
          {
            const std::uint64_t value =
                (remainder[index] + (Mask + 1U)) - At(divisor, index) - borrow;
            remainder[index] = static_cast<std::uint32_t>(value & Mask);
            borrow           = 1U - (value >> 32U);
          }
          quotient[place / 32] |= 1U << (place % 32);
        }
      }
      return Trimmed(quotient);
    }

    /**
     *  The sum of two signed numbers at one precision.
     *
     *  @param  lhs One.
     *  @param  rhs The other.
     *  @return     The sum; zero is not negative.
     */
    auto Sum(const Fixed &lhs, const Fixed &rhs) -> Fixed
    {
      if (lhs.Negative == rhs.Negative)
      {
        return {.Words = Add(lhs.Words, rhs.Words), .Negative = lhs.Negative};
      }
      const int  order  = Compare(lhs.Words, rhs.Words);
      const bool larger = order >= 0;
      Fixed sum = {.Words =
                       larger ? Difference(lhs.Words, rhs.Words) : Difference(rhs.Words, lhs.Words),
                   .Negative = larger ? lhs.Negative : rhs.Negative};
      sum.Negative = sum.Negative && !sum.Words.empty();
      return sum;
    }

    /**
     *  A signed number the other way.
     *
     *  @param  value The number.
     *  @return       Its negative.
     */
    auto Negated(Fixed value) -> Fixed
    {
      value.Negative = !value.Negative && !value.Words.empty();
      return value;
    }

    /**
     *  The product of two signed numbers at one precision, truncated.
     *
     *  @param  lhs  One.
     *  @param  rhs  The other.
     *  @param  bits The precision.
     *  @return      The product.
     */
    auto Product(const Fixed &lhs, const Fixed &rhs, int bits) -> Fixed
    {
      Fixed product    = {.Words    = Lowered(Multiply(lhs.Words, rhs.Words), bits),
                          .Negative = lhs.Negative != rhs.Negative};
      product.Negative = product.Negative && !product.Words.empty();
      return product;
    }

    /**
     *  A double at a precision, truncated toward zero: exact when bits holds its lowest bit.
     *
     *  @param  value The double, finite.
     *  @param  bits  The precision.
     *  @return       The number.
     */
    auto Of(double value, int bits) -> Fixed
    {
      int          exponent = 0;
      const double fraction = std::frexp(std::abs(value), &exponent);
      const auto   whole    = static_cast<std::uint64_t>(std::ldexp(fraction, 53));
      const Words  words    = Trimmed(
          {static_cast<std::uint32_t>(whole & Mask), static_cast<std::uint32_t>(whole >> 32U)});
      const int shift = exponent - 53 + bits;
      return {.Words    = shift >= 0 ? Raised(words, shift) : Lowered(words, -shift),
              .Negative = value < 0.0};
    }

    /**
     *  1 at a precision.
     *
     *  @param  bits The precision.
     *  @return      2^bits as words.
     */
    auto One(int bits) -> Words
    {
      return Raised({1U}, bits);
    }

    /**
     *  atan(1 / n) by its series, to within the terms it takes of 2^-bits.
     *
     *  @param  inverse n, from 2 to 65535.
     *  @param  bits    The precision.
     *  @return         The arctangent.
     */
    auto Inverse(std::uint32_t inverse, int bits) -> Fixed
    {
      Words         power = Over(One(bits), inverse);
      Fixed         sum;
      std::uint32_t odd      = 1;
      bool          negative = false;
      while (!power.empty())
      {
        sum       = Sum(sum, {.Words = Over(power, odd), .Negative = negative});
        power     = Over(power, inverse * inverse);
        odd      += 2;
        negative  = !negative;
      }
      return sum;
    }

  } // namespace

  auto Pi(int bits) -> Fixed
  {
    constexpr int Held = 4096;
    if (bits <= Held)
    {
      /**
       *  floor(pi 2^4096) truncated, within 1 unit.
       */
      return {.Words = Lowered(Words(Table::PiWords.begin(), Table::PiWords.end()), Held - bits)};
    }
    const Fixed fifth = Inverse(5, bits);
    const Fixed other = Inverse(239, bits);
    return Sum({.Words = Times(fifth.Words, 16)},
               {.Words = Times(other.Words, 4), .Negative = true});
  }

  namespace
  {

    /**
     *  The cosine and sine of radians within a quarter turn and a hair of 0, by Taylor's series,
     *  each term a third of the one before or less; within 2^(19 - bits) for bits to 2^16.
     *
     *  @param  angle The radians.
     *  @param  bits  The precision.
     *  @return       The cosine, then the sine.
     */
    auto Trigonometric(const Fixed &angle, int bits) -> List<Fixed>
    {
      const Fixed square = Product(angle, angle, bits);
      Fixed       cosine = {.Words = One(bits)};
      Fixed       sine   = angle;
      Fixed       even   = cosine;
      Fixed       odd    = angle;
      /**
       *  The odd term never outlasts the even: past the first, each is the even one's times |angle|
       *  / (2 step + 1) or less, below 1 for an angle within 3 radians.
       */
      for (std::uint32_t step = 1; !even.Words.empty(); ++step)
      {
        even       = Negated(Product(even, square, bits));
        even.Words = Over(even.Words, ((2 * step) - 1) * (2 * step));
        odd        = Negated(Product(odd, square, bits));
        odd.Words  = Over(odd.Words, (2 * step) * ((2 * step) + 1));
        cosine     = Sum(cosine, even);
        sine       = Sum(sine, odd);
      }
      return {cosine, sine};
    }

    /**
     *  A cosine and a sine turned by quarter turns, and the sine's sign flipped where asked.
     *
     *  @param  turned  The cosine, then the sine.
     *  @param  turns   Counterclockwise quarter turns, 0 to 3.
     *  @param  flipped Whether the sine changes sign after.
     *  @return         The cosine, then the sine, turned.
     */
    auto Quartered(const List<Fixed> &turned, int turns, bool flipped) -> List<Fixed>
    {
      const auto  turn   = static_cast<unsigned>(turns);
      const bool  odd    = (turn & 1U) != 0U;
      const Fixed along  = odd ? turned[1] : turned[0];
      const Fixed upward = odd ? turned[0] : turned[1];
      const bool  left   = ((turn + 1U) & 2U) != 0U;
      const bool  below  = ((turn & 2U) != 0U) != flipped;
      return {left ? Negated(along) : along, below ? Negated(upward) : upward};
    }

    /**
     *  The place of the lowest bit of any of some doubles, as a power of two: 2^-place is it.
     *
     *  @param  values The doubles.
     *  @return        The place, 0 when every one is whole or zero.
     */
    auto Lowest(std::span<const double> values) -> int
    {
      int lowest = 0;
      for (const double value : values)
      {
        int          exponent = 0;
        const double fraction = std::frexp(value, &exponent);
        const auto   whole    = static_cast<std::uint64_t>(std::abs(std::ldexp(fraction, 53)));
        /**
         *  A zero has 64 zeros and exponent 0, which places it at -11, below any other.
         */
        lowest = std::max(lowest, 53 - exponent - std::countr_zero(whole));
      }
      return lowest;
    }

    /**
     *  A sum of doubles as a signed whole number of 2^-place, exactly.
     *
     *  @param  values The doubles.
     *  @param  place  The place, at least each double's lowest bit's.
     *  @return        The number.
     */
    auto Whole(std::span<const double> values, int place) -> Fixed
    {
      Fixed sum;
      for (const double value : values)
      {
        sum = Sum(sum, Of(value, place));
      }
      return sum;
    }
  } // namespace

  auto Nearest(Twofold::Two pair, double error, const Computed &value) -> double
  {
    const double high = pair.High;
    const double low  = pair.Low;
    if (const std::optional<double> settled = Turning::Settled(high, low, error))
    {
      return *settled;
    }
    const double rounded = high + low;
    const double size    = std::abs(rounded);
    if (!(size >= 0x1p-1022))
    {
      /**
       *  On the subnormal grid: the value truncated to it, from far more bits than it holds, is the
       *  point below it or the one below that; the halfway point above settles which.
       */
      constexpr int Bits  = 2048;
      const Fixed   held  = value(Bits);
      const Words   steps = Lowered(held.Words, Bits - Exact);
      const double  sign  = held.Negative ? -1.0 : 1.0;
      const double  base =
          sign * std::ldexp(static_cast<double>(At(steps, 0) | (At(steps, 1) << 32U)), -Exact);
      const double step = sign * 0x1p-1074;
      return Side(value, base, step) * sign > 0 ? base + step : base;
    }
    /**
     *  Near a halfway point: the one the residual points to, past which the neighbour is nearer.
     */
    const double residual = (high - rounded) + low;
    const auto   bits     = std::bit_cast<std::uint64_t>(size);
    const double unit     = std::bit_cast<double>(bits & 0x7ff0000000000000U) * 0x1p-52;
    const bool   finer    = static_cast<bool>(static_cast<unsigned>((bits & 0xfffffffffffffU) == 0U)
                                              & static_cast<unsigned>(residual * rounded < 0.0)
                                              & static_cast<unsigned>(size > 0x1p-1022));
    const double step     = std::copysign(finer ? unit / 2.0 : unit, residual);
    return Side(value, rounded, step) * step > 0.0 ? rounded + step : rounded;
  }

  auto Slowed() -> std::uint64_t &
  {
    thread_local std::uint64_t count = 0;
    return count;
  }

  auto Side(const Computed &value, double base, double step) -> int
  {
    ++Slowed();
    int exponent = 0;
    static_cast<void>(std::frexp(step, &exponent));
    /**
     *  The halfway point, base + step / 2, is exact from 2 - exponent bits on; at the bits that
     *  first hold it the value's error is above it, so the loop always doubles at least once.
     */
    for (int bits = std::max(64, 4 - exponent);; bits *= 2)
    {
      const Fixed stride = Of(step, bits);
      const Fixed middle =
          Sum(Of(base, bits), {.Words = Lowered(stride.Words, 1), .Negative = stride.Negative});
      const Fixed distance = Sum(value(bits), Negated(middle));
      if (Compare(distance.Words, {1U << 21U}) > 0)
      {
        return distance.Negative ? -1 : 1;
      }
    }
  }

  auto Degree(Twofold::Two degrees, int quarter, int bits) -> List<Fixed>
  {
    /**
     *  Within 45 + 1/16 degrees: pi's error, 2^16 units, times a quarter, and Taylor's 2^19.
     */
    const Fixed held  = Sum(Of(degrees.High, bits), Of(degrees.Low, bits));
    Fixed       angle = Product(held, Pi(bits), bits);
    angle.Words       = Over(angle.Words, 180);
    return Quartered(Trigonometric(angle, bits), quarter, false);
  }

  auto Radian(double radians, int bits) -> List<Fixed>
  {
    /**
     *  Exactly: n = rint(r / (pi / 2)) and r - n pi / 2, pi to as many more bits as r has whole
     *  ones and 8 besides, so n's error times pi's stays below 2^-bits.
     */
    int exponent = 0;
    static_cast<void>(std::frexp(radians, &exponent));
    const int   extra   = std::max(0, exponent) + 8;
    const int   wide    = bits + extra;
    const Fixed size    = Of(std::abs(radians), wide);
    const Fixed quarter = {.Words = Lowered(Pi(wide).Words, 1)};
    /**
     *  Clamped, so the cast is defined where the guess goes unused.
     */
    const auto guess = static_cast<std::uint64_t>(
        std::rint(std::min(std::abs(radians), 0x1p50) * 0x1.45f306dc9c883p-1));
    const Words count = std::abs(radians) < 0x1p50
                          ? Trimmed({static_cast<std::uint32_t>(guess & Mask),
                                     static_cast<std::uint32_t>(guess >> 32U)})
                          : Divided(Add(size.Words, Lowered(quarter.Words, 1)), quarter.Words);
    const Fixed rest  = Sum(size, Negated({.Words = Multiply(count, quarter.Words)}));
    const Fixed angle = {.Words = Lowered(rest.Words, extra), .Negative = rest.Negative};
    /**
     *  The quarter turns, then the sign: cos is even, sin odd.
     */
    return Quartered(Trigonometric(angle, bits), static_cast<int>(At(count, 0) & 3U),
                     radians < 0.0);
  }

  auto Arc(std::span<const double> y, std::span<const double> x, int bits) -> Fixed
  {
    const int   place = std::max(Lowest(y), Lowest(x));
    const Fixed rise  = Whole(y, place);
    const Fixed run   = Whole(x, place);
    const bool  steep = Compare(rise.Words, run.Words) > 0;
    const Words over  = steep ? run.Words : rise.Words;
    const Words under = steep ? rise.Words : run.Words;

    /**
     *  Euler's series: atan t = u (1 + 2/3 w + 2 4/(3 5) w^2 + ...), t = over / under, u = t / (1 +
     *  t^2) and w = t^2 / (1 + t^2), at most 1/2, so each term is half the one before or less.
     */
    const Words square = Add(Multiply(over, over), Multiply(under, under));
    const Fixed weight = {.Words = Divided(Raised(Multiply(over, over), bits), square)};
    Fixed       term   = {.Words = Divided(Raised(Multiply(over, under), bits), square)};
    Fixed       arc    = term;
    for (std::uint32_t step = 1; !term.Words.empty(); ++step)
    {
      term       = Product(term, weight, bits);
      term.Words = Over(Times(term.Words, 2 * step), (2 * step) + 1);
      arc        = Sum(arc, term);
    }

    /**
     *  Unfold the octant: pi / 2 or pi, plus or minus the octant's angle, then the sign of y.
     */
    const Fixed half     = Pi(bits);
    const Fixed sideways = run.Negative ? half : Fixed{};
    const Fixed base     = steep ? Fixed{.Words = Lowered(half.Words, 1)} : sideways;
    const Fixed heading  = Sum(base, steep == run.Negative ? arc : Negated(arc));
    return rise.Negative ? Negated(heading) : heading;
  }

  auto NearestDegree(Way way, Twofold::Two degrees, int quarter) -> std::array<double, 2>
  {
    /**
     *  Captured by reference, so std::function holds them in place, not in memory it asks for.
     */
    return {Nearest(way.Cosine, Settling, [&degrees, quarter](int bits) -> Fixed
                    { return Degree(degrees, quarter, bits)[0]; }),
            Nearest(way.Sine, Settling, [&degrees, quarter](int bits) -> Fixed
                    { return Degree(degrees, quarter, bits)[1]; })};
  }

  auto NearestRadian(Way way, double radians) -> std::array<double, 2>
  {
    return {Nearest(way.Cosine, Settling,
                    [radians](int bits) -> Fixed { return Radian(radians, bits)[0]; }),
            Nearest(way.Sine, Settling,
                    [radians](int bits) -> Fixed { return Radian(radians, bits)[1]; })};
  }

  auto NearestSine(double near, Twofold::Two degrees) -> double
  {
    return Nearest({.High = near, .Low = 0.0}, 0x1p-50,
                   [degrees](int bits) -> Fixed { return Degree(degrees, 0, bits)[1]; });
  }

  auto NearestArc(Twofold::Two pair, std::span<const double> y, std::span<const double> x) -> double
  {
    return Nearest(pair, Settling, [&y, &x](int bits) -> Fixed { return Arc(y, x, bits); });
  }
} // namespace Geoxy::Precise
