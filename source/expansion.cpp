/***************************************************************************************************
 *  @file      source/expansion.cpp
 *  @brief     Numbers held exactly as whole numbers of 64-bit limbs times a power of two, and the
 *             double nearest a quotient that may hold a square root.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "expansion.hpp"

#include "real.hpp"
#include "scratch.hpp"
#include "twofold.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <utility>

namespace Geoxy::Expansion
{
  namespace
  {
    using Twofold::Add;
    using Twofold::Difference;
    using Twofold::Divide;
    using Twofold::Double;
    using Twofold::Multiply;
    using Twofold::Root;

    /**
     *  The word a limb's sign fills a limb above it with: all ones below zero, else none.
     *
     *  @param  limb The largest limb.
     *  @return      The word.
     */
    auto SignWord(std::uint64_t limb) -> std::uint64_t
    {
      return std::uint64_t{0} - (limb >> 63U);
    }

    /**
     *  A number's limbs as a sum reads them: each at its place, zero below the least, the sign's
     *  word above the largest.
     */
    struct LimbRun
    {
      std::span<const std::uint64_t> Limbs = {}; /**< The limbs, from the least.            */
      std::int32_t                   Low   = 0;  /**< The least limb's place.               */
      std::uint64_t                  Above = 0;  /**< The word above the largest, its sign. */
    };

    /**
     *  A number's limbs lined up.
     *
     *  @param  number The number.
     *  @return        Its limbs.
     */
    auto RunOf(const Number &number) -> LimbRun
    {
      const std::span<const std::uint64_t> limbs = number.Parts.Limbs();
      return {.Limbs = limbs,
              .Low   = number.Parts.Low(),
              .Above = limbs.empty() ? 0 : SignWord(limbs.back())};
    }

    /**
     *  A number's limb at a place.
     *
     *  @param  lined The number's limbs.
     *  @param  place The place.
     *  @return       The limb there.
     */
    auto LimbAt(const LimbRun &lined, std::int32_t place) -> std::uint64_t
    {
      const std::int64_t index = static_cast<std::int64_t>(place) - lined.Low;
      if (index < 0)
      {
        return 0;
      }
      const auto offset = static_cast<std::size_t>(index);
      return offset < lined.Limbs.size() ? lined.Limbs[offset] : lined.Above;
    }

    /**
     *  A sum or a difference of two numbers, limb by limb with a carry, over the places either
     *  holds and one more for the carry: two's complement wraps the rest away.
     *
     *  @param  first  One number.
     *  @param  second The other.
     *  @param  minus  Whether the second is taken from the first: its limbs inverted, one carried
     *                 in.
     *  @return        The sum or the difference.
     */
    auto Totalled(const Number &first, const Number &second, bool minus) -> Number
    {
      const LimbRun one   = RunOf(first);
      const LimbRun other = RunOf(second);
      const auto    ends  = [](const LimbRun &lined) -> std::int32_t
      {
        return lined.Low + static_cast<std::int32_t>(lined.Limbs.size());
      };
      const std::int32_t             low  = std::min(one.Low, other.Low);
      const std::int32_t             high = std::max(ends(one), ends(other)) + 1;
      Number                         out;
      const std::span<std::uint64_t> room  = out.Parts.Room(static_cast<std::size_t>(high - low));
      const std::uint64_t            flip  = minus ? ~std::uint64_t{0} : 0;
      std::uint64_t                  carry = flip & 1U;
      for (std::int32_t place = low; place < high; ++place)
      {
        const std::uint64_t mine                    = LimbAt(one, place);
        const std::uint64_t theirs                  = LimbAt(other, place) ^ flip;
        const std::uint64_t partial                 = mine + theirs;
        const std::uint64_t total                   = partial + carry;
        carry                                       = static_cast<std::uint64_t>(partial < mine)
                                                    + static_cast<std::uint64_t>(total < partial);
        room[static_cast<std::size_t>(place - low)] = total;
      }
      out.Parts.Settle(low);
      return out;
    }

    /**
     *  The product of two limbs.
     */
    struct Doubleword
    {
      std::uint64_t High = 0; /**< The upper 64 bits. */
      std::uint64_t Low  = 0; /**< The lower 64 bits. */
    };

    /**
     *  The product of two limbs: one instruction where the compiler has 128-bit integers, four
     *  products of halves elsewhere.
     *
     *  @param  first  One limb.
     *  @param  second The other.
     *  @return        Their product.
     */
    auto WordProduct(std::uint64_t first, std::uint64_t second) -> Doubleword
    {
#ifdef __SIZEOF_INT128__
      __extension__ using Unsigned = unsigned __int128;

      const Unsigned product = static_cast<Unsigned>(first) * second;
      return {.High = static_cast<std::uint64_t>(product >> 64U),
              .Low  = static_cast<std::uint64_t>(product)};
#else
      constexpr std::uint64_t Half = 0xFFFF'FFFFU;

      const std::uint64_t lowLow   = (first & Half) * (second & Half);
      const std::uint64_t lowHigh  = (first & Half) * (second >> 32U);
      const std::uint64_t highLow  = (first >> 32U) * (second & Half);
      const std::uint64_t highHigh = (first >> 32U) * (second >> 32U);
      const std::uint64_t middle   = (lowLow >> 32U) + (lowHigh & Half) + (highLow & Half);
      return {.High = highHigh + (lowHigh >> 32U) + (highLow >> 32U) + (middle >> 32U),
              .Low  = (middle << 32U) + (lowLow & Half)};
#endif
    }

    /**
     *  Limbs negated in place, two's complement: inverted, one added.
     *
     *  @param limbs The limbs.
     */
    void Complemented(std::span<std::uint64_t> limbs)
    {
      std::uint64_t carry = 1;
      for (std::uint64_t &limb : limbs)
      {
        const std::uint64_t total = ~limb + carry;
        carry                     = static_cast<std::uint64_t>(total < carry);
        limb                      = total;
      }
    }

    /**
     *  Limbs less others, in place, with a borrow; what is borrowed past the last dropped, as two's
     *  complement wraps it.
     *
     *  @param limbs The limbs, at least as many as the others.
     *  @param taken The limbs taken from them.
     */
    void Lessened(std::span<std::uint64_t> limbs, std::span<const std::uint64_t> taken)
    {
      std::uint64_t borrow = 0;
      for (std::size_t index = 0; index < taken.size(); ++index)
      {
        const std::uint64_t limb    = limbs[index];
        const std::uint64_t partial = limb - taken[index];
        const std::uint64_t total   = partial - borrow;
        borrow                      = static_cast<std::uint64_t>(limb < taken[index])
                                    + static_cast<std::uint64_t>(partial < borrow);
        limbs[index]                = total;
      }
    }

    /**
     *  A number's magnitude read bit by bit: a negative number's limbs inverted and one added,
     *  which carries out of the least limb alone, as no number's least limb is zero.
     */
    struct Reading
    {
      std::span<const std::uint64_t> Limbs = {}; /**< The number's limbs.             */
      std::uint64_t                  Mask  = 0;  /**< All ones for a negative number. */

      /**
       *  A limb of the magnitude, none past the number's.
       *
       *  @param  index Which, from the least.
       *  @return       The limb.
       */
      [[nodiscard]] auto At(std::size_t index) const -> std::uint64_t
      {
        if (index >= Limbs.size())
        {
          return 0;
        }
        return (Limbs[index] ^ Mask) + (index == 0 ? (Mask & 1U) : 0U);
      }

      /**
       *  The 64 bits of the magnitude from a bit up, none below the least limb nor past the
       *  largest.
       *
       *  @param  bit The bit, from the least limb's least.
       *  @return     The bits.
       */
      [[nodiscard]] auto Window(std::int64_t bit) const -> std::uint64_t
      {
        if (bit <= -64)
        {
          return 0;
        }
        if (bit < 0)
        {
          return At(0) << static_cast<unsigned>(-bit);
        }
        const auto          limb  = static_cast<std::size_t>(bit / 64);
        const auto          shift = static_cast<unsigned>(bit % 64);
        const std::uint64_t upper = shift == 0 ? 0 : At(limb + 1) << (64U - shift);
        return (At(limb) >> shift) | upper;
      }

      /**
       *  Whether any bit of the magnitude lies below a bit: the least limb, never zero, does below
       *  any bit past it.
       *
       *  @param  bit The bit, from the least limb's least.
       *  @return     Whether one does.
       */
      [[nodiscard]] auto Under(std::int64_t bit) const -> bool
      {
        if (bit <= 0)
        {
          return false;
        }
        if (bit >= 64)
        {
          return true;
        }
        return (At(0) & ((std::uint64_t{1} << static_cast<unsigned>(bit)) - 1U)) != 0;
      }
    };

    /**
     *  A double from a significand and the place of its least bit, an infinity past the largest.
     *
     *  @param  significand The significand: below 2^53, at least 2^52 unless subnormal.
     *  @param  place       Its least bit's place, -1074 for a subnormal.
     *  @return             The double, positive.
     */
    auto Composed(std::uint64_t significand, std::int64_t place) -> double
    {
      constexpr std::uint64_t Fraction = (std::uint64_t{1} << 52U) - 1U;
      if ((significand >> 52U) == 0)
      {
        return std::bit_cast<double>(significand);
      }
      const std::int64_t field = place + 1075;
      if (field >= 2047)
      {
        return std::numeric_limits<double>::infinity();
      }
      return std::bit_cast<double>((static_cast<std::uint64_t>(field) << 52U)
                                   | (significand & Fraction));
    }

    /**
     *  The double nearest a number, ties to even, and the double nearest what it leaves within
     *  2^-116 of the number: the 53 bits from the leading one down, fewer below the normal doubles,
     *  the bit below them and whether any lies below that; and the 64 bits below them for the rest,
     *  less a whole step where the double was rounded up.
     *
     *  @param  number The number.
     *  @return        The two doubles, the second zero where the first is all the number or an
     *                 infinity.
     */
    auto Rounding(const Number &number) -> Double
    {
      const std::span<const std::uint64_t> limbs = number.Parts.Limbs();
      if (limbs.empty())
      {
        return {};
      }
      const Reading reading = {.Limbs = limbs, .Mask = SignWord(limbs.back())};
      std::size_t   top     = limbs.size() - 1;
      while (reading.At(top) == 0)
      {
        --top;
      }
      const std::int64_t low = 64 * static_cast<std::int64_t>(number.Parts.Low());
      const std::int64_t first =
          low + (64 * static_cast<std::int64_t>(top)) + 63 - std::countl_zero(reading.At(top));
      const std::int64_t least       = std::max<std::int64_t>(first - 52, -1074);
      const std::int64_t shift       = least - low;
      std::uint64_t      significand = reading.Window(shift) & ((std::uint64_t{1} << 53U) - 1U);
      const bool         half        = (reading.Window(shift - 1) & 1U) != 0;
      const bool         upward  = half && (reading.Under(shift - 1) || (significand & 1U) != 0);
      significand               += upward ? 1U : 0U;
      std::int64_t place         = least;
      if ((significand >> 53U) != 0)
      {
        significand >>= 1U;
        ++place;
      }
      const double size = Composed(significand, place);
      const double sign = reading.Mask != 0 ? -1.0 : 1.0;
      if (!std::isfinite(size))
      {
        return {.High = sign * size, .Low = 0.0};
      }
      const std::uint64_t window = reading.Window(shift - 64);
      const double        left =
          upward ? -(static_cast<double>(~window) + 1.0) : static_cast<double>(window);
      return {.High = sign * size, .Low = sign * std::ldexp(left, static_cast<int>(least - 64))};
    }

    /**
     *  Where two numbers' halfway points lie: the double a step above a double, as an ordered whole
     *  number, so that neighbouring doubles are neighbouring numbers.
     *
     *  @param  value The double.
     *  @return       Its place in the order of doubles.
     */
    auto Order(double value) -> std::int64_t
    {
      const auto bits = std::bit_cast<std::uint64_t>(value);
      const auto size = static_cast<std::int64_t>(bits & ~(std::uint64_t{1} << 63U));
      return (bits >> 63U) == 0 ? size : -size;
    }

    /**
     *  The double at a place in the order of doubles.
     *
     *  @param  place The place.
     *  @return       The double.
     */
    auto Unordered(std::int64_t place) -> double
    {
      if (place >= 0)
      {
        return std::bit_cast<double>(place);
      }
      return std::bit_cast<double>(static_cast<std::uint64_t>(-place) | (std::uint64_t{1} << 63U));
    }

    /**
     *  A quotient made ready to compare twice it against many numbers: what does not change between
     *  them, worked out once. Twice the quotient is compared, so that the halfway points between
     *  doubles are sums of two doubles, exactly.
     */
    struct Ready
    {
      const Quotient *Parts  = nullptr; /**< The quotient.                            */
      int             Side   = 1;       /**< Which root.                              */
      Number          Twice  = {};      /**< 2 top.                                   */
      Number          Rooted = {};      /**< (2 bottom)^2 root, for squared compares. */
      int             Root   = 0;       /**< The root's sign.                         */
      int             Other  = 0;       /**< The root's share's sign.                 */
      int             Below  = 0;       /**< The divisor's sign.                      */
    };

    /**
     *  A quotient made ready to compare.
     *
     *  @param  parts The quotient.
     *  @param  side  Which root.
     *  @return       It, ready.
     */
    auto Readied(const Quotient &parts, int side) -> Ready
    {
      return {.Parts  = &parts,
              .Side   = side,
              .Twice  = Scaled(parts.Top, 2.0),
              .Rooted = Scaled(Product(Product(parts.Bottom, parts.Bottom), parts.Root), 4.0),
              .Root   = Sign(parts.Root),
              .Other  = side * Sign(parts.Bottom),
              .Below  = Sign(parts.Below)};
    }

    /**
     *  The sign of twice (top + side bottom sqrt(root)) / below less a number, exactly.
     *
     *  @param  ready The quotient, ready.
     *  @param  mark  The number.
     *  @return       1, -1 or 0.
     */
    auto Compared(const Ready &ready, const Number &mark) -> int
    {
      const Number rest  = Difference(ready.Twice, Product(mark, ready.Parts->Below));
      const int    first = Sign(rest);
      int          sign  = first;
      if (ready.Root != 0 && ready.Other != 0)
      {
        if (first == 0 || first == ready.Other)
        {
          sign = ready.Other;
        }
        else
        {
          const int larger = Sign(Difference(Product(rest, rest), ready.Rooted));
          sign             = larger > 0 ? first : 0;
          if (larger < 0)
          {
            sign = ready.Other;
          }
        }
      }
      return sign * ready.Below;
    }

    /**
     *  Twice the halfway point between the double at a place and the next above: their sum,
     *  exactly. Past the largest double either way, where one of the two is an infinity, 2^1024
     *  stands for it, the halfway point rounding goes by: twice it is twice the largest double and
     *  2^971.
     *
     *  @param  place The place: of a double, or of an infinity.
     *  @return       Twice the halfway point.
     */
    auto Halfway(std::int64_t place) -> Number
    {
      const double here = Unordered(place);
      const double next = Unordered(place + 1);
      if (std::isfinite(here) && std::isfinite(next))
      {
        return Sum(Of(here), Of(next));
      }
      Number past = Sum(Scaled(Of(std::numeric_limits<double>::max()), 2.0), Of(0x1p971));
      if (!std::isfinite(next))
      {
        return past;
      }
      return Difference(Number{}, past);
    }

    /**
     *  A quotient worked out in double-double arithmetic: its value, within 2^-96 of its terms'
     *  size, and that bound.
     */
    struct Estimate
    {
      double Value = 0.0; /**< The value's leading part.     */
      double Low   = 0.0; /**< What the leading part leaves. */
      double Bound = 0.0; /**< How far the value may be off. */
    };

    /**
     *  A quotient worked out in double-double arithmetic, each of its numbers read as the double
     *  nearest it and the double nearest what that leaves.
     *
     *  @param  parts The quotient.
     *  @param  side  Which root.
     *  @return       The estimate; its bound infinite where it does not hold.
     */
    auto Estimated(const Quotient &parts, int side) -> Estimate
    {
      const Double top      = Rounding(parts.Top);
      const Double radicand = Rounding(parts.Root);
      const Double bottom   = Rounding(parts.Bottom);
      const Double share    = Multiply(Multiply(bottom, Root(radicand)),
                                       {.High = static_cast<double>(side), .Low = 0.0});
      const Double below    = Rounding(parts.Below);
      /**
       *  A number no double holds, below the least subnormal's half: no bound, so that every double
       *  is in question.
       */
      const auto vanished = [](const Number &number, const Double &read) -> bool
      {
        return !number.Parts.Empty() && Zero(read.High);
      };
      if (vanished(parts.Top, top)
          || vanished(parts.Root, radicand)
          || vanished(parts.Bottom, bottom)
          || vanished(parts.Below, below))
      {
        return {.Value = 0.0, .Low = 0.0, .Bound = std::numeric_limits<double>::infinity()};
      }
      /**
       *  Where the top and the root's share nearly cancel, as for the nearer root of a quadratic
       *  whose terms are large, their sum is worked as its conjugate's quotient: (t + s) = (t^2 -
       *  s^2) / (t - s), the numerator exact, as s^2 = bottom^2 root, the denominator free of
       *  cancellation; so the estimate keeps its digits however small the value: each part within
       *  2^-100 or so, the whole within 2^-88 of the value, where it is finite and well above the
       *  subnormals: a value past the doubles' range, or lost below it, leaves the estimate to the
       *  plain way. A numerator of zero is a value of zero, exactly.
       */
      const bool cancels =
          (top.High > 0.0 && share.High < 0.0) || (top.High < 0.0 && share.High > 0.0);
      if (cancels)
      {
        const Number gap    = Difference(Product(parts.Top, parts.Top),
                                         Product(Product(parts.Bottom, parts.Bottom), parts.Root));
        const double larger = std::max(std::abs(top.High), std::abs(share.High));
        const bool   inside = larger > 0x1p-300 && larger < 0x1p300;
        if (inside && Sign(gap) == 0)
        {
          /**
           *  t^2 = s^2 with the two of opposite signs, their squares well inside the doubles'
           *  range: t + s is zero, exactly.
           */
          return {.Value = 0.0, .Low = 0.0, .Bound = 0.0};
        }
        /**
         *  A numerator or a divisor near the subnormals, as for a line far shorter than the square
         *  root of the least double, or terms whose squares fall there though the value does not,
         *  where a rounding is a step of the least double and no share of its size, leaves the
         *  estimate to the plain way too; so does a value past the doubles' range or lost below.
         */
        const Double rest    = Rounding(gap);
        const Double divisor = Multiply(Difference(top, share), below);
        const bool   clear   = std::abs(rest.High) > 0x1p-900 && std::abs(divisor.High) > 0x1p-900;
        const Double value   = clear ? Divide(rest, divisor) : Double{};
        const double bound   = 0x1p-88 * std::abs(value.High);
        const double size    = std::abs(value.High);
        if (std::isfinite(size) && size > 0x1p-900)
        {
          return {.Value = value.High, .Low = value.Low, .Bound = bound};
        }
      }
      const Double value = Divide(Add(top, share), below);
      const double scale = (std::abs(top.High) + std::abs(share.High)) / std::abs(below.High);
      const double bound = 0x1p-96 * (scale + std::abs(value.High));
      if (!std::isfinite(value.High) || !std::isfinite(bound))
      {
        return {.Value = (top.High + share.High) / below.High,
                .Low   = 0.0,
                .Bound = std::numeric_limits<double>::infinity()};
      }
      /**
       *  A term or the value near the subnormals, where a product's rounding is a step of the least
       *  double and no share of its size: no bound, so that every double is in question.
       */
      double least = std::numeric_limits<double>::infinity();
      for (const double term :
           {top.High, radicand.High, bottom.High, share.High, below.High, value.High})
      {
        least = Zero(term) ? least : std::min(least, std::abs(term));
      }
      if (least < 0x1p-900)
      {
        return {.Value = value.High,
                .Low   = value.Low,
                .Bound = std::numeric_limits<double>::infinity()};
      }
      return {.Value = value.High, .Low = value.Low, .Bound = bound};
    }

    /**
     *  The double nearest an estimate, where no halfway point between doubles lies within its
     *  bound, the nearer halfway point taken as half the smaller step either side, as at a power of
     *  two the step below is half the step above; else nothing.
     *
     *  @param  estimate The estimate.
     *  @return          The double nearest, where settled.
     */
    auto Settled(const Estimate &estimate) -> std::optional<double>
    {
      const double value = estimate.Value;
      const double gap   = Spacing(value);
      if (gap - (2.0 * std::abs(estimate.Low)) > 2.0 * estimate.Bound)
      {
        return estimate.Value;
      }
      return std::nullopt;
    }
  } // namespace

  Parts::Parts(const Parts &other) : m_Low(other.m_Low), m_Count(other.m_Count)
  {
    const std::span<const std::uint64_t> limbs = other.Limbs();
    if (m_Count <= Held)
    {
      std::ranges::copy(limbs, m_Held.begin());
    }
    else
    {
      m_Heap.assign(limbs.begin(), limbs.end());
    }
  }

  Parts::Parts(Parts &&other) noexcept : m_Low(other.m_Low), m_Count(other.m_Count)
  {
    if (m_Count <= Held)
    {
      std::ranges::copy(other.Limbs(), m_Held.begin());
    }
    else
    {
      m_Heap        = std::move(other.m_Heap);
      other.m_Count = 0;
    }
  }

  auto Parts::operator=(const Parts &other) -> Parts &
  {
    if (this != &other)
    {
      const std::span<const std::uint64_t> limbs = other.Limbs();
      if (other.m_Count <= Held)
      {
        std::ranges::copy(limbs, m_Held.begin());
      }
      else
      {
        m_Heap.assign(limbs.begin(), limbs.end());
      }
      m_Low   = other.m_Low;
      m_Count = other.m_Count;
    }
    return *this;
  }

  auto Parts::operator=(Parts &&other) noexcept -> Parts &
  {
    if (this != &other)
    {
      m_Low   = other.m_Low;
      m_Count = other.m_Count;
      if (m_Count <= Held)
      {
        std::ranges::copy(other.Limbs(), m_Held.begin());
      }
      else
      {
        m_Heap        = std::move(other.m_Heap);
        other.m_Count = 0;
      }
    }
    return *this;
  }

  auto Parts::Room(std::size_t count) -> std::span<std::uint64_t>
  {
    m_Count = static_cast<std::uint32_t>(count);
    if (count <= Held)
    {
      return std::span{m_Held}.first(count);
    }
    m_Heap.resize(count);
    return std::span{m_Heap};
  }

  void Parts::Settle(std::int32_t low)
  {
    const std::span<std::uint64_t> room =
        m_Count <= Held ? std::span{m_Held}.first(m_Count) : std::span{m_Heap}.first(m_Count);
    std::size_t start = 0;
    while (start < room.size() && room[start] == 0)
    {
      ++start;
    }
    std::size_t end = room.size();
    while (end - start > 1 && room[end - 1] == SignWord(room[end - 2]))
    {
      --end;
    }
    const std::span<std::uint64_t> kept = room.subspan(start, end - start);
    if (room.size() > Held && kept.size() <= Held)
    {
      std::ranges::copy(kept, m_Held.begin());
      m_Heap = {};
    }
    else if (start > 0)
    {
      std::ranges::copy(kept, room.begin());
    }
    m_Count = static_cast<std::uint32_t>(kept.size());
    m_Low   = kept.empty() ? 0 : low + static_cast<std::int32_t>(start);
  }

  auto Of(double value) -> Number
  {
    Number out;
    if (Zero(value))
    {
      return out;
    }
    const auto          bits        = std::bit_cast<std::uint64_t>(value);
    const auto          field       = static_cast<std::int32_t>((bits >> 52U) & 0x7FFU);
    const std::uint64_t fraction    = bits & ((std::uint64_t{1} << 52U) - 1U);
    const std::uint64_t significand = field == 0 ? fraction : fraction | (std::uint64_t{1} << 52U);
    const std::int32_t  place       = field == 0 ? -1074 : field - 1075;
    /**
     *  The least bit's place is 64 limb + shift: the limb floor(place / 64), from -17 at -1074.
     */
    const std::int32_t             limb  = ((place + 1088) / 64) - 17;
    const auto                     shift = static_cast<unsigned>(place - (64 * limb));
    const std::span<std::uint64_t> room  = out.Parts.Room(2);
    room[0]                              = significand << shift;
    room[1]                              = shift == 0 ? 0 : significand >> (64U - shift);
    if ((bits >> 63U) != 0)
    {
      Complemented(room);
    }
    out.Parts.Settle(limb);
    return out;
  }

  auto Sum(const Number &first, const Number &second) -> Number
  {
    if (first.Parts.Empty())
    {
      return second;
    }
    if (second.Parts.Empty())
    {
      return first;
    }
    return Totalled(first, second, false);
  }

  auto Difference(const Number &first, const Number &second) -> Number
  {
    if (second.Parts.Empty())
    {
      return first;
    }
    return Totalled(first, second, true);
  }

  auto Scaled(const Number &number, double factor) -> Number
  {
    return Product(number, Of(factor));
  }

  auto Product(const Number &first, const Number &second) -> Number
  {
    Number out;
    if (first.Parts.Empty() || second.Parts.Empty())
    {
      return out;
    }
    /**
     *  The limbs multiplied as unsigned, then each negative factor's sign weight taken back: a
     *  two's complement factor of n limbs is its unsigned reading less 2^(64 n) where it is
     *  negative, so the product is the unsigned one less the other factor at the first's top place
     *  for each, the product of the two weights past the room; the product needs no more room than
     *  both.
     */
    const std::span<const std::uint64_t> one   = first.Parts.Limbs();
    const std::span<const std::uint64_t> other = second.Parts.Limbs();
    const std::span<std::uint64_t>       room  = out.Parts.Room(one.size() + other.size());
    std::ranges::fill(room, 0U);
    for (std::size_t mine = 0; mine < one.size(); ++mine)
    {
      std::uint64_t carry = 0;
      for (std::size_t theirs = 0; theirs < other.size(); ++theirs)
      {
        const Doubleword    part  = WordProduct(one[mine], other[theirs]);
        std::uint64_t      &limb  = room[mine + theirs];
        const std::uint64_t low   = part.Low + limb;
        const std::uint64_t total = low + carry;
        carry                     = part.High
                                  + static_cast<std::uint64_t>(low < limb)
                                  + static_cast<std::uint64_t>(total < low);
        limb                      = total;
      }
      room[mine + other.size()] = carry;
    }
    if (SignWord(one.back()) != 0)
    {
      Lessened(room.subspan(one.size()), other);
    }
    if (SignWord(other.back()) != 0)
    {
      Lessened(room.subspan(other.size()), one);
    }
    out.Parts.Settle(first.Parts.Low() + second.Parts.Low());
    return out;
  }

  auto Sign(const Number &number) -> int
  {
    const std::span<const std::uint64_t> limbs = number.Parts.Limbs();
    if (limbs.empty())
    {
      return 0;
    }
    return SignWord(limbs.back()) != 0 ? -1 : 1;
  }

  auto Approximate(const Number &number) -> double
  {
    return Rounding(number).High;
  }

  auto Rounded(const Quotient &parts, int side) -> double
  {
    const Estimate estimate = Estimated(parts, side);
    const bool     bounded  = std::isfinite(estimate.Bound);
    if (bounded)
    {
      if (const std::optional<double> settled = Settled(estimate))
      {
        return *settled;
      }
    }
    const Ready ready = Readied(parts, side);
    if (bounded && Compared(ready, Scaled(Of(estimate.Value), 2.0)) == 0)
    {
      /**
       *  A crossing that falls on a double, as at 0 or at a vertex, often does: no search.
       */
      return estimate.Value;
    }
    /**
     *  The least place whose halfway point above lies at or above the number, halved down to one
     *  place, each step one exact comparison. The value lies within the estimate's bound, and the
     *  estimate within half a step of its leading double, either way: the place is between the
     *  doubles two steps past both, where the value is past the halfway point below and short of
     *  the one above. Where the estimate has no bound, as past the doubles' range or near the
     *  subnormals, every place is in question, the infinities' too, each standing for every number
     *  past the largest double's halfway point its way: some 64 comparisons. The halving takes the
     *  middle of two places without their difference, which may pass the largest whole number.
     */
    constexpr std::int64_t Lowest  = -0x7FF0'0000'0000'0000;
    constexpr std::int64_t Highest =  0x7FF0'0000'0000'0000;
    const auto             above   = [&ready](std::int64_t place) -> bool
    {
      return Compared(ready, Halfway(place)) <= 0;
    };
    std::int64_t low  = Lowest - 1;
    std::int64_t high = Highest;
    if (bounded)
    {
      low  = std::max(Order(estimate.Value - (2.0 * estimate.Bound)) - 2, Lowest - 1);
      high = std::min(Order(estimate.Value + (2.0 * estimate.Bound)) + 2, Highest);
    }
    while (low < high - 1)
    {
      const std::int64_t middle    = std::midpoint(low, high);
      (above(middle) ? high : low) = middle;
    }
    if (Compared(ready, Halfway(high)) == 0)
    {
      const double here = Unordered(high);
      const double next = Unordered(high + 1);
      const auto   even = [](double value) -> bool
      {
        return (std::bit_cast<std::uint64_t>(value) & 1U) == 0;
      };
      return even(here) ? here : next;
    }
    return Unordered(high);
  }
} // namespace Geoxy::Expansion
