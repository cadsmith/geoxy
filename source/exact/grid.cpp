/***************************************************************************************************
 *  @file      source/exact/grid.cpp
 *  @brief     Both tests' determinants in 256 bits of integer.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/grid.hpp"

#include "exact/triangle.hpp"

#include "geoxy/point.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>

namespace Geoxy::Exact
{
  namespace
  {
    /**
     *  The most bits the coordinates may span: each then below 2^61, the determinant below 2^252.
     */
    constexpr int Span = 61;

    /**
     *  The most bits the orientation test's coordinates may span in two words: each then below
     *  2^125, the determinant below 2^253.
     */
    constexpr int WideSpan = 125;

    /**
     *  A signed integer of 256 bits, two's complement, its lowest word first.
     */
    using Wide = std::array<std::uint64_t, 4>;

    /**
     *  A double's bits, read as a whole number of its lowest set bit's place.
     */
    struct Bits
    {
      std::uint64_t Significand = 0;     /**< The 53 bit significand, its hidden bit set. */
      int           Exponent    = 0;     /**< The place of the significand's lowest bit.  */
      bool          Negative    = false; /**< Whether the double is below zero.           */
    };

    /**
     *  A normal double's significand and place.
     *
     *  @param  value The double, normal.
     *  @return       Its bits: value = +-Significand 2^Exponent.
     */
    auto Read(double value) -> Bits
    {
      constexpr std::uint64_t Fraction = (std::uint64_t{1} << 52U) - 1;

      const auto bits   = std::bit_cast<std::uint64_t>(value);
      const auto biased = static_cast<int>((bits >> 52U) & 0x7FFU);
      return {.Significand = (bits & Fraction) + (std::uint64_t{1} << 52U),
              .Exponent    = biased - 1075,
              .Negative    = (bits >> 63U) != 0};
    }

    /**
     *  The lowest place set and the highest place of some doubles, each normal or zero.
     */
    struct Spread
    {
      int Lowest  = std::numeric_limits<int>::max(); /**< The lowest set bit's place. */
      int Highest = std::numeric_limits<int>::min(); /**< The highest bit's place.    */
    };

    /**
     *  A wide integer's negation, modulo 2^256.
     *
     *  @param  wide The integer.
     *  @return      -wide.
     */
    auto Negated(const Wide &wide) -> Wide
    {
      Wide          negated = {};
      std::uint64_t carry   = 1;
      for (std::size_t index = 0; index < wide.size(); ++index)
      {
        negated.at(index) = ~wide.at(index) + carry;
        carry             = static_cast<std::uint64_t>(carry != 0 && negated.at(index) == 0);
      }
      return negated;
    }

    /**
     *  A 128 bit magnitude as a signed wide integer.
     *
     *  @param  magnitude The magnitude.
     *  @param  negative  Whether it is taken below zero.
     *  @return           +-magnitude.
     */
    auto Signed(Product magnitude, bool negative) -> Wide
    {
      const Wide wide = {magnitude.Low, magnitude.High, 0, 0};
      return negative ? Negated(wide) : wide;
    }

    /**
     *  Adds one wide integer to another, modulo 2^256.
     *
     *  @param sum  The integer added to.
     *  @param term The integer added.
     */
    void Add(Wide &sum, const Wide &term)
    {
      std::uint64_t carry = 0;
      for (std::size_t index = 0; index < sum.size(); ++index)
      {
        const std::uint64_t partial = sum.at(index) + term.at(index);
        const std::uint64_t total   = partial + carry;
        carry         = static_cast<std::uint64_t>(partial < term.at(index) || total < partial);
        sum.at(index) = total;
      }
    }

    /**
     *  Whether a wide integer is below zero.
     *
     *  @param  wide The integer.
     *  @return      Its top bit.
     */
    auto Negative(const Wide &wide) -> bool
    {
      return (wide.back() >> 63U) != 0;
    }

    /**
     *  The magnitude of a wide integer below 2^128 in magnitude.
     *
     *  @param  wide The integer.
     *  @return      Its magnitude, two words.
     */
    auto Magnitude(Wide wide) -> Product
    {
      if (Negative(wide))
      {
        wide = Negated(wide);
      }
      return {.High = wide.at(1), .Low = wide.at(0)};
    }

    /**
     *  The product of two magnitudes below 2^128, exactly.
     *
     *  @param  first  One magnitude.
     *  @param  second The other.
     *  @return        Their product, as a nonnegative wide integer.
     */
    auto Times(Product first, Product second) -> Wide
    {
      const Product lowLow   = Multiply(first.Low, second.Low);
      const Product lowHigh  = Multiply(first.Low, second.High);
      const Product highLow  = Multiply(first.High, second.Low);
      const Product highHigh = Multiply(first.High, second.High);

      Wide product = {lowLow.Low, lowLow.High, highHigh.Low, highHigh.High};
      Add(product, {0, lowHigh.Low, lowHigh.High, 0});
      Add(product, {0, highLow.Low, highLow.High, 0});
      return product;
    }

    /**
     *  Doubles read once: each one's bits, zero for zero, and where their set bits lie.
     *
     *  @tparam CountT How many.
     */
    template <std::size_t CountT>
    struct Readings
    {
      std::array<Bits, CountT> Each  = {}; /**< Each double's bits.   */
      Spread                   Where = {}; /**< Their set bits' span. */
    };

    /**
     *  Reads doubles that are each normal or zero.
     *
     *  @tparam CountT How many.
     *  @param  values The doubles.
     *  @return        Their readings; none where one is not normal or zero.
     */
    template <std::size_t CountT>
    auto ReadAll(const std::array<double, CountT> &values) -> std::optional<Readings<CountT>>
    {
      Readings<CountT> readings;
      for (std::size_t index = 0; index < CountT; ++index)
      {
        const double        value  = values.at(index);
        const auto          raw    = std::bit_cast<std::uint64_t>(value);
        const std::uint64_t biased = (raw >> 52U) & 0x7FFU;
        if ((raw << 1U) == 0)
        {
          continue;
        }
        if (biased == 0 || biased == 0x7FFU)
        {
          return std::nullopt;
        }
        const Bits bits         = Read(value);
        readings.Each.at(index) = bits;
        readings.Where.Lowest =
            std::min(readings.Where.Lowest, bits.Exponent + std::countr_zero(bits.Significand));
        readings.Where.Highest = std::max(readings.Where.Highest, bits.Exponent + 52);
      }
      return readings;
    }

    /**
     *  A signed integer of 128 bits, two's complement: a coordinate or a difference of two, in grid
     *  units, for the orientation test.
     */
    using Pair = Product;

    /**
     *  A signed 128 bit integer negated or kept, modulo 2^128: ~x + 1 is -x.
     *
     *  @param  value  The integer.
     *  @param  negate Whether it is negated.
     *  @return        -value or value.
     */
    auto Negate(Pair value, bool negate) -> Pair
    {
      if (!negate)
      {
        return value;
      }
      return {.High = ~value.High + static_cast<std::uint64_t>(value.Low == 0),
              .Low  = ~value.Low + 1};
    }

    /**
     *  A double, read, as a whole number of grid units.
     *
     *  @pre           The double is a whole number of units below 2^125.
     *  @param  bits   The double's bits, zero for zero.
     *  @param  lowest The grid, a place at or below the double's lowest set bit.
     *  @return        value 2^-lowest.
     */
    auto Units(Bits bits, int lowest) -> Pair
    {
      if (bits.Significand == 0)
      {
        return {};
      }
      const int shift = bits.Exponent - lowest;
      Pair      units;
      if (shift < 0)
      {
        units.Low = bits.Significand >> static_cast<unsigned>(-shift);
      }
      else if (shift < 64)
      {
        units.Low  = bits.Significand << static_cast<unsigned>(shift);
        units.High = (bits.Significand >> 1U) >> static_cast<unsigned>(63 - shift);
      }
      else
      {
        units.High = bits.Significand << static_cast<unsigned>(shift - 64);
      }
      return Negate(units, bits.Negative);
    }

    /**
     *  The difference of two signed 128 bit integers, modulo 2^128.
     *
     *  @param  end   The integer arrived at.
     *  @param  start The integer left.
     *  @return       end - start.
     */
    auto Minus(Pair end, Pair start) -> Pair
    {
      return {.High = end.High - start.High - static_cast<std::uint64_t>(end.Low < start.Low),
              .Low  = end.Low - start.Low};
    }

    /**
     *  Whether a signed 128 bit integer is below zero.
     *
     *  @param  value The integer.
     *  @return       Its top bit.
     */
    auto Below(Pair value) -> bool
    {
      return (value.High >> 63U) != 0;
    }

    /**
     *  The magnitude of a signed 128 bit integer above -2^127.
     *
     *  @param  value The integer.
     *  @return       |value|.
     */
    auto Absolute(Pair value) -> Product
    {
      return Negate(value, Below(value));
    }

    /**
     *  How two nonnegative wide integers compare.
     *
     *  @param  first  One.
     *  @param  second The other.
     *  @return        1, 0 or -1 as the first is larger, the same or smaller.
     */
    auto Compare(const Wide &first, const Wide &second) -> int
    {
      for (std::size_t index = first.size(); index > 0; --index)
      {
        if (first.at(index - 1) != second.at(index - 1))
        {
          return first.at(index - 1) > second.at(index - 1) ? 1 : -1;
        }
      }
      return 0;
    }

    /**
     *  The sign of a product of two signed 128 bit integers.
     *
     *  @param  magnitude The product's magnitude.
     *  @param  first     One factor.
     *  @param  second    The other.
     *  @return           1, 0 or -1.
     */
    auto SignOf(const Wide &magnitude, Pair first, Pair second) -> int
    {
      if (std::ranges::all_of(magnitude, [](std::uint64_t word) -> bool { return word == 0; }))
      {
        return 0;
      }
      return Below(first) != Below(second) ? -1 : 1;
    }

    /**
     *  A corner's difference from the point, in grid units, for the circle test: each below 2^62.
     */
    struct Delta
    {
      std::int64_t X = 0; /**< Along the x axis. */
      std::int64_t Y = 0; /**< Along the y axis. */
    };

    /**
     *  A double, read, as a whole number of grid units in one word, for the circle test.
     *
     *  @pre           The double is a whole number of units below 2^61.
     *  @param  bits   The double's bits, zero for zero.
     *  @param  lowest The grid, a place at or below the double's lowest set bit.
     *  @return        value 2^-lowest.
     */
    auto Count(Bits bits, int lowest) -> std::int64_t
    {
      if (bits.Significand == 0)
      {
        return 0;
      }
      const int           shift = bits.Exponent - lowest;
      const std::uint64_t whole = shift >= 0 ? bits.Significand << static_cast<unsigned>(shift)
                                             : bits.Significand >> static_cast<unsigned>(-shift);
      const auto          count = static_cast<std::int64_t>(whole);
      return bits.Negative ? -count : count;
    }

    /**
     *  A corner's term of the determinant: its lift times its minor, exactly.
     *
     *  @param  own    The corner's difference from the point, in grid units.
     *  @param  next   The corner after it.
     *  @param  before The corner before it.
     *  @return        (x^2 + y^2)(next.x before.y - next.y before.x).
     */
    auto Term(const Delta &own, const Delta &next, const Delta &before) -> Wide
    {
      const auto word = [](std::int64_t value) -> std::uint64_t
      {
        return value < 0 ? std::uint64_t{0} - static_cast<std::uint64_t>(value)
                         : static_cast<std::uint64_t>(value);
      };

      Wide lift = Signed(Multiply(word(own.X), word(own.X)), false);
      Add(lift, Signed(Multiply(word(own.Y), word(own.Y)), false));

      Wide minor = Signed(Multiply(word(next.X), word(before.Y)), (next.X < 0) != (before.Y < 0));
      Add(minor, Signed(Multiply(word(next.Y), word(before.X)), (next.Y < 0) == (before.X < 0)));

      const Wide product = Times(Magnitude(lift), Magnitude(minor));
      return Negative(minor) ? Negated(product) : product;
    }

    /**
     *  The sign of a wide integer.
     *
     *  @param  wide The integer.
     *  @return      1, 0 or -1.
     */
    auto SignOf(const Wide &wide) -> int
    {
      if (Negative(wide))
      {
        return -1;
      }
      return std::ranges::any_of(wide, [](std::uint64_t word) -> bool { return word != 0; }) ? 1
                                                                                             : 0;
    }
  } // namespace

  auto Multiply(std::uint64_t first, std::uint64_t second) -> Product
  {
#ifdef __SIZEOF_INT128__
    __extension__ using Unsigned = unsigned __int128;

    const Unsigned product = static_cast<Unsigned>(first) * second;
    return {.High = static_cast<std::uint64_t>(product >> 64U),
            .Low  = static_cast<std::uint64_t>(product)};
#else
    return MultiplyHalves(first, second);
#endif
  }

  /**
   *  Every coordinate is read as a significand and the place of its lowest bit; the lowest place of
   *  all is the grid. Where the highest bit lies at most Span places above it, each coordinate is a
   *  whole number of grid units below 2^61 and each difference one below 2^62, exact in a word.
   */
  auto OnGrid(const Triangle &circle, Point point) -> std::optional<int>
  {
    const std::array<double, 8> values = {circle.First.X,  circle.First.Y, circle.Second.X,
                                          circle.Second.Y, circle.Third.X, circle.Third.Y,
                                          point.X,         point.Y};

    const std::optional<Readings<8>> readings = ReadAll(values);
    if (!readings)
    {
      return std::nullopt;
    }
    const int lowest = readings->Where.Lowest;
    if (readings->Where.Highest < lowest)
    {
      return 0;
    }
    if (readings->Where.Highest - lowest >= Span)
    {
      return std::nullopt;
    }

    const auto count = [&readings, lowest](std::size_t index) -> std::int64_t
    {
      return Count(readings->Each.at(index), lowest);
    };
    const auto corner = [&count](std::size_t index) -> Delta
    {
      return {.X = count(index) - count(6), .Y = count(index + 1) - count(7)};
    };
    const Delta a = corner(0);
    const Delta b = corner(2);
    const Delta c = corner(4);

    Wide determinant = Term(a, b, c);
    Add(determinant, Term(b, c, a));
    Add(determinant, Term(c, a, b));
    return SignOf(determinant);
  }

  /**
   *  As OnGrid, with coordinates below 2^125 units in two words: differences below 2^126, products
   *  below 2^252 and the determinant below 2^253, exact in 256 bits.
   */
  auto OrientGrid(const Triangle &corners) -> std::optional<int>
  {
    const std::array<double, 6> values = {corners.First.X,  corners.First.Y, corners.Second.X,
                                          corners.Second.Y, corners.Third.X, corners.Third.Y};

    const std::optional<Readings<6>> readings = ReadAll(values);
    if (!readings)
    {
      return std::nullopt;
    }
    const int lowest = readings->Where.Lowest;
    if (readings->Where.Highest < lowest)
    {
      return 0;
    }
    if (readings->Where.Highest - lowest >= WideSpan)
    {
      return std::nullopt;
    }

    const auto units = [&readings, lowest](std::size_t index) -> Pair
    {
      return Units(readings->Each.at(index), lowest);
    };
    const Pair firstX  = Minus(units(0), units(4));
    const Pair firstY  = Minus(units(1), units(5));
    const Pair secondX = Minus(units(2), units(4));
    const Pair secondY = Minus(units(3), units(5));

    const Wide left      = Times(Absolute(firstX), Absolute(secondY));
    const Wide right     = Times(Absolute(firstY), Absolute(secondX));
    const int  leftSign  = SignOf(left, firstX, secondY);
    const int  rightSign = SignOf(right, firstY, secondX);
    if (leftSign != rightSign)
    {
      return leftSign > rightSign ? 1 : -1;
    }
    return leftSign * Compare(left, right);
  }
} // namespace Geoxy::Exact
