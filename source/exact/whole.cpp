/***************************************************************************************************
 *  @file      source/exact/whole.cpp
 *  @brief     Signed whole numbers of a few thousand bits, the two determinants in them, and where
 *             two lines cross.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/whole.hpp"

#include "exact/grid.hpp"
#include "exact/triangle.hpp"

#include "geoxy/point.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>

namespace Geoxy::Exact
{
  namespace
  {
    /**
     *  The words a coordinate or a difference of two takes: below 2^2099.
     */
    constexpr std::size_t Coordinate = 33;

    /**
     *  A signed whole number, its magnitude in words of 64 bits, lowest first.
     *
     *  @tparam CapacityT The most words its magnitude takes.
     */
    template <std::size_t CapacityT>
    struct Whole
    {
      std::array<std::uint64_t, CapacityT> Words    = {};    /**< The magnitude, lowest first. */
      std::size_t                          Length   = 0;     /**< The words in use, 0 for 0.   */
      bool                                 Negative = false; /**< Whether it is below zero.    */

      /**
       *  The words in use.
       *
       *  @return The magnitude, its top word not zero.
       */
      [[nodiscard]] auto Magnitude() const -> std::span<const std::uint64_t>
      {
        return std::span{Words}.first(Length);
      }

      /**
       *  The sign.
       *
       *  @return 1, 0 or -1.
       */
      [[nodiscard]] auto Sign() const -> int
      {
        if (Length == 0)
        {
          return 0;
        }
        return Negative ? -1 : 1;
      }

      /**
       *  Drops the zero words at the top.
       */
      void Trim()
      {
        while (Length > 0 && Words.at(Length - 1) == 0)
        {
          --Length;
        }
      }
    };

    /**
     *  How two magnitudes compare.
     *
     *  @param  first  One magnitude, trimmed.
     *  @param  second The other, trimmed.
     *  @return        1, 0 or -1 as the first is larger, the same or smaller.
     */
    auto Ranked(std::span<const std::uint64_t> first, std::span<const std::uint64_t> second) -> int
    {
      if (first.size() != second.size())
      {
        return first.size() > second.size() ? 1 : -1;
      }
      for (std::size_t index = first.size(); index > 0; --index)
      {
        if (first[index - 1] != second[index - 1])
        {
          return first[index - 1] > second[index - 1] ? 1 : -1;
        }
      }
      return 0;
    }

    /**
     *  The sum of two magnitudes.
     *
     *  @pre           The sum has room for the longer and a carry.
     *  @param  first  One magnitude.
     *  @param  second The other.
     *  @param  sum    Where the sum goes, zeroed.
     *  @return        How many words it takes.
     */
    auto AddMagnitudes(std::span<const std::uint64_t> first, std::span<const std::uint64_t> second,
                       std::span<std::uint64_t> sum) -> std::size_t
    {
      const std::size_t longer = std::max(first.size(), second.size());
      std::uint64_t     carry  = 0;
      for (std::size_t index = 0; index < longer; ++index)
      {
        const std::uint64_t one     = index < first.size() ? first[index] : 0;
        const std::uint64_t two     = index < second.size() ? second[index] : 0;
        const std::uint64_t partial = one + two;
        const std::uint64_t total   = partial + carry;
        carry                       = static_cast<std::uint64_t>(partial < one || total < partial);
        sum[index]                  = total;
      }
      if (carry == 0)
      {
        return longer;
      }
      sum[longer] = carry;
      return longer + 1;
    }

    /**
     *  The difference of two magnitudes, the larger first.
     *
     *  @pre               The first is at least the second, and the difference has room for it.
     *  @param  larger     The magnitude subtracted from.
     *  @param  smaller    The magnitude subtracted.
     *  @param  difference Where the difference goes.
     *  @return            How many words the first takes.
     */
    auto SubtractMagnitudes(std::span<const std::uint64_t> larger,
                            std::span<const std::uint64_t> smaller,
                            std::span<std::uint64_t>       difference) -> std::size_t
    {
      std::uint64_t borrow = 0;
      for (std::size_t index = 0; index < larger.size(); ++index)
      {
        const std::uint64_t taken   = index < smaller.size() ? smaller[index] : 0;
        const std::uint64_t partial = larger[index] - taken;
        const std::uint64_t total   = partial - borrow;
        borrow            = static_cast<std::uint64_t>(larger[index] < taken || partial < borrow);
        difference[index] = total;
      }
      return larger.size();
    }

    /**
     *  The sum of two signed whole numbers.
     *
     *  @tparam CapacityT The capacity of the sum, room for either and a carry.
     *  @tparam FirstT    The capacity of one.
     *  @tparam SecondT   The capacity of the other.
     *  @param  first     One number.
     *  @param  second    The other.
     *  @param  negate    Whether the second is subtracted instead.
     *  @return           The sum.
     */
    template <std::size_t CapacityT, std::size_t FirstT, std::size_t SecondT>
    auto Plus(const Whole<FirstT> &first, const Whole<SecondT> &second, bool negate)
        -> Whole<CapacityT>
    {
      const bool       secondNegative = second.Negative != negate;
      Whole<CapacityT> sum;
      if (first.Negative == secondNegative)
      {
        sum.Length   = AddMagnitudes(first.Magnitude(), second.Magnitude(), sum.Words);
        sum.Negative = first.Negative;
      }
      else if (Ranked(first.Magnitude(), second.Magnitude()) >= 0)
      {
        sum.Length   = SubtractMagnitudes(first.Magnitude(), second.Magnitude(), sum.Words);
        sum.Negative = first.Negative;
      }
      else
      {
        sum.Length   = SubtractMagnitudes(second.Magnitude(), first.Magnitude(), sum.Words);
        sum.Negative = secondNegative;
      }
      sum.Trim();
      return sum;
    }

    /**
     *  The product of two signed whole numbers, schoolbook.
     *
     *  @tparam FirstT  The capacity of one.
     *  @tparam SecondT The capacity of the other.
     *  @param  first   One number.
     *  @param  second  The other.
     *  @return         The product, with room for both.
     */
    template <std::size_t FirstT, std::size_t SecondT>
    auto Multiple(const Whole<FirstT> &first, const Whole<SecondT> &second)
        -> Whole<FirstT + SecondT>
    {
      Whole<FirstT + SecondT> product;
      for (std::size_t one = 0; one < first.Length; ++one)
      {
        std::uint64_t carry = 0;
        for (std::size_t two = 0; two < second.Length; ++two)
        {
          const Product       part  = Multiply(first.Words.at(one), second.Words.at(two));
          std::uint64_t      &word  = product.Words.at(one + two);
          const std::uint64_t low   = part.Low + word;
          const std::uint64_t total = low + carry;
          carry                     = part.High
                                    + static_cast<std::uint64_t>(low < word)
                                    + static_cast<std::uint64_t>(total < low);
          word                      = total;
        }
        product.Words.at(one + second.Length) = carry;
      }
      product.Length   = first.Length + second.Length;
      product.Negative = first.Negative != second.Negative;
      product.Trim();
      return product;
    }

    /**
     *  A finite double read as a significand and the place of its lowest bit.
     */
    struct Place
    {
      std::uint64_t Significand = 0;     /**< Its significand, the hidden bit set if normal. */
      int           Exponent    = 0;     /**< The place of the significand's lowest bit.     */
      bool          Negative    = false; /**< Whether the double is below zero.              */
    };

    /**
     *  A finite double's significand and place, subnormal ones too.
     *
     *  @param  value The double, finite.
     *  @return       Its place: value = +-Significand 2^Exponent.
     */
    auto Unpacked(double value) -> Place
    {
      constexpr std::uint64_t Fraction = (std::uint64_t{1} << 52U) - 1;

      const auto bits     = std::bit_cast<std::uint64_t>(value);
      const auto biased   = static_cast<int>((bits >> 52U) & 0x7FFU);
      const bool negative = (bits >> 63U) != 0;
      if (biased == 0)
      {
        return {.Significand = bits & Fraction, .Exponent = -1074, .Negative = negative};
      }
      return {.Significand = (bits & Fraction) + (std::uint64_t{1} << 52U),
              .Exponent    = biased - 1075,
              .Negative    = negative};
    }

    /**
     *  The lowest place any of some finite doubles has a bit set in.
     *
     *  @param  values The doubles.
     *  @return        The place, or the largest int when all are zero.
     */
    auto Lowest(std::span<const double> values) -> int
    {
      int lowest = std::numeric_limits<int>::max();
      for (const double value : values)
      {
        const Place place = Unpacked(value);
        if (place.Significand != 0)
        {
          lowest = std::min(lowest, place.Exponent + std::countr_zero(place.Significand));
        }
      }
      return lowest;
    }

    /**
     *  A finite double as a whole number of grid units.
     *
     *  @param  value  The double.
     *  @param  lowest The grid, a place at or below the double's lowest set bit.
     *  @return        value 2^-lowest.
     */
    auto Quantized(double value, int lowest) -> Whole<Coordinate>
    {
      const Place       place = Unpacked(value);
      Whole<Coordinate> units;
      units.Negative = place.Negative;
      if (place.Significand == 0)
      {
        return units;
      }
      const int shift = place.Exponent - lowest;
      if (shift < 0)
      {
        units.Words.at(0) = place.Significand >> static_cast<unsigned>(-shift);
      }
      else
      {
        const auto word      = static_cast<std::size_t>(shift) / 64;
        const auto bit       = static_cast<unsigned>(shift) % 64;
        units.Words.at(word) = place.Significand << bit;
        if (bit != 0)
        {
          units.Words.at(word + 1) = place.Significand >> (64U - bit);
        }
      }
      units.Length = units.Words.size();
      units.Trim();
      return units;
    }

    /**
     *  Whether every double is finite.
     *
     *  @param  values The doubles.
     *  @return        Whether none is NaN or infinite.
     */
    auto Finite(std::span<const double> values) -> bool
    {
      return std::ranges::all_of(values, [](double value) -> bool { return std::isfinite(value); });
    }

    /**
     *  The difference of two coordinates in grid units.
     *
     *  @param  end   The coordinate arrived at.
     *  @param  start The coordinate left.
     *  @return       end - start.
     */
    auto Deducted(const Whole<Coordinate> &end, const Whole<Coordinate> &start) -> Whole<Coordinate>
    {
      return Plus<Coordinate>(end, start, true);
    }

    /**
     *  A corner's difference from the point, in grid units.
     */
    struct Offset
    {
      Whole<Coordinate> X; /**< Along the x axis. */
      Whole<Coordinate> Y; /**< Along the y axis. */
    };

    /**
     *  The cross product of two differences, exactly.
     *
     *  @param  first  The difference turned from.
     *  @param  second The difference turned to.
     *  @return        first.X second.Y - first.Y second.X.
     */
    auto Cross(const Offset &first, const Offset &second) -> Whole<2 * Coordinate>
    {
      return Plus<2 * Coordinate>(Multiple(first.X, second.Y), Multiple(first.Y, second.X), true);
    }

    /**
     *  A corner's term of the circle determinant, exactly.
     *
     *  @param  own    The corner's difference.
     *  @param  next   The corner after it.
     *  @param  before The corner before it.
     *  @return        (x^2 + y^2) times the cross product of the other two.
     */
    auto Summand(const Offset &own, const Offset &next, const Offset &before)
        -> Whole<4 * Coordinate>
    {
      const Whole<2 * Coordinate> lift =
          Plus<2 * Coordinate>(Multiple(own.X, own.X), Multiple(own.Y, own.Y), false);
      return Multiple(lift, Cross(next, before));
    }

    /**
     *  The words a crossing's coordinate's numerator takes: a coordinate times a cross product, and
     *  a way's times another, each below 2^6298, their sum below 2^6299.
     */
    constexpr std::size_t Numerator = (3 * Coordinate);

    /**
     *  The words that numerator takes on the grid of 2^-1075, the grid's place moved up by as much
     *  as 2046 bits: below 2^8345.
     */
    constexpr std::size_t GridWords = Numerator + 33;

    /**
     *  The words a halfway point's double takes on the grid of 2^-1075 times a cross product: below
     *  2^2099 times 2^4199.
     */
    constexpr std::size_t Halved = (3 * Coordinate) + 1;

    /**
     *  How many bits a magnitude takes, past its top word's leading zeros.
     *
     *  @tparam CapacityT The capacity.
     *  @param  whole     The magnitude, not zero.
     *  @return           Its width.
     */
    template <std::size_t CapacityT>
    auto Bitwidth(const Whole<CapacityT> &whole) -> int
    {
      const std::uint64_t top = whole.Words.at(whole.Length - 1);
      return (64 * static_cast<int>(whole.Length - 1)) + (64 - std::countl_zero(top));
    }

    /**
     *  A magnitude's leading 64 bits, as a double: the magnitude is that times 2^(width - 64),
     *  within 2^-63 of it before the conversion rounds.
     *
     *  @tparam CapacityT The capacity.
     *  @param  whole     The magnitude, not zero.
     *  @return           Its leading bits.
     */
    template <std::size_t CapacityT>
    auto Topmost(const Whole<CapacityT> &whole) -> double
    {
      const std::size_t   top   = whole.Length - 1;
      const auto          spare = static_cast<unsigned>(std::countl_zero(whole.Words.at(top)));
      const std::uint64_t below = top > 0 ? whole.Words.at(top - 1) : 0;
      const std::uint64_t bits  = (whole.Words.at(top) << spare) | ((below >> 1U) >> (63U - spare));
      return static_cast<double>(bits);
    }

    /**
     *  A magnitude moved up by some bits.
     *
     *  @tparam CapacityT The capacity of the result, room for the magnitude so moved.
     *  @tparam FromT     The capacity of the magnitude.
     *  @param  whole     The magnitude.
     *  @param  bits      How far, not below 0.
     *  @return           The magnitude times 2^bits.
     */
    template <std::size_t CapacityT, std::size_t FromT>
    auto Raised(const Whole<FromT> &whole, int bits) -> Whole<CapacityT>
    {
      const auto       words = static_cast<std::size_t>(bits) / 64;
      const auto       bit   = static_cast<unsigned>(bits) % 64;
      Whole<CapacityT> raised;
      for (std::size_t index = 0; index < whole.Length; ++index)
      {
        const std::uint64_t word            = whole.Words.at(index);
        raised.Words.at(index + words)     |= word << bit;
        raised.Words.at(index + words + 1)  = (word >> 1U) >> (63U - bit);
      }
      raised.Length = whole.Length + words + 1;
      raised.Trim();
      return raised;
    }

    /**
     *  A crossing's coordinate as a quotient of magnitudes: the numerator on the grid of 2^-1075,
     *  and the cross product of the ways below it.
     */
    struct Quotient
    {
      Whole<GridWords>      Above = {}; /**< The numerator's magnitude on 2^-1075's grid. */
      Whole<2 * Coordinate> Below = {}; /**< The cross product's magnitude.               */
    };

    /**
     *  How the halfway point above a double lies from a quotient: the double and the next as whole
     *  numbers of 2^-1074, their sum times the quotient's denominator against its numerator, past
     *  the largest double the next its 2^1024.
     *
     *  @param  quotient The quotient.
     *  @param  place    The double's place among the doubles not below zero, its bits.
     *  @return          1, 0 or -1 as the halfway point lies above the quotient, at it or below.
     */
    auto Midway(const Quotient &quotient, std::uint64_t place) -> int
    {
      const Whole<Coordinate>     low  = Quantized(std::bit_cast<double>(place), -1074);
      const Whole<Coordinate>     high = Quantized(std::bit_cast<double>(place + 1), -1074);
      const Whole<Coordinate + 1> sum  = Plus<Coordinate + 1>(low, high, false);
      const Whole<Halved>         mark = Multiple(sum, quotient.Below);
      return Ranked(mark.Magnitude(), quotient.Above.Magnitude());
    }

    /**
     *  The double nearest a quotient not below zero, ties to even, an infinity past the largest:
     *  the least place whose halfway point above lies at or above it, found between places eight
     *  either side of a rough quotient of the magnitudes' leading bits, which lies within six of
     *  the nearest, as their own roundings are three steps of a double at most, six in the binade
     *  below.
     *
     *  @param  quotient The quotient.
     *  @param  rough    The quotient, roughly: within 2^-50 of it, rounded once where subnormal.
     *  @return          The double nearest.
     */
    auto Closest(const Quotient &quotient, double rough) -> double
    {
      constexpr auto Endless =
          std::bit_cast<std::uint64_t>(std::numeric_limits<double>::infinity());
      constexpr std::uint64_t Margin = 8;

      const auto start = std::bit_cast<std::uint64_t>(rough);
      /**
       *  The least place's halfway point lies above the quotient: the one below it, if any, is the
       *  low end, and the high end's lies above.
       */
      std::uint64_t low  = start > Margin ? start - Margin : 0;
      std::uint64_t high = std::min(start + Margin, Endless);
      if (Midway(quotient, low) >= 0)
      {
        high = low;
      }
      while (high - low > 1)
      {
        const std::uint64_t middle                   = low + ((high - low) / 2);
        (Midway(quotient, middle) >= 0 ? high : low) = middle;
      }
      const bool tie = high < Endless && Midway(quotient, high) == 0;
      return std::bit_cast<double>(high + static_cast<std::uint64_t>(tie && (high & 1U) != 0));
    }
  } // namespace

  auto OrientWhole(const Triangle &corners) -> int
  {
    const std::array<double, 6> values = {corners.First.X,  corners.First.Y, corners.Second.X,
                                          corners.Second.Y, corners.Third.X, corners.Third.Y};
    if (!Finite(values))
    {
      return 0;
    }
    const int lowest = Lowest(values);
    if (lowest == std::numeric_limits<int>::max())
    {
      return 0;
    }
    const Whole<Coordinate> thirdX = Quantized(corners.Third.X, lowest);
    const Whole<Coordinate> thirdY = Quantized(corners.Third.Y, lowest);
    const Offset            first  = {.X = Deducted(Quantized(corners.First.X, lowest), thirdX),
                                      .Y = Deducted(Quantized(corners.First.Y, lowest), thirdY)};
    const Offset            second = {.X = Deducted(Quantized(corners.Second.X, lowest), thirdX),
                                      .Y = Deducted(Quantized(corners.Second.Y, lowest), thirdY)};
    return Cross(first, second).Sign();
  }

  auto CircleWhole(const Triangle &circle, Point point) -> int
  {
    const std::array<double, 8> values = {circle.First.X,  circle.First.Y, circle.Second.X,
                                          circle.Second.Y, circle.Third.X, circle.Third.Y,
                                          point.X,         point.Y};
    if (!Finite(values))
    {
      return 0;
    }
    const int lowest = Lowest(values);
    if (lowest == std::numeric_limits<int>::max())
    {
      return 0;
    }
    const Whole<Coordinate> pointX = Quantized(point.X, lowest);
    const Whole<Coordinate> pointY = Quantized(point.Y, lowest);
    const auto              offset = [&](Point corner) -> Offset
    {
      return {.X = Deducted(Quantized(corner.X, lowest), pointX),
              .Y = Deducted(Quantized(corner.Y, lowest), pointY)};
    };
    const Offset a = offset(circle.First);
    const Offset b = offset(circle.Second);
    const Offset c = offset(circle.Third);

    const Whole<4 * Coordinate> sum =
        Plus<4 * Coordinate>(Summand(a, b, c), Summand(b, c, a), false);
    return Plus<4 * Coordinate>(sum, Summand(c, a, b), false).Sign();
  }

  auto CrossWhole(const std::array<Straightway, 2> &ways) -> std::optional<Point>
  {
    const std::array<double, 8> values = {ways[0].Through.X, ways[0].Through.Y, ways[0].Toward.X,
                                          ways[0].Toward.Y,  ways[1].Through.X, ways[1].Through.Y,
                                          ways[1].Toward.X,  ways[1].Toward.Y};
    if (!Finite(values))
    {
      return std::nullopt;
    }
    const int lowest = Lowest(values);
    if (lowest == std::numeric_limits<int>::max())
    {
      return std::nullopt;
    }
    const auto offset = [lowest](const Straightway &way) -> std::array<Offset, 2>
    {
      const Offset through = {.X = Quantized(way.Through.X, lowest),
                              .Y = Quantized(way.Through.Y, lowest)};
      const Offset toward  = {.X = Quantized(way.Toward.X, lowest),
                              .Y = Quantized(way.Toward.Y, lowest)};
      if (way.Way)
      {
        return {through, toward};
      }
      return {through,
              Offset{.X = Deducted(toward.X, through.X), .Y = Deducted(toward.Y, through.Y)}};
    };
    const std::array<Offset, 2> one   = offset(ways[0]);
    const std::array<Offset, 2> other = offset(ways[1]);
    const Whole<2 * Coordinate> turn  = Cross(one[1], other[1]);
    if (turn.Sign() == 0)
    {
      return std::nullopt;
    }
    const Offset gap = {.X = Deducted(other[0].X, one[0].X), .Y = Deducted(other[0].Y, one[0].Y)};
    const Whole<2 * Coordinate> along = Cross(gap, other[1]);
    Whole<2 * Coordinate>       below = turn;
    below.Negative                    = false;
    /**
     *  Each coordinate (P t + d a) / t on the grid, t the turn and a how far along: its numerator
     *  moved up onto 2^-1075's grid, and its magnitudes' leading bits giving it roughly.
     */
    const auto axis = [&](const Whole<Coordinate> &start, const Whole<Coordinate> &step) -> double
    {
      const Whole<Numerator> top =
          Plus<Numerator>(Multiple(start, turn), Multiple(step, along), false);
      if (top.Sign() == 0)
      {
        return turn.Negative ? -0.0 : 0.0;
      }
      const Quotient quotient = {.Above = Raised<GridWords>(top, lowest + 1075), .Below = below};
      const double   rough =
          std::ldexp(Topmost(top) / Topmost(below), Bitwidth(top) - Bitwidth(below) + lowest);
      const double nearest = Closest(quotient, rough);
      return top.Negative != turn.Negative ? -nearest : nearest;
    };
    return Point{.X = axis(one[0].X, one[1].X), .Y = axis(one[0].Y, one[1].Y)};
  }
} // namespace Geoxy::Exact
