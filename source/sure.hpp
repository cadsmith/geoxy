/***************************************************************************************************
 *  @file      source/sure.hpp
 *  @brief     Numbers in double-double that carry a bound on how far they may be from the value
 *             they stand for, and the certificate that rounds one to the double nearest that value.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SURE_HPP
#define GEOXY_SURE_HPP

#include "real.hpp"
#include "twofold.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>

namespace Geoxy::Sure
{
  /**
   *  How much an operation's own rounding may add, as a share of its operands' sizes, for a sum:
   *  the sloppy sum of two double-doubles is off by under 6 u^2 of them, u = 2^-53, and this is
   *  2^-100, some ten times that.
   */
  inline constexpr double SumGrain = 0x1p-100;

  /**
   *  For a product: the double-double product is off by under 9 u^2 of its operands' product, and
   *  this is 2^-100, some seven times that.
   */
  inline constexpr double ProductGrain = 0x1p-100;

  /**
   *  For a quotient or a square root: each is off by under 16 u^2 of its result, and this is 2^-98,
   *  sixteen times that.
   */
  inline constexpr double QuotientGrain = 0x1p-98;

  /**
   *  What a bound is grown by once worked out, so that its own few roundings, each under 2^-53 of
   *  it, can never shrink it below what it stands for.
   */
  inline constexpr double Grow = 1.0 + 0x1p-40;

  /**
   *  What a bound is shrunk by where it bounds from below, for the same reason.
   */
  inline constexpr double Shrink = 1.0 - 0x1p-40;

  /**
   *  An absolute part every bound worked out carries: what roundings below the normal doubles may
   *  lose, under 2^-1070 each, many times over. A value as small as this is never certain.
   */
  inline constexpr double Tiny = 0x1p-1000;

  /**
   *  A number in double-double and how far it may be from the exact value it stands for: |High +
   *  Low - exact| <= Bound. Its two doubles are normalized, |Low| no more than half a step of a
   *  double at High.
   */
  struct Number
  {
    Twofold::Double Value = {};  /**< The number.                            */
    double          Bound = 0.0; /**< How far from the exact value, at most. */
  };

  /**
   *  A double, exactly.
   *
   *  @param  value The double.
   *  @return       It.
   */
  inline auto Exact(double value) -> Number
  {
    return {.Value = {.High = value, .Low = 0.0}, .Bound = 0.0};
  }

  /**
   *  The difference of two doubles, exactly, by TwoSum, which is exact whatever their size.
   *
   *  @param  head The double.
   *  @param  tail The double taken from it.
   *  @return      The difference.
   */
  inline auto Apart(double head, double tail) -> Number
  {
    const Twofold::Two difference = Twofold::TwoSum(head, -tail);
    return {.Value = {.High = difference.High, .Low = difference.Low}, .Bound = 0.0};
  }

  /**
   *  The product of two doubles, by TwoProduct: exact unless it falls among the doubles below the
   *  normal ones, where it may lose what Tiny covers.
   *
   *  @param  first  One double.
   *  @param  second The other.
   *  @return        The product.
   */
  inline auto Product(double first, double second) -> Number
  {
    const Twofold::Two product = Twofold::TwoProduct(first, second);
    const double       bound   = std::abs(product.High) < 0x1p-960 ? Tiny : 0.0;
    return {.Value = {.High = product.High, .Low = product.Low}, .Bound = bound};
  }

  /**
   *  The size of a number, |High| + |Low|.
   *
   *  @param  number The number.
   *  @return        Its size.
   */
  inline auto Size(const Number &number) -> double
  {
    return std::abs(number.Value.High) + std::abs(number.Value.Low);
  }

  /**
   *  A number with its sign changed, exactly.
   *
   *  @param  number The number.
   *  @return        Its negative.
   */
  inline auto operator-(const Number &number) -> Number
  {
    return {.Value = Twofold::Negated(number.Value), .Bound = number.Bound};
  }

  /**
   *  The sum of two numbers: their values summed in double-double, their bounds added to the sum's
   *  own rounding.
   *
   *  @param  lhs One number.
   *  @param  rhs The other.
   *  @return     The sum.
   */
  inline auto operator+(const Number &lhs, const Number &rhs) -> Number
  {
    const Twofold::Double value = Twofold::Add(lhs.Value, rhs.Value);
    const double          bound =
        (lhs.Bound + rhs.Bound + (SumGrain * (Size(lhs) + Size(rhs))) + Tiny) * Grow;
    return {.Value = value, .Bound = bound};
  }

  /**
   *  The difference of two numbers.
   *
   *  @param  lhs The number.
   *  @param  rhs The number taken from it.
   *  @return     The difference.
   */
  inline auto operator-(const Number &lhs, const Number &rhs) -> Number
  {
    return lhs + (-rhs);
  }

  /**
   *  The product of two numbers: |x y - x' y'| <= |x'| e_y + |y'| e_x + e_x e_y for values x', y'
   *  within e_x and e_y of x and y, and the product's own rounding.
   *
   *  @param  lhs One number.
   *  @param  rhs The other.
   *  @return     The product.
   */
  inline auto operator*(const Number &lhs, const Number &rhs) -> Number
  {
    const Twofold::Double value = Twofold::Multiply(lhs.Value, rhs.Value);
    const double          left  = Size(lhs);
    const double          right = Size(rhs);
    const double          bound = ((left * rhs.Bound)
                                   + (right * lhs.Bound)
                                   + (lhs.Bound * rhs.Bound)
                                   + (ProductGrain * left * right)
                                   + Tiny)
                                * Grow;
    return {.Value = value, .Bound = bound};
  }

  /**
   *  The least a number's magnitude may be: its size less its bound, below zero where the bound
   *  reaches zero.
   *
   *  @param  number The number.
   *  @return        The least.
   */
  inline auto Least(const Number &number) -> double
  {
    return (std::abs(number.Value.High) - std::abs(number.Value.Low) - number.Bound) * Shrink;
  }

  /**
   *  The most a number's magnitude may be.
   *
   *  @param  number The number.
   *  @return        The most.
   */
  inline auto Most(const Number &number) -> double
  {
    return (Size(number) + number.Bound) * Grow;
  }

  /**
   *  The quotient of two numbers: |x / y - x' / y'| <= (e_x + |x' / y'| e_y) / |y| for values
   *  within e_x and e_y, |y| at least the denominator's least, and the quotient's own rounding.
   *  Where the denominator's bound reaches zero the bound is infinite, or not a number: no decision
   *  a kernel makes from it is then sure, and Rounded settles nothing.
   *
   *  @param  top    The numerator.
   *  @param  bottom The denominator.
   *  @return        The quotient.
   */
  inline auto Quotient(const Number &top, const Number &bottom) -> Number
  {
    const double          least = Least(bottom);
    const Twofold::Double value = Twofold::Divide(top.Value, bottom.Value);
    const double          size  = (std::abs(value.High) + std::abs(value.Low)) * Grow;
    const double          bound =
        least > 0.0
            ? (((top.Bound + (size * bottom.Bound)) / least) + (QuotientGrain * size) + Tiny) * Grow
            : std::numeric_limits<double>::infinity();
    return {.Value = value, .Bound = bound};
  }

  /**
   *  The square root of a number: |sqrt(x) - sqrt(x')| <= |x - x'| / (sqrt(x) + sqrt(x')), no more
   *  than e_x / sqrt(m) for m the least either may be, and the root's own rounding. Where the
   *  number's bound reaches zero, or it is below zero, the bound is infinite or not a number, as a
   *  quotient's is.
   *
   *  @param  number The number.
   *  @return        The root.
   */
  inline auto Root(const Number &number) -> Number
  {
    const double          least = (number.Value.High + number.Value.Low - number.Bound) * Shrink;
    const Twofold::Double value = Twofold::Root(number.Value);
    const double          size  = std::abs(value.High) + std::abs(value.Low);
    const double          bound =
        least > 0.0
            ? ((number.Bound / std::sqrt(least)) + (QuotientGrain * size) + Tiny) * Grow * Grow
            : std::numeric_limits<double>::infinity();
    return {.Value = value, .Bound = bound};
  }

  /**
   *  The double nearest the exact value a number stands for, where its bound settles it: with the
   *  number r + e, r a double and e what is left, the exact value lies within b of it, and r is the
   *  nearest double to every value from r + e - b to r + e + b exactly when both of those round to
   *  r, as round to nearest is monotone. The hardware rounds both, ties to even as the answer must;
   *  e + b and e - b are each rounded away from r + e by the growth of b, so the values tested lie
   *  past the true ends. A value exactly halfway between two doubles is never settled, as one of
   *  its ends rounds past it.
   *
   *  @param  number The number.
   *  @return        The double nearest its exact value; nothing where the bound does not settle it,
   *                 or it is not finite.
   */
  inline auto Rounded(const Number &number) -> std::optional<double>
  {
    const double rounded = number.Value.High + number.Value.Low;
    const double rest    = (number.Value.High - rounded) + number.Value.Low;
    const double bound   = (number.Bound * Grow) + (std::abs(rest) * 0x1p-50);
    if (std::isfinite(rounded)
        && Exactly(rounded + (rest + bound), rounded)
        && Exactly(rounded + (rest - bound), rounded))
    {
      return rounded;
    }
    return std::nullopt;
  }

  /**
   *  The doubles nearest the exact values some numbers stand for, where their bounds settle every
   *  one, as Rounded settles each.
   *
   *  @tparam CountT  How many.
   *  @param  numbers The numbers.
   *  @return         The doubles; nothing where any is not settled.
   */
  template <std::size_t CountT>
  auto Rounded(const std::array<Number, CountT> &numbers)
      -> std::optional<std::array<double, CountT>>
  {
    std::array<double, CountT> out = {};
    for (std::size_t index = 0; index < CountT; ++index)
    {
      const std::optional<double> rounded = Rounded(numbers.at(index));
      if (!rounded)
      {
        return std::nullopt;
      }
      out.at(index) = *rounded;
    }
    return out;
  }
} // namespace Geoxy::Sure

#endif /* GEOXY_SURE_HPP */
