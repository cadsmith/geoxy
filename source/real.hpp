/***************************************************************************************************
 *  @file      source/real.hpp
 *  @brief     Doubles compared exactly, and runs of them put in order.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_REAL_HPP
#define GEOXY_REAL_HPP

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>

namespace Geoxy
{
  /**
   *  The positive infinity of doubles.
   */
  inline constexpr double Infinite = std::numeric_limits<double>::infinity();

  /**
   *  Whether a double is exactly a value, where exactness is what is meant.
   *
   *  @param  value   The double.
   *  @param  exactly The value.
   *  @return         Whether they are equal.
   */
  inline auto Exactly(double value, double exactly) -> bool
  {
    return std::equal_to<double>{}(value, exactly);
  }

  /**
   *  Whether a double is neither above nor below zero.
   *
   *  @param  value The double.
   *  @return       Whether it is zero of either sign, or NaN, which has no sign to go by.
   */
  inline auto Zero(double value) -> bool
  {
    return !(value < 0.0) && !(value > 0.0);
  }

  /**
   *  The smaller step from a double to its neighbours: its unit in the last place, half that at a
   *  power of two, where the neighbour below is nearer, read from its bits, as nextafter would give
   *  it at a call's cost; the subnormal step at zero and below the normals; not a number for an
   *  infinity or not a number.
   *
   *  @param  value The double.
   *  @return       The step.
   */
  inline auto Spacing(double value) -> double
  {
    constexpr std::uint64_t Fraction = (std::uint64_t{1} << 52U) - 1U;
    const std::uint64_t     bits = std::bit_cast<std::uint64_t>(value) & ~(std::uint64_t{1} << 63U);
    const std::uint64_t     field = bits >> 52U;
    if (field == 0x7FF)
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    if (field <= 1)
    {
      return 0x1p-1074;
    }
    const double unit = field > 52 ? std::bit_cast<double>((field - 52) << 52U)
                                   : std::bit_cast<double>(std::uint64_t{1} << (field - 1));
    return (bits & Fraction) == 0 ? unit / 2 : unit;
  }

  /**
   *  The step from a magnitude up to the next double, read from its bits, as nextafter toward
   *  infinity would give it at a call's cost: its unit in the last place, the subnormal step at
   *  zero, an infinity past the largest double, not a number for an infinity or not a number.
   *
   *  @param  magnitude The magnitude, not below zero.
   *  @return           The step.
   */
  inline auto Above(double magnitude) -> double
  {
    return std::bit_cast<double>(std::bit_cast<std::uint64_t>(magnitude) + 1U) - magnitude;
  }

  /**
   *  A run of doubles sorted, each value kept once.
   *
   *  @tparam ValuesT The list's kind.
   *  @param  values  The run.
   */
  template <typename ValuesT>
  void Ordered(ValuesT &values)
  {
    std::ranges::sort(values);
    const auto repeated = std::ranges::unique(values);
    values.erase(repeated.begin(), repeated.end());
  }

  /**
   *  A run of doubles sorted, values a few steps of a double apart kept once, the least of them:
   *  one root found twice, from either side of a seam, comes out once.
   *
   *  @tparam ValuesT The list's kind.
   *  @param  values  The run.
   */
  template <typename ValuesT>
  void Gathered(ValuesT &values)
  {
    std::ranges::sort(values);
    const auto close = [](double lhs, double rhs) -> bool
    {
      const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
      return rhs - lhs <= 64.0 * std::numeric_limits<double>::epsilon() * scale;
    };
    const auto repeated = std::ranges::unique(values, close);
    values.erase(repeated.begin(), repeated.end());
  }
} // namespace Geoxy

#endif /* GEOXY_REAL_HPP */
