/***************************************************************************************************
 *  @file      tests/geoxy/exact/same.hpp
 *  @brief     Whether the components of an expansion are the original's.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TESTS_EXACT_SAME_HPP
#define GEOXY_TESTS_EXACT_SAME_HPP

#include "exact/two.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>

namespace Geoxy::Exact
{
  /**
   *  The bits of a double.
   */
  inline auto Bits(double value) -> std::uint64_t
  {
    return std::bit_cast<std::uint64_t>(value);
  }

  /**
   *  Whether a component is the original's: the same bits, but for the sign of a zero, since a
   *  fused product's zero roundoff, or a square's, is +0 where Shewchuk's split can give -0
   *  (exact/two.hpp).
   */
  inline auto Same(double ours, double original) -> bool
  {
    if (IsZero(ours) && IsZero(original))
    {
      return true;
    }
    return Bits(ours) == Bits(original);
  }

  /**
   *  Whether two expansions are the same, component by component.
   */
  inline auto Same(std::span<const double> ours, std::span<const double> original) -> bool
  {
    if (ours.size() != original.size())
    {
      return false;
    }
    for (std::size_t index = 0; index < ours.size(); ++index)
    {
      if (!Same(ours[index], original[index]))
      {
        return false;
      }
    }
    return true;
  }
} // namespace Geoxy::Exact

#endif /* GEOXY_TESTS_EXACT_SAME_HPP */
