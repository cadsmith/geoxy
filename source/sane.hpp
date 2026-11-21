/***************************************************************************************************
 *  @file      source/sane.hpp
 *  @brief     A tolerance made fit to work with, as every meeting takes it first.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SANE_HPP
#define GEOXY_SANE_HPP

#include "geoxy/tolerance.hpp"

#include <cmath>
#include <limits>

namespace Geoxy
{
  /**
   *  A tolerance made fit to work with: a length or an angle that is not a number, or below 0, is
   *  0; one past every double is the largest double.
   *
   *  A length that is not a number fails every test of a distance against it, so a meeting within
   *  it found nothing and said nothing of it; one without end has every point near every other, and
   *  joining every hit to every other cost the square of their count, each join a search for the
   *  nearest point, without bound on curves whose numbers are not finite. Taken so, the first is
   *  met as exactly as the curves' numbers allow, and the second within the largest length a double
   *  holds.
   *
   *  @param  tolerance The tolerance asked for.
   *  @return           The tolerance to meet within.
   */
  inline auto Sane(const Tolerance &tolerance) -> Tolerance
  {
    constexpr double Most = std::numeric_limits<double>::max();
    const auto       fit  = [](double value) -> double
    {
      return value >= 0.0 ? std::fmin(value, Most) : 0.0;
    };
    return {.Length = fit(tolerance.Length), .Radians = fit(tolerance.Radians)};
  }
} // namespace Geoxy

#endif /* GEOXY_SANE_HPP */
