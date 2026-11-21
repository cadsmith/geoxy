/***************************************************************************************************
 *  @file      source/swift/casteljau.hpp
 *  @brief     A polynomial held by its Bernstein coefficients worked out at a place by de
 *             Casteljau's algorithm, in doubles, or in double-double with a bound, as the swift
 *             kernels for Bézier curves take it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_CASTELJAU_HPP
#define GEOXY_SWIFT_CASTELJAU_HPP

#include "sure.hpp"

#include <array>
#include <cstddef>

namespace Geoxy::Swift
{
  /**
   *  The highest degree the kernels take: a cubic's.
   */
  inline constexpr std::size_t Highest = 3;

  /**
   *  The highest degree of a spline's span the kernels of a conic take: twice it, the degree of the
   *  conic's equation on the span, is sixteen.
   */
  inline constexpr std::size_t Broadest = 8;

  /**
   *  Bernstein coefficients, as many as the degree and one.
   */
  using Coefficients = std::array<double, Highest + 1>;

  /**
   *  Bernstein coefficients in double-double, each with its bound.
   */
  using Bounded = std::array<Sure::Number, Highest + 1>;

  /**
   *  A polynomial's value by de Casteljau's algorithm in doubles: each round the convex combination
   *  of two neighbours, c (1 - t) + c' t, so that every value stays among the coefficients' and a
   *  round is off by no more than three steps of a double at their largest, and what the
   *  coefficients were off by (tests/proof/twin.py's casteljau).
   *
   *  @tparam CountT       How many coefficients it holds room for.
   *  @param  coefficients Its Bernstein coefficients.
   *  @param  degree       Its degree.
   *  @param  along        Where, from 0 to 1.
   *  @return              The value.
   */
  template <std::size_t CountT>
  auto Valued(std::array<double, CountT> coefficients, std::size_t degree, double along) -> double
  {
    const double rest = 1.0 - along;
    for (std::size_t round = 1; round <= degree; ++round)
    {
      for (std::size_t index = 0; index + round <= degree; ++index)
      {
        coefficients.at(index) =
            (coefficients.at(index) * rest) + (coefficients.at(index + 1) * along);
      }
    }
    return coefficients[0];
  }

  /**
   *  A polynomial's value by de Casteljau's algorithm in double-double, at a place on the grid of
   *  2^-53, so that 1 - t is a double too and every weight is exact.
   *
   *  @tparam CountT       How many coefficients it holds room for.
   *  @param  coefficients Its Bernstein coefficients, with their bounds.
   *  @param  degree       Its degree.
   *  @param  along        Where, on the grid.
   *  @return              The value, with its bound.
   */
  template <std::size_t CountT>
  auto Valued(std::array<Sure::Number, CountT> coefficients, std::size_t degree, double along)
      -> Sure::Number
  {
    const Sure::Number ahead = Sure::Exact(along);
    const Sure::Number rest  = Sure::Exact(1.0 - along);
    for (std::size_t round = 1; round <= degree; ++round)
    {
      for (std::size_t index = 0; index + round <= degree; ++index)
      {
        coefficients.at(index) =
            (coefficients.at(index) * rest) + (coefficients.at(index + 1) * ahead);
      }
    }
    return coefficients[0];
  }
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_CASTELJAU_HPP */
