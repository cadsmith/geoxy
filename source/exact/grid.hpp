/***************************************************************************************************
 *  @file      source/exact/grid.hpp
 *  @brief     Both tests' determinants in integers, where every coordinate lies on one grid.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_EXACT_GRID_HPP
#define GEOXY_EXACT_GRID_HPP

#include "exact/triangle.hpp"

#include "geoxy/point.hpp"

#include <cstdint>
#include <optional>

/**
 *  A double is a whole number of its last bit's place. Where every coordinate of the four points
 *  spans at most 61 bits, from the highest bit of the largest down to the lowest set bit of any,
 *  each is a whole number of that lowest place below 2^61, its differences below 2^62, and the
 *  determinant a sum of three products of four of them: below 2^252, exact in 256 bits of integer.
 *  Survey and CAD coordinates snapped to a grid span far less, and so do integer lattices, the
 *  places where the circle test meets exact zeros most.
 *
 *  The orientation test's determinant is a difference of two products of two: there each coordinate
 *  may span 125 bits, two words, its differences below 2^126 and the determinant below 2^253. Any
 *  doubles of like magnitude span 53 bits; ones whose magnitudes part by a factor up to 2^72 still
 *  fit.
 *
 *  tests/proof/integers.v proves both, and MultiplyHalves, step for step: each sign the sign of the
 *  determinant of the doubles as IEEE 754 reads their bits.
 */

namespace Geoxy::Exact
{
  /**
   *  A product of two 64 bit words.
   */
  struct Product
  {
    std::uint64_t High = 0; /**< The upper 64 bits. */
    std::uint64_t Low  = 0; /**< The lower 64 bits. */
  };

  /**
   *  The product of two words, from their 32 bit halves: every compiler's.
   *
   *  @param  first  One word.
   *  @param  second The other.
   *  @return        Their 128 bit product.
   */
  constexpr auto MultiplyHalves(std::uint64_t first, std::uint64_t second) -> Product
  {
    constexpr std::uint64_t Half = 0xFFFF'FFFFU;

    const std::uint64_t lowLow   = (first & Half) * (second & Half);
    const std::uint64_t lowHigh  = (first & Half) * (second >> 32U);
    const std::uint64_t highLow  = (first >> 32U) * (second & Half);
    const std::uint64_t highHigh = (first >> 32U) * (second >> 32U);
    const std::uint64_t middle   = (lowLow >> 32U) + (lowHigh & Half) + (highLow & Half);
    return {.High = highHigh + (lowHigh >> 32U) + (highLow >> 32U) + (middle >> 32U),
            .Low  = (middle << 32U) + (lowLow & Half)};
  }

  /**
   *  The product of two words: one instruction where the compiler has 128 bit integers,
   *  MultiplyHalves elsewhere.
   *
   *  @param  first  One word.
   *  @param  second The other.
   *  @return        Their 128 bit product.
   */
  auto Multiply(std::uint64_t first, std::uint64_t second) -> Product;

  /**
   *  The sign of the circle test's determinant in integers, where the coordinates lie on one grid.
   *
   *  @param  circle The corners the circle passes through.
   *  @param  point  The point asked about.
   *  @return        1, 0 or -1, as InCircle; none where a coordinate is not normal or zero, or the
   *                 coordinates span more than 61 bits.
   */
  auto OnGrid(const Triangle &circle, Point point) -> std::optional<int>;

  /**
   *  The sign of the orientation test's determinant in integers, where the coordinates lie on one
   *  grid.
   *
   *  @param  corners The triangle.
   *  @return         1, 0 or -1, as Orient; none where a coordinate is not normal or zero, or the
   *                  coordinates span more than 125 bits.
   */
  auto OrientGrid(const Triangle &corners) -> std::optional<int>;
} // namespace Geoxy::Exact

#endif /* GEOXY_EXACT_GRID_HPP */
