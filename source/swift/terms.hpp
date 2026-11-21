/***************************************************************************************************
 *  @file      source/swift/terms.hpp
 *  @brief     The counts and the terms of the polynomials a curve's kernel is worked out in.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_TERMS_HPP
#define GEOXY_SWIFT_TERMS_HPP

#include "sure.hpp"
#include "swift/casteljau.hpp"

#include <array>
#include <cstddef>

namespace Geoxy::Swift
{
  using Sure::Grow;

  using Sure::Number;

  using Sure::Shrink;

  /**
   *  A piece's numbers, one for each control point a room holds.
   *
   *  @tparam ValueT The numbers' kind.
   *  @tparam CountT The room: the piece's degree and one at the most.
   */
  template <typename ValueT, std::size_t CountT>
  using Poled = std::array<ValueT, CountT>;

  /**
   *  How many Bernstein coefficients a conic's equation on a piece of a room holds: twice the
   *  room's degree and one.
   *
   *  @tparam CountT The piece's room.
   */
  template <std::size_t CountT>
  inline constexpr std::size_t Ranks = (2 * (CountT - 1)) + 1;

  /**
   *  A Bézier curve's room, a quadratic's or a cubic's: its piece, and a quarter of an arc or an
   *  ellipse, are met in so many cells.
   */
  inline constexpr std::size_t Thin = Highest + 1;

  /**
   *  A spline's span's room, of any degree to Broadest.
   */
  inline constexpr std::size_t Ample = Broadest + 1;

  /**
   *  The highest degree of a conic's equation on a quadratic or cubic piece, whose isolation and
   *  polish keep to so many coefficients.
   */
  inline constexpr std::size_t Narrow = 2 * Highest;

  /**
   *  The equation's Bernstein coefficients on a piece of a room.
   *
   *  @tparam CountT The piece's room.
   */
  template <std::size_t CountT>
  using Terms = std::array<double, Ranks<CountT>>;

  /**
   *  A binomial coefficient, n choose k, exact for every degree the kernel takes.
   *
   *  @param  count How many, n.
   *  @param  taken How many taken, k.
   *  @return       The coefficient.
   */
  constexpr auto Choose(std::size_t count, std::size_t taken) -> double
  {
    double out = 1.0;
    for (std::size_t index = 0; index < taken; ++index)
    {
      out = out * static_cast<double>(count - index) / static_cast<double>(index + 1);
    }
    return out;
  }

  /**
   *  A direction turned by some quarters of a turn, exactly.
   *
   *  @tparam ValueT   The coordinates' kind.
   *  @param  way      The direction's coordinates.
   *  @param  quarters How many quarters, counterclockwise.
   *  @return          The coordinates turned.
   */
  template <typename ValueT>
  auto Quartered(std::array<ValueT, 2> way, std::size_t quarters) -> std::array<ValueT, 2>
  {
    for (std::size_t turn = 0; turn < quarters; ++turn)
    {
      way = {-way[1], way[0]};
    }
    return way;
  }
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_TERMS_HPP */
