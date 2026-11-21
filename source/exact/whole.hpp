/***************************************************************************************************
 *  @file      source/exact/whole.hpp
 *  @brief     Both determinants, and where two lines cross, in whole numbers of any size: exact for
 *             every finite double.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_EXACT_WHOLE_HPP
#define GEOXY_EXACT_WHOLE_HPP

#include "exact/triangle.hpp"

#include "geoxy/point.hpp"

#include <array>
#include <optional>

/**
 *  Every finite double is a whole number of 2^-1074, below 2^2098 of them: so are the coordinates
 *  of any points, their differences below 2^2099, the orientation determinant below 2^4199 and the
 *  circle determinant below 2^8400. The two are summed here in words of 64 bits, sign and
 *  magnitude, on the grid of the lowest set bit among the coordinates.
 *
 *  It is the predicates' last resort, for places the adaptive stages' arithmetic cannot hold and
 *  the grids of exact/grid.hpp do not reach. Slow beside them, and never wrong: tests/proof/
 *  integers.v proves it step for step, each sign the sign of the determinant of the doubles as IEEE
 *  754 reads their bits.
 *
 *  Where two lines cross is worked out on the same grid: each way's cross product with the other
 *  below 2^4199, each coordinate's numerator below 2^6299, the coordinate their quotient, rounded
 *  by comparing it with the halfway points between doubles in whole numbers of 2^-1075.
 */

namespace Geoxy::Exact
{
  /**
   *  Which way three points turn, exactly, wherever they lie.
   *
   *  @param  corners The triangle.
   *  @return         The sign of Orient's determinant; 0 where a coordinate is NaN or infinite.
   */
  auto OrientWhole(const Triangle &corners) -> int;

  /**
   *  Where a point lies from the circle through three corners, exactly, wherever they lie.
   *
   *  @param  circle The corners the circle passes through.
   *  @param  point  The point asked about.
   *  @return        The sign of InCircle's determinant; 0 where a coordinate is NaN or infinite.
   */
  auto CircleWhole(const Triangle &circle, Point point) -> int;

  /**
   *  A straight carrier as CrossWhole reads it: a point it runs through, and where it heads from
   *  there, a second point it runs through or its way.
   */
  struct Straightway
  {
    Point Through = {};    /**< A point it runs through.                                   */
    Point Toward  = {};    /**< A second point it runs through, or its way from the first. */
    bool  Way     = false; /**< Whether Toward is the way, not a second point.             */
  };

  /**
   *  Where two straight carriers cross, exactly, wherever they lie: worked out in whole numbers of
   *  the lowest place any of their numbers sets, each coordinate the double nearest the true one,
   *  ties to even, an infinity past the largest double, and a zero of the sign of the ways' cross
   *  product, as a zero over it has. The last resort of a crossing whose numbers part too far for
   *  sums of doubles to hold their products.
   *
   *  @param  ways The two carriers.
   *  @return      The crossing; nothing where they are parallel, or a number is not finite.
   */
  auto CrossWhole(const std::array<Straightway, 2> &ways) -> std::optional<Point>;
} // namespace Geoxy::Exact

#endif /* GEOXY_EXACT_WHOLE_HPP */
