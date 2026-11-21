/***************************************************************************************************
 *  @file      source/exact/orient.hpp
 *  @brief     The exact orientation test, its first stage inline.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_EXACT_ORIENT_HPP
#define GEOXY_EXACT_ORIENT_HPP

#include "exact/bounds.hpp"
#include "exact/triangle.hpp"
#include "exact/two.hpp"

#include "geoxy/point.hpp"

#include <cmath>

/**
 *  Nearly every call ends in the first stage, a few operations, so the call costs more than the
 *  test: the first stage is inline, and the library's callers compile it with their own loops. The
 *  later stages stay out of line.
 *
 *  The first stage is one test, |det| > bound (|left + right| + u_N), as Ozaki et al. write it and
 *  Agafonkin's robust-predicates does without u_N. Shewchuk branches on the signs of the two
 *  products first, which random input mispredicts. When the signs agree, |left + right| is his
 *  magnitude; when they differ, |det| >= max(|left|, |right|) >= |left + right| decides, as his
 *  branches do, unless both products are within a few u_N of zero. So the one test decides the
 *  cases his decide, and none wrongly, whatever the doubles (Decided, below).
 *
 *  The header is the library's own, compiled with contraction off: inline in a caller built
 *  otherwise, its products could fuse and break the bound.
 */

namespace Geoxy::Exact
{
  /**
   *  The first stage's approximation of the determinant, and the magnitude its error is relative
   *  to.
   */
  struct Approximation
  {
    double Determinant = 0.0; /**< The determinant, rounded.                      */
    double Magnitude   = 0.0; /**< The sum of the magnitudes of its two products. */
  };

  /**
   *  The orientation test's first stage.
   *
   *  @param  corners The triangle.
   *  @return         The determinant rounded, and its magnitude.
   */
  inline auto Approximate(const Triangle &corners) -> Approximation
  {
    const double left  = (corners.First.X - corners.Third.X) * (corners.Second.Y - corners.Third.Y);
    const double right = (corners.First.Y - corners.Third.Y) * (corners.Second.X - corners.Third.X);
    return {.Determinant = left - right, .Magnitude = std::abs(left + right)};
  }

  /**
   *  Whether the first stage's sign is beyond doubt, for every double.
   *
   *  Ozaki, Bünger, Ogita, Oishi and Rump, "Simple floating-point filters for the two-dimensional
   *  orientation problem" (2016), (3.31): |det| > theta (|left + right| + u_N) never passes a wrong
   *  sign, through underflow and overflow alike, since an infinity or a NaN fails the comparison.
   *  Their theta, 3u - 94906242u^2, takes 55 bits; Shewchuk's (3 + 16u)u is larger, and every step
   *  of their proof holds for it. The test is strict, so a zero determinant is never decided here,
   *  where a product's underflow could have made it. tests/proof/first.py checks it with Z3.
   *
   *  @param  approximation The first stage's approximation.
   *  @return               Whether it is decided.
   */
  inline auto Decided(const Approximation &approximation) -> bool
  {
    return std::abs(approximation.Determinant) > OrientBoundA * (approximation.Magnitude + Normal);
  }

  /**
   *  The orientation test's stages past the first.
   *
   *  Where the differences of the corners are exact, one fused comparison; what it leaves open,
   *  whole numbers (exact/grid.hpp, exact/whole.hpp).
   *
   *  @param  corners The triangle.
   *  @return         1, 0 or -1.
   */
  auto OrientAdapt(const Triangle &corners) -> int;

  /**
   *  Which way three points turn, exactly: the sign of twice the triangle's signed area.
   *
   *  Shewchuk's first stage, Ozaki's guard against underflow, then one fused comparison or whole
   *  numbers: the sign alone, exact for every finite double.
   *
   *  @param  first  The first corner.
   *  @param  second The second corner.
   *  @param  third  The corner whose side of the line from the first to the second is asked.
   *  @return        1 counterclockwise, -1 clockwise, 0 collinear; 0 too where a coordinate is NaN
   *                 or infinite.
   */
  [[nodiscard]] inline auto Orient(Point first, Point second, Point third) -> int
  {
    const Triangle      corners       = {.First = first, .Second = second, .Third = third};
    const Approximation approximation = Approximate(corners);
    if (Decided(approximation)) [[likely]]
    {
      return Sign(approximation.Determinant);
    }
    return OrientAdapt(corners);
  }
} // namespace Geoxy::Exact

#endif /* GEOXY_EXACT_ORIENT_HPP */
