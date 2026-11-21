/***************************************************************************************************
 *  @file      source/exact/orient.cpp
 *  @brief     The exact orientation test's stages past the first.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/orient.hpp"

#include "exact/grid.hpp"
#include "exact/roundoff.hpp"
#include "exact/triangle.hpp"
#include "exact/two.hpp"
#include "exact/whole.hpp"

#include "geoxy/vector.hpp"

#include <cmath>
#include <optional>

namespace Geoxy::Exact
{
  namespace
  {
#if GEOXY_EXACT_FMA
    /**
     *  Whether the product of two displacements' components has its roundoff exact as a double: no
     *  overflow, and the product at least 2^-968, so its 106 bits stay above 2^-1074; or a factor
     *  zero.
     *
     *  @pre           Both factors are finite, as an exact difference is.
     *  @param  first  One factor.
     *  @param  second The other.
     *  @return        Whether fma(first, second, -first second) is the roundoff exactly.
     */
    auto Safe(double first, double second) -> bool
    {
      const double product = std::abs(first * second);
      if (product >= 0x1.0p-968 && product <= 0x1.0p1020)
      {
        return true;
      }
      return std::fpclassify(first) == FP_ZERO || std::fpclassify(second) == FP_ZERO;
    }

    /**
     *  The sign of a cross product of displacements that are themselves exact, by one fused
     *  comparison.
     *
     *  With p = c d rounded, e = fma(c, d, -p) is its roundoff exactly, so c d = p + e, and a b - c
     *  d = x - e with x = a b - p. f = fma(a, b, -p) is x rounded once. Rounding is monotone and e
     *  is a double, so x <= e forces f <= e: f > e proves the cross product positive, and f < e
     *  negative. Only f = e leaves the sign open, as every zero cross product does
     *  (tests/proof/rounding.v, cross_sign).
     *
     *  @param  first  The displacement turned from, exact.
     *  @param  second The displacement turned to, exact.
     *  @return        1 or -1; 0 where the sign is open.
     */
    auto CrossSign(Vector first, Vector second) -> int
    {
      if (!Safe(first.X, second.Y) || !Safe(first.Y, second.X))
      {
        return 0;
      }
      const double product  = first.Y * second.X;
      const double roundoff = std::fma(first.Y, second.X, -product);
      const double rest     = std::fma(first.X, second.Y, -product);
      if (rest > roundoff)
      {
        return 1;
      }
      if (rest < roundoff)
      {
        return -1;
      }
      return 0;
    }
#endif
  } // namespace

  /**
   *  Where the four differences are exact and their products' roundoffs too, wherever the
   *  coordinates lie, the determinant is a cross product of doubles, one fused comparison away:
   *  points near one another have exact differences, and so do integers below 2^53. What that
   *  leaves open is summed exactly in whole numbers: two words of them where the coordinates span
   *  at most 125 bits (exact/grid.hpp), as many as the coordinates take elsewhere
   *  (exact/whole.hpp).
   */
  auto OrientAdapt(const Triangle &corners) -> int
  {
#if GEOXY_EXACT_FMA
    const Vector first  = corners.First - corners.Third;
    const Vector second = corners.Second - corners.Third;
    if (Unrounded(Roundoff(corners.First, corners.Third))
        && Unrounded(Roundoff(corners.Second, corners.Third)))
    {
      if (const int sign = CrossSign(first, second); sign != 0)
      {
        return sign;
      }
    }
#endif
    if (const std::optional<int> sign = OrientGrid(corners))
    {
      return *sign;
    }
    return OrientWhole(corners);
  }
} // namespace Geoxy::Exact
