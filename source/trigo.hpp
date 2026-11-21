/***************************************************************************************************
 *  @file      source/trigo.hpp
 *  @brief     Geoxy's own sines, cosines, tangents and arc tangents, within about 0.56 steps of a
 *             double and the same bits on every machine, with degrees reduced exactly.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

/**
 *  Ported from trigo, the author's library of the same name, to what the library calls. Each
 *  function is worked from sums, products, quotients, square roots and fused multiply-adds alone,
 *  every one rounded as IEEE 754 says, so a machine that rounds them so gives the same bits: arm64
 *  and x86-64 alike, where the C library's functions differ from vendor to vendor. A machine with
 *  no fused multiply-add in its hardware gets std::fma from its C library, which rounds as it must
 *  but slowly. Direction and Heading keep their own correctly rounded routines.
 */

#ifndef GEOXY_TRIGO_HPP
#define GEOXY_TRIGO_HPP

#include "geoxy/vector.hpp"

namespace Geoxy::Trigo
{
  /**
   *  The sine.
   *
   *  @param  radians The angle.
   *  @return         Its sine, within 0.56 steps of a double; not a number for an infinite angle.
   */
  [[nodiscard]] auto Sin(double radians) -> double;

  /**
   *  The cosine.
   *
   *  @param  radians The angle.
   *  @return         Its cosine, within 0.56 steps of a double; not a number for an infinite angle.
   */
  [[nodiscard]] auto Cos(double radians) -> double;

  /**
   *  The tangent.
   *
   *  @param  radians The angle.
   *  @return         Its tangent, within 0.59 steps of a double; not a number for an infinite
   *                  angle.
   */
  [[nodiscard]] auto Tan(double radians) -> double;

  /**
   *  The cosine and the sine together, for about the price of one.
   *
   *  @param  radians The angle.
   *  @return         The cosine and the sine, each as Cos and Sin give it.
   */
  [[nodiscard]] auto Way(double radians) -> Vector;

  /**
   *  The cosine and the sine of an angle in degrees, reduced by whole quarter turns exactly: 90
   *  degrees is (0, 1) and 30 has a sine of 1/2, which no sine of the radians can give.
   *
   *  @param  degrees The angle.
   *  @return         The cosine and the sine, within 0.56 steps of a double; +0 where either is
   *                  zero but the sine of a negative angle, which is -0.
   */
  [[nodiscard]] auto WayDegrees(double degrees) -> Vector;

  /**
   *  The arc tangent.
   *
   *  @param  tangent The tangent.
   *  @return         Its angle, in [-pi / 2, pi / 2], within 0.51 steps of a double.
   */
  [[nodiscard]] auto Atan(double tangent) -> double;

  /**
   *  The angle of a displacement, as IEEE 754 sets it for zeros and infinities.
   *
   *  @param  rise The rise.
   *  @param  run  The run.
   *  @return      Its angle, in [-pi, pi], within 0.51 steps of a double.
   */
  [[nodiscard]] auto Atan2(double rise, double run) -> double;

  /**
   *  The arc sine.
   *
   *  @param  sine The sine.
   *  @return      Its angle, in [-pi / 2, pi / 2], within 0.53 steps of a double; not a number past
   *               1 either way.
   */
  [[nodiscard]] auto Asin(double sine) -> double;

  /**
   *  The arc cosine.
   *
   *  @param  cosine The cosine.
   *  @return        Its angle, in [0, pi], within 0.52 steps of a double; not a number past 1
   *                 either way.
   */
  [[nodiscard]] auto Acos(double cosine) -> double;

  /**
   *  The length of a displacement, worked from its squares as a pair and one Newton step, scaled by
   *  a power of two, so that no square overflows or underflows.
   *
   *  @param  run  The run.
   *  @param  rise The rise.
   *  @return      sqrt(run^2 + rise^2), infinite when either is, even beside a NaN.
   */
  [[nodiscard]] auto Hypot(double run, double rise) -> double;
  /**
   *  The binary logarithm, exact at powers of two.
   *
   *  @param  value The value.
   *  @return       log2 of it, within about 0.6 steps of a double; -infinity at 0, not a number
   *                below it.
   */
  [[nodiscard]] auto Log2(double value) -> double;

  /**
   *  2 to a power, exact at whole powers.
   *
   *  @param  power The power.
   *  @return       2^power, within about 0.8 steps of a double.
   */
  [[nodiscard]] auto Exp2(double power) -> double;

  /**
   *  A positive base to a power, as 2 to the power times the base's binary logarithm: off by the
   *  product's rounding times its size, a few steps of a double where it is near 10.
   *
   *  @param  base  The base, not negative.
   *  @param  power The power.
   *  @return       base^power.
   */
  [[nodiscard]] auto Pow(double base, double power) -> double;

  /**
   *  The inverse hyperbolic sine, ln(x + sqrt(x^2 + 1)), the sum worked as a pair.
   *
   *  @param  value x.
   *  @return       Its inverse hyperbolic sine, within about 0.6 steps of a double.
   */
  [[nodiscard]] auto Asinh(double value) -> double;

  /**
   *  The cube root: Newton's steps from a guess of the bits, the last from the cube's exact
   *  remainder.
   *
   *  @param  value The value.
   *  @return       Its cube root, within about 0.51 steps of a double.
   */
  [[nodiscard]] auto Cbrt(double value) -> double;
} // namespace Geoxy::Trigo

#endif /* GEOXY_TRIGO_HPP */
