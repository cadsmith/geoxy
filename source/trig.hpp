/***************************************************************************************************
 *  @file      source/trig.hpp
 *  @brief     The sines, arc tangents, lengths, logarithms and roots the library calls outside its
 *             correctly rounded Direction and Heading: its own, source/trigo, or the C library's.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

/**
 *  Built with GEOXY_TRIGO, every function here is source/trigo's: within 0.56 steps of a double,
 *  degrees reduced exactly, and the same bits on every machine. Without it, each is the <cmath>
 *  call it stands for, as the C library at hand computes it.
 */

#ifndef GEOXY_TRIG_HPP
#define GEOXY_TRIG_HPP

#include "geoxy/angle.hpp"
#include "geoxy/vector.hpp"

#ifdef GEOXY_TRIGO
#  include "trigo.hpp"
#else
#  include <cmath>
#endif

namespace Geoxy::Trig
{
  /**
   *  The sine.
   *
   *  @param  radians The angle.
   *  @return         Its sine.
   */
  [[nodiscard]] inline auto Sin(double radians) -> double
  {
#ifdef GEOXY_TRIGO
    return Trigo::Sin(radians);
#else
    return std::sin(radians);
#endif
  }

  /**
   *  The cosine.
   *
   *  @param  radians The angle.
   *  @return         Its cosine.
   */
  [[nodiscard]] inline auto Cos(double radians) -> double
  {
#ifdef GEOXY_TRIGO
    return Trigo::Cos(radians);
#else
    return std::cos(radians);
#endif
  }

  /**
   *  The tangent.
   *
   *  @param  radians The angle.
   *  @return         Its tangent.
   */
  [[nodiscard]] inline auto Tan(double radians) -> double
  {
#ifdef GEOXY_TRIGO
    return Trigo::Tan(radians);
#else
    return std::tan(radians);
#endif
  }

  /**
   *  The arc tangent.
   *
   *  @param  value The tangent.
   *  @return       Its angle, in [-pi / 2, pi / 2].
   */
  [[nodiscard]] inline auto Atan(double value) -> double
  {
#ifdef GEOXY_TRIGO
    return Trigo::Atan(value);
#else
    return std::atan(value);
#endif
  }

  /**
   *  The angle of a displacement.
   *
   *  @param  y The rise.
   *  @param  x The run.
   *  @return   Its angle, in [-pi, pi].
   */
  [[nodiscard]] inline auto Atan2(double y, double x) -> double
  {
#ifdef GEOXY_TRIGO
    return Trigo::Atan2(y, x);
#else
    return std::atan2(y, x);
#endif
  }

  /**
   *  The arc sine.
   *
   *  @param  value The sine.
   *  @return       Its angle, in [-pi / 2, pi / 2].
   */
  [[nodiscard]] inline auto Asin(double value) -> double
  {
#ifdef GEOXY_TRIGO
    return Trigo::Asin(value);
#else
    return std::asin(value);
#endif
  }

  /**
   *  The arc cosine.
   *
   *  @param  value The cosine.
   *  @return       Its angle, in [0, pi].
   */
  [[nodiscard]] inline auto Acos(double value) -> double
  {
#ifdef GEOXY_TRIGO
    return Trigo::Acos(value);
#else
    return std::acos(value);
#endif
  }

  /**
   *  The length of a displacement.
   *
   *  @param  x The run.
   *  @param  y The rise.
   *  @return   sqrt(x^2 + y^2), infinite when either is.
   */
  [[nodiscard]] inline auto Hypot(double x, double y) -> double
  {
#ifdef GEOXY_TRIGO
    return Trigo::Hypot(x, y);
#else
    return std::hypot(x, y);
#endif
  }

  /**
   *  The binary logarithm.
   *
   *  @param  value The value.
   *  @return       log2 of it.
   */
  [[nodiscard]] inline auto Log2(double value) -> double
  {
#ifdef GEOXY_TRIGO
    return Trigo::Log2(value);
#else
    return std::log2(value);
#endif
  }

  /**
   *  2 to a power.
   *
   *  @param  power The power.
   *  @return       2^power.
   */
  [[nodiscard]] inline auto Exp2(double power) -> double
  {
#ifdef GEOXY_TRIGO
    return Trigo::Exp2(power);
#else
    return std::exp2(power);
#endif
  }

  /**
   *  A base, not negative, to a power.
   *
   *  @param  base  The base.
   *  @param  power The power.
   *  @return       base^power.
   */
  [[nodiscard]] inline auto Pow(double base, double power) -> double
  {
#ifdef GEOXY_TRIGO
    return Trigo::Pow(base, power);
#else
    return std::pow(base, power);
#endif
  }

  /**
   *  The inverse hyperbolic sine.
   *
   *  @param  value The value.
   *  @return       asinh of it.
   */
  [[nodiscard]] inline auto Asinh(double value) -> double
  {
#ifdef GEOXY_TRIGO
    return Trigo::Asinh(value);
#else
    return std::asinh(value);
#endif
  }

  /**
   *  The cube root.
   *
   *  @param  value The value.
   *  @return       Its cube root.
   */
  [[nodiscard]] inline auto Cbrt(double value) -> double
  {
#ifdef GEOXY_TRIGO
    return Trigo::Cbrt(value);
#else
    return std::cbrt(value);
#endif
  }

  /**
   *  The cosine and sine of an angle in radians, as a unit vector.
   *
   *  @param  radians The angle.
   *  @return         The cosine and the sine.
   */
  [[nodiscard]] inline auto Way(double radians) -> Vector
  {
#ifdef GEOXY_TRIGO
    return Trigo::Way(radians);
#else
    return {.X = std::cos(radians), .Y = std::sin(radians)};
#endif
  }

  /**
   *  The cosine and sine of an angle, not certified as Direction is: of its degrees when it is held
   *  in degrees and source/trigo is built, so 90 degrees is on the axis; else of its radians.
   *
   *  @param  angle The angle.
   *  @return       The cosine and the sine.
   */
  [[nodiscard]] inline auto Way(Angle angle) -> Vector
  {
#ifdef GEOXY_TRIGO
    if (angle.HeldInDegrees())
    {
      return Trigo::WayDegrees(angle.InDegrees());
    }
#endif
    return Way(angle.InRadians());
  }

  /**
   *  The cosine and sine at the end of a sweep, not certified as Direction is: of the degrees
   *  summed when both are held in degrees and source/trigo is built, so a sweep from 0 to 90
   *  degrees ends on the axis; else of the radians summed.
   *
   *  @param  start The start.
   *  @param  sweep The sweep.
   *  @return       The cosine and the sine.
   */
  [[nodiscard]] inline auto Way(Angle start, Angle sweep) -> Vector
  {
#ifdef GEOXY_TRIGO
    if (start.HeldInDegrees() && sweep.HeldInDegrees())
    {
      return Trigo::WayDegrees(start.InDegrees() + sweep.InDegrees());
    }
#endif
    return Way(start.InRadians() + sweep.InRadians());
  }
} // namespace Geoxy::Trig

#endif /* GEOXY_TRIG_HPP */
