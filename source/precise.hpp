/***************************************************************************************************
 *  @file      source/precise.hpp
 *  @brief     Angles to any precision, for the rare answer the fast path cannot round: sines,
 *             cosines and arctangents in whole-number fixed point, and which side of a halfway
 *             point between doubles they fall on.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_PRECISE_HPP
#define GEOXY_PRECISE_HPP

#include "scratch.hpp"
#include "twofold.hpp"

#include <array>
#include <cstdint>
#include <functional>
#include <span>

/**
 *  Real numbers held as whole numbers of 2^-bits, in 32-bit words, as precise as asked: the slow
 *  path behind correctly rounded angles. Each value is computed to within 2^(20 - bits) of the
 *  exact one, and Side doubles bits until the halfway point is cleared, which it always is: no
 *  sine, cosine or arctangent of a double that the fast path hands on is a halfway point between
 *  doubles (Lindemann's theorem in radians, Niven's in degrees).
 */
namespace Geoxy::Precise
{
  /**
   *  A real number m 2^-bits, bits given by its use: m in 32-bit words, least significant first.
   */
  struct Fixed
  {
    List<std::uint32_t> Words    = {};    /**< m's words, least significant first. */
    bool                Negative = false; /**< Whether the number is below zero.   */
  };

  /**
   *  A computation to any precision: the value to within 2^(20 - bits) for the bits asked.
   */
  using Computed = std::function<Fixed(int bits)>;

  /**
   *  pi by a table of 4096 bits, past them by Machin's formula, 16 atan(1/5) - 4 atan(1/239).
   *
   *  @param  bits The precision, to 2^13.
   *  @return      pi, within 2^(16 - bits).
   */
  auto Pi(int bits) -> Fixed;

  /**
   *  How many times Side has run on this thread: how often the slow path is taken, for the tests
   *  and the measures.
   *
   *  @return The count.
   */
  auto Slowed() -> std::uint64_t &;

  /**
   *  The double nearest a value known as a pair of doubles: from the pair alone where its error
   *  cannot carry it across a halfway point, else by Side on the value computed to any precision.
   *
   *  @param  pair  The pair, its small part at most an ulp of its big part.
   *  @param  error How far the pair may be from the value, relatively.
   *  @param  value The value, to any precision.
   *  @return       The double nearest the value, subnormal ones included.
   */
  auto Nearest(Twofold::Two pair, double error, const Computed &value) -> double;

  /**
   *  Which side of a halfway point between doubles a computed value lies on.
   *
   *  @param  value The value, computed to any precision.
   *  @param  base  The double the halfway point starts from.
   *  @param  step  The step to the next double, signed, a power of two.
   *  @return       1 when the value is above base + step / 2, -1 when below; never 0.
   */
  auto Side(const Computed &value, double base, double step) -> int;

  /**
   *  The cosine and sine of degrees given as a sum of two doubles, within a quarter turn and a hair
   *  of 0.
   *
   *  @param  degrees The degrees, an exact sum of the pair.
   *  @param  quarter Quarter turns to add, 0 to 3.
   *  @param  bits    The precision.
   *  @return         The cosine, then the sine, of the degrees and the quarter turns.
   */
  auto Degree(Twofold::Two degrees, int quarter, int bits) -> List<Fixed>;

  /**
   *  The cosine and sine of radians, any double, reduced by quarter turns exactly first.
   *
   *  @param  radians The radians, finite.
   *  @param  bits    The precision.
   *  @return         The cosine, then the sine.
   */
  auto Radian(double radians, int bits) -> List<Fixed>;

  /**
   *  The heading in radians, in (-pi, pi], of a displacement whose components are exact sums of
   *  doubles.
   *
   *  @param  y    Along the y axis, as doubles that sum to it.
   *  @param  x    Along the x axis, likewise; not both 0.
   *  @param  bits The precision.
   *  @return      The heading; a -0 component counts as +0.
   */
  auto Arc(std::span<const double> y, std::span<const double> x, int bits) -> Fixed;

  /**
   *  A cosine and a sine, each a pair of doubles within 2^-68 of the exact value, relatively.
   */
  struct Way
  {
    Twofold::Two Cosine = {}; /**< The cosine. */
    Twofold::Two Sine   = {}; /**< The sine.   */
  };

  /**
   *  The nearest doubles to the cosine and sine of degrees, from their pairs where those settle,
   *  else to any precision: the slow way, which the fast path calls in this file to keep its own
   *  frame in registers.
   *
   *  @param  way     The pairs.
   *  @param  degrees The degrees, as for Degree.
   *  @param  quarter Quarter turns to add, 0 to 3.
   *  @return         The cosine, then the sine.
   */
  auto NearestDegree(Way way, Twofold::Two degrees, int quarter) -> std::array<double, 2>;

  /**
   *  The nearest doubles to the cosine and sine of radians, likewise.
   *
   *  @param  way     The pairs.
   *  @param  radians The radians, finite.
   *  @return         The cosine, then the sine.
   */
  auto NearestRadian(Way way, double radians) -> std::array<double, 2>;

  /**
   *  The nearest double to the sine of tiny degrees held as a pair, likewise.
   *
   *  @param  near    The sine, within 2^-50 of the exact one, relatively.
   *  @param  degrees The degrees, below 2^-500.
   *  @return         The sine.
   */
  auto NearestSine(double near, Twofold::Two degrees) -> double;

  /**
   *  The nearest double to a heading, likewise.
   *
   *  @param  pair The heading as a pair, within 2^-68 of the exact one, relatively.
   *  @param  y    As for Arc.
   *  @param  x    As for Arc.
   *  @return      The heading.
   */
  auto NearestArc(Twofold::Two pair, std::span<const double> y, std::span<const double> x)
      -> double;
} // namespace Geoxy::Precise

#endif /* GEOXY_PRECISE_HPP */
