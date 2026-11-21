/***************************************************************************************************
 *  @file      source/aberth.hpp
 *  @brief     The complex roots of a polynomial in Bernstein form, by the Aberth-Ehrlich iteration.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_ABERTH_HPP
#define GEOXY_ABERTH_HPP

#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <cstddef>
#include <span>

namespace Geoxy::Aberth
{
  /**
   *  The most roots a polynomial here has: those of a rational piece's speed, of degree 2 MaxDegree
   *  - 2, and a few more.
   */
  inline constexpr std::size_t Most = 2 * static_cast<std::size_t>(MaxDegree);

  /**
   *  How many rounds a polynomial's roots take at most where nothing else is asked: Aberth's
   *  iteration settles simple roots in a dozen, and a root of multiplicity m in some 50 m at worst.
   */
  inline constexpr int Rounds = 128;

  /**
   *  A polynomial's roots, each a complex number held as a vector: its x the real part, its y the
   *  imaginary.
   */
  struct Roots
  {
    std::size_t              Count  = 0;  /**< How many.                   */
    std::array<Vector, Most> Values = {}; /**< The roots, the first Count. */
  };

  /**
   *  The roots of a polynomial in Bernstein form over 0 to 1, with complex coefficients.
   *
   *  One of degree 1 or 2 is solved in closed form, in w = t / (1 - t), where it is (1 - t)^m times
   *  a polynomial whose coefficients are the Bernstein ones times binomials, no difference of them
   *  taken: by the stable quadratic formula, and each root w taken back as w / (1 + w).
   *
   *  Any other is solved by Aberth and Ehrlich's iteration, which moves every approximation at once
   *  by Newton's step corrected for the others, z -= N / (1 - N sum 1 / (z - w)), N = p(z) / p'(z),
   *  and converges cubically to simple roots. It starts on circles about the middle, 1 / 2, whose
   *  radii the upper convex hull of the points (k, log |c_k|) gives for the coefficients c_k in
   *  powers of 2 t - 1, Bini's choice, so that roots of sizes far apart are each started near their
   *  own. It evaluates p and p' in Bernstein form, Horner's rule in t / (1 - t) where |t| is at
   *  most |1 - t|, and in (1 - t) / t elsewhere, each ratio at most 1: so nothing overflows, and
   *  near an end p is told by the coefficients there, to their own digits, where a polynomial whose
   *  coefficients range far, as a weight's beside a heavy one does, has lost them in powers. An
   *  approximation is settled once p is within its own rounding there, 4 m steps of a double of the
   *  sum of the terms' lengths, or once its step is within a step of a double of it.
   *
   *  What it gives are roots of the coefficients as they are, to within their own rounding, and
   *  where an approximation has not settled once the rounds are spent, as far as it got. A root
   *  without end, where the polynomial's degree in t is less than its form's, is far off: 2^500 in
   *  closed form, and for a higher degree as many fewer as the powers at the top that are 0, or
   *  where they are 0 to within their rounding alone, an approximation gone off far.
   *
   *  @param  coefficients The Bernstein coefficients, two at least and at most Most + 1, not all 0.
   *  @param  rounds       How many rounds the approximations may take.
   *  @return              The roots, as many as the form's degree at most.
   */
  auto Solved(std::span<const Vector> coefficients, int rounds) -> Roots;
} // namespace Geoxy::Aberth

#endif /* GEOXY_ABERTH_HPP */
