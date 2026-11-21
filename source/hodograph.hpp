/***************************************************************************************************
 *  @file      source/hodograph.hpp
 *  @brief     The length of a polynomial curve from its hodograph: in closed form where the
 *             hodograph is straight, by Gauss-Legendre rules whose error is bounded from its roots
 *             where it is a parabola.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_HODOGRAPH_HPP
#define GEOXY_HODOGRAPH_HPP

#include "quadrature.hpp"

#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <cstddef>

namespace Geoxy
{
  /**
   *  A polynomial curve's derivative, its hodograph, in Bernstein form over the curve's parameter
   *  from 0 to 1, a factor held apart so that it multiplies the length once: a cubic Bézier curve's
   *  is 3 times the differences of its control points, a quadratic's 2 times.
   */
  struct Hodograph
  {
    std::size_t           Degree = 0;   /**< Its degree, at most 2.                    */
    std::array<Vector, 3> Points = {};  /**< Its control points, the first Degree + 1. */
    double                Factor = 1.0; /**< What they are to be multiplied by.        */
  };

  /**
   *  The most halvings a stretch's length by Measured may take, whatever its effort allows.
   */
  inline constexpr std::size_t Halved = 252;

  /**
   *  A quadratic Bézier curve's hodograph: the differences of its control points, times 2.
   *
   *  @param  quadratic The curve.
   *  @return           The hodograph, of degree 1.
   */
  auto Hodographed(const Bezier::Quadratic &quadratic) -> Hodograph;

  /**
   *  A cubic Bézier curve's hodograph: the differences of its control points, times 3.
   *
   *  @param  cubic The curve.
   *  @return       The hodograph, of degree 2.
   */
  auto Hodographed(const Bezier::Cubic &cubic) -> Hodograph;

  /**
   *  The length a polynomial curve runs over a stretch of its parameter: the integral of its
   *  hodograph's length, past 0 and 1 too.
   *
   *  A hodograph of degree 0 runs at one speed. One of degree 1, h(t) = h0 + a t, has the speed |a|
   *  sqrt((t - p)^2 + q^2), whose integral is in closed form, both differences in it rationalised
   *  so that nothing cancels (see Linear).
   *
   *  One of degree 2 is a complex quadratic h(t) = x'(t) + i y'(t), and the speed sqrt(h(t)
   *  conj(h)(t)) is analytic but at h's two roots and their conjugates. The stretch is cut where a
   *  root's real part falls, there the speed may turn back or have a cusp, and each piece is halved
   *  while the Bernstein ellipse about it that reaches the nearest root is thinner than rho = 3. On
   *  a piece of half width w, with |f| at most M on the ellipse of rho, Gauss-Legendre with n
   *  points is off by at most (64 / 15) M w rho^(-2 n) / (rho^2 - 1) (Trefethen, Approximation
   *  Theory and Approximation Practice, theorem 19.3); M is the leading coefficient's length times
   *  each root's distance from the piece's middle plus the ellipse's half axis. Each piece takes
   *  the fewest points, in pairs, that keep that within three quarters of its share of the error,
   *  the share its width's of the whole.
   *
   *  A root whose imaginary part q is so small that sqrt((t - p)^2 + q^2) differs from |t - p| by
   *  little over a piece is taken as real there: the speed is then |t - p| times what is analytic,
   *  itself analytic on the piece, as the cut keeps p off its inside. It differs from the speed by
   *  q^2 / (2 |t - p|) at most, so by at most the other root's reach times q^2 (log(1 + 2 w / q) +
   *  1) in the integral, and by 568 times that over the rule, no node of 40 points or fewer lying
   *  nearer an end than 0.00176 w: a quarter of the share bounds both.
   *
   *  The error asked is at least 2^-55 of the length, of a lower bound on it: the polygon through
   *  the points at the quarters of the stretch, each side the hodograph's exact integral by the
   *  rule of 2 points. The pieces are summed with their roundoffs. What is left is rounding in the
   *  speeds, a few steps of a double; and the roots, found by the stable quadratic formula, are
   *  taken as they are found: their own rounding is not in the bound.
   *
   *  @param  hodograph The hodograph, no coordinate of it not a number.
   *  @param  stretch   The stretch, its start at most its end.
   *  @param  effort    How far off the length may be, and how many halvings, at most Halved, may be
   *                    spent; once spent, each piece left takes the most points.
   *  @return           The length: 0 for a hodograph of nothing or a stretch of no width; without
   *                    end for any other where the stretch's width or a coordinate of the hodograph
   *                    is.
   */
  auto Measured(const Hodograph &hodograph, Range stretch, Effort effort) -> double;
} // namespace Geoxy

#endif /* GEOXY_HODOGRAPH_HPP */
