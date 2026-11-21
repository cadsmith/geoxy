/***************************************************************************************************
 *  @file      source/crossings/bows.hpp
 *  @brief     Where two polylines' bulges' circles cross, from the bulges' own numbers: their
 *             radical line met with the smaller circle, in double-double with error bounds.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CROSSINGS_BOWS_HPP
#define GEOXY_CROSSINGS_BOWS_HPP

#include "sure.hpp"

#include "geoxy/point.hpp"

#include <array>
#include <optional>

namespace Geoxy::Crossings
{
  /**
   *  A polyline's bulge: its two vertices and its bulge, its circle the exact rational one through
   *  both.
   */
  struct Belly
  {
    Point  Start = {};  /**< Its first vertex, P0.             */
    Point  End   = {};  /**< Its second vertex, P1, not P0.    */
    double Bulge = 0.0; /**< Its bulge, b, not under Straight. */
  };

  /**
   *  Where two bulges' circles cross, each crossing's coordinates in double-double with a bound on
   *  how far each is from the exact one.
   */
  struct Crossbow
  {
    std::array<std::array<Sure::Number, 2>, 2> At = {}; /**< The two crossings, x and y. */
  };

  /**
   *  Whether the second of two bulges has the smaller circle, its radius worked out roughly:
   *  Crossbows meets their radical line with the smaller, where its root is clearest.
   *
   *  @param  bellies The two bulges.
   *  @return         Whether it does; not where they are as large.
   */
  auto Narrower(const std::array<Belly, 2> &bellies) -> bool;

  /**
   *  Where two bulges' circles cross near a point, the double nearest: Crossbows's radical line met
   *  with the first circle in plain doubles, in the frame of the first bulge's start, for a guess
   *  at the crossing nearer the point; then one Newton step from the guess on both circles'
   *  quadratics, each worked out exactly but for sure.hpp's bounds. A quadratic is its value at the
   *  guess, its gradient there, and 4 b times the step's square, nothing more: so the step's error
   *  is bounded by the step's own size squared, and the crossing settled where that and the sure
   *  numbers' bounds leave no halfway point near, and where no other crossing lies nearer the
   *  point.
   *
   *  @param  bellies The two bulges, the smaller circle's first, as Narrower tells.
   *  @param  near    The point.
   *  @return         The crossing nearer the point, the double nearest; nothing where it is not
   *                  settled so.
   */
  auto Hopped(const std::array<Belly, 2> &bellies, Point near) -> std::optional<Point>;

  /**
   *  Where two bulges' circles cross, worked from the bulges' own numbers with no center and no
   *  radius.
   *
   *  A bulge from P0 to P1 of bulge b has its circle the points X where 4 b (X - P0).(X - P1) - 2
   *  (1 - b^2) (X - P0).n is zero, n the chord turned a quarter to the left: with Y = X - P0 and c
   *  the chord, 4 b |Y|^2 = Y.q for q = 4 b c + 2 (1 - b^2) n, a circle through Y = 0 of center q /
   *  (8 b). The other bulge, from Q0 of chord d and bulge e, s = Q0 - P0, is 4 e |Y - s|^2 = (Y
   *  - s).p for p = 4 e d + 2 (1 - e^2) n_d. Taking e times the first from b times the second
   *    leaves their radical line, A.Y = B, A = 8 b e s - e q + b p and B = 4 b e |s|^2 + b s.p; A
   *    is 8 b e times the step between the centers. The line meets the first circle where Y = (8 b
   *    B A + (A'.q +- sqrt(D)) A') / (8 b |A|^2), A' the step turned a quarter, D = |A|^2 |q|^2 -
   *    (A.q - 8 b B)^2: D over |A|^2 |q|^2 is the square of the half chord over the first's radius,
   *    which for the smaller circle is well apart from nothing wherever the circles cross at any
   *    slant.
   *
   *  @param  bellies The two bulges, the one met with the radical line first: the smaller, as
   *                  Narrower tells.
   *  @return         Both crossings, the root added first; nothing where the radicand's bound
   *                  reaches zero, as where the circles do not cross. A number past every double,
   *                  or a quotient over a divisor whose bound reaches zero, has a bound not finite,
   *                  which settles nothing.
   */
  auto Crossbows(const std::array<Belly, 2> &bellies) -> std::optional<Crossbow>;
} // namespace Geoxy::Crossings

#endif /* GEOXY_CROSSINGS_BOWS_HPP */
