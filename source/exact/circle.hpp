/***************************************************************************************************
 *  @file      source/exact/circle.hpp
 *  @brief     The exact circle test, its first stage inline.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_EXACT_CIRCLE_HPP
#define GEOXY_EXACT_CIRCLE_HPP

#include "exact/bounds.hpp"
#include "exact/triangle.hpp"
#include "exact/two.hpp"

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <cmath>

/**
 *  The first stage is inline for the reason the orientation test's is (exact/orient.hpp).
 *  Shewchuk's is already a single test, |det| > bound * permanent, kept with one term added for
 *  underflow. A product that underflows is off by up to half the smallest double, u_S, not by a
 *  part of itself, and the error reaches the determinant times the factors after it, each at most
 *  twice the largest lift L: in all a corner's products and its permanent's are off by at most 15
 *  u_S L + 4 u_S, which 32 u_S (S + 1) covers, S the sum of the lifts, at least L. Overflow leaves
 *  the permanent infinite or NaN, since rounding keeps it at least the determinant's magnitude term
 *  by term, and the strict comparison fails on either. So the test decides no wrong sign for any
 *  doubles. tests/proof/first.py checks it with Z3.
 */

namespace Geoxy::Exact
{
  /**
   *  The circle test's first stage approximation of the determinant, and the permanent its error is
   *  relative to.
   */
  struct Lifted
  {
    double Determinant = 0.0; /**< The determinant, rounded.               */
    double Permanent   = 0.0; /**< The sum of the magnitudes of its terms. */
    double Lifts       = 0.0; /**< The sum of the three lifts.             */
  };

  /**
   *  The circle test's first stage.
   *
   *  @param  circle The corners the circle passes through.
   *  @param  point  The point asked about.
   *  @return        The determinant of the lifted points rounded, and its permanent.
   */
  inline auto Lift(const Triangle &circle, Point point) -> Lifted
  {
    const Vector a = circle.First - point;
    const Vector b = circle.Second - point;
    const Vector c = circle.Third - point;

    const double bxcy  = b.X * c.Y;
    const double cxby  = c.X * b.Y;
    const double alift = (a.X * a.X) + (a.Y * a.Y);
    const double cxay  = c.X * a.Y;
    const double axcy  = a.X * c.Y;
    const double blift = (b.X * b.X) + (b.Y * b.Y);
    const double axby  = a.X * b.Y;
    const double bxay  = b.X * a.Y;
    const double clift = (c.X * c.X) + (c.Y * c.Y);

    return {
        .Determinant = (alift * (bxcy - cxby)) + (blift * (cxay - axcy)) + (clift * (axby - bxay)),
        .Permanent   = ((std::abs(bxcy) + std::abs(cxby)) * alift)
                     + ((std::abs(cxay) + std::abs(axcy)) * blift)
                     + ((std::abs(axby) + std::abs(bxay)) * clift),
        .Lifts       = (alift + blift) + clift,
    };
  }

  /**
   *  The circle test's stages past the first: a compensated second stage, then Shewchuk's
   *  incircleadapt's third and last (exact/circle.cpp). Where a coordinate is neither zero nor
   *  between 2^-150 and 2^150, past what their products hold, whole numbers instead: on a grid
   *  where the coordinates fit one (exact/grid.hpp), else of any size (exact/whole.hpp).
   *
   *  @param  circle    The corners the circle passes through.
   *  @param  point     The point asked about.
   *  @param  permanent The first stage's permanent.
   *  @return           1, 0 or -1.
   */
  auto CircleAdapt(const Triangle &circle, Point point, double permanent) -> int;

  /**
   *  Where a point lies from the circle through three corners, exactly: the sign of the determinant
   *  of the points lifted onto a paraboloid.
   *
   *  Shewchuk's adaptive incircle, its second stage compensated, its result the sign alone, exact
   *  for every finite double: past his limits, products that would underflow or overflow, the sign
   *  is summed in whole numbers.
   *
   *  @param  circle The corners the circle passes through.
   *  @param  point  The point asked about.
   *  @return        1 inside the circle when its corners turn counterclockwise, -1 outside, 0 on
   *                 it; the signs swap when they turn clockwise; 0 where a coordinate is NaN or
   *                 infinite.
   */
  [[nodiscard]] inline auto InCircle(const Triangle &circle, Point point) -> int
  {
    const Lifted lifted = Lift(circle, point);
    if (std::abs(lifted.Determinant)
        > (CircleBoundA * lifted.Permanent) + (CircleBoundUnderflow * (lifted.Lifts + 1.0)))
        [[likely]]
    {
      return Sign(lifted.Determinant);
    }
    return CircleAdapt(circle, point, lifted.Permanent);
  }
} // namespace Geoxy::Exact

#endif /* GEOXY_EXACT_CIRCLE_HPP */
