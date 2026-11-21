/***************************************************************************************************
 *  @file      source/exact/circle.cpp
 *  @brief     The circle test's stages past the first.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/circle.hpp"

#include "exact/bounds.hpp"
#include "exact/grid.hpp"
#include "exact/roundoff.hpp"
#include "exact/triangle.hpp"
#include "exact/two.hpp"
#include "exact/whole.hpp"

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>

/**
 *  The second stage is compensated (Compensated, below, with its proof): the determinant of the
 *  rounded differences to within a few u^2 of the permanent, where Shewchuk's incircleadapt sums it
 *  exactly. The third is his first order correction by the differences' roundoffs, its bound
 *  widened by the estimate's error.
 *
 *  What they leave undecided is summed exactly, in whole numbers: in 256 bits where the coordinates
 *  lie on one grid (exact/grid.hpp), as snapped and integer places do, else in as many words as
 *  they take (exact/whole.hpp).
 */

namespace Geoxy::Exact
{
  namespace
  {
    /**
     *  What the later stages know of one corner of the circle, seen from the point.
     */
    struct Corner
    {
      Vector Difference = {}; /**< The corner minus the point, rounded.    */
      Vector Roundoff   = {}; /**< What rounding took from the difference. */
    };

    /**
     *  A corner and the two after it, in turn.
     */
    struct Turn
    {
      const Corner *Own    = nullptr; /**< The corner whose terms are summed. */
      const Corner *Next   = nullptr; /**< The corner after it.               */
      const Corner *Before = nullptr; /**< The corner before it.              */
    };

    /**
     *  A corner's term of the third stage: the first order correction by the differences'
     *  roundoffs.
     *
     *  @param  turn The corner and the two after it.
     *  @return      Its term, rounded as Shewchuk rounds it.
     */
    auto Correction(const Turn &turn) -> double
    {
      const Vector a     = turn.Own->Difference;
      const Vector b     = turn.Next->Difference;
      const Vector c     = turn.Before->Difference;
      const Vector aTail = turn.Own->Roundoff;
      const Vector bTail = turn.Next->Roundoff;
      const Vector cTail = turn.Before->Roundoff;
      return (((a.X * a.X) + (a.Y * a.Y))
              * (((b.X * cTail.Y) + (c.Y * bTail.X)) - ((b.Y * cTail.X) + (c.X * bTail.Y))))
           + (2.0 * ((a.X * aTail.X) + (a.Y * aTail.Y)) * ((b.X * c.Y) - (b.Y * c.X)));
    }

    /**
     *  A value as a double and a far smaller correction, which sum to it but for a few u^2 of it.
     */
    struct Doubled
    {
      double High = 0.0; /**< The value, rounded.          */
      double Low  = 0.0; /**< What rounding took, rounded. */
    };

    /**
     *  A difference's lift, compensated: x^2 + y^2 to within 4u^2(1 + 6u) of it.
     *
     *  @param  difference The difference.
     *  @return            The lift, its low part at most 2u(1 + 5u) of it.
     */
    auto LiftNear(Vector difference) -> Doubled
    {
      const Two x    = Square(difference.X);
      const Two y    = Square(difference.Y);
      const Two high = TwoSum(x.Big, y.Big);
      return {.High = high.Big, .Low = (x.Small + y.Small) + high.Small};
    }

    /**
     *  The cross product of two differences, compensated: next.X before.Y - next.Y before.X to
     *  within 4u^2(1 + 6u) of N = |next.X before.Y| + |next.Y before.X|.
     *
     *  @param  next   The difference turned from.
     *  @param  before The difference turned to.
     *  @return        The minor, its low part at most 2u(1 + 5u) N.
     */
    auto MinorNear(Vector next, Vector before) -> Doubled
    {
      const Two first  = TwoProduct(next.X, before.Y);
      const Two second = TwoProduct(next.Y, before.X);
      const Two high   = TwoDiff(first.Big, second.Big);
      return {.High = high.Big, .Low = (first.Small - second.Small) + high.Small};
    }

    /**
     *  A lift times a minor, compensated: the high parts' product exactly, the cross terms added to
     *  its roundoff, the low parts' product dropped. No fused product: the estimate is the same
     *  double on every platform.
     *
     *  @param  lift  The lift.
     *  @param  minor The minor.
     *  @return       Their product, to within 25u^2(1 + 10u) L N.
     */
    auto TermNear(Doubled lift, Doubled minor) -> Doubled
    {
      const Two high = TwoProduct(lift.High, minor.High);
      return {.High = high.Big,
              .Low  = high.Small + ((lift.High * minor.Low) + (lift.Low * minor.High))};
    }

    /**
     *  The determinant of the lifted rounded differences, compensated: the second stage.
     *
     *  Write u = 2^-53, a corner's rounded difference from the point (x, y), the other two b and c,
     *  its lift L = x^2 + y^2, its minor M = bx cy - by cx and N = |bx cy| + |by cx|. P' = sum L N
     *  is within (1 + u)^7 of the permanent P the first stage rounds, and the precondition keeps
     *  every product here from underflowing or overflowing.
     *
     *  Lift: x^2 = h1 + l1, y^2 = h2 + l2 and h1 + h2 = H + l3 exactly, the l's at most 2u(1 + 2u)
     *  L in all, so Lo = (l1 + l2) + l3, rounded twice, is off by at most 4u^2(1 + 6u) L and is
     *  itself at most 2u(1 + 5u) L. The minor alike, against N. Term: the high parts' product is p
     *  + e exactly; Lo Mlo is dropped, the lift's and minor's errors carried, and e + (H Mlo + Lo
     *  Mh) rounded four times: at most (4 + 4 + 4 + 8 + 5)u^2(1 + 10u) L N off. Sum: the three p
     *  sum exactly to a double and two roundoffs of at most 2u(1 + 9u) P'; the five low parts, at
     *  most 7u(1 + 11u) P', are rounded three deep, at most 21u^2(1 + 15u) P' off. So the estimate
     *  E is within u |E| + 47u^2 P of D, the determinant of the rounded differences.
     *
     *  Where the differences are exact, D is the determinant, and |E| >= 64u^2 P decides it, since
     *  (1 - u)^2 64 > 47: CircleBoundCompensated. Where they are not, each exact lift and minor is
     *  within (2u + u^2) of the rounded one's, so the determinant is within ((1 + u)^4 - 1) P' <=
     *  (4u + 40u^2) P of D, and |E| >= (4 + 128u)u P decides it, since (1 - u)^2 (4 + 128u) > 4 +
     *  87u: CircleBoundB. The third stage adds Shewchuk's first order corrections: the terms of
     *  second order are at most 7u^2(1 + 8u) P and the corrections, nine roundings deep over at
     *  most 4u P', at most 36u^2(1 + 10u) P' off, so with E's error the result C is within u |E| +
     *  u |C| + 91u^2 P, which (44 + 64)u^2 P + 3u |E|, the third stage's bound with
     *  CircleBoundCompensated added, covers. tests/proof/circle.py checks every step with Z3.
     *
     *  @param  a The first corner.
     *  @param  b The second.
     *  @param  c The third.
     *  @return   E.
     */
    auto Compensated(const Corner &a, const Corner &b, const Corner &c) -> double
    {
      const Doubled first = TermNear(LiftNear(a.Difference), MinorNear(b.Difference, c.Difference));
      const Doubled second =
          TermNear(LiftNear(b.Difference), MinorNear(c.Difference, a.Difference));
      const Doubled third = TermNear(LiftNear(c.Difference), MinorNear(a.Difference, b.Difference));
      const Two     once  = TwoSum(first.High, second.High);
      const Two     twice = TwoSum(once.Big, third.High);
      const double  low   = ((first.Low + second.Low) + third.Low) + (once.Small + twice.Small);
      return twice.Big + low;
    }

    /**
     *  Whether every coordinate of the four points is zero or from 2^-150 up to 2^150 in magnitude,
     *  where no product of four differences of coordinates, or of their roundoffs, underflows or
     *  overflows.
     *
     *  @param  circle The corners the circle passes through.
     *  @param  point  The point asked about.
     *  @return        Whether the stages' arithmetic holds its products.
     */
    auto Held(const Triangle &circle, Point point) -> bool
    {
      constexpr unsigned Range = 150;
      return std::max({Rung(circle.First.X, Range), Rung(circle.First.Y, Range),
                       Rung(circle.Second.X, Range), Rung(circle.Second.Y, Range),
                       Rung(circle.Third.X, Range), Rung(circle.Third.Y, Range),
                       Rung(point.X, Range), Rung(point.Y, Range)})
           < std::uint64_t{2} * Range;
    }

    /**
     *  The circle test's sign, exactly: in 256 bits where the coordinates lie on one grid, else in
     *  as many words as they take.
     *
     *  @param  circle The corners the circle passes through.
     *  @param  point  The point asked about.
     *  @return        1, 0 or -1.
     */
    auto CircleExactly(const Triangle &circle, Point point) -> int
    {
      if (const std::optional<int> sign = OnGrid(circle, point))
      {
        return *sign;
      }
      return CircleWhole(circle, point);
    }
  } // namespace

  auto CircleAdapt(const Triangle &circle, Point point, double permanent) -> int
  {
    if (!Held(circle, point)) [[unlikely]]
    {
      return CircleExactly(circle, point);
    }
    const Corner a = {.Difference = circle.First - point,
                      .Roundoff   = Roundoff(circle.First, point)};
    const Corner b = {.Difference = circle.Second - point,
                      .Roundoff   = Roundoff(circle.Second, point)};
    const Corner c = {.Difference = circle.Third - point,
                      .Roundoff   = Roundoff(circle.Third, point)};

    double determinant = Compensated(a, b, c);
    if (IsZero(a.Roundoff) && IsZero(b.Roundoff) && IsZero(c.Roundoff))
    {
      if (Decided(determinant, CircleBoundCompensated * permanent))
      {
        return Sign(determinant);
      }
      return CircleExactly(circle, point);
    }
    if (Decided(determinant, CircleBoundB * permanent))
    {
      return Sign(determinant);
    }

    const Turn first  = {.Own = &a, .Next = &b, .Before = &c};
    const Turn second = {.Own = &b, .Next = &c, .Before = &a};
    const Turn third  = {.Own = &c, .Next = &a, .Before = &b};

    const double bound  = ((CircleBoundC + CircleBoundCompensated) * permanent)
                        + (ResultBound * std::abs(determinant));
    determinant        += (Correction(first) + Correction(second)) + Correction(third);
    if (Decided(determinant, bound))
    {
      return Sign(determinant);
    }
    return CircleExactly(circle, point);
  }
} // namespace Geoxy::Exact
