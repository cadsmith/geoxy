/***************************************************************************************************
 *  @file      source/crossings/bows.cpp
 *  @brief     Where two polylines' bulges' circles cross, from the bulges' own numbers: their
 *             radical line met with the smaller circle, in double-double with error bounds.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "crossings/bows.hpp"

#include "sure.hpp"

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>

namespace Geoxy::Crossings
{
  namespace
  {
    using Sure::Number;

    /**
     *  A pair of numbers: a point's or a step's coordinates.
     */
    using Duo = std::array<Number, 2>;

    /**
     *  A number times a power of two, exactly, its bound with it.
     *
     *  @param  number The number.
     *  @param  factor The power of two.
     *  @return        The number scaled.
     */
    auto Doubled(const Number &number, double factor) -> Number
    {
      return {.Value = {.High = number.Value.High * factor, .Low = number.Value.Low * factor},
              .Bound = number.Bound * factor};
    }

    /**
     *  The dot product of two pairs.
     *
     *  @param  first  One.
     *  @param  second The other.
     *  @return        The product.
     */
    auto Inner(const Duo &first, const Duo &second) -> Number
    {
      return (first[0] * second[0]) + (first[1] * second[1]);
    }

    /**
     *  The step from one point to another, exactly.
     *
     *  @param  head The point reached.
     *  @param  tail The point left.
     *  @return      head - tail.
     */
    auto Stepped(Point head, Point tail) -> Duo
    {
      return {Sure::Apart(head.X, tail.X), Sure::Apart(head.Y, tail.Y)};
    }

    /**
     *  A bulge's circle's q, 4 b c + 2 (1 - b^2) n, for its chord c and bulge b, n the chord turned
     *  a quarter to the left: the circle through Y = 0 where 4 b |Y|^2 = Y.q.
     *
     *  @param  chord The chord.
     *  @param  bulge The bulge.
     *  @return       q.
     */
    auto Pulled(const Duo &chord, double bulge) -> Duo
    {
      const Number four = Sure::Exact(4.0 * bulge);
      const Number rest = Doubled(Sure::Exact(1.0) - Sure::Product(bulge, bulge), 2.0);
      return {(four * chord[0]) - (rest * chord[1]), (four * chord[1]) + (rest * chord[0])};
    }

    /**
     *  The radius of a bulge's circle, roughly, |c| (1 + b^2) / (4 |b|), to take the smaller circle
     *  first.
     *
     *  @param  belly The bulge.
     *  @return       The radius.
     */
    auto Girthed(const Belly &belly) -> double
    {
      const double bulge = belly.Bulge;
      return Length(belly.End - belly.Start) * ((1.0 + (bulge * bulge)) / (4.0 * std::abs(bulge)));
    }

    /**
     *  Where two bulges' circles cross, roughly: Crossbows's radical line met with the first circle
     *  in plain doubles, in the frame of the first bulge's start, and of the two crossings the one
     *  nearer a point, where they lie clearly apart in distance from it.
     *
     *  @param  bellies The two bulges, the smaller circle's first.
     *  @param  toward  The point, in the frame.
     *  @return         The crossing, in the frame; nothing where the circles do not clearly cross,
     *                  or the two lie nearly as near.
     */
    auto Guessed(const std::array<Belly, 2> &bellies, Point toward) -> std::optional<Point>
    {
      const Point  origin = bellies[0].Start;
      const double bend   = bellies[0].Bulge;
      const double wide   = bellies[1].Bulge;
      const Vector chord  = bellies[0].End - origin;
      const Vector gap    = bellies[1].Start - origin;
      const Vector other  = bellies[1].End - bellies[1].Start;
      const double rest   = 2.0 * (1.0 - (bend * bend));
      const double far    = 2.0 * (1.0 - (wide * wide));
      const Vector q      = (chord * (4.0 * bend)) + (Perpendicular(chord) * rest);
      const Vector p      = (other * (4.0 * wide)) + (Perpendicular(other) * far);
      const double both   = 8.0 * bend * wide;
      const Vector a      = (gap * both) - (q * wide) + (p * bend);
      const double lead   = (0.5 * both * Dot(gap, gap)) + (bend * Dot(gap, p));
      const double span   = Dot(a, a);
      const double lean   = Dot(a, q) - (8.0 * bend * lead);
      const double disc   = (span * Dot(q, q)) - (lean * lean);
      if (!(disc > 0.0))
      {
        return std::nullopt;
      }
      const Vector turn  = Perpendicular(a);
      const double along = Dot(turn, q);
      const double root  = std::sqrt(disc);
      const double below = 8.0 * bend * span;
      const double foot  = 8.0 * bend * lead;
      const Point  plus  = Point{} + (((a * foot) + (turn * (along + root))) / below);
      const Point  minus = Point{} + (((a * foot) + (turn * (along - root))) / below);
      const double one   = Distance(plus, toward);
      const double two   = Distance(minus, toward);
      if (!(std::abs(one - two) > 0x1p-20 * (one + two)))
      {
        return std::nullopt;
      }
      return one < two ? plus : minus;
    }

    /**
     *  A bulge's quadratic at a point and its gradient there, each a sure number.
     */
    struct Leveled
    {
      Number Value = {}; /**< G(Y).              */
      Duo    Slope = {}; /**< G's gradient at Y. */
    };

    /**
     *  A bulge's quadratic G(Y) = 4 b (Y - P0).(Y - P1) - 2 (1 - b^2) (Y - P0).n at a point of a
     *  frame, and its gradient there, 4 b ((Y - P0) + (Y - P1)) - 2 (1 - b^2) n, in sure numbers: G
     *  at Y + D is G(Y) + G'(Y).D + 4 b |D|^2, exactly.
     *
     *  @param  belly  The bulge.
     *  @param  origin The frame's origin.
     *  @param  local  The point, in the frame.
     *  @return        G and its gradient there.
     */
    auto Level(const Belly &belly, Point origin, Point local) -> Leveled
    {
      const Duo    chord = Stepped(belly.End, belly.Start);
      const Duo    start = Stepped(belly.Start, origin);
      const Duo    end   = Stepped(belly.End, origin);
      const Duo    from  = {Sure::Exact(local.X) - start[0], Sure::Exact(local.Y) - start[1]};
      const Duo    onto  = {Sure::Exact(local.X) - end[0], Sure::Exact(local.Y) - end[1]};
      const Duo    turn  = {-chord[1], chord[0]};
      const double bulge = belly.Bulge;
      const Number four  = Sure::Exact(4.0 * bulge);
      const Number rest  = Doubled(Sure::Exact(1.0) - Sure::Product(bulge, bulge), 2.0);
      return {.Value = (four * Inner(from, onto)) - (rest * Inner(from, turn)),
              .Slope = {(four * (from[0] + onto[0])) - (rest * turn[0]),
                        (four * (from[1] + onto[1])) - (rest * turn[1])}};
    }
  } // namespace

  auto Narrower(const std::array<Belly, 2> &bellies) -> bool
  {
    return Girthed(bellies[1]) < Girthed(bellies[0]);
  }

  auto Hopped(const std::array<Belly, 2> &bellies, Point near) -> std::optional<Point>
  {
    const Point                origin = bellies[0].Start;
    const Point                toward = {.X = near.X - origin.X, .Y = near.Y - origin.Y};
    const std::optional<Point> guess  = Guessed(bellies, toward);
    if (!guess)
    {
      return std::nullopt;
    }
    const Leveled one = Level(bellies[0], origin, *guess);
    const Leveled two = Level(bellies[1], origin, *guess);
    const Number  det = (one.Slope[0] * two.Slope[1]) - (one.Slope[1] * two.Slope[0]);
    const Duo step = {Sure::Quotient((two.Value * one.Slope[1]) - (one.Value * two.Slope[1]), det),
                      Sure::Quotient((one.Value * two.Slope[0]) - (two.Value * one.Slope[0]), det)};
    const double least = Sure::Least(det);
    const double rows  = std::max(Sure::Most(two.Slope[1]) + Sure::Most(one.Slope[1]),
                                  Sure::Most(two.Slope[0]) + Sure::Most(one.Slope[0]));
    const double bend  = std::max(std::abs(bellies[0].Bulge), std::abs(bellies[1].Bulge));
    const double reach = std::max(Sure::Most(step[0]), Sure::Most(step[1]));
    /**
     *  The circles' crossings are the fixed points of D = d - J^-1 (4 b1, 4 b2) |D|^2, d the step
     *  worked out, within reach of it, and J's inverse at most rows over the determinant's least in
     *  its largest row: where pull, that times the larger bulge, times 64 reach is no more than 1,
     *  there is one within 2 reach of the guess, within 32 pull reach^2 of the step, and no other
     *  within 1 / (16 pull) of the guess (bows.py's quadratic). A determinant whose bound reaches
     *  nothing passes no such test.
     */
    if (!(64.0 * rows * bend * reach * Sure::Grow <= least))
    {
      return std::nullopt;
    }
    const double pull = rows / least * bend * Sure::Grow;
    /**
     *  The point, gap from the guess, is then nearer that crossing than the other by more than
     *  rounding moves either (bows.py's reach).
     */
    const double size = std::abs(origin.X)
                      + std::abs(origin.Y)
                      + std::abs(guess->X)
                      + std::abs(guess->Y)
                      + std::abs(near.X)
                      + std::abs(near.Y);
    const double gap  = (Distance(*guess, toward) * Sure::Grow) + (0x1p-50 * size);
    if (!((2.0 * gap) + (3.0 * reach) + (0x1p-50 * size) < Sure::Shrink / (16.0 * pull)))
    {
      return std::nullopt;
    }
    const double          quadratic = 32.0 * pull * reach * reach * Sure::Grow;
    std::array<double, 2> out       = {};
    for (std::size_t axis = 0; axis < 2; ++axis)
    {
      const double start = axis == 0 ? origin.X : origin.Y;
      const double guide = axis == 0 ? guess->X : guess->Y;
      Number       value = (Sure::Exact(start) + Sure::Exact(guide)) + step.at(axis);
      value.Bound        = (value.Bound + quadratic) * Sure::Grow;
      const std::optional<double> rounded = Sure::Rounded(value);
      if (!rounded)
      {
        return std::nullopt;
      }
      out.at(axis) = *rounded;
    }
    return Point{.X = out[0], .Y = out[1]};
  }

  auto Crossbows(const std::array<Belly, 2> &bellies) -> std::optional<Crossbow>
  {
    const Belly &small = bellies[0];
    const Belly &large = bellies[1];
    const double b     = small.Bulge;
    const double wide  = large.Bulge;
    const Duo    gap   = Stepped(large.Start, small.Start);
    const Duo    q     = Pulled(Stepped(small.End, small.Start), b);
    const Duo    p     = Pulled(Stepped(large.End, large.Start), wide);
    const Number both  = Doubled(Sure::Product(b, wide), 8.0);
    const Number low   = Sure::Exact(b);
    const Number high  = Sure::Exact(wide);
    const Duo    a     = {(both * gap[0]) - (high * q[0]) + (low * p[0]),
                          (both * gap[1]) - (high * q[1]) + (low * p[1])};
    const Number lead  = (Doubled(both, 0.5) * Inner(gap, gap)) + (low * Inner(gap, p));
    const Duo    turn  = {-a[1], a[0]};
    const Number span  = Inner(a, a);
    const Number along = Inner(turn, q);
    const Number eight = Doubled(low, 8.0);
    const Number lean  = Inner(a, q) - (eight * lead);
    const Number disc  = (span * Inner(q, q)) - (lean * lean);
    if (!((disc.Value.High - std::abs(disc.Value.Low)) * Sure::Shrink > disc.Bound))
    {
      return std::nullopt;
    }
    const Number root  = Sure::Root(disc);
    const Number below = eight * span;
    const Number foot  = eight * lead;
    Crossbow     out;
    for (std::size_t sign = 0; sign < 2; ++sign)
    {
      const Number share = sign == 0 ? along + root : along - root;
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        const Number top         = (foot * a.at(axis)) + (share * turn.at(axis));
        const double start       = axis == 0 ? small.Start.X : small.Start.Y;
        out.At.at(sign).at(axis) = Sure::Exact(start) + Sure::Quotient(top, below);
      }
    }
    return out;
  }
} // namespace Geoxy::Crossings
