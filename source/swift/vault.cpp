/***************************************************************************************************
 *  @file      source/swift/vault.cpp
 *  @brief     A quadratic Bézier curve met with a curve: its implicit equation along each piece of
 *             the other.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/vault.hpp"

#include "basis.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "sure.hpp"
#include "swift/arcade.hpp"
#include "swift/casteljau.hpp"
#include "swift/circuit.hpp"
#include "swift/lane.hpp"
#include "swift/piece.hpp"
#include "swift/roots.hpp"
#include "swift/terms.hpp"
#include "trig.hpp"
#include "twofold.hpp"

#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <variant>

namespace Geoxy::Swift
{
  namespace
  {
    /**
     *  A quadratic Bézier curve's parabola. Each side is the exact difference of two doubles, D
     *  their cross product in double-double; |B'| = 2 |(1 - s) e2 + s e0| is no more than 2
     *  max(|e0|, |e2|) from 0 to 1, and its least is 2 |D| / |e2 - e0|, the distance from nothing
     *  of the line through e2 and e0 times two. A step of its place moves it by no more than the
     *  thorough way's Strict charges a Bézier curve's, four times twice its longest control step's
     *  coordinate times a step of the doubles at 2.
     *
     *  @param  quadratic The curve.
     *  @return           Its parabola; nothing where a number lies outside the range the kernels
     *                    take or D may be as small as Narrowest, a curve too flat to be one.
     */
    auto Vaulting(const Bezier::Quadratic &quadratic) -> std::optional<Vault>
    {
      const Point start   = quadratic.Start;
      const Point control = quadratic.Control;
      const Point end     = quadratic.End;
      if (!Held({start.X, start.Y, control.X, control.Y, end.X, end.Y}))
      {
        return std::nullopt;
      }
      const auto side = [](Point head, Point tail) -> std::array<Number, 2>
      {
        return {Sure::Apart(head.X, tail.X), Sure::Apart(head.Y, tail.Y)};
      };
      Vault vault = {.Hub   = control,
                     .Sides = {side(end, control), side(end, start), side(control, start)}};
      for (std::size_t which = 0; which < vault.Sides.size(); ++which)
      {
        const std::array<Number, 2> &exact = vault.Sides.at(which);
        vault.Ways.at(which)               = {.X = exact[0].Value.High, .Y = exact[1].Value.High};
        vault.Slips.at(which) =
            std::max(std::abs(exact[0].Value.Low), std::abs(exact[1].Value.Low));
        vault.Long.at(which) =
            (Length(vault.Ways.at(which)) + (2.0 * vault.Slips.at(which))) * Grow;
      }
      const std::array<Number, 2> &after  = vault.Sides[0];
      const std::array<Number, 2> &before = vault.Sides[2];
      vault.Area                          = (before[0] * after[1]) - (before[1] * after[0]);
      const double least                  = Sure::Least(vault.Area);
      if (!(least > Narrowest))
      {
        return std::nullopt;
      }
      const double axis =
          Trig::Hypot(Sure::Most(before[0] - after[0]), Sure::Most(before[1] - after[1])) * Grow;
      double pace = 0.0;
      for (const Vector step : {control - start, end - control})
      {
        pace = std::max({pace, 2.0 * std::abs(step.X), 2.0 * std::abs(step.Y)});
      }
      vault.Bend =
          ((2.0 * vault.Long[1] * vault.Long[1]) + (8.0 * vault.Long[0] * vault.Long[2])) * Grow;
      vault.Fast   = 2.0 * std::max(vault.Long[0], vault.Long[2]) * Grow;
      vault.Slow   = 2.0 * least / axis * Shrink;
      vault.Stride = 4.0 * pace * Above(2.0);
      vault.Hull   = Grown(Grown(Box{.Min = start, .Max = start}, control), end);
      return vault;
    }
  } // namespace

  auto Springing(const Vault &vault, const Piece<Ample> &piece, double reach) -> Equation<Ample>
  {
    const std::size_t                   degree = piece.Degree;
    const double                        area   = vault.Area.Value.High;
    const double                        drift  = std::abs(vault.Area.Value.Low) + vault.Area.Bound;
    std::array<Poled<double, Ample>, 3> lines  = {};
    std::array<Poled<double, Ample>, 3> offs   = {};
    double                              most   = 0.0;
    double                              heavy  = 0.0;
    for (std::size_t index = 0; index <= degree; ++index)
    {
      const Vector point  = piece.Moved.at(index);
      const Vector slack  = piece.Offs.at(index);
      const double weight = piece.Weights.at(index);
      for (std::size_t which = 0; which < vault.Ways.size(); ++which)
      {
        const Vector way   = vault.Ways.at(which);
        const double cross = (way.X * point.Y) - (way.Y * point.X);
        const double off =
            (std::abs(way.X) * slack.Y)
            + (std::abs(way.Y) * slack.X)
            + (vault.Slips.at(which) * (std::abs(point.X) + std::abs(point.Y) + slack.X + slack.Y))
            + (4.0 * Unit * (std::abs(way.X * point.Y) + std::abs(way.Y * point.X)));
        const bool weighed        = which == 1;
        lines.at(which).at(index) = weighed ? (area * weight) - cross : cross;
        offs.at(which).at(index)  = (weighed ? off
                                                   + (4.0 * Unit * std::abs(area * weight))
                                                   + (drift * (weight + piece.Heft))
                                                   + (std::abs(area) * piece.Heft)
                                             : off)
                                  * Grow;
      }
      const auto sized = [&lines, &offs, index](std::size_t which) -> double
      {
        return std::abs(lines.at(which).at(index)) + offs.at(which).at(index);
      };
      most  = std::max(most, ((2.0 * sized(1) * vault.Long[1])
                              + (4.0 * sized(0) * vault.Long[2])
                              + (4.0 * sized(2) * vault.Long[0]))
                                 / (weight - piece.Heft));
      heavy = std::max(heavy, weight + piece.Heft);
    }
    Equation<Ample> equation = {.Degree = 2 * degree};
    for (std::size_t power = 0; power <= equation.Degree; ++power)
    {
      double sum  = 0.0;
      double size = 0.0;
      double off  = 0.0;
      for (std::size_t one = power > degree ? power - degree : 0; one <= std::min(power, degree);
           ++one)
      {
        const std::size_t two       = power - one;
        const double      times     = Choose(degree, one) * Choose(degree, two);
        const double      lead      = lines[1].at(one);
        const double      next      = lines[1].at(two);
        const double      first     = lines[0].at(one);
        const double      last      = lines[2].at(two);
        const double      leadOff   = offs[1].at(one);
        const double      nextOff   = offs[1].at(two);
        const double      firstOff  = offs[0].at(one);
        const double      lastOff   = offs[2].at(two);
        sum                        += times * ((lead * next) - (4.0 * first * last));
        size += times * (std::abs(lead * next) + std::abs(4.0 * first * last));
        off  += times
              * ((std::abs(lead) * nextOff)
                 + (std::abs(next) * leadOff)
                 + (leadOff * nextOff)
                 + (4.0
                    * ((std::abs(first) * lastOff)
                       + (std::abs(last) * firstOff)
                       + (firstOff * lastOff))));
      }
      const double below         = Choose(equation.Degree, power);
      const double value         = sum / below;
      equation.Whole.C.at(power) = value;
      equation.First =
          std::max(equation.First, (((off + (16.0 * Unit * size)) / below) * Grow) + Sure::Tiny);
      equation.Largest = std::max(equation.Largest, std::abs(value));
    }
    const Terms<Ample> &values = equation.Whole.C;
    for (std::size_t power = 0; power + 2 <= equation.Degree; ++power)
    {
      equation.Turn =
          std::max(equation.Turn,
                   std::abs(values.at(power + 2) - (2.0 * values.at(power + 1)) + values.at(power))
                       + (8.0 * equation.First)
                       + (8.0 * Unit * equation.Largest));
    }
    const auto order  = static_cast<double>(equation.Degree);
    equation.Turn    *= order * (order - 1.0) * (1.0 + 0x1p-48);
    equation.Heavy    = heavy * heavy;
    equation.Reach    = reach * equation.Heavy * (most + (0.5 * vault.Bend * reach)) * Grow;
    return equation;
  }

  auto Keystone(const Vault &vault, const std::array<Number, 3> &point) -> Number
  {
    const Number gapX  = point[0] - (Sure::Exact(vault.Hub.X) * point[2]);
    const Number gapY  = point[1] - (Sure::Exact(vault.Hub.Y) * point[2]);
    const auto   cross = [&vault, &gapX, &gapY](std::size_t which) -> Number
    {
      return (vault.Sides.at(which)[0] * gapY) - (vault.Sides.at(which)[1] * gapX);
    };
    const Number middle = (vault.Area * point[2]) - cross(1);
    return (middle * middle) - (Sure::Exact(4.0) * cross(0) * cross(2));
  }

  auto Sprung(const Arcade<Vault> &arcade, const Root &root, List<Ribbed> &found) -> bool
  {
    const Vault                 &vault    = *arcade.Carrier;
    const Piece<Ample>          &piece    = *arcade.Part;
    const Equation<Ample>       &equation = arcade.Equation;
    const Honed                  honed    = Honing(arcade, root);
    const Number                &place    = honed.Place;
    const std::array<Number, 2> &crossing = honed.Crossing;
    const double                 near     = honed.Near;
    const Number                 gapX     = crossing[0] - Sure::Exact(vault.Hub.X);
    const Number                 gapY     = crossing[1] - Sure::Exact(vault.Hub.Y);
    const Number                 axisX    = vault.Sides[2][0] - vault.Sides[0][0];
    const Number                 axisY    = vault.Sides[2][1] - vault.Sides[0][1];
    const Number                 bent =
        Sure::Exact(0.5)
        + Sure::Quotient((axisX * gapY) - (axisY * gapX), Sure::Exact(2.0) * vault.Area);
    const double margin =
        ((near / vault.Slow * Grow) + std::abs(bent.Value.Low) + bent.Bound) * Grow;
    const Side side = Within({.Start = 0.0, .End = 1.0}, bent.Value.High, margin);
    if (side != Side::In)
    {
      return side == Side::Out;
    }
    const double sine = root.Slope
                      / (equation.Heavy * 2.0 * Sure::Most(vault.Area) * vault.Fast * piece.Speed)
                      * Shrink;
    if (!Steep(sine, arcade.Tolerance))
    {
      return false;
    }
    const std::optional<std::array<double, 4>> rounded =
        Sure::Rounded(std::array{crossing[0], crossing[1], Lodged(piece, place), bent});
    if (rounded)
    {
      found.push_back({.At    = {.X = (*rounded)[0], .Y = (*rounded)[1]},
                       .Along = (*rounded)[2],
                       .Bent  = (*rounded)[3]});
    }
    return rounded.has_value();
  }

  auto Clearing(const Vault &vault, Point point, double clear) -> bool
  {
    const Number gapX  = Sure::Apart(point.X, vault.Hub.X);
    const Number gapY  = Sure::Apart(point.Y, vault.Hub.Y);
    const auto   cross = [&vault, &gapX, &gapY](std::size_t which) -> Number
    {
      return (vault.Sides.at(which)[0] * gapY) - (vault.Sides.at(which)[1] * gapX);
    };
    const Number first  = cross(0);
    const Number last   = cross(2);
    const Number middle = vault.Area - cross(1);
    const Number level  = (middle * middle) - (Sure::Exact(4.0) * first * last);
    const double most   = ((2.0 * Sure::Most(middle) * vault.Long[1])
                           + (4.0 * Sure::Most(first) * vault.Long[2])
                           + (4.0 * Sure::Most(last) * vault.Long[0]))
                        * Grow;
    return Sure::Least(level) > ((most * clear) + (0.5 * vault.Bend * clear * clear)) * Grow;
  }

  auto Vaulted(const Curves &curves, const Tolerance &tolerance, Intersections &into)
      -> std::optional<bool>
  {
    const std::size_t which = std::holds_alternative<Bezier::Quadratic>(curves[0].get()) ? 0 : 1;
    const auto       *quadratic = std::get_if<Bezier::Quadratic>(&curves.at(which).get());
    const auto       *spline    = std::get_if<Spline>(&curves.at(1 - which).get());
    if (quadratic == nullptr
        || spline == nullptr
        || !(Rational(*spline) || std::cmp_greater(spline->Degree, Highest)))
    {
      return std::nullopt;
    }
    const std::optional<Vault> vault = Vaulting(*quadratic);
    return vault
        && Trussed(*vault, {.Pair = &curves, .Side = 1 - which, .Tolerance = tolerance}, into);
  }
} // namespace Geoxy::Swift
