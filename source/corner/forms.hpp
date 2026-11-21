/***************************************************************************************************
 *  @file      source/corner/forms.hpp
 *  @brief     Where a circle of a radius touches two lines, a line and a circle, or two circles, as
 *             forms of nested square roots of the curves' own numbers; worked in doubles, in sure
 *             double-doubles, or exactly, from one set of formulas.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CORNER_FORMS_HPP
#define GEOXY_CORNER_FORMS_HPP

#include "corner/surd.hpp"
#include "expansion.hpp"
#include "real.hpp"
#include "sure.hpp"
#include "twofold.hpp"

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace Geoxy::Corner
{
  /**
   *  A double as a number of a kind.
   *
   *  @tparam NumberT The kind.
   *  @param  value   The double.
   *  @return         It, exactly.
   */
  template <typename NumberT>
  auto Held(double value) -> NumberT;

  /**
   *  A double as a double.
   *
   *  @param  value The double.
   *  @return       It.
   */
  template <>
  inline auto Held<double>(double value) -> double
  {
    return value;
  }

  /**
   *  A double as a sure number.
   *
   *  @param  value The double.
   *  @return       It, exactly.
   */
  template <>
  inline auto Held<Sure::Number>(double value) -> Sure::Number
  {
    return Sure::Exact(value);
  }

  /**
   *  A double as an exact number.
   *
   *  @param  value The double.
   *  @return       It.
   */
  template <>
  inline auto Held<Big>(double value) -> Big
  {
    return {.Value = Expansion::Of(value)};
  }

  /**
   *  A number held in two doubles as a number of a kind.
   *
   *  @tparam NumberT The kind.
   *  @param  high    Its leading part.
   *  @param  low     What the leading part leaves, within half a step of it.
   *  @return         Their sum, exactly where the kind holds it.
   */
  template <typename NumberT>
  auto Twin(double high, double low) -> NumberT;

  /**
   *  Two doubles' sum as a double.
   *
   *  @param  high The leading part.
   *  @param  low  The rest.
   *  @return      Their sum, rounded.
   */
  template <>
  inline auto Twin<double>(double high, double low) -> double
  {
    return high + low;
  }

  /**
   *  Two doubles' sum as a sure number.
   *
   *  @param  high The leading part.
   *  @param  low  The rest, within half a step of it.
   *  @return      Their sum, exactly.
   */
  template <>
  inline auto Twin<Sure::Number>(double high, double low) -> Sure::Number
  {
    return {.Value = {.High = high, .Low = low}, .Bound = 0.0};
  }

  /**
   *  Two doubles' sum as an exact number.
   *
   *  @param  high The leading part.
   *  @param  low  The rest, within half a step of it.
   *  @return      Their sum.
   */
  template <>
  inline auto Twin<Big>(double high, double low) -> Big
  {
    return {.Value = Expansion::Sum(Expansion::Of(low), Expansion::Of(high))};
  }

  /**
   *  The difference of two doubles as a number of a kind.
   *
   *  @tparam NumberT The kind.
   *  @param  head    The double.
   *  @param  tail    The double taken from it.
   *  @return         The difference, exactly where the kind holds it.
   */
  template <typename NumberT>
  auto Apart(double head, double tail) -> NumberT
  {
    const Twofold::Two difference = Twofold::TwoSum(head, -tail);
    return Twin<NumberT>(difference.High, difference.Low);
  }

  /**
   *  A number, or its negative.
   *
   *  @tparam NumberT The kind.
   *  @param  value   The number.
   *  @param  side    1 for the number, -1 for its negative.
   *  @return         It.
   */
  template <typename NumberT>
  auto Flip(const NumberT &value, int side) -> NumberT
  {
    return side > 0 ? value : -value;
  }

  /**
   *  A line a fillet touches, taken whole: a point on it and its way, held exactly in two doubles,
   *  as a line's end less its start leaves it.
   */
  struct Straight
  {
    Point  Anchor = {}; /**< A point it passes through. */
    Vector High   = {}; /**< Its way, rounded.          */
    Vector Low    = {}; /**< What rounding left out.    */
  };

  /**
   *  A circle a fillet touches.
   */
  struct Round
  {
    Point  Center = {};  /**< Its center.                */
    double Radius = 0.0; /**< Its radius, above nothing. */
  };

  /**
   *  Which circle of a radius is asked for: the side of each curve it lies on and which root.
   */
  struct Ask
  {
    double Radius = 0.0; /**< The radius, not below nothing. */
    int    First  = 1;   /**< The side of the first: of a line, 1 left of its way and -1 right; of
                              a circle, 1 outside it and -1 within or around it.                */
    int Second = 1;      /**< The side of the second, likewise.               */
    int Root   = 1;      /**< Which root, 1 or -1, where a circle is touched. */
  };

  /**
   *  A point as two forms.
   *
   *  @tparam NumberT The numbers.
   */
  template <typename NumberT>
  struct Spot
  {
    Form<NumberT> X = {}; /**< Its x. */
    Form<NumberT> Y = {}; /**< Its y. */
  };

  /**
   *  A circle touching two curves, as forms: its center, where it touches each, and what must not
   *  be below nothing for it to be at all.
   *
   *  @tparam NumberT The numbers.
   */
  template <typename NumberT>
  struct Touching
  {
    Spot<NumberT> Center   = {}; /**< Its center.                          */
    Spot<NumberT> First    = {}; /**< Where it touches the first curve.    */
    Spot<NumberT> Second   = {}; /**< Where it touches the second.         */
    Form<NumberT> Radicand = {}; /**< Not below nothing where it is found. */
  };

  /**
   *  A form from the roots it shares with others and its own three factors and denominator.
   *
   *  @tparam NumberT The numbers.
   *  @param  shared  The roots' arguments: Root, Base and Bend.
   *  @param  terms   Whole, Inner and Outer.
   *  @param  below   The denominator.
   *  @return         The form.
   */
  template <typename NumberT>
  auto Formed(const Form<NumberT> &shared, const std::array<NumberT, 3> &terms,
              const NumberT &below) -> Form<NumberT>
  {
    Form<NumberT> form = shared;
    form.Whole         = std::get<0>(terms);
    form.Inner         = std::get<1>(terms);
    form.Outer         = std::get<2>(terms);
    form.Below         = below;
    return form;
  }

  /**
   *  The circle of a radius touching two lines, on the sides asked, and their corner where the
   *  radius is nothing.
   *
   *  With d1, d2 their ways, n1, n2 their squares, c = d1 x d2, and X their corner, X c = P1 c -
   *  (d2 x (P2 - P1)) d1; the center is X + (r / c) (s1 sqrt(n1) d2 - s2 sqrt(n2) d1), the touches
   *  X + (r / c) d1 (s1 (d1 . d2) / sqrt(n1) - s2 sqrt(n2)) and X + (r / c) d2 (s1 sqrt(n1) - s2
   *  (d1 . d2) / sqrt(n2)), each over one denominator.
   *
   *  @tparam NumberT The numbers.
   *  @param  first   One line.
   *  @param  second  The other, not parallel to it.
   *  @param  ask     The radius and the side of each.
   *  @return         The circle.
   */
  template <typename NumberT>
  auto LineLine(const Straight &first, const Straight &second, const Ask &ask) -> Touching<NumberT>
  {
    const auto          dx1         = Twin<NumberT>(first.High.X, first.Low.X);
    const auto          dy1         = Twin<NumberT>(first.High.Y, first.Low.Y);
    const auto          dx2         = Twin<NumberT>(second.High.X, second.Low.X);
    const auto          dy2         = Twin<NumberT>(second.High.Y, second.Low.Y);
    const auto          gapX        = Apart<NumberT>(second.Anchor.X, first.Anchor.X);
    const auto          gapY        = Apart<NumberT>(second.Anchor.Y, first.Anchor.Y);
    const NumberT       lead        = (dx1 * dx1) + (dy1 * dy1);
    const NumberT       trail       = (dx2 * dx2) + (dy2 * dy2);
    const NumberT       across      = (dx1 * dy2) - (dy1 * dx2);
    const NumberT       dot         = (dx1 * dx2) + (dy1 * dy2);
    const NumberT       turned      = (dx2 * gapY) - (dy2 * gapX);
    const NumberT       cornerX     = (Held<NumberT>(first.Anchor.X) * across) - (turned * dx1);
    const NumberT       cornerY     = (Held<NumberT>(first.Anchor.Y) * across) - (turned * dy1);
    const auto          radius      = Held<NumberT>(ask.Radius);
    const NumberT       near        = Flip(radius, ask.First);
    const NumberT       far         = Flip(radius, -ask.Second);
    const auto          one         = Held<NumberT>(1.0);
    const Form<NumberT> shared      = {.Root = lead, .Base = trail};
    const NumberT       firstBelow  = across * lead;
    const NumberT       secondBelow = across * trail;
    Touching<NumberT>   touching;
    touching.Center = {.X = Formed(shared, {cornerX, near * dx2, far * dx1}, across),
                       .Y = Formed(shared, {cornerY, near * dy2, far * dy1}, across)};
    touching.First  = {
        .X = Formed(shared, {cornerX * lead, near * dx1 * dot, far * dx1 * lead}, firstBelow),
        .Y = Formed(shared, {cornerY * lead, near * dy1 * dot, far * dy1 * lead}, firstBelow)};
    touching.Second = {
        .X = Formed(shared, {cornerX * trail, near * dx2 * trail, far * dx2 * dot}, secondBelow),
        .Y = Formed(shared, {cornerY * trail, near * dy2 * trail, far * dy2 * dot}, secondBelow)};
    touching.Radicand = {.Whole = one, .Below = one};
    return touching;
  }

  /**
   *  The circle of a radius touching a line and a circle, on the sides asked, by the root asked;
   *  where the radius is nothing, a point where they meet.
   *
   *  With d the line's way and n its square, w the circle's center less the line's point, e = d x
   *  w, f = d . w, s r the radius on the line's side and rho the circle's radius grown by the
   *  radius outside it and shrunk by it within: the center is P + ((f + root sqrt(Q)) / n) d + s r
   *  perp(d) / sqrt(n), Q = n (rho^2 - r^2) - e^2 + 2 s r e sqrt(n); the line is touched at its
   *  foot, and the circle at C + (O - C) R / rho.
   *
   *  @tparam NumberT The numbers.
   *  @param  line    The line.
   *  @param  circle  The circle.
   *  @param  ask     The radius, the side of the line, the side of the circle and the root.
   *  @return         The circle, the line's touch first; Q its radicand.
   */
  template <typename NumberT>
  auto LineCircle(const Straight &line, const Round &circle, const Ask &ask) -> Touching<NumberT>
  {
    const auto          dx     = Twin<NumberT>(line.High.X, line.Low.X);
    const auto          dy     = Twin<NumberT>(line.High.Y, line.Low.Y);
    const auto          wayX   = Apart<NumberT>(circle.Center.X, line.Anchor.X);
    const auto          wayY   = Apart<NumberT>(circle.Center.Y, line.Anchor.Y);
    const NumberT       square = (dx * dx) + (dy * dy);
    const NumberT       turned = (dx * wayY) - (dy * wayX);
    const NumberT       ahead  = (dx * wayX) + (dy * wayY);
    const auto          radius = Held<NumberT>(ask.Radius);
    const auto          big    = Held<NumberT>(circle.Radius);
    const NumberT       reach  = big + Flip(radius, ask.Second);
    const NumberT       side   = Flip(radius, ask.First);
    const NumberT       base = (square * ((reach * reach) - (radius * radius))) - (turned * turned);
    const auto          bend = Held<NumberT>(2.0) * side * turned;
    const NumberT       wholeX  = (Held<NumberT>(line.Anchor.X) * square) + (dx * ahead);
    const NumberT       wholeY  = (Held<NumberT>(line.Anchor.Y) * square) + (dy * ahead);
    const NumberT       outerX  = Flip(dx, ask.Root);
    const NumberT       outerY  = Flip(dy, ask.Root);
    const NumberT       innerX  = -(side * dy);
    const NumberT       innerY  = side * dx;
    const NumberT       zero    = {};
    const auto          one     = Held<NumberT>(1.0);
    const NumberT       below   = square * reach;
    const NumberT       awayX   = (dx * ahead) - (wayX * square);
    const NumberT       awayY   = (dy * ahead) - (wayY * square);
    const NumberT       centerX = (Held<NumberT>(circle.Center.X) * below) + (big * awayX);
    const NumberT       centerY = (Held<NumberT>(circle.Center.Y) * below) + (big * awayY);
    const Form<NumberT> shared  = {.Root = square, .Base = base, .Bend = bend};
    Touching<NumberT>   touching;
    touching.First    = {.X = Formed(shared, {wholeX, zero, outerX}, square),
                         .Y = Formed(shared, {wholeY, zero, outerY}, square)};
    touching.Center   = {.X = Formed(shared, {wholeX, innerX, outerX}, square),
                         .Y = Formed(shared, {wholeY, innerY, outerY}, square)};
    touching.Second   = {.X = Formed(shared, {centerX, big * innerX, big * outerX}, below),
                         .Y = Formed(shared, {centerY, big * innerY, big * outerY}, below)};
    touching.Radicand = {.Whole = base, .Inner = bend, .Root = square, .Below = one};
    return touching;
  }

  /**
   *  The circle of a radius touching two circles, on the sides asked, by the root asked; where the
   *  radius is nothing, a point where they meet.
   *
   *  With w the second center less the first, n its square, and rho1, rho2 each radius grown by the
   *  radius outside it and shrunk by it within, the center is C1 + (M w + root sqrt(Q) perp(w)) /
   *  (2 n), M = n + rho1^2 - rho2^2 and Q = 4 n rho1^2 - M^2; each circle is touched at C + (O - C)
   *  R / rho.
   *
   *  @tparam NumberT The numbers.
   *  @param  first   One circle.
   *  @param  second  The other.
   *  @param  ask     The radius, the side of each and the root.
   *  @return         The circle; Q its radicand.
   */
  template <typename NumberT>
  auto CircleCircle(const Round &first, const Round &second, const Ask &ask) -> Touching<NumberT>
  {
    const auto          wayX   = Apart<NumberT>(second.Center.X, first.Center.X);
    const auto          wayY   = Apart<NumberT>(second.Center.Y, first.Center.Y);
    const NumberT       square = (wayX * wayX) + (wayY * wayY);
    const auto          radius = Held<NumberT>(ask.Radius);
    const auto          near   = Held<NumberT>(first.Radius);
    const auto          far    = Held<NumberT>(second.Radius);
    const NumberT       reach1 = near + Flip(radius, ask.First);
    const NumberT       reach2 = far + Flip(radius, ask.Second);
    const NumberT       middle = square + (reach1 * reach1) - (reach2 * reach2);
    const auto          twice  = Held<NumberT>(2.0) * square;
    const NumberT       root = (Held<NumberT>(4.0) * square * reach1 * reach1) - (middle * middle);
    const NumberT       innerX      = Flip(wayY, -ask.Root);
    const NumberT       innerY      = Flip(wayX, ask.Root);
    const NumberT       zero        = {};
    const auto          one         = Held<NumberT>(1.0);
    const NumberT       back        = middle - twice;
    const NumberT       firstBelow  = twice * reach1;
    const NumberT       secondBelow = twice * reach2;
    const auto          startX      = Held<NumberT>(first.Center.X);
    const auto          startY      = Held<NumberT>(first.Center.Y);
    const auto          endX        = Held<NumberT>(second.Center.X);
    const auto          endY        = Held<NumberT>(second.Center.Y);
    const Form<NumberT> shared      = {.Root = root};
    Touching<NumberT>   touching;
    touching.Center = {
        .X = Formed(shared, {(twice * startX) + (middle * wayX), innerX, zero}, twice),
        .Y = Formed(shared, {(twice * startY) + (middle * wayY), innerY, zero}, twice)};
    touching.First = {
        .X = Formed(shared, {(firstBelow * startX) + (near * middle * wayX), near * innerX, zero},
                    firstBelow),
        .Y = Formed(shared, {(firstBelow * startY) + (near * middle * wayY), near * innerY, zero},
                    firstBelow)};
    touching.Second = {
        .X = Formed(shared, {(secondBelow * endX) + (far * back * wayX), far * innerX, zero},
                    secondBelow),
        .Y = Formed(shared, {(secondBelow * endY) + (far * back * wayY), far * innerY, zero},
                    secondBelow)};
    touching.Radicand = {.Whole = root, .Below = one};
    return touching;
  }

  /**
   *  A form's value in doubles, roots of what rounds below nothing taken as nothing: near enough to
   *  choose among the answers by.
   *
   *  @param  form The form.
   *  @return      Its value, rounded.
   */
  inline auto Approximate(const Form<double> &form) -> double
  {
    const double root  = std::sqrt(std::max(form.Root, 0.0));
    const double outer = std::sqrt(std::max(form.Base + (form.Bend * root), 0.0));
    return (form.Whole + (form.Inner * root) + (form.Outer * outer)) / form.Below;
  }

  /**
   *  Whether a sure number is nothing, exactly.
   *
   *  @param  number The number.
   *  @return        Whether its value and its bound are both nothing.
   */
  inline auto Naught(const Sure::Number &number) -> bool
  {
    return Zero(number.Value.High) && Zero(number.Bound);
  }

  /**
   *  A form's value in sure double-doubles, with how far it may be from the exact one: a term whose
   *  factor is nothing exactly is left out, so its root is not asked for.
   *
   *  @param  form The form.
   *  @return      Its value.
   */
  inline auto Surely(const Form<Sure::Number> &form) -> Sure::Number
  {
    Sure::Number top = form.Whole;
    if (!Naught(form.Inner))
    {
      top = top + (form.Inner * Sure::Root(form.Root));
    }
    if (!Naught(form.Outer))
    {
      top = top + (form.Outer * Sure::Root(form.Base + (form.Bend * Sure::Root(form.Root))));
    }
    return Sure::Quotient(top, form.Below);
  }
} // namespace Geoxy::Corner

#endif /* GEOXY_CORNER_FORMS_HPP */
