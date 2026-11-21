/***************************************************************************************************
 *  @file      source/corner/circle.cpp
 *  @brief     The circles of a radius touching two lines or circles, found in doubles, told to be
 *             by the sign of their radicand, and settled to the nearest doubles.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "corner/circle.hpp"

#include "corner/forms.hpp"
#include "corner/surd.hpp"
#include "expansion.hpp"
#include "real.hpp"
#include "sure.hpp"

#include "geoxy/point.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy::Corner
{
  namespace
  {
    /**
     *  A circle's sides swapped with a line's, for a circle met before a line.
     *
     *  @param  ask The ask, the circle's side first.
     *  @return     The ask, the line's side first.
     */
    auto Swapped(const Ask &ask) -> Ask
    {
      return {.Radius = ask.Radius, .First = ask.Second, .Second = ask.First, .Root = ask.Root};
    }

    /**
     *  The circle of an ask as forms of a kind of numbers, the touches in the problem's order.
     *
     *  @tparam NumberT The numbers.
     *  @param  problem The problem.
     *  @param  ask     Which circle.
     *  @return         Its forms.
     */
    template <typename NumberT>
    auto Made(const Problem &problem, const Ask &ask) -> Touching<NumberT>
    {
      const Carrier &first  = problem.Carriers[0];
      const Carrier &second = problem.Carriers[1];
      if (first.Linear && second.Linear)
      {
        return LineLine<NumberT>(first.Line, second.Line, ask);
      }
      if (first.Linear)
      {
        return LineCircle<NumberT>(first.Line, second.Circle, ask);
      }
      if (second.Linear)
      {
        Touching<NumberT> touching = LineCircle<NumberT>(second.Line, first.Circle, Swapped(ask));
        std::swap(touching.First, touching.Second);
        return touching;
      }
      return CircleCircle<NumberT>(first.Circle, second.Circle, ask);
    }

    /**
     *  What must not be below nothing for a circle of a line or circle and a circle to be: the
     *  radicand LineCircle and CircleCircle work out, alone.
     *
     *  @tparam NumberT The numbers.
     *  @param  problem The problem, a circle in it.
     *  @param  ask     Which circle.
     *  @return         The radicand, as a form of one root.
     */
    template <typename NumberT>
    auto Under(const Problem &problem, const Ask &ask) -> Form<NumberT>
    {
      const Carrier &first  = problem.Carriers[0];
      const Carrier &second = problem.Carriers[1];
      const auto     one    = Held<NumberT>(1.0);
      const auto     radius = Held<NumberT>(ask.Radius);
      if (first.Linear || second.Linear)
      {
        const bool      ahead  = first.Linear;
        const Straight &line   = ahead ? first.Line : second.Line;
        const Round    &round  = ahead ? second.Circle : first.Circle;
        const auto      dx     = Twin<NumberT>(line.High.X, line.Low.X);
        const auto      dy     = Twin<NumberT>(line.High.Y, line.Low.Y);
        const auto      wayX   = Apart<NumberT>(round.Center.X, line.Anchor.X);
        const auto      wayY   = Apart<NumberT>(round.Center.Y, line.Anchor.Y);
        const NumberT   square = (dx * dx) + (dy * dy);
        const NumberT   turned = (dx * wayY) - (dy * wayX);
        const NumberT   reach =
            Held<NumberT>(round.Radius) + Flip(radius, ahead ? ask.Second : ask.First);
        const NumberT side = Flip(radius, ahead ? ask.First : ask.Second);
        return {.Whole = (square * ((reach * reach) - (radius * radius))) - (turned * turned),
                .Inner = Held<NumberT>(2.0) * side * turned,
                .Root  = square,
                .Below = one};
      }
      const auto    wayX   = Apart<NumberT>(second.Circle.Center.X, first.Circle.Center.X);
      const auto    wayY   = Apart<NumberT>(second.Circle.Center.Y, first.Circle.Center.Y);
      const NumberT square = (wayX * wayX) + (wayY * wayY);
      const auto    reach1 = Held<NumberT>(first.Circle.Radius) + Flip(radius, ask.First);
      const auto    reach2 = Held<NumberT>(second.Circle.Radius) + Flip(radius, ask.Second);
      const NumberT middle = square + (reach1 * reach1) - (reach2 * reach2);
      return {.Whole = (Held<NumberT>(4.0) * square * reach1 * reach1) - (middle * middle),
              .Below = one};
    }

    /**
     *  The sign of a form's value, worked in sure arithmetic and, where its bound leaves it open,
     *  exactly.
     *
     *  @param  problem The problem.
     *  @param  ask     Which circle's radicand.
     *  @return         1, -1 or 0.
     */
    auto Radicand(const Problem &problem, const Ask &ask) -> int
    {
      const Carrier &first  = problem.Carriers[0];
      const Carrier &second = problem.Carriers[1];
      if (first.Linear && second.Linear)
      {
        return 1;
      }
      const Sure::Number value = Surely(Under<Sure::Number>(problem, ask));
      if (Sure::Least(value) > 0.0)
      {
        return value.Value.High > 0.0 ? 1 : -1;
      }
      const Form<Big> exact = Under<Big>(problem, ask);
      return Signed(exact.Whole, exact.Inner, exact.Root);
    }

    /**
     *  A number scaled by a power of two.
     *
     *  @param  value The number.
     *  @param  power The power.
     *  @return       value 2^power.
     */
    auto Raised(double value, int power) -> double
    {
      return std::ldexp(value, power);
    }

    /**
     *  A point scaled by a power of two.
     *
     *  @param  point The point.
     *  @param  power The power.
     *  @return       It, scaled.
     */
    auto Raised(Point point, int power) -> Point
    {
      return {.X = Raised(point.X, power), .Y = Raised(point.Y, power)};
    }

    /**
     *  A circle's forms as six coordinates, center first.
     *
     *  @tparam NumberT  The numbers.
     *  @param  touching The circle.
     *  @return          The center's, the first touch's and the second's, x before y.
     */
    template <typename NumberT>
    auto Listed(const Touching<NumberT> &touching) -> std::array<Form<NumberT>, 6>
    {
      return {touching.Center.X, touching.Center.Y, touching.First.X,
              touching.First.Y,  touching.Second.X, touching.Second.Y};
    }

    /**
     *  Six coordinates as a circle's center and touches, scaled back.
     *
     *  @param  values The coordinates, as Listed lists them.
     *  @param  power  The power of two to scale them by.
     *  @return        The circle.
     */
    auto Touched(std::span<const double, 6> values, int power) -> Touch
    {
      return {.Center  = Raised(Point{.X = values[0], .Y = values[1]}, power),
              .Touches = {Raised(Point{.X = values[2], .Y = values[3]}, power),
                          Raised(Point{.X = values[4], .Y = values[5]}, power)}};
    }

    /**
     *  The doubles nearest the values of forms: each certified in sure arithmetic where its bound
     *  settles it, the rest found by exact comparisons on the forms made exactly, once, where one
     *  is not.
     *
     *  @tparam CountT How many.
     *  @tparam MakeT  What makes the forms exactly.
     *  @param  fast   The forms in sure numbers.
     *  @param  exact  Makes them exactly, as an array of CountT.
     *  @return        The doubles nearest.
     */
    template <std::size_t CountT, typename MakeT>
    auto Pinned(const std::array<Form<Sure::Number>, CountT> &fast, const MakeT &exact)
        -> std::array<double, CountT>
    {
      std::array<double, CountT>       values = {};
      std::array<Sure::Number, CountT> sure   = {};
      bool                             open   = false;
      for (std::size_t index = 0; index < CountT; ++index)
      {
        sure.at(index)                      = Surely(fast.at(index));
        const std::optional<double> rounded = Sure::Rounded(sure.at(index));
        values.at(index)                    = rounded.value_or(0.0);
        open                                = open || !rounded;
      }
      if (open)
      {
        const std::array<Form<Big>, CountT> sharp = exact();
        for (std::size_t index = 0; index < CountT; ++index)
        {
          const Sure::Number &number = sure.at(index);
          if (!Sure::Rounded(number))
          {
            values.at(index) =
                Pinpointed(sharp.at(index), number.Value.High + number.Value.Low, number.Bound);
          }
        }
      }
      return values;
    }

    /**
     *  The power of two that scales numbers so that the largest is near 1.
     *
     *  @param  numbers The numbers.
     *  @return         The power the largest's magnitude has: it is in [2^(p - 1), 2^p).
     */
    auto Power(std::initializer_list<double> numbers) -> int
    {
      double largest = 0.0;
      for (const double number : numbers)
      {
        largest = std::max(largest, std::abs(number));
      }
      int power = 0;
      std::frexp(largest, &power);
      return power;
    }

    /**
     *  A line scaled by a power of two.
     *
     *  @param  line  The line.
     *  @param  power The power.
     *  @return       It, scaled.
     */
    auto Raised(const Straight &line, int power) -> Straight
    {
      return {.Anchor = Raised(line.Anchor, power),
              .High   = {.X = Raised(line.High.X, power), .Y = Raised(line.High.Y, power)},
              .Low    = {.X = Raised(line.Low.X, power), .Y = Raised(line.Low.Y, power)}};
    }

    /**
     *  The half circle from a point to a line parallel to its own: the point's foot, the midpoint
     *  between, and the radius, as forms.
     *
     *  With d the line's way, n its square and k = d x (F - P), the foot is F - (k / n) perp(d),
     *  the center halfway, and the radius sqrt(k^2 n) / (2 n).
     *
     *  @tparam NumberT The numbers.
     *  @param  line    The line.
     *  @param  from    The point.
     *  @return         The center's, the foot's coordinates and the radius.
     */
    template <typename NumberT>
    auto Spanning(const Straight &line, Point from) -> std::array<Form<NumberT>, 5>
    {
      const auto    dx     = Twin<NumberT>(line.High.X, line.Low.X);
      const auto    dy     = Twin<NumberT>(line.High.Y, line.Low.Y);
      const auto    gapX   = Apart<NumberT>(from.X, line.Anchor.X);
      const auto    gapY   = Apart<NumberT>(from.Y, line.Anchor.Y);
      const NumberT square = (dx * dx) + (dy * dy);
      const NumberT turned = (dx * gapY) - (dy * gapX);
      const auto    twice  = Held<NumberT>(2.0) * square;
      const auto    startX = Held<NumberT>(from.X);
      const auto    startY = Held<NumberT>(from.Y);
      return {Form<NumberT>{.Whole = (twice * startX) + (turned * dy), .Below = twice},
              Form<NumberT>{.Whole = (twice * startY) - (turned * dx), .Below = twice},
              Form<NumberT>{.Whole = (square * startX) + (turned * dy), .Below = square},
              Form<NumberT>{.Whole = (square * startY) - (turned * dx), .Below = square},
              Form<NumberT>{
                  .Inner = Held<NumberT>(1.0), .Root = turned * turned * square, .Below = twice}};
    }

    /**
     *  Where CHAMFER cuts two lines, as forms: the first cut's coordinates, then the second's.
     *
     *  @tparam NumberT The numbers.
     *  @param  first   One line, scaled.
     *  @param  second  The other, scaled.
     *  @param  cut     The distances, scaled, the sides kept and the angle's sine and cosine.
     *  @return         The four forms, and the sine rule's denominator, which must be above
     *                  nothing.
     */
    template <typename NumberT>
    auto Cutting(const Straight &first, const Straight &second, const Cutback &cut)
        -> std::array<Form<NumberT>, 5>
    {
      const auto    dx1     = Twin<NumberT>(first.High.X, first.Low.X);
      const auto    dy1     = Twin<NumberT>(first.High.Y, first.Low.Y);
      const auto    dx2     = Twin<NumberT>(second.High.X, second.Low.X);
      const auto    dy2     = Twin<NumberT>(second.High.Y, second.Low.Y);
      const auto    gapX    = Apart<NumberT>(second.Anchor.X, first.Anchor.X);
      const auto    gapY    = Apart<NumberT>(second.Anchor.Y, first.Anchor.Y);
      const NumberT lead    = (dx1 * dx1) + (dy1 * dy1);
      const NumberT trail   = (dx2 * dx2) + (dy2 * dy2);
      const NumberT across  = (dx1 * dy2) - (dy1 * dx2);
      const NumberT turned  = (dx2 * gapY) - (dy2 * gapX);
      const NumberT cornerX = (Held<NumberT>(first.Anchor.X) * across) - (turned * dx1);
      const NumberT cornerY = (Held<NumberT>(first.Anchor.Y) * across) - (turned * dy1);
      const auto    near    = Held<NumberT>(cut.Distances[0] * cut.Keeps[0]);
      const auto    one     = Held<NumberT>(1.0);
      const NumberT below   = across * lead;
      std::array<Form<NumberT>, 5> forms = {
          Form<NumberT>{
              .Whole = cornerX * lead, .Inner = near * dx1 * across, .Root = lead, .Below = below},
          Form<NumberT>{
              .Whole = cornerY * lead, .Inner = near * dy1 * across, .Root = lead, .Below = below},
          Form<NumberT>{.Whole = cornerX * trail,
                        .Inner = Held<NumberT>(cut.Distances[1] * cut.Keeps[1]) * dx2 * across,
                        .Root  = trail,
                        .Below = across * trail},
          Form<NumberT>{.Whole = cornerY * trail,
                        .Inner = Held<NumberT>(cut.Distances[1] * cut.Keeps[1]) * dy2 * across,
                        .Root  = trail,
                        .Below = across * trail},
          Form<NumberT>{.Whole = one, .Below = one}};
      if (cut.Slant)
      {
        const auto    sine   = Held<NumberT>((*cut.Slant)[0]);
        const auto    cosine = Held<NumberT>((*cut.Slant)[1]);
        const NumberT dot    = (dx1 * dx2) + (dy1 * dy2);
        const NumberT apart  = Flip(dot, -cut.Keeps[0] * cut.Keeps[1]);
        const NumberT size   = Flip(across, cut.Across);
        const NumberT width  = (sine * apart) + (cosine * size);
        const auto    reach  = Held<NumberT>(cut.Distances[0] * cut.Keeps[1]) * sine * across;
        forms[2]             = {
            .Whole = cornerX * width, .Inner = reach * dx2, .Root = lead, .Below = across * width};
        forms[3] = {
            .Whole = cornerY * width, .Inner = reach * dy2, .Root = lead, .Below = across * width};
        forms[4] = {.Whole = width, .Below = one};
      }
      return forms;
    }
  } // namespace

  auto Posed(const std::array<Carrier, 2> &carriers, double radius) -> Problem
  {
    int power = Power({radius});
    for (const Carrier &carrier : carriers)
    {
      const Straight &line   = carrier.Line;
      const Round    &circle = carrier.Circle;
      power = std::max(power, Power({line.Anchor.X, line.Anchor.Y, line.High.X, line.High.Y,
                                     circle.Center.X, circle.Center.Y, circle.Radius}));
    }
    Problem problem = {.Carriers = carriers, .Radius = Raised(radius, -power), .Power = power};
    for (Carrier &carrier : problem.Carriers)
    {
      carrier.Line          = Raised(carrier.Line, -power);
      carrier.Circle.Center = Raised(carrier.Circle.Center, -power);
      carrier.Circle.Radius = Raised(carrier.Circle.Radius, -power);
    }
    return problem;
  }

  auto Asks(const Problem &problem) -> std::vector<Ask>
  {
    const bool       lines = problem.Carriers[0].Linear && problem.Carriers[1].Linear;
    const bool       none  = !(problem.Radius > 0.0);
    std::vector<Ask> asks;
    const std::array signs = {1, -1};
    for (const int first : signs)
    {
      for (const int second : signs)
      {
        for (const int root : signs)
        {
          const bool repeated = (none && (first < 0 || second < 0)) || (lines && root < 0);
          if (!repeated)
          {
            asks.push_back(
                {.Radius = problem.Radius, .First = first, .Second = second, .Root = root});
          }
        }
      }
    }
    return asks;
  }

  auto Approximated(const Problem &problem, const Ask &ask) -> std::optional<Touch>
  {
    const std::array sides = {ask.First, ask.Second};
    for (std::size_t index = 0; index < 2; ++index)
    {
      const Carrier &carrier = problem.Carriers.at(index);
      if (!carrier.Linear && sides.at(index) < 0 && Exactly(carrier.Circle.Radius, problem.Radius))
      {
        return std::nullopt;
      }
    }
    const Carrier &first  = problem.Carriers[0];
    const Carrier &second = problem.Carriers[1];
    if (!first.Linear && !second.Linear && first.Circle.Center == second.Circle.Center)
    {
      return std::nullopt;
    }
    const int sign = Radicand(problem, ask);
    if (sign < 0 || (sign == 0 && ask.Root < 0))
    {
      return std::nullopt;
    }
    const std::array<Form<double>, 6> forms  = Listed(Made<double>(problem, ask));
    std::array<double, 6>             values = {};
    std::ranges::transform(forms, values.begin(),
                           [](const Form<double> &form) -> double { return Approximate(form); });
    return Touched(values, problem.Power);
  }

  auto Settled(const Problem &problem, const Ask &ask) -> Touch
  {
    const std::array<double, 6> values = Pinned(Listed(Made<Sure::Number>(problem, ask)),
                                                [&problem, &ask] -> std::array<Form<Big>, 6>
                                                { return Listed(Made<Big>(problem, ask)); });
    return Touched(values, problem.Power);
  }

  auto Bridged(const Straight &line, Point from) -> Span
  {
    const int power =
        Power({line.Anchor.X, line.Anchor.Y, line.High.X, line.High.Y, from.X, from.Y});
    const Straight              scaled = Raised(line, -power);
    const Point                 start  = Raised(from, -power);
    const std::array<double, 5> values =
        Pinned(Spanning<Sure::Number>(scaled, start), [&scaled, start] -> std::array<Form<Big>, 5>
               { return Spanning<Big>(scaled, start); });
    return {.Touch  = {.Center  = Raised(Point{.X = values[0], .Y = values[1]}, power),
                       .Touches = {from, Raised(Point{.X = values[2], .Y = values[3]}, power)}},
            .Radius = Raised(values[4], power)};
  }

  auto Bevelled(const Straight &first, const Straight &second, const Cutback &cut)
      -> std::optional<std::array<Point, 2>>
  {
    const int power =
        std::max({Power({first.Anchor.X, first.Anchor.Y, first.High.X, first.High.Y,
                         cut.Distances[0], cut.Distances[1]}),
                  Power({second.Anchor.X, second.Anchor.Y, second.High.X, second.High.Y})});
    const Straight one    = Raised(first, -power);
    const Straight other  = Raised(second, -power);
    Cutback        scaled = cut;
    scaled.Distances      = {Raised(cut.Distances[0], -power), Raised(cut.Distances[1], -power)};
    const Big dx1         = Twin<Big>(one.High.X, one.Low.X);
    const Big dy1         = Twin<Big>(one.High.Y, one.Low.Y);
    const Big dx2         = Twin<Big>(other.High.X, other.Low.X);
    const Big dy2         = Twin<Big>(other.High.Y, other.Low.Y);
    scaled.Across         = Expansion::Sign(((dx1 * dy2) - (dy1 * dx2)).Value);
    const std::array<Form<Sure::Number>, 5> fast = Cutting<Sure::Number>(one, other, scaled);
    if (cut.Slant)
    {
      const Sure::Number width = Surely(fast[4]);
      int                sign  = width.Value.High > 0.0 ? 1 : -1;
      if (!(Sure::Least(width) > 0.0))
      {
        sign = Expansion::Sign(Cutting<Big>(one, other, scaled)[4].Whole.Value);
      }
      if (sign <= 0)
      {
        return std::nullopt;
      }
    }
    const std::array<Form<Sure::Number>, 4> points = {fast[0], fast[1], fast[2], fast[3]};
    const std::array<double, 4> values = Pinned(points,
                                                [&one, &other, &scaled] -> std::array<Form<Big>, 4>
                                                {
                                                  const std::array<Form<Big>, 5> exact =
                                                      Cutting<Big>(one, other, scaled);
                                                  return {exact[0], exact[1], exact[2], exact[3]};
                                                });
    return std::array{Raised(Point{.X = values[0], .Y = values[1]}, power),
                      Raised(Point{.X = values[2], .Y = values[3]}, power)};
  }

  auto Parallel(const Straight &first, const Straight &second) -> bool
  {
    const Big dx1 = Twin<Big>(first.High.X, first.Low.X);
    const Big dy1 = Twin<Big>(first.High.Y, first.Low.Y);
    const Big dx2 = Twin<Big>(second.High.X, second.Low.X);
    const Big dy2 = Twin<Big>(second.High.Y, second.Low.Y);
    return Expansion::Sign(((dx1 * dy2) - (dy1 * dx2)).Value) == 0;
  }
} // namespace Geoxy::Corner
