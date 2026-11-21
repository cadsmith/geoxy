/***************************************************************************************************
 *  @file      source/hodograph.cpp
 *  @brief     The length of a polynomial curve from its hodograph: in closed form where the
 *             hodograph is straight, by Gauss-Legendre rules whose error is bounded from its roots
 *             where it is a parabola.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "hodograph.hpp"

#include "legendre.hpp"
#include "quadrature.hpp"
#include "real.hpp"
#include "trig.hpp"
#include "twofold.hpp"

#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>

namespace Geoxy
{
  namespace
  {
    /**
     *  The thinnest Bernstein ellipse a piece is integrated over: thinner, it is halved.
     */
    constexpr double RhoMin = 3.0;

    /**
     *  How many pieces a stretch is held in at once: its three first, and one more per halving.
     */
    constexpr std::size_t Capacity = Halved + 4;

    /**
     *  A complex number's product with another, each held as a vector: its x the real part, its y
     *  the imaginary.
     *
     *  @param  lhs One number.
     *  @param  rhs The other.
     *  @return     Their product.
     */
    auto Times(Vector lhs, Vector rhs) -> Vector
    {
      return {.X = (lhs.X * rhs.X) - (lhs.Y * rhs.Y), .Y = (lhs.X * rhs.Y) + (lhs.Y * rhs.X)};
    }

    /**
     *  A complex number's quotient by another, both first scaled by the power of two that brings
     *  the divisor near 1, so that its length squared neither overflows nor underflows.
     *
     *  @param  lhs The dividend.
     *  @param  rhs The divisor, not 0.
     *  @return     The quotient.
     */
    auto Over(Vector lhs, Vector rhs) -> Vector
    {
      const int    power = std::ilogb(Largest(rhs));
      const Vector top   = {.X = std::ldexp(lhs.X, -power), .Y = std::ldexp(lhs.Y, -power)};
      const Vector under = {.X = std::ldexp(rhs.X, -power), .Y = std::ldexp(rhs.Y, -power)};
      return Times(top, {.X = under.X, .Y = -under.Y}) / Dot(under, under);
    }

    /**
     *  A complex number's principal square root.
     *
     *  @param  value The number.
     *  @return       Its root, the real part not negative.
     */
    auto Root(Vector value) -> Vector
    {
      const double size = Length(value);
      return {.X = std::sqrt(std::max(0.0, (size + value.X) / 2.0)),
              .Y = std::copysign(std::sqrt(std::max(0.0, (size - value.X) / 2.0)), value.Y)};
    }

    /**
     *  A complex quadratic's power form, s t^2 + l t + c.
     */
    struct Power
    {
      Vector Square   = {}; /**< The coefficient of t^2. */
      Vector Linear   = {}; /**< The coefficient of t.   */
      Vector Constant = {}; /**< The constant.           */
    };

    /**
     *  The roots of a complex quadratic, by the stable formula: the root of the larger magnitude
     *  from -(l + sqrt(l^2 - 4 s c)) / 2 over s, the sign of the root of the discriminant that of
     *  l, so that nothing cancels, and the other from c over the same; a root at 0 twice where l
     *  and c are both 0.
     *
     *  @param  power The quadratic, its square's coefficient not 0.
     *  @return       Its two roots.
     */
    auto Solutions(const Power &power) -> std::array<Vector, 2>
    {
      Vector root =
          Root(Times(power.Linear, power.Linear) - (Times(power.Square, power.Constant) * 4.0));
      if (Dot(power.Linear, root) < 0.0)
      {
        root = -root;
      }
      const Vector half  = (power.Linear + root) * -0.5;
      const Vector first = Over(half, power.Square);
      return {first, Zero(Largest(half)) ? first : Over(power.Constant, half)};
    }

    /**
     *  A hodograph of degree 2 at a parameter, by de Casteljau's algorithm.
     *
     *  @param  points Its control points.
     *  @param  along  The parameter, past 0 and 1 too.
     *  @return        The velocity there.
     */
    auto Velocity(const std::array<Vector, 3> &points, double along) -> Vector
    {
      const double rest  = 1.0 - along;
      const Vector start = (points[0] * rest) + (points[1] * along);
      const Vector end   = (points[1] * rest) + (points[2] * along);
      return (start * rest) + (end * along);
    }

    /**
     *  A lower bound on the length a hodograph of degree 2 runs over a stretch: the polygon through
     *  its points at the stretch's quarters, each side the hodograph's integral over its quarter,
     *  which the rule of 2 points gives exactly for a polynomial of degree 3, less 2^-40 of it for
     *  the rounding of the sides.
     *
     *  @param  points  Its control points.
     *  @param  stretch The stretch.
     *  @return         The lower bound.
     */
    auto Lower(const std::array<Vector, 3> &points, Range stretch) -> double
    {
      const double node    = Legendre::Nodes[0];
      const double quarter = (stretch.End - stretch.Start) / 4.0;
      const double half    = quarter / 2.0;
      double       total   = 0.0;
      for (const double offset : {0.5, 1.5, 2.5, 3.5})
      {
        const double middle  = stretch.Start + (quarter * offset);
        total               += half
                             * Length(Velocity(points, middle - (half * node))
                                      + Velocity(points, middle + (half * node)));
      }
      return total * (1.0 - 0x1p-40);
    }

    /**
     *  The length a hodograph that is a straight line runs over a stretch, in closed form.
     *
     *  For h(t) = h0 + a t the speed is |a| sqrt(u^2 + q^2), u = t - p, p the parameter nearest the
     *  root and q its distance from the parameter's line, both over |a|. Its integral from u0 to u1
     *  is |a| / 2 (u1 s1 - u0 s0 + q^2 (asinh(u1 / q) - asinh(u0 / q))), s = sqrt(u^2 + q^2). Where
     *  u0 and u1 have one sign, both differences cancel, and are rationalised:
     *
     *  - u1 s1 - u0 s0 = (u1 - u0) (u1 + u0) (u1^2 + u0^2 + q^2) / (u1 s1 + u0 s0);
     *  - asinh(u1 / q) - asinh(u0 / q) = asinh((u1 - u0) (u1 + u0) / (u1 s0 + u0 s1));
     *
     *  u1 - u0 the stretch's width, as exact as it is given; where they differ in sign, each term
     *  is a sum. It is worked out in D = a . h, X = |a x h| and S = |h| at the ends, u = D / |a|^2,
     *  q = X / |a|^2 and s = S / |a|, so that no quotient rounds before the last. A hodograph whose
     *  change over the stretch is within 2^-30 of its speed at the middle runs at that speed, off
     *  by 2^-60 / 12 of the length at most.
     *
     *  @param  start   The hodograph at 0.
     *  @param  end     The hodograph at 1, not both of them nothing.
     *  @param  stretch The stretch, its start below its end.
     *  @return         The length.
     */
    auto Linear(Vector start, Vector end, Range stretch) -> double
    {
      constexpr double Tiny  = 0x1p-400;
      const int        power = -(std::ilogb(std::max(Largest(start), Largest(end))) + 1);
      const Vector     low   = {.X = std::ldexp(start.X, power), .Y = std::ldexp(start.Y, power)};
      const Vector     high  = {.X = std::ldexp(end.X, power), .Y = std::ldexp(end.Y, power)};
      const auto       velocity = [low, high](double along) -> Vector
      {
        return (low * (1.0 - along)) + (high * along);
      };
      const Vector slope  = high - low;
      const double width  = stretch.End - stretch.Start;
      const double square = Dot(slope, slope);
      const double rate   = std::sqrt(square);
      const double speed  = Length(velocity(stretch.Start + (width / 2.0)));
      if (!(rate * width > 0x1p-30 * speed))
      {
        return std::ldexp(width * speed, -power);
      }
      const Vector opening   = velocity(stretch.Start);
      const Vector closing   = velocity(stretch.End);
      const double cross     = std::abs(Cross(slope, opening));
      double       near      = Dot(slope, opening);
      double       far       = Dot(slope, closing);
      double       nearSpeed = Length(opening);
      double       farSpeed  = Length(closing);
      if (!(far > 0.0))
      {
        const double turned = near;
        near                = -far;
        far                 = -turned;
        std::swap(nearSpeed, farSpeed);
      }
      const double curl  = cross * cross / (2.0 * square * rate);
      double       whole = 0.0;
      if (!(near < 0.0))
      {
        const double sum = far + near;
        whole            = width
                         * sum
                         * (((far * far) + (near * near)) + (cross * cross))
                         / (2.0 * square * ((far * farSpeed) + (near * nearSpeed)));
        if (cross > Tiny)
        {
          whole += curl * Trig::Asinh(rate * width * sum / ((far * nearSpeed) + (near * farSpeed)));
        }
      }
      else
      {
        whole = ((far * farSpeed) - (near * nearSpeed)) / (2.0 * square);
        if (cross > Tiny)
        {
          whole += curl * (Trig::Asinh(far / cross) + Trig::Asinh(-near / cross));
        }
      }
      return std::ldexp(whole, -power);
    }

    /**
     *  A hodograph of degree 2, ready to integrate: its control points, scaled, its leading
     *  coefficient's length on that scale, and its roots.
     */
    struct Parabola
    {
      std::array<Vector, 3> Points = {};  /**< The control points, scaled.       */
      double                Lead   = 0.0; /**< The leading coefficient's length. */
      std::array<Vector, 2> Roots  = {};  /**< The roots, as complex numbers.    */
    };

    /**
     *  The distance of a point from the origin, by the square root of the sum of squares: past
     *  1e154 it is without end, which counts a root that far as unreachable, as it is.
     *
     *  @param  x Its first coordinate.
     *  @param  y Its second.
     *  @return   The distance.
     */
    auto Modulus(double x, double y) -> double
    {
      return std::sqrt((x * x) + (y * y));
    }

    /**
     *  How thin the Bernstein ellipse about a piece is, and where to cut the piece if it is too
     *  thin.
     */
    struct Slimness
    {
      double Rho   = 0.0; /**< The ellipse's rho.                    */
      double Split = 0.0; /**< Where to cut the piece if it is thin. */
    };

    /**
     *  How thin the Bernstein ellipse about a piece is that reaches the nearest root counted: rho =
     *  a + sqrt(a^2 - 1), a the sum of the root's distances from the piece's ends over its width.
     *
     *  A root is not counted where it is taken as real: its imaginary part so small, against the
     *  piece's share of the error, that the speed is |t - p| times what is analytic, but for a
     *  quarter of the share (see Measured), log(1 + 2 w / q) taken at its most, 43, q being at
     *  least 2^-60 w.
     *
     *  The cut is 1.875 times the nearest root's distance from the end it lies beyond, from that
     *  end: a piece that long has rho at least 3 whichever way the root lies from the end, so the
     *  pieces grow away from a root by a factor of nearly 3. It is kept an eighth of the piece in
     *  from either end, so that a piece far wider than its root's distance shrinks by 8 at a time.
     *
     *  @param  parabola The hodograph.
     *  @param  piece    The piece, no root's real part inside it.
     *  @param  share    The piece's share of the error.
     *  @return          The ellipse's rho, without end when no root is counted; and the cut.
     */
    auto Thickness(const Parabola &parabola, Range piece, double share) -> Slimness
    {
      constexpr double Smallest = std::numeric_limits<double>::min();
      const double     width    = piece.End - piece.Start;
      const double     half     = width / 2.0;
      const double     middle   = piece.Start + half;
      Slimness         reach    = {.Rho = Infinite, .Split = middle};
      for (std::size_t index = 0; index < 2; ++index)
      {
        const Vector root  = parabola.Roots.at(index);
        const Vector other = parabola.Roots.at(1 - index);
        const double lift  = std::max({std::abs(root.Y), half * 0x1p-60, Smallest});
        const double cost =
            parabola.Lead * (Modulus(middle - other.X, other.Y) + half) * lift * lift * 683.0;
        const double before = Modulus(piece.Start - root.X, root.Y);
        const double after  = Modulus(piece.End - root.X, root.Y);
        const double focal  = (before + after) / std::max(width, Smallest);
        const double rho    = focal + std::sqrt(std::max(0.0, (focal * focal) - 1.0));
        if (cost > share / 4.0 && rho < reach.Rho)
        {
          const double cut =
              root.X <= piece.Start ? piece.Start + (1.875 * before) : piece.End - (1.875 * after);
          reach = {.Rho = rho,
                   .Split =
                       std::clamp(cut, piece.Start + (width / 8.0), piece.End - (width / 8.0))};
        }
      }
      return reach;
    }

    /**
     *  What a piece is integrated to: how thin its ellipse is, and its share of the error.
     */
    struct Allowance
    {
      double Rho   = 0.0; /**< The ellipse's rho.          */
      double Share = 0.0; /**< The piece's share of error. */
    };

    /**
     *  A piece's length by the fewest pairs of Gauss-Legendre points whose error bound is within
     *  three quarters of its share, or 40 points. A piece whose ellipse is thinner than RhoMin, its
     *  halvings spent, takes the bound of one 1.25 thick, which is no bound for it.
     *
     *  The rule of n points is exact to degree 2 n - 1, so its error is (64 / 15) M w rho^(2 - 2 n)
     *  / (rho^2 - 1) at most (Trefethen's theorem 19.3, whose n + 1 points make 2 n here). Each
     *  node is placed from the piece's nearer end, and the sum times the half width held exactly:
     *  nodes about a rounded middle, a rounded half width apart, would stretch the piece by a step
     *  of a double at its ends and add the speed there times that step.
     *
     *  @param  parabola  The hodograph.
     *  @param  piece     The piece.
     *  @param  allowance Its ellipse and its share of the error.
     *  @return           The length.
     */
    auto GaussLegendre(const Parabola &parabola, Range piece, Allowance allowance)
        -> Twofold::Double
    {
      const double rho    = std::clamp(allowance.Rho, 1.25, 64.0);
      const double half   = (piece.End - piece.Start) / 2.0;
      const double middle = piece.Start + half;
      const double reach  = half * (rho + (1.0 / rho)) / 2.0;
      double       most   = parabola.Lead;
      for (const Vector root : parabola.Roots)
      {
        most *= Modulus(middle - root.X, root.Y) + reach;
      }
      const double fall  = 1.0 / (rho * rho * rho * rho);
      double       bound = 64.0 / 15.0 * most * half / ((rho * rho) - 1.0) / (rho * rho);
      std::size_t  pairs = 1;
      while (bound > 0.75 * allowance.Share && pairs < Legendre::Pairs)
      {
        bound *= fall;
        ++pairs;
      }
      const Legendre::Rule rule = Legendre::Ruled(pairs);
      Twofold::Double      sum  = {};
      for (std::size_t index = 0; index < pairs; ++index)
      {
        const double inset  = half * (1.0 - rule.Nodes[index]);
        const Vector before = Velocity(parabola.Points, piece.Start + inset);
        const Vector after  = Velocity(parabola.Points, piece.End - inset);
        sum = Twofold::Add(sum, Twofold::Lifted(rule.Weights[index]
                                                * (std::sqrt(Dot(before, before))
                                                   + std::sqrt(Dot(after, after)))));
      }
      const Twofold::Two across = Twofold::TwoSum(piece.End, -piece.Start);
      return Twofold::Multiply(sum, {.High = across.High / 2.0, .Low = across.Low / 2.0});
    }

    /**
     *  The length a hodograph of degree 2 runs over a stretch: cut at its roots' real parts, each
     *  piece halved while its ellipse is thinner than RhoMin, each then by the fewest points its
     *  bound allows, the lengths summed with their roundoffs.
     *
     *  @param  parabola The hodograph.
     *  @param  stretch  The stretch, its start below its end.
     *  @param  effort   The error, on the hodograph's scale, and the halvings to spend.
     *  @return          The length, on that scale, as a double and its roundoff.
     */
    auto Curved(const Parabola &parabola, Range stretch, Effort effort) -> Twofold::Double
    {
      const double width = stretch.End - stretch.Start;
      const double target =
          std::max(effort.Within.Error, 0x1p-55 * Lower(parabola.Points, stretch));
      std::array<double, 4> cuts = {stretch.Start, stretch.End, stretch.Start, stretch.Start};
      for (std::size_t index = 0; index < 2; ++index)
      {
        const double real  = parabola.Roots.at(index).X;
        cuts.at(2 + index) = stretch.Start < real && real < stretch.End ? real : stretch.Start;
      }
      std::ranges::sort(cuts);
      std::array<Range, Capacity> pending{};
      std::size_t                 count = 0;
      for (std::size_t index = cuts.size() - 1; index > 0; --index)
      {
        if (cuts.at(index - 1) < cuts.at(index))
        {
          pending.at(count++) = {.Start = cuts.at(index - 1), .End = cuts.at(index)};
        }
      }
      const std::size_t budget = std::min(effort.Halvings, Halved);
      std::size_t       spent  = 0;
      Twofold::Double   total  = {};
      while (count > 0)
      {
        const Range    piece = pending.at(--count);
        const double   share = target * ((piece.End - piece.Start) / width);
        const Slimness reach = Thickness(parabola, piece, share);
        if (reach.Rho < RhoMin && spent < budget)
        {
          pending.at(count++) = {.Start = reach.Split, .End = piece.End};
          pending.at(count++) = {.Start = piece.Start, .End = reach.Split};
          ++spent;
          continue;
        }
        total =
            Twofold::Add(total, GaussLegendre(parabola, piece, {.Rho = reach.Rho, .Share = share}));
      }
      return total;
    }

    /**
     *  A hodograph's control points times a power of two.
     *
     *  @param  hodograph The hodograph.
     *  @param  power     The power.
     *  @return           The points, scaled.
     */
    auto Magnified(const Hodograph &hodograph, int power) -> std::array<Vector, 3>
    {
      std::array<Vector, 3> points{};
      for (std::size_t index = 0; index < points.size(); ++index)
      {
        const Vector point = hodograph.Points.at(index);
        points.at(index)   = {.X = std::ldexp(point.X, power), .Y = std::ldexp(point.Y, power)};
      }
      return points;
    }

    /**
     *  The length a hodograph of degree 2 runs over a stretch.
     *
     *  Its roots come from its points scaled near 1; its speeds from them scaled again by the
     *  square of the power of two past the stretch's farthest end, up to 2^500, so that no square
     *  of a velocity on it overflows. Where its leading coefficient is 0 it is straight.
     *
     *  @param  hodograph The hodograph, of degree 2.
     *  @param  stretch   The stretch, finite, its start below its end.
     *  @param  effort    The error and the halvings to spend.
     *  @return           The length.
     */
    auto Parabolic(const Hodograph &hodograph, Range stretch, Effort effort) -> double
    {
      double largest = 0.0;
      for (const Vector point : hodograph.Points)
      {
        largest = std::max(largest, Largest(point));
      }
      const int                   unit   = -(std::ilogb(largest) + 1);
      const std::array<Vector, 3> points = Magnified(hodograph, unit);
      const Power power = {.Square   = (points[2] - points[1]) - (points[1] - points[0]),
                           .Linear   = (points[1] - points[0]) * 2.0,
                           .Constant = points[0]};
      if (Zero(Largest(power.Square)))
      {
        return hodograph.Factor * Linear(hodograph.Points[0], hodograph.Points[2], stretch);
      }
      const double   farthest = std::max({1.0, std::abs(stretch.Start), std::abs(stretch.End)});
      const int      far      = -2 * std::min(std::ilogb(farthest) + 1, 500);
      const Parabola parabola = {.Points = Magnified(hodograph, unit + far),
                                 .Lead   = std::ldexp(Length(power.Square), far),
                                 .Roots  = Solutions(power)};
      const Effort   scaled   = {
          .Within   = {.Error = std::ldexp(effort.Within.Error / hodograph.Factor, unit + far)},
          .Halvings = effort.Halvings};
      const Twofold::Double length =
          Twofold::Multiply(Curved(parabola, stretch, scaled), Twofold::Lifted(hodograph.Factor));
      return std::ldexp(length.High + length.Low, -(unit + far));
    }
  } // namespace

  auto Hodographed(const Bezier::Quadratic &quadratic) -> Hodograph
  {
    return {.Degree = 1,
            .Points = {quadratic.Control - quadratic.Start, quadratic.End - quadratic.Control, {}},
            .Factor = 2.0};
  }

  auto Hodographed(const Bezier::Cubic &cubic) -> Hodograph
  {
    return {.Degree = 2,
            .Points = {cubic.Out - cubic.Start, cubic.In - cubic.Out, cubic.End - cubic.In},
            .Factor = 3.0};
  }

  auto Measured(const Hodograph &hodograph, Range stretch, Effort effort) -> double
  {
    double largest = 0.0;
    for (const Vector point : std::span{hodograph.Points}.first(hodograph.Degree + 1))
    {
      largest = std::max(largest, Largest(point));
    }
    const double width = stretch.End - stretch.Start;
    if (!(largest > 0.0) || !(width > 0.0))
    {
      return 0.0;
    }
    if (!std::isfinite(width) || !std::isfinite(largest))
    {
      return Infinite;
    }
    if (hodograph.Degree == 0)
    {
      return hodograph.Factor * width * Length(hodograph.Points[0]);
    }
    if (hodograph.Degree == 1)
    {
      return hodograph.Factor * Linear(hodograph.Points[0], hodograph.Points[1], stretch);
    }
    return Parabolic(hodograph, stretch, effort);
  }
} // namespace Geoxy
