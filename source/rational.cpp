/***************************************************************************************************
 *  @file      source/rational.cpp
 *  @brief     A spline's length, area and moments strip by strip: its spans cut where their weights
 *             crowd an end, each strip's integrals by Gauss-Legendre rules whose error is bounded
 *             from the roots of its speed and of its weight.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "rational.hpp"

#include "aberth.hpp"
#include "basis.hpp"
#include "eased.hpp"
#include "hodograph.hpp"
#include "kernel/piece.hpp"
#include "legendre.hpp"
#include "quadrature.hpp"
#include "real.hpp"
#include "trig.hpp"
#include "twofold.hpp"

#include "geoxy/at.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace Geoxy
{
  namespace
  {
    constexpr double Smallest = std::numeric_limits<double>::min();

    /**
     *  The thinnest Bernstein ellipse a piece is integrated over: thinner, it is halved.
     */
    constexpr double RhoFloor = 3.0;

    /**
     *  How near an end of a piece, in its own parameter, a root of its weight may lie and be told
     *  apart from the end to some 28 bits; nearer, the piece is cut that far from the end.
     */
    constexpr double Close = 0x1p-24;

    /**
     *  How many times a span is cut at most on its way to strips: Close to the 42nd is 2^-1008, as
     *  narrow as a strip may be and its width in its span still a double of all its digits.
     */
    constexpr int MostCuts = 42;

    /**
     *  How many pieces a stretch is held in at once: one between each two cuts, and one more per
     *  halving.
     */
    constexpr std::size_t StackDepth = Halved + Kernel::Terms + 4;

    /**
     *  The most coefficients a polynomial here has, and so the rows of binomials it needs.
     */
    constexpr std::size_t Rows = Aberth::Most + 1;

    /**
     *  The binomial coefficients of every degree under Rows, each row as the multiplicative formula
     *  gives it: degree choose 0 to degree choose degree.
     */
    constexpr std::array<std::array<double, Rows>, Rows> Binomials =
        [] -> std::array<std::array<double, Rows>, Rows>
    {
      std::array<std::array<double, Rows>, Rows> table = {};
      for (std::size_t degree = 0; degree < Rows; ++degree)
      {
        std::array<double, Rows> &row = table.at(degree);
        row.at(0)                     = 1.0;
        for (std::size_t index = 1; index <= degree; ++index)
        {
          row.at(index) = row.at(index - 1)
                        * static_cast<double>(degree - index + 1)
                        / static_cast<double>(index);
        }
      }
      return table;
    }();

    /**
     *  A parameter and 1 less it, each as exactly as it can be told: a node near 1 is placed by how
     *  far short of 1 it is, which a double near 1 tells to a step of 2^-53 and no finer.
     */
    struct Abscissa
    {
      double Along = 0.0; /**< The parameter t. */
      double Rest  = 1.0; /**< 1 - t.           */
    };

    /**
     *  The ratio a polynomial in Bernstein form is summed in at a parameter: t / (1 - t) where t is
     *  nearer 0 than 1, and (1 - t) / t elsewhere, at most 1 in size.
     *
     *  @param  spot The parameter.
     *  @return      The ratio.
     */
    auto Leanness(Abscissa spot) -> double
    {
      return std::abs(spot.Along) <= std::abs(spot.Rest) ? spot.Along / spot.Rest
                                                         : spot.Rest / spot.Along;
    }

    /**
     *  A polynomial in Bernstein form at a parameter, times (1 + r) to its degree, r the ratio
     *  Leanness gives: Horner's rule in r, the coefficients taken from the end the parameter is
     *  nearer, so that nothing overflows however far past 0 and 1 the parameter lies. Divided by (1
     *  + r) to its degree, it is the polynomial at r / (1 + r), or at 1 / (1 + r) where the sum is
     *  taken from 1, exactly but for the sum's rounding: at a parameter a step of a double from the
     *  one asked, as a node's is.
     *
     *  @tparam ValueT What a coefficient is: a double or a vector.
     *  @param  values The coefficients, one at least.
     *  @param  spot   The parameter.
     *  @return        The polynomial there, times (1 + r) to its degree.
     */
    template <typename ValueT>
    auto Hornered(std::span<const ValueT> values, Abscissa spot) -> ValueT
    {
      const std::size_t               degree  = values.size() - 1;
      const bool                      forward = std::abs(spot.Along) <= std::abs(spot.Rest);
      const double                    ratio   = Leanness(spot);
      const std::array<double, Rows> &row     = Binomials.at(degree);
      ValueT                          sum     = values[forward ? degree : 0] * row.at(degree);
      for (std::size_t index = degree; index-- > 0;)
      {
        sum = (sum * ratio) + (values[forward ? index : degree - index] * row.at(index));
      }
      return sum;
    }

    /**
     *  (1 + r) to a power at a parameter, r the ratio Leanness gives, in double-double and rounded
     *  once: what Hornered's sums are divided by, to their degrees, for the polynomials' values.
     *  Taken from r itself, it is the base of the very parameter the sums were taken at, whose
     *  rounding is not grown by the power, as the power of a 1 - t rounded apart from t would be.
     *
     *  @param  spot  The parameter.
     *  @param  power The power, negative too.
     *  @return       (1 + r) to the power.
     */
    auto Powered(Abscissa spot, int power) -> double
    {
      const Twofold::Two sum    = Twofold::TwoSum(1.0, Leanness(spot));
      Twofold::Double    factor = {.High = sum.High, .Low = sum.Low};
      Twofold::Double    result = Twofold::Lifted(1.0);
      auto               left   = static_cast<unsigned>(power < 0 ? -power : power);
      while (left > 0)
      {
        if ((left & 1U) != 0)
        {
          result = Twofold::Multiply(result, factor);
        }
        factor   = Twofold::Multiply(factor, factor);
        left   >>= 1U;
      }
      const Twofold::Double grown =
          power < 0 ? Twofold::Divide(Twofold::Lifted(1.0), result) : result;
      return grown.High + grown.Low;
    }

    /**
     *  The distance of a point from the origin, by the square root of the sum of squares: past
     *  1e154 it is without end, which counts a root that far as unreachable, as it is.
     *
     *  @param  x Its first coordinate.
     *  @param  y Its second.
     *  @return   The distance.
     */
    auto Norm(double x, double y) -> double
    {
      return std::sqrt((x * x) + (y * y));
    }

    /**
     *  The log2 of the product of some roots' distances from a point on the real line, carried as a
     *  fraction and a power of two so that it neither overflows nor underflows.
     *
     *  @param  roots The roots.
     *  @param  point The point.
     *  @return       The log2 of the product.
     */
    auto Distanced(const Aberth::Roots &roots, double point) -> double
    {
      double fraction = 1.0;
      int    power    = 0;
      for (const Vector root : std::span{roots.Values}.first(roots.Count))
      {
        int exponent  = 0;
        fraction      = std::frexp(fraction * Norm(point - root.X, root.Y), &exponent);
        power        += exponent;
      }
      return static_cast<double>(power) + Trig::Log2(fraction);
    }

    /**
     *  The log2 of a bound from below on a polynomial over a piece's Bernstein ellipse of a half
     *  axis: its leading coefficient's length times the product of each root's own half axis less
     *  that one, at least the smallest double.
     *
     *  @param  singular The polynomial's roots and leading coefficient.
     *  @param  piece    The piece.
     *  @param  reach    The ellipse's half axis, at least half the piece's width.
     *  @return          The log2 of the bound.
     */
    auto Minorant(const Singular &singular, Range piece, double reach) -> double
    {
      double fraction = 1.0;
      int    power    = 0;
      for (const Vector root : std::span{singular.Roots.Values}.first(singular.Roots.Count))
      {
        const double axis =
            (Norm(piece.Start - root.X, root.Y) + Norm(piece.End - root.X, root.Y)) / 2.0;
        int exponent  = 0;
        fraction      = std::frexp(fraction * std::max(axis - reach, Smallest), &exponent);
        power        += exponent;
      }
      return singular.Lead + static_cast<double>(power) + Trig::Log2(fraction);
    }

    /**
     *  The log2 of a bound on a strip's speed over a piece's Bernstein ellipse of a half axis, R:
     *  on it |N| is at most |a| times the product of each of its roots' distance from the piece's
     *  middle plus R, and |W| at least |b| times the product of each of its roots' own half axis
     *  less R, at least the smallest double; and the speed is |N| / W^2 times the factor.
     *
     *  @param  pace  The speed.
     *  @param  piece The piece.
     *  @param  reach The ellipse's half axis R, at least half the piece's width.
     *  @return       The log2 of the bound.
     */
    auto Most(const Pace &pace, Range piece, double reach) -> double
    {
      const double middle   = piece.Start + ((piece.End - piece.Start) / 2.0);
      double       fraction = pace.Factor;
      int          power    = 0;
      int          exponent = 0;
      for (const Vector root : std::span{pace.Turns.Roots.Values}.first(pace.Turns.Roots.Count))
      {
        fraction  = std::frexp(fraction * (Norm(middle - root.X, root.Y) + reach), &exponent);
        power    += exponent;
      }
      for (const Vector root : std::span{pace.Poles.Roots.Values}.first(pace.Poles.Roots.Count))
      {
        const double gap = std::max(
            ((Norm(piece.Start - root.X, root.Y) + Norm(piece.End - root.X, root.Y)) / 2.0) - reach,
            Smallest);
        fraction  = std::frexp(fraction / gap, &exponent);
        power    += exponent;
        fraction  = std::frexp(fraction / gap, &exponent);
        power    += exponent;
      }
      return pace.Turns.Lead
           - (2.0 * pace.Poles.Lead)
           + static_cast<double>(power)
           + Trig::Log2(fraction);
    }

    /**
     *  A polynomial's roots, from its Bernstein form over 0 to 1, and its leading coefficient's
     *  length as its value at the one of 0, 1/4, 1/2, 3/4 and 1 where it is largest makes it: there
     *  the product of the distances to the roots, times it, is the polynomial's length, so that the
     *  bounds Distanced and Minorant give from the roots hold the polynomial as it is near the
     *  piece, even where its degree is less than its form's and an approximation has gone off far.
     *
     *  @param  values The coefficients, complex, one at least.
     *  @return        The roots in t and the log2 of the leading coefficient's length: without end
     *                 below, and no roots, for the polynomial 0.
     */
    auto Singularities(std::span<const Vector> values) -> Singular
    {
      const std::size_t degree   = values.size() - 1;
      Singular          singular = {.Roots = {}, .Lead = -Infinite};
      if (!std::ranges::any_of(values, [](Vector value) -> bool { return !Zero(Largest(value)); }))
      {
        return singular;
      }
      if (degree > 0)
      {
        singular.Roots = Aberth::Solved(values, Aberth::Rounds);
      }
      double best = 0.0;
      for (const double along : {0.0, 0.25, 0.5, 0.75, 1.0})
      {
        const Abscissa spot = {.Along = along, .Rest = 1.0 - along};
        const double   value =
            Length(Hornered(values, spot)) * Powered(spot, -static_cast<int>(degree));
        if (value > best)
        {
          best          = value;
          singular.Lead = Trig::Log2(value) - Distanced(singular.Roots, along);
        }
      }
      return singular;
    }

    /**
     *  A piece's weight's roots.
     *
     *  @param  piece The piece.
     *  @return       Its weight's roots in its own parameter, and the log2 of the leading
     *                coefficient's length.
     */
    auto WeightRoots(const Kernel::Piece &piece) -> Singular
    {
      std::array<Vector, Rows> weights{};
      for (std::size_t index = 0; index <= piece.Degree; ++index)
      {
        weights.at(index) = {.X = piece.Poles.at(index).W, .Y = 0.0};
      }
      return Singularities(std::span<const Vector>{weights}.first(piece.Degree + 1));
    }

    /**
     *  How thin the Bernstein ellipse about a piece is that reaches a root: rho = a + sqrt(a^2 -
     *  1), a the sum of the root's distances from the piece's ends over its width.
     *
     *  @param  piece The piece.
     *  @param  root  The root.
     *  @return       The ellipse's rho.
     */
    auto Thinness(Range piece, Vector root) -> double
    {
      const double focal = (Norm(piece.Start - root.X, root.Y) + Norm(piece.End - root.X, root.Y))
                         / std::max(piece.End - piece.Start, Smallest);
      return focal + std::sqrt(std::max(0.0, (focal * focal) - 1.0));
    }

    /**
     *  Where to cut a piece thin about a root: 1.875 times the root's distance from the end it lies
     *  beyond, from that end, where a piece that long has rho at least 3 whichever way the root
     *  lies, kept an eighth of the piece in from either end.
     *
     *  @param  piece The piece.
     *  @param  root  The root, its real part not inside the piece.
     *  @return       The cut.
     */
    auto Cleft(Range piece, Vector root) -> double
    {
      const double width = piece.End - piece.Start;
      const double cut   = root.X <= piece.Start
                             ? piece.Start + (1.875 * Norm(piece.Start - root.X, root.Y))
                             : piece.End - (1.875 * Norm(piece.End - root.X, root.Y));
      return std::clamp(cut, piece.Start + (width / 8.0), piece.End - (width / 8.0));
    }

    /**
     *  How thin the ellipses about a piece are, and where to cut it if it is too thin.
     */
    struct Girth
    {
      double Rho   = Infinite; /**< The thinnest of any root counted.     */
      double Pole  = Infinite; /**< The thinnest of the weight's roots.   */
      double Split = 0.0;      /**< Where to cut the piece if it is thin. */
    };

    /**
     *  How thin the ellipses about a piece are that reach the weight's roots.
     *
     *  @param  poles The weight's roots.
     *  @param  piece The piece.
     *  @return       The thinnest, and where to cut.
     */
    auto Thinnest(const Singular &poles, Range piece) -> Girth
    {
      Girth reach = {.Rho   = Infinite,
                     .Pole  = Infinite,
                     .Split = piece.Start + ((piece.End - piece.Start) / 2.0)};
      for (const Vector root : std::span{poles.Roots.Values}.first(poles.Roots.Count))
      {
        const double rho = Thinness(piece, root);
        if (rho < reach.Rho)
        {
          reach = {.Rho = rho, .Pole = rho, .Split = Cleft(piece, root)};
        }
      }
      return reach;
    }

    /**
     *  How thin the ellipses about a piece are that reach the roots of the speed counted, and of
     *  its weight.
     *
     *  A root of the speed's numerator is not counted where it is taken as real: its imaginary part
     *  q so small, against the piece's share of the error, that the speed is |t - p| times what is
     *  analytic, but for a quarter of the share (see Measured), the rest of the speed bounded on
     *  the piece by the other roots and the weight's, log(1 + 2 h / q) taken at its most, 43, q
     *  being at least 2^-60 h, and the rule's 568 added: 683 times q^2 that bound.
     *
     *  @param  pace  The speed.
     *  @param  piece The piece, no root's real part inside it.
     *  @param  share The piece's share of the error.
     *  @return       The thinnest, and where to cut.
     */
    auto Girthed(const Pace &pace, Range piece, double share) -> Girth
    {
      const double half   = (piece.End - piece.Start) / 2.0;
      const double middle = piece.Start + half;
      Girth        reach  = Thinnest(pace.Poles, piece);
      const double rest   = Most(pace, piece, half);
      const double floor  = Trig::Exp2(Trig::Log2(share / 4.0 / 683.0) - rest);
      for (const Vector root : std::span{pace.Turns.Roots.Values}.first(pace.Turns.Roots.Count))
      {
        const double lift = std::max({std::abs(root.Y), half * 0x1p-60, Smallest});
        const double rho  = Thinness(piece, root);
        if (lift * lift > floor * (Norm(middle - root.X, root.Y) + half) && rho < reach.Rho)
        {
          reach.Rho   = rho;
          reach.Split = Cleft(piece, root);
        }
      }
      return reach;
    }

    /**
     *  The fewest pairs of Gauss-Legendre points whose bound on a piece is within three quarters of
     *  its share of the error: (64 / 15) M h rho^(2 - 2 n) / (rho^2 - 1) for n points.
     *
     *  @param  most  The log2 of M h.
     *  @param  rho   The ellipse's rho, from 1.25 to 64.
     *  @param  share The log2 of the share.
     *  @return       The pairs, from 1 to Legendre::Pairs; the most where nothing bounds the error.
     */
    auto Needed(double most, double rho, double share) -> double
    {
      const double fall = Trig::Log2(rho);
      const double first =
          (most - share) - Trig::Log2(((rho * rho) - 1.0) * rho * rho * 0.75 * 15.0 / 64.0);
      const double pairs = std::ceil(1.0 + (first / (4.0 * fall)));
      return std::fmax(1.0, std::fmin(pairs, static_cast<double>(Legendre::Pairs)));
    }

    /**
     *  The rhos a piece may be integrated over: none past the nearest root counted, whose rho
     *  itself is one where that root is the speed's, and the powers 9 / 10, 3 / 4 and 1 / 2 of the
     *  nearest of the weight's, which M grows without end toward; each from 1.25 to 64.
     *
     *  @param  reach How thin the ellipses about the piece are.
     *  @return       The rhos.
     */
    auto Candidates(Girth reach) -> std::array<double, 4>
    {
      const auto toward = [reach](double power) -> double
      {
        return std::min(Trig::Pow(reach.Pole, power), reach.Rho);
      };
      std::array<double, 4> rhos = {reach.Rho < reach.Pole ? reach.Rho : toward(0.9), toward(0.9),
                                    toward(0.75), toward(0.5)};
      for (double &rho : rhos)
      {
        rho = std::clamp(rho, 1.25, 64.0);
      }
      return rhos;
    }

    /**
     *  What a piece is integrated to: how thin its ellipses are, and its share of the error.
     */
    struct Allotment
    {
      Girth  Thick = {};  /**< How thin its ellipses are. */
      double Share = 0.0; /**< Its share of the error.    */
    };

    /**
     *  Integrals over a piece by a rule of pairs of Gauss-Legendre points, each node placed from
     *  the piece's nearer end, its distance from 1 too, and each sum times the half width held
     *  exactly: nodes about a rounded middle, a rounded half width apart, would stretch the piece
     *  by a step of a double at its ends.
     *
     *  @tparam FunctionT What gives the functions' values at a parameter, an array of doubles.
     *  @param  function  The functions.
     *  @param  piece     The piece.
     *  @param  pairs     How many pairs of points.
     *  @return           The integrals, each with its roundoff.
     */
    template <typename FunctionT>
    auto Ruled(const FunctionT &function, Range piece, std::size_t pairs)
    {
      constexpr std::size_t              Count = decltype(function(Abscissa{})){}.size();
      const double                       half  = (piece.End - piece.Start) / 2.0;
      const double                       start = 1.0 - piece.Start;
      const double                       end   = 1.0 - piece.End;
      const Legendre::Rule               rule  = Legendre::Ruled(pairs);
      std::array<Twofold::Double, Count> sums{};
      for (std::size_t index = 0; index < pairs; ++index)
      {
        const double inset  = half * (1.0 - rule.Nodes[index]);
        const auto   before = function({.Along = piece.Start + inset, .Rest = start - inset});
        const auto   after  = function({.Along = piece.End - inset, .Rest = end + inset});
        for (std::size_t which = 0; which < Count; ++which)
        {
          sums.at(which) =
              Twofold::Add(sums.at(which), Twofold::Lifted(rule.Weights[index]
                                                           * (before.at(which) + after.at(which))));
        }
      }
      const Twofold::Two across = Twofold::TwoSum(piece.End, -piece.Start);
      for (Twofold::Double &sum : sums)
      {
        sum = Twofold::Multiply(sum, {.High = across.High / 2.0, .Low = across.Low / 2.0});
      }
      return sums;
    }

    /**
     *  A strip's speed at a parameter.
     *
     *  @param  pace The speed.
     *  @param  spot The parameter.
     *  @return      |N| / W^2 there.
     */
    auto Speed(const Pace &pace, Abscissa spot) -> double
    {
      const Vector top = Hornered(std::span<const Vector>{pace.Top}.first(pace.Degree + 1), spot);
      const double bottom =
          Hornered(std::span<const double>{pace.Bottom}.first(pace.Weighing + 1), spot);
      const int power = static_cast<int>(pace.Degree) - (2 * static_cast<int>(pace.Weighing));
      return Length(top) / (bottom * bottom) * Powered(spot, -power) * pace.Factor;
    }

    /**
     *  A strip's point at a parameter, less its offset.
     *
     *  @param  pace The speed, which holds the strip's control points.
     *  @param  spot The parameter.
     *  @return      The point, less the offset.
     */
    auto Plotted(const Pace &pace, Abscissa spot) -> Vector
    {
      const std::size_t count = pace.Order + 1;
      return Hornered(std::span<const Vector>{pace.Places}.first(count), spot)
           / Hornered(std::span<const double>{pace.Masses}.first(count), spot);
    }

    /**
     *  A lower bound on the length a strip runs over a stretch: the polygon through its points at
     *  the stretch's ends and between them, less 2^-40 of it for the rounding of the points.
     *
     *  @param  pace    The speed.
     *  @param  stretch The stretch.
     *  @param  sides   How many sides the polygon has.
     *  @return         The lower bound.
     */
    auto Inscribed(const Pace &pace, Range stretch, int sides) -> double
    {
      const double step  = (stretch.End - stretch.Start) / static_cast<double>(sides);
      const double rest  = 1.0 - stretch.Start;
      Vector       last  = Plotted(pace, {.Along = stretch.Start, .Rest = rest});
      double       total = 0.0;
      for (int side = 1; side <= sides; ++side)
      {
        const double inset  = step * static_cast<double>(side);
        const Vector next   = Plotted(pace, {.Along = stretch.Start + inset, .Rest = rest - inset});
        total              += Length(next - last);
        last                = next;
      }
      return std::ldexp(total * (1.0 - 0x1p-40), pace.Size);
    }

    /**
     *  A piece's length by the fewest pairs of points its bound allows, over whichever of the rhos
     *  open to it takes the fewest.
     *
     *  @param  pace      The speed.
     *  @param  piece     The piece.
     *  @param  allowance How thin its ellipses are, and its share of the error.
     *  @return           The length.
     */
    auto Lengthened(const Pace &pace, Range piece, Allotment allowance) -> Twofold::Double
    {
      const double half  = (piece.End - piece.Start) / 2.0;
      const double share = Trig::Log2(allowance.Share);
      auto         pairs = static_cast<double>(Legendre::Pairs);
      double       last  = 0.0;
      for (const double rho : Candidates(allowance.Thick))
      {
        if (!Exactly(rho, last))
        {
          const double most =
              Most(pace, piece, half * (rho + (1.0 / rho)) / 2.0) + Trig::Log2(half);
          pairs = std::min(pairs, Needed(most, rho, share));
        }
        last = rho;
      }
      return Ruled([&pace](Abscissa spot) -> std::array<double, 1> { return {Speed(pace, spot)}; },
                   piece, static_cast<std::size_t>(pairs))[0];
    }

    /**
     *  A stretch cut where the real parts of roots near it fall inside it, the pieces in order from
     *  its end, as a stack takes them: a root whose ellipse about the stretch is thinner than
     *  RhoFloor, which the stretch would be halved about anyway, and so is a root that may be taken
     *  as real, the speed's turn back or cusp, or the weight's peak; one farther off needs no cut.
     *
     *  @param  roots   The roots, of the speed's numerator and of the weight.
     *  @param  stretch The stretch.
     *  @param  pending Where the pieces go.
     *  @return         How many.
     */
    auto Splayed(std::array<const Singular *, 2> roots, Range stretch,
                 std::array<Range, StackDepth> &pending) -> std::size_t
    {
      std::array<double, (2 * Rows) + 2> cuts{};
      std::size_t                        count = 0;
      cuts.at(count++)                         = stretch.Start;
      cuts.at(count++)                         = stretch.End;
      for (const Singular *singular : roots)
      {
        for (const Vector root : std::span{singular->Roots.Values}.first(singular->Roots.Count))
        {
          if (stretch.Start < root.X && root.X < stretch.End && Thinness(stretch, root) < RhoFloor)
          {
            cuts.at(count++) = root.X;
          }
        }
      }
      std::ranges::sort(std::span{cuts}.first(count));
      std::size_t pieces = 0;
      for (std::size_t index = count - 1; index > 0; --index)
      {
        if (cuts.at(index - 1) < cuts.at(index))
        {
          pending.at(pieces++) = {.Start = cuts.at(index - 1), .End = cuts.at(index)};
        }
      }
      return pieces;
    }

    /**
     *  A span's Bézier control points, homogeneous, each pole taken less an offset: the blossoms at
     *  the span's ends, as Kernel::Section takes them.
     *
     *  @param  spline The spline, valid.
     *  @param  span   The index of the knot that starts the span, the span not empty.
     *  @param  offset The point the poles are taken less.
     *  @return        The piece.
     */
    auto Lifted(const Spline &spline, std::size_t span, Point offset) -> Kernel::Piece
    {
      const Basis basis = Spanned(spline);
      const int   power = Leveling(spline, span);
      const auto  poles = [&spline, offset, power](std::size_t index) -> Homogeneous
      {
        const double weight =
            spline.Weights.empty() ? 1.0 : std::ldexp(spline.Weights[index], power);
        const Vector away = spline.Poles[index] - offset;
        return {.X = away.X * weight, .Y = away.Y * weight, .W = weight};
      };
      Kernel::Piece                 piece = {.Degree = basis.Degree, .Poles = {}};
      std::array<double, MaxDegree> values{};
      const std::span<double>       blossom = std::span{values}.first(basis.Degree);
      for (std::size_t index = 0; index <= basis.Degree; ++index)
      {
        std::ranges::fill(blossom, spline.Knots[span]);
        std::ranges::fill(blossom.last(index), spline.Knots[span + 1]);
        piece.Poles.at(index) = Blossom(basis, poles, {.Span = span, .Values = blossom});
      }
      return piece;
    }

    /**
     *  A polynomial span's hodograph's points, over its own parameter from 0 to 1, but for the
     *  degree they are to be multiplied by.
     *
     *  Its derivative is the spline of a degree less whose control points are d (P_(i+1) - P_i) /
     *  (u_(i+d+1) - u_(i+1)), on its knots but the first and the last; over a span of width w, its
     *  hodograph in the span's own parameter is w times that, whose Bézier points are its blossoms
     *  at the span's ends. Each control point is the difference of two poles times w over a width
     *  the span lies within: none is divided by nothing, and none of the differences is grown.
     *
     *  @param  spline The spline, valid.
     *  @param  span   The index of the knot that starts the span, the span not empty.
     *  @return        The points.
     */
    auto Derivated(const Spline &spline, std::size_t span) -> std::array<Vector, MaxDegree>
    {
      const auto                    degree = static_cast<std::size_t>(spline.Degree);
      const std::span<const double> knots  = spline.Knots;
      const Basis                   slopes = {.Degree = degree - 1,
                                              .Knots  = knots.subspan(1, knots.size() - 2),
                                              .Count  = spline.Poles.size() - 1};
      const double                  start  = knots[span];
      const double                  end    = knots[span + 1];
      const auto                    steps  = [&spline, knots, degree,
                                              width = end - start](std::size_t index) -> Homogeneous
      {
        const Vector step = (spline.Poles[index + 1] - spline.Poles[index])
                          * (width / (knots[index + degree + 1] - knots[index + 1]));
        return {.X = step.X, .Y = step.Y, .W = 1.0};
      };
      std::array<Vector, MaxDegree> points{};
      std::array<double, MaxDegree> values{};
      for (std::size_t point = 0; point < degree; ++point)
      {
        for (std::size_t index = 0; index + 1 < degree; ++index)
        {
          values.at(index) = index + 1 + point < degree ? start : end;
        }
        const Homogeneous blossom = Blossom(
            slopes, steps, {.Span = span - 1, .Values = std::span{values}.first(degree - 1)});
        points.at(point) = {.X = blossom.X, .Y = blossom.Y};
      }
      return points;
    }

    /**
     *  A piece's control points times the power of two that brings its heaviest weight to 1 or a
     *  little more, unless that would leave its lightest less than a double of all its digits, as
     *  Leveling brings a span's.
     *
     *  @param  piece The piece.
     *  @return       The piece, scaled.
     */
    auto Lightened(Kernel::Piece piece) -> Kernel::Piece
    {
      constexpr int Least    = -1021;
      double        heaviest = 0.0;
      double        lightest = std::numeric_limits<double>::max();
      for (const Homogeneous &pole : Kernel::Poles(piece))
      {
        heaviest = std::max(heaviest, pole.W);
        lightest = std::min(lightest, pole.W);
      }
      const int power = std::max(-std::ilogb(heaviest), Least - std::ilogb(lightest));
      for (Homogeneous &pole : std::span{piece.Poles}.first(piece.Degree + 1))
      {
        pole = {.X = std::ldexp(pole.X, power),
                .Y = std::ldexp(pole.Y, power),
                .W = std::ldexp(pole.W, power)};
      }
      return piece;
    }

    /**
     *  A piece of a span on its way to being cut into strips.
     */
    struct Stripping
    {
      Kernel::Piece   Shape = {}; /**< The piece, its heaviest weight near 1.   */
      Twofold::Double Start = {}; /**< How far past the span's start it starts. */
      Twofold::Double Rest  = {}; /**< How far short of the span's end it ends. */
      Twofold::Double Width = {}; /**< How much of the span it is.              */
      int             Depth = 0;  /**< How many cuts made it.                   */
    };

    /**
     *  Where to cut a piece whose weight has a root too near an end: Close from that end.
     *
     *  @param  poles The weight's roots.
     *  @return       The cut, from 0 to 1; 0 where no root is that near.
     */
    auto Crowded(const Singular &poles) -> double
    {
      double cut = 0.0;
      for (const Vector root : std::span{poles.Roots.Values}.first(poles.Roots.Count))
      {
        cut = Norm(root.X, root.Y) < Close ? Close : cut;
        cut = Norm(1.0 - root.X, root.Y) < Close ? 1.0 - Close : cut;
      }
      return cut;
    }

    /**
     *  A span's piece cut into strips, in order: cut Close from an end its weight has a root too
     *  near, and its parts again, until no root is, or until cut as often as doubles allow.
     *
     *  @param whole  The span's piece, its heaviest weight near 1.
     *  @param strip  What every strip of the span shares: its offset and its span.
     *  @param strips Where the strips go.
     */
    void Shredded(const Kernel::Piece &whole, const Strip &strip, std::vector<Strip> &strips)
    {
      std::vector<Stripping> pending = {Stripping{
          .Shape = whole, .Start = {}, .Rest = {}, .Width = Twofold::Lifted(1.0), .Depth = 0}};
      while (!pending.empty())
      {
        const Stripping top = pending.back();
        pending.pop_back();
        const Singular poles = WeightRoots(top.Shape);
        const double   cut   = Crowded(poles);
        if (Zero(cut) || top.Depth == MostCuts)
        {
          Strip even       = strip;
          even.Shape       = top.Shape;
          even.Poles       = poles;
          even.Where.Start = top.Start;
          even.Where.Rest  = top.Rest;
          even.Where.Width = top.Width.High + top.Width.Low;
          strips.push_back(even);
          continue;
        }
        const std::array<Kernel::Piece, 2> parts = Kernel::Split(top.Shape, cut);
        const Twofold::Double before = Twofold::Multiply(top.Width, Twofold::Lifted(cut));
        const Twofold::Double after  = Twofold::Difference(top.Width, before);
        pending.push_back({.Shape = Lightened(parts[1]),
                           .Start = Twofold::Add(top.Start, before),
                           .Rest  = top.Rest,
                           .Width = after,
                           .Depth = top.Depth + 1});
        pending.push_back({.Shape = Lightened(parts[0]),
                           .Start = top.Start,
                           .Rest  = Twofold::Add(top.Rest, after),
                           .Width = before,
                           .Depth = top.Depth + 1});
      }
    }

    /**
     *  A piece's coordinates brought near 1 by a power of two, its weights left as they are, so
     *  that no product of two of them overflows however far out the curve lies.
     *
     *  @param  piece The piece.
     *  @param  size  Where the power is put: the log2 of its largest coordinate, which every length
     *                of the piece is to be taken times 2 to.
     *  @return       The piece, its coordinates scaled.
     */
    auto Shrunken(Kernel::Piece piece, int &size) -> Kernel::Piece
    {
      double largest = Smallest;
      for (const Homogeneous &pole : Kernel::Poles(piece))
      {
        largest = std::max({largest, std::abs(pole.X), std::abs(pole.Y)});
      }
      size = std::ilogb(largest);
      for (Homogeneous &pole : std::span{piece.Poles}.first(piece.Degree + 1))
      {
        pole = {.X = std::ldexp(pole.X, -size), .Y = std::ldexp(pole.Y, -size), .W = pole.W};
      }
      return piece;
    }

    /**
     *  What a strip's area and moment rates are made of, about a point: X x X' for X its
     *  homogeneous coordinates less the point's times the weight, its coordinates less its offset,
     *  and its weight, all in Bernstein form, the first two each times a power of two of its own.
     */
    struct Fan
    {
      std::size_t                       Order   = 0;  /**< The strip's degree n.               */
      std::array<double, Rows>          Turning = {}; /**< X x X', of degree 2 n - 2.          */
      std::array<Vector, MaxDegree + 1> Places  = {}; /**< Its coordinates less its offset.    */
      std::array<double, MaxDegree + 1> Masses  = {}; /**< Its weight.                         */
      int                               Size    = 0;  /**< The log2 the coordinates are under. */
      int                               Turned  = 0;  /**< The log2 X x X' is under.           */
    };

    /**
     *  A coordinate of a strip's control point less a point's times its weight, in double-double: y
     *  + d w, the product and the sum exact but for the low part's own rounding.
     *
     *  @param  coordinate The control point's coordinate, less the strip's offset's times its
     *                     weight.
     *  @param  shift      The offset less the point, along the same axis.
     *  @param  weight     The control point's weight.
     *  @return            The coordinate less the point's times the weight.
     */
    auto Recentered(double coordinate, double shift, double weight) -> Twofold::Double
    {
      const Twofold::Two product = Twofold::TwoProduct(shift, weight);
      const Twofold::Two sum     = Twofold::TwoSum(coordinate, product.High);
      return {.High = sum.High, .Low = sum.Low + product.Low};
    }

    /**
     *  What a strip's area and moment rates are made of, about a point.
     *
     *  With Y the strip's coordinates less its offset times the weight and d the offset less the
     *  point, X = Y + d W, and X x X' is the pairs of X's control points, X_k x X_l, by the pairing
     *  identity Kernel::Derived works W X' - W' X out by. Each X_k is worked out in double-double,
     *  and each product and sum, and each coefficient rounded once: where the fan is thin, X x X'
     *  is far smaller than X and X', and the coefficients keep the digits the control points give
     *  it. Y and d are first brought near 1 by one power of two, the largest of them and of d times
     *  a weight, so that nothing overflows where the strip is far smaller than its distance from
     *  the point.
     *
     *  @param  strip The strip.
     *  @param  point The point.
     *  @return       The rates' parts.
     */
    auto Fanning(const Strip &strip, Point point) -> Fan
    {
      const std::size_t   order = strip.Shape.Degree;
      int                 size  = 0;
      const Kernel::Piece shape = Shrunken(strip.Shape, size);
      const Vector        shift = strip.Offset - point;
      double              most  = Smallest;
      for (const Homogeneous &pole : Kernel::Poles(strip.Shape))
      {
        most = std::max({most, std::abs(pole.X), std::abs(pole.Y), Largest(shift) * pole.W});
      }
      const int    turned = std::ilogb(most);
      const Vector moved  = {.X = std::ldexp(shift.X, -turned), .Y = std::ldexp(shift.Y, -turned)};
      Fan          fan    = {.Order   = order,
                             .Turning = {},
                             .Places  = {},
                             .Masses  = {},
                             .Size    = size,
                             .Turned  = 2 * turned};
      std::array<std::array<Twofold::Double, 2>, MaxDegree + 1> lifted{};
      for (std::size_t index = 0; index <= order; ++index)
      {
        const Homogeneous &pole = strip.Shape.Poles.at(index);
        fan.Places.at(index)    = {.X = shape.Poles.at(index).X, .Y = shape.Poles.at(index).Y};
        fan.Masses.at(index)    = pole.W;
        lifted.at(index)        = {Recentered(std::ldexp(pole.X, -turned), moved.X, pole.W),
                                   Recentered(std::ldexp(pole.Y, -turned), moved.Y, pole.W)};
      }
      const std::size_t                 twice  = (2 * order) - 2;
      const std::array<double, Rows>   &single = Binomials.at(order);
      std::array<Twofold::Double, Rows> sums{};
      for (std::size_t first = 0; first < order; ++first)
      {
        for (std::size_t second = first + 1; second <= order; ++second)
        {
          const auto &low             = lifted.at(first);
          const auto &high            = lifted.at(second);
          const auto  cross           = Twofold::Difference(Twofold::Multiply(low[0], high[1]),
                                                            Twofold::Multiply(low[1], high[0]));
          sums.at(first + second - 1) = Twofold::Add(
              sums.at(first + second - 1),
              Twofold::Multiply(cross, Twofold::Lifted(static_cast<double>(second - first)
                                                       * single.at(first)
                                                       * single.at(second))));
        }
      }
      for (std::size_t index = 0; index <= twice; ++index)
      {
        const Twofold::Double value =
            Twofold::Divide(sums.at(index), Twofold::Lifted(Binomials.at(twice).at(index)));
        fan.Turning.at(index) = value.High + value.Low;
      }
      return fan;
    }

    /**
     *  A strip's area and moment rates at a parameter: X x X' / W^2, and Y / W times that, Y its
     *  coordinates less its offset.
     *
     *  @param  fan  The rates' parts.
     *  @param  spot The parameter.
     *  @return      The area rate, and the moment rate along x and along y.
     */
    auto Rates(const Fan &fan, Abscissa spot) -> std::array<double, 3>
    {
      const std::size_t count = fan.Order + 1;
      const double      turning =
          Hornered(std::span<const double>{fan.Turning}.first((2 * fan.Order) - 1), spot);
      const double mass  = Hornered(std::span<const double>{fan.Masses}.first(count), spot);
      const Vector place = Hornered(std::span<const Vector>{fan.Places}.first(count), spot) / mass;
      const double area  = turning / (mass * mass) * Powered(spot, 2);
      return {area, place.X * area, place.Y * area};
    }

    /**
     *  The sizes of what a piece of a strip encloses: each rate at its length, by the rule of 2
     *  points.
     *
     *  @param  fan   The rates' parts.
     *  @param  piece The piece.
     *  @return       The area's size and the moment's.
     */
    auto Scale(const Fan &fan, Range piece) -> std::array<double, 2>
    {
      const std::array<Twofold::Double, 2> sizes = Ruled(
          [&fan](Abscissa spot) -> std::array<double, 2>
          {
            const std::array<double, 3> rates = Rates(fan, spot);
            return {std::abs(rates[0]), Length(Vector{.X = rates[1], .Y = rates[2]})};
          },
          piece, 1);
      return {sizes[0].High / 2.0, sizes[1].High / 3.0};
    }

    /**
     *  What a piece's area and moment are integrated to: the piece, how thin its ellipses are, and
     *  the area's share of the error and the moment's.
     */
    struct Shared
    {
      Range                 Piece  = {}; /**< The piece.                     */
      Girth                 Thick  = {}; /**< How thin its ellipses are.     */
      std::array<double, 2> Shares = {}; /**< The area's share and moment's. */
    };

    /**
     *  The fewest pairs of points whose bounds on a piece's area and moment rates keep each within
     *  its share.
     *
     *  On the piece's ellipse of rho, a polynomial of degree d in Bernstein form is at most its
     *  largest coefficient times (|z| + |1 - z|)^d, at most max(1, |2 c - 1|) + 2 R there, c the
     *  piece's middle and R the ellipse's half axis; the weight is at least as Minorant bounds it.
     *
     *  @param  fan    The rates' parts.
     *  @param  poles  The weight's roots.
     *  @param  shared The piece, how thin its ellipses are, and the shares.
     *  @return        The pairs.
     */
    auto Enough(const Fan &fan, const Singular &poles, const Shared &shared) -> std::size_t
    {
      const std::size_t count   = fan.Order + 1;
      const Range       piece   = shared.Piece;
      const double      half    = (piece.End - piece.Start) / 2.0;
      const double      middle  = piece.Start + half;
      double            turning = 0.0;
      double            place   = 0.0;
      for (const double value : std::span{fan.Turning}.first((2 * fan.Order) - 1))
      {
        turning = std::max(turning, std::abs(value));
      }
      for (const Vector value : std::span{fan.Places}.first(count))
      {
        place = std::max(place, Length(value));
      }
      auto                        pairs  = static_cast<double>(Legendre::Pairs);
      const std::array<double, 2> shares = {Trig::Log2(shared.Shares[0]),
                                            Trig::Log2(shared.Shares[1])};
      for (const double rho : Candidates(shared.Thick))
      {
        const double axis = half * (rho + (1.0 / rho)) / 2.0;
        const double spread =
            Trig::Log2(std::max(1.0, std::abs((2.0 * middle) - 1.0)) + (2.0 * axis));
        const double light = Minorant(poles, piece, axis);
        const double area  = Trig::Log2(turning * half)
                           + (static_cast<double>((2 * fan.Order) - 2) * spread)
                           - (2.0 * light);
        const double reach =
            Trig::Exp2(Trig::Log2(place) + (static_cast<double>(fan.Order) * spread) - light);
        pairs = std::min(pairs, std::max(Needed(area, rho, shares[0]),
                                         Needed(area + Trig::Log2(reach), rho, shares[1])));
      }
      return static_cast<std::size_t>(pairs);
    }
  } // namespace

  auto Whole(const Spline &spline, std::size_t span) -> Strip
  {
    const auto                    degree = static_cast<std::size_t>(spline.Degree);
    const std::span<const double> weights =
        spline.Weights.empty() ? std::span<const double>{}
                               : std::span{spline.Weights}.subspan(span - degree, degree + 1);
    const bool plain = std::ranges::all_of(weights, [&weights](double weight) -> bool
                                           { return Exactly(weight, weights.front()); });
    Point      offset =
        At(spline, spline.Knots[span] + ((spline.Knots[span + 1] - spline.Knots[span]) / 2.0));
    for (std::size_t index = span - degree; index <= span; ++index)
    {
      offset = std::isfinite(Length(spline.Poles[index] - offset)) ? offset : Point{};
    }
    Strip strip = {.Shape  = Lightened(Lifted(spline, span, offset)),
                   .Offset = offset,
                   .Where  = {.Span  = span,
                              .From  = spline.Knots[span],
                              .To    = spline.Knots[span + 1],
                              .Start = {},
                              .Rest  = {},
                              .Width = 1.0},
                   .Steps  = {},
                   .Poles  = {.Roots = {}, .Lead = 0.0},
                   .Plain  = plain};
    if (plain)
    {
      strip.Steps = Derivated(spline, span);
    }
    else
    {
      strip.Poles = WeightRoots(strip.Shape);
    }
    return strip;
  }

  auto Stripped(const Spline &spline) -> std::vector<Strip>
  {
    std::vector<Strip> strips;
    const auto         degree = static_cast<std::size_t>(spline.Degree);
    for (std::size_t span = degree; span < spline.Poles.size(); ++span)
    {
      if (!(spline.Knots[span] < spline.Knots[span + 1]))
      {
        continue;
      }
      const Strip strip = Whole(spline, span);
      if (strip.Plain)
      {
        strips.push_back(strip);
        continue;
      }
      Shredded(strip.Shape, strip, strips);
    }
    return strips;
  }

  auto Local(const Slice &where, double along) -> double
  {
    const double span = where.To - where.From;
    if (where.Start.High < 0.5)
    {
      const double ahead = (along - where.From) / span;
      return ((ahead - where.Start.High) - where.Start.Low) / where.Width;
    }
    const double behind = (where.To - along) / span;
    return 1.0 - (((behind - where.Rest.High) - where.Rest.Low) / where.Width);
  }

  auto Ahead(const Slice &where, double along) -> bool
  {
    const double span = where.To - where.From;
    if (where.Start.High < 0.5)
    {
      const Twofold::Double reach = Twofold::Multiply(where.Start, Twofold::Lifted(span));
      return (along - where.From) - reach.High <= reach.Low;
    }
    const Twofold::Double reach = Twofold::Multiply(
        Twofold::Difference(Twofold::Lifted(1.0), where.Start), Twofold::Lifted(span));
    return (where.To - along) - reach.High >= reach.Low;
  }

  auto Past(const Slice &where, double along) -> bool
  {
    const double span = where.To - where.From;
    if (where.Rest.High > 0.5)
    {
      const Twofold::Double reach = Twofold::Multiply(
          Twofold::Difference(Twofold::Lifted(1.0), where.Rest), Twofold::Lifted(span));
      return (along - where.From) - reach.High >= reach.Low;
    }
    const Twofold::Double reach = Twofold::Multiply(where.Rest, Twofold::Lifted(span));
    return (where.To - along) - reach.High <= reach.Low;
  }

  auto Paced(const Strip &strip) -> Pace
  {
    const std::size_t   order = strip.Shape.Degree;
    int                 size  = 0;
    const Kernel::Piece shape = Shrunken(strip.Shape, size);
    Pace                pace  = {.Degree   = 0,
                                 .Top      = {},
                                 .Weighing = 0,
                                 .Bottom   = {1.0},
                                 .Factor   = 1.0,
                                 .Order    = order,
                                 .Places   = {},
                                 .Masses   = {},
                                 .Size     = 0,
                                 .Turns    = {},
                                 .Poles    = strip.Poles};
    for (std::size_t index = 0; index <= order; ++index)
    {
      const Homogeneous &pole = shape.Poles.at(index);
      pace.Places.at(index)   = {.X = pole.X, .Y = pole.Y};
      pace.Masses.at(index)   = pole.W;
    }
    pace.Size = size;
    if (strip.Plain)
    {
      pace.Degree = order - 1;
      for (std::size_t index = 0; index < order; ++index)
      {
        const Vector step  = strip.Steps.at(index);
        pace.Top.at(index) = {.X = std::ldexp(step.X, -size), .Y = std::ldexp(step.Y, -size)};
      }
      pace.Factor = std::ldexp(static_cast<double>(order), size);
    }
    else
    {
      const Kernel::Hodograph speed = Kernel::Derived(shape);
      pace.Degree                   = speed.Degree;
      pace.Weighing                 = order;
      for (std::size_t index = 0; index <= speed.Degree; ++index)
      {
        pace.Top.at(index) = {.X = speed.Across.at(index), .Y = speed.Rising.at(index)};
      }
      std::ranges::copy(std::span{pace.Masses}.first(order + 1), pace.Bottom.begin());
      pace.Factor = std::ldexp(1.0, size);
    }
    const std::span<Vector> top     = std::span{pace.Top}.first(pace.Degree + 1);
    double                  largest = Smallest;
    for (const Vector value : top)
    {
      largest = std::max(largest, Largest(value));
    }
    const int power = -std::ilogb(largest);
    for (Vector &value : top)
    {
      value = {.X = std::ldexp(value.X, power), .Y = std::ldexp(value.Y, power)};
    }
    pace.Factor = std::ldexp(pace.Factor, -power);
    pace.Turns  = Singularities(top);
    return pace;
  }

  auto Measured(const Pace &pace, Range stretch, Effort effort) -> double
  {
    const double width = stretch.End - stretch.Start;
    if (!(pace.Turns.Lead > -Infinite) || !(width > 0.0))
    {
      return 0.0;
    }
    if (!std::isfinite(width))
    {
      return Infinite;
    }
    const double target = std::max(effort.Within.Error, 0x1p-55 * Inscribed(pace, stretch, 4));
    std::array<Range, StackDepth> pending{};
    std::size_t                   count  = Splayed({&pace.Turns, &pace.Poles}, stretch, pending);
    const std::size_t             budget = std::min(effort.Halvings, Halved);
    std::size_t                   spent  = 0;
    Twofold::Double               total  = {};
    while (count > 0)
    {
      const Range  piece = pending.at(--count);
      const double share = std::max(target * ((piece.End - piece.Start) / width),
                                    0x1p-55 * Inscribed(pace, piece, 1));
      const Girth  reach = Girthed(pace, piece, share);
      if (reach.Rho < RhoFloor && spent < budget)
      {
        pending.at(count++) = {.Start = reach.Split, .End = piece.End};
        pending.at(count++) = {.Start = piece.Start, .End = reach.Split};
        ++spent;
        continue;
      }
      total = Twofold::Add(total, Lengthened(pace, piece, {.Thick = reach, .Share = share}));
    }
    return total.High + total.Low;
  }

  auto Walked(const Pace &pace, double distance, double tolerance) -> double
  {
    constexpr int Steps  = 100;
    const Effort  effort = {.Within = {.Error = tolerance / 2.0}};
    const double  whole  = Measured(pace, {.Start = 0.0, .End = 1.0}, effort);
    double        low    = 0.0;
    double        high   = 1.0;
    double        along  = std::clamp(distance / std::max(whole, Smallest), 0.0, 1.0);
    for (int step = 0; step < Steps; ++step)
    {
      const double off = Measured(pace, {.Start = 0.0, .End = along}, effort) - distance;
      if (!(std::abs(off) > tolerance))
      {
        break;
      }
      low                = off < 0.0 ? along : low;
      high               = off < 0.0 ? high : along;
      const double speed = Speed(pace, {.Along = along, .Rest = 1.0 - along});
      const double next  = along - (off / std::max(speed, Smallest));
      along              = low < next && next < high ? next : low + ((high - low) / 2.0);
    }
    return along;
  }

  auto Enclosed(const Strip &strip, const Inquiry &asked) -> Kernel::Enclosure
  {
    const Fan                   fan    = Fanning(strip, asked.Origin);
    const std::array<double, 2> scale  = Scale(fan, {.Start = 0.0, .End = 1.0});
    const std::array<double, 2> target = {
        std::max(std::ldexp(asked.Error, -fan.Turned), 0x1p-55 * scale[0]),
        std::max(std::ldexp(asked.Error, -fan.Turned - fan.Size), 0x1p-55 * scale[1])};
    const Singular                none = {};
    std::array<Range, StackDepth> pending{};
    std::size_t count = Splayed({&none, &strip.Poles}, {.Start = 0.0, .End = 1.0}, pending);
    std::size_t spent = 0;
    std::array<Twofold::Double, 3> total{};
    while (count > 0)
    {
      const Range piece = pending.at(--count);
      const Girth reach = Thinnest(strip.Poles, piece);
      if (reach.Rho < RhoFloor && spent < std::min(asked.Halvings, Halved))
      {
        pending.at(count++) = {.Start = reach.Split, .End = piece.End};
        pending.at(count++) = {.Start = piece.Start, .End = reach.Split};
        ++spent;
        continue;
      }
      const double                width = piece.End - piece.Start;
      const std::array<double, 2> sized = Scale(fan, piece);
      const std::size_t           pairs = Enough(
          fan, strip.Poles,
          {.Piece  = piece,
           .Thick  = reach,
           .Shares = {std::max(target[0] * width, 0x1p-55 * sized[0]),
                      asked.Moments ? std::max(target[1] * width, 0x1p-55 * sized[1]) : Infinite}});
      const std::array<Twofold::Double, 3> sums =
          Ruled([&fan](Abscissa spot) -> std::array<double, 3> { return Rates(fan, spot); }, piece,
                pairs);
      for (std::size_t which = 0; which < 3; ++which)
      {
        total.at(which) = Twofold::Add(total.at(which), sums.at(which));
      }
    }
    const double area = std::ldexp((total[0].High + total[0].Low) / 2.0, fan.Turned);
    const Vector own  = {.X = std::ldexp(total[1].High + total[1].Low, fan.Turned + fan.Size),
                         .Y = std::ldexp(total[2].High + total[2].Low, fan.Turned + fan.Size)};
    return {.Area   = area,
            .Moment = (own / 3.0) + ((strip.Offset - asked.Origin) * (2.0 * area / 3.0))};
  }
} // namespace Geoxy
