/***************************************************************************************************
 *  @file      source/lengths/ellipse.cpp
 *  @brief     The length along an ellipse, by Carlson's elliptic integrals: each whole turn of it,
 *             and what is left.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "lengths/ellipse.hpp"

#include "degrees.hpp"
#include "lengths/pairs.hpp"
#include "radians.hpp"
#include "real.hpp"
#include "trig.hpp"
#include "turning.hpp"
#include "twofold.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>

namespace Geoxy::Lengths
{
  namespace
  {
    /**
     *  A number as a pair.
     *
     *  @param  value The number.
     *  @return       The pair.
     */
    auto PairOf(double value) -> Double
    {
      return {.High = value, .Low = 0.0};
    }

    /**
     *  A pair times a power of two, exactly.
     *
     *  @param  pair  The pair.
     *  @param  power The power of two.
     *  @return       The product.
     */
    auto ByPower(Double pair, double power) -> Double
    {
      return {.High = pair.High * power, .Low = pair.Low * power};
    }

    /**
     *  The sum of two pairs, left as a high part and what it leaves: within 2^-104 of the sum of
     *  their sizes, the low parts' sum being rounded once.
     *
     *  @param  lhs The first pair.
     *  @param  rhs The second.
     *  @return     The sum.
     */
    auto PairSum(Double lhs, Double rhs) -> Double
    {
      const Twofold::Two lead = Twofold::TwoSum(lhs.High, rhs.High);
      return {.High = lead.High, .Low = lead.Low + (lhs.Low + rhs.Low)};
    }

    /**
     *  The product of two pairs, left as a high part and what it leaves: within some 2^-103 of
     *  itself, the high parts' product exact as a pair and their cross terms rounded.
     *
     *  @param  lhs The first pair.
     *  @param  rhs The second.
     *  @return     The product.
     */
    auto PairProduct(Double lhs, Double rhs) -> Double
    {
      const Twofold::Two lead = Twofold::TwoProduct(lhs.High, rhs.High);
      return {.High = lead.High, .Low = lead.Low + ((lhs.High * rhs.Low) + (lhs.Low * rhs.High))};
    }

    /**
     *  The quotient of two pairs, by the high part's quotient and one correction through the
     *  remainder, which a fused multiply-add finds exactly: within some 2^-102 of itself.
     *
     *  @param  top    The dividend.
     *  @param  bottom The divisor, not nothing.
     *  @return        The quotient.
     */
    auto PairQuotient(Double top, Double bottom) -> Double
    {
      const double inverse = 1.0 / bottom.High;
      const double high    = top.High * inverse;
      const double rest    = std::fma(-high, bottom.High, top.High);
      return {.High = high, .Low = ((rest + top.Low) - (high * bottom.Low)) * inverse};
    }

    /**
     *  The square root of a pair, by one step of Newton's method from its high part's, the
     *  remainder found exactly by a fused multiply-add: within some 2^-103 of itself; nothing for
     *  nothing.
     *
     *  @param  value The pair, not below nothing.
     *  @return       The root.
     */
    auto PairRoot(Double value) -> Double
    {
      const double high = std::sqrt(value.High);
      const double rest = std::fma(-high, high, value.High);
      return {.High = high, .Low = (rest + value.Low) * (0.5 / std::max(high, 0x1p-1000))};
    }

    /**
     *  Carlson's symmetric elliptic integrals RF(x, y, 1) and RD(x, y, 1), as pairs.
     */
    struct Integrals
    {
      Double First  = {}; /**< RF(x, y, 1), of the first kind.  */
      Double Second = {}; /**< RD(x, y, 1), of the second kind. */
    };

    /**
     *  The seventh order series that ends RF's duplication, less its leading 1, in doubles: DLMF
     *  19.36.1 at deviations X and Y from the mean, Z = -X - Y.
     *
     *  @param  x The first deviation, 1 - x / A.
     *  @param  y The second.
     *  @return   The series less 1.
     */
    auto FirstSeries(double x, double y) -> double
    {
      const double z      = -x - y;
      const double second = (x * y) - (z * z);
      const double third  = x * y * z;
      return (-second / 10.0)
           + (third / 14.0)
           + (second * second / 24.0)
           - (3.0 * second * third / 44.0)
           - (5.0 * second * second * second / 208.0)
           + (3.0 * third * third / 104.0)
           + (second * second * third / 16.0);
    }

    /**
     *  The seventh order series that ends RD's duplication, less its leading 1, in doubles: DLMF
     *  19.36.2 at deviations X and Y from the mean, Z = -(X + Y) / 3.
     *
     *  @param  x The first deviation, 1 - x / A.
     *  @param  y The second.
     *  @return   The series less 1.
     */
    auto SecondSeries(double x, double y) -> double
    {
      const double z      = -(x + y) / 3.0;
      const double square = z * z;
      const double second = (x * y) - (6.0 * square);
      const double third  = ((3.0 * x * y) - (8.0 * square)) * z;
      const double fourth = 3.0 * ((x * y) - square) * square;
      const double fifth  = x * y * square * z;
      return (-3.0 * second / 14.0)
           + (third / 6.0)
           + (9.0 * second * second / 88.0)
           - (3.0 * fourth / 22.0)
           - (9.0 * second * third / 52.0)
           + (3.0 * fifth / 26.0)
           - (second * second * second / 16.0)
           + (3.0 * third * third / 40.0)
           + (3.0 * second * fourth / 20.0)
           + (45.0 * second * second * third / 272.0)
           - (9.0 * ((third * fourth) + (second * fifth)) / 68.0);
    }

    /**
     *  How far a pair lies from a mean, as a share of the mean, in a double.
     *
     *  @param  mean  The mean.
     *  @param  value The pair.
     *  @return       1 - value / mean.
     */
    auto Deviation(Double mean, Double value) -> double
    {
      const Double apart = PairSum(mean, Twofold::Negated(value));
      return (apart.High + apart.Low) / mean.High;
    }

    /**
     *  A series's value, 1 and what it adds, as a pair.
     *
     *  @param  rest What it adds to 1.
     *  @return      The pair.
     */
    auto OnePlus(double rest) -> Double
    {
      const Twofold::Two sum = Twofold::FastTwoSum(1.0, rest);
      return {.High = sum.High, .Low = sum.Low};
    }

    /**
     *  RF(x, y, 1) and RD(x, y, 1) in double-double, by one run of Carlson's duplication for both.
     *
     *  The two take the same arguments, and each step moves them alike a quarter of the way toward
     *  their means, x_(n+1) = (x_n + lambda_n) / 4 for lambda_n = sqrt(x_n y_n) + sqrt(y_n z_n) +
     *  sqrt(z_n x_n); RD sums 4^-n / (sqrt z_n (z_n + lambda_n)) besides. A mean A_n of either is
     *  A_0's moved alike, so every argument lies within 4^-n |A_0 - x_0| of it; once that is 2^-9
     *  of the mean or less, the seventh order series of DLMF 19.36.1 and 19.36.2 ends each, off by
     *  at most 0.016 and 0.08 times the eighth power of the deviation (as mpmath finds it), so by
     *  2^-75 of itself at most. Its terms past the leading 1 are below 2^-15 and are worked out in
     *  doubles, to 2^-68 of themselves; the rest is in pairs, each step to about 2^-102.
     *
     *  @pre      x and y are not negative, and not both zero.
     *  @param  x The first argument, as a pair.
     *  @param  y The second.
     *  @return   RF(x, y, 1) and RD(x, y, 1).
     */
    auto Duplicated(Double x, Double y) -> Integrals
    {
      const double first   = (x.High + y.High + 1.0) / 3.0;
      const double second  = (x.High + y.High + 3.0) / 5.0;
      const double reach   = 0x1p9
                           * std::max({std::abs(first - x.High), std::abs(first - y.High),
                                       std::abs(first - 1.0), std::abs(second - x.High),
                                       std::abs(second - y.High), std::abs(second - 1.0)});
      Double       p       = x;
      Double       q       = y;
      Double       r       = PairOf(1.0);
      Double       sum     = {};
      double       meanOne = first;
      double       meanTwo = second;
      double       scale   = 1.0;
      while (reach * scale >= std::min(meanOne, meanTwo))
      {
        const Double rootP = PairRoot(p);
        const Double rootQ = PairRoot(q);
        const Double rootR = PairRoot(r);
        const Double step  = PairSum(PairSum(PairProduct(rootP, rootQ), PairProduct(rootQ, rootR)),
                                     PairProduct(rootR, rootP));
        sum      = PairSum(sum, PairQuotient(PairOf(scale), PairProduct(rootR, PairSum(r, step))));
        p        = ByPower(PairSum(p, step), 0.25);
        q        = ByPower(PairSum(q, step), 0.25);
        r        = ByPower(PairSum(r, step), 0.25);
        meanOne  = (meanOne + step.High) / 4.0;
        meanTwo  = (meanTwo + step.High) / 4.0;
        scale   /= 4.0;
      }
      const Double both = PairSum(p, q);
      const Double one  = PairQuotient(PairSum(both, r), PairOf(3.0));
      const Double two  = PairQuotient(PairSum(both, PairProduct(r, PairOf(3.0))), PairOf(5.0));
      return {.First  = PairQuotient(OnePlus(FirstSeries(Deviation(one, p), Deviation(one, q))),
                                     PairRoot(one)),
              .Second = PairSum(
                  PairProduct(sum, PairOf(3.0)),
                  PairQuotient(
                      ByPower(OnePlus(SecondSeries(Deviation(two, p), Deviation(two, q))), scale),
                      PairProduct(two, PairRoot(two))))};
    }

    /**
     *  An ellipse's shape as its length's integral takes it: the square of its speed at its unit
     *  circle's angle t is a² (1 - m sin²(t - phi)), a its longest semi-axis and phi the angle of
     *  the unit circle that goes there, m = 1 - b² / a², b its shortest.
     *
     *  Everything is divided by a power of two near the map's size, and worked out in pairs from
     *  the map's products, each exact as a pair: b² as det² / a², so that 1 - m keeps its precision
     *  however flat the ellipse.
     */
    struct Elongation
    {
      Double Major  = {};  /**< a², divided twice.                      */
      Double Flat   = {};  /**< m.                                      */
      Double Third  = {};  /**< m / 3.                                  */
      Double Rest   = {};  /**< 1 - m, b² / a².                         */
      Double Cosine = {};  /**< cos phi.                                */
      Double Sine   = {};  /**< sin phi.                                */
      double Phase  = 0.0; /**< phi in degrees, near.                   */
      double Scale  = 1.0; /**< The power of two the map is divided by. */
    };

    /**
     *  The sum of two products of doubles, as a pair: each product exact as a pair, the sum to
     *  about 2^-104 of the products' sizes.
     *
     *  @param  first  The first two doubles.
     *  @param  second The second two.
     *  @return        The sum.
     */
    auto Products(Vector first, Vector second) -> Double
    {
      const Twofold::Two lead  = Twofold::TwoProduct(first.X, first.Y);
      const Twofold::Two trail = Twofold::TwoProduct(second.X, second.Y);
      return Twofold::Add({.High = lead.High, .Low = lead.Low},
                          {.High = trail.High, .Low = trail.Low});
    }

    /**
     *  An ellipse's shape for its length.
     *
     *  With x and y the map's axes, the square of its speed is (x·x + y·y) / 2 + ((y·y - x·x) / 2)
     *  cos 2t - (x·y) sin 2t, so a² is the mean of x·x and y·y plus the length of (c, s) = ((y·y
     *  - x·x) / 2, -x·y), and 2 phi is the angle of (c, s): phi's cosine and sine come by the
     *    half-angle formulas, the larger of the two by a root, the other by a quotient. A circle
     *    has no such angle, and takes phi = 0.
     *
     *  @param  map The ellipse's map.
     *  @return     The shape; a² alone, not above 0, for a map of nothing or one not finite.
     */
    auto Elongated(const Matrix &map) -> Elongation
    {
      const double most =
          std::max({std::abs(map.X.X), std::abs(map.X.Y), std::abs(map.Y.X), std::abs(map.Y.Y)});
      int power = 0;
      std::frexp(most, &power);
      const double scale  = std::ldexp(1.0, power);
      const Vector x      = map.X / scale;
      const Vector y      = map.Y / scale;
      const Double across = Products({.X = x.X, .Y = x.X}, {.X = x.Y, .Y = x.Y});
      const Double upward = Products({.X = y.X, .Y = y.X}, {.X = y.Y, .Y = y.Y});
      const Double cross  = Products({.X = x.X, .Y = y.X}, {.X = x.Y, .Y = y.Y});
      const Double area   = Products({.X = x.X, .Y = y.Y}, {.X = -x.Y, .Y = y.X});
      const Double cosine = ByPower(Twofold::Difference(upward, across), 0.5);
      const Double sine   = Twofold::Negated(cross);
      const Double swing =
          PairRoot(Twofold::Add(PairProduct(cosine, cosine), PairProduct(sine, sine)));
      const Double major = Twofold::Add(ByPower(Twofold::Add(across, upward), 0.5), swing);
      if (!(major.High > 0.0))
      {
        return {.Major = major};
      }
      const Double flat         = PairQuotient(ByPower(swing, 2.0), major);
      Double       turnedCosine = PairOf(1.0);
      Double       turnedSine   = {};
      if (swing.High > 0.0)
      {
        const Double twice = ByPower(swing, 2.0);
        if (cosine.High >= 0.0)
        {
          turnedCosine = PairRoot(PairQuotient(Twofold::Add(swing, cosine), twice));
          turnedSine   = PairQuotient(sine, PairProduct(twice, turnedCosine));
        }
        else
        {
          const Double root = PairRoot(PairQuotient(Twofold::Difference(swing, cosine), twice));
          turnedSine        = std::signbit(sine.High) ? Twofold::Negated(root) : root;
          turnedCosine      = PairQuotient(sine, PairProduct(twice, turnedSine));
        }
      }
      return {.Major  = major,
              .Flat   = flat,
              .Third  = PairQuotient(flat, PairOf(3.0)),
              .Rest   = PairQuotient(PairProduct(area, area), PairProduct(major, major)),
              .Cosine = turnedCosine,
              .Sine   = turnedSine,
              .Phase  = Trig::Atan2(turnedSine.High, turnedCosine.High) * (180.0 / Pi),
              .Scale  = scale};
    }

    /**
     *  The cosine and sine of an angle, as pairs.
     */
    struct Cosines
    {
      Double Cosine = {}; /**< The cosine. */
      Double Sine   = {}; /**< The sine.   */
    };

    /**
     *  The cosine and sine of degrees held as a pair.
     *
     *  @param  degrees The degrees, finite.
     *  @return         Their cosine and sine, within 2^-68 of themselves, as Turning::CosinePair
     *                  gives them.
     */
    auto Spun(Double degrees) -> Cosines
    {
      const std::array<Twofold::Two, 2> way =
          Turning::CosinePair({.High = degrees.High, .Low = degrees.Low});
      return {.Cosine = {.High = way[0].High, .Low = way[0].Low},
              .Sine   = {.High = way[1].High, .Low = way[1].Low}};
    }

    /**
     *  An angle turned by another: the cosine and sine of their sum.
     *
     *  @param  way The angle.
     *  @param  by  The angle it is turned by.
     *  @return     The sum's cosine and sine.
     */
    auto Rotated(const Cosines &way, const Cosines &by) -> Cosines
    {
      return {.Cosine = PairSum(PairProduct(way.Cosine, by.Cosine),
                                Twofold::Negated(PairProduct(way.Sine, by.Sine))),
              .Sine = PairSum(PairProduct(way.Sine, by.Cosine), PairProduct(way.Cosine, by.Sine))};
    }

    /**
     *  An angle of the integral of the second kind as whole half turns k and the sine and cosine of
     *  what is left, v = u - k pi, within a quarter turn of nothing: E(u) = 2 k E(m) + E(v).
     */
    struct Remnant
    {
      double Turns  = 0.0; /**< k.                        */
      Double Sine   = {};  /**< sin v.                    */
      Double Cosine = {};  /**< cos v, not below nothing. */
    };

    /**
     *  An angle of the integral of the second kind, its half turns k first u / pi rounded from the
     *  doubles, then moved by one where the cosine of what is left, worked out in pairs, is below
     *  nothing: where u is so near a half turn and a quarter that the doubles cannot tell which way
     *  it lies.
     *
     *  @param  turns u / pi, near.
     *  @param  way   The cosine and sine of u.
     *  @return       Its half turns and what is left.
     */
    auto HalfTurned(double turns, const Cosines &way) -> Remnant
    {
      double half   = std::rint(turns);
      Double sine   = way.Sine;
      Double cosine = way.Cosine;
      if (!Zero(std::fmod(half, 2.0)))
      {
        sine   = Twofold::Negated(sine);
        cosine = Twofold::Negated(cosine);
      }
      if (cosine.High < 0.0)
      {
        half   += sine.High > 0.0 ? 1.0 : -1.0;
        sine    = Twofold::Negated(sine);
        cosine  = Twofold::Negated(cosine);
      }
      return {.Turns = half, .Sine = sine, .Cosine = cosine};
    }

    /**
     *  The incomplete integral of the second kind within a quarter turn of nothing, E(v) = sin v
     *  RF(cos² v, 1 - m sin² v, 1) - (m / 3) sin³ v RD(cos² v, 1 - m sin² v, 1), 1 - m sin² v taken
     *  as cos² v + (1 - m) sin² v, so that it keeps its precision near the end of a flat ellipse's
     *  major axis.
     *
     *  @param  oval The ellipse's shape, of some area.
     *  @param  left The sine and cosine of v.
     *  @return      E(v).
     */
    auto Incomplete(const Elongation &oval, const Remnant &left) -> Double
    {
      const Double    square = PairProduct(left.Cosine, left.Cosine);
      const Double    sines  = PairProduct(left.Sine, left.Sine);
      const Integrals both   = Duplicated(square, PairSum(square, PairProduct(oval.Rest, sines)));
      return PairSum(PairProduct(left.Sine, both.First),
                     Twofold::Negated(PairProduct(
                         PairProduct(oval.Third, PairProduct(left.Sine, sines)), both.Second)));
    }

    /**
     *  The complete integral of the second kind, by Gauss's arithmetic-geometric mean: from a_0 =
     *  1, g_0 = sqrt(1 - m) and c_0 = sqrt(m), a_(n+1) = (a_n + g_n) / 2, g_(n+1) = sqrt(a_n g_n)
     *  and c_(n+1) = c_n² / (4 a_(n+1)), and E(m) = pi / (2 a_n) (1 - sum 2^(n-1) c_n²), DLMF
     *  19.8.1 and 19.8.6, in pairs. Each c is about the square of the last over four, so once one
     *  is below 2^-40 the terms left are below 2^-80 of the sum's 1, and a_n is the mean to as
     *  near. A map of no area has m = 1 and E(m) = 1.
     *
     *  @param  oval The ellipse's shape.
     *  @return      E(m).
     */
    auto Complete(const Elongation &oval) -> Double
    {
      if (!(oval.Rest.High > 0.0))
      {
        return PairOf(1.0);
      }
      Double mean  = PairOf(1.0);
      Double root  = PairRoot(oval.Rest);
      Double gap   = PairRoot(oval.Flat);
      Double sum   = ByPower(oval.Flat, 0.5);
      double power = 1.0;
      while (gap.High >= 0x1p-40)
      {
        const Double next  = ByPower(PairSum(mean, root), 0.5);
        root               = PairRoot(PairProduct(mean, root));
        gap                = PairQuotient(PairProduct(gap, gap), ByPower(next, 4.0));
        mean               = next;
        sum                = PairSum(sum, ByPower(PairProduct(gap, gap), power));
        power             *= 2.0;
      }
      return PairQuotient(PairProduct({.High = Table::HalfPi.High, .Low = Table::HalfPi.Low},
                                      PairSum(PairOf(1.0), Twofold::Negated(sum))),
                          mean);
    }

    /**
     *  What E(v1) - E(v0) comes to for v0 and v1 within a quarter turn of nothing, by the addition
     *  theorem (DLMF 19.11.1 and 19.11.2): the angle w at which E(w) - m sin v1 sin v0 sin w is it,
     *  and that last term.
     */
    struct Junction
    {
      Remnant Left   = {}; /**< w, as half turns and what is left. */
      Double  Excess = {}; /**< m sin v1 sin v0 sin w.             */
    };

    /**
     *  The angle of the addition theorem: sin w = (sin v1 cos v0 D0 - sin v0 cos v1 D1) / g and cos
     *  w = (cos v1 cos v0 + sin v1 sin v0 D1 D0) / g, for D = sqrt(1 - m sin² v) and g = 1 - m sin²
     *  v0 sin² v1. w lies within a half turn of nothing, F(w) being F(v1) - F(v0); so E(v1) - E(v0)
     *  takes one integral, not two, and a short stretch is no difference of two near ones.
     *
     *  v1 is v0 turned by the stretch's sweep d, and a half turn where their half turns differ by
     *  an odd number: so sin(v1 - v0) is sin d, and sin² v1 - sin² v0 is sin d (sin d cos 2 v0 +
     *  cos d sin 2 v0). The numerator of sin w is D0 sin(v1 - v0) + m sin v0 cos v1 (sin² v1 - sin²
     *  v0) / (D0 + D1), which leaves nothing to cancel however short the stretch, and g is cos² v0
     *  + sin² v0 cos² v1 + (1 - m) sin² v0 sin² v1, which keeps its precision where both ends lie
     *  near the end of a flat ellipse's major axis.
     *
     *  @param  oval The ellipse's shape, of some area.
     *  @param  ends Where the stretch starts and where it ends.
     *  @param  step The cosine and sine of d, a half turn on where the ends' half turns differ by
     *               an odd number.
     *  @return      The angle and the excess.
     */
    auto Bridged(const Elongation &oval, const std::array<Remnant, 2> &ends, const Cosines &step)
        -> Junction
    {
      const Remnant &first       = ends[0];
      const Remnant &last        = ends[1];
      const Double   startCosine = PairProduct(first.Cosine, first.Cosine);
      const Double   startSine   = PairProduct(first.Sine, first.Sine);
      const Double   endCosine   = PairProduct(last.Cosine, last.Cosine);
      const Double   endSine     = PairProduct(last.Sine, last.Sine);
      const Double   startRoot = PairRoot(PairSum(startCosine, PairProduct(oval.Rest, startSine)));
      const Double   endRoot   = PairRoot(PairSum(endCosine, PairProduct(oval.Rest, endSine)));
      const Double   both      = PairProduct(first.Sine, last.Sine);
      const Double   gap       = PairSum(PairSum(startCosine, PairProduct(startSine, endCosine)),
                                         PairProduct(oval.Rest, PairProduct(startSine, endSine)));
      const Double   squares   = PairProduct(
          step.Sine,
          PairSum(PairProduct(step.Sine, PairSum(startCosine, Twofold::Negated(startSine))),
                  PairProduct(step.Cosine, ByPower(PairProduct(first.Sine, first.Cosine), 2.0))));
      const Double rise = PairSum(
          PairProduct(startRoot, step.Sine),
          PairQuotient(
              PairProduct(PairProduct(oval.Flat, PairProduct(first.Sine, last.Cosine)), squares),
              PairSum(startRoot, endRoot)));
      const Double sine = PairQuotient(rise, gap);
      const Double cosine =
          PairQuotient(PairSum(PairProduct(last.Cosine, first.Cosine),
                               PairProduct(both, PairProduct(endRoot, startRoot))),
                       gap);
      return {.Left   = HalfTurned(0.0, {.Cosine = cosine, .Sine = sine}),
              .Excess = PairProduct(oval.Flat, PairProduct(both, sine))};
    }

    /**
     *  An angle as degrees in a pair, whole turns taken off where they do not count.
     *
     *  @param  angle The angle, finite.
     *  @return       Held in degrees, its degrees; held in radians, its degrees within 2^-104,
     *                whole turns taken off, as Turning::Reduced gives them.
     */
    auto Started(Angle angle) -> Double
    {
      if (angle.HeldInDegrees())
      {
        return PairOf(angle.InDegrees());
      }
      const Twofold::Two degrees = Turning::Reduced(angle.InRadians());
      return {.High = degrees.High, .Low = degrees.Low};
    }

    /**
     *  A share of a sweep in degrees, as a pair, every whole turn of it kept.
     *
     *  @param  sweep The sweep.
     *  @param  share The share, as a pair.
     *  @return       The sweep times the share, in degrees, within about 2^-104 of itself.
     */
    auto SweptDegrees(Angle sweep, Double share) -> Double
    {
      if (sweep.HeldInDegrees())
      {
        return PairProduct(PairOf(sweep.InDegrees()), share);
      }
      return PairProduct(PairProduct(PairOf(sweep.InRadians()), share),
                         {.High = Table::ToDegrees.High, .Low = Table::ToDegrees.Low});
    }
  } // namespace

  auto Arclength(const Ellipse &ellipse, Range stretch) -> double
  {
    const auto held = [](Angle angle) -> double
    {
      return angle.HeldInDegrees() ? angle.InDegrees() : angle.InRadians();
    };
    const Elongation oval = Elongated(ellipse.Map);
    if (!(oval.Major.High > 0.0))
    {
      return 0.0;
    }
    if (!std::isfinite(held(ellipse.Start) + held(ellipse.Sweep)))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const Twofold::Two share = Twofold::TwoSum(stretch.End, -stretch.Start);
    const Double       from =
        Twofold::Add(Started(ellipse.Start), SweptDegrees(ellipse.Sweep, PairOf(stretch.Start)));
    const Double span = SweptDegrees(ellipse.Sweep, {.High = share.High, .Low = share.Low});
    if (!std::isfinite(from.High + span.High))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const Cosines begun =
        Rotated(Spun(from), {.Cosine = oval.Cosine, .Sine = Twofold::Negated(oval.Sine)});
    const Cosines swing = Spun(span);
    const Remnant first = HalfTurned((from.High - oval.Phase) / 180.0, begun);
    const Remnant last =
        HalfTurned((from.High + span.High - oval.Phase) / 180.0, Rotated(begun, swing));
    double        turns = last.Turns - first.Turns;
    const Cosines step  = Zero(std::fmod(turns, 2.0))
                            ? swing
                            : Cosines{.Cosine = Twofold::Negated(swing.Cosine),
                                      .Sine   = Twofold::Negated(swing.Sine)};
    Double        apart = PairSum(PairProduct(first.Sine, PairSum(step.Cosine, PairOf(-1.0))),
                                  PairProduct(first.Cosine, step.Sine));
    if (oval.Rest.High > 0.0)
    {
      const Junction bridge  = Bridged(oval, {first, last}, step);
      turns                 += bridge.Left.Turns;
      apart = PairSum(Incomplete(oval, bridge.Left), Twofold::Negated(bridge.Excess));
    }
    if (!Zero(turns))
    {
      apart = PairSum(apart, PairProduct(Complete(oval), PairOf(2.0 * turns)));
    }
    const Double length = PairProduct(PairRoot(oval.Major), apart);
    return std::abs(length.High + length.Low) * oval.Scale;
  }
} // namespace Geoxy::Lengths
