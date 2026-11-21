/***************************************************************************************************
 *  @file      source/swift/strand.cpp
 *  @brief     A line, a ray or an endless line met with a spline of any degree, rational or not:
 *             the line's equation on the spline held as a spline of its own, by a coefficient per
 *             pole, each span's Bernstein coefficients its blossoms, its roots isolated in doubles
 *             and each polished by one double-double Newton step, the point and places proved the
 *             doubles nearest the true ones.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "real.hpp"
#include "straight.hpp"
#include "sure.hpp"
#include "swift/casteljau.hpp"
#include "swift/lane.hpp"
#include "swift/pieced.hpp"
#include "swift/roots.hpp"
#include "trig.hpp"
#include "twofold.hpp"

#include "geoxy/at.hpp"
#include "geoxy/box.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>

namespace Geoxy::Swift
{
  namespace
  {
    using Sure::Number;

    /**
     *  Room for a span's coefficients: a spline of degree fifteen at the most is met so; one of
     *  higher degree is left to the thorough way.
     */
    constexpr std::size_t Strands = 16;

    /**
     *  How far apart, in reaches, a span's poles and a line's ends may lie and the span still be
     *  met, as Strung keeps a section: the span lies in its poles' hull.
     */
    constexpr double Skirted = 128.0;

    /**
     *  The least a fraction of knot differences may be and not be nothing, as Blossomed holds it: a
     *  fraction below the normal doubles is rounded to a step of them, not to a share of itself.
     */
    constexpr double Wispiest = 0x1p-1000;

    /**
     *  A span of the spline as the kernel weighs it: where it lies and what bounds it.
     */
    struct Hank
    {
      const Spline *Shape  = nullptr; /**< The spline.                                  */
      std::size_t   Knot   = 0;       /**< The index of the knot that starts the span.  */
      std::size_t   Degree = 0;       /**< The spline's degree.                         */
      bool          Rated  = false;   /**< Whether its weights are weighed.             */
      double        Light  = 1.0;     /**< The least weight of the span's poles.        */
      double        Heavy  = 1.0;     /**< The most.                                    */
      Box           Hull   = {};      /**< The box of the span's poles, which holds it. */
    };

    /**
     *  A pole's weight, 1 where the spline's weights are not weighed.
     *
     *  @param  hank  The span.
     *  @param  index The pole.
     *  @return       Its weight.
     */
    auto Laden(const Hank &hank, std::size_t index) -> double
    {
      return hank.Rated ? hank.Shape->Weights[index] : 1.0;
    }

    /**
     *  The line's equation at a pole, times its weight: w (d x (P - S)), the differences exact and
     *  the products in double-double. A spline's point, as the sum of its poles times their weights
     *  and basis functions over its weight, gives the line's equation times its weight as the same
     *  sum of these: one spline function, of the spline's own knots and degree.
     *
     *  @param  lane  The straight curve.
     *  @param  hank  The span.
     *  @param  index The pole.
     *  @return       The coefficient, with its bound.
     */
    auto Weft(const Lane &lane, const Hank &hank, std::size_t index) -> Number
    {
      const Point  pole  = hank.Shape->Poles[index];
      const Number cross = (lane.Way[0] * Sure::Apart(pole.Y, lane.Start.Y))
                         - (lane.Way[1] * Sure::Apart(pole.X, lane.Start.X));
      return hank.Rated ? Sure::Exact(Laden(hank, index)) * cross : cross;
    }

    /**
     *  A span's Bernstein coefficient: the blossom of the span's coefficients at its start so many
     *  times, then at its end, by de Boor's rounds in doubles, each a convex mix of two cells by a
     *  fraction of knot differences (tests/proof/section.py: within the cells' error and ten steps
     *  of a double at their largest a round).
     *
     *  @param  hank  The span.
     *  @param  cells The span's coefficients, rounded, its degree and one.
     *  @param  ends  How many rounds take the span's end, the last ones.
     *  @return       The coefficient; nothing where a fraction falls below the normal doubles.
     */
    auto Plaited(const Hank &hank, std::array<double, Strands> cells, std::size_t ends)
        -> std::optional<double>
    {
      const std::size_t             degree = hank.Degree;
      const std::size_t             span   = hank.Knot;
      const std::span<const double> knots  = hank.Shape->Knots;
      for (std::size_t round = 1; round <= degree; ++round)
      {
        const double value = round + ends > degree ? knots[span + 1] : knots[span];
        for (std::size_t index = degree; index >= round; --index)
        {
          const double start    = knots[span - degree + index];
          const double end      = knots[span + 1 + index - round];
          const double ahead    = value - start;
          const double fraction = ahead / (end - start);
          if (ahead > 0.0 && !(fraction >= Wispiest))
          {
            return std::nullopt;
          }
          cells.at(index) = (cells.at(index - 1) * (1.0 - fraction)) + (cells.at(index) * fraction);
        }
      }
      return cells.at(degree);
    }

    /**
     *  De Boor's algorithm at a place of a span, in double-double, every rounding bounded: the
     *  place is u = a + w t on the span from a to b, w = b - a, and u less a knot k is w t + (a -
     *  k), each part exact or bounded, so that a place a sliver from the span's start keeps its
     *  digits. Each fraction is then (u - start) / (end - start), its rest 1 less it, and the
     *  place's own bound is carried through every mix.
     *
     *  @tparam CountT How many values a cell holds.
     *  @param  hank   The span.
     *  @param  cells  The cells: the span's poles' values, its degree and one.
     *  @param  along  The place, from 0 to 1, with its bound.
     *  @return        The values at the place.
     */
    template <std::size_t CountT>
    auto Warped(const Hank &hank, std::array<std::array<Number, CountT>, Strands> cells,
                const Number &along) -> std::array<Number, CountT>
    {
      const std::size_t             degree = hank.Degree;
      const std::size_t             span   = hank.Knot;
      const std::span<const double> knots  = hank.Shape->Knots;
      const double                  start  = knots[span];
      const Number                  ahead  = Sure::Apart(knots[span + 1], start) * along;
      for (std::size_t round = 1; round <= degree; ++round)
      {
        for (std::size_t index = degree; index >= round; --index)
        {
          const double from  = knots[span - degree + index];
          const double until = knots[span + 1 + index - round];
          const Number fraction =
              Sure::Quotient(ahead + Sure::Apart(start, from), Sure::Apart(until, from));
          const Number rest = Sure::Exact(1.0) - fraction;
          for (std::size_t axis = 0; axis < CountT; ++axis)
          {
            cells.at(index).at(axis) =
                (cells.at(index - 1).at(axis) * rest) + (cells.at(index).at(axis) * fraction);
          }
        }
      }
      return cells.at(degree);
    }

    /**
     *  The most the spline's derivative by the span's own place, from 0 to 1, may be on the span. A
     *  polynomial span's derivative is the spline of its derivative poles d (P_i+1 - P_i) /
     *  (u_i+d+1
     *  - u_i+1), times the span's width: no longer than the longest of them. A rational one's is
     *    (N' - p W') / W, N' - p W' = d sum (w_k+1 (Q_k+1 - p) - w_k (Q_k - p)) B_k for its Bézier
     *    points Q_k and weights w_k, each Q_k in its poles' hull and each w_k among their weights:
     *    no more than 2 d w_max times the hull's diagonal, over w_min. Each bound is grown past the
     *    roundings that work it out.
     *
     *  @param  hank The span.
     *  @return      The most |p'| may be.
     */
    auto Swiftest(const Hank &hank) -> double
    {
      const Spline &spline = *hank.Shape;
      const auto    order  = static_cast<double>(hank.Degree);
      const double  width  = spline.Knots[hank.Knot + 1] - spline.Knots[hank.Knot];
      if (hank.Rated)
      {
        const Box   &hull   = hank.Hull;
        const double across = Trig::Hypot(hull.Max.X - hull.Min.X, hull.Max.Y - hull.Min.Y);
        return 2.0 * order * hank.Heavy * across / hank.Light * (1.0 + 0x1p-48);
      }
      double most = 0.0;
      for (std::size_t index = hank.Knot - hank.Degree; index < hank.Knot; ++index)
      {
        const Vector step = spline.Poles[index + 1] - spline.Poles[index];
        const double over = spline.Knots[index + hank.Degree + 1] - spline.Knots[index + 1];
        most              = std::max(most, Length(step) / over);
      }
      return order * most * width * (1.0 + 0x1p-48);
    }

    /**
     *  The most the spline's derivative by its own parameter may be on the span, as Derivative
     *  works it out, and twice that: Pace's bound for a section, here from the span's poles.
     *
     *  @param  hank The span.
     *  @return      The pace.
     */
    auto Cadenced(const Hank &hank) -> double
    {
      const Spline &spline = *hank.Shape;
      const double  width  = spline.Knots[hank.Knot + 1] - spline.Knots[hank.Knot];
      return 4.0 * Swiftest(hank) / width;
    }

    /**
     *  The line's equation on a span, as the kernel weighs it: each pole's coefficient in
     *  double-double, for its true values; the span's Bernstein coefficients in doubles, within a
     *  bound of the true ones, and the most any of them is.
     */
    struct Loom
    {
      std::array<std::array<Number, 1>, Strands> Wefts   = {};  /**< Each pole's coefficient. */
      Stretch<Strands>                           Whole   = {};  /**< The span's coefficients. */
      double                                     Bound   = 0.0; /**< How far each may be off. */
      double                                     Largest = 0.0; /**< The most any of them is. */
      bool Far = false; /**< Whether the span is far from the line: its coefficients, which its
                             equation is a convex mix of, all past the reach one way. */
    };

    /**
     *  The line's equation on a span: its coefficients, rounded, each within its double-double's
     *  low part and bound of the exact one; their Bernstein coefficients on the span within that
     *  and 11 d u of their largest (tests/proof/section.py, as Spanned's slack for a section).
     *
     *  A span whose coefficients all lie past the reach one way, beyond their bound, is far: the
     *  equation on it is their mix by the basis functions, which are not below nothing and sum to
     *  1, and its Bernstein coefficients are not worked out.
     *
     *  @param  lane The straight curve.
     *  @param  hank The span.
     *  @param  cut  The reach the equation must pass, T.
     *  @return      The equation; nothing where a fraction of the knots falls below the normal
     *               doubles.
     */
    auto Loomed(const Lane &lane, const Hank &hank, double cut) -> std::optional<Loom>
    {
      const std::size_t           degree  = hank.Degree;
      const std::size_t           base    = hank.Knot - degree;
      Loom                        loom    = {};
      std::array<double, Strands> rounded = {};
      double                      error   = 0.0;
      double                      size    = 0.0;
      for (std::size_t index = 0; index <= degree; ++index)
      {
        const Number weft       = Weft(lane, hank, base + index);
        loom.Wefts.at(index)[0] = weft;
        rounded.at(index)       = weft.Value.High;
        error                   = std::max(error, std::abs(weft.Value.Low) + weft.Bound);
        size                    = std::max(size, std::abs(weft.Value.High));
      }
      const double clear = (cut + error) * Sure::Grow;
      const auto   span  = std::span{rounded}.first(degree + 1);
      loom.Far =
          std::ranges::all_of(span, [clear](double value) -> bool { return value > clear; })
          || std::ranges::all_of(span, [clear](double value) -> bool { return value < -clear; });
      if (loom.Far)
      {
        return loom;
      }
      for (std::size_t index = 0; index <= degree; ++index)
      {
        const std::optional<double> coefficient = Plaited(hank, rounded, index);
        if (!coefficient)
        {
          return std::nullopt;
        }
        loom.Whole.C.at(index) = *coefficient;
        loom.Largest           = std::max(loom.Largest, std::abs(*coefficient));
      }
      const auto order = static_cast<double>(degree);
      loom.Bound =
          ((error + (11.0 * order * Unit * (size + error))) * Sure::Grow * Sure::Grow) + Sure::Tiny;
      return loom;
    }

    /**
     *  What every root of a span's equation is proved with: the curves, the equation and its
     *  bounds, the spline's poles in homogeneous coordinates, and how fast the span may run.
     */
    struct Weave
    {
      const Geoxy::Swift::Lane                  *Lane      = nullptr; /**< The straight curve.   */
      const Hank                                *Span      = nullptr; /**< The span.             */
      const Loom                                *Equation  = nullptr; /**< Its equation.         */
      Geoxy::Tolerance                           Tolerance = {};      /**< The tolerance.        */
      std::array<std::array<Number, 3>, Strands> Poles     = {};      /**< w P and w, exactly.   */
      std::array<double, Strands>                Slopes    = {};      /**< f' 's coefficients.   */
      double                                     Cut       = 0.0; /**< T, the equation's reach. */
      double                                     Turn      = 0.0; /**< The most |f''| is.    */
      double                                     Slack     = 0.0; /**< f' 's error in doubles. */
      double                                     Speed     = 0.0; /**< The most |p'| is.     */
    };

    /**
     *  A span's equation's bounds for Newton's step and the poles it is proved on: f'' as the
     *  line's kernel bounds it, n (n - 1) times its coefficients' largest second difference with
     *  what their bound and rounding add; f' in doubles within twice n times each coefficient's
     *  bound, and the roundings of the differences and of de Casteljau's n - 1 rounds, three steps
     *  of their largest, 2 n largest, a round (tests/proof/strand.py).
     *
     *  @param weave The curves, the span and its equation, filled in.
     */
    void Woven(Weave &weave)
    {
      const Hank  &hank    = *weave.Span;
      const Loom  &loom    = *weave.Equation;
      const auto   degree  = hank.Degree;
      const auto   order   = static_cast<double>(degree);
      const auto  &whole   = loom.Whole.C;
      const double largest = loom.Largest;
      for (std::size_t index = 0; index + 2 <= degree; ++index)
      {
        weave.Turn =
            std::max(weave.Turn,
                     std::abs(whole.at(index + 2) - (2.0 * whole.at(index + 1)) + whole.at(index))
                         + (8.0 * loom.Bound)
                         + (8.0 * Unit * largest));
      }
      weave.Turn *= order * (order - 1.0) * (1.0 + 0x1p-48);
      for (std::size_t index = 0; index < degree; ++index)
      {
        weave.Slopes.at(index) = order * (whole.at(index + 1) - whole.at(index));
      }
      weave.Slack =
          (((2.0 * order * loom.Bound) + (((6.0 * order * order) + (4.0 * order)) * Unit * largest))
           * Sure::Grow)
          + std::numeric_limits<double>::min();
      weave.Speed          = Swiftest(hank);
      const Spline &spline = *hank.Shape;
      for (std::size_t index = 0; index <= degree; ++index)
      {
        const std::size_t which  = hank.Knot - degree + index;
        const Point       pole   = spline.Poles[which];
        const double      weight = Laden(hank, which);
        weave.Poles.at(index) =
            hank.Rated ? std::array{Sure::Product(weight, pole.X), Sure::Product(weight, pole.Y),
                                    Sure::Exact(weight)}
                       : std::array{Sure::Exact(pole.X), Sure::Exact(pole.Y), Sure::Exact(1.0)};
      }
    }

    /**
     *  The spline's point at a place of a span, the place's bound carried: a rational span's
     *  homogeneous point over its weight, a polynomial one's two coordinates alone.
     *
     *  @param  weave The span and its poles.
     *  @param  place The place, from 0 to 1, with its bound.
     *  @return       The point's coordinates.
     */
    auto Pointed(const Weave &weave, const Number &place) -> std::array<Number, 2>
    {
      const Hank &hank = *weave.Span;
      if (hank.Rated)
      {
        const std::array<Number, 3> lifted = Warped<3>(hank, weave.Poles, place);
        return {Sure::Quotient(lifted[0], lifted[2]), Sure::Quotient(lifted[1], lifted[2])};
      }
      std::array<std::array<Number, 2>, Strands> poles = {};
      for (std::size_t index = 0; index <= hank.Degree; ++index)
      {
        poles.at(index) = {weave.Poles.at(index)[0], weave.Poles.at(index)[1]};
      }
      return Warped<2>(hank, poles, place);
    }

    /**
     *  A root of a span's equation proved a crossing and placed: polished by one double-double
     *  Newton step, its point the spline's at the place, the place's bound carried, and its places
     *  the doubles nearest the true ones, as the kernel of a line and a Bézier curve proves them.
     *
     *  @param  weave What the root is proved with.
     *  @param  root  The root's stretch.
     *  @param  into  Where the crossing goes, where it lies on the straight curve.
     *  @return       Whether it was sure.
     */
    auto Tied(const Weave &weave, const Root &root, Spots &into) -> bool
    {
      const Lane      &lane      = *weave.Lane;
      const Hank      &hank      = *weave.Span;
      const Loom      &loom      = *weave.Equation;
      const Tolerance &tolerance = weave.Tolerance;
      if (!Steep(root.Slope / (lane.Speed * (1.0 + 0x1p-50) * weave.Speed * hank.Heavy), tolerance))
      {
        return false;
      }
      const double along = hank.Degree <= Highest
                             ? Guessed<Strands>(root, loom.Whole.C, hank.Degree)
                             : Guessed<Strands, 2 * Polishing>(root, loom.Whole.C, hank.Degree);
      const Number value = Warped<1>(hank, loom.Wefts, Sure::Exact(along))[0];
      const double slope = Steady(Valued(weave.Slopes, hank.Degree - 1, along));
      /**
       *  The mean value theorem, as the line's kernel takes it: the root t* lies within R =
       *  (|f(t0)| + E) / g of t0, g the least |f'| on the stretch, and the step t1 = t0 - f0 /
       *  f'(t0) is within [E + 2u |f0| + F (slack + |f''| R) / g] / |f'(t0)| of it, F the most
       *  |f(t0)| may be.
       */
      const double most  = (Sure::Size(value) + value.Bound) * Sure::Grow;
      const double far   = most / root.Slope * Sure::Grow;
      const double step  = (value.Value.High + value.Value.Low) / slope;
      const double off   = ((((value.Bound + (2.0 * Unit * Sure::Size(value)))
                              + (most * (weave.Slack + (weave.Turn * far)) / root.Slope))
                             / std::abs(slope))
                            + (2.0 * Unit * std::abs(step)))
                         * Sure::Grow;
      const Number place = Number{
          .Value = Twofold::Add(Twofold::Lifted(along), Twofold::Lifted(-step)), .Bound = off};
      /**
       *  The crossing's point: the spline's at the place, over its weight where it is rational.
       */
      const std::array<Number, 2> point = Pointed(weave, place);
      const Number                x     = point[0];
      const Number                y     = point[1];
      const Number                dx    = lane.Way[0];
      const Number                dy    = lane.Way[1];
      const Number line = Sure::Quotient(((x - Sure::Exact(lane.Start.X)) * dx)
                                             + ((y - Sure::Exact(lane.Start.Y)) * dy),
                                         (dx * dx) + (dy * dy));
      /**
       *  A line's end within the reach of the curve lies within R' = T / g times the speed, and the
       *  reach, of a crossing: the crossing then is that far from the line's end along it, or it is
       *  not sure.
       */
      const double near =
          (Reach(tolerance) + (weave.Speed * weave.Cut / root.Slope)) / (lane.Speed * Sure::Shrink);
      const double margin =
          near + (line.Bound * Sure::Grow) + (0x1p-50 * std::abs(line.Value.High));
      const Side where = Within(lane.Along, line.Value.High, margin);
      if (where != Side::In)
      {
        return where == Side::Out;
      }
      const Spline &spline  = *hank.Shape;
      const double  start   = spline.Knots[hank.Knot];
      const double  end     = spline.Knots[hank.Knot + 1];
      const Number  curve   = Sure::Exact(start) + (Sure::Apart(end, start) * place);
      const auto    settled = Sure::Rounded(std::array{x, y, Led(lane, line), curve});
      if (!settled)
      {
        return false;
      }
      const auto [pointX, pointY, lineT, curveT] = *settled;
      /**
       *  The thorough way holds the crossing's candidates to the tolerance alone where twice what a
       *  step of each place moves its curve is less, as Strung asks of a section.
       */
      const double moved    = 2.0 * lane.Speed * Above((2.0 * std::abs(lineT)) + 1.0);
      const double farthest = std::max(std::abs(start), std::abs(end));
      const double coarse   = 2.0 * (moved + (Cadenced(hank) * Above(2.0 * farthest)));
      if (!(2.0 * coarse <= tolerance.Length))
      {
        return false;
      }
      into.push_back({.At = {.X = pointX, .Y = pointY}, .Along = {lineT, curveT}, .Own = true});
      return true;
    }

    /**
     *  The roots of the line's equation on a span, each a crossing proved and placed, as the kernel
     *  of a line and a Bézier curve proves and places them, on the span's own numbers.
     *
     *  @param  weave The straight curve, the span and the tolerance; the rest is filled in.
     *  @param  into  Where the crossings go.
     *  @return       Whether it was sure.
     */
    auto Hanked(Weave &weave, Spots &into) -> bool
    {
      const Lane      &lane      = *weave.Lane;
      const Hank      &hank      = *weave.Span;
      const Tolerance &tolerance = weave.Tolerance;
      /**
       *  The line's equation is the curve's distance from the carrier times |d| and the weight W:
       *  where |f| passes T = reach |d| w_max, the curve is farther than the reach.
       */
      const double cut = Reach(tolerance) * lane.Speed * (1.0 + 0x1p-50) * hank.Heavy * Sure::Grow;
      const std::optional<Loom> loom = Loomed(lane, hank, cut);
      if (!loom || loom->Far)
      {
        return loom.has_value();
      }
      const std::optional<Roots<Strands>> roots =
          Isolated(loom->Whole, {loom->Bound, loom->Largest, cut}, hank.Degree);
      if (!roots || roots->Count == 0)
      {
        return roots.has_value();
      }
      const Spline &spline = *hank.Shape;
      const Point   from   = At(spline, spline.Knots[hank.Knot]);
      const Run     run    = {.Start = lane.Start, .Way = lane.Rough, .Along = lane.Along};
      const Run     chord  = {.Start = from,
                              .Way   = At(spline, spline.Knots[hank.Knot + 1]) - from,
                              .Along = {.Start = 0.0, .End = 1.0}};
      if (Dot(chord.Way, chord.Way) > 0.0 && !Unaligned(run, chord, tolerance))
      {
        return false;
      }
      weave.Equation = &*loom;
      weave.Cut      = cut;
      weave.Turn     = 0.0;
      Woven(weave);
      return std::ranges::all_of(std::span{roots->At}.first(roots->Count),
                                 [&weave, &into](const Root &root) -> bool
                                 { return Tied(weave, root, into); });
    }
  } // namespace

  auto Skeined(const Lane &lane, const Spline &spline, const Tolerance &tolerance)
      -> std::optional<Spots>
  {
    const auto degree = static_cast<std::size_t>(spline.Degree);
    if (degree < 2 || degree >= Strands || !Valid(spline))
    {
      return std::nullopt;
    }
    const bool rated = !spline.Weights.empty();
    const bool held =
        std::ranges::all_of(spline.Poles, [](Point pole) -> bool { return Held({pole.X, pole.Y}); })
        && std::ranges::all_of(spline.Knots, [](double knot) -> bool { return Held({knot}); })
        && std::ranges::all_of(spline.Weights, [](double weight) -> bool
                               { return weight > Narrowest && Held({weight}); });
    if (!held)
    {
      return std::nullopt;
    }
    const bool   bounded = lane.Count == 2;
    const Box    reach   = Grown(Grown(Box{}, lane.Ends[0]), lane.Ends[1]);
    const double clear   = Skirted * Reach(tolerance);
    Spots        spots;
    for (std::size_t span = degree; span < spline.Poles.size(); ++span)
    {
      if (!(spline.Knots[span] < spline.Knots[span + 1]))
      {
        continue;
      }
      Hank hank = {.Shape  = &spline,
                   .Knot   = span,
                   .Degree = degree,
                   .Rated  = rated,
                   .Light  = std::numeric_limits<double>::infinity(),
                   .Heavy  = 0.0,
                   .Hull   = {.Min = spline.Poles[span], .Max = spline.Poles[span]}};
      for (std::size_t index = span - degree; index <= span; ++index)
      {
        hank.Hull  = Grown(hank.Hull, spline.Poles[index]);
        hank.Light = std::min(hank.Light, Laden(hank, index));
        hank.Heavy = std::max(hank.Heavy, Laden(hank, index));
      }
      const Box   &hull = hank.Hull;
      const double gap  = std::max({hull.Min.X - reach.Max.X, reach.Min.X - hull.Max.X,
                                    hull.Min.Y - reach.Max.Y, reach.Min.Y - hull.Max.Y});
      if (bounded && gap > clear)
      {
        continue;
      }
      Weave weave = {.Lane = &lane, .Span = &hank, .Equation = nullptr, .Tolerance = tolerance};
      if (!Hanked(weave, spots))
      {
        return std::nullopt;
      }
    }
    return spots;
  }
} // namespace Geoxy::Swift
