/***************************************************************************************************
 *  @file      source/swift/braid.cpp
 *  @brief     A cubic Bézier curve met with a curve: its implicit equation, a Bézout determinant,
 *             along each piece of the other.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/braid.hpp"

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

#include "geoxy/bezier/cubic.hpp"
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
     *  How high a spline's span's degree may be for a cubic Bézier curve's equation on it, of three
     *  times its degree: five, so fifteen, within the sixteen the isolation's running bound and the
     *  slope's are proved to (tests/proof/strand.py) and an ample equation holds.
     */
    constexpr std::size_t Tightest = 5;

    /**
     *  How large an entry of a cubic's matrix may be on a span, so that no product of three and no
     *  sum of such products, nor what bounds them, runs past the doubles: 2^300.
     */
    constexpr double Thickest = 0x1p300;

    /**
     *  A cubic Bézier curve's implicit cubic. Each side is the exact difference of two doubles, a
     *  side's multiple and each A_ij in double-double, each S and K rounded and its slip kept. A
     *  step of its place moves it by no more than the thorough way's Strict charges a Bézier
     *  curve's, four times three times its longest control step's coordinate times a step of the
     *  doubles at 2.
     *
     *  @param  cubic The curve.
     *  @return       Its cubic; nothing where a number lies outside the range the kernels take.
     */
    auto Braiding(const Bezier::Cubic &cubic) -> std::optional<Braid>
    {
      const std::array<Point, 4> poles = {cubic.Start, cubic.Out, cubic.In, cubic.End};
      if (!std::ranges::all_of(poles, [](Point pole) -> bool { return Held({pole.X, pole.Y}); }))
      {
        return std::nullopt;
      }
      using Pair     = std::array<Number, 2>;
      const auto gap = [&poles](std::size_t head, std::size_t tail) -> Pair
      {
        return {Sure::Apart(poles.at(head).X, poles.at(tail).X),
                Sure::Apart(poles.at(head).Y, poles.at(tail).Y)};
      };
      const auto cross = [](const Pair &one, const Pair &other) -> Number
      {
        return (one[0] * other[1]) - (one[1] * other[0]);
      };
      const auto times = [](double factor, const Pair &side) -> Pair
      {
        return {Sure::Exact(factor) * side[0], Sure::Exact(factor) * side[1]};
      };
      const Pair   first  = gap(1, 0);
      const Pair   second = gap(2, 0);
      const Pair   third  = gap(3, 0);
      const Pair   inner  = times(9.0, gap(2, 1));
      const Number outer  = cross(third, second);
      const Number across = cross(third, first);
      const Number middle = cross(second, first);
      Braid braid = {.Hub   = cubic.Start,
                     .Sides = {times(3.0, gap(3, 2)), times(3.0, gap(3, 1)), third,
                               Pair{third[0] + inner[0], third[1] + inner[1]}, times(3.0, second),
                               times(3.0, first)},
                     .Areas = {Sure::Exact(3.0) * outer, Sure::Exact(3.0) * across,
                               Sure::Exact(0.0), Sure::Exact(9.0) * middle, Sure::Exact(0.0),
                               Sure::Exact(0.0)}};
      for (std::size_t which = 0; which < Fibres; ++which)
      {
        const Pair   &side     = braid.Sides.at(which);
        const Number &area     = braid.Areas.at(which);
        const Vector  way      = {.X = side[0].Value.High, .Y = side[1].Value.High};
        braid.Ways.at(which)   = way;
        braid.Slips.at(which)  = std::max(std::abs(side[0].Value.Low) + side[0].Bound,
                                          std::abs(side[1].Value.Low) + side[1].Bound)
                               * Grow;
        braid.Long.at(which)   = (Length(way) + (2.0 * braid.Slips.at(which))) * Grow;
        braid.Levels.at(which) = area.Value.High;
        braid.Drifts.at(which) = (std::abs(area.Value.Low) + area.Bound) * Grow;
      }
      const std::array<Number, 3> lambdas = {Sure::Exact(3.0) * middle, -across, outer};
      for (std::size_t form = 0; form < 2; ++form)
      {
        const std::array<std::size_t, 3> rows =
            form == 0 ? std::array<std::size_t, 3>{0, 1, 2} : std::array<std::size_t, 3>{1, 3, 4};
        Pair   turn = {Sure::Exact(0.0), Sure::Exact(0.0)};
        Number base = Sure::Exact(0.0);
        for (std::size_t index = 0; index < rows.size(); ++index)
        {
          const Number &lambda = lambdas.at(index);
          const auto    entry  = rows.at(index);
          turn                 = {turn[0] + (lambda * braid.Sides.at(entry)[0]),
                                  turn[1] + (lambda * braid.Sides.at(entry)[1])};
          base                 = base + (lambda * braid.Areas.at(entry));
        }
        braid.Turns.at(form) = turn;
        braid.Bases.at(form) = base;
      }
      double pace = 0.0;
      for (std::size_t index = 0; index + 1 < poles.size(); ++index)
      {
        const Vector step = poles.at(index + 1) - poles.at(index);
        pace              = std::max({pace, 3.0 * std::abs(step.X), 3.0 * std::abs(step.Y)});
      }
      braid.Stride = 4.0 * pace * Above(2.0);
      Box hull     = {.Min = poles[0], .Max = poles[0]};
      for (const Point pole : poles)
      {
        hull = Grown(hull, pole);
      }
      braid.Hull = hull;
      return braid;
    }

    /**
     *  A cubic's entries on a span: each a polynomial of the span's degree whose Bernstein
     *  coefficient at control point i is the entry at w_i (P_i - P0) and w_i, how far each may be
     *  off, and the most each entry's size may be at a point of the span.
     */
    struct Cords
    {
      std::array<Poled<double, Ample>, Fibres> Values = {}; /**< The coefficients.        */
      std::array<Poled<double, Ample>, Fibres> Offs   = {}; /**< How far each may be off. */
      std::array<double, Fibres>               Most   = {}; /**< Each entry's most size.  */
    };

    /**
     *  A cubic's entries on a span moved to P0. Each is K w - S x G for the moved point G and its
     *  weight w, worked out as the parabola's lines are (arch.py's rounding, apart and carried):
     *  within 4 u of the products' sizes past S's slip times the point, S times the point's slack
     *  and, where K is not nothing, K's drift times the weight and K times the heft; and a step of
     *  the doubles below the normal ones, which Tiny covers, where a product falls among them. At a
     *  point of the span the entry is a convex mix of its coefficients over the weights, so no
     *  larger than the most of their sizes over the lightest weights.
     *
     *  @param  braid The cubic.
     *  @param  piece The span, moved to P0.
     *  @return       The entries; nothing where one may be past Thickest.
     */
    auto Corded(const Braid &braid, const Piece<Ample> &piece) -> std::optional<Cords>
    {
      Cords cords = {};
      for (std::size_t index = 0; index <= piece.Degree; ++index)
      {
        const Vector point  = piece.Moved.at(index);
        const Vector slack  = piece.Offs.at(index);
        const double weight = piece.Weights.at(index);
        for (std::size_t which = 0; which < Fibres; ++which)
        {
          const Vector way   = braid.Ways.at(which);
          const double level = braid.Levels.at(which);
          const double cross = (way.X * point.Y) - (way.Y * point.X);
          const bool   based = which == 0 || which == 1 || which == 3;
          const double off =
              ((std::abs(way.X) * slack.Y)
               + (std::abs(way.Y) * slack.X)
               + (braid.Slips.at(which)
                  * (std::abs(point.X) + std::abs(point.Y) + slack.X + slack.Y))
               + (4.0 * Unit * (std::abs(way.X * point.Y) + std::abs(way.Y * point.X)))
               + (based ? (4.0 * Unit * std::abs(level * weight))
                              + (braid.Drifts.at(which) * (weight + piece.Heft))
                              + (std::abs(level) * piece.Heft)
                        : 0.0)
               + Sure::Tiny)
              * Grow;
          const double value               = based ? (level * weight) - cross : -cross;
          cords.Values.at(which).at(index) = value;
          cords.Offs.at(which).at(index)   = off;
          cords.Most.at(which) = std::max(cords.Most.at(which),
                                          (std::abs(value) + off) / (weight - piece.Heft) * Grow);
          if (!(std::abs(value) + off <= Thickest))
          {
            return std::nullopt;
          }
        }
      }
      return cords;
    }

    /**
     *  A minor of a cubic's matrix on a span, x y - z v for four entries, of twice the span's
     *  degree n, held by its Bernstein coefficients times C(2 n, k): c_k the sum over i + j = k of
     *  C(n, i) C(n, j) (x_i y_j - z_i v_j), and how far each may be off.
     */
    struct Twill
    {
      Terms<Ample> Values = {}; /**< The coefficients, times C(2 n, k). */
      Terms<Ample> Offs   = {}; /**< How far each may be off.           */
    };

    /**
     *  A minor of a cubic's matrix on a span. Each term two products and a difference times the
     *  binomial weights' product, whole and exact, and the running sum of no more than n + 1 terms,
     *  are within ring.py's 16 u of the terms' sizes for every degree to Broadest (its term,
     *  running and wide), past what the entries may be off times each other (its product), and a
     *  step of the doubles below the normal ones a product, which Tiny covers.
     *
     *  @param  cords  The entries.
     *  @param  degree The span's degree, n.
     *  @param  which  The entries x, y, z and v.
     *  @return        The minor.
     */
    auto Twilled(const Cords &cords, std::size_t degree, std::array<std::size_t, 4> which) -> Twill
    {
      const auto &[left, right, other, last] = which;
      const auto value = [&cords](std::size_t entry, std::size_t index) -> double
      {
        return cords.Values.at(entry).at(index);
      };
      const auto off = [&cords](std::size_t entry, std::size_t index) -> double
      {
        return cords.Offs.at(entry).at(index);
      };
      Twill twill = {};
      for (std::size_t power = 0; power <= 2 * degree; ++power)
      {
        double sum  = 0.0;
        double size = 0.0;
        double slip = 0.0;
        for (std::size_t one = power > degree ? power - degree : 0; one <= std::min(power, degree);
             ++one)
        {
          const std::size_t two    = power - one;
          const double      times  = Choose(degree, one) * Choose(degree, two);
          const double      x      = value(left, one);
          const double      y      = value(right, two);
          const double      z      = value(other, one);
          const double      v      = value(last, two);
          sum                     += times * ((x * y) - (z * v));
          size                    += times * (std::abs(x * y) + std::abs(z * v));
          slip                    += times
                                   * ((std::abs(x) * off(right, two))
                                      + (std::abs(y) * off(left, one))
                                      + (off(left, one) * off(right, two))
                                      + (std::abs(z) * off(last, two))
                                      + (std::abs(v) * off(other, one))
                                      + (off(other, one) * off(last, two)));
        }
        twill.Values.at(power) = sum;
        twill.Offs.at(power)   = (slip + (16.0 * Unit * size) + Sure::Tiny) * Grow;
      }
      return twill;
    }

    /**
     *  The most f(q) - f(p) may be where each entry moves by no more than its side's length times r
     *  from its value at p, by the order of r: f's monomials a d g, a e^2, b^2 g, 2 b c e and c^2 d
     *  each expanded, every product bounded by its factors' bounds (braid.py's expanded, scaled,
     *  monomial and gathered).
     *
     *  @param  braid The cubic.
     *  @param  sizes The most each entry's size may be at p.
     *  @return       The terms in r, r^2 and r^3.
     */
    auto Teased(const Braid &braid, const std::array<double, Fibres> &sizes)
        -> std::array<double, 3>
    {
      constexpr std::array<std::array<std::size_t, 3>, 5> Monomials = {
          {{0, 3, 5}, {0, 4, 4}, {1, 1, 5}, {1, 2, 4}, {2, 2, 3}}};
      constexpr std::array<double, 5>   Counts = {1.0, 1.0, 1.0, 2.0, 1.0};
      const std::array<double, Fibres> &way    = braid.Long;
      std::array<double, 3>             orders = {};
      for (std::size_t index = 0; index < Monomials.size(); ++index)
      {
        const std::array<std::size_t, 3> &monomial  = Monomials.at(index);
        const std::size_t                 x         = monomial[0];
        const std::size_t                 y         = monomial[1];
        const std::size_t                 z         = monomial[2];
        const double                      count     = Counts.at(index);
        orders[0]                                  += count
                                                    * ((way.at(x) * sizes.at(y) * sizes.at(z))
                                                       + (sizes.at(x) * way.at(y) * sizes.at(z))
                                                       + (sizes.at(x) * sizes.at(y) * way.at(z)));
        orders[1]                                  += count
                                                    * ((way.at(x) * way.at(y) * sizes.at(z))
                                                       + (way.at(x) * sizes.at(y) * way.at(z))
                                                       + (sizes.at(x) * way.at(y) * way.at(z)));
        orders[2]                                  += count * way.at(x) * way.at(y) * way.at(z);
      }
      return {orders[0] * Grow, orders[1] * Grow, orders[2] * Grow};
    }

    /**
     *  A cubic's entries at a homogeneous point (N, W), in double-double: K W - S x (N - P0 W).
     *
     *  @param  braid The cubic.
     *  @param  point The point's numerator's coordinates and its weight.
     *  @return       The entries.
     */
    auto Carded(const Braid &braid, const std::array<Number, 3> &point)
        -> std::array<Number, Fibres>
    {
      const Number               gapX    = point[0] - (Sure::Exact(braid.Hub.X) * point[2]);
      const Number               gapY    = point[1] - (Sure::Exact(braid.Hub.Y) * point[2]);
      std::array<Number, Fibres> entries = {};
      for (std::size_t which = 0; which < Fibres; ++which)
      {
        const std::array<Number, 2> &side = braid.Sides.at(which);
        entries.at(which) =
            (braid.Areas.at(which) * point[2]) - ((side[0] * gapY) - (side[1] * gapX));
      }
      return entries;
    }

    /**
     *  A cubic's matrix's determinant, a (d g - e^2) - b (b g - c e) + c (b e - c d), in
     *  double-double.
     *
     *  @param  entries The entries a, b, c, d, e and g.
     *  @return         The determinant.
     */
    auto Tufted(const std::array<Number, Fibres> &entries) -> Number
    {
      const auto minor = [&entries](std::array<std::size_t, 4> which) -> Number
      {
        return (entries.at(which[0]) * entries.at(which[1]))
             - (entries.at(which[2]) * entries.at(which[3]));
      };
      return (entries[0] * minor({3, 5, 4, 4}))
           - (entries[1] * minor({1, 5, 2, 4}))
           + (entries[2] * minor({1, 4, 2, 3}));
    }

    /**
     *  The most |grad f| may be where the entries' sizes are no more than some: the sum over the
     *  entries of |df / dE| times |S|, the partial derivatives d g - e^2, 2 (c e - b g), 2 (b e - c
     *  d), a g - c^2, 2 (b c - a e) and a d - b^2 bounded by their terms' sizes (braid.py's
     *  partials).
     *
     *  @param  braid The cubic.
     *  @param  sizes The most each entry's size may be.
     *  @return       The most |grad f| may be.
     */
    auto Pitched(const Braid &braid, const std::array<double, Fibres> &sizes) -> double
    {
      const auto pair = [&sizes](std::array<std::size_t, 4> which) -> double
      {
        return (sizes.at(which[0]) * sizes.at(which[1]))
             + (sizes.at(which[2]) * sizes.at(which[3]));
      };
      const std::array<double, Fibres> partials = {
          pair({3, 5, 4, 4}), pair({1, 5, 2, 4}) * 2.0, pair({1, 4, 2, 3}) * 2.0,
          pair({0, 5, 2, 2}), pair({0, 4, 1, 2}) * 2.0, pair({0, 3, 1, 1})};
      double sum = 0.0;
      for (std::size_t which = 0; which < Fibres; ++which)
      {
        sum += partials.at(which) * braid.Long.at(which);
      }
      return sum * Grow;
    }
  } // namespace

  auto Springing(const Braid &braid, const Piece<Ample> &piece, double reach)
      -> std::optional<Equation<Ample>>
  {
    const std::size_t          degree = piece.Degree;
    const std::optional<Cords> cords  = Corded(braid, piece);
    if (!cords)
    {
      return std::nullopt;
    }
    const std::array<Twill, 3> minors   = {Twilled(*cords, degree, {3, 5, 4, 4}),
                                           Twilled(*cords, degree, {1, 5, 2, 4}),
                                           Twilled(*cords, degree, {1, 4, 2, 3})};
    Equation<Ample>            equation = {.Degree = 3 * degree};
    for (std::size_t power = 0; power <= equation.Degree; ++power)
    {
      double sum  = 0.0;
      double size = 0.0;
      double slip = 0.0;
      for (std::size_t one = power > 2 * degree ? power - (2 * degree) : 0;
           one <= std::min(power, degree); ++one)
      {
        const std::size_t     two   = power - one;
        const double          times = Choose(degree, one);
        std::array<double, 3> terms = {};
        for (std::size_t row = 0; row < minors.size(); ++row)
        {
          const double entry  = cords->Values.at(row).at(one);
          const double moved  = cords->Offs.at(row).at(one);
          const double minor  = minors.at(row).Values.at(two);
          const double drift  = minors.at(row).Offs.at(two);
          terms.at(row)       = entry * minor;
          size               += times * std::abs(entry * minor);
          slip += times * ((std::abs(entry) * drift) + (std::abs(minor) * moved) + (moved * drift));
        }
        sum += times * ((terms[0] - terms[1]) + terms[2]);
      }
      const double below         = Choose(equation.Degree, power);
      const double value         = sum / below;
      equation.Whole.C.at(power) = value;
      equation.First =
          std::max(equation.First, (((slip + (16.0 * Unit * size)) / below) * Grow) + Sure::Tiny);
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
    double heavy      = 0.0;
    for (std::size_t index = 0; index <= degree; ++index)
    {
      heavy = std::max(heavy, piece.Weights.at(index) + piece.Heft);
    }
    const std::array<double, 3> orders = Teased(braid, cords->Most);
    equation.Heavy                     = heavy * heavy * heavy;
    equation.Reach =
        reach * equation.Heavy * (orders[0] + (reach * (orders[1] + (reach * orders[2])))) * Grow;
    return equation;
  }

  auto Keystone(const Braid &braid, const std::array<Number, 3> &point) -> Number
  {
    return Tufted(Carded(braid, point));
  }

  auto Sprung(const Arcade<Braid> &arcade, const Root &root, List<Ribbed> &found) -> bool
  {
    const Braid                 &braid    = *arcade.Carrier;
    const Piece<Ample>          &piece    = *arcade.Part;
    const Honed                  honed    = Honing(arcade, root);
    const std::array<Number, 2> &crossing = honed.Crossing;
    const std::array<Number, 3>  point    = {crossing[0], crossing[1], Sure::Exact(1.0)};
    const Number                 gapX     = crossing[0] - Sure::Exact(braid.Hub.X);
    const Number                 gapY     = crossing[1] - Sure::Exact(braid.Hub.Y);
    const auto                   form     = [&braid, &gapX, &gapY](std::size_t which) -> Number
    {
      const std::array<Number, 2> &turn = braid.Turns.at(which);
      return braid.Bases.at(which) - ((turn[0] * gapY) - (turn[1] * gapX));
    };
    const auto moving = [&honed](const std::array<Number, 2> &turn) -> double
    {
      return Trig::Hypot(Sure::Most(turn[0]), Sure::Most(turn[1])) * honed.Near * Grow;
    };
    const Number                ahead  = form(0);
    const Number                behind = form(1);
    const Number                across = behind - ahead;
    const std::array<Number, 2> slope  = {braid.Turns[1][0] - braid.Turns[0][0],
                                          braid.Turns[1][1] - braid.Turns[0][1]};
    const Number                about =
        Sure::Quotient({.Value = behind.Value, .Bound = behind.Bound + moving(braid.Turns[1])},
                       {.Value = across.Value, .Bound = across.Bound + moving(slope)});
    const Side side = Within({.Start = 0.0, .End = 1.0}, about.Value.High,
                             (std::abs(about.Value.Low) + about.Bound) * Grow);
    if (side != Side::In)
    {
      return side == Side::Out;
    }
    const std::array<Number, Fibres> entries = Carded(braid, point);
    std::array<double, Fibres>       sizes   = {};
    std::ranges::transform(entries, sizes.begin(),
                           [](const Number &entry) -> double { return Sure::Most(entry); });
    const double sine =
        root.Slope / (arcade.Equation.Heavy * Pitched(braid, sizes) * piece.Speed) * Shrink;
    if (!Steep(sine, arcade.Tolerance))
    {
      return false;
    }
    const std::optional<std::array<double, 4>> rounded = Sure::Rounded(std::array{
        crossing[0], crossing[1], Lodged(piece, honed.Place), Sure::Quotient(behind, across)});
    if (rounded)
    {
      found.push_back({.At    = {.X = (*rounded)[0], .Y = (*rounded)[1]},
                       .Along = (*rounded)[2],
                       .Bent  = (*rounded)[3]});
    }
    return rounded.has_value();
  }

  auto Clearing(const Braid &braid, Point point, double clear) -> bool
  {
    const std::array<Number, Fibres> entries =
        Carded(braid, {Sure::Exact(point.X), Sure::Exact(point.Y), Sure::Exact(1.0)});
    std::array<double, Fibres> sizes = {};
    std::ranges::transform(entries, sizes.begin(),
                           [](const Number &entry) -> double { return Sure::Most(entry); });
    const std::array<double, 3> orders = Teased(braid, sizes);
    return Sure::Least(Tufted(entries))
         > clear * (orders[0] + (clear * (orders[1] + (clear * orders[2])))) * Grow;
  }

  auto Braided(const Curves &curves, const Tolerance &tolerance, Intersections &into)
      -> std::optional<bool>
  {
    const std::size_t which  = std::holds_alternative<Bezier::Cubic>(curves[0].get()) ? 0 : 1;
    const auto       *cubic  = std::get_if<Bezier::Cubic>(&curves.at(which).get());
    const auto       *spline = std::get_if<Spline>(&curves.at(1 - which).get());
    if (cubic == nullptr
        || spline == nullptr
        || !(Rational(*spline) || std::cmp_greater(spline->Degree, Highest)))
    {
      return std::nullopt;
    }
    if (std::cmp_greater(spline->Degree, Tightest))
    {
      return false;
    }
    const std::optional<Braid> braid = Braiding(*cubic);
    return braid
        && Trussed(*braid, {.Pair = &curves, .Side = 1 - which, .Tolerance = tolerance}, into);
  }
} // namespace Geoxy::Swift
