/***************************************************************************************************
 *  @file      source/kernel/piece.cpp
 *  @brief     Rational Bézier pieces: cut, bounded, judged flat, and made from curves.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "kernel/piece.hpp"

#include "basis.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "trig.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace Geoxy::Kernel
{
  namespace
  {
    /**
     *  How many cells a span of a low degree takes room for: a degree below it is worked out in
     *  that many, so that what is copied and cleared is no more than it uses.
     */
    constexpr std::size_t Lowly = 8;

    /**
     *  A spline's span as a piece's control points, each the blossom at the span's start d - k
     *  times, then at its end k times: de Boor's triangle whose first d - k rounds take the start.
     *  Those rounds are one triangle's, the start's, worked out once: after r of them the cells
     *  from r on are the blossom's first r rounds, to the bit, and the rounds at the end go on from
     *  a copy of them, each mix taken as Blossom takes it.
     *
     *  @tparam CountT How many cells to take room for: more than the degree.
     *  @param  spline The spline, valid.
     *  @param  span   The index of the knot that starts the span, which is not empty.
     *  @param  piece  The piece, its degree set, whose control points are filled in.
     */
    template <std::size_t CountT>
    void Sectional(const Spline &spline, std::size_t span, Piece &piece)
    {
      const std::size_t               degree = piece.Degree;
      const int                       power  = Leveling(spline, span);
      const std::vector<double>      &knots  = spline.Knots;
      std::array<Homogeneous, CountT> cells  = {};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        cells.at(index) = Pole(spline, span - degree + index, power);
      }
      const auto round = [&knots, span, degree](std::array<Homogeneous, CountT> &work,
                                                std::size_t which, double along) -> void
      {
        for (std::size_t index = degree; index >= which; --index)
        {
          const double start = knots[span - degree + index];
          const double end   = knots[span + 1 + index - which];
          work.at(index) = Mix(work.at(index - 1), work.at(index), (along - start) / (end - start));
        }
      };
      for (std::size_t taken = 0; taken <= degree; ++taken)
      {
        std::array<Homogeneous, CountT> work = cells;
        for (std::size_t which = taken + 1; which <= degree; ++which)
        {
          round(work, which, knots[span + 1]);
        }
        piece.Poles.at(degree - taken) = work.at(degree);
        if (taken < degree)
        {
          round(cells, taken + 1, knots[span]);
        }
      }
    }

    /**
     *  The distance from a point to a segment.
     *
     *  @param  point The point.
     *  @param  start The segment's start.
     *  @param  end   The segment's end.
     *  @return       The distance to the nearest point of the segment.
     */
    auto Away(Point point, Point start, Point end) -> double
    {
      const Vector chord  = end - start;
      const double square = Dot(chord, chord);
      const double along =
          square > 0.0 ? std::clamp(Dot(point - start, chord) / square, 0.0, 1.0) : 0.0;
      return Distance(point, Lerp(start, end, along));
    }

    /**
     *  How far a box reaches past another.
     *
     *  @param  part  The box.
     *  @param  whole The box it should be within.
     *  @return       The most it reaches past it along either axis, or 0.
     */
    auto Excess(const Box &part, const Box &whole) -> double
    {
      return std::max({whole.Min.X - part.Min.X, part.Max.X - whole.Max.X, whole.Min.Y - part.Min.Y,
                       part.Max.Y - whole.Max.Y, 0.0});
    }

    /**
     *  A polynomial's coefficients in Bernstein form, the first degree + 1 used.
     */
    using Coefficients = std::array<double, Terms>;

    /**
     *  A polynomial in Bernstein form: its degree and its coefficients.
     */
    struct Series
    {
      std::size_t  Degree = 0;  /**< The degree.       */
      Coefficients Values = {}; /**< The coefficients. */
    };

    /**
     *  The binomial coefficients of every degree under Terms, worked out once, each row as the
     *  multiplicative formula gives it: degree choose 0 to degree choose degree.
     */
    constexpr std::array<Coefficients, Terms> Binomials = [] -> std::array<Coefficients, Terms>
    {
      std::array<Coefficients, Terms> table = {};
      for (std::size_t degree = 0; degree < Terms; ++degree)
      {
        Coefficients &row = table.at(degree);
        row.at(0)         = 1.0;
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
     *  The product of two polynomials in Bernstein form.
     *
     *  @param  lhs One polynomial.
     *  @param  rhs The other; the degrees together less than Terms.
     *  @return     Their product, of the two degrees together.
     */
    auto operator*(const Series &lhs, const Series &rhs) -> Series
    {
      Series              product = {.Degree = lhs.Degree + rhs.Degree, .Values = {}};
      const Coefficients &left    = Binomials.at(lhs.Degree);
      const Coefficients &right   = Binomials.at(rhs.Degree);
      const Coefficients &both    = Binomials.at(product.Degree);
      for (std::size_t first = 0; first <= lhs.Degree; ++first)
      {
        for (std::size_t second = 0; second <= rhs.Degree; ++second)
        {
          product.Values.at(first + second) +=
              left.at(first) * right.at(second) * lhs.Values.at(first) * rhs.Values.at(second);
        }
      }
      for (std::size_t index = 0; index <= product.Degree; ++index)
      {
        product.Values.at(index) /= both.at(index);
      }
      return product;
    }

    /**
     *  The sum of two polynomials in Bernstein form of one degree.
     *
     *  @param  lhs One polynomial.
     *  @param  rhs The other, of the same degree.
     *  @return     Their sum.
     */
    auto operator+(const Series &lhs, const Series &rhs) -> Series
    {
      Series sum = lhs;
      for (std::size_t index = 0; index <= sum.Degree; ++index)
      {
        sum.Values.at(index) += rhs.Values.at(index);
      }
      return sum;
    }

    /**
     *  The difference of two polynomials in Bernstein form of one degree.
     *
     *  @param  lhs One polynomial.
     *  @param  rhs The other, of the same degree.
     *  @return     The first less the second.
     */
    auto operator-(const Series &lhs, const Series &rhs) -> Series
    {
      Series difference = lhs;
      for (std::size_t index = 0; index <= difference.Degree; ++index)
      {
        difference.Values.at(index) -= rhs.Values.at(index);
      }
      return difference;
    }

    /**
     *  A polynomial's coefficients made positive: worked through the same products, a bound on each
     *  coefficient of the result's size, to tell rounding from a value.
     *
     *  @param  series The polynomial.
     *  @return        Its coefficients' magnitudes.
     */
    auto Sized(const Series &series) -> Series
    {
      Series sized = series;
      for (double &value : sized.Values)
      {
        value = std::abs(value);
      }
      return sized;
    }

    /**
     *  The integral of a polynomial in Bernstein form over 0 to 1: the mean of its coefficients.
     *
     *  @param  series The polynomial.
     *  @return        The integral.
     */
    auto Integral(const Series &series) -> double
    {
      double sum = 0.0;
      for (const double value : std::span{series.Values}.first(series.Degree + 1))
      {
        sum += value;
      }
      return sum / static_cast<double>(series.Degree + 1);
    }

    /**
     *  The derivative of a polynomial in Bernstein form.
     *
     *  @pre           The degree is 1 at least.
     *  @param  series The polynomial.
     *  @return        Its derivative, of one degree less.
     */
    auto Slope(const Series &series) -> Series
    {
      Series     slope = {.Degree = series.Degree - 1, .Values = {}};
      const auto scale = static_cast<double>(series.Degree);
      for (std::size_t index = 0; index < series.Degree; ++index)
      {
        slope.Values.at(index) = scale * (series.Values.at(index + 1) - series.Values.at(index));
      }
      return slope;
    }

    /**
     *  The polynomial that pairs of a piece's control values make together, each pair of a lesser
     *  index with a greater giving one value.
     *
     *  B(n, i) B'(n, j) - B'(n, i) B(n, j) is (j - i) C(n, i) C(n, j) / C(2 n - 2, i + j - 1) times
     *  B(2 n - 2, i + j - 1) for i less than j, so one polynomial times another's derivative, less
     *  its own derivative times the other, is the pairs' differences so scaled and summed.
     *
     *  @tparam PairT  What gives a pair's value.
     *  @pre           The degree is 1 at least.
     *  @param  degree The degree n of what is paired.
     *  @param  pair   The value of the pair i, j, i less than j.
     *  @return        The polynomial of the pairs, of degree 2 n - 2.
     */
    template <typename PairT>
    auto Paired(std::size_t degree, const PairT &pair) -> Series
    {
      Series              paired = {.Degree = (2 * degree) - 2, .Values = {}};
      const Coefficients &single = Binomials.at(degree);
      const Coefficients &both   = Binomials.at(paired.Degree);
      for (std::size_t first = 0; first < degree; ++first)
      {
        for (std::size_t second = first + 1; second <= degree; ++second)
        {
          const double scale =
              static_cast<double>(second - first) * single.at(first) * single.at(second);
          paired.Values.at(first + second - 1) += scale * pair(first, second);
        }
      }
      for (std::size_t index = 0; index <= paired.Degree; ++index)
      {
        paired.Values.at(index) /= both.at(index);
      }
      return paired;
    }

    /**
     *  One polynomial times another's derivative, less its own derivative times the other: the
     *  numerator of the derivative of the second over the first.
     *
     *  Worked out pair by pair, a pair of small coefficients keeps its digits beside a pair of
     *  large ones. The same from the polynomials' slopes loses them: each product then holds the
     *  large pairs' terms, which cancel and leave their rounding where the small pair's stood, as
     *  at the end of a piece whose weights there are light beside the rest.
     *
     *  @pre        The degree is 1 at least.
     *  @param  lhs The first polynomial.
     *  @param  rhs The second, of the same degree.
     *  @return     lhs rhs' - lhs' rhs, of degree 2 n - 2.
     */
    auto Wronskian(const Series &lhs, const Series &rhs) -> Series
    {
      return Paired(lhs.Degree,
                    [&lhs, &rhs](std::size_t first, std::size_t second) -> double
                    {
                      return (lhs.Values.at(first) * rhs.Values.at(second))
                           - (lhs.Values.at(second) * rhs.Values.at(first));
                    });
    }

    /**
     *  A bound on each coefficient of two polynomials' Wronskian: the same pairs, their products
     *  made positive and added.
     *
     *  @pre        The degree is 1 at least.
     *  @param  lhs The first polynomial.
     *  @param  rhs The second, of the same degree.
     *  @return     The bound, of degree 2 n - 2.
     */
    auto Bounding(const Series &lhs, const Series &rhs) -> Series
    {
      return Paired(lhs.Degree,
                    [&lhs, &rhs](std::size_t first, std::size_t second) -> double
                    {
                      return std::abs(lhs.Values.at(first) * rhs.Values.at(second))
                           + std::abs(lhs.Values.at(second) * rhs.Values.at(first));
                    });
    }

    /**
     *  How a span's Bézier control points are made of its poles: how much of each pole each control
     *  point holds, and the poles' weights.
     */
    struct Sharing
    {
      /**
       *  A row of shares: a share a pole, and room for the two more a blossom tells at a time.
       */
      using Row = std::array<double, MaxDegree + 3>;

      std::size_t                       Degree  = 0;  /**< The degree.                         */
      std::size_t                       First   = 0;  /**< The index of the span's first pole. */
      std::array<Row, MaxDegree + 1>    Shares  = {}; /**< A row a control point.              */
      std::array<double, MaxDegree + 1> Weights = {}; /**< The poles' weights, scaled.         */
    };

    /**
     *  How a span's Bézier control points are made of its poles.
     *
     *  Each control point is the blossom of the poles at the span's ends, so its shares are the
     *  blossom of poles told apart, three at a time, one a coordinate. The weights are all times
     *  one power of two that sets the heaviest as far above 1 as the lightest is below, so that a
     *  product of two neither overflows nor comes to nothing.
     *
     *  @pre           The spline is valid and the span is not empty.
     *  @param  spline The spline.
     *  @param  span   The index of the knot that starts the span.
     *  @return        The shares and the weights.
     */
    auto Shared(const Spline &spline, std::size_t span) -> Sharing
    {
      const Basis basis   = Spanned(spline);
      Sharing     sharing = {
          .Degree = basis.Degree, .First = span - basis.Degree, .Shares = {}, .Weights = {}};
      std::array<double, MaxDegree> values{};
      const std::span<double>       blossom = std::span{values}.first(basis.Degree);
      for (std::size_t row = 0; row <= basis.Degree; ++row)
      {
        std::ranges::fill(blossom, spline.Knots[span]);
        std::ranges::fill(blossom.last(row), spline.Knots[span + 1]);
        for (std::size_t pole = 0; pole <= basis.Degree; pole += 3)
        {
          const std::size_t which = sharing.First + pole;
          const auto        told  = [which](std::size_t index) -> Homogeneous
          {
            return {.X = index == which ? 1.0 : 0.0,
                    .Y = index == which + 1 ? 1.0 : 0.0,
                    .W = index == which + 2 ? 1.0 : 0.0};
          };
          const Homogeneous shares = Blossom(basis, told, {.Span = span, .Values = blossom});
          Sharing::Row     &into   = sharing.Shares.at(row);
          into.at(pole)            = shares.X;
          into.at(pole + 1)        = shares.Y;
          into.at(pole + 2)        = shares.W;
        }
      }
      std::ranges::fill(sharing.Weights, 1.0);
      if (!spline.Weights.empty())
      {
        const auto about = std::span{spline.Weights}.subspan(sharing.First, basis.Degree + 1);
        const auto [lightest, heaviest] = std::ranges::minmax(about);
        const int power                 = -(std::ilogb(heaviest) + std::ilogb(lightest)) / 2;
        std::ranges::transform(about, sharing.Weights.begin(), [power](double weight) -> double
                               { return std::ldexp(weight, power); });
      }
      return sharing;
    }

    /**
     *  The polynomial the pairs of a span's poles make together, each pair weighed by its weights
     *  and by how much of the two each pair of Bézier control points holds.
     *
     *  A pair of control points is every pair of poles, by a minor of the shares, never negative as
     *  the blossom's are: what a pair of poles gives is so carried to the control points with
     *  nothing cancelling but what the poles themselves cancel. A control point, a sum of poles
     *  rounded, has lost a light pole beside a heavy one; the poles' pairs have not.
     *
     *  @tparam PairT   What gives a pair of poles' value.
     *  @param  sharing How the span's control points are made of its poles.
     *  @param  pair    The value of the poles near and far, near the lesser index, without their
     *                  weights.
     *  @return         The polynomial, of degree 2 n - 2.
     */
    template <typename PairT>
    auto Polar(const Sharing &sharing, const PairT &pair) -> Series
    {
      const std::size_t                                            degree = sharing.Degree;
      std::array<std::array<double, MaxDegree + 1>, MaxDegree + 1> weighed{};
      for (std::size_t near = 0; near < degree; ++near)
      {
        for (std::size_t far = near + 1; far <= degree; ++far)
        {
          weighed.at(near).at(far) =
              sharing.Weights.at(near) * sharing.Weights.at(far) * pair(near, far);
        }
      }
      return Paired(degree,
                    [&sharing, &weighed, degree](std::size_t lower, std::size_t upper) -> double
                    {
                      const auto &low  = sharing.Shares.at(lower);
                      const auto &high = sharing.Shares.at(upper);
                      double      sum  = 0.0;
                      for (std::size_t near = 0; near < degree; ++near)
                      {
                        for (std::size_t far = near + 1; far <= degree; ++far)
                        {
                          sum += ((low.at(near) * high.at(far)) - (low.at(far) * high.at(near)))
                               * weighed.at(near).at(far);
                        }
                      }
                      return sum;
                    });
    }

    /**
     *  A stretch of the parameter, from its start to its end.
     */
    struct Stretch
    {
      double Start = 0.0; /**< Where it starts. */
      double End   = 1.0; /**< Where it ends.   */
    };

    /**
     *  Whether two values are of opposite signs, neither zero. Their product tells as much only
     *  where it does not underflow: of two values under 2^-538 it is zero, and tells nothing.
     *
     *  @param  lhs One value.
     *  @param  rhs The other.
     *  @return     Whether one is positive and the other negative.
     */
    auto Opposed(double lhs, double rhs) -> bool
    {
      return (lhs < 0.0 && rhs > 0.0) || (lhs > 0.0 && rhs < 0.0);
    }

    /**
     *  Whether a piece found nothing and has nothing to tell of its ends: one of no length, or one
     *  the function is zero all along. A chain's function turns across it as across a seam.
     *
     *  @param  solution What the piece found.
     *  @return          Whether it has no roots and both its end values are zero.
     */
    auto Blank(const Solution &solution) -> bool
    {
      return solution.Roots.empty() && Zero(solution.Start) && Zero(solution.End);
    }

    /**
     *  How many times a run of coefficients changes sign, zeros skipped: at least as many as the
     *  roots of the polynomial within, and as many as them when it is 0 or 1.
     *
     *  @param  values The coefficients.
     *  @return        The count.
     */
    auto Changes(std::span<const double> values) -> std::size_t
    {
      std::size_t changes = 0;
      double      last    = 0.0;
      for (const double value : values)
      {
        if (Zero(value))
        {
          continue;
        }
        if (!Zero(last) && (value > 0.0) != (last > 0.0))
        {
          ++changes;
        }
        last = value;
      }
      return changes;
    }

    /**
     *  A polynomial's value and slope at a parameter.
     */
    struct Valued
    {
      double Value = 0.0; /**< The value.          */
      double Slope = 0.0; /**< Its rate of change. */
    };

    /**
     *  A polynomial's value and slope, by de Casteljau's algorithm: the slope from the last two
     *  values but one round, times the degree.
     *
     *  @param  values The coefficients, two at least.
     *  @param  along  The parameter.
     *  @param  work   Room for as many values, overwritten.
     *  @return        The value and the slope.
     */
    auto Value(std::span<const double> values, double along, std::span<double> work) -> Valued
    {
      const std::size_t degree = values.size() - 1;
      std::ranges::copy(values, work.begin());
      for (std::size_t round = 1; round < degree; ++round)
      {
        for (std::size_t index = 0; index + round <= degree; ++index)
        {
          work[index] += along * (work[index + 1] - work[index]);
        }
      }
      const double slope = static_cast<double>(degree) * (work[1] - work[0]);
      return {.Value = work[0] + (along * (work[1] - work[0])), .Slope = slope};
    }

    /**
     *  A part's coefficients halved, by de Casteljau's algorithm: those over its first half and
     *  those over its second.
     *
     *  @param cells  The part's coefficients, overwritten.
     *  @param before The first half's, as many.
     *  @param after  The second half's, as many.
     */
    void HalvedCells(std::span<double> cells, std::span<double> before, std::span<double> after)
    {
      const std::size_t degree = cells.size() - 1;
      before[0]                = cells[0];
      after[degree]            = cells[degree];
      for (std::size_t round = 1; round <= degree; ++round)
      {
        for (std::size_t index = 0; index + round <= degree; ++index)
        {
          cells[index] = (cells[index] + cells[index + 1]) / 2.0;
        }
        before[round]         = cells[0];
        after[degree - round] = cells[degree - round];
      }
    }

    /**
     *  The next place to try for a root, by Newton's method from where the polynomial was just
     *  worked out: one step of a double toward the root where the step is smaller, halfway across
     *  the bracket where the step would leave it or the slope is zero.
     *
     *  @param  bracket The stretch the root is known in, its ends where the sign differs.
     *  @param  middle  Where the polynomial was worked out.
     *  @param  valued  Its value and slope there.
     *  @return         The place.
     */
    auto Guessed(Stretch bracket, double middle, Valued valued) -> double
    {
      if (Zero(valued.Slope))
      {
        return bracket.Start + ((bracket.End - bracket.Start) / 2.0);
      }
      double guess = middle - (valued.Value / valued.Slope);
      if (Exactly(guess, middle))
      {
        guess = std::nextafter(middle, Exactly(middle, bracket.End) ? bracket.Start : bracket.End);
      }
      const bool inside = bracket.Start < guess && guess < bracket.End;
      if (!inside)
      {
        guess = bracket.Start + ((bracket.End - bracket.Start) / 2.0);
      }
      return guess;
    }

    /**
     *  The root of a polynomial within a stretch whose ends differ in sign: Newton's method kept in
     *  the bracket, until the bracket's ends are neighbouring doubles, as bisection would end.
     *
     *  Newton's method is trusted for Guided steps and no more. Where the slope overflows, its step
     *  is nothing and each guess is one double on from the last, which would walk the whole bracket
     *  a double at a time; past Guided steps the bracket is halved instead, and ends where it would
     *  have.
     *
     *  @param  values  The polynomial's coefficients.
     *  @param  stretch The stretch.
     *  @param  work    Room for as many values as the coefficients, overwritten; its first, on
     *                  entry, the sign at the stretch's start.
     *  @return         The root, to the last bit.
     */
    auto Bisected(std::span<const double> values, Stretch stretch, std::span<double> work) -> double
    {
      constexpr int Guided  = 64;
      const bool    rising  = work[0] < 0.0;
      Stretch       bracket = stretch;
      double        middle  = bracket.Start + ((bracket.End - bracket.Start) / 2.0);
      for (int step = 0; bracket.Start < middle && middle < bracket.End; ++step)
      {
        const Valued valued = Value(values, middle, work);
        if (Zero(valued.Value))
        {
          return middle;
        }
        if ((valued.Value > 0.0) == rising)
        {
          bracket.End = middle;
        }
        else
        {
          bracket.Start = middle;
        }
        middle = step < Guided ? Guessed(bracket, middle, valued)
                               : bracket.Start + ((bracket.End - bracket.Start) / 2.0);
      }
      return bracket.Start + ((bracket.End - bracket.Start) / 2.0);
    }

    /**
     *  A point turned about the origin by whole quarter turns, exactly.
     *
     *  @param  point    The point.
     *  @param  quarters How many quarter turns, counterclockwise.
     *  @return          The point turned.
     */
    auto Turned(Point point, std::size_t quarters) -> Point
    {
      for (std::size_t turn = 0; turn < quarters; ++turn)
      {
        point = {.X = -point.Y, .Y = point.X};
      }
      return point;
    }

    /**
     *  A piece seen from a point: its homogeneous coordinates less the point's times the weight,
     *  and its weight, as polynomials.
     */
    struct Relative
    {
      Series Across = {}; /**< Along x, less the point's. */
      Series Rising = {}; /**< Along y, less the point's. */
      Series Weight = {}; /**< The weight.                */
    };

    /**
     *  A piece seen from a point.
     *
     *  @param  piece The piece.
     *  @param  point The point.
     *  @return       Its coordinates less the point's, and its weight.
     */
    auto Toward(const Piece &piece, Point point) -> Relative
    {
      Relative   from  = {.Across = {.Degree = piece.Degree, .Values = {}},
                          .Rising = {.Degree = piece.Degree, .Values = {}},
                          .Weight = {.Degree = piece.Degree, .Values = {}}};
      const auto poles = Poles(piece);
      for (std::size_t index = 0; index <= piece.Degree; ++index)
      {
        from.Across.Values.at(index) = poles[index].X - (point.X * poles[index].W);
        from.Rising.Values.at(index) = poles[index].Y - (point.Y * poles[index].W);
        from.Weight.Values.at(index) = poles[index].W;
      }
      return from;
    }

    /**
     *  Whether a piece has one weight throughout: a polynomial piece, whatever the weight.
     *
     *  @param  piece The piece.
     *  @return       Whether every control point's weight is the first's, and that positive.
     */
    auto Uniform(const Piece &piece) -> bool
    {
      const double weight = Poles(piece).front().W;
      return weight > 0.0
          && std::ranges::all_of(Poles(piece), [weight](const Homogeneous &pole) -> bool
                                 { return Exactly(pole.W, weight); });
    }

    /**
     *  A piece seen from a point, its coefficients made positive.
     *
     *  @param  from The piece seen from the point.
     *  @return      Its coefficients' magnitudes.
     */
    auto Sized(const Relative &from) -> Relative
    {
      return {
          .Across = Sized(from.Across), .Rising = Sized(from.Rising), .Weight = Sized(from.Weight)};
    }

    /**
     *  Two polynomials scaled by one power of two, the largest of their coefficients brought near a
     *  power of two.
     *
     *  A piece cut from beside a heavy weight has weights far under 1 and a derivative's numerator
     *  far under those: their product underflows where either alone does not. Scaled so, exactly
     *  and each pair by one positive number, the product keeps its signs and its roots.
     *
     *  @param first  One polynomial, scaled in place.
     *  @param second The other, scaled by as much.
     *  @param toward The power of two the largest coefficient is brought near.
     */
    void Leveled(Series &first, Series &second, int toward)
    {
      double most = 0.0;
      for (const Series *series : {&first, &second})
      {
        for (const double value : std::span{series->Values}.first(series->Degree + 1))
        {
          most = std::max(most, std::abs(value));
        }
      }
      if (!(most > 0.0) || !std::isfinite(most))
      {
        return;
      }
      const int power = toward - std::ilogb(most);
      for (Series *series : {&first, &second})
      {
        for (double &value : std::span{series->Values}.first(series->Degree + 1))
        {
          value = std::ldexp(value, power);
        }
      }
    }

    /**
     *  A piece seen from a point and its derivative's numerator, each scaled to near 1.
     */
    struct Trodden
    {
      Relative From   = {}; /**< The piece seen from the point. */
      Series   Across = {}; /**< The numerator along x.         */
      Series   Rising = {}; /**< The numerator along y.         */
    };

    /**
     *  A piece and its derivative's numerator made ready to be multiplied.
     *
     *  @param  tread The piece and the numerator.
     *  @param  point The point the piece is seen from.
     *  @return       Both, each scaled by a power of two of its own.
     */
    auto Readied(const Tread &tread, Point point) -> Trodden
    {
      Trodden ready = {.From   = Toward(tread.Shape, point),
                       .Across = {.Degree = tread.Speed.Degree, .Values = tread.Speed.Across},
                       .Rising = {.Degree = tread.Speed.Degree, .Values = tread.Speed.Rising}};
      /**
       *  The numerator's coefficients may lie as far apart as the weights' squares: its largest is
       *  set high, so that its least, at a piece's light end, is still a number.
       */
      constexpr int High = 800;
      Leveled(ready.From.Across, ready.From.Rising, 0);
      Leveled(ready.Across, ready.Rising, High);
      return ready;
    }

    /**
     *  Whether a polynomial worked out from a piece is zero to within its rounding: every
     *  coefficient within a few steps of a double of its bound.
     *
     *  @param  value The polynomial.
     *  @param  bound A bound on each coefficient's size.
     *  @return       Whether it is.
     */
    auto Negligible(const Series &value, const Series &bound) -> bool
    {
      const double noise = 64.0 * std::numeric_limits<double>::epsilon();
      bool         flat  = true;
      for (std::size_t index = 0; index <= value.Degree; ++index)
      {
        flat = flat && std::abs(value.Values.at(index)) <= noise * bound.Values.at(index);
      }
      return flat;
    }

    /**
     *  The roots of a polynomial worked out from a piece, unless it is zero to within its rounding.
     *
     *  @param  value  The polynomial.
     *  @param  bound  The same worked through magnitudes: a bound on each coefficient's size.
     *  @param  budget How many halvings the roots may take.
     *  @return        Its roots and its values at the ends; none, and both values 0, when every
     *                 coefficient is within its rounding of 0.
     */
    auto Unless(const Series &value, const Series &bound, std::size_t budget) -> Solution
    {
      if (Negligible(value, bound))
      {
        return {};
      }
      const std::span<const double> values = std::span{value.Values}.first(value.Degree + 1);
      return {.Roots = Roots(values, budget), .Start = values.front(), .End = values.back()};
    }

    /**
     *  The zeros of a polynomial worked out from a piece, and where it turns, unless it is zero to
     *  within its rounding.
     *
     *  @param  value  The polynomial.
     *  @param  bound  A bound on each coefficient's size.
     *  @param  budget How many halvings the roots may take.
     *  @return        Its zeros, its end values, and its derivative's zeros.
     */
    auto Turned(const Series &value, const Series &bound, std::size_t budget) -> Zeros
    {
      if (Negligible(value, bound))
      {
        return {.Crossing = {}, .Turns = {}, .On = true};
      }
      Zeros zeros = {.Crossing = Unless(value, bound, budget), .Turns = {}};
      if (value.Degree >= 2)
      {
        const Series slope = Slope(value);
        zeros.Turns        = Roots(std::span{slope.Values}.first(slope.Degree + 1), budget);
      }
      return zeros;
    }
  } // namespace

  auto Split(const Piece &piece, double along) -> std::array<Piece, 2>
  {
    std::array<Piece, 2> halves = {Piece{.Degree = piece.Degree, .Poles = {}},
                                   Piece{.Degree = piece.Degree, .Poles = {}}};
    const std::size_t    degree = piece.Degree;

    std::array<Homogeneous, MaxDegree + 1> work = piece.Poles;
    const std::span<Homogeneous>           cells{work};
    const std::span<Homogeneous>           before{halves[0].Poles};
    const std::span<Homogeneous>           after{halves[1].Poles};
    before[0]     = cells[0];
    after[degree] = cells[degree];
    for (std::size_t round = 1; round <= degree; ++round)
    {
      for (std::size_t index = 0; index + round <= degree; ++index)
      {
        cells[index] = Mix(cells[index], cells[index + 1], along);
      }
      before[round]         = cells[0];
      after[degree - round] = cells[degree - round];
    }
    return halves;
  }

  auto Sampled(const Piece &piece, double along) -> Sample
  {
    const std::size_t                      degree = piece.Degree;
    std::array<Homogeneous, MaxDegree + 1> work   = piece.Poles;
    const std::span<Homogeneous>           cells{work};
    for (std::size_t round = 1; round < degree; ++round)
    {
      for (std::size_t index = 0; index + round <= degree; ++index)
      {
        cells[index] = Mix(cells[index], cells[index + 1], along);
      }
    }
    const Homogeneous here   = Mix(cells[0], cells[1], along);
    const auto        scale  = static_cast<double>(degree);
    const Homogeneous rate   = {.X = scale * (cells[1].X - cells[0].X),
                                .Y = scale * (cells[1].Y - cells[0].Y),
                                .W = scale * (cells[1].W - cells[0].W)};
    const double      square = here.W * here.W;
    return {.At  = Projected(here),
            .Way = {.X = ((rate.X * here.W) - (here.X * rate.W)) / square,
                    .Y = ((rate.Y * here.W) - (here.Y * rate.W)) / square}};
  }

  auto Reached(const Piece &piece, double along) -> Point
  {
    const std::size_t                      degree = piece.Degree;
    std::array<Homogeneous, MaxDegree + 1> work   = piece.Poles;
    const std::span<Homogeneous>           cells{work};
    for (std::size_t round = 1; round <= degree; ++round)
    {
      for (std::size_t index = 0; index + round <= degree; ++index)
      {
        cells[index] = Mix(cells[index], cells[index + 1], along);
      }
    }
    return Projected(cells[0]);
  }

  auto Hull(const Piece &piece) -> Box
  {
    Box box;
    for (const Homogeneous &pole : Poles(piece))
    {
      box = Grown(box, Projected(pole));
    }
    return box;
  }

  auto Flat(const Piece &piece, double tolerance) -> bool
  {
    const Point start = Start(piece);
    const Point end   = End(piece);
    return std::ranges::all_of(Poles(piece), [&](const Homogeneous &pole) -> bool
                               { return Away(Projected(pole), start, end) <= tolerance; });
  }

  auto Bound(const Piece &piece, Box box, std::size_t budget) -> Box
  {
    List<Piece> pending = {piece};
    while (!pending.empty())
    {
      const Piece top = pending.back();
      pending.pop_back();
      const Box hull = Hull(top);
      if (Contains(box, hull.Min) && Contains(box, hull.Max))
      {
        continue;
      }
      const double allowed = 4.0 * std::max(Grain(hull.Min), Grain(hull.Max));
      if (budget == 0 || !(Excess(hull, box) > allowed))
      {
        box = Joined(box, hull);
        continue;
      }
      --budget;
      const std::array<Piece, 2> halves = Split(top, 0.5);
      box                               = Grown(box, End(halves[0]));
      pending.push_back(halves[1]);
      pending.push_back(halves[0]);
    }
    return box;
  }

  auto Roots(std::span<const double> values, std::size_t budget) -> List<double>
  {
    const std::size_t count = values.size();

    List<double> roots;
    if (Zero(values.front()))
    {
      roots.push_back(0.0);
    }
    if (Zero(values.back()))
    {
      roots.push_back(1.0);
    }

    /**
     *  The parts yet to look into, their coefficients one part after another, as many each as the
     *  polynomial's; and room for the part in hand and its two halves.
     */
    List<double>            pending(values.begin(), values.end());
    List<Stretch>           stretches = {Stretch{}};
    List<double>            work(3 * count);
    const std::span<double> part   = std::span{work}.first(count);
    const std::span<double> before = std::span{work}.subspan(count, count);
    const std::span<double> after  = std::span{work}.last(count);
    while (!stretches.empty())
    {
      const Stretch top = stretches.back();
      stretches.pop_back();
      std::ranges::copy(std::span{pending}.last(count), part.begin());
      pending.resize(pending.size() - count);
      const std::size_t changes = Changes(part);
      if (changes == 0)
      {
        continue;
      }
      if ((changes == 1 && !Zero(part.front()) && !Zero(part.back())) || budget == 0)
      {
        if (Opposed(part.front(), part.back()))
        {
          before[0] = part.front();
          roots.push_back(Bisected(values, top, before));
        }
        continue;
      }
      --budget;
      const double middle = top.Start + ((top.End - top.Start) / 2.0);
      HalvedCells(part, before, after);
      if (Zero(after[0]))
      {
        roots.push_back(middle);
      }
      pending.insert(pending.end(), after.begin(), after.end());
      stretches.push_back({.Start = middle, .End = top.End});
      pending.insert(pending.end(), before.begin(), before.end());
      stretches.push_back({.Start = top.Start, .End = middle});
    }
    Ordered(roots);
    return roots;
  }

  auto Chained(std::span<const Solution> solutions, bool closed) -> List<Spot>
  {
    const std::size_t count = solutions.size();
    List<Spot>        spots;
    for (std::size_t index = 0; index < count; ++index)
    {
      const Solution &solution = solutions[index];
      for (const double root : solution.Roots)
      {
        spots.push_back({.Index = index, .Along = root});
      }
      if (Blank(solution))
      {
        continue;
      }
      std::size_t next = index + 1;
      while (next < count && Blank(solutions[next]))
      {
        ++next;
      }
      if (next == count && closed)
      {
        next = 0;
        while (Blank(solutions[next]))
        {
          ++next;
        }
      }
      if (next < count && Opposed(solution.End, solutions[next].Start))
      {
        spots.push_back({.Index = next, .Along = 0.0});
      }
    }
    return spots;
  }

  auto Roots(const Piece &piece, Level level, std::size_t budget) -> Solution
  {
    Coefficients            values = {};
    const std::span<double> run    = std::span{values}.first(piece.Degree + 1);
    const auto              poles  = Poles(piece);
    for (std::size_t index = 0; index <= piece.Degree; ++index)
    {
      run[index] = Coordinate(poles[index], level.Along) - (level.Value * poles[index].W);
    }
    return {.Roots = Roots(run, budget), .Start = run.front(), .End = run.back()};
  }

  auto Perpendiculars(const Piece &piece, Point point, std::size_t budget) -> Solution
  {
    const Relative from  = Toward(piece, point);
    const Relative sized = Sized(from);
    if (Uniform(piece))
    {
      /**
       *  One weight throughout: W' is zero and W a constant, so the numerator is W times U U' + V
       *  V', of degree 2 n - 1, with the same roots.
       */
      const Series turning =
          (from.Across * Slope(from.Across)) + (from.Rising * Slope(from.Rising));
      const Series bound =
          (sized.Across * Sized(Slope(from.Across))) + (sized.Rising * Sized(Slope(from.Rising)));
      return Unless(turning, bound, budget);
    }
    /**
     *  W (U U' + V V') - W' (U^2 + V^2) is U (W U' - W' U) + V (W V' - W' V), and W U' - W' U is W
     *  X' - W' X whatever the point: the numerator of the piece's own derivative, pair by pair.
     */
    const Relative whole   = Toward(piece, Point{});
    const Relative large   = Sized(whole);
    const Series   turning = (from.Across * Wronskian(whole.Weight, whole.Across))
                           + (from.Rising * Wronskian(whole.Weight, whole.Rising));
    const Series   bound   = (sized.Across * Bounding(large.Weight, large.Across))
                           + (sized.Rising * Bounding(large.Weight, large.Rising));
    return Unless(turning, bound, budget);
  }

  auto Derived(const Piece &piece) -> Hodograph
  {
    const Relative whole  = Toward(piece, Point{});
    const Series   across = Wronskian(whole.Weight, whole.Across);
    const Series   rising = Wronskian(whole.Weight, whole.Rising);
    return {.Degree = across.Degree, .Across = across.Values, .Rising = rising.Values};
  }

  auto Halved(const Hodograph &speed) -> std::array<Hodograph, 2>
  {
    const std::size_t        count  = speed.Degree + 1;
    std::array<Hodograph, 2> halves = {
        Hodograph{.Degree = speed.Degree, .Across = {}, .Rising = {}},
        Hodograph{.Degree = speed.Degree, .Across = {}, .Rising = {}}};
    Coefficients across = speed.Across;
    Coefficients rising = speed.Rising;
    HalvedCells(std::span{across}.first(count), std::span{halves[0].Across}.first(count),
                std::span{halves[1].Across}.first(count));
    HalvedCells(std::span{rising}.first(count), std::span{halves[0].Rising}.first(count),
                std::span{halves[1].Rising}.first(count));
    return halves;
  }

  auto Perpendiculars(const Tread &tread, Point point, std::size_t budget) -> Solution
  {
    const Trodden  ready = Readied(tread, point);
    const Relative sized = Sized(ready.From);
    return Unless((ready.From.Across * ready.Across) + (ready.From.Rising * ready.Rising),
                  (sized.Across * Sized(ready.Across)) + (sized.Rising * Sized(ready.Rising)),
                  budget);
  }

  auto Touches(const Tread &tread, Point point, std::size_t budget) -> Solution
  {
    const Trodden  ready = Readied(tread, point);
    const Relative sized = Sized(ready.From);
    return Unless((ready.From.Across * ready.Rising) - (ready.From.Rising * ready.Across),
                  (sized.Across * Sized(ready.Rising)) + (sized.Rising * Sized(ready.Across)),
                  budget);
  }

  auto Closest(const Piece &piece, Point point) -> double
  {
    List<double> alongs = Perpendiculars(piece, point, Budget).Roots;
    alongs.push_back(0.0);
    alongs.push_back(1.0);
    double best  = 0.0;
    double least = Distance(point, Start(piece));
    for (const double along : alongs)
    {
      const double away = Distance(point, Sampled(piece, along).At);
      if (away < least)
      {
        best  = along;
        least = away;
      }
    }
    return best;
  }

  auto Touches(const Piece &piece, Point point, std::size_t budget) -> Solution
  {
    const Relative from  = Toward(piece, point);
    const Relative sized = Sized(from);
    if (Uniform(piece))
    {
      const Series turning =
          (from.Across * Slope(from.Rising)) - (from.Rising * Slope(from.Across));
      const Series bound =
          (sized.Across * Sized(Slope(from.Rising))) + (sized.Rising * Sized(Slope(from.Across)));
      return Unless(turning, bound, budget);
    }
    return Unless(Wronskian(from.Across, from.Rising), Bounding(sized.Across, sized.Rising),
                  budget);
  }

  auto Zeroed(const Piece &piece, const Carrier &carrier, std::size_t budget) -> Zeros
  {
    const Relative from  = Toward(piece, carrier.Through);
    const Vector   way   = carrier.Way;
    Series         value = from.Across;
    Series         bound = from.Across;
    for (std::size_t index = 0; index <= piece.Degree; ++index)
    {
      const double across    = from.Across.Values.at(index);
      const double rising    = from.Rising.Values.at(index);
      value.Values.at(index) = (across * way.Y) - (rising * way.X);
      bound.Values.at(index) = std::abs(across * way.Y) + std::abs(rising * way.X);
    }
    return Turned(value, bound, budget);
  }

  auto Zeroed(const Piece &piece, const Matrix &unit, std::size_t budget) -> Zeros
  {
    Series     across = {.Degree = piece.Degree, .Values = {}};
    Series     rising = across;
    Series     weight = across;
    const auto poles  = Poles(piece);
    for (std::size_t index = 0; index <= piece.Degree; ++index)
    {
      const Homogeneous &pole = poles[index];
      across.Values.at(index) =
          (unit.X.X * pole.X) + (unit.Y.X * pole.Y) + (unit.Origin.X * pole.W);
      rising.Values.at(index) =
          (unit.X.Y * pole.X) + (unit.Y.Y * pole.Y) + (unit.Origin.Y * pole.W);
      weight.Values.at(index) = pole.W;
    }
    const Series value = ((across * across) + (rising * rising)) - (weight * weight);
    const Series bound =
        (Sized(across) * Sized(across)) + (Sized(rising) * Sized(rising)) + (weight * weight);
    return Turned(value, bound, budget);
  }

  auto Enclosed(const Piece &piece, Point point) -> Enclosure
  {
    Series     across = {.Degree = piece.Degree, .Values = {}};
    Series     rising = across;
    const auto poles  = Poles(piece);
    for (std::size_t index = 0; index <= piece.Degree; ++index)
    {
      const Point pole        = Projected(poles[index]);
      across.Values.at(index) = pole.X - point.X;
      rising.Values.at(index) = pole.Y - point.Y;
    }
    const Series turning = (across * Slope(rising)) - (rising * Slope(across));
    return {.Area = Integral(turning) / 2.0,
            .Moment =
                Vector{.X = Integral(across * turning), .Y = Integral(rising * turning)} / 3.0};
  }

  auto Quarter(const Matrix &map, std::size_t index) -> Piece
  {
    const Point  start = Turned({.X = 1, .Y = 0}, index);
    const Point  bend  = Turned({.X = 1, .Y = 1}, index);
    const Point  end   = Turned({.X = 0, .Y = 1}, index);
    const double half  = std::sqrt(0.5);
    Piece        piece = Polynomial(std::array{map * start, map * bend, map * end});
    piece.Poles[1]     = {.X = piece.Poles[1].X * half, .Y = piece.Poles[1].Y * half, .W = half};
    return piece;
  }

  auto Quartered(std::size_t index, double along) -> double
  {
    const double quarter = Pi / 2.0;
    if (Exactly(along, 1.0))
    {
      return quarter * static_cast<double>(index + 1);
    }
    const double start = quarter * static_cast<double>(index);
    const double half  = std::sqrt(0.5);
    return start + (2.0 * Trig::Atan(along * half / (1.0 - along + (along * half))));
  }

  auto Across(Range span, double along) -> double
  {
    if (Exactly(along, 1.0))
    {
      return span.End;
    }
    return span.Start + (along * (span.End - span.Start));
  }

  auto Straight(Point start, Point end) -> Piece
  {
    Piece piece    = {.Degree = 1, .Poles = {}};
    piece.Poles[0] = {.X = start.X, .Y = start.Y, .W = 1.0};
    piece.Poles[1] = {.X = end.X, .Y = end.Y, .W = 1.0};
    return piece;
  }

  auto Polynomial(std::span<const Point> points) -> Piece
  {
    Piece                        piece = {.Degree = points.size() - 1, .Poles = {}};
    const std::span<Homogeneous> poles{piece.Poles};
    for (std::size_t index = 0; index < points.size(); ++index)
    {
      poles[index] = {.X = points[index].X, .Y = points[index].Y, .W = 1.0};
    }
    return piece;
  }

  auto Polynomial(const Bezier::Quadratic &quadratic) -> Piece
  {
    const std::array points = {quadratic.Start, quadratic.Control, quadratic.End};
    return Polynomial(points);
  }

  auto Polynomial(const Bezier::Cubic &cubic) -> Piece
  {
    const std::array points = {cubic.Start, cubic.Out, cubic.In, cubic.End};
    return Polynomial(points);
  }

  auto Section(const Spline &spline, std::size_t span) -> Piece
  {
    Piece piece = {.Degree = static_cast<std::size_t>(spline.Degree), .Poles = {}};
    if (piece.Degree < Lowly)
    {
      Sectional<Lowly>(spline, span, piece);
    }
    else
    {
      Sectional<MaxDegree + 1>(spline, span, piece);
    }
    return piece;
  }

  auto Derived(const Spline &spline, std::size_t span) -> Hodograph
  {
    const Sharing sharing = Shared(spline, span);
    const auto    paired  = [&spline, &sharing](Axis axis) -> Series
    {
      return Polar(sharing,
                   [&spline, &sharing, axis](std::size_t near, std::size_t far) -> double
                   {
                     return Coordinate(spline.Poles[sharing.First + far], axis)
                          - Coordinate(spline.Poles[sharing.First + near], axis);
                   });
    };
    const Series across = paired(Axis::X);
    const Series rising = paired(Axis::Y);
    return {.Degree = across.Degree, .Across = across.Values, .Rising = rising.Values};
  }

  auto Touches(const Spline &spline, std::size_t span, Point point) -> Solution
  {
    const Sharing sharing = Shared(spline, span);
    const auto    legs    = [&spline, &sharing, point](std::size_t near,
                                                       std::size_t far) -> std::array<double, 2>
    {
      const Vector away  = spline.Poles[sharing.First + near] - point;
      const Vector along = spline.Poles[sharing.First + far] - spline.Poles[sharing.First + near];
      return std::array<double, 2>{away.X * along.Y, away.Y * along.X};
    };
    const Series turning = Polar(sharing,
                                 [&legs](std::size_t near, std::size_t far) -> double
                                 {
                                   const std::array<double, 2> products = legs(near, far);
                                   return products[0] - products[1];
                                 });
    const Series bound   = Polar(sharing,
                                 [&legs](std::size_t near, std::size_t far) -> double
                                 {
                                 const std::array<double, 2> products = legs(near, far);
                                 return std::abs(products[0]) + std::abs(products[1]);
                                 });
    return Unless(turning, bound, Budget);
  }
} // namespace Geoxy::Kernel
