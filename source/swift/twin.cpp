/***************************************************************************************************
 *  @file      source/swift/twin.cpp
 *  @brief     Two quadratic or cubic Bézier curves met: their crossings isolated in doubles in a
 *             frame moved to the first curve's start, each proved alone in its box by Krawczyk's
 *             operator, polished by one double-double Newton step that Krawczyk's operator on a
 *             tiny box bounds, its point and places proved the doubles nearest the true ones; and
 *             about each the curves kept from anything the thorough way's tolerance would decide.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/twin.hpp"

#include "scratch.hpp"
#include "sure.hpp"
#include "swift/casteljau.hpp"
#include "swift/lane.hpp"
#include "swift/roots.hpp"
#include "swift/section.hpp"
#include "trig.hpp"
#include "twofold.hpp"

#include "geoxy/domain.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <utility>

namespace Geoxy::Swift
{
  namespace
  {
    using Sure::Grow;
    using Sure::Number;
    using Sure::Shrink;
    using Sure::Tiny;

    /**
     *  How many pairs of pieces the isolation halves before it gives up: a pair of random cubics
     *  takes a few dozen, the hardest it is sure of a few hundred, and two curves that touch or lie
     *  along each other would take them all.
     */
    constexpr int Spending = 1024;

    /**
     *  How many halvings one curve's piece may take: its ends and its middle stay doubles exactly,
     *  and a piece of 2^-32 of the curve is finer than any crossing the kernel is sure of needs.
     */
    constexpr int Finest = 32;

    /**
     *  From how many halvings on a pair is tried by Krawczyk's operator, however large.
     */
    constexpr int Trying = 6;

    /**
     *  From how many halvings on a pair left open is tried again in a box centred where Newton's
     *  method goes: a crossing on the edge between two pieces is proved in neither.
     */
    constexpr int Centring = 16;

    /**
     *  How many such boxes the isolation tries at the most: the first few settle a crossing on an
     *  edge or along a shallow strip; about a touch, every one fails.
     */
    constexpr int Recentring = 64;

    /**
     *  The most crossings the kernel is sure of: eight, one fewer than two cubics may have, so that
     *  nine places of the first curve evenly spaced never each lie by a crossing of their own.
     */
    constexpr std::size_t Most = 8;

    /**
     *  How far apart the isolation keeps what does not cross, in reaches: so far that about a
     *  crossing steeper than a degree or so it alone keeps the rest of the second curve from the
     *  stretches of the first the thorough way may weigh for lying along; a shallower one asks a
     *  ball about it.
     */
    constexpr double Wide = 64.0;

    /**
     *  Two numbers, one for each curve, the first's first: places, radii, sizes.
     */
    using Places = std::array<double, 2>;

    /**
     *  A quadratic or cubic Bézier curve as the kernel takes it: its control points as given and
     *  moved to the frame, its hodograph's, and bounds on its derivatives on its whole stretch.
     */
    struct Bent
    {
      std::array<Coefficients, 2> Given  = {};  /**< Its control points' x and y, as given.     */
      std::array<Coefficients, 2> Moved  = {};  /**< Them moved to the frame, each rounded.     */
      std::array<Coefficients, 2> Way    = {};  /**< The hodograph's, n (P_i+1 - P_i), rounded. */
      std::size_t                 Degree = 0;   /**< Its degree, two or three.                  */
      double                      Speed  = 0.0; /**< The most |B'| may be.                      */
      double                      Bend   = 0.0; /**< The most |B''| may be.                     */
      double                      Turn   = 0.0; /**< The most |B'''| may be.                    */
      double                      Slip   = 0.0; /**< The most a coordinate of B' worked out in
                                                     doubles is off.                           */
      const Section *Source = nullptr;          /**< The piece, as given. */
    };

    /**
     *  The two curves and what holds for both: the frame they are moved to, and how far what is
     *  worked out in it may be off.
     */
    struct Twin
    {
      std::array<Bent, 2> Curves = {};  /**< The curves, the first's first.               */
      double              Scale  = 0.0; /**< The largest size of a moved coordinate.      */
      double              Shift  = 0.0; /**< The most a moved coordinate is off: u Scale. */
      double              Miss   = 0.0; /**< The most a coordinate of the gap between the
                                             curves worked out in the frame is off.             */
      double Grain = 0.0;               /**< 64 u times the largest size of a given coordinate:
                                             the most the thorough way's own working of a point
                                             or a distance is off.                              */
    };

    /**
     *  A piece of a curve, in the frame: its control points, its stretch of the curve's parameter,
     *  and how far its control points may be from the true piece's.
     */
    struct Slice
    {
      std::array<Coefficients, 2> At    = {};  /**< Its control points' x and y.  */
      double                      Lo    = 0.0; /**< Where its stretch starts.     */
      double                      Hi    = 1.0; /**< Where it ends.                */
      double                      Off   = 0.0; /**< The most a coordinate is off. */
      int                         Depth = 0;   /**< How many halvings made it.    */
    };

    /**
     *  A piece of each curve, as the isolation holds them.
     */
    struct Couple
    {
      std::array<Slice, 2> Slices = {}; /**< The pieces, the first curve's first. */
      int                  Depth  = 0;  /**< How many halvings made them.         */
    };

    /**
     *  A box of the two places for Krawczyk's operator, and what it takes: the box, the
     *  derivatives' ranges on it, and the gap between the curves at its middle as worked out.
     */
    struct Trial
    {
      Places               Center = {}; /**< The box's middle.                           */
      Places               Radius = {}; /**< Its half widths.                            */
      std::array<Range, 4> Slopes = {}; /**< A'x, A'y, B'x and B'y on it.                */
      Places               Gap    = {}; /**< A - B at the middle, worked out.            */
      Places               Miss   = {}; /**< How far each coordinate of that may be off. */
    };

    /**
     *  What Krawczyk's operator found of a box.
     */
    struct Verdict
    {
      bool   One    = false; /**< Whether it proves one root in the box, only one. */
      Places Offset = {};    /**< -Y F at the middle, worked out: the Newton step. */
      Places Reach  = {};    /**< How far each way the operator's image runs from the
                                  middle moved by the offset: every root in the box
                                  lies in it.                                          */
      Places Unit = {};      /**< How far from the root a place of the box may lie
                                  where the curves come within one of each other.      */
    };

    /**
     *  A box a root is proved alone in.
     */
    struct Cage
    {
      Places Center = {}; /**< The box's middle.                                */
      Places Radius = {}; /**< Its half widths.                                 */
      Places Root   = {}; /**< The middle moved by the Newton step, worked out. */
      Places Spread = {}; /**< How far from that the root may lie, each way.    */
      Places Unit   = {}; /**< How far from the root a place of the box may lie, where the
                               curves come within one of each other.                     */
    };

    /**
     *  The boxes roots were proved in.
     */
    using Cages = List<Cage>;

    /**
     *  A root, and the boxes it was proved in.
     */
    struct Rootbox
    {
      List<std::size_t> Boxes = {}; /**< Which of the boxes found.          */
      Places            Unit  = {}; /**< The most of their bands, per unit. */
    };

    /**
     *  The roots.
     */
    using Rootboxes = List<Rootbox>;

    /**
     *  A crossing proved: its point and places, the doubles nearest the true ones, where its places
     *  lie, and how far the curves may come from each other about it.
     */
    struct Settled
    {
      Point  At     = {}; /**< The point.                        */
      Places Along  = {}; /**< Its places.                       */
      Places Center = {}; /**< The places, to within the spread. */
      Places Spread = {}; /**< How far each may be.              */
      Places Unit   = {}; /**< How far from them a place may lie where the curves come
                               within one of each other.                                 */
      Range Hull = {};    /**< The second curve's places its boxes span, and once weighed
                               the stretch of them its pieces about it lie in.            */
      std::array<Range, 2> Bands = {}; /**< Its band on each curve, where the curves come
                                            within the tolerance.                        */
      Places Runs  = {};               /**< How far each curve may run on its band. */
      double Reach = 0.0;              /**< R: how far a piece of the second the thorough
                                            way weighs against a stretch of the first's
                                            band may lie from that stretch's middle.     */
      double Margin = 0.0;             /**< Above zero where the crossing is steep and the
                                            second curve runs one way about it.          */
    };

    /**
     *  A disk the derivative of a curve keeps to on a stretch, its values and the Bernstein
     *  coefficients of every part of the stretch alike.
     */
    struct Disk
    {
      Vector Center = {};  /**< Its middle. */
      double Radius = 0.0; /**< Its radius. */
    };

    /**
     *  The lengths the checks about the crossings weigh.
     */
    struct Levels
    {
      double Length = 0.0; /**< The tolerance's. */
      double Within = 0.0; /**< The most the curves may be apart where the thorough way
                                finds them within the tolerance.                          */
      double Twice = 0.0;  /**< Likewise for twice the tolerance.           */
      double Grain = 0.0;  /**< The most its own working of a point is off. */
      double Sine  = 0.0;  /**< The sine of the tolerance's angle.          */
    };

    /**
     *  A curve as the kernel takes it. The hodograph's control points n (P_i+1 - P_i) are each
     *  within 2 n u D of the exact ones a coordinate, D the largest size of a step's coordinate:
     *  the step rounded once, then the product. Then |B'| is no more than the largest of theirs,
     *  |B''| than n - 1 times their largest step, and a cubic's |B'''| than twice their second
     *  difference, each with what those roundings add; and B' worked out from them in doubles is
     *  off by their own error and three steps of a double at their largest a round of de
     *  Casteljau's algorithm.
     *
     *  Control points within a slack s of the true ones a coordinate move each of the hodograph's
     *  by 2 n s a coordinate, so by sqrt 2 times that in length, its steps by twice that and its
     *  second differences by four times: each bound grows by as much, and the slip by 2 n s.
     *
     *  @param  section The piece, its control points as given.
     *  @param  origin  The frame's origin.
     *  @return         The curve.
     */
    auto Made(const Section &section, Point origin) -> Bent
    {
      const std::span<const Point> poles = std::span{section.Poles}.first(section.Degree + 1);
      Bent                         bent;
      bent.Degree      = section.Degree;
      bent.Source      = &section;
      const auto order = static_cast<double>(bent.Degree);
      for (std::size_t index = 0; index <= bent.Degree; ++index)
      {
        const Point pole        = poles[index];
        bent.Given[0].at(index) = pole.X;
        bent.Given[1].at(index) = pole.Y;
        bent.Moved[0].at(index) = pole.X - origin.X;
        bent.Moved[1].at(index) = pole.Y - origin.Y;
      }
      std::array<Vector, Highest> ways  = {};
      double                      step  = 0.0;
      double                      way   = 0.0;
      double                      speed = 0.0;
      for (std::size_t index = 0; index < bent.Degree; ++index)
      {
        const Vector difference = poles[index + 1] - poles[index];
        ways.at(index)          = difference * order;
        bent.Way[0].at(index)   = ways.at(index).X;
        bent.Way[1].at(index)   = ways.at(index).Y;
        step                    = std::max({step, std::abs(difference.X), std::abs(difference.Y)});
        way   = std::max({way, std::abs(ways.at(index).X), std::abs(ways.at(index).Y)});
        speed = std::max(speed, Length(ways.at(index)));
      }
      double bend = 0.0;
      for (std::size_t index = 0; index + 2 <= bent.Degree; ++index)
      {
        bend = std::max(bend, Length(ways.at(index + 1) - ways.at(index)));
      }
      const double off   = 2.0 * order * Unit * step;
      const double moved = 2.0 * order * section.Slack;
      const double apart = std::numbers::sqrt2 * moved * Grow;
      bent.Speed         = (speed + (2.0 * off) + apart) * Grow;
      bent.Bend = (order - 1.0) * (bend + (4.0 * off) + (4.0 * Unit * way) + (2.0 * apart)) * Grow;
      bent.Turn = bent.Degree == Highest ? 2.0
                                               * (Length((ways[2] - ways[1]) - (ways[1] - ways[0]))
                                                  + (8.0 * off)
                                                  + (16.0 * Unit * way)
                                                  + (4.0 * apart))
                                               * Grow
                                         : 0.0;
      bent.Slip = ((off + (3.0 * (order - 1.0) * Unit * (way + off)) + moved) * Grow) + Tiny;
      return bent;
    }

    /**
     *  How far a curve's control points moved to the frame may be from its true piece's: the move's
     *  rounding, and the slack of its control points as given.
     *
     *  @param  twin  The curves.
     *  @param  which Which curve.
     *  @return       The most a coordinate is off.
     */
    auto Offset(const Twin &twin, std::size_t which) -> double
    {
      return twin.Shift + twin.Curves.at(which).Source->Slack;
    }

    /**
     *  A curve's derivative at a place, in doubles, off by its slip a coordinate.
     *
     *  @param  bent  The curve.
     *  @param  along The place, from 0 to 1.
     *  @return       The derivative.
     */
    auto Slope(const Bent &bent, double along) -> Vector
    {
      return {.X = Valued(bent.Way[0], bent.Degree - 1, along),
              .Y = Valued(bent.Way[1], bent.Degree - 1, along)};
    }

    /**
     *  The gap A(s) - B(t) between the curves, in the frame, in doubles: off by the twin's miss a
     *  coordinate.
     *
     *  @param  twin  The curves.
     *  @param  place The places s and t, each from 0 to 1.
     *  @return       The gap.
     */
    auto Mismatch(const Twin &twin, Places place) -> Vector
    {
      const Bent &first  = twin.Curves[0];
      const Bent &second = twin.Curves[1];
      return {.X = Valued(first.Moved[0], first.Degree, place[0])
                 - Valued(second.Moved[0], second.Degree, place[1]),
              .Y = Valued(first.Moved[1], first.Degree, place[0])
                 - Valued(second.Moved[1], second.Degree, place[1])};
    }

    /**
     *  A range's middle and the most either end is from it.
     *
     *  @param  range The range.
     *  @return       The middle and the radius, which covers the range whatever the rounding.
     */
    auto Centred(Range range) -> std::array<double, 2>
    {
      const double middle = 0.5 * (range.Start + range.End);
      return {middle, std::max(range.End - middle, middle - range.Start) * Grow};
    }

    /**
     *  Krawczyk's operator on a box X about m, K = m - Y F(m) + (I - Y J(X)) (X - m), J(X) the
     *  Jacobian's ranges on the box and Y the inverse of their middles: K inside the box's inside
     *  proves the box holds one root and only one, and every Jacobian in J(X) invertible.
     *
     *  The offset -Y F(m) is worked out directly, not as m - Y F(m), whose rounding would swamp a
     *  box finer than a step of a double; it is off by |Y| times F's miss and four steps of a
     *  double of its terms. Each entry of |I - Y J| is no more than the entry at the middles,
     *  worked out and off by four steps of a double of 1 and its terms, and |Y| times the ranges'
     *  radii.
     *
     *  Where it proves one root, the band: for a place x of the box where |F(x)| <= l, F(x) = M (x
     *  - x*) for some M in J(X), its rows the gradients somewhere on the segment, so x - x* = (I -
     *    C)^-1 Y F(x) with C = I - Y M, |C| <= Ĉ; Ĉ r < r for the box's radii r, so Ĉ's spectral
     *    radius is below one and |x - x*| <= (I - Ĉ)^-1 |Y F| <= l (I - Ĉ)^-1 g, g each row of Y's
     *    length. The unit band is that for l = 1.
     *
     *  @param  trial The box and what it takes.
     *  @return       What the operator proves, its offset and reach, and the unit band.
     */
    auto Krawczyk(const Trial &trial) -> Verdict
    {
      Verdict      verdict;
      const auto   a   = Centred(trial.Slopes[0]);
      const auto   c   = Centred(trial.Slopes[1]);
      const auto   b   = Centred(trial.Slopes[2]);
      const auto   d   = Centred(trial.Slopes[3]);
      const double det = (b[0] * c[0]) - (a[0] * d[0]);
      /**
       *  A singular middle gives a Y past the doubles, whose image settles nothing: the least
       *  normal double of the determinant's sign keeps the division off zero.
       */
      const double safe =
          std::copysign(std::max(std::abs(det), std::numeric_limits<double>::min()), det);
      const double y11     = -d[0] / safe;
      const double y12     = b[0] / safe;
      const double y21     = -c[0] / safe;
      const double y22     = a[0] / safe;
      const double gapX    = trial.Gap[0];
      const double gapY    = trial.Gap[1];
      const double offsetS = -((y11 * gapX) + (y12 * gapY));
      const double offsetT = -((y21 * gapX) + (y22 * gapY));
      const double slopS   = (((std::abs(y11) * trial.Miss[0]) + (std::abs(y12) * trial.Miss[1]))
                              + (4.0 * Unit * (std::abs(y11 * gapX) + std::abs(y12 * gapY))))
                           * Grow;
      const double slopT   = (((std::abs(y21) * trial.Miss[0]) + (std::abs(y22) * trial.Miss[1]))
                              + (4.0 * Unit * (std::abs(y21 * gapX) + std::abs(y22 * gapY))))
                           * Grow;
      const auto   entry   = [](double middle, std::array<double, 2> terms, double spread) -> double
      {
        return (std::abs(middle + terms[0] + terms[1])
                + spread
                + (4.0 * Unit * (std::abs(middle) + std::abs(terms[0]) + std::abs(terms[1]))))
             * Grow;
      };
      const double c11 = entry(1.0, {-(y11 * a[0]), -(y12 * c[0])},
                               (std::abs(y11) * a[1]) + (std::abs(y12) * c[1]));
      const double c12 =
          entry(0.0, {y11 * b[0], y12 * d[0]}, (std::abs(y11) * b[1]) + (std::abs(y12) * d[1]));
      const double c21 = entry(0.0, {-(y21 * a[0]), -(y22 * c[0])},
                               (std::abs(y21) * a[1]) + (std::abs(y22) * c[1]));
      const double c22 =
          entry(1.0, {y21 * b[0], y22 * d[0]}, (std::abs(y21) * b[1]) + (std::abs(y22) * d[1]));
      const double reachS = (slopS + (c11 * trial.Radius[0]) + (c12 * trial.Radius[1])) * Grow;
      const double reachT = (slopT + (c21 * trial.Radius[0]) + (c22 * trial.Radius[1])) * Grow;
      /**
       *  Inside, Ĉ r < r makes (1 - c11) (1 - c22) > c12 c21; where rounding leaves it no more, the
       *  band runs past the doubles and settles nothing.
       */
      const double rest =
          std::max(((1.0 - c11) * (1.0 - c22) * Shrink * Shrink) - (c12 * c21 * Grow),
                   std::numeric_limits<double>::min());
      const double rowS = Trig::Hypot(y11, y12) * Grow;
      const double rowT = Trig::Hypot(y21, y22) * Grow;
      verdict.Offset    = {offsetS, offsetT};
      verdict.Reach     = {reachS, reachT};
      verdict.One       = (std::abs(offsetS) + reachS) * Grow < trial.Radius[0]
                       && (std::abs(offsetT) + reachT) * Grow < trial.Radius[1];
      verdict.Unit      = {((((1.0 - c22) * rowS) + (c12 * rowT)) * Grow / (rest * Shrink)) * Grow,
                           (((c21 * rowS) + ((1.0 - c11) * rowT)) * Grow / (rest * Shrink)) * Grow};
      return verdict;
    }

    /**
     *  A box a root is proved in, from the trial and the verdict: the root's place within its
     *  spread of the middle moved by the step, that sum's rounding added.
     *
     *  @param  trial   The box.
     *  @param  verdict What Krawczyk's operator proved.
     *  @return         The box found.
     */
    auto Proved(const Trial &trial, const Verdict &verdict) -> Cage
    {
      Cage found = {.Center = trial.Center,
                    .Radius = trial.Radius,
                    .Root   = {},
                    .Spread = {},
                    .Unit   = verdict.Unit};
      for (std::size_t which = 0; which < 2; ++which)
      {
        const double center  = trial.Center.at(which);
        const double offset  = verdict.Offset.at(which);
        found.Root.at(which) = center + offset;
        found.Spread.at(which) =
            (verdict.Reach.at(which) + (2.0 * Unit * (std::abs(center) + std::abs(offset)))) * Grow;
      }
      return found;
    }

    /**
     *  A piece halved by de Casteljau's algorithm at its middle: each new control point an average
     *  of averages, each rounded once, a step of a double at the frame's scale or less, so each
     *  half's control points are within the piece's error and the grain of the true half's.
     *
     *  @param  slice  The piece.
     *  @param  degree Its curve's degree.
     *  @param  grain  How far a halving may move a control point: n u times the frame's scale.
     *  @return        The two halves, the first's first.
     */
    auto Halved(const Slice &slice, std::size_t degree, double grain) -> std::array<Slice, 2>
    {
      std::array<Slice, 2> halves = {slice, slice};
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        Coefficients work                = slice.At.at(axis);
        halves[0].At.at(axis)[0]         = work[0];
        halves[1].At.at(axis).at(degree) = work.at(degree);
        for (std::size_t round = 1; round <= degree; ++round)
        {
          for (std::size_t index = 0; index + round <= degree; ++index)
          {
            work.at(index) = (work.at(index) + work.at(index + 1)) * 0.5;
          }
          halves[0].At.at(axis).at(round)          = work[0];
          halves[1].At.at(axis).at(degree - round) = work.at(degree - round);
        }
      }
      const double middle = slice.Lo + ((slice.Hi - slice.Lo) * 0.5);
      halves[0].Hi        = middle;
      halves[1].Lo        = middle;
      for (Slice &half : halves)
      {
        half.Off   = ((slice.Off + grain) * Grow) + Tiny;
        half.Depth = slice.Depth + 1;
      }
      return halves;
    }

    /**
     *  The box a piece's control points span: its least and most x, and its least and most y.
     *
     *  @param  slice  The piece.
     *  @param  degree Its curve's degree.
     *  @return        The box.
     */
    auto SliceBox(const Slice &slice, std::size_t degree) -> std::array<double, 4>
    {
      const auto abscissas = std::span{slice.At[0]}.first(degree + 1);
      const auto ordinates = std::span{slice.At[1]}.first(degree + 1);
      return {std::ranges::min(abscissas), std::ranges::max(abscissas), std::ranges::min(ordinates),
              std::ranges::max(ordinates)};
    }

    /**
     *  One piece's fat line against another's control points: whether the other lies wholly on one
     *  side of the band the first's control points span about its chord, farther than a length.
     *  With w the chord, every cross product w x (X - P0) worked out is off by four steps of a
     *  double of its terms, each coordinate of X - P0 under twice the frame's scale; the true
     *  pieces lie within their errors of their control points, which moves a cross product by at
     *  most (|wx| + |wy|) times the error; and a gap g in cross products is a gap g / |w| in
     *  length.
     *
     *  @param  pieces The first piece and the other, each with its degree.
     *  @param  scale  The frame's scale.
     *  @param  length The length.
     *  @return        Whether the other lies that far off.
     */
    auto Fat(std::array<std::pair<const Slice *, std::size_t>, 2> pieces, double scale,
             double length) -> bool
    {
      const auto [line, lineDegree]   = pieces[0];
      const auto [other, otherDegree] = pieces[1];
      const Vector start              = {.X = line->At[0][0], .Y = line->At[1][0]};
      const Vector way                = {.X = line->At[0].at(lineDegree) - start.X,
                                         .Y = line->At[1].at(lineDegree) - start.Y};
      const double chord              = Length(way);
      if (!(chord > 0.0))
      {
        return false;
      }
      const auto across = [&start, &way](const Slice &slice, std::size_t index) -> double
      {
        return Cross(way, Vector{.X = slice.At[0].at(index) - start.X,
                                 .Y = slice.At[1].at(index) - start.Y});
      };
      double low  = 0.0;
      double high = 0.0;
      for (std::size_t index = 1; index < lineDegree; ++index)
      {
        low  = std::min(low, across(*line, index));
        high = std::max(high, across(*line, index));
      }
      double otherLow  = std::numeric_limits<double>::infinity();
      double otherHigh = -std::numeric_limits<double>::infinity();
      for (std::size_t index = 0; index <= otherDegree; ++index)
      {
        otherLow  = std::min(otherLow, across(*other, index));
        otherHigh = std::max(otherHigh, across(*other, index));
      }
      const double sides = std::abs(way.X) + std::abs(way.Y);
      const double slack = ((length * chord * (1.0 + 0x1p-50))
                            + ((line->Off + other->Off) * sides)
                            + (16.0 * Unit * sides * scale * Grow))
                         * Grow;
      return otherLow - high > slack || low - otherHigh > slack;
    }

    /**
     *  Whether two pieces are certainly farther apart than a length: their boxes, padded by their
     *  errors, or either's fat line against the other.
     *
     *  @param  couple The pieces.
     *  @param  twin   The curves.
     *  @param  span   The boxes the pieces' control points span, and the length.
     *  @return        Whether they are.
     */
    auto Far(const Couple &couple, const Twin &twin,
             const std::pair<std::array<std::array<double, 4>, 2>, double> &span) -> bool
    {
      const auto &[boxes, length]   = span;
      const Slice      &one         = couple.Slices[0];
      const Slice      &other       = couple.Slices[1];
      const std::size_t oneDegree   = twin.Curves[0].Degree;
      const std::size_t otherDegree = twin.Curves[1].Degree;
      const auto       &box         = boxes[0];
      const auto       &beside      = boxes[1];
      const double      clear       = (length + one.Off + other.Off) * Grow;
      if (beside[0] - box[1]    > clear
          || box[0] - beside[1] > clear
          || beside[2] - box[3] > clear
          || box[2] - beside[3] > clear)
      {
        return true;
      }
      return Fat({std::pair{&one, oneDegree}, std::pair{&other, otherDegree}}, twin.Scale, length)
          || Fat({std::pair{&other, otherDegree}, std::pair{&one, oneDegree}}, twin.Scale, length);
    }

    /**
     *  The ranges a piece's derivative keeps to, in its curve's own parameter, x's and y's: the
     *  hull of n (P_i+1 - P_i) / w, w the stretch's width, a power of two; each difference rounded
     *  once and the product once, and the true piece's control points within its error of these,
     *  which moves a coefficient by 2 n times the error over w.
     *
     *  @param  slice  The piece.
     *  @param  degree Its curve's degree.
     *  @return        The ranges.
     */
    auto Speeds(const Slice &slice, std::size_t degree) -> std::array<Range, 2>
    {
      const double         width = slice.Hi - slice.Lo;
      const auto           order = static_cast<double>(degree);
      std::array<Range, 2> out   = {};
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        double low     = std::numeric_limits<double>::infinity();
        double high    = -std::numeric_limits<double>::infinity();
        double largest = 0.0;
        for (std::size_t index = 0; index < degree; ++index)
        {
          const double step  = slice.At.at(axis).at(index + 1) - slice.At.at(axis).at(index);
          const double value = order * step / width;
          low                = std::min(low, value);
          high               = std::max(high, value);
          largest            = std::max(largest, std::abs(step));
        }
        const double error = ((2.0 * order * (slice.Off + (Unit * largest))) / width) * Grow;
        const double room = (error + (2.0 * Unit * std::max(std::abs(low), std::abs(high)))) * Grow;
        out.at(axis)      = {.Start = low - room, .End = high + room};
      }
      return out;
    }

    /**
     *  The curves' derivatives' ranges on a box about a place, from their derivatives there and
     *  their bends: within the bend times the radius of the derivative at the middle, which is off
     *  by the slip.
     *
     *  @param  twin   The curves.
     *  @param  center The box's middle.
     *  @param  radius Its half widths.
     *  @return        A'x, A'y, B'x and B'y's ranges.
     */
    auto Around(const Twin &twin, Places center, Places radius) -> std::array<Range, 4>
    {
      std::array<Range, 4> out = {};
      for (std::size_t which = 0; which < 2; ++which)
      {
        const Bent  &bent       = twin.Curves.at(which);
        const Vector slope      = Slope(bent, center.at(which));
        const double room       = ((bent.Bend * radius.at(which)) + bent.Slip) * Grow;
        const double sizes      = 2.0 * Unit * std::max(std::abs(slope.X), std::abs(slope.Y));
        const double wider      = (room + sizes) * Grow;
        out.at(2 * which)       = {.Start = slope.X - wider, .End = slope.X + wider};
        out.at((2 * which) + 1) = {.Start = slope.Y - wider, .End = slope.Y + wider};
      }
      return out;
    }

    /**
     *  A few steps of Newton's method on the gap between the curves, in doubles, in the frame;
     *  stopped where a place leaves its curve's stretch, the Jacobian is singular, or the steps are
     *  lost.
     *
     *  @param  twin  The curves.
     *  @param  place Where to start.
     *  @param  steps How many steps at the most.
     *  @return       Where it went.
     */
    auto Approached(const Twin &twin, Places place, int steps) -> Places
    {
      for (int step = 0; step < steps; ++step)
      {
        if (std::max({-place[0], place[0] - 1.0, -place[1], place[1] - 1.0}) > 0.0)
        {
          break;
        }
        const Vector gap   = Mismatch(twin, place);
        const Vector one   = Slope(twin.Curves[0], place[0]);
        const Vector other = Slope(twin.Curves[1], place[1]);
        const double det   = (other.X * one.Y) - (one.X * other.Y);
        if (!(std::abs(det) > 0.0))
        {
          break;
        }
        const double moveS = ((other.X * gap.Y) - (other.Y * gap.X)) / det;
        const double moveT = ((one.X * gap.Y) - (one.Y * gap.X)) / det;
        place              = {place[0] - moveS, place[1] - moveT};
        if (std::abs(moveS) <= 0x1p-52 && std::abs(moveT) <= 0x1p-52)
        {
          break;
        }
      }
      return place;
    }

    /**
     *  A pair of pieces tried by Krawczyk's operator on their stretches, the derivatives' ranges
     *  their hulls'.
     *
     *  @param  twin   The curves.
     *  @param  couple The pieces.
     *  @return        The trial.
     */
    auto Tried(const Twin &twin, const Couple &couple) -> Trial
    {
      const Slice &one    = couple.Slices[0];
      const Slice &other  = couple.Slices[1];
      const Places center = {one.Lo + ((one.Hi - one.Lo) * 0.5),
                             other.Lo + ((other.Hi - other.Lo) * 0.5)};
      const auto   first  = Speeds(one, twin.Curves[0].Degree);
      const auto   second = Speeds(other, twin.Curves[1].Degree);
      const Vector gap    = Mismatch(twin, center);
      return {.Center = center,
              .Radius = {(one.Hi - one.Lo) * 0.5, (other.Hi - other.Lo) * 0.5},
              .Slopes = {first[0], first[1], second[0], second[1]},
              .Gap    = {gap.X, gap.Y},
              .Miss   = {twin.Miss, twin.Miss}};
    }

    /**
     *  What a box centred where Newton's method goes found.
     */
    struct Aimed
    {
      std::optional<Cage> Box      = {};    /**< The box proved, where it is. */
      bool                Touching = false; /**< Whether the curves run within 2^-20 of
                                                 parallel where Newton's method went.   */
    };

    /**
     *  A pair of pieces the operator left unproved tried again in a box centred where Newton's
     *  method goes from the operator's image, wide enough to hold the pieces' box whole and as wide
     *  again as the pieces: a root on the edge between two pieces is proved in neither, and about a
     *  shallow crossing the curves stay near each other along a strip of pairs, each holding no
     *  root, that only a box about the root holds. Where Newton's method ends with the curves
     *  within 2^-20 of parallel, they touch there, or cross too shallowly for the kernel ever to be
     *  sure, and the isolation gives up.
     *
     *  @param  twin    The curves.
     *  @param  cell    The pieces' box, as tried.
     *  @param  verdict What the operator found of it.
     *  @return         The box proved, where it is, and whether the curves touch.
     */
    auto Recentred(const Twin &twin, const Trial &cell, const Verdict &verdict) -> Aimed
    {
      const Places start = {cell.Center[0] + verdict.Offset[0], cell.Center[1] + verdict.Offset[1]};
      const Places place = Approached(twin, start, 3);
      Places       radius = {};
      for (std::size_t which = 0; which < 2; ++which)
      {
        const double spot = place.at(which);
        radius.at(which) =
            (std::abs(spot - cell.Center.at(which)) * Grow) + (2.0 * cell.Radius.at(which));
        if (!(spot >= radius.at(which)) || !(spot <= 1.0 - radius.at(which)))
        {
          return {};
        }
      }
      const Vector one   = Slope(twin.Curves[0], place[0]);
      const Vector other = Slope(twin.Curves[1], place[1]);
      if (std::abs(Cross(one, other)) < 0x1p-20 * Length(one) * Length(other))
      {
        return {.Box = std::nullopt, .Touching = true};
      }
      const Vector  gap   = Mismatch(twin, place);
      const Trial   trial = {.Center = place,
                             .Radius = radius,
                             .Slopes = Around(twin, place, radius),
                             .Gap    = {gap.X, gap.Y},
                             .Miss   = {twin.Miss, twin.Miss}};
      const Verdict proof = Krawczyk(trial);
      return {.Box      = proof.One ? std::optional<Cage>{Proved(trial, proof)} : std::nullopt,
              .Touching = false};
    }

    /**
     *  A pair of pieces tried by Krawczyk's operator on their box, and where it is left open from
     *  enough halvings on, while tries are left, in a box centred where Newton's method goes.
     *
     *  @param  twin   The curves.
     *  @param  couple The pieces.
     *  @param  tries  How many centred boxes are left to try, one taken where one is tried.
     *  @return        The box proved, where it is, and whether the curves touch.
     */
    auto Proving(const Twin &twin, const Couple &couple, int &tries) -> Aimed
    {
      const Trial   trial   = Tried(twin, couple);
      const Verdict verdict = Krawczyk(trial);
      if (verdict.One)
      {
        return {.Box = Proved(trial, verdict), .Touching = false};
      }
      return couple.Depth >= Centring && tries-- > 0 ? Recentred(twin, trial, verdict) : Aimed{};
    }

    /**
     *  Whether a pair of pieces lies wholly in a box already proved: its stretches within the
     *  box's.
     *
     *  @param  couple The pieces.
     *  @param  founds The boxes.
     *  @return        Whether it does.
     */
    auto Covered(const Couple &couple, const Cages &founds) -> bool
    {
      const Slice &one    = couple.Slices[0];
      const Slice &other  = couple.Slices[1];
      const Places center = {one.Lo + ((one.Hi - one.Lo) * 0.5),
                             other.Lo + ((other.Hi - other.Lo) * 0.5)};
      const Places radius = {(one.Hi - one.Lo) * 0.5, (other.Hi - other.Lo) * 0.5};
      return std::ranges::any_of(
          founds,
          [&center, &radius](const Cage &found) -> bool
          {
            return (std::abs(center[0] - found.Center[0]) + radius[0]) * Grow <= found.Radius[0]
                && (std::abs(center[1] - found.Center[1]) + radius[1]) * Grow <= found.Radius[1];
          });
    }

    /**
     *  The roots of the gap between the curves, each proved alone in a box by Krawczyk's operator:
     *  pairs of pieces halved, the larger first, dropped where certainly farther apart than the
     *  length, tried by the operator once small or deep enough. Every place of the two curves'
     *  stretches lies in a pair dropped or in a box proved, so wherever the curves come within the
     *  length of each other, it is in a box proved.
     *
     *  @param  twin   The curves.
     *  @param  length How far apart what is dropped must be.
     *  @return        The boxes; nothing where the budget runs out, a piece grows too fine, or
     *                 Newton's method ends where the curves touch.
     */
    auto BoxedRoots(const Twin &twin, double length) -> std::optional<Cages>
    {
      const std::array<double, 2> grains = {
          static_cast<double>(twin.Curves[0].Degree) * Unit * twin.Scale * (1.0 + 0x1p-30),
          static_cast<double>(twin.Curves[1].Degree) * Unit * twin.Scale * (1.0 + 0x1p-30)};
      List<Couple> stack;
      stack.push_back(Couple{.Slices = {Slice{.At    = twin.Curves[0].Moved,
                                              .Lo    = 0.0,
                                              .Hi    = 1.0,
                                              .Off   = Offset(twin, 0),
                                              .Depth = 0},
                                        Slice{.At    = twin.Curves[1].Moved,
                                              .Lo    = 0.0,
                                              .Hi    = 1.0,
                                              .Off   = Offset(twin, 1),
                                              .Depth = 0}},
                             .Depth  = 0});
      Cages founds;
      int   budget = Spending;
      int   tries  = Recentring;
      while (!stack.empty())
      {
        const Couple couple = stack.back();
        stack.pop_back();
        const std::array<std::array<double, 4>, 2> boxes = {
            SliceBox(couple.Slices[0], twin.Curves[0].Degree),
            SliceBox(couple.Slices[1], twin.Curves[1].Degree)};
        if (Far(couple, twin, {boxes, length}) || Covered(couple, founds))
        {
          continue;
        }
        const auto  &one   = boxes[0];
        const auto  &other = boxes[1];
        const Places spans = {std::max(one[1] - one[0], one[3] - one[2]),
                              std::max(other[1] - other[0], other[3] - other[2])};
        if (std::max(spans[0], spans[1]) < 0.25 * twin.Scale || couple.Depth > Trying)
        {
          const Aimed aimed = Proving(twin, couple, tries);
          if (aimed.Touching)
          {
            return std::nullopt;
          }
          if (aimed.Box)
          {
            founds.push_back(*aimed.Box);
            continue;
          }
        }
        const std::size_t which = spans[0] >= spans[1] ? 0 : 1;
        if (std::min(budget--, Finest - couple.Slices.at(which).Depth) <= 0)
        {
          return std::nullopt;
        }
        const auto halves =
            Halved(couple.Slices.at(which), twin.Curves.at(which).Degree, grains.at(which));
        for (std::size_t half = 2; half-- > 0;)
        {
          Couple next           = couple;
          next.Slices.at(which) = halves.at(half);
          next.Depth            = couple.Depth + 1;
          stack.push_back(next);
        }
      }
      return founds;
    }

    /**
     *  Whether the root one box holds certainly lies in another: its spread about its place inside
     *  the other's inside.
     *
     *  @param  one   The box whose root.
     *  @param  other The other box.
     *  @return       Whether it does.
     */
    auto Holds(const Cage &one, const Cage &other) -> bool
    {
      return (std::abs(one.Root[0] - other.Center[0]) + one.Spread[0]) * Grow < other.Radius[0]
          && (std::abs(one.Root[1] - other.Center[1]) + one.Spread[1]) * Grow < other.Radius[1];
    }

    /**
     *  The roots the boxes hold: a box whose root lies in one before it, or that holds the root of
     *  one before it, holds that one's root, as each holds one only; else its root is taken as one
     *  not seen before. Two boxes of one root not told so give two crossings at one point, which
     *  are never far enough apart to be sure of.
     *
     *  @param  founds The boxes.
     *  @return        The roots.
     */
    auto Coalesced(const Cages &founds) -> Rootboxes
    {
      List<std::size_t> label;
      Rootboxes         roots;
      for (std::size_t index = 0; index < founds.size(); ++index)
      {
        const Cage                &found = founds[index];
        std::optional<std::size_t> same;
        for (std::size_t before = 0; before < index; ++before)
        {
          if (Holds(found, founds[before]) || Holds(founds[before], found))
          {
            same = label[before];
          }
        }
        if (!same)
        {
          same = roots.size();
          roots.push_back({});
        }
        label.push_back(*same);
        Rootbox &root = roots[*same];
        root.Boxes.push_back(index);
        root.Unit = {std::max(root.Unit[0], found.Unit[0]), std::max(root.Unit[1], found.Unit[1])};
      }
      return roots;
    }

    /**
     *  A stretch about a place, rounded outward: every value within the radius of the center lies
     *  in it.
     *
     *  @param  center The center, from about 0 to about 1.
     *  @param  radius The radius.
     *  @return        The stretch.
     */
    auto About(double center, double radius) -> Range
    {
      const double reach = (radius * Grow) + 0x1p-52;
      return {.Start = center - reach, .End = center + reach};
    }

    /**
     *  A root proved and placed: polished by Newton's method in doubles, set on the grid of 2^-53
     *  so that 1 - s and 1 - t are doubles too, the gap there worked out in double-double from the
     *  curves' own numbers; a tiny box about it, twice the step and what the gap's bound may move
     *  it, proved to hold one root by Krawczyk's operator, its derivatives' ranges from the
     *  derivatives there and the bends; that box inside one the root was proved alone in, so its
     *  root is this one. The operator's image then holds the root: it is the step's place, within
     *  the image's reach, and the point by Taylor from the grid place on the first curve, A(s0) +
     *  A'(s0) D + A''(xi) D^2 / 2 for D = s* - s0, which the step is within the reach of; A'(s0)
     *  worked out is off by the slip. Each of the four is then the double nearest its value where
     *  the certificate settles it.
     *
     *  @param  twin   The curves.
     *  @param  root   The root.
     *  @param  founds The boxes it was proved in, among others.
     *  @return        The crossing; nothing where it is not settled.
     */
    auto Certified(const Twin &twin, const Rootbox &root, const Cages &founds)
        -> std::optional<Settled>
    {
      const Bent              &first  = twin.Curves[0];
      const Bent              &second = twin.Curves[1];
      const List<std::size_t> &boxes  = root.Boxes;
      const Cage              &cell   = founds.at(boxes.front());
      const Places             near   = Approached(twin, cell.Root, 8);
      const Places grid = {std::clamp(std::round(near[0] * 0x1p53) * 0x1p-53, 0.0, 1.0),
                           std::clamp(std::round(near[1] * 0x1p53) * 0x1p-53, 0.0, 1.0)};
      const std::array<Number, 2> here  = Exactly(*first.Source, grid[0]);
      const std::array<Number, 2> there = Exactly(*second.Source, grid[1]);
      const std::array<Number, 2> gap   = {here[0] - there[0], here[1] - there[1]};
      const Places                miss  = {(gap[0].Bound + std::abs(gap[0].Value.Low)) * Grow,
                                           (gap[1].Bound + std::abs(gap[1].Value.Low)) * Grow};
      const Places                value = {gap[0].Value.High, gap[1].Value.High};
      const Vector                one   = Slope(first, grid[0]);
      const Vector                other = Slope(second, grid[1]);
      const double                det   = (other.X * one.Y) - (one.X * other.Y);
      const double                safe =
          std::copysign(std::max(std::abs(det), std::numeric_limits<double>::min()), det);
      const double  y11    = -other.Y / safe;
      const double  y12    = other.X / safe;
      const double  y21    = -one.Y / safe;
      const double  y22    = one.X / safe;
      const double  most   = std::max(miss[0], miss[1]);
      const Places  radius = {(4.0 * std::abs((y11 * value[0]) + (y12 * value[1])))
                                  + (8.0 * (std::abs(y11) + std::abs(y12)) * most)
                                  + 0x1p-70,
                              (4.0 * std::abs((y21 * value[0]) + (y22 * value[1])))
                                  + (8.0 * (std::abs(y21) + std::abs(y22)) * most)
                                  + 0x1p-70};
      const Trial   trial  = {.Center = grid,
                              .Radius = radius,
                              .Slopes = Around(twin, grid, radius),
                              .Gap    = value,
                              .Miss   = miss};
      const Verdict proof  = Krawczyk(trial);
      const Cage    tiny   = Proved(trial, proof);
      double        out    = std::numeric_limits<double>::infinity();
      for (const std::size_t box : boxes)
      {
        const Cage &around = founds.at(box);
        out                = std::min(
            out, std::max(((std::abs(tiny.Center[0] - around.Center[0]) + tiny.Radius[0]) * Grow)
                              - around.Radius[0],
                          ((std::abs(tiny.Center[1] - around.Center[1]) + tiny.Radius[1]) * Grow)
                              - around.Radius[1]));
      }
      /**
       *  A tiny box not proved to hold a root, or not inside a box of this one, leaves every bound
       *  past the doubles: the certificate settles nothing.
       */
      const double lost =
          std::max(static_cast<double>(!proof.One), static_cast<double>(!(out < 0.0))) * 0x1p1000;
      const auto placed = [&twin, &grid, &proof, lost](std::size_t which) -> Number
      {
        const Number along = {.Value = Twofold::Add(Twofold::Lifted(grid.at(which)),
                                                    Twofold::Lifted(proof.Offset.at(which))),
                              .Bound = proof.Reach.at(which) + lost};
        return Spread(*twin.Curves.at(which).Source, along);
      };
      const double move  = proof.Offset[0];
      const double reach = proof.Reach[0];
      const double tail =
          (((first.Speed * reach)
            + (first.Slip * std::abs(move))
            + (first.Bend * (std::abs(move) + reach) * (std::abs(move) + reach) * 0.5))
           * Grow)
          + lost;
      const auto pointed = [&here, move, tail](std::size_t which, double slope) -> Number
      {
        const Number lean = Sure::Product(slope, move);
        return here.at(which) + Number{.Value = lean.Value, .Bound = lean.Bound + tail};
      };
      const auto rounded =
          Sure::Rounded(std::array{pointed(0, one.X), pointed(1, one.Y), placed(0), placed(1)});
      if (!rounded)
      {
        return std::nullopt;
      }
      const auto [pointX, pointY, placeS, placeT] = *rounded;
      Range hull                                  = {.Start = 1.0, .End = 0.0};
      for (const std::size_t box : boxes)
      {
        const Range around = About(founds.at(box).Center[1], founds.at(box).Radius[1]);
        hull = {.Start = std::min(hull.Start, around.Start), .End = std::max(hull.End, around.End)};
      }
      return Settled{.At     = {.X = pointX, .Y = pointY},
                     .Along  = {placeS, placeT},
                     .Center = tiny.Root,
                     .Spread = tiny.Spread,
                     .Unit   = root.Unit,
                     .Hull   = hull};
    }

    /**
     *  The stretch of a curve's places about a crossing where the curves may come within a length
     *  of each other: every place of the crossing's boxes where they do lies in it.
     *
     *  @param  settled The crossing.
     *  @param  length  The length, no more than the isolation kept apart.
     *  @param  which   Which curve's.
     *  @return         The stretch.
     */
    auto Band(const Settled &settled, double length, std::size_t which) -> Range
    {
      return About(settled.Center.at(which),
                   settled.Spread.at(which) + (length * settled.Unit.at(which)));
    }

    /**
     *  A stretch cut to a curve's.
     *
     *  @param  range The stretch.
     *  @return       It within 0 to 1.
     */
    auto CutRange(Range range) -> Range
    {
      return {.Start = std::max(range.Start, 0.0), .End = std::min(range.End, 1.0)};
    }

    /**
     *  How far apart two stretches are.
     *
     *  @param  one   One.
     *  @param  other The other.
     *  @return       The gap between them: above zero where they do not meet.
     */
    auto Gulf(Range one, Range other) -> double
    {
      return std::max(other.Start - one.End, one.Start - other.End);
    }

    /**
     *  The disk a curve's derivative keeps to on a stretch within its own: about the derivative at
     *  the middle, which is off by the slip, the bend times the half width h and a cubic's turn
     *  times h^2 / 2. A value there is within the bend times h of the middle's; a Bernstein
     *  coefficient of a part about a place m of half width h' is the derivative at an end of it,
     *  or, for a cubic, at m less h'^2 B''' / 2, so within as much of the middle's.
     *
     *  @param  bent    The curve.
     *  @param  stretch The stretch, within 0 to 1.
     *  @return         The disk.
     */
    auto Hodograph(const Bent &bent, Range stretch) -> Disk
    {
      const double middle = stretch.Start + ((stretch.End - stretch.Start) / 2.0);
      const double half   = (((stretch.End - stretch.Start) / 2.0) * Grow) + 0x1p-52;
      return {.Center = Slope(bent, middle),
              .Radius = ((bent.Bend * half) + (bent.Turn * half * half * 0.5) + (2.0 * bent.Slip))
                      * Grow};
    }

    /**
     *  How far a curve may run on a stretch: its width times the most its derivative's disk
     *  reaches.
     *
     *  @param  disk    The disk.
     *  @param  stretch The stretch.
     *  @return         The length, at the most.
     */
    auto Run(const Disk &disk, Range stretch) -> double
    {
      return (Length(disk.Center) + disk.Radius) * (stretch.End - stretch.Start) * Grow * Grow;
    }

    /**
     *  The least sine of the angle between a vector of one disk and one of another, every pair
     *  turning the same way: |a0 x b0| less |a0| rb + |b0| ra + ra rb and the cross product's
     *  rounding, over (|a0| + ra) (|b0| + rb).
     *
     *  @param  one   One disk.
     *  @param  other The other.
     *  @return       The sine, at the least; zero where the pairs may turn both ways, not a number
     *                where neither disk keeps from the origin.
     */
    auto LeastSine(const Disk &one, const Disk &other) -> double
    {
      const Vector a     = one.Center;
      const Vector b     = other.Center;
      const double sizeA = Length(a);
      const double sizeB = Length(b);
      const double terms = std::abs(a.X * b.Y) + std::abs(a.Y * b.X);
      const double least =
          (std::abs(Cross(a, b))
           - (((sizeA * other.Radius) + (sizeB * one.Radius) + (one.Radius * other.Radius)) * Grow)
           - (4.0 * Unit * terms))
          * Shrink;
      const double most = (sizeA + one.Radius) * (sizeB + other.Radius) * Grow;
      return std::max(least, 0.0) / std::max(most, std::numeric_limits<double>::min()) * Shrink;
    }

    /**
     *  A crossing weighed: its bands on both curves, where the curves come within the tolerance,
     *  how far each curve runs on its band, and how steep the crossing is about them.
     *
     *  The thorough way may weigh a stretch of the first curve within its band against a piece of
     *  the second for lying along it. Such a stretch, of a curve halved and never whole, has a
     *  parent at most twice as long and at most the band's width beyond it each way, which the
     *  thorough way halved as the larger of its pair; so the piece spans no more than the parent,
     *  and has the stretch's middle m within the tolerance of its box. Every point of the piece
     *  then lies within R = sqrt 2 (tolerance + the parent's span) of A(m). Where no point of the
     *  second outside the hull of the crossing's boxes, cut to its band of R, comes that near,
     *  which the isolation shows where R is no more than what it kept apart and no other crossing's
     *  band of R meets this one's, the piece lies in that hull; the second curve runs there within
     *  its derivative's disk, of radius under its middle's length, so its points as far as places
     *  apart run apart along that middle at its length less the radius; and some place of the band
     *  on the second is within the tolerance of A(m). So the piece lies within (R + tolerance) /
     *  (that length less the radius) of the band.
     *
     *  The crossing is steep where the derivatives' disks on the first's band and on that stretch
     *  of the second turn apart by more than twice the tolerance's angle: then every hit the
     *  thorough way finds in the bands passes the same way, and every pair of control steps of the
     *  stretch and the piece parts, so it does not weigh them.
     *
     *  @param  twin     The curves.
     *  @param  settled  The crossing, its places proved.
     *  @param  measures The lengths, and how far apart the isolation kept what does not cross.
     *  @return          The crossing, weighed.
     */
    auto Weigh(const Twin &twin, Settled settled, const std::pair<Levels, double> &measures)
        -> Settled
    {
      const auto &[levels, kept] = measures;
      const Bent  &first         = twin.Curves[0];
      const Bent  &second        = twin.Curves[1];
      const Range  bandS         = Band(settled, levels.Within, 0);
      const Range  bandT         = Band(settled, levels.Within, 1);
      const Disk   disk          = Hodograph(first, CutRange(bandS));
      const double width         = bandS.End - bandS.Start;
      const Disk   wide =
          Hodograph(first, CutRange({.Start = bandS.Start - width, .End = bandS.End + width}));
      const double speed = (Length(wide.Center) + wide.Radius) * Grow;
      const double reach =
          ((std::numbers::sqrt2 * Grow * (levels.Length + (2.0 * width * speed * Grow)))
           + levels.Grain)
          * Grow;
      const Range  near   = Band(settled, std::min(reach, kept), 1);
      const Range  hull   = CutRange({.Start = std::max(near.Start, settled.Hull.Start),
                                      .End   = std::min(near.End, settled.Hull.End)});
      const Disk   around = Hodograph(second, hull);
      const double slack  = (Length(around.Center) - around.Radius) * Shrink;
      const double slide =
          ((reach + levels.Within) / std::max(slack, std::numeric_limits<double>::min())) * Grow;
      const Disk piece =
          Hodograph(second, CutRange({.Start = bandT.Start - slide, .End = bandT.End + slide}));
      settled.Bands  = {bandS, bandT};
      settled.Runs   = {Run(disk, CutRange(bandS)),
                        Run(Hodograph(second, CutRange(bandT)), CutRange(bandT))};
      settled.Reach  = reach;
      settled.Hull   = hull;
      settled.Margin = std::min(LeastSine(disk, piece) - ((2.0 * levels.Sine) + 0x1p-30), slack);
      return settled;
    }

    /**
     *  Whether a crossing is clear of everything about it the thorough way would decide: steep, and
     *  the second curve running one way, as weighed; its bands within both curves' stretches, so
     *  that no end of either is within the tolerance of the other; and every end of either curve
     *  farther than the tolerance from any point of its band.
     *
     *  @param  twin    The curves.
     *  @param  settled The crossing, weighed.
     *  @param  levels  The lengths.
     *  @return         Whether it is.
     */
    auto Clear(const Twin &twin, const Settled &settled, const Levels &levels) -> bool
    {
      const auto &[bandS, bandT] = settled.Bands;
      const double inside = std::min({bandS.Start, 1.0 - bandS.End, bandT.Start, 1.0 - bandT.End});
      double       ends   = std::numeric_limits<double>::infinity();
      for (std::size_t which = 0; which < 2; ++which)
      {
        const Bent  &bent = twin.Curves.at(which);
        const double far  = (levels.Within + settled.Runs.at(which) + levels.Grain) * Grow;
        for (const std::size_t end : {std::size_t{0}, bent.Degree})
        {
          const Point pole = {.X = bent.Given[0].at(end), .Y = bent.Given[1].at(end)};
          ends             = std::min(ends, Distance(settled.At, pole) - far);
        }
      }
      return std::min({settled.Margin, inside, ends}) > 0.0;
    }

    /**
     *  Whether every two crossings stay two hits, as the thorough way joins them, and the bands
     *  keep the first curve from lying along the second: no more crossings than the kernel is sure
     *  of; each two points farther apart than the tolerance, each hit's point being within its
     *  band's run of its crossing; the first curve halfway between them, the short way round as
     *  well, in no band of twice the tolerance, so farther than that from the second; their bands
     *  on the first apart by more than any band is wide, so that nine places evenly spaced lie
     *  either in one band or each in one of its own; and, where the reach R that weighing either
     *  found is no more than the isolation kept apart, its band of R not meeting the other's band.
     *
     *  @param  settled  The crossings, weighed.
     *  @param  measures The lengths, and how far apart the isolation kept what does not cross.
     *  @return          Whether they do.
     */
    auto Unmerged(std::span<const Settled> settled, const std::pair<Levels, double> &measures)
        -> bool
    {
      const auto &[levels, kept] = measures;
      const double many   = static_cast<double>(Most) + 0.5 - static_cast<double>(settled.size());
      double       widest = 0.0;
      for (const Settled &one : settled)
      {
        widest = std::max(widest, one.Bands[0].End - one.Bands[0].Start);
      }
      for (std::size_t index = 0; index < settled.size(); ++index)
      {
        for (std::size_t other = index + 1; other < settled.size(); ++other)
        {
          const Settled &one       = settled[index];
          const Settled &another   = settled[other];
          const Range    bandOne   = one.Bands[0];
          const Range    bandOther = another.Bands[0];
          const double   far       = (levels.Within
                                      + std::max(one.Runs[0], one.Runs[1])
                                      + std::max(another.Runs[0], another.Runs[1])
                                      + levels.Grain)
                                   * Grow;
          const Range    half      = {.Start = ((bandOne.Start + bandOther.Start) / 2.0) - 0x1p-52,
                                      .End   = ((bandOne.End + bandOther.End) / 2.0) + 0x1p-52};
          double         lying     = std::numeric_limits<double>::infinity();
          for (const Settled &any : settled)
          {
            const Range twice = Band(any, levels.Twice, 0);
            for (const double shift : {-0.5, 0.0, 0.5})
            {
              lying = std::min(lying,
                               Gulf({.Start = half.Start + shift, .End = half.End + shift}, twice));
            }
          }
          const auto covered = [kept](const Settled &any) -> double
          {
            return std::min(any.Reach, kept) * static_cast<double>(any.Reach <= kept);
          };
          const double crossed = std::min(Gulf(bandOne, Band(another, covered(another), 0)),
                                          Gulf(bandOther, Band(one, covered(one), 0)));
          if (!(std::min({Distance(one.At, another.At) - far,
                          (Gulf(bandOne, bandOther) * Shrink) - (widest * Grow), crossed, lying,
                          many})
                > 0.0))
          {
            return false;
          }
        }
      }
      return true;
    }

    /**
     *  Whether no point of the second curve outside a stretch of its places comes within a ball:
     *  pieces of it halved, dropped where their box lies wholly outside the ball or their places
     *  within the stretch.
     *
     *  @param  twin The curves.
     *  @param  ball The ball's center as given and its radius.
     *  @param  hull The stretch.
     *  @return      Whether none does; not where it is not shown in a few hundred halvings.
     */
    auto Alone(const Twin &twin, const std::pair<Point, double> &ball, Range hull) -> bool
    {
      const Bent  &second = twin.Curves[1];
      const Point  origin = {.X = twin.Curves[0].Given[0][0], .Y = twin.Curves[0].Given[1][0]};
      const Vector center = ball.first - origin;
      const double radius = ((ball.second + (4.0 * twin.Shift)) * Grow) + Tiny;
      const double grain = static_cast<double>(second.Degree) * Unit * twin.Scale * (1.0 + 0x1p-30);
      List<Slice>  stack;
      int          budget = Spending / 4;
      stack.push_back(
          Slice{.At = second.Moved, .Lo = 0.0, .Hi = 1.0, .Off = Offset(twin, 1), .Depth = 0});
      while (!stack.empty())
      {
        const Slice slice = stack.back();
        stack.pop_back();
        const auto   box   = SliceBox(slice, second.Degree);
        const double aside = std::max({box[0] - center.X, center.X - box[1], 0.0});
        const double above = std::max({box[2] - center.Y, center.Y - box[3], 0.0});
        const double away =
            std::max(std::min(slice.Lo - hull.Start, hull.End - slice.Hi),
                     (Trig::Hypot(aside, above) * Shrink) - ((radius + (2.0 * slice.Off)) * Grow));
        if (away >= 0.0)
        {
          continue;
        }
        if (--budget < 0)
        {
          return false;
        }
        const auto halves = Halved(slice, second.Degree, grain);
        stack.push_back(halves[1]);
        stack.push_back(halves[0]);
      }
      return true;
    }

    /**
     *  Whether no stretch of the first curve lies along the second as the thorough way finds one:
     *  nine places of it evenly spaced all within the tolerance of the second, more than four
     *  tolerances long, and not parting from the piece of the second it is weighed against.
     *  Unmerged keeps such a stretch within one band, which runs no farther than its run: where
     *  that is no more than four tolerances, the stretch is not long enough. Else, as weighing a
     *  crossing shows, the piece lies in the hull of the crossing's boxes, where they part, if no
     *  point of the second outside the hull comes within R of A(m): the isolation shows it where R
     *  is no more than it kept apart, as Unmerged keeps every other crossing's band of R from this
     *  one's; else a ball about the crossing of R and the run holds no point of the second outside
     *  it.
     *
     *  @param  twin     The curves.
     *  @param  settled  The crossings, weighed.
     *  @param  measures The lengths, and how far apart the isolation kept what does not cross.
     *  @return          Whether none does.
     */
    auto Lain(const Twin &twin, std::span<const Settled> settled,
              const std::pair<Levels, double> &measures) -> bool
    {
      const auto &[levels, kept] = measures;
      return std::ranges::all_of(
          settled,
          [&twin, &levels, kept](const Settled &one) -> bool
          {
            return one.Runs[0] <= 4.0 * levels.Length * Shrink
                || one.Reach <= kept
                || Alone(twin, {one.At, one.Reach + one.Runs[0] + levels.Grain}, one.Hull);
          });
    }

    /**
     *  Whether two curves meet where the kernel can never be sure of them, as the thorough way
     *  would decide: an end of one within a length of an end of the other, as consecutive pieces of
     *  a path are; or the two within the length all along, their control points, the first's
     *  against the second's one way round or the other, each that near the other's, so that every
     *  point of one is that near the point of the other at the same place, or at the same place
     *  from the other end, and the isolation could never keep them apart.
     *
     *  @param  twin   The curves.
     *  @param  length The length.
     *  @return        Whether they may.
     */
    auto Abutting(const Twin &twin, double length) -> bool
    {
      const Bent       &first  = twin.Curves[0];
      const Bent       &second = twin.Curves[1];
      const std::size_t degree = std::min(first.Degree, second.Degree);
      double            ahead  = 0.0;
      double            back   = 0.0;
      double            ends   = std::numeric_limits<double>::infinity();
      for (std::size_t index = 0; index <= degree; ++index)
      {
        for (std::size_t axis = 0; axis < 2; ++axis)
        {
          const double mine = first.Given.at(axis).at(index);
          ahead             = std::max(ahead, std::abs(mine - second.Given.at(axis).at(index)));
          back = std::max(back, std::abs(mine - second.Given.at(axis).at(second.Degree - index)));
        }
      }
      for (const std::size_t one : {std::size_t{0}, first.Degree})
      {
        for (const std::size_t other : {std::size_t{0}, second.Degree})
        {
          ends = std::min(ends,
                          std::max(std::abs(first.Given[0].at(one) - second.Given[0].at(other)),
                                   std::abs(first.Given[1].at(one) - second.Given[1].at(other))));
        }
      }
      return std::min({ahead, back, ends}) <= length;
    }

    /**
     *  Whether a curve is straight to within an angle as the thorough way finds it, or may be:
     *  every control step within the angle of its chord, and running its way; its chord of some
     *  length.
     *
     *  @param  bent    The curve.
     *  @param  radians The angle.
     *  @return         Whether it may be.
     */
    auto Straight(const Bent &bent, double radians) -> bool
    {
      const auto pole = [&bent](std::size_t index) -> Point
      {
        return {.X = bent.Given[0].at(index), .Y = bent.Given[1].at(index)};
      };
      const Vector way  = pole(bent.Degree) - pole(0);
      const double sine = Trig::Sin(radians);
      double       off  = -Dot(way, way);
      for (std::size_t index = 1; index <= bent.Degree; ++index)
      {
        const Vector step  = pole(index) - pole(index - 1);
        const double sizes = Length(step) * Length(way);
        off = std::max({off, std::abs(Cross(step, way)) - ((sine * sizes * Grow) + Tiny),
                        -(Dot(step, way) + (0x1p-50 * sizes) + Tiny)});
      }
      return !(off > 0.0);
    }
  } // namespace

  auto Twined(const std::array<const Section *, 2> &sections, const Tolerance &tolerance)
      -> std::optional<Twains>
  {
    if (!std::ranges::all_of(sections,
                             [](const Section *section) -> bool
                             {
                               return std::ranges::all_of(
                                   std::span{section->Poles}.first(section->Degree + 1),
                                   [](Point pole) -> bool { return Held({pole.X, pole.Y}); });
                             }))
    {
      return std::nullopt;
    }
    const Point origin = sections[0]->Poles.front();
    Twin        twin   = {.Curves = {Made(*sections[0], origin), Made(*sections[1], origin)},
                          .Scale  = 0.0,
                          .Shift  = 0.0,
                          .Miss   = 0.0,
                          .Grain  = 0.0};
    double      size   = 0.0;
    for (const Bent &bent : twin.Curves)
    {
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        for (std::size_t index = 0; index <= bent.Degree; ++index)
        {
          twin.Scale = std::max(twin.Scale, std::abs(bent.Moved.at(axis).at(index)));
          size       = std::max(size, std::abs(bent.Given.at(axis).at(index)));
        }
      }
    }
    /**
     *  A moved coordinate is rounded once, so off by u Scale at most; a gap worked out by de
     *  Casteljau's algorithm on both curves, three steps of a double at the scale a round, and
     *  their difference rounded once; and each curve's control points as given off the true ones by
     *  their slack, which the gap and the thorough way's own working carry as well.
     */
    const auto   degrees = static_cast<double>(twin.Curves[0].Degree + twin.Curves[1].Degree);
    const double slacks  = sections[0]->Slack + sections[1]->Slack;
    twin.Shift           = (Unit * twin.Scale * Grow) + Tiny;
    twin.Miss =
        (((((3.0 * degrees) + 2.0) * Unit * twin.Scale) + (2.0 * twin.Shift) + slacks) * Grow)
        + Tiny;
    twin.Grain          = (64.0 * Unit * size) + slacks;
    const Levels levels = {.Length = tolerance.Length,
                           .Within = (tolerance.Length + twin.Grain) * Grow,
                           .Twice  = ((2.0 * tolerance.Length) + twin.Grain) * Grow,
                           .Grain  = twin.Grain,
                           .Sine   = Trig::Sin(tolerance.Radians)};
    const double kept   = Wide * Reach(tolerance);
    if (Abutting(twin, kept))
    {
      return std::nullopt;
    }
    const std::optional<Cages> founds = BoxedRoots(twin, kept);
    if (!founds
        || (Straight(twin.Curves[0], tolerance.Radians)
            && Straight(twin.Curves[1], tolerance.Radians)))
    {
      return std::nullopt;
    }
    List<Settled> settled;
    for (const Rootbox &root : Coalesced(*founds))
    {
      const std::optional<Settled> one = Certified(twin, root, *founds);
      if (!one)
      {
        return std::nullopt;
      }
      settled.push_back(Weigh(twin, *one, {levels, kept}));
    }
    if (!std::ranges::all_of(settled, [&twin, &levels](const Settled &one) -> bool
                             { return Clear(twin, one, levels); })
        || !Unmerged(settled, {levels, kept})
        || !Lain(twin, settled, {levels, kept}))
    {
      return std::nullopt;
    }
    Twains twains;
    for (const Settled &one : settled)
    {
      twains.push_back({.At     = one.At,
                        .Along  = one.Along,
                        .Near   = one.Bands,
                        .Twice  = {Band(one, levels.Twice, 0), Band(one, levels.Twice, 1)},
                        .Runs   = std::max(one.Runs[0], one.Runs[1]),
                        .Within = (levels.Within + levels.Grain) * Grow});
    }
    return twains;
  }
} // namespace Geoxy::Swift
