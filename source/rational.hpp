/***************************************************************************************************
 *  @file      source/rational.hpp
 *  @brief     A spline's length, area and moments strip by strip: its spans cut where their weights
 *             are uneven, each strip's integrals by Gauss-Legendre rules whose error is bounded
 *             from the roots of its speed and of its weight.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_RATIONAL_HPP
#define GEOXY_RATIONAL_HPP

#include "aberth.hpp"
#include "eased.hpp"
#include "hodograph.hpp"
#include "kernel/piece.hpp"
#include "quadrature.hpp"

#include "geoxy/domain.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace Geoxy
{
  /**
   *  A polynomial's roots in a strip's parameter, and its leading coefficient's size: what bounds
   *  it off and on the real line, as |a| times the product of the distances to its roots.
   */
  struct Singular
  {
    Aberth::Roots Roots = {};  /**< The roots, in the strip's parameter.          */
    double        Lead  = 0.0; /**< The log2 of the leading coefficient's length. */
  };

  /**
   *  A piece of a spline's span, ready to be measured in its own parameter from 0 to 1: a span
   *  whole, or a part of one whose weight has a root crowding an end of it.
   *
   *  Its control points are homogeneous, taken less the point of the curve at its span's middle
   *  times their weights: a span's points are its curve's to every digit, not lost in its poles'
   *  distance from the origin, and where a heavy weight holds the curve at its pole for most of the
   *  span, as near that pole as the curve is. Its weights are times a power of two that brings the
   *  heaviest near 1, unless that would leave the lightest less than a double of all its digits. A
   *  polynomial span's hodograph is worked out apart, from the differences of its poles, as the
   *  derivative's blossoms, so that a span far smaller than its poles' spread keeps every digit of
   *  how fast it runs.
   */
  struct Strip
  {
    Kernel::Piece                 Shape  = {};    /**< Its control points.                   */
    Point                         Offset = {};    /**< What they are taken less.             */
    Slice                         Where  = {};    /**< Where it lies on the spline.          */
    std::array<Vector, MaxDegree> Steps  = {};    /**< A polynomial span's hodograph.        */
    Singular                      Poles  = {};    /**< Its weight's roots; none where plain. */
    bool                          Plain  = false; /**< Whether its weights are all one.      */
  };

  /**
   *  A span of a valid spline as one strip, uncut: what a stretch run on past the spline's ends is
   *  measured on, the span's own parameter telling a place far past its end where a strip cut
   *  beside it would not.
   *
   *  @param  spline The spline, valid.
   *  @param  span   The index of the knot that starts the span, the span not empty.
   *  @return        The strip.
   */
  auto Whole(const Spline &spline, std::size_t span) -> Strip;

  /**
   *  A valid spline's strips, in order: each span that is not empty, cut 2^-24 of it from an end
   *  its weight has a root nearer than that, and its parts again, until none is, or until cut 42
   *  times. A root that near an end is told apart from it to no more than its own distance's digits
   *  less 24; cut so, every strip's weight's roots lie a fair part of it away from it, so that its
   *  integrals are smooth in its own parameter, and a node near either end of it is placed to the
   *  digits its distance from that end has. A rush past the curve within 2^-1008 of a span's end,
   *  beside weights 2^1008 apart in it, is beyond that: its strip is taken as it is, and measured
   *  short of where a double tells it.
   *
   *  @param  spline The spline, valid.
   *  @return        Its strips.
   */
  auto Stripped(const Spline &spline) -> std::vector<Strip>;

  /**
   *  A strip's parameter at a spline's parameter, told from the nearer end of its span, so that a
   *  strip a sliver wide beside the span's end is told to every digit the spline's parameter has.
   *
   *  @param  where Where the strip lies on the spline.
   *  @param  along The spline's parameter.
   *  @return       The strip's parameter, before 0 and past 1 too.
   */
  auto Local(const Slice &where, double along) -> double;

  /**
   *  Whether a spline's parameter lies at or before where a strip starts, told from the nearer end
   *  of its span in double-double, so that a strip whose start rounds to its span's end in the
   *  spline's parameter is still told after the span's start and before its end.
   *
   *  @param  where Where the strip lies on the spline.
   *  @param  along The spline's parameter.
   *  @return       Whether it does.
   */
  auto Ahead(const Slice &where, double along) -> bool;

  /**
   *  Whether a spline's parameter lies at or past where a strip ends, told as Ahead tells it.
   *
   *  @param  where Where the strip lies on the spline.
   *  @param  along The spline's parameter.
   *  @return       Whether it does.
   */
  auto Past(const Slice &where, double along) -> bool;

  /**
   *  A strip's speed, ready to be integrated: |N| / W^2 times a factor, N the numerator of its
   *  derivative, W X' - W' X for X and W its homogeneous coordinates and weight, of degree 2 n - 2,
   *  or for a polynomial strip its hodograph, of degree n - 1, and W 1; both in Bernstein form, N
   *  scaled near 1 and its coordinates' scale in the factor, and the roots of each; and the strip's
   *  control points for its points, their coordinates scaled likewise.
   */
  struct Pace
  {
    std::size_t                       Degree   = 0;   /**< N's degree.                         */
    std::array<Vector, Kernel::Terms> Top      = {};  /**< N's coefficients.                   */
    std::size_t                       Weighing = 0;   /**< W's degree: 0 where it is 1.        */
    std::array<double, MaxDegree + 1> Bottom   = {};  /**< W's coefficients.                   */
    double                            Factor   = 1.0; /**< What |N| / W^2 is times.            */
    std::size_t                       Order    = 0;   /**< The strip's degree.                 */
    std::array<Vector, MaxDegree + 1> Places   = {};  /**< Its coordinates, weighted, scaled.  */
    std::array<double, MaxDegree + 1> Masses   = {};  /**< Its weights.                        */
    int                               Size     = 0;   /**< The log2 its coordinates are under. */
    Singular                          Turns    = {};  /**< N's roots.                          */
    Singular                          Poles    = {};  /**< W's roots.                          */
  };

  /**
   *  A strip's speed made ready to be integrated.
   *
   *  @param  strip The strip.
   *  @return       Its speed.
   */
  auto Paced(const Strip &strip) -> Pace;

  /**
   *  The length a strip runs over a stretch of its parameter, past 0 and 1 too.
   *
   *  The speed |N| / W^2 is analytic but at N's roots and their conjugates, where it has branch
   *  points, and W's, where it has poles. The stretch is cut where the real part falls of a root
   *  whose Bernstein ellipse about the stretch is thinner than rho = 3, there the speed may turn
   *  back or have a cusp or peak, and each piece is halved while the ellipse about it that reaches
   *  the nearest root counted is thinner than 3. On a piece of half width h, with |f| at most M on
   *  the ellipse of rho, Gauss-Legendre with n points is off by at most (64 / 15) M h rho^(2 - 2 n)
   *  / (rho^2 - 1) (Trefethen, Approximation Theory and Approximation Practice, theorem 19.3). On
   *  that ellipse, of half axis R = h (rho + 1 / rho) / 2, |N| is at most |a| times the product of
   *  each root's distance from the piece's middle plus R, and |W| at least |b| times the product of
   *  each root's half axis A less R: a root lies on the ellipse about the piece of half axis A, the
   *  half sum of its distances from the piece's ends, and every point of a smaller one is at least
   *  A - R from it, by the triangle inequality on the two sums of distances to the foci. M is
   *  finite on the ellipse through N's nearest root, and grows without end toward W's: rho is the
   *  nearest root's where it is N's, and the powers 9 / 10, 3 / 4 and 1 / 2 of W's nearest's, none
   *  past N's, whichever takes the fewest points. Each piece takes the fewest points, in pairs,
   *  that keep the bound within three quarters of its share of the error: its width's share of the
   *  whole, or 2^-55 of its chord where that is more, since beside a pole a narrow piece carries
   *  much of the length.
   *
   *  A root of N whose imaginary part q is so small that sqrt((t - p)^2 + q^2) differs from |t - p|
   *  by little over a piece is taken as real there, as Measured takes a cubic's: the speed is then
   *  |t - p| times what is analytic, the cut keeping p off the piece's inside. It differs from the
   *  speed by q^2 / (2 |t - p|) times the rest of the speed at most, the rest bounded on the piece
   *  as M is on its ellipse with R = h, so by at most that times q^2 (log(1 + 2 h / q) + 1) in the
   *  integral, and 568 times that over the rule: a quarter of the share bounds both.
   *
   *  The error asked is at least 2^-55 of a lower bound on the length: the polygon through the
   *  points at the stretch's quarters, less 2^-40 of it for their rounding. Each speed is summed in
   *  Bernstein form in the ratio of the nearer of t and 1 - t to the other, and divided by one plus
   *  that ratio to its power, which so belongs to the very parameter the sum was taken at: no step
   *  of 1 - t is grown by the degree. The pieces are summed with their roundoffs; what is left is
   *  rounding in the speeds, a few steps of a double, and the roots, found by Aberth's iteration
   *  and the leading coefficients from them and the polynomial's value at the largest of three
   *  points, are taken as they are found: their own rounding is not in the bound.
   *
   *  @param  pace    The strip's speed.
   *  @param  stretch The stretch, its start at most its end.
   *  @param  effort  How far off the length may be, and how many halvings, at most Halved, may be
   *                  spent; once spent, each piece left takes the most points.
   *  @return         The length: 0 for a strip that does not move or a stretch of no width, and
   *                  without end for any other whose width is.
   */
  auto Measured(const Pace &pace, Range stretch, Effort effort) -> double;

  /**
   *  Where on a strip the length from its start is a distance: Newton's method on Measured, its
   *  derivative the speed, kept within a bracket the length is known to cross, and halving the
   *  bracket where a step would leave it.
   *
   *  @param  pace      The strip's speed.
   *  @param  distance  The distance, from 0 to the strip's length.
   *  @param  tolerance How far off the length at the answer may be.
   *  @return           The strip's parameter there, from 0 to 1.
   */
  auto Walked(const Pace &pace, double distance, double tolerance) -> double;

  /**
   *  What a strip's enclosure is asked about: the point, how far off it may be, and whether its
   *  moment is asked too.
   */
  struct Inquiry
  {
    Point       Origin   = {};     /**< The point the area and moment are taken about.  */
    double      Error    = 0.0;    /**< How far off the area and each moment may be.    */
    bool        Moments  = false;  /**< Whether the moment is asked too.                */
    std::size_t Halvings = Halved; /**< How many halvings may be spent, Halved at most. */
  };

  /**
   *  What a strip encloses with a point: the area of the fan from the point and its first moment
   *  about the point, as Green's theorem gives them.
   *
   *  With X and W the strip's homogeneous coordinates less the point's times the weight, and its
   *  weight, the fan's area grows at X x X' / W^2 / 2, X x X' a polynomial of degree 2 n - 2,
   *  worked out as Y x Y' + d x (W Y' - W' Y) for Y the coordinates less the strip's offset and d
   *  the offset less the point, the two scaled apart so that neither overflows; its moment grows at
   *  Y / W times that, a third as fast, and d times two thirds of the area. Both are analytic but
   *  at W's roots, poles, and integrated as Measured integrates the speed, but cut at the poles'
   *  real parts alone, and each numerator bounded on the ellipse by its largest Bernstein
   *  coefficient times (|z| + |1 - z|)^d, at most max(1, |2 c - 1|) + 2 R there, c the piece's
   *  middle: B_i(z) taken at its length, the binomial theorem sums them to that power. The error
   *  asked is at least 2^-55 of what each would be with its rate taken at its length, as the rule
   *  of 2 points gives it, on the whole strip and on each piece.
   *
   *  @param  strip The strip.
   *  @param  asked The point, how far off the area and moment may be, whether the moment is asked,
   *                and how many halvings may be spent.
   *  @return       The area and the moment, taken about the point; the moment as far as the area's
   *                points tell it, where it is not asked.
   */
  auto Enclosed(const Strip &strip, const Inquiry &asked) -> Kernel::Enclosure;
} // namespace Geoxy

#endif /* GEOXY_RATIONAL_HPP */
