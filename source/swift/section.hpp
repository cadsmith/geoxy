/***************************************************************************************************
 *  @file      source/swift/section.hpp
 *  @brief     A quadratic or cubic piece of a curve as the swift kernels take it: its control
 *             points rounded, how far they may be from the true ones, the stretch of the curve's
 *             parameter it runs over, and its true points in double-double with a bound.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_SECTION_HPP
#define GEOXY_SWIFT_SECTION_HPP

#include "scratch.hpp"
#include "sure.hpp"
#include "swift/casteljau.hpp"

#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"

#include <array>
#include <cstddef>
#include <optional>

namespace Geoxy::Swift
{
  /**
   *  As many weights of 1 as a piece has room for.
   *
   *  @tparam CountT How many.
   */
  template <std::size_t CountT>
  inline constexpr std::array<double, CountT> Evened = [] -> std::array<double, CountT>
  {
    std::array<double, CountT> ones = {};
    ones.fill(1.0);
    return ones;
  }();

  /**
   *  A polynomial or rational piece of a curve: a Bézier curve whole, or a span of a spline, its
   *  Bézier control points the span's blossom at its ends.
   *
   *  A span's control points are worked out in doubles, each within the slack of the true one a
   *  coordinate, so that the piece they make is within the slack of the true piece everywhere: the
   *  kernels isolate crossings on them, charging the slack. A crossing is proved on the true piece,
   *  its points worked out from the spline's own poles and knots in double-double.
   *
   *  @tparam CountT Room for its control points: its degree and one at the most.
   */
  template <std::size_t CountT>
  struct Sectional
  {
    std::array<Point, CountT> Poles  = {};               /**< The control points, rounded. */
    std::size_t               Degree = 0;                /**< The degree, from two.        */
    double                    Slack  = 0.0;              /**< The most a control point's
                                                              coordinate is off the true one's.   */
    Range Span = {.Start = 0.0, .End = 1.0};             /**< The curve's
                                           parameter over the piece.           */
    Box           Hull   = {};                           /**< A box the true piece cannot leave. */
    const Spline *Source = nullptr;                      /**< The spline it is a span of; none
                                                              where its control points are exact. */
    std::size_t Knot = 0;                                /**< The index of the knot that starts the
                                                              span.                               */
    std::array<double, CountT> Weights = Evened<CountT>; /**< A rational span's
                                                      weights, rounded, within the slack.  */
    bool Rated = false;                                  /**< Whether it is a rational span's: its
                                                              control points are then w_i P_i, its
                                                              weights the w_i, each within the
                                                              slack.                              */
    Point Base = {}; /**< The point its control points are taken from: they are the true ones
                          less it, w_i (P_i - B) where rational, so that a span far from nothing
                          keeps its digits; nothing for a narrow section's, and a Bézier
                          curve's.                                                          */
  };

  /**
   *  A quadratic or cubic polynomial piece, as the kernels of Bézier curves take it.
   */
  using Section = Sectional<Highest + 1>;

  /**
   *  A piece of degree two to Broadest, polynomial or rational, as the kernel of a conic takes it.
   */
  using Broad = Sectional<Broadest + 1>;

  /**
   *  A quadratic Bézier curve as a section: its control points exact, over its whole parameter.
   *
   *  @param  quadratic The curve.
   *  @return           The section.
   */
  auto Sectioned(const Bezier::Quadratic &quadratic) -> Section;

  /**
   *  A cubic Bézier curve as a section.
   *
   *  @param  cubic The curve.
   *  @return       The section.
   */
  auto Sectioned(const Bezier::Cubic &cubic) -> Section;

  /**
   *  A polynomial spline cut into sections, one per span that is not empty, in order.
   *
   *  Each control point is the span's blossom at its ends, de Boor's rounds each a mix of two cells
   *  by a fraction of knot differences, all in doubles: every exact fraction lies from 0 to 1,
   *  every cell is a convex combination of the span's poles, and a round adds at most ten steps u
   *  of a double at the largest size M of their coordinates to what the cells were off by, which it
   *  carries as it was (tests/proof/section.py): within 11 n u M of the true one after n rounds.
   *
   *  @param  spline The spline.
   *  @return        Its sections; nothing where it is not valid, its weights differ, its degree is
   *                 not two or three, a number lies outside the range the kernels take, or a
   *                 fraction of its knots falls among the doubles below the normal ones.
   */
  auto Sectioned(const Spline &spline) -> std::optional<List<Section>>;

  /**
   *  A spline of degree two to Broadest cut into broad sections, polynomial or rational: a
   *  polynomial spline's spans as Sectioned cuts them, of any such degree; a rational one's spans
   *  each the blossom of its homogeneous poles (w P, w), the control points w_i P_i and the weights
   *  w_i, all within the slack 11 n u M of the true ones, M their largest coordinate, as de Boor's
   *  rounds carry the products' first rounding as they are (tests/proof/section.py, for any n). The
   *  hull holds the true control points, each w_i P_i / w_i, the true weights above zero.
   *
   *  @param  spline The spline.
   *  @param  near   The box outside which no span is cut: a span whose poles' box misses it, which
   *                 the span lies in, is left out.
   *  @return        Its sections; nothing where it is not valid, its degree is below two or past
   *                 Broadest, a number lies outside the range the kernels take, a rational one's
   *                 weight is not above 2^-200, a fraction of its knots falls among the doubles
   *                 below the normal ones, or a span's weight may be nothing.
   */
  auto Weighed(const Spline &spline, const Box &near) -> std::optional<List<Broad>>;

  /**
   *  A rational section's true point in homogeneous coordinates at a place, a number with its
   *  bound: de Boor's algorithm on the spline's own poles, each w P and w, at the parameter the
   *  place stands for, every rounding and the place's bound carried.
   *
   *  @pre            The section is a rational spline's span.
   *  @param  section The section.
   *  @param  along   The place, from 0 to 1, with its bound.
   *  @return         The point's numerator's x and y and its weight.
   */
  auto Homogeneously(const Broad &section, const Sure::Number &along)
      -> std::array<Sure::Number, 3>;

  /**
   *  The most a rational section's derivative by its own place may be, from 0 to 1: (N' - p W') /
   *  W, N' - p W' = n sum (w_i+1 (P_i+1 - p) - w_i (P_i - p)) B_i, no more than 2 n w_max times the
   *  hull's diagonal, over w_min, the weights within the slack.
   *
   *  @pre            The section is a rational spline's span.
   *  @param  section The section.
   *  @return         The most |p'| may be.
   */
  auto Fastest(const Broad &section) -> double;

  /**
   *  A section's true point at a place on the grid of 2^-53, so that 1 - t is a double: a Bézier
   *  curve's by de Casteljau's algorithm on its control points, a spline's span's by de Boor's on
   *  the spline's own poles at the parameter the place stands for, both in double-double, every
   *  rounding bounded.
   *
   *  @param  section The section.
   *  @param  along   The place, from 0 to 1, on the grid.
   *  @return         The point's x and y, each with its bound.
   */
  auto Exactly(const Section &section, double along) -> std::array<Sure::Number, 2>;

  /**
   *  A broad section's true point at a place on the grid of 2^-53, as Exactly gives a section's.
   *
   *  @param  section The broad section.
   *  @param  along   The place, from 0 to 1, on the grid.
   *  @return         The point's x and y, each with its bound.
   */
  auto Exactly(const Broad &section, double along) -> std::array<Sure::Number, 2>;

  /**
   *  The curve's parameter at a place of a section: the place itself for a Bézier curve, the span's
   *  start and its width times the place for a spline's.
   *
   *  @param  section The section.
   *  @param  along   The place, with its bound.
   *  @return         The curve's parameter, with its bound.
   */
  auto Spread(const Section &section, const Sure::Number &along) -> Sure::Number;

  /**
   *  The curve's parameter at a place of a broad section, as Spread gives a section's.
   *
   *  @param  section The broad section.
   *  @param  along   The place, with its bound.
   *  @return         The curve's parameter, with its bound.
   */
  auto Spread(const Broad &section, const Sure::Number &along) -> Sure::Number;
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_SECTION_HPP */
