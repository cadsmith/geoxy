/***************************************************************************************************
 *  @file      source/lengths/reach.hpp
 *  @brief     The place a distance along a curve: the stretch it falls in bracketed, and the place
 *             in it found.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_LENGTHS_REACH_HPP
#define GEOXY_LENGTHS_REACH_HPP

#include "lengths/pairs.hpp"
#include "rational.hpp"
#include "walk.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/length.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace Geoxy::Lengths
{
  /**
   *  Whether a stretch of a spline's parameter reaches where the spline, run on past an end, stops
   *  running: where its weight comes to nothing, the curve off to no point there and its length to
   *  that place without end. Told by the first root of the end span's weight, and by the weight at
   *  the stretch's own ends, which is not above nothing past such a place.
   *
   *  @param  spline  The spline, valid.
   *  @param  stretch The stretch.
   *  @return         Whether it does: never for a spline of no weights, nor within its domain.
   */
  auto Stopped(const Spline &spline, Range stretch) -> bool;

  /**
   *  A place of a curve from its point at a parameter: no place, not a number, where the parameter
   *  is none or without end, as for a distance the curve never gets to, or one that is no distance.
   *
   *  @param  point The curve's point at the parameter.
   *  @param  along The parameter.
   *  @return       The place.
   */
  auto Situated(Point point, double along) -> Place;

  /**
   *  What a search for a length goes by: the length from a curve's start to a parameter, its speed
   *  there, and its domain.
   */
  struct Gauge
  {
    std::function<double(double)> Measure = {}; /**< The length from the domain's start. */
    std::function<double(double)> Speed   = {}; /**< The speed at a parameter.           */
    Range                         Domain  = {}; /**< The domain, finite.                 */
  };

  /**
   *  What a search for a length on a curve goes by.
   *
   *  @tparam CurveT  The curve's kind.
   *  @param  curve   The curve.
   *  @param  measure The length from its domain's start to a parameter.
   *  @return         The gauge.
   */
  template <typename CurveT>
  auto Gauged(const CurveT &curve, const std::function<double(double)> &measure) -> Gauge
  {
    return {.Measure = measure,
            .Speed  = [&curve](double along) -> double { return Length(Derivative(curve, along)); },
            .Domain = Domain(curve)};
  }

  /**
   *  What a search for a length on a polynomial spline goes by: its legs and each one's length from
   *  the domain's start to where it starts, found once, so that each measure integrates one leg,
   *  the first run on before the domain and the last past it.
   *
   *  @param  spline   The spline, valid, of degree 3 at most and of one weight or none.
   *  @param  accuracy How far off a measure may be: half for the lengths to the legs, half for the
   *                   leg it falls in.
   *  @return          The gauge.
   */
  auto Coursed(const Spline &spline, Accuracy accuracy) -> Gauge;

  /**
   *  A spline's strips, where each starts on its parameter, and the length from the domain's start
   *  to each, the last the whole length.
   */
  struct Course
  {
    std::vector<Strip>  Strips = {}; /**< The strips.                                  */
    std::vector<double> Starts = {}; /**< Where each starts on the spline's parameter. */
    std::vector<double> Before = {}; /**< The length to each, and the whole length.    */
  };

  /**
   *  A spline's strips and the lengths to each, found once.
   *
   *  @param  spline The spline, valid.
   *  @param  error  How far off the whole length may be, each strip's share its share of the
   *                 domain.
   *  @return        The course.
   */
  auto Charted(const Spline &spline, double error) -> Course;

  /**
   *  What a search for a length on a spline past its domain goes by, strip by strip: each measure
   *  integrates one strip, the first run on before the domain and the last past it, and is without
   *  end across where the spline, run on, stops.
   *
   *  @param  spline The spline, valid.
   *  @param  course Its strips and the lengths to them.
   *  @param  error  How far off a measure of the strip it falls in may be.
   *  @return        The gauge.
   */
  auto Strewn(const Spline &spline, Course course, double error) -> Gauge;

  /**
   *  Where a length measured from a curve's start reaches a distance: Newton's method on the
   *  length, its derivative the speed, kept within a bracket the length is known to cross, and
   *  halving the bracket where a step would leave it.
   *
   *  @param  gauge     What it goes by.
   *  @param  distance  The distance, negative before the start.
   *  @param  tolerance How far off the length at the answer may be.
   *  @return           The parameter; not a number where the curve, run on, never gets as far as
   *                    the distance, so that there is no such place.
   */
  auto Reached(const Gauge &gauge, double distance, double tolerance) -> double;

  /**
   *  How far off a length found by searching may be: the accuracy asked, or a few steps of a double
   *  at the distance.
   *
   *  @param  distance The distance.
   *  @param  accuracy The accuracy asked.
   *  @return          The tolerance.
   */
  auto Allowed(double distance, double accuracy) -> double;

  /**
   *  A polyline segment's length in double-double, within 2^-66 of itself.
   *
   *  A straight one's is |b - a|: the difference exact (Lean::Apart), its square within 2^-100
   *  (Lean::Dot) and the root within 2^-102 (Lean::Root). A bulge's is its radius times its sweep,
   *  c (1 + b^2) / (4 |b|) times 4 |atan b|: the chord's length c so, 1 + b^2 exactly, and atan b /
   *  b within 2^-68 (ArcTangent, Lean::Over), all positive, multiplied in double-double. Where the
   *  chord's square would leave the doubles, past 2^900 or under 2^-900, it is the segment's own
   *  Length, rounded.
   *
   *  @param  step The segment.
   *  @return      Its length.
   */
  auto Lengthed(const Step &step) -> Double;
} // namespace Geoxy::Lengths

#endif /* GEOXY_LENGTHS_REACH_HPP */
