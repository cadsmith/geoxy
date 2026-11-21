/***************************************************************************************************
 *  @file      source/lengths/spans.hpp
 *  @brief     The length along a spline or a Bézier curve, span by span, and along a polyline, leg
 *             by leg.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_LENGTHS_SPANS_HPP
#define GEOXY_LENGTHS_SPANS_HPP

#include "hodograph.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <cstddef>
#include <initializer_list>
#include <vector>

namespace Geoxy::Lengths
{
  /**
   *  How long a straight run that never ends is.
   *
   *  @param  way Its way.
   *  @return     Infinity, or 0 for no way.
   */
  auto Forever(Vector way) -> double;

  /**
   *  The length of a stretch of a curve that runs at one speed.
   *
   *  @param  stretch The stretch of its parameter.
   *  @param  speed   Its speed, the length of each 1 of its parameter.
   *  @return         The length, negative for a stretch that runs backward.
   */
  auto Evenly(Range stretch, double speed) -> double;

  /**
   *  The parameter at a distance along a curve that runs at one speed.
   *
   *  @param  start    Where the distance is measured from.
   *  @param  speed    Its speed.
   *  @param  distance The distance.
   *  @return          The parameter; the start for a curve of no speed.
   */
  auto Evenly(double start, double speed, double distance) -> double;

  /**
   *  A stretch run forward: its start at most its end.
   *
   *  @param  stretch The stretch, either way round.
   *  @return         The stretch from its lower end to its higher.
   */
  auto Forward(Range stretch) -> Range;

  /**
   *  Whether every coordinate of some points is finite, told without a branch per coordinate: a sum
   *  of each times 0 is not a number exactly where one is not finite.
   *
   *  @param  points The points.
   *  @return        Whether they are all finite.
   */
  auto AllFinite(std::initializer_list<Point> points) -> bool;

  /**
   *  A span of a polynomial spline of degree 3 at most: its hodograph in the span's own parameter
   *  from 0 to 1, and the span.
   */
  struct Haul
  {
    Hodograph Shape = {}; /**< The hodograph.                      */
    Range     Span  = {}; /**< The span of the spline's parameter. */
  };

  /**
   *  Whether a spline is a polynomial one of degree 3 at most, whose length its legs give: one with
   *  no weights, or with every weight the same, is polynomial.
   *
   *  @param  spline The spline.
   *  @return        Whether it is.
   */
  auto Unweighted(const Spline &spline) -> bool;

  /**
   *  A valid polynomial spline's legs, one per span that is not empty.
   *
   *  Its derivative is the spline of a degree less whose control points are d (P_(i+1) - P_i) /
   *  (u_(i+d+1) - u_(i+1)), on its knots but the first and the last; over a span u_k to u_(k+1), of
   *  width w, its hodograph in the span's own parameter is w times that, whose Bézier points are
   *  its blossoms at the span's ends, u_k repeated d - 1 - j times and u_(k+1) j times. The span
   *  lies within every u_(i+1) to u_(i+d+1) it reads, so each control point is the difference of
   *  two poles times w over that width, at most 1: none is divided by nothing, and none of the
   *  differences is grown. The pole differences are taken first, so that a curve far from the
   *  origin keeps its precision.
   *
   *  @param  spline The spline, valid, of degree 3 at most and of one weight or none.
   *  @return        Its legs, in order.
   */
  auto Hauls(const Spline &spline) -> std::vector<Haul>;

  /**
   *  The length of a stretch of a leg.
   *
   *  @param  leg     The leg.
   *  @param  stretch The stretch of the spline's parameter, run forward, past the span's ends too.
   *  @param  error   How far off the length may be.
   *  @return         The length.
   */
  auto Traced(const Haul &leg, Range stretch, double error) -> double;

  /**
   *  The length of a stretch of a polynomial spline, leg by leg, the first run on before the domain
   *  and the last past it, each leg's share of the accuracy its share of the stretch, the lengths
   *  summed with their roundoffs.
   *
   *  @param  legs     The spline's legs.
   *  @param  stretch  The stretch, run forward.
   *  @param  accuracy How far off the length may be.
   *  @return          The length.
   */
  auto Legwise(const std::vector<Haul> &legs, Range stretch, Accuracy accuracy) -> double;

  /**
   *  The length of a stretch of a spline strip by strip, the lengths summed with their roundoffs. A
   *  strip the stretch covers whole is measured from 0 to 1, told so by its ends in double-double,
   *  so that no strip's end is moved by a rounding, nor a strip lost whose ends round to one
   *  parameter. What runs on before the domain or past it is measured on its end span whole, whose
   *  own parameter tells a place far off where a strip cut beside the end would not.
   *
   *  @param  spline  The spline, valid.
   *  @param  stretch The stretch, run forward.
   *  @param  error   How far off the length may be.
   *  @return         The length.
   */
  auto Spanwise(const Spline &spline, Range stretch, double error) -> double;

  /**
   *  A length given the sign of the way a stretch runs.
   *
   *  @param  stretch The stretch.
   *  @param  length  The length, run forward.
   *  @return         The length, negative for a stretch that runs backward.
   */
  auto Signed(Range stretch, double length) -> double;
} // namespace Geoxy::Lengths

#endif /* GEOXY_LENGTHS_SPANS_HPP */
