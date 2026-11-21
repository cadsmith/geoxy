/***************************************************************************************************
 *  @file      source/near/feet.hpp
 *  @brief     The feet of the perpendiculars from a point to a curve: straight runs and arcs in
 *             closed form, every other curve by the kernel's feet on its pieces, listed in order.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_NEAR_FEET_HPP
#define GEOXY_NEAR_FEET_HPP

#include "scratch.hpp"
#include "walk.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <cstddef>
#include <span>

namespace Geoxy::Near
{
  /**
   *  A few places put in order along the curve, those at one parameter kept in the order given: an
   *  insertion sort, which takes no memory.
   *
   *  @param places The places.
   */
  void Lineup(std::span<Place> places);

  /**
   *  The point of a circle on the way from its center to a point, or on the way back: c + R v / |v|
   *  for v the point less the center and R the radius, signed for the way back, rounded once at
   *  each coordinate.
   *
   *  v is exact (Apart), |v|^2 within 2^-100 of itself (Dot), |v| within 2^-102 (Root), R / |v|
   *  within 2^-100 (Over), and each coordinate before its rounding within 2^-99 (|c| + |R|)
   *  (Moved): the nearest double to the exact point but a hair from halfway between two. Where |v|
   *  is past 2^450 or under 2^-450, whose square would leave the doubles, it is worked in plain
   *  doubles, within a few steps of a double of |R|.
   *
   *  @param  center The circle's center.
   *  @param  radius Its radius, signed: the way back where negative.
   *  @param  point  The point, not the center.
   *  @return        The circle's point.
   */
  auto Dropped(Point center, double radius, Point point) -> Point;

  /**
   *  A polyline segment's foot of a point, worked in lean double-double from the segment's start
   *  and rounded once at each coordinate: the double nearest the exact foot but a hair from
   *  halfway, every step within 2^-98 of its terms.
   *
   *  A straight segment's is the exact projection a + d (w . d) / (d . d), d = b - a and w = p - a
   *  exactly. A bulge's is its circle's point toward the point, or away from it, worked along its
   *  chord as Bowed works its points, with no center rounded: the center m = d / 2 + k d turned a
   *  quarter, k = (1 - b^2) / (4 b), and the radius |m|, the start being on the circle: the foot m
   *  + v sqrt(|m|^2 / |v|^2) for v = w - m. Where a way's coordinate is past 2^450, or one not
   *  nothing is under 2^-450, its square would leave the doubles: there it is At's point.
   *
   *  @param  step  The segment, and At's parameter of the foot along it.
   *  @param  point The point.
   *  @param  side  For a bulge, 1 for the point toward the point, -1 away from it.
   *  @return       The foot.
   */
  auto Trod(const Step &step, Point point, double side) -> Point;

  /**
   *  A line's feet, as a list drawing on the thread's pool.
   *
   *  @param  line  The curve.
   *  @param  point The point.
   *  @return       The feet, in order along the curve.
   */
  auto Listed(const Line &line, Point point) -> List<Place>;

  /**
   *  A ray's feet, as a list drawing on the thread's pool.
   *
   *  @param  ray   The curve.
   *  @param  point The point.
   *  @return       The feet, in order along the curve.
   */
  auto Listed(const Ray &ray, Point point) -> List<Place>;

  /**
   *  An endless line's feet, as a list drawing on the thread's pool.
   *
   *  @param  xline The curve.
   *  @param  point The point.
   *  @return       The feet, in order along the curve.
   */
  auto Listed(const Xline &xline, Point point) -> List<Place>;

  /**
   *  An arc's feet, as a list drawing on the thread's pool.
   *
   *  @param  arc   The curve.
   *  @param  point The point.
   *  @return       The feet, in order along the curve.
   */
  auto Listed(const Arc &arc, Point point) -> List<Place>;

  /**
   *  An ellipse's feet, as a list drawing on the thread's pool.
   *
   *  @param  ellipse The curve.
   *  @param  point   The point.
   *  @return         The feet, in order along the curve.
   */
  auto Listed(const Ellipse &ellipse, Point point) -> List<Place>;

  /**
   *  A polyline's feet, as a list drawing on the thread's pool.
   *
   *  @param  polyline The curve.
   *  @param  point    The point.
   *  @return          The feet, in order along the curve.
   */
  auto Listed(const Polyline &polyline, Point point) -> List<Place>;

  /**
   *  A spline's feet, as a list drawing on the thread's pool.
   *
   *  @param  spline The curve.
   *  @param  point  The point.
   *  @return        The feet, in order along the curve.
   */
  auto Listed(const Spline &spline, Point point) -> List<Place>;

  /**
   *  A quadratic Bézier curve's feet, as a list drawing on the thread's pool.
   *
   *  @param  quadratic The curve.
   *  @param  point     The point.
   *  @return           The feet, in order along the curve.
   */
  auto Listed(const Bezier::Quadratic &quadratic, Point point) -> List<Place>;

  /**
   *  A cubic Bézier curve's feet, as a list drawing on the thread's pool.
   *
   *  @param  cubic The curve.
   *  @param  point The point.
   *  @return       The feet, in order along the curve.
   */
  auto Listed(const Bezier::Cubic &cubic, Point point) -> List<Place>;

  /**
   *  The place nearest a point among some: of those as near, the first along the curve, then the
   *  first given, as a stable sort by parameter and the first strict minimum after it would choose,
   *  with nothing sorted. Where the first along the curve is at no number of distance, none is
   *  nearer, and it is the one.
   *
   *  @pre           There is a place at least.
   *  @param  point  The point.
   *  @param  places The places.
   *  @return        The nearest.
   */
  auto Least(Point point, std::span<const Place> places) -> Place;

  /**
   *  The place of a polyline nearest a point, every segment solved: its feet chained and kept once,
   *  and its vertices, weighed by Least.
   *
   *  @param  polyline The polyline.
   *  @param  point    The point.
   *  @return          The place.
   */
  auto Everywhere(const Polyline &polyline, Point point) -> Place;
} // namespace Geoxy::Near

#endif /* GEOXY_NEAR_FEET_HPP */
