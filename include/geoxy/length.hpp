/***************************************************************************************************
 *  @file      geoxy/length.hpp
 *  @brief     The length of a curve or a stretch of it, and the place at a length: exact where a
 *             closed form has it, to an accuracy where not.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_LENGTH_HPP
#define GEOXY_LENGTH_HPP

#include "geoxy/accuracy.hpp"
#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/place.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

namespace Geoxy
{
  /**
   *  The length of a line.
   *
   *  @param  line The line.
   *  @return      The distance between its ends.
   */
  [[nodiscard]] GEOXY_API auto Length(const Line &line) -> double;

  /**
   *  The length of a ray.
   *
   *  @param  ray The ray.
   *  @return     Infinity, or 0 for a ray of no way.
   */
  [[nodiscard]] GEOXY_API auto Length(const Ray &ray) -> double;

  /**
   *  The length of an xline.
   *
   *  @param  xline The xline.
   *  @return       Infinity, or 0 for an xline of no way.
   */
  [[nodiscard]] GEOXY_API auto Length(const Xline &xline) -> double;

  /**
   *  The length of an arc.
   *
   *  @param  arc The arc.
   *  @return     Its radius times its sweep, both taken positive.
   */
  [[nodiscard]] GEOXY_API auto Length(const Arc &arc) -> double;

  /**
   *  The length of an ellipse, by Carlson's elliptic integrals: exact but for rounding.
   *
   *  Its speed squared is A + R cos(2 (angle - c)), so its length is the larger radius times the
   *  difference of the incomplete elliptic integral of the second kind at its ends' angles, less c,
   *  of parameter 2 R / (A + R); worked out by Carlson's duplication to the last bits, every half
   *  turn of it a complete integral.
   *
   *  @param  ellipse The ellipse.
   *  @return         The length, every turn it sweeps counted.
   */
  [[nodiscard]] GEOXY_API auto Length(const Ellipse &ellipse) -> double;

  /**
   *  The length of a polyline.
   *
   *  @param  polyline The polyline.
   *  @return          The sum of its segments' lengths, a bulged one's as its arc's.
   */
  [[nodiscard]] GEOXY_API auto Length(const Polyline &polyline) -> double;

  /**
   *  The length of a spline.
   *
   *  A polynomial spline of degree 3 at most, one with no weights or with every weight the same, is
   *  measured span by span from its hodograph, as a Bézier curve of its degree is. Any other, one
   *  with weights or of degree 4 or more, is measured strip by strip: each span as its rational
   *  Bézier piece, cut where a root of its weight lies within 2^-24 of an end of it in its own
   *  parameter, so that a heavy weight's rush past the curve is told to every digit; and each
   *  strip's speed integrated by Gauss-Legendre rules whose error is bounded from the roots of the
   *  speed's numerator and of the weight, cut at the turns and cusps, halved toward the near roots
   *  and given the fewest points the bound allows.
   *
   *  Either way the length is within 2^-55 of itself, or the accuracy where that is coarser, and a
   *  few steps of a double for rounding: within 4 steps of 30-digit references on circles and
   *  ellipses written as NURBS, weights 1 to 1e16 and 1e200 apart, degrees 1 to 12, cusps and
   *  curves that run back along themselves.
   *
   *  @param  spline   The spline.
   *  @param  accuracy How far the length may be off.
   *  @return          The length; 0 for a spline that is not valid.
   */
  [[nodiscard]] GEOXY_API auto Length(const Spline &spline, Accuracy accuracy) -> double;

  /**
   *  The length of a quadratic Bézier curve, in closed form.
   *
   *  Its speed is |a| sqrt((t - p)^2 + q^2), whose integral is worked out with both its differences
   *  rationalised, so that nothing cancels where the curve turns back: within a few steps of a
   *  double, whatever the accuracy asked.
   *
   *  @param  quadratic The curve.
   *  @param  accuracy  How far the length may be off; it is never further off than rounding.
   *  @return           The length; not a number for a curve with a coordinate that is none or
   *                    without end.
   */
  [[nodiscard]] GEOXY_API auto Length(const Bezier::Quadratic &quadratic, Accuracy accuracy)
      -> double;

  /**
   *  The length of a cubic Bézier curve, by Gauss-Legendre rules whose error is bounded.
   *
   *  Its hodograph is a complex quadratic, and its speed is analytic but at the hodograph's two
   *  roots and their conjugates: the parameter is cut where a root's real part lies, where the
   *  curve may turn back or have a cusp, each piece is halved while a root lies near it, and each
   *  takes the fewest points whose bound, by the root's Bernstein ellipse, is within its share of
   *  2^-55 of the length, or of the accuracy where that is coarser. The length is then within a few
   *  steps of a double of the true one, cusps and turns back too.
   *
   *  @param  cubic    The curve.
   *  @param  accuracy How far the length may be off.
   *  @return          The length; not a number for a curve with a coordinate that is none or
   *                   without end.
   */
  [[nodiscard]] GEOXY_API auto Length(const Bezier::Cubic &cubic, Accuracy accuracy) -> double;

  /**
   *  The length of any curve.
   *
   *  @param  curve    The curve.
   *  @param  accuracy How far the length may be off, for a kind that is approximated.
   *  @return          The length, as its own kind gives it.
   */
  [[nodiscard]] GEOXY_API auto Length(const Curve &curve, Accuracy accuracy) -> double;

  /**
   *  The length of a stretch of a line.
   *
   *  @param  line    The line.
   *  @param  stretch The stretch of its parameter, past its ends too.
   *  @return         The length, negative for a stretch that runs backward.
   */
  [[nodiscard]] GEOXY_API auto Length(const Line &line, Range stretch) -> double;

  /**
   *  The length of a stretch of a ray.
   *
   *  @param  ray     The ray.
   *  @param  stretch The stretch of its parameter, finite.
   *  @return         The length, negative for a stretch that runs backward.
   */
  [[nodiscard]] GEOXY_API auto Length(const Ray &ray, Range stretch) -> double;

  /**
   *  The length of a stretch of an xline.
   *
   *  @param  xline   The xline.
   *  @param  stretch The stretch of its parameter, finite.
   *  @return         The length, negative for a stretch that runs backward.
   */
  [[nodiscard]] GEOXY_API auto Length(const Xline &xline, Range stretch) -> double;

  /**
   *  The length of a stretch of an arc.
   *
   *  @param  arc     The arc.
   *  @param  stretch The stretch of its parameter, past its ends too.
   *  @return         The length, negative for a stretch that runs backward.
   */
  [[nodiscard]] GEOXY_API auto Length(const Arc &arc, Range stretch) -> double;

  /**
   *  The length of a stretch of an ellipse, by Carlson's integrals.
   *
   *  @param  ellipse The ellipse.
   *  @param  stretch The stretch of its parameter, past its ends too.
   *  @return         The length, negative for a stretch that runs backward.
   */
  [[nodiscard]] GEOXY_API auto Length(const Ellipse &ellipse, Range stretch) -> double;

  /**
   *  The length of a stretch of a polyline.
   *
   *  @param  polyline The polyline.
   *  @param  stretch  The stretch of its parameter, its first and last segments running on past its
   *                   ends.
   *  @return          The length, negative for a stretch that runs backward; 0 for fewer than two
   *                   vertices; not a number for a stretch with an end that is not one.
   */
  [[nodiscard]] GEOXY_API auto Length(const Polyline &polyline, Range stretch) -> double;

  /**
   *  The length of a stretch of a spline, its spans' each as the whole spline's is, cut at its
   *  knots.
   *
   *  @param  spline   The spline.
   *  @param  stretch  The stretch of its parameter, past its ends too.
   *  @param  accuracy How far the length may be off.
   *  @return          The length, negative for a stretch that runs backward; 0 for a spline that is
   *                   not valid; without end for a stretch that reaches where a spline with
   *                   weights, run on, stops running, its weight come to nothing; not a number for
   *                   a stretch with an end that is not one.
   */
  [[nodiscard]] GEOXY_API auto Length(const Spline &spline, Range stretch, Accuracy accuracy)
      -> double;

  /**
   *  The length of a stretch of a quadratic Bézier curve, in closed form.
   *
   *  @param  quadratic The curve.
   *  @param  stretch   The stretch of its parameter, past its ends too.
   *  @param  accuracy  How far the length may be off.
   *  @return           The length, negative for a stretch that runs backward.
   */
  [[nodiscard]] GEOXY_API auto Length(const Bezier::Quadratic &quadratic, Range stretch,
                                      Accuracy accuracy) -> double;

  /**
   *  The length of a stretch of a cubic Bézier curve, by Gauss-Legendre rules whose error is
   *  bounded.
   *
   *  @param  cubic    The curve.
   *  @param  stretch  The stretch of its parameter, past its ends too.
   *  @param  accuracy How far the length may be off.
   *  @return          The length, negative for a stretch that runs backward.
   */
  [[nodiscard]] GEOXY_API auto Length(const Bezier::Cubic &cubic, Range stretch, Accuracy accuracy)
      -> double;

  /**
   *  The length of a stretch of any curve.
   *
   *  @param  curve    The curve.
   *  @param  stretch  The stretch of its parameter.
   *  @param  accuracy How far the length may be off, for a kind that is approximated.
   *  @return          The length, as its own kind gives it.
   */
  [[nodiscard]] GEOXY_API auto Length(const Curve &curve, Range stretch, Accuracy accuracy)
      -> double;

  /**
   *  The place at a distance along a line from its start.
   *
   *  @param  line     The line.
   *  @param  distance The distance, negative before its start, past its length beyond its end.
   *  @return          The place; for a line of no length its start at a distance of nothing, and no
   *                   place, every number of it not a number, at any other, as at a distance that
   *                   is no number or without end.
   */
  [[nodiscard]] GEOXY_API auto PlaceAt(const Line &line, double distance) -> Place;

  /**
   *  The place at a distance along a ray from its start.
   *
   *  @param  ray      The ray.
   *  @param  distance The distance, negative before its start.
   *  @return          The place; for a ray of no way its start at a distance of nothing, and no
   *                   place, as a line's of no length, at any other.
   */
  [[nodiscard]] GEOXY_API auto PlaceAt(const Ray &ray, double distance) -> Place;

  /**
   *  The place at a distance along an xline from the point it is given by.
   *
   *  @param  xline    The xline.
   *  @param  distance The distance, negative behind the point.
   *  @return          The place; for an xline of no way its point at a distance of nothing, and no
   *                   place, as a line's of no length, at any other.
   */
  [[nodiscard]] GEOXY_API auto PlaceAt(const Xline &xline, double distance) -> Place;

  /**
   *  The place at a distance along an arc from its start.
   *
   *  @param  arc      The arc.
   *  @param  distance The distance, negative before its start, past its length beyond its end.
   *  @return          The place; for an arc of no length its start at a distance of nothing, and no
   *                   place, as a line's of no length, at any other.
   */
  [[nodiscard]] GEOXY_API auto PlaceAt(const Arc &arc, double distance) -> Place;

  /**
   *  The place at a distance along an ellipse from its start, by Newton's method on its length.
   *
   *  @param  ellipse  The ellipse.
   *  @param  distance The distance, negative before its start, past its length beyond its end.
   *  @return          The place, its length from the start the distance to a few steps of a double;
   *                   no place, as a line's of no length, for an ellipse that is a point.
   */
  [[nodiscard]] GEOXY_API auto PlaceAt(const Ellipse &ellipse, double distance) -> Place;

  /**
   *  The place at a distance along a polyline from its start.
   *
   *  @param  polyline The polyline.
   *  @param  distance The distance, negative before its start along its first segment, past its
   *                   length beyond its end along its last.
   *  @return          The place; the origin at 0 for no vertices; for one vertex, or a last segment
   *                   of no length, the vertex at a distance of nothing and no place, as a line's
   *                   of no length, at any other.
   */
  [[nodiscard]] GEOXY_API auto PlaceAt(const Polyline &polyline, double distance) -> Place;

  /**
   *  The place at a distance along a spline from its start, by Newton's method on its length: a
   *  polynomial one of degree 3 at most has its spans' lengths found once, and each step integrates
   *  the one span it falls in; any other has its strips' lengths found once, and the place within
   *  its domain is sought in the one strip it falls in, in that strip's own parameter, as finely as
   *  the spline's parameter can be told there.
   *
   *  Past its ends the spline runs on as its end spans do, for as far as doubles go. One with
   *  weights stops running where its weight comes to nothing, off to no point there and infinitely
   *  long to that place: every distance on that side has its place short of there, a double from it
   *  at the farthest.
   *
   *  @param  spline   The spline.
   *  @param  distance The distance, negative before its start, past its length beyond its end.
   *  @param  accuracy How far the length from the start to the place may be off; no finer than the
   *                   spline's speed is known, which near where it stops running is not fine.
   *  @return          The place; the origin at 0 for a spline that is not valid. No place, every
   *                   number of it not a number, where the spline never gets as far as the
   *                   distance: one that is a point, or one with weights that runs on toward a
   *                   point it nears for ever, its length coming to an end.
   */
  [[nodiscard]] GEOXY_API auto PlaceAt(const Spline &spline, double distance, Accuracy accuracy)
      -> Place;

  /**
   *  The place at a distance along a quadratic Bézier curve from its start.
   *
   *  @param  quadratic The curve.
   *  @param  distance  The distance, negative before its start, past its length beyond its end, as
   *                    far as doubles go.
   *  @param  accuracy  How far the length from the start to the place may be off.
   *  @return           The place; none, every number of it not a number, at any distance but
   *                    nothing for a curve that is a point.
   */
  [[nodiscard]] GEOXY_API auto PlaceAt(const Bezier::Quadratic &quadratic, double distance,
                                       Accuracy accuracy) -> Place;

  /**
   *  The place at a distance along a cubic Bézier curve from its start.
   *
   *  @param  cubic    The curve.
   *  @param  distance The distance, negative before its start, past its length beyond its end.
   *  @param  accuracy How far the length from the start to the place may be off.
   *  @return          The place; none, every number of it not a number, at any distance but nothing
   *                   for a curve that is a point.
   */
  [[nodiscard]] GEOXY_API auto PlaceAt(const Bezier::Cubic &cubic, double distance,
                                       Accuracy accuracy) -> Place;

  /**
   *  The place at a distance along any curve from its start.
   *
   *  @param  curve    The curve.
   *  @param  distance The distance.
   *  @param  accuracy How far the length to the place may be off, for a kind that is searched.
   *  @return          The place, as its own kind gives it.
   */
  [[nodiscard]] GEOXY_API auto PlaceAt(const Curve &curve, double distance, Accuracy accuracy)
      -> Place;
} // namespace Geoxy

#endif /* GEOXY_LENGTH_HPP */
