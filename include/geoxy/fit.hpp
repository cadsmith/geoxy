/***************************************************************************************************
 *  @file      geoxy/fit.hpp
 *  @brief     Splines made: through points, from a polyline's vertices, with knots inserted or a
 *             degree raised; and a spline made a polyline again.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_FIT_HPP
#define GEOXY_FIT_HPP

#include "geoxy/accuracy.hpp"
#include "geoxy/api.hpp"
#include "geoxy/error.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <cstdint>
#include <expected>
#include <optional>
#include <span>

namespace Geoxy
{
  /**
   *  How a fit spaces its knots, as AutoCAD's SPLINE asks: each step from one point's knot to the
   *  next is the chord between them, its square root, or one.
   */
  enum class Knotting : std::uint8_t
  {
    Chord,   /**< The chord's length: unit speed along the chords, AutoCAD's own rule. */
    Root,    /**< The square root of the chord's length, centripetal.                  */
    Uniform, /**< One for every chord.                                                 */
  };

  /**
   *  How a spline is fit through points: its knot rule, a tangent at either end if one is asked,
   *  and whether it closes.
   */
  struct Fitting
  {
    Knotting              Knots  = Knotting::Chord; /**< How the knots are spaced.              */
    std::optional<Vector> Start  = {};              /**< The way it leaves its first point.     */
    std::optional<Vector> End    = {};              /**< The way it arrives at its last point.  */
    bool                  Closed = false;           /**< Whether it runs on from last to first. */
  };

  /**
   *  The cubic spline through points, as AutoCAD's SPLINE fits one.
   *
   *  Its knots are the points' parameters, unscaled: the running sum of the chords, of their square
   *  roots, or of ones, from 0 at the first point. Open, it is clamped, so that n + 1 points make n
   *  + 3 poles, and its first and last poles are its first and last points exactly. An end with a
   *  tangent leaves along it at the speed of the chord next to it over that chord's knot step, so
   *  at unit speed for chord knots; an end without one is natural, its second derivative nothing.
   *  Closed, it is periodic: its knots run on by the loop's length either side, its first three
   *  poles repeated at its end, and a last point the same as the first is the first.
   *
   *  The poles solve a tridiagonal system, cyclic when closed, about the first point, and once more
   *  for what the first answer leaves, worked out in twice the precision: the spline meets every
   *  point to a few steps of a double at its size.
   *
   *  @param  points  The points, in order.
   *  @param  fitting The knot rule, the end tangents, and whether it closes.
   *  @return         The spline; Degenerate for fewer than two points, or three when closed, two
   *                  points in a row the same or so near that their knots are, or a tangent of no
   *                  length; NotFinite for a point or a tangent that is not finite, or knots or
   *                  poles past the doubles.
   */
  [[nodiscard]] GEOXY_API auto Fit(std::span<const Point> points, const Fitting &fitting)
      -> std::expected<Spline, Error>;

  /**
   *  A spline with a knot inserted, the same curve: every new pole a blossom of the span that holds
   *  it, weights and all.
   *
   *  @param  spline The spline.
   *  @param  knot   The knot, within its domain.
   *  @param  times  How many times to insert it; none leaves the spline as it is.
   *  @return        The spline with the knot; none for a spline that is not valid, a knot outside
   *                 the domain or not finite, a negative count, or a knot that would then be
   *                 repeated more than the degree.
   */
  [[nodiscard]] GEOXY_API auto Inserted(const Spline &spline, double knot, int times)
      -> std::optional<Spline>;

  /**
   *  A spline of a higher degree, the same curve: each pole the mean of the blossoms with one value
   *  left out, every knot within the domain repeated once more for each degree, clamped at the
   *  domain's ends.
   *
   *  @param  spline The spline.
   *  @param  times  How many degrees to raise it by.
   *  @return        The spline; none for a spline that is not valid, a negative count, or a degree
   *                 past MaxDegree.
   */
  [[nodiscard]] GEOXY_API auto Elevated(const Spline &spline, int times) -> std::optional<Spline>;

  /**
   *  The spline whose control frame is a polyline's vertices, as PEDIT Spline makes it: uniform
   *  knots, clamped at its ends when open, periodic when closed, the bulges left out.
   *
   *  @param  polyline The polyline.
   *  @param  degree   The degree: 2 for AutoCAD's quadratic, 3 for its cubic.
   *  @return          The spline; none for a degree outside 1 to MaxDegree, or for a polyline of no
   *                   more vertices than the degree.
   */
  [[nodiscard]] GEOXY_API auto Framed(const Polyline &polyline, int degree)
      -> std::optional<Spline>;

  /**
   *  A spline as a polyline of chords, each within the accuracy of it, as SPLINEDIT converts one.
   *
   *  @param  spline   The spline.
   *  @param  accuracy How far the chords may stray from it.
   *  @return          The polyline, closed where the spline ends where it starts; no vertices for a
   *                   spline that is not valid.
   */
  [[nodiscard]] GEOXY_API auto Explode(const Spline &spline, Accuracy accuracy) -> Polyline;
} // namespace Geoxy

#endif /* GEOXY_FIT_HPP */
