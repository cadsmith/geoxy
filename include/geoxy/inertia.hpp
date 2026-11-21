/***************************************************************************************************
 *  @file      geoxy/inertia.hpp
 *  @brief     The second moments of what a curve encloses, about its centroid, and their principal
 *             axes, as MASSPROP gives them.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_INERTIA_HPP
#define GEOXY_INERTIA_HPP

#include "geoxy/accuracy.hpp"
#include "geoxy/api.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <optional>

namespace Geoxy
{
  /**
   *  What a curve encloses, weighed: its area, its centroid, and its second moments about the
   *  centroid, each signed as the area is, negative for a curve that runs clockwise.
   */
  struct Moments
  {
    double Area     = 0.0; /**< The signed area.                                       */
    Point  Centroid = {};  /**< The centroid.                                          */
    double AboutX   = 0.0; /**< The moment about the x axis through it: of (y - cy)^2. */
    double AboutY   = 0.0; /**< The moment about the y axis through it: of (x - cx)^2. */
    double Product  = 0.0; /**< The product of inertia: of (x - cx) (y - cy).          */
  };

  /**
   *  The principal moments about a centroid, and the axis of the least.
   */
  struct Principal
  {
    double Least = 0.0; /**< The least moment about an axis through the centroid. */
    double Most  = 0.0; /**< The most, about the axis square to it.               */
    Vector Axis  = {};  /**< The way of the axis of the least, a unit vector.     */
  };

  /**
   *  The area, centroid and second moments about it of what a curve encloses, closed by the line
   *  from its end to its start as Area closes it.
   *
   *  The straight part of a polyline is summed by the polygon's sums about its first vertex. Every
   *  other part, a bulge's, an arc's, an ellipse's, a Bézier curve's or a spline span's, is taken
   *  with its chord, in the frame of the chord, where Green's theorem gives the chord nothing: each
   *  moment is the integral along the part alone of a form with a factor of the part's distance
   *  from its chord, by Gauss-Legendre rules halved until two agree, so nothing cancels however
   *  flat the part. All of it is summed about the first point, and moved to the centroid once.
   *
   *  @param  curve    The curve.
   *  @param  accuracy How far the area may be off, and each moment as far times the curve's size to
   *                   its power; at or below 0, as near as two rules agreeing to 2^-44 allows.
   *  @return          The moments; none for a curve that encloses no area, or never ends.
   */
  [[nodiscard]] GEOXY_API auto Inertia(const Curve &curve, Accuracy accuracy)
      -> std::optional<Moments>;

  /**
   *  The principal moments of second moments about a centroid: the least and most moments about an
   *  axis through it, and the way of the axis of the least, along which the area spreads most.
   *
   *  @param  moments The moments, those of a clockwise curve taken negated.
   *  @return         The principal moments and the axis.
   */
  [[nodiscard]] GEOXY_API auto Axes(const Moments &moments) -> Principal;
} // namespace Geoxy

#endif /* GEOXY_INERTIA_HPP */
