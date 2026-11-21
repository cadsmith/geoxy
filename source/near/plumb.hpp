/***************************************************************************************************
 *  @file      source/near/plumb.hpp
 *  @brief     The place of a spline nearest a point: its spans arched, sounded and pruned, and the
 *             feet on what is left weighed.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_NEAR_PLUMB_HPP
#define GEOXY_NEAR_PLUMB_HPP

#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <cstddef>
#include <optional>

namespace Geoxy::Near
{
  /**
   *  The place of a polynomial spline nearest a point, by branch and bound, as Least would choose
   *  it from every span's feet and every knot: the knots weighed first, then the spans taken in
   *  order of how near their boxes come, each looked into by halving until it is shown to lie
   *  farther than the nearest place so far, as every point of it At gives does, or solved for its
   *  feet as the whole spline would be, its feet weighed; until a span's box lies that far.
   *
   *  @param  spline The spline, valid and polynomial.
   *  @param  point  The point.
   *  @return        The place; nothing where a number lies past Widest, the degree is MostPoles or
   *                 more, or a foot lies by a knot.
   */
  auto Pruned(const Spline &spline, Point point) -> std::optional<Place>;

  /**
   *  The place of a rational spline nearest a point, by branch and bound, as Least would choose it
   *  from the feet Listed finds on every span and every knot: the knots weighed first, then the
   *  spans in order of how near their boxes come, each solved for its feet on its piece and its
   *  derivative's numerator as Listed solves it, until a span's box lies farther than the nearest
   *  place so far, as every point of it At gives does.
   *
   *  @param  spline The spline, valid and rational.
   *  @param  point  The point.
   *  @return        The place; nothing where a number lies past Widest, the weights lie more than
   *                 2^40 apart, the degree is MostPoles or more, a span is not even, where a paving
   *                 stands for the spline, or a foot lies by a knot, where Listed finds and keeps
   *                 its feet otherwise.
   */
  auto Rationally(const Spline &spline, Point point) -> std::optional<Place>;

  /**
   *  How much nearer a point the place of a circle at one angle from the way to the point is than
   *  the place at another, at the least: for angles a <= b from that way, within a half turn, the
   *  squares of the distances differ by 4 d r (sin²(b / 2) - sin²(a / 2)) = 4 d r sin((b - a) / 2)
   *  sin((a + b) / 2), d the point's distance from the center and r the radius, and the distances,
   *  each at most d + r, by at least that over 2 (d + r); sin x is 2 x / pi or more from nothing to
   *  a quarter turn.
   *
   *  @param  sizes d and r, both above nothing.
   *  @param  near  a, the nearer place's angle from the way, from nothing to a half turn.
   *  @param  far   b, the farther's, from a to a half turn.
   *  @return       The least difference of the distances.
   */
  auto CloserBy(Vector sizes, double near, double far) -> double;
} // namespace Geoxy::Near

#endif /* GEOXY_NEAR_PLUMB_HPP */
