/***************************************************************************************************
 *  @file      source/tangency.hpp
 *  @brief     Whether a line and a circle, or two circles, touch exactly, and where: decided by
 *             their own short exact sums, the circle's quadratic of the second degree, not a
 *             conic's of the fourth; and where such two near tangent come nearest, rounded once.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TANGENCY_HPP
#define GEOXY_TANGENCY_HPP

#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/point.hpp"

#include <array>
#include <optional>

namespace Geoxy::Crossings
{
  /**
   *  What the short way decided about two curves' carriers touching.
   */
  struct Touching
  {
    bool                 Decided = false; /**< Whether it decided. */
    std::optional<Point> At      = {};    /**< Where they touch, rounded; none where they do
                                               not.                                         */
  };

  /**
   *  Whether two curves' carriers touch exactly, and the double nearest where: a line, a ray or an
   *  endless line and an arc's circle, the circle's own quadratic |P + t d - O|^2 = r^2 with its
   *  two roots one, (d x (P - O))^2 = r^2 |d|^2, touching at the foot of the center; two arcs'
   *  circles, their centers as far apart as the radii's sum or difference s, |g|^2 = s^2, touching
   *  at c1 + g r1 / s. Each sum is held exactly in a few doubles; the point is the double nearest
   *  the true one where its double-double quotient settles it, as nearly always. The same point the
   *  exact crossing finds, more cheaply.
   *
   *  @param  curves The two curves.
   *  @return        Whether they touch and where; undecided for any other curves, where a straight
   *                 curve has no length or two circles one center, or where the point does not
   *                 settle, for the exact crossing to decide.
   */
  auto Tangency(const std::array<const Curve *, 2> &curves) -> Touching;

  /**
   *  Where two curves' carriers cross near a hit, a line, a ray or an endless line and an arc's
   *  circle, or two arcs' circles, by their own short exact sums, as the exact crossing finds it:
   *  each coordinate of each root the double nearest it where its double-double estimate settles
   *  it, the root's share and the rest summed through their conjugate where they cancel; the root
   *  nearer the hit, the first where they are as near.
   *
   *  @param  curves The two curves.
   *  @param  near   Where the hit is.
   *  @return        The crossing; undecided for any other curves, where the carriers do not cross,
   *                 or a coordinate does not settle, for the exact crossing to decide.
   */
  auto Shortly(const std::array<const Curve *, 2> &curves, Point near) -> Touching;

  /**
   *  Where two circles near tangent come nearest, or cross deepest: halfway between each one's
   *  point on the line of their centers, c1 + w1 r1 u and c2 + w2 r2 u for u the way from the first
   *  center to the second and each w 1 or -1, which is (c1 + c2) / 2 + k g / (2 |g|), g = c2 - c1,
   *  k = w1 r1 + w2 r2: each coordinate ((c1 + c2) |g|^2 + k g sqrt(|g|^2)) / (2 |g|^2), worked out
   *  in expansions from the arcs' own numbers and rounded once, the double nearest it. So the same
   *  point to the bit whichever circle is first, the ways and the way between swapping with them,
   *  and turned, mirrored or scaled by a power of two with the plane, as the point two bulges of
   *  the same circles touch at is; worked out in doubles from the centers' gap, it was as much as
   *  381 steps of a double off where the centers lie near each other.
   *
   *  @pre         The centers differ.
   *  @param  arcs The two arcs, whose circles.
   *  @param  ways w1 and w2: which way along u each circle's point lies from its center.
   *  @return      The point.
   */
  auto Nearing(const std::array<const Arc *, 2> &arcs, const std::array<double, 2> &ways) -> Point;

  /**
   *  Where a line, a ray or an endless line near tangent to a circle comes nearest it, or crosses
   *  deepest: halfway between the foot of the center on the carrier and the circle's point toward
   *  it. For the carrier P + t d and g = P - c, the foot lies at c + w / |d|^2 and the circle's
   *  point at c + r w / (|d| |d x g|), w = g |d|^2 - d (g . d): each coordinate of the point
   *  halfway, ((2 c |d|^2 + w) |d x g| + r w sqrt(|d|^2)) / (2 |d|^2 |d x g|), worked out in
   *  expansions from the curves' own numbers, a line's way the exact difference of its ends, and
   *  rounded once, the double nearest it, as the point where the line meets the same circle as a
   *  bulge is.
   *
   *  @pre             The straight curve is a line, a ray or an endless line of some length, whose
   *                   carrier misses the center.
   *  @param  straight The line, ray or endless line.
   *  @param  arc      The arc, whose circle.
   *  @return          The point.
   */
  auto Nearing(const Curve &straight, const Arc &arc) -> Point;
} // namespace Geoxy::Crossings

#endif /* GEOXY_TANGENCY_HPP */
