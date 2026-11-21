/***************************************************************************************************
 *  @file      source/corner/smooth.hpp
 *  @brief     The circles of a radius touching two curves one of which is an ellipse, a spline or a
 *             Bézier curve, found on their offsets by Newton's method; and where such curves meet.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CORNER_SMOOTH_HPP
#define GEOXY_CORNER_SMOOTH_HPP

#include "corner/circle.hpp"
#include "corner/guide.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/place.hpp"

#include <array>
#include <optional>
#include <vector>

namespace Geoxy::Corner
{
  /**
   *  The circles of a radius touching two curves, at least one of them smooth, or where the two
   *  meet for the radius 0.
   *
   *  A line is taken whole and an arc as its circle; an ellipse as its whole ellipse; a spline or a
   *  Bézier curve over its domain. The radius 0 asks where they meet, as Intersect finds it.
   *  Otherwise each curve's offsets by the radius either side are laid as chords through points
   *  every few steps of its parameter, a line's and a circle's exactly; where two offsets' chords
   *  cross, the center of a circle touching both lies near, and Newton's method on the two offsets'
   *  points, their derivatives the curves' times one less the radius by the curvature, finds it to
   *  a few steps of a double.
   *
   *  @param  first  One curve's guide.
   *  @param  second The other's.
   *  @param  radius The radius, not below nothing.
   *  @return        The circles or the corners, each touch with its parameter; none settled again
   *                 exactly.
   */
  auto Smoothly(const Guide &first, const Guide &second, double radius) -> std::vector<Found>;

  /**
   *  Two smooth curves' offsets: the curves, and each one's signed distance, positive to its left.
   */
  struct Offsets
  {
    std::array<const Guide *, 2> Guides    = {}; /**< The curves' guides.  */
    std::array<double, 2>        Distances = {}; /**< Each one's distance. */
  };

  /**
   *  Where two smooth curves' offsets meet, from a guess, by Newton's method on c1 + s1 r n1 = c2 +
   *  s2 r n2, each offset's derivative c' (1 - s r k) for k its curve's signed curvature; a step
   *  across offsets that run alike divided by the least normal double, not by nothing.
   *
   *  @param  offsets The curves and their distances.
   *  @param  guess   The parameters to start from.
   *  @return         The parameters where they meet, the last step under 2^-50 of the size they
   *                  were stepped from and the offsets' points within 2^-30 of the first curve's
   *                  size at the guess; none after 200 steps that do not settle, as where a curve
   *                  has no way.
   */
  auto Polished(const Offsets &offsets, std::array<double, 2> guess)
      -> std::optional<std::array<double, 2>>;

  /**
   *  A curve as Intersect meets it for a corner: a line's or a circle's carrier whole, an ellipse
   *  whole, a spline or a Bézier curve as it is.
   *
   *  @param  guide The curve's guide.
   *  @return       The curve to meet.
   */
  auto Carried(const Guide &guide) -> Curve;

  /**
   *  A place's parameter on a curve, from where Intersect found it on what Carried gave.
   *
   *  @param  guide The curve's guide.
   *  @param  place Where on what was met.
   *  @return       The parameter on the curve: an ellipse's scaled from its whole turn.
   */
  auto Cornered(const Guide &guide, const Place &place) -> double;
} // namespace Geoxy::Corner

#endif /* GEOXY_CORNER_SMOOTH_HPP */
