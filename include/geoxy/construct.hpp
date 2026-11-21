/***************************************************************************************************
 *  @file      geoxy/construct.hpp
 *  @brief     Circles and lines built from what they pass through or touch: CIRCLE's 2P, 3P, TTR
 *             and TTT, and LINE's tangents.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CONSTRUCT_HPP
#define GEOXY_CONSTRUCT_HPP

#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/segment.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace Geoxy
{
  /**
   *  How a circle built tangent to a curve lies against it.
   */
  enum class Qualifier : std::uint8_t
  {
    Left,      /**< Left of a line, as it runs from its start to its end. */
    Right,     /**< Right of a line.                                      */
    Outside,   /**< Outside a circle, the two touching from without.      */
    Enclosing, /**< Around a circle, holding it, touching it from within. */
    Enclosed,  /**< Within a circle, held by it, touching it from within. */
  };

  /**
   *  Where a circle built tangent to a curve touches it, and how it lies against it.
   */
  struct Grazing
  {
    Point     At   = {};                 /**< Where it touches the curve, on it. */
    Qualifier Kind = Qualifier::Outside; /**< How it lies against it.            */
  };

  /**
   *  A circle tangent to curves, with where it touches each and how it lies against it.
   *
   *  @tparam CountT How many curves it touches.
   */
  template <std::size_t CountT>
  struct Tangential
  {
    Arc                         Circle  = {}; /**< The circle, a whole turn from the x axis.    */
    std::array<Grazing, CountT> Touches = {}; /**< Where it touches each curve, in their order. */
  };

  /**
   *  The circle on a diameter, CIRCLE's two points.
   *
   *  @param  first  One end of the diameter.
   *  @param  second The other.
   *  @return        The circle about their midpoint through both; of no radius where they are one
   *                 point.
   */
  [[nodiscard]] GEOXY_API auto Circle(Point first, Point second) -> Arc;

  /**
   *  The circle through three points, CIRCLE's three points.
   *
   *  @param  first  One point.
   *  @param  second Another.
   *  @param  third  The third.
   *  @return        The circle of Arc::Through's center and radius; none where it has none.
   */
  [[nodiscard]] GEOXY_API auto Circle(Point first, Point second, Point third) -> std::optional<Arc>;

  /**
   *  Every circle of a radius tangent to two curves, CIRCLE's tangent, tangent, radius.
   *
   *  A line is taken whole, as the line it lies on, and an arc as its whole circle; whether a touch
   *  falls on the curve itself is the caller's to ask of the point given. Each circle's center
   *  solves Apollonius's equations, one for each curve and the side it lies on: a line's offset by
   *  the radius, a circle's grown or shrunk by it, the two circles' difference a line.
   *
   *  @param  first  One curve.
   *  @param  second The other.
   *  @param  radius The radius, above nothing.
   *  @return        Every circle, up to four for two lines and eight for a line and a circle or two
   *                 circles, in no order but that of their qualifiers. None for a line of no
   *                 length, for parallel lines, which have none or a whole line of them, or for a
   *                 radius not above nothing; a circle of no radius is a point to pass through.
   */
  [[nodiscard]] GEOXY_API auto TangentCircles(const Segment &first, const Segment &second,
                                              double radius) -> std::vector<Tangential<2>>;

  /**
   *  Every circle tangent to three curves, CIRCLE's tangent, tangent, tangent: Apollonius's
   *  problem.
   *
   *  Each line is taken whole and each arc as its whole circle. For each side of each curve the
   *  equations are linear in the center and the radius but for one circle's, which the others',
   *  less it, leave linear: two unknowns then follow the third, the one the equations are best
   *  conditioned to leave free, and one quadratic gives it.
   *
   *  @param  first  One curve.
   *  @param  second Another.
   *  @param  third  The third.
   *  @return        Every circle, up to four for three lines and eight otherwise, each with where
   *                 it touches the three and how it lies against them. None for a line of no
   *                 length, and none for a choice of sides whose equations do not fix one: three
   *                 lines of which two are parallel give two circles, not four.
   */
  [[nodiscard]] GEOXY_API auto TangentCircles(const Segment &first, const Segment &second,
                                              const Segment &third) -> std::vector<Tangential<3>>;

  /**
   *  The lines from a point tangent to a circle, LINE's tangent snap.
   *
   *  @param  from   The point.
   *  @param  circle The circle, an arc taken as its whole circle.
   *  @return        The lines from the point to where they touch it, the one on the circle's left
   *                 as seen from the point first: two from outside; one of no length, at the point,
   *                 from a point on it; none from within, or from the center of a circle of no
   *                 radius.
   */
  [[nodiscard]] GEOXY_API auto TangentLines(Point from, const Arc &circle) -> std::vector<Line>;

  /**
   *  The lines tangent to two circles, LINE's tangent to tangent.
   *
   *  A line touches both where the circles' centers lie at its distance of their radii from it, on
   *  one side for an outer tangent and on either side for an inner one: its unit normal n has n ·
   *  (c2 - c1) = r2 - r1 or -(r1 + r2), and turns either way from the centers' line by as much.
   *
   *  @param  first  One circle, an arc taken as its whole circle.
   *  @param  second The other.
   *  @return        The two outer tangents, then the two inner ones, each from where it touches the
   *                 first circle to where it touches the second; a pair that coincides, as the
   *                 inner ones of circles that touch, given twice; a pair that does not exist, as
   *                 the inner ones of circles that cross, left out. None for concentric circles.
   */
  [[nodiscard]] GEOXY_API auto TangentLines(const Arc &first, const Arc &second)
      -> std::vector<Line>;
} // namespace Geoxy

#endif /* GEOXY_CONSTRUCT_HPP */
