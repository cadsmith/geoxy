/***************************************************************************************************
 *  @file      source/corner/circle.hpp
 *  @brief     The circles of a radius touching two lines or circles: each found in doubles to
 *             choose among them, then the one chosen worked out as the doubles nearest its center
 *             and its touches.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CORNER_CIRCLE_HPP
#define GEOXY_CORNER_CIRCLE_HPP

#include "corner/forms.hpp"

#include "geoxy/point.hpp"

#include <array>
#include <optional>
#include <vector>

namespace Geoxy::Corner
{
  /**
   *  A line or a circle a fillet touches.
   */
  struct Carrier
  {
    bool     Linear = true; /**< Whether it is a line.  */
    Straight Line   = {};   /**< The line, where it is. */
    Round    Circle = {};   /**< The circle, where not. */
  };

  /**
   *  Two lines or circles and a radius, scaled by a power of two so that their largest number is
   *  near 1, where no product the forms take overflows or is lost below the doubles.
   */
  struct Problem
  {
    std::array<Carrier, 2> Carriers = {};  /**< The two, scaled.                      */
    double                 Radius   = 0.0; /**< The radius, scaled.                   */
    int                    Power    = 0;   /**< The power of two they were scaled by. */
  };

  /**
   *  A circle touching two curves: its center and where it touches each.
   */
  struct Touch
  {
    Point                Center  = {}; /**< The center.                      */
    std::array<Point, 2> Touches = {}; /**< Where it touches each, in order. */
  };

  /**
   *  A circle of the radius touching two curves, or a corner where they meet, as found: where, the
   *  parameter of each touch on its curve, and how to work it out again exactly.
   */
  struct Found
  {
    Geoxy::Corner::Touch  Touch  = {};   /**< The center and the touches.                 */
    std::array<double, 2> Alongs = {};   /**< Each touch's parameter on its curve.        */
    Geoxy::Corner::Ask    Ask    = {};   /**< Which circle, for lines and circles.        */
    bool                  Exact  = true; /**< Whether Settled works it out again exactly. */
  };

  /**
   *  Two lines or circles and a radius, scaled.
   *
   *  @param  carriers The two, as given.
   *  @param  radius   The radius.
   *  @return          The problem.
   */
  auto Posed(const std::array<Carrier, 2> &carriers, double radius) -> Problem;

  /**
   *  Every way a circle of the problem's radius may lie against its two curves: the sides of each,
   *  and the root, the radius 0 asking only for the roots.
   *
   *  @param  problem The problem.
   *  @return         The asks, each a circle where its radicand is not below nothing.
   */
  auto Asks(const Problem &problem) -> std::vector<Ask>;

  /**
   *  A circle of the problem, in doubles: near enough to choose among them by.
   *
   *  @param  problem The problem.
   *  @param  ask     Which circle.
   *  @return         It, as given; none where its radicand is below nothing, as sure arithmetic or
   *                  exact tells it, or where it would divide by nothing.
   */
  auto Approximated(const Problem &problem, const Ask &ask) -> std::optional<Touch>;

  /**
   *  A circle of the problem, every coordinate the double nearest the exact one: certified in sure
   *  double-doubles, or found by exact comparisons where no bound settles it.
   *
   *  @pre            The circle is: Approximated gives it.
   *  @param  problem The problem.
   *  @param  ask     Which circle.
   *  @return         It, as given.
   */
  auto Settled(const Problem &problem, const Ask &ask) -> Touch;

  /**
   *  A circle of a radius found with it: the half circle joining a point to a line.
   */
  struct Span
  {
    Geoxy::Corner::Touch Touch  = {};  /**< The center, the point and its foot. */
    double               Radius = 0.0; /**< The radius.                         */
  };

  /**
   *  The half circle from a point to a line, as FILLET joins parallel lines: the point's foot on
   *  the line, the center halfway, and half their distance, each the double nearest.
   *
   *  @param  line The line.
   *  @param  from The point.
   *  @return      The half circle's center, the point and the foot, and its radius.
   */
  auto Bridged(const Straight &line, Point from) -> Span;

  /**
   *  How CHAMFER cuts two lines back from their corner: the distance along each, or along the first
   *  and the angle's sine and cosine, and the side of the corner each keeps.
   */
  struct Cutback
  {
    std::array<double, 2>                Distances = {}; /**< Along the first and the second. */
    std::array<int, 2>                   Keeps     = {}; /**< The side kept of each, 1 or -1. */
    std::optional<std::array<double, 2>> Slant     = {}; /**< The angle's sine and cosine; where
                                                              given, the second distance is not
                                                              read.                          */
    int Across = 1;                                      /**< The sign of the ways' cross
                                                              product, worked out within.    */
  };

  /**
   *  Where CHAMFER cuts two lines, each the double nearest the exact point.
   *
   *  With X the corner, d1, d2 the ways, n1, n2 their squares and c = d1 x d2, the first is cut at
   *  X + k1 e1 d1 / sqrt(n1); the second at X + k2 e2 d2 / sqrt(n2), or with an angle a, at X + k2
   *  d2 e1 sin(a) sqrt(n1) / (sin(a) (-k1 k2 d1 . d2) + cos(a) |c|), the triangle's sine rule.
   *
   *  @param  first  One line.
   *  @param  second The other, not parallel to it.
   *  @param  cut    The distances or the distance and the angle, and the sides kept.
   *  @return        The two cuts; none where the angle's line never meets the second's part kept.
   */
  auto Bevelled(const Straight &first, const Straight &second, const Cutback &cut)
      -> std::optional<std::array<Point, 2>>;

  /**
   *  Whether two lines are parallel, exactly.
   *
   *  @param  first  One.
   *  @param  second The other.
   *  @return        Whether their ways' cross product is nothing.
   */
  auto Parallel(const Straight &first, const Straight &second) -> bool;
} // namespace Geoxy::Corner

#endif /* GEOXY_CORNER_CIRCLE_HPP */
