/***************************************************************************************************
 *  @file      geoxy/arc.hpp
 *  @brief     Arcs of circles, and circles as arcs of a whole turn.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_ARC_HPP
#define GEOXY_ARC_HPP

#include "geoxy/angle.hpp"
#include "geoxy/api.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <functional>
#include <optional>

namespace Geoxy
{
  /**
   *  An arc of a circle: a start and a signed sweep about a center.
   *
   *  A positive sweep runs counterclockwise and a negative one clockwise, so an arc is never worked
   *  out from two ends: a sweep of a whole turn is a circle, and one of none is a point. Held as
   *  DXF's ARC is not, with two ends and no direction, so that no rounding of an end turns a whole
   *  circle into nothing.
   */
  struct Arc
  {
    Point  Center = {};  /**< The center.                                        */
    double Radius = 0.0; /**< The radius.                                        */
    Angle  Start  = {};  /**< Where it starts, from the x axis counterclockwise. */
    Angle  Sweep  = {};  /**< How far it runs, counterclockwise when positive.   */

    /**
     *  The arc from a start through a point to an end, ARC's three points.
     *
     *  Its center is the three's circumcenter, worked out from the start, so its radius is the
     *  start's distance from it; it runs the way the three turn, as Orient tells it exactly, and
     *  its angles are in radians. Its ends are where At puts them, within a few steps of a double
     *  of the points given.
     *
     *  @param  start   Where it starts.
     *  @param  through A point it passes between its ends.
     *  @param  end     Where it ends.
     *  @return         The arc; none where the three are on one line or two are one point, where
     *                  the cross product of their two legs from the start rounds to nothing, or
     *                  where its center is past every double.
     */
    [[nodiscard]] GEOXY_API static auto Through(Point start, Point through, Point end)
        -> std::optional<Arc>;

    /**
     *  The arc from a start about a center to an end's direction, ARC's start, center and end.
     *
     *  Counterclockwise, as AutoCAD draws it: its radius is the start's distance from the center,
     *  and it runs until the direction of the end, which need not lie on it. An end on the start's
     *  own direction makes a whole circle.
     *
     *  @param  start  Where it starts.
     *  @param  center The center.
     *  @param  end    A point in the direction it ends at.
     *  @return        The arc, its angles in radians; none where the start or the end is the
     *                 center.
     */
    [[nodiscard]] GEOXY_API static auto Around(Point start, Point center, Point end)
        -> std::optional<Arc>;

    /**
     *  The arc from a start to an end of a radius, ARC's start, end and radius.
     *
     *  Counterclockwise, as AutoCAD draws it: the minor arc for a positive radius, its center on
     *  the left of the way from the start to the end, and the major arc for a negative one.
     *
     *  @param  start  Where it starts.
     *  @param  end    Where it ends.
     *  @param  radius The radius, signed.
     *  @return        The arc, its radius the magnitude given and its angles in radians; none where
     *                 the start and end are one point or further apart than twice the radius.
     */
    [[nodiscard]] GEOXY_API static auto Spanning(Point start, Point end, double radius)
        -> std::optional<Arc>;

    /**
     *  The arc from a start, leaving it along a way, to an end, ARC's start, end and direction.
     *
     *  @param  start Where it starts.
     *  @param  end   Where it ends.
     *  @param  way   The way it leaves its start, of any length.
     *  @return       The arc, counterclockwise where the end is left of the way and clockwise where
     *                it is right, its angles in radians; none for a way of no length, where the end
     *                lies on the line along the way, the cross product of the way and the chord
     *                rounding to nothing, or where its center is past every double.
     */
    [[nodiscard]] GEOXY_API static auto Leaving(Point start, Point end, Vector way)
        -> std::optional<Arc>;

    /**
     *  Whether two arcs are equal exactly.
     *
     *  @param  lhs One arc.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Arc &lhs, const Arc &rhs) -> bool
    {
      return lhs.Center == rhs.Center
          && std::equal_to<double>{}(lhs.Radius, rhs.Radius)
          && lhs.Start == rhs.Start
          && lhs.Sweep == rhs.Sweep;
    }
  };

  /**
   *  A circle: the arc of a whole turn counterclockwise from the x axis.
   *
   *  @param  center The center.
   *  @param  radius The radius.
   *  @return        The arc.
   */
  [[nodiscard]] constexpr auto Circle(Point center, double radius) -> Arc
  {
    return {.Center = center, .Radius = radius, .Start = {}, .Sweep = Turn};
  }
} // namespace Geoxy

#endif /* GEOXY_ARC_HPP */
