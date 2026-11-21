/***************************************************************************************************
 *  @file      source/corner/choose.hpp
 *  @brief     Which circle of a radius rounds a corner of two curves, and which side of each is
 *             kept.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CORNER_CHOOSE_HPP
#define GEOXY_CORNER_CHOOSE_HPP

#include "corner/circle.hpp"
#include "corner/guide.hpp"

#include "geoxy/point.hpp"

#include <array>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy::Corner
{
  /**
   *  Two curves' guides, in order.
   */
  using Guides = std::array<const Guide *, 2>;

  /**
   *  Which side of a corner a curve keeps: the side its pick lies on, the corner's parameter taken
   *  as the one nearest the curve, where it repeats.
   */
  struct Side
  {
    double Corner = 0.0; /**< The corner's parameter.                            */
    int    Keep   = 0;   /**< 1 past it, -1 before it, 0 for a pick right on it. */
  };

  /**
   *  A circle and how the curves are kept from its touches.
   */
  struct Choice
  {
    Geoxy::Corner::Found Found = {};   /**< The circle, its parameters those the parts run to. */
    std::array<int, 2>   Keeps = {};   /**< The side kept of each touch.                       */
    bool                 Ahead = true; /**< Whether the arc runs counterclockwise.             */
    double               Score = 0.0;  /**< How far its touches are from what they are near.   */
    double               Sweep = 0.0;  /**< How far the arc turns, either way.                 */
  };

  /**
   *  The circles of a radius touching two curves, or their corners for the radius 0.
   *
   *  @param  guides The curves' guides.
   *  @param  radius The radius.
   *  @return        Each, with its touches' parameters.
   */
  auto Circles(const Guides &guides, double radius) -> std::vector<Found>;

  /**
   *  The sign of a number.
   *
   *  @param  value The number.
   *  @return       1, -1, or 0 for nothing.
   */
  auto Sign(double value) -> int;

  /**
   *  Which side of a corner a curve keeps.
   *
   *  @param  guide  The curve's guide.
   *  @param  corner The corner's parameter, below one period where it repeats.
   *  @return        The side: on a closed curve the nearer way round to the pick; on an arc the
   *                 corner past the nearer end where it is off the arc.
   */
  auto Sided(const Guide &guide, double corner) -> Side;

  /**
   *  The sweep of an arc from one touch to another about a center, the way it runs.
   *
   *  @param  center  The center.
   *  @param  touches Where it starts and where it ends.
   *  @param  ahead   Whether it runs counterclockwise.
   *  @return         In radians, in (0, 2 pi] counterclockwise and [-2 pi, 0) clockwise.
   */
  auto Swept(Point center, const std::array<Point, 2> &touches, bool ahead) -> double;

  /**
   *  Which way the arc of a circle runs, where it runs smoothly from the first curve's part kept
   *  onto the second's: arriving at its first touch along the first curve, toward it from the part
   *  kept, and leaving its second along the second, toward the part kept there.
   *
   *  @param  guides The curves' guides.
   *  @param  found  The circle, its parameters those the parts run to.
   *  @param  keeps  The side kept of each touch.
   *  @return        Whether it runs counterclockwise; none where it would turn back at a touch.
   */
  auto Flowing(const Guides &guides, const Found &found, const std::array<int, 2> &keeps)
      -> std::optional<bool>;

  /**
   *  Whether one choice is better than another: nearer what its touches should be near, then of the
   *  shorter arc.
   *
   *  @param  one   A choice.
   *  @param  other The best so far, if any.
   *  @return       Whether one is better.
   */
  auto Better(const Choice &one, const std::optional<Choice> &other) -> bool;

  /**
   *  A circle, where it runs smoothly from one curve onto the other with the sides kept, as a
   *  choice.
   *
   *  @param  guides The curves' guides.
   *  @param  found  The circle, its parameters those the parts run to.
   *  @param  keeps  The side kept of each touch.
   *  @return        The choice, scored by nothing yet; none where it does not run smoothly.
   */
  auto Chosen(const Guides &guides, const Found &found, const std::array<int, 2> &keeps)
      -> std::optional<Choice>;

  /**
   *  A corner of two curves, and every corner they meet at: those past it bound the parts kept.
   */
  struct Nook
  {
    const Found           *At  = nullptr; /**< The corner rounded.          */
    std::span<const Found> All = {};      /**< Every corner, it among them. */
  };

  /**
   *  The circle in the angle between the parts of two curves kept from their corner, nearest it.
   *
   *  The parts kept run from the corner the way each pick lies, up to the next corner where the
   *  curves themselves meet that way, not their run ons, past which the angle is another's.
   *
   *  @param  guides  The curves' guides.
   *  @param  nook    The corner, and every corner.
   *  @param  circles The circles of the radius touching both.
   *  @return         The circle whose touches lie on the parts kept and whose arc runs smoothly
   *                  from one onto the other, nearest the corner; none where none does.
   */
  auto Nooked(const Guides &guides, const Nook &nook, const std::vector<Found> &circles)
      -> std::optional<Choice>;

  /**
   *  The ways a curve may be kept from a touch: the part from it the pick lies on, either way round
   *  where its parameter repeats.
   *
   *  @param  guide The curve's guide.
   *  @param  along The touch's parameter.
   *  @return       Each way, as the touch's parameter the part runs from or to and the side kept.
   */
  auto Options(const Guide &guide, double along) -> std::vector<std::pair<double, int>>;

  /**
   *  The circle touching two curves that never meet, nearest the picks, whose arc runs smoothly
   *  from the part of one kept onto the other's.
   *
   *  @param  guides  The curves' guides.
   *  @param  picks   Where each was picked.
   *  @param  circles The circles of the radius touching both.
   *  @return         The circle; none where no arc runs smoothly.
   */
  auto Bridging(const Guides &guides, const std::array<Point, 2> &picks,
                const std::vector<Found> &circles) -> std::optional<Choice>;

  /**
   *  The corner of two curves a fillet rounds: one where they meet on both curves themselves, not
   *  on what they run on to, before one on either's run on, and of those the nearest their picks.
   *
   *  @param  guides  The curves' guides.
   *  @param  corners The corners.
   *  @param  picks   The picks.
   *  @return         The corner; none for no corner.
   */
  auto Closest(const Guides &guides, const std::vector<Found> &corners,
               const std::array<Point, 2> &picks) -> std::optional<Found>;

  /**
   *  A circle chosen, settled: its center and touches the doubles nearest the exact ones where both
   *  curves are lines or circles, and each touch's parameter worked out again from its point, the
   *  one nearest the parameter chosen.
   *
   *  @param  guides The curves' guides.
   *  @param  found  The circle, as chosen.
   *  @param  radius The radius it was found for.
   *  @return        It, settled.
   */
  auto Placed(const Guides &guides, const Found &found, double radius) -> Found;
} // namespace Geoxy::Corner

#endif /* GEOXY_CORNER_CHOOSE_HPP */
