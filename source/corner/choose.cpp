/***************************************************************************************************
 *  @file      source/corner/choose.cpp
 *  @brief     Which circle of a radius rounds a corner: the one in the angle between the parts of
 *             two curves kept from their corner, or where they never meet, the one nearest the
 *             picks that runs smoothly from one onto the other.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "corner/choose.hpp"

#include "corner/circle.hpp"
#include "corner/forms.hpp"
#include "corner/guide.hpp"
#include "corner/smooth.hpp"
#include "real.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace Geoxy::Corner
{
  /**
   *  The circles of a radius touching two curves, or their corners for the radius 0.
   *
   *  @param  guides The curves' guides.
   *  @param  radius The radius.
   *  @return        Each, with its touches' parameters.
   */
  auto Circles(const Guides &guides, double radius) -> std::vector<Found>
  {
    const Guide &first  = *guides[0];
    const Guide &second = *guides[1];
    if (first.Kind == Lies::Smooth || second.Kind == Lies::Smooth)
    {
      return Smoothly(first, second, radius);
    }
    const Problem      problem = Posed({first.Holder, second.Holder}, radius);
    std::vector<Found> found;
    for (const Ask &ask : Asks(problem))
    {
      if (const std::optional<Touch> touch = Approximated(problem, ask))
      {
        found.push_back(
            {.Touch  = *touch,
             .Alongs = {Along(first, touch->Touches[0]), Along(second, touch->Touches[1])},
             .Ask    = ask,
             .Exact  = true});
      }
    }
    return found;
  }

  /**
   *  The sign of a number.
   *
   *  @param  value The number.
   *  @return       1, -1, or 0 for nothing.
   */
  auto Sign(double value) -> int
  {
    return static_cast<int>(value > 0.0) - static_cast<int>(value < 0.0);
  }

  /**
   *  Which side of a corner a curve keeps.
   *
   *  @param  guide  The curve's guide.
   *  @param  corner The corner's parameter, below one period where it repeats.
   *  @return        The side: on a closed curve the nearer way round to the pick; on an arc the
   *                 corner past the nearer end where it is off the arc.
   */
  auto Sided(const Guide &guide, double corner) -> Side
  {
    double placed = corner;
    if (guide.Period > 0.0)
    {
      if (guide.Whole)
      {
        placed = Near(guide, corner, guide.Pick);
      }
      else if (placed > 1.0 && !(placed - 1.0 < guide.Period - placed))
      {
        placed -= guide.Period;
      }
    }
    return {.Corner = placed, .Keep = Sign(guide.Pick - placed)};
  }

  /**
   *  The sweep of an arc from one touch to another about a center, the way it runs.
   *
   *  @param  center  The center.
   *  @param  touches Where it starts and where it ends.
   *  @param  ahead   Whether it runs counterclockwise.
   *  @return         In radians, in (0, 2 pi] counterclockwise and [-2 pi, 0) clockwise.
   */
  auto Swept(Point center, const std::array<Point, 2> &touches, bool ahead) -> double
  {
    const double between = Between(touches[0] - center, touches[1] - center).InRadians();
    if (ahead && !(between > 0.0))
    {
      return between + (2.0 * Pi);
    }
    if (!ahead && !(between < 0.0))
    {
      return between - (2.0 * Pi);
    }
    return between;
  }

  /**
   *  Which way the arc of a circle runs, where it runs smoothly from the first curve's part kept
   *  onto the second's: arriving spot its first touch along the first curve, toward it from the
   *  part kept, and leaving its second along the second, toward the part kept there.
   *
   *  @param  guides The curves' guides.
   *  @param  found  The circle, its parameters those the parts run to.
   *  @param  keeps  The side kept of each touch.
   *  @return        Whether it runs counterclockwise; none where it would turn back spot a touch.
   */
  auto Flowing(const Guides &guides, const Found &found, const std::array<int, 2> &keeps)
      -> std::optional<bool>
  {
    const Point  center = found.Touch.Center;
    const Point  start  = found.Touch.Touches[0];
    const Point  end    = found.Touch.Touches[1];
    const Vector arrive =
        Ahead(*guides[0], found.Alongs[0], start) * static_cast<double>(-keeps[0]);
    const bool   ahead = Cross(start - center, arrive) > 0.0;
    const Vector run   = Perpendicular(end - center) * (ahead ? 1.0 : -1.0);
    const Vector leave = Ahead(*guides[1], found.Alongs[1], end) * static_cast<double>(keeps[1]);
    if (!(Dot(run, leave) > 0.0))
    {
      return std::nullopt;
    }
    return ahead;
  }

  /**
   *  Whether one choice is better than another: nearer what its touches should be near, then of the
   *  shorter arc.
   *
   *  @param  one   A choice.
   *  @param  other The best so far, if any.
   *  @return       Whether one is better.
   */
  auto Better(const Choice &one, const std::optional<Choice> &other) -> bool
  {
    if (!other)
    {
      return true;
    }
    if (one.Score < other->Score)
    {
      return true;
    }
    return Exactly(one.Score, other->Score) && std::abs(one.Sweep) < std::abs(other->Sweep);
  }

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
      -> std::optional<Choice>
  {
    return Flowing(guides, found, keeps)
        .transform(
            [&found, &keeps](bool ahead) -> Choice
            {
              return {.Found = found,
                      .Keeps = keeps,
                      .Ahead = ahead,
                      .Score = 0.0,
                      .Sweep = Swept(found.Touch.Center, found.Touch.Touches, ahead)};
            });
  }

  /**
   *  The circle in the angle between the parts of two curves kept from their corner, nearest it.
   *
   *  @param  guides  The curves' guides.
   *  @param  nook    The corner.
   *  @param  circles The circles of the radius touching both.
   *  @return         The circle whose touches lie on the parts kept and whose arc runs smoothly
   *                  from one onto the other, nearest the corner; none where none does.
   */
  auto Nooked(const Guides &guides, const Nook &nook, const std::vector<Found> &circles)
      -> std::optional<Choice>
  {
    const Found              &corner = *nook.At;
    const std::array<Side, 2> sides  = {Sided(*guides[0], corner.Alongs[0]),
                                        Sided(*guides[1], corner.Alongs[1])};
    const Point               spot   = corner.Touch.Touches[0];
    std::array<double, 2>     limits = {Infinite, Infinite};
    for (std::size_t index = 0; index < 2; ++index)
    {
      const Guide &guide = *guides.at(index);
      const Side  &side  = sides.at(index);
      for (const Found &other : nook.All)
      {
        double ahead = (other.Alongs.at(index) - side.Corner) * side.Keep;
        if (guide.Period > 0.0)
        {
          ahead -= std::floor(ahead / guide.Period) * guide.Period;
        }
        const bool met = Within(*guides[0], other.Alongs[0]) && Within(*guides[1], other.Alongs[1]);
        if (met && ahead > 0x1p-40)
        {
          limits.at(index) = std::min(limits.at(index), ahead);
        }
      }
    }
    std::optional<Choice> best;
    for (const Found &circle : circles)
    {
      Found              placed = circle;
      std::array<int, 2> keeps  = {};
      bool               inside = true;
      for (std::size_t index = 0; index < 2; ++index)
      {
        const Side &side  = sides.at(index);
        double     &along = placed.Alongs.at(index);
        along             = Near(*guides.at(index), along, side.Corner);
        const double past = along - side.Corner;
        const int    keep = side.Keep != 0 ? side.Keep : Sign(past);
        keeps.at(index)   = keep;
        inside            = inside && keep * past >= 0.0 && keep * past <= limits.at(index);
      }
      if (!inside)
      {
        continue;
      }
      if (std::optional<Choice> choice = Chosen(guides, placed, keeps))
      {
        choice->Score =
            Distance(placed.Touch.Touches[0], spot) + Distance(placed.Touch.Touches[1], spot);
        if (Better(*choice, best))
        {
          best = choice;
        }
      }
    }
    return best;
  }

  /**
   *  The ways a curve may be kept from a touch: the part from it the pick lies on, either way round
   *  where its parameter repeats.
   *
   *  @param  guide The curve's guide.
   *  @param  along The touch's parameter.
   *  @return       Each way, as the touch's parameter the part runs from or to and the side kept.
   */
  auto Options(const Guide &guide, double along) -> std::vector<std::pair<double, int>>
  {
    std::vector<std::pair<double, int>> options;
    if (guide.Period > 0.0)
    {
      double below = Near(guide, along, guide.Pick);
      if (below > guide.Pick)
      {
        below -= guide.Period;
      }
      const double above = below + guide.Period;
      if (guide.Whole || below > 1.0 - guide.Period)
      {
        options.emplace_back(below, 1);
      }
      if (guide.Whole || above < guide.Period)
      {
        options.emplace_back(above, -1);
      }
      return options;
    }
    if (guide.Pick >= along)
    {
      options.emplace_back(along, 1);
    }
    if (guide.Pick <= along)
    {
      options.emplace_back(along, -1);
    }
    return options;
  }

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
                const std::vector<Found> &circles) -> std::optional<Choice>
  {
    std::optional<Choice> best;
    for (const Found &circle : circles)
    {
      for (const auto &[near, keepNear] : Options(*guides[0], circle.Alongs[0]))
      {
        for (const auto &[far, keepFar] : Options(*guides[1], circle.Alongs[1]))
        {
          Found placed  = circle;
          placed.Alongs = {near, far};
          if (std::optional<Choice> choice = Chosen(guides, placed, {keepNear, keepFar}))
          {
            choice->Score = Distance(circle.Touch.Touches[0], picks[0])
                          + Distance(circle.Touch.Touches[1], picks[1]);
            if (Better(*choice, best))
            {
              best = choice;
            }
          }
        }
      }
    }
    return best;
  }

  /**
   *  The corner of two curves nearest their picks.
   *
   *  @param  guides  The curves' guides.
   *  @param  corners The corners.
   *  @param  picks   The picks.
   *  @return         The nearest; none for no corner.
   */
  auto Closest(const Guides &guides, const std::vector<Found> &corners,
               const std::array<Point, 2> &picks) -> std::optional<Found>
  {
    std::optional<Found> best;
    int                  fewest = 0;
    double               least  = 0.0;
    for (const Found &corner : corners)
    {
      const Point  spot     = corner.Touch.Touches[0];
      const double distance = Distance(spot, picks[0]) + Distance(spot, picks[1]);
      const int    off      = static_cast<int>(!Within(*guides[0], corner.Alongs[0]))
                            + static_cast<int>(!Within(*guides[1], corner.Alongs[1]));
      if (!best || off < fewest || (off == fewest && distance < least))
      {
        best   = corner;
        fewest = off;
        least  = distance;
      }
    }
    return best;
  }

  auto Placed(const Guides &guides, const Found &found, double radius) -> Found
  {
    if (!found.Exact)
    {
      return found;
    }
    Found placed = found;
    placed.Touch = Settled(Posed({guides[0]->Holder, guides[1]->Holder}, radius), found.Ask);
    for (std::size_t index = 0; index < 2; ++index)
    {
      const Guide &guide = *guides.at(index);
      placed.Alongs.at(index) =
          Near(guide, Along(guide, placed.Touch.Touches.at(index)), found.Alongs.at(index));
    }
    return placed;
  }
} // namespace Geoxy::Corner
