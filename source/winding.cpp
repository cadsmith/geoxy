/***************************************************************************************************
 *  @file      source/winding.cpp
 *  @brief     How many times a curve run between two points crosses the ray from a point: lines by
 *             the exact orientation, arcs by quarters and the exact distance from their center,
 *             other curves halved until no part can wind about the point.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "winding.hpp"

#include "expansion.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/length.hpp"
#include "geoxy/line.hpp"
#include "geoxy/orient.hpp"
#include "geoxy/point.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  The most turns of an arc counted.
     */
    constexpr double MostArcTurns = 0x1p24;

    /**
     *  How many times a traverse is halved at most, past which its part is taken as its line.
     */
    constexpr int MostHalvings = 64;

    /**
     *  How many times the line from one point to another crosses the ray from a point, by Sunday's
     *  rule: upward where it starts at or below the point's height and ends above it with the point
     *  on its left, downward the other way round with the point on its right.
     *
     *  @param  start Where it starts.
     *  @param  end   Where it ends.
     *  @param  point The point.
     *  @return       1, -1 or 0.
     */
    auto Sunday(Point start, Point end, Point point) -> std::int64_t
    {
      if (start.Y <= point.Y)
      {
        return point.Y < end.Y && Orient(start, end, point) == Orientation::Counterclockwise ? 1
                                                                                             : 0;
      }
      return end.Y <= point.Y && Orient(start, end, point) == Orientation::Clockwise ? -1 : 0;
    }

    /**
     *  Where a point lies against an arc's circle, exactly.
     *
     *  @param  arc   The arc.
     *  @param  point The point.
     *  @return       The sign of the square of its distance from the center less the radius's: -1
     *                inside, 0 on it, 1 outside.
     */
    auto Inward(const Arc &arc, Point point) -> int
    {
      const double across = std::abs(point.X - arc.Center.X) + std::abs(point.Y - arc.Center.Y);
      const double radius = std::abs(arc.Radius);
      /**
       *  Each difference is within a step of a double of the exact one: a sum of the two under half
       *  the radius puts the point inside, its distance at most that sum; one over twice the radius
       *  outside, its distance at least that sum over the square root of 2.
       */
      if (across < radius / 2)
      {
        return -1;
      }
      if (across > 2 * radius)
      {
        return 1;
      }
      using namespace Expansion;
      const Number dx = Difference(Of(point.X), Of(arc.Center.X));
      const Number dy = Difference(Of(point.Y), Of(arc.Center.Y));
      const Number r  = Of(arc.Radius);
      return Sign(Difference(Sum(Product(dx, dx), Product(dy, dy)), Product(r, r)));
    }

    /**
     *  A part of an arc that rises or falls all along, within one half of its circle.
     */
    struct Climb
    {
      Point From  = {};    /**< Where it starts.                                     */
      Point To    = {};    /**< Where it ends.                                       */
      bool  Right = false; /**< Whether it lies in the half right of the center's x. */
    };

    /**
     *  How many times a part of an arc that rises or falls all along crosses the ray from a point.
     *
     *  Where it spans the point's height, the circle's point at that height in its half is right of
     *  the point where the point is left of the center, or, in the right half, inside the circle;
     *  in the left half, where the point is left of the center and outside the circle.
     *
     *  @param  arc   The arc.
     *  @param  rim   The part.
     *  @param  point The point.
     *  @return       1, -1 or 0.
     */
    auto Climbed(const Arc &arc, const Climb &rim, Point point) -> std::int64_t
    {
      const bool rises = rim.From.Y <= point.Y && point.Y < rim.To.Y;
      const bool falls = rim.To.Y <= point.Y && point.Y < rim.From.Y;
      if (!rises && !falls)
      {
        return 0;
      }
      const bool before = point.X < arc.Center.X;
      const bool beyond =
          rim.Right ? before || Inward(arc, point) < 0 : before && Inward(arc, point) > 0;
      if (!beyond)
      {
        return 0;
      }
      return rises ? 1 : -1;
    }

    /**
     *  The point of an arc's circle at a whole number of quarter turns, as At works it out: the
     *  center and the radius along an axis.
     *
     *  @param  arc     The arc.
     *  @param  quarter The quarter turns, from the x axis counterclockwise.
     *  @return         The point.
     */
    auto Axial(const Arc &arc, double quarter) -> Point
    {
      const double turned = quarter - (4.0 * std::floor(quarter / 4.0));
      const Point  center = arc.Center;
      if (turned < 1.0)
      {
        return {.X = center.X + arc.Radius, .Y = center.Y};
      }
      if (turned < 2.0)
      {
        return {.X = center.X, .Y = center.Y + arc.Radius};
      }
      if (turned < 3.0)
      {
        return {.X = center.X - arc.Radius, .Y = center.Y};
      }
      return {.X = center.X, .Y = center.Y - arc.Radius};
    }

    /**
     *  Whether the quarter of a circle between two quarter turns lies right of its center's x.
     *
     *  @param  arc     The arc, whose radius may be negative and turn its points half about.
     *  @param  quarter The quarter turns its part starts after.
     *  @return         Whether it does.
     */
    auto Easterly(const Arc &arc, double quarter) -> bool
    {
      const double turned = quarter - (4.0 * std::floor(quarter / 4.0));
      const bool   right  = turned < 1.0 || turned >= 3.0;
      return right != (arc.Radius < 0.0);
    }

    /**
     *  How many times an arc of less than a turn, from one point to another, crosses the ray from a
     *  point: cut at each quarter turn it passes, from its start, into parts that rise or fall all
     *  along.
     *
     *  @param  arc   The arc, its start in degrees reduced below a turn and its sweep below one.
     *  @param  ends  Where it starts and ends.
     *  @param  point The point.
     *  @return       The crossings, signed.
     */
    auto Quadrants(const Arc &arc, const Traverse &ends, Point point) -> std::int64_t
    {
      const double start = arc.Start.InDegrees() / 90.0;
      const double sweep = arc.Sweep.InDegrees() / 90.0;
      const double step  = sweep > 0.0 ? 1.0 : -1.0;
      double       next  = sweep > 0.0 ? std::floor(start) + 1.0 : std::ceil(start) - 1.0;
      Climb        rim   = {.From = ends.From, .Right = Easterly(arc, next - std::max(step, 0.0))};
      std::int64_t count = 0;
      while (step * (start + sweep - next) > 0.0)
      {
        rim.To  = Axial(arc, next);
        count  += Climbed(arc, rim, point);
        rim     = {.From = rim.To, .Right = Easterly(arc, next + std::min(step, 0.0))};
        next   += step;
      }
      rim.To = ends.To;
      return count + Climbed(arc, rim, point);
    }

    /**
     *  How many times an arc crosses the ray from a point: once for each whole turn where the point
     *  is inside its circle, and the rest of its sweep by quarters.
     *
     *  @param  arc   The arc.
     *  @param  ends  Where it starts and ends.
     *  @param  point The point.
     *  @return       The crossings, signed.
     */
    auto Circling(const Arc &arc, const Traverse &ends, Point point) -> std::int64_t
    {
      const double sweep = arc.Sweep.InDegrees();
      if (!(std::abs(sweep) > 0.0))
      {
        return Sunday(ends.From, ends.To, point);
      }
      /**
       *  The remainder is exact, and so is what it leaves: a whole number of turns.
       */
      const double left = std::fmod(sweep, 360.0);
      const double turns =
          std::fmax(-MostArcTurns, std::fmin(MostArcTurns, (sweep - left) / 360.0));
      const auto         whole  = static_cast<std::int64_t>(turns);
      const Arc          rest   = {.Center = arc.Center,
                                   .Radius = arc.Radius,
                                   .Start  = Degrees(std::fmod(arc.Start.InDegrees(), 360.0)),
                                   .Sweep  = Degrees(left)};
      const std::int64_t inside = whole != 0 && Inward(arc, point) < 0 ? whole : 0;
      return inside + Quadrants(rest, ends, point);
    }

    /**
     *  A stretch of a curve's parameter and the points that stand for its ends.
     */
    struct Stint
    {
      Range Along = {}; /**< The stretch of the parameter. */
      Point From  = {}; /**< The point at its start.       */
      Point To    = {}; /**< The point at its end.         */
    };

    /**
     *  A curve's crossings with the ray from a point, part by part.
     */
    class Halver
    {
    public:
      /**
       *  The crossings of a curve with the ray from a point.
       *
       *  @param curve The curve.
       *  @param point The point.
       */
      Halver(const Curve &curve, Point point) : m_Curve(&curve), m_Point(point) {}

      /**
       *  How many times a stretch crosses the ray: as the line between its ends where it is shorter
       *  than half the distance from the point to its start, as it then lies in a disc about its
       *  start that does not hold the point and sweeps about it as that line does, or where its
       *  length is not a number; else each half.
       *
       *  @param  stretch The stretch.
       *  @param  depth   How many more times it may be halved.
       *  @return         The crossings, signed.
       */
      [[nodiscard]] auto Count(const Stint &stretch, int depth) const -> std::int64_t
      {
        const double reach  = Distance(stretch.From, m_Point);
        const double length = Length(*m_Curve, stretch.Along, Accuracy{.Error = reach / 8.0});
        if (depth == 0 || !(length >= reach / 2.0))
        {
          return Sunday(stretch.From, stretch.To, m_Point);
        }
        const double middle = (stretch.Along.Start / 2.0) + (stretch.Along.End / 2.0);
        const Point  split  = At(*m_Curve, middle);
        return Count({.Along = {.Start = stretch.Along.Start, .End = middle},
                      .From  = stretch.From,
                      .To    = split},
                     depth - 1)
             + Count({.Along = {.Start = middle, .End = stretch.Along.End},
                      .From  = split,
                      .To    = stretch.To},
                     depth - 1);
      }

    private:
      const Curve *m_Curve; /**< The curve. */
      Point        m_Point; /**< The point. */
    };
  } // namespace

  auto Ascents(const Traverse &traverse, Point point) -> std::int64_t
  {
    const Curve &curve = *traverse.Shape;
    if (std::holds_alternative<Line>(curve))
    {
      return Sunday(traverse.From, traverse.To, point);
    }
    if (const auto *arc = std::get_if<Arc>(&curve))
    {
      return Circling(*arc, traverse, point);
    }
    const Halver swept{curve, point};
    return swept.Count({.Along = Domain(curve), .From = traverse.From, .To = traverse.To},
                       MostHalvings);
  }
} // namespace Geoxy
