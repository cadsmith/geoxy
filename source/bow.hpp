/***************************************************************************************************
 *  @file      source/bow.hpp
 *  @brief     A bulged segment worked out from its start along its chord, not from its center.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BOW_HPP
#define GEOXY_BOW_HPP

#include "trig.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/point.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <cmath>

namespace Geoxy
{
  /**
   *  The bulge at and over which no arc is remote: its circle's center and radius are then within
   *  4096 times the size of its ends, wherever they lie.
   */
  inline constexpr double Slight = 0x1p-13;

  /**
   *  A bulged segment as its start, its chord and its sweep.
   *
   *  A small bulge puts the center far off: at 5e6 a chord of 1000 bulged 1e-6 has its center 2.5e8
   *  away, whose rounding is many times a pick's width. Worked from the start along the chord,
   *  every point keeps the precision of the vertices.
   */
  struct Bight
  {
    Point  Start = {};  /**< Where it starts.                                     */
    Vector Chord = {};  /**< From its start to its end.                           */
    double Sweep = 0.0; /**< Its sweep in radians, four atan(bulge).              */
    double Half  = 0.0; /**< The sine of half its sweep, 2 / (bulge + 1 / bulge). */
  };

  /**
   *  A bulged segment's circle.
   */
  struct Hoop
  {
    Point  Center = {};  /**< Its center. */
    double Radius = 0.0; /**< Its radius. */
  };

  /**
   *  The circle of the bow from one vertex to the next by a bulge, as every arc made from a bulge
   *  has it: the center off the chord's midpoint by (1 - b^2) / (4 b) of the chord turned a
   *  quarter, the radius (1 + b^2) / (4 |b|) of the chord's length.
   *
   *  @param  start Where it starts.
   *  @param  end   Where it ends, not its start.
   *  @param  bulge The bulge, not under Straight in magnitude.
   *  @return       The circle.
   */
  inline auto Hooped(Point start, Point end, double bulge) -> Hoop
  {
    const Vector chord  = end - start;
    const double square = bulge * bulge;
    return {.Center =
                Midpoint(start, end) + (Perpendicular(chord) * ((1.0 - square) / (4.0 * bulge))),
            .Radius = Length(chord) * ((1.0 + square) / (4.0 * std::abs(bulge)))};
  }

  /**
   *  Whether a bulged segment is remote: its circle, the center's larger coordinate grown by the
   *  radius, more than 4096 times the size of its ends' numbers. What is worked out through such a
   *  center is rounded by thousands of steps of a double at the ends, where the arc is: a remote
   *  arc is worked out from its ends and its bulge. An arc that is small beside how far it lies
   *  from the origin is not remote however slight its bulge, its center being told as well as its
   *  ends are.
   *
   *  @param  start Where it starts.
   *  @param  end   Where it ends.
   *  @param  bulge Its bulge.
   *  @return       Whether it is; never for a bulge under Straight, of Slight or over, or a segment
   *                of no length.
   */
  inline auto Remote(Point start, Point end, double bulge) -> bool
  {
    if (!(std::abs(bulge) >= Straight) || !(std::abs(bulge) < Slight) || start == end)
    {
      return false;
    }
    const Hoop   hoop = Hooped(start, end, bulge);
    const double size =
        std::max({std::abs(start.X), std::abs(start.Y), std::abs(end.X), std::abs(end.Y)});
    return !(std::max(std::abs(hoop.Center.X), std::abs(hoop.Center.Y)) + hoop.Radius
             <= 4096.0 * size);
  }

  /**
   *  A displacement turned by an angle.
   *
   *  @param  vector The displacement.
   *  @param  angle  The angle, counterclockwise.
   *  @return        The displacement turned, exactly at every quarter turn.
   */
  inline auto Turned(Vector vector, Angle angle) -> Vector
  {
    const Vector way = Direction(angle);
    return {.X = (way.X * vector.X) - (way.Y * vector.Y),
            .Y = (way.Y * vector.X) + (way.X * vector.Y)};
  }

  /**
   *  The bow from one vertex to the next by a bulge.
   *
   *  The sine of half its sweep, sin(2 atan b), is 2 b / (1 + b^2), worked out from the bulge with
   *  no trigonometry: near a whole turn the sweep's own rounding would move it, and the radius c /
   *  (2 sin(s / 2)) with it, by hundreds of steps of a double.
   *
   *  @param  start Where it starts.
   *  @param  end   Where it ends.
   *  @param  bulge The bulge, not under Straight in magnitude.
   *  @return       The bow.
   */
  inline auto Bowed(Point start, Point end, double bulge) -> Bight
  {
    return {.Start = start,
            .Chord = end - start,
            .Sweep = 4.0 * Trig::Atan(bulge),
            .Half  = 2.0 / (bulge + (1.0 / bulge))};
  }

  /**
   *  The point a part of a bow's sweep along it: its start, plus the chord to that point, which is
   *  the whole chord turned by half the sweep left over and scaled by the ratio of the half sines.
   *
   *  @param  bow   The bow.
   *  @param  swept How much of its sweep, in radians, its own sign.
   *  @return       The point.
   */
  inline auto Along(const Bight &bow, double swept) -> Point
  {
    const double scale = Trig::Sin(swept / 2.0) / bow.Half;
    return bow.Start + Turned(bow.Chord * scale, Radians((swept - bow.Sweep) / 2.0));
  }

  /**
   *  A bow's derivative with respect to the fraction of its sweep.
   *
   *  @param  bow   The bow.
   *  @param  swept How much of its sweep, in radians.
   *  @return       The tangent, its length the radius times the sweep.
   */
  inline auto Tangent(const Bight &bow, double swept) -> Vector
  {
    const double speed = bow.Sweep / (2.0 * bow.Half);
    return Turned(bow.Chord * speed, Radians(swept - (bow.Sweep / 2.0)));
  }
} // namespace Geoxy

#endif /* GEOXY_BOW_HPP */
