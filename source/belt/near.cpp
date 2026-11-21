/***************************************************************************************************
 *  @file      source/belt/near.cpp
 *  @brief     How near a point lies to a shape, and how the shape winds about it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "belt/near.hpp"

#include "belt/piece.hpp"

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace Geoxy::Belt
{
  namespace
  {
    /**
     *  Whether an arc rises as it leaves its start: on the circle's right where it runs
     *  counterclockwise, on its left where clockwise. One that starts at the circle's top or bottom
     *  is past it.
     *
     *  @param  arc The arc.
     *  @return     Whether it does.
     */
    auto Rising(const Piece &arc) -> bool
    {
      const Vector start = arc.Start - arc.Center;
      const double aside = arc.Left ? start.X : -start.X;
      return aside > 0.0 || (!(aside < 0.0) && start.Y < 0.0);
    }

    /**
     *  How a stretch of an arc that only rises or only falls winds about a point: once where it
     *  passes the point's height to the point's right, counterclockwise as it rises.
     *
     *  @param  arc    The arc.
     *  @param  levels The heights the stretch starts and ends at.
     *  @param  point  The point.
     *  @return        1, -1 or 0.
     */
    auto Passing(const Piece &arc, const std::array<double, 2> &levels, Point point) -> int
    {
      if ((levels[0] <= point.Y) == (levels[1] <= point.Y))
      {
        return 0;
      }
      const bool   rises = levels[1] > levels[0];
      const double rise  = point.Y - arc.Center.Y;
      const double reach = std::sqrt(std::max(0.0, (arc.Radius - rise) * (arc.Radius + rise)));
      const double where = rises == arc.Left ? arc.Center.X + reach : arc.Center.X - reach;
      if (!(where > point.X))
      {
        return 0;
      }
      return rises ? 1 : -1;
    }
  } // namespace

  auto Chorded(const std::array<Point, 2> &ends, Point point) -> int
  {
    if ((ends[0].Y <= point.Y) == (ends[1].Y <= point.Y))
    {
      return 0;
    }
    const bool   rises = ends[1].Y > ends[0].Y;
    const double side  = Cross(ends[1] - ends[0], point - ends[0]);
    if (rises ? side > 0.0 : side < 0.0)
    {
      return rises ? 1 : -1;
    }
    return 0;
  }

  auto Rounded(const Piece &arc, Point point) -> int
  {
    if (arc.Shallow)
    {
      /**
       *  A nearly straight arc winds as its chord does, and once more about a point between the
       *  two: on the side of the chord the arc bows to, and within its circle.
       */
      const Vector chord = arc.End - arc.Start;
      const double side  = Cross(chord, point - arc.Start);
      const bool   bowed = arc.Left ? side < 0.0 : side > 0.0;
      const int    lens  = bowed && Outside(arc, point) < 0.0 ? 1 : 0;
      return Chorded({arc.Start, arc.End}, point) + (arc.Left ? lens : -lens);
    }
    const Vector start  = arc.Start - arc.Center;
    const Vector upward = {.X = 0.0, .Y = 1.0};
    bool         rises  = Rising(arc);
    const double span   = Quartered(start, arc.End - arc.Center, arc.Left);
    const double first  = Quartered(start, rises ? upward : -upward, arc.Left);
    int          count  = 0;
    double       level  = arc.Start.Y;
    for (const double turn : {first, first + 2.0})
    {
      if (turn >= span)
      {
        break;
      }
      const double ended  = arc.Center.Y + (rises ? arc.Radius : -arc.Radius);
      count              += Passing(arc, {level, ended}, point);
      level               = ended;
      rises               = !rises;
    }
    return count + Passing(arc, {level, arc.End.Y}, point);
  }
} // namespace Geoxy::Belt
