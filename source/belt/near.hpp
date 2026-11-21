/***************************************************************************************************
 *  @file      source/belt/near.hpp
 *  @brief     How near a point lies to a shape, and how the shape winds about it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BELT_NEAR_HPP
#define GEOXY_BELT_NEAR_HPP

#include "belt/piece.hpp"
#include "belt/work.hpp"
#include "scratch.hpp"

#include "geoxy/box.hpp"
#include "geoxy/offset.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>

namespace Geoxy::Belt
{
  /**
   *  Whether a point lies outside a box.
   *
   *  @param  box   The box.
   *  @param  point The point.
   *  @return       Whether it does.
   */
  inline auto Beyond(const Box &box, Point point) -> bool
  {
    return point.X < box.Min.X || point.X > box.Max.X || point.Y < box.Min.Y || point.Y > box.Max.Y;
  }

  /**
   *  The shape's pieces a point may lie within the distance of.
   *
   *  @param  work  The shape.
   *  @param  point The point.
   *  @return       Their numbers.
   */
  inline auto Around(const Work &work, Point point) -> std::span<const std::size_t>
  {
    return work.Shape.size() <= Handful ? std::span<const std::size_t>{work.Everyone}
                                        : work.Near.Around(point);
  }

  /**
   *  Whether a piece of the shape is left out when a curve's point is measured: one the curve lies
   *  the distance from by its making; or an arc that meets the vertex the curve hinges on, where
   *  its nearest point is that vertex, which the curve lies the distance from likewise.
   *
   *  @param  work  The shape.
   *  @param  curve The curve, and the point of it measured.
   *  @param  piece The piece.
   *  @return       Whether it is.
   */
  inline auto Spared(const Work &work, const std::pair<const Piece *, Point> &curve,
                     std::size_t piece) -> bool
  {
    const Piece &from  = *curve.first;
    const Point  point = curve.second;
    if (piece == from.Spared[0] || piece == from.Spared[1])
    {
      return true;
    }
    if (piece != from.Beside[0] && piece != from.Beside[1])
    {
      return false;
    }
    const Piece &arc   = work.Shape[piece];
    const bool   first = Long(from.Pivot - arc.Start) < Long(from.Pivot - arc.End);
    const Point  other = first ? arc.End : arc.Start;
    return !Facing(arc, point) && Long(point - from.Pivot) <= Long(point - other);
  }

  /**
   *  Whether a point lies past a flat end of a piece of the shape: on or beyond the line square to
   *  an open shape's first piece at its start, or to its last at its end, where the ends are not
   *  round. What the piece sweeps ends at that line, so the point is not near it, however near the
   *  end itself.
   *
   *  @param  work  The shape.
   *  @param  piece The piece.
   *  @param  point The point.
   *  @return       Whether it does.
   */
  inline auto Past(const Work &work, std::size_t piece, Point point) -> bool
  {
    const bool first = piece == 0;
    const bool last  = piece + 1 == work.Shape.size();
    if (work.Closed || work.Asked.Ends == Cap::Round || (!first && !last))
    {
      return false;
    }
    const Piece &shape = work.Shape[piece];
    if (shape.Round)
    {
      /**
       *  The arc's own center is as far from every point of it, and past neither end.
       */
      if (point == shape.Center)
      {
        return false;
      }
      /**
       *  Within the arc's sweep by more than a rounding of the angle, the point is not past either
       *  end; else it is past the end it is nearer.
       */
      const double along  = Along(shape, point);
      const bool   starts = Long(point - shape.Start) < Long(point - shape.End);
      return (along <= 1e-12 || along >= 1.0 - 1e-12) && (starts ? first : last);
    }
    /**
     *  A point on the line across the end, as near as doubles tell, is past it: that line is the
     *  band's edge there, not within it.
     */
    const Vector way   = shape.End - shape.Start;
    const double slack = work.Reach * Long(way);
    return (first && Dot(point - shape.Start, way) <= slack)
        || (last && Dot(point - shape.End, way) >= -slack);
  }

  /**
   *  How far a curve's point past a flat end of a piece is from what the band still holds of that
   *  piece there: the vertex at its other end, where the shape goes on from it.
   *
   *  @param  work  The shape.
   *  @param  curve The curve, or none, and the point of it measured.
   *  @param  piece The piece, the shape's first or its last.
   *  @return       The distance to that vertex; infinite for a shape of one piece, flat at both
   *                ends, and for a curve that is a circle about that very vertex, the distance from
   *                it by its making.
   */
  inline auto Hub(const Work &work, const std::pair<const Piece *, Point> &curve, std::size_t piece)
      -> double
  {
    const Piece &shape = work.Shape[piece];
    const Point  hub   = piece == 0 ? shape.End : shape.Start;
    const Piece *from  = curve.first;
    const bool   about =
        from->Round && from->Source == None && Long(from->Center - hub) <= work.Reach;
    if (work.Shape.size() == 1 || about)
    {
      return std::numeric_limits<double>::infinity();
    }
    return Long(curve.second - hub);
  }

  /**
   *  How near the shape comes to a point of a curve, the pieces the curve is spared left out.
   *
   *  @param  work  The shape.
   *  @param  from  The curve.
   *  @param  point The point.
   *  @return       The distance to the nearest of the others, or to one of them nearer than the
   *                distance by more than a billionth of the shape's size; infinite where there is
   *                none.
   */
  inline auto Nearest(const Work &work, const Piece &from, Point point) -> double
  {
    const Piece *curve = &from;
    const bool   flat  = !work.Closed && work.Asked.Ends != Cap::Round;
    /**
     *  Nearer than the distance by more than a billionth of the shape's size is inside the band,
     *  however much nearer: no more is measured.
     */
    const double inside = work.Asked.Distance - (work.Scale * 1e-9);
    double       best   = std::numeric_limits<double>::infinity();
    const auto   tried  = [&](std::size_t piece) -> bool
    {
      const bool beside = piece == from.Beside[0] || piece == from.Beside[1];
      if (Beyond(work.Pads[piece], point)
          || piece == from.Spared[0]
          || piece == from.Spared[1]
          || (beside && Spared(work, {curve, point}, piece)))
      {
        return false;
      }
      best = std::min(best, flat && Past(work, piece, point) ? Hub(work, {curve, point}, piece)
                                                             : Apart(work.Shape[piece], point));
      return best < inside;
    };
    /**
     *  The pieces before and after the curve's own are the likeliest to hold it within the
     *  distance, where it is: they are measured first.
     */
    const std::size_t count = work.Shape.size();
    if (from.Source != None
        && count > 2
        && (tried((from.Source + 1) % count) || tried((from.Source + count - 1) % count)))
    {
      return best;
    }
    for (const std::size_t piece : Around(work, point))
    {
      if (tried(piece))
      {
        break;
      }
    }
    return best;
  }

  /**
   *  How a line winds about a point, by a ray from the point to its right: once where its two ends
   *  lie either side of the ray's level, an end on the level counting as below, and it passes to
   *  the point's right; counterclockwise as it rises.
   *
   *  @param  ends  The line's start and end.
   *  @param  point The point.
   *  @return       1, -1 or 0.
   */
  auto Chorded(const std::array<Point, 2> &ends, Point point) -> int;

  /**
   *  How an arc winds about a point, by a ray from the point to its right: the arc rises or falls
   *  between its ends and the circle's top and bottom it passes, each at the most once, and each
   *  such stretch is counted where it passes the point's height to the point's right. No angle is
   *  worked out.
   *
   *  @param  arc   The arc.
   *  @param  point The point.
   *  @return       The crossings counted, counterclockwise above 0.
   */
  auto Rounded(const Piece &arc, Point point) -> int;

  /**
   *  How many times a loop of pieces winds about a point: the ray from it toward +x, each crossing
   *  signed.
   *
   *  @param  loop  The pieces, each ending where the next starts.
   *  @param  point The point.
   *  @return       The winding number.
   */
  inline auto Wound(std::span<const Piece> loop, Point point) -> int
  {
    int count = 0;
    for (const Piece &piece : loop)
    {
      if (piece.Round)
      {
        count += Rounded(piece, point);
        continue;
      }
      count += Chorded({piece.Start, piece.End}, point);
    }
    return count;
  }
} // namespace Geoxy::Belt

#endif /* GEOXY_BELT_NEAR_HPP */
