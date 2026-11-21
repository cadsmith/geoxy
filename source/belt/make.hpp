/***************************************************************************************************
 *  @file      source/belt/make.hpp
 *  @brief     The curves a band's edge can run along, made: each piece moved to its sides, and the
 *             carriers they lie on.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BELT_MAKE_HPP
#define GEOXY_BELT_MAKE_HPP

#include "belt/carrier.hpp"
#include "belt/piece.hpp"
#include "belt/work.hpp"
#include "scratch.hpp"
#include "twofold.hpp"

#include "geoxy/offset.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy::Belt
{
  /**
   *  A direction of unit length.
   *
   *  @param  way A way, not of no length.
   *  @return     The same way, as long as 1.
   */
  inline auto Aimed(Vector way) -> Vector
  {
    return way / Long(way);
  }

  /**
   *  The shape's pieces from the polyline's vertices: lines and arcs, one from each vertex to the
   *  next.
   *
   *  @param work The polyline; its pieces are set.
   */
  void Shaped(Work &work);

  /**
   *  A piece moved to one side by the distance.
   *
   *  @param  piece The piece.
   *  @param  by    How far, to its left when positive.
   *  @return       The piece moved; an arc moved through its center with a radius of 0 or less.
   */
  inline auto Sided(const Piece &piece, double by) -> Piece
  {
    Piece moved = piece;
    if (!piece.Round)
    {
      /**
       *  The unit normal first, then the distance: a level or an upright line's is 0 and 1 exactly,
       *  and the line moves by just the distance.
       */
      const Vector way  = piece.End - piece.Start;
      const double span = Long(way);
      const Vector side = Vector{.X = -way.Y / span, .Y = way.X / span} * by;
      moved.Start       = piece.Start + side;
      moved.End         = piece.End + side;
      return moved;
    }
    moved.Radius =
        (by > 0.0) == piece.Left ? piece.Radius - std::abs(by) : piece.Radius + std::abs(by);
    if (piece.Shallow)
    {
      /**
       *  A nearly straight arc's ends are moved square to the way it heads at each.
       */
      const std::array<Vector, 2> ways = Tangents(piece);
      const double                span = Long(ways[0]);
      moved.Start = piece.Start + (Vector{.X = -ways[0].Y / span, .Y = ways[0].X / span} * by);
      moved.End   = piece.End + (Vector{.X = -ways[1].Y / span, .Y = ways[1].X / span} * by);
      return moved;
    }
    if (moved.Radius > 0.0)
    {
      /**
       *  Each end moved along its own radius, from where it is: a center far off is known less
       *  finely than the ends, and what is measured from it would be as coarse.
       */
      /**
       *  By the distance itself over the radius: the new radius less the old is the distance
       *  rounded at the radius's size, a hundred steps of a double off for an arc of thousands.
       */
      const double step = ((by > 0.0) == piece.Left ? -std::abs(by) : std::abs(by)) / piece.Radius;
      moved.Start       = piece.Start + ((piece.Start - piece.Center) * step);
      moved.End         = piece.End + ((piece.End - piece.Center) * step);
    }
    return moved;
  }

  /**
   *  A curve put among the curves.
   *
   *  @param  work  The curves.
   *  @param  curve The curve.
   *  @param  way   How far to its piece's left it was moved, 0 for a circle about a vertex.
   *  @return       Its number.
   */
  inline auto Put(Work &work, const Piece &curve, double way) -> std::size_t
  {
    work.Curves.push_back(curve);
    work.Ways.push_back(way);
    return work.Curves.size() - 1;
  }

  /**
   *  Every piece moved to its right and to its left, but an arc moved through its center: the one
   *  moved right as it runs, the one moved left run the other way, each with the band on its left.
   *
   *  @param work The shape; the curves are added.
   */
  void Moved(Work &work);

  /**
   *  How an end of a curve is made, told as the curve is drawn.
   *
   *  @param work The curves; the making is kept.
   *  @param end  The curve, and whether its end, not its start.
   *  @param by   The making.
   */
  inline void Tell(Work &work, const std::pair<std::size_t, bool> &end, const Recipe &by)
  {
    if (work.Told.size() <= end.first)
    {
      work.Told.resize(end.first + 1);
    }
    work.Told[end.first].at(end.second ? 1 : 0) = by;
  }

  /**
   *  A circle about a vertex or an end: where it is, the pieces beside it and the curves it joins.
   */
  struct Hinge
  {
    Point                 About = {};                 /**< The vertex. */
    std::array<Vector, 2> Sides = {};                 /**< The unit directions it runs from
                                                           and to, counterclockwise.          */
    std::array<std::size_t, 2> Beside = {None, None}; /**< The pieces beside it. */
    std::array<std::size_t, 2> Joined = {None, None}; /**< The curves it joins.  */
    std::array<std::size_t, 2> Ended  = {None, None}; /**< The piece each of its two ends is
                                                           moved from.                       */
    std::array<std::size_t, 2> Tails = {};            /**< For each: 1 where the vertex is
                                                           that piece's end, 0 its start.     */
    std::int8_t Flank = 0; /**< Which side of the shape it is on: 1 its right, -1 its left, 0
                                neither, as round an end.                                    */
  };

  /**
   *  An arc of a circle of the distance about a vertex, put among the curves.
   *
   *  A line beside it lies the distance from all of it by its making, and is not measured against
   *  it; an arc beside it may curl nearer, and is.
   *
   *  @param work  The shape; the curve is added.
   *  @param hinge The circle.
   */
  void Hinged(Work &work, const Hinge &hinge);

  /**
   *  The half circles beyond an open shape's two ends, each from one side's curve to the other's.
   *
   *  @param work The shape; the curves are added.
   */
  void Capped(Work &work);

  /**
   *  The box of the shape, grown by the pad, and the sizes that go by it.
   *
   *  @param work The shape; its box, scale, reach and pad are set.
   */
  void Sized(Work &work);

  /**
   *  The carrier of a curve that is no line across a corner or an end, which keep the ones they
   *  were drawn by: a circle about a vertex of the distance's radius, or its piece's carrier moved
   *  by the distance.
   *
   *  @param  work  The curves and the pieces' carriers.
   *  @param  curve The curve.
   *  @return       Its carrier.
   */
  inline auto Carrying(const Work &work, std::size_t curve) -> Carrier
  {
    const Piece &piece = work.Curves[curve];
    if (piece.Source != None)
    {
      return Shifted(work.WideShape[piece.Source], work.Shape[piece.Source].Left, work.Ways[curve]);
    }
    return {.Round  = true,
            .X      = Twofold::Lifted(piece.Center.X),
            .Y      = Twofold::Lifted(piece.Center.Y),
            .Radius = Twofold::Lifted(work.Asked.Distance)};
  }

  /**
   *  A piece's carrier, made when first asked.
   *
   *  @param  work  The shape.
   *  @param  piece The piece.
   *  @return       Its carrier.
   */
  inline auto Broad(Work &work, std::size_t piece) -> const Carrier &
  {
    if (work.Carried.empty())
    {
      work.WideShape.resize(work.Shape.size());
      work.Carried.assign(work.Shape.size(), 0);
    }
    if (work.Carried[piece] == 0)
    {
      work.Carried[piece]   = 1;
      work.WideShape[piece] = Carried(work.Shape[piece], work.Vertices[piece].Bulge);
    }
    return work.WideShape[piece];
  }

  /**
   *  A curve's carrier, made when first asked; its piece's carrier with it.
   *
   *  @param  work  The curves.
   *  @param  curve The curve.
   *  @return       Its carrier.
   */
  inline auto Widened(Work &work, std::size_t curve) -> const Carrier &
  {
    if (work.Curves[curve].Source != None)
    {
      Broad(work, work.Curves[curve].Source);
    }
    if (work.Made.size() < work.Curves.size())
    {
      work.WideCurves.resize(work.Curves.size());
      work.Made.resize(work.Curves.size(), 0);
    }
    if (work.Made[curve] == 0)
    {
      work.Made[curve]       = 1;
      work.WideCurves[curve] = Carrying(work, curve);
    }
    return work.WideCurves[curve];
  }

  /**
   *  A curve given the carrier it was drawn by, where lifting it from its rounded ends would make
   *  another, a step of a double off.
   *
   *  @param work  The curves; the curve's carrier is set.
   *  @param curve The curve.
   *  @param wide  Its carrier.
   */
  inline void Borne(Work &work, std::size_t curve, const Carrier &wide)
  {
    work.WideCurves.resize(work.Curves.size());
    work.Made.resize(work.Curves.size(), 0);
    work.Made[curve]       = 1;
    work.WideCurves[curve] = wide;
  }

  /**
   *  Whether a point of a curve's line or circle falls within the curve, or within the reach of an
   *  end of it.
   *
   *  @param  work  The curves.
   *  @param  curve The curve.
   *  @param  point The point.
   *  @return       Whether it does.
   */
  inline auto Within(const Work &work, std::size_t curve, Point point) -> bool
  {
    const Piece &piece = work.Curves[curve];
    const double along = Along(piece, point);
    const double hair  = work.Reach / std::max(Long(piece), work.Reach);
    return along >= -hair && along <= 1.0 + hair;
  }

  /**
   *  Where a point of a curve's line or circle lies along the curve, and whether it falls within it
   *  or within the reach of an end of it. A point of an arc's circle the arc does not sweep to, and
   *  not by an end of it, is told so by the ways to its ends, with no angle worked out.
   *
   *  @param  work  The curves.
   *  @param  curve The curve.
   *  @param  point The point.
   *  @return       The place, and whether it falls within; the place is not worked out where it
   *                does not.
   */
  inline auto Lying(const Work &work, std::size_t curve, Point point) -> std::pair<double, bool>
  {
    const Piece &piece = work.Curves[curve];
    if (piece.Round
        && !Facing(piece, point)
        && Dot(point - piece.Start, point - piece.Start) > 4.0 * work.Reach * work.Reach
        && Dot(point - piece.End, point - piece.End) > 4.0 * work.Reach * work.Reach)
    {
      return {0.0, false};
    }
    const double along = Along(piece, point);
    /**
     *  Within it outright, as most are: how long a hair is along it is not worked out.
     */
    if (along >= 0.0 && along <= 1.0)
    {
      return {along, true};
    }
    const double hair = work.Reach / std::max(Long(piece), work.Reach);
    return {along, along >= -hair && along <= 1.0 + hair};
  }

  /**
   *  Whether two curves run along one line or one circle.
   *
   *  @param  work   The reach.
   *  @param  first  One curve.
   *  @param  second The other, a line if the first is and an arc if it is.
   *  @return        Whether they do.
   */
  inline auto Together(const Work &work, const Piece &first, const Piece &second) -> bool
  {
    const double reach = work.Reach;
    if (first.Round)
    {
      return std::abs(first.Center.X - second.Center.X) <= reach
          && std::abs(first.Center.Y - second.Center.Y) <= reach
          && std::abs(first.Radius - second.Radius)     <= reach;
    }
    const Vector way   = first.End - first.Start;
    const double slack = reach * Long(way);
    return std::abs(Cross(way, second.Start - first.Start)) <= slack
        && std::abs(Cross(way, second.End - first.Start))   <= slack;
  }

  /**
   *  Whether two carriers are one line or one circle, as near as two doubles tell: two lines that
   *  run one way, the one through the other; two circles of one center and one radius.
   *
   *  @param  work   The shape's size.
   *  @param  first  One carrier.
   *  @param  second The other, a line if the first is and a circle if it is.
   *  @return        Whether they are.
   */
  inline auto Alike(const Work &work, const Carrier &first, const Carrier &second) -> bool
  {
    const Twofold::Double offX  = Twofold::Difference(second.X, first.X);
    const Twofold::Double offY  = Twofold::Difference(second.Y, first.Y);
    const double          slack = 1e-26 * work.Scale;
    if (first.Round)
    {
      const double far =
          std::max({std::abs(offX.High), std::abs(offY.High),
                    std::abs(Twofold::Difference(second.Radius, first.Radius).High)});
      return far <= slack;
    }
    const Twofold::Double turn  = Twofold::Difference(Twofold::Multiply(first.WayX, second.WayY),
                                                      Twofold::Multiply(first.WayY, second.WayX));
    const Twofold::Double aside = Twofold::Difference(Twofold::Multiply(first.WayX, offY),
                                                      Twofold::Multiply(first.WayY, offX));
    const double          along = std::abs(first.WayX.High) + std::abs(first.WayY.High);
    const double          next  = std::abs(second.WayX.High) + std::abs(second.WayY.High);
    return std::abs(turn.High) <= 1e-26 * along * next && std::abs(aside.High) <= slack * along;
  }
} // namespace Geoxy::Belt

#endif /* GEOXY_BELT_MAKE_HPP */
