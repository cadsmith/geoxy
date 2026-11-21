/***************************************************************************************************
 *  @file      source/belt/corner.cpp
 *  @brief     The corners of a band's edge: the round about a vertex, its wedge and its mitre.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "belt/corner.hpp"

#include "belt/carrier.hpp"
#include "belt/cut.hpp"
#include "belt/make.hpp"
#include "belt/near.hpp"
#include "belt/piece.hpp"
#include "belt/work.hpp"
#include "scratch.hpp"
#include "twofold.hpp"

#include "geoxy/box.hpp"
#include "geoxy/offset.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy::Belt
{
  namespace
  {
    /**
     *  A corner the shape turns at, and its outer side.
     */
    struct Corner
    {
      std::size_t Before = 0;     /**< The piece that ends at it.                   */
      std::size_t After  = 0;     /**< The piece that starts at it.                 */
      bool        Right  = false; /**< Whether its outer side is the shape's right. */
      Vector      Into   = {};    /**< The heading into it, of unit length.         */
    };

    /**
     *  How far a piece runs on past its end, along its line or round its circle, to a point of it.
     *
     *  @param  piece The piece.
     *  @param  point The point, on its line or circle.
     *  @param  reach Within which the point is the piece's end.
     *  @return       The length of the run; below 0 where the piece does not come to the point: a
     *                line whose end lies past it, an arc that comes round to its own start first.
     */
    auto Run(const Piece &piece, Point point, double reach) -> double
    {
      if (!piece.Round)
      {
        const double past = (Along(piece, point) - 1.0) * Long(piece);
        return past >= -reach ? std::max(past, 0.0) : -1.0;
      }
      if (piece.Shallow)
      {
        /**
         *  Along a nearly straight arc a run is as long as the way straight to its end.
         */
        const Vector ahead = Tangents(piece)[1];
        const double past  = Dot(point - piece.End, ahead) / Long(ahead);
        return past >= -reach ? std::max(past, 0.0) : -1.0;
      }
      const Vector end   = piece.End - piece.Center;
      const double turn  = Swung(end, point - piece.Center, piece.Left);
      const double sweep = Swung(piece.Start - piece.Center, end, piece.Left);
      return turn < Whole - sweep ? turn * piece.Radius : -1.0;
    }

    /**
     *  A stretch of a piece's line or circle from one point of it to another, the way the piece
     *  runs.
     *
     *  @param  piece The piece.
     *  @param  from  Where the stretch starts.
     *  @param  onto  Where it ends.
     *  @return       The stretch.
     */
    auto Stretched(Piece piece, Point from, Point onto) -> Piece
    {
      piece.Start = from;
      piece.End   = onto;
      return piece.Round ? Arced(piece) : piece;
    }

    /**
     *  The two pieces at a corner moved to its outer side, as they run.
     *
     *  @param  work   The shape.
     *  @param  corner The corner.
     *  @return        The piece before it and the piece after it, moved.
     */
    auto Outer(const Work &work, const Corner &corner) -> std::array<Piece, 2>
    {
      const double by = corner.Right ? -work.Asked.Distance : work.Asked.Distance;
      return {Sided(work.Shape[corner.Before], by), Sided(work.Shape[corner.After], by)};
    }

    /**
     *  The point a corner's two moved pieces meet at when run on: of the points their lines or
     *  circles share that the first comes to past its end and the second before its start, the one
     *  nearest the corner.
     *
     *  @param  work   The shape and the curves.
     *  @param  corner The corner.
     *  @param  sides  The two moved pieces.
     *  @return        The point and whether there is one.
     */
    auto Meeting(Work &work, const Corner &corner, const std::array<Piece, 2> &sides)
        -> std::pair<Point, bool>
    {
      Hits hits = Met(sides[0], sides[1]);
      if (hits.Close)
      {
        const List<std::size_t> &moved = corner.Right ? work.Rights : work.Lefts;
        const Carrier            wide  = Widened(work, moved[corner.Before]);
        hits                           = Met(wide, Widened(work, moved[corner.After]));
      }
      const Point about = work.Shape[corner.After].Start;
      Point       best  = {};
      double      least = std::numeric_limits<double>::infinity();
      for (const Point point : std::span{hits.At}.first(hits.Count))
      {
        const double apart = Long(point - about);
        if (Run(sides[0], point, work.Reach)              >= 0.0
            && Run(Reversed(sides[1]), point, work.Reach) >= 0.0
            && apart                                      < least)
        {
          least = apart;
          best  = point;
        }
      }
      return {best, least < std::numeric_limits<double>::infinity()};
    }

    /**
     *  The line across a corner: as doubles draw it, and in two doubles a number.
     */
    struct Across
    {
      Piece   Drawn = {}; /**< The line, a stretch of it about its foot. */
      Carrier Wide  = {}; /**< Its carrier.                              */
    };

    /**
     *  The first point of a line that a piece comes to, run on past its end. Where the piece's
     *  circle comes near touching the line, the two are met in two doubles a number; and where they
     *  are then found to miss by less than their numbers' rounding, they touch.
     *
     *  @param  work   The reach and the curves.
     *  @param  moved  The piece, and the number of the curve along its line or circle.
     *  @param  across The line.
     *  @return        The point and whether there is one.
     */
    auto Reached(Work &work, const std::pair<Piece, std::size_t> &moved, const Across &across)
        -> std::pair<Point, bool>
    {
      const Piece &piece = moved.first;
      /**
       *  Always in two doubles: the line is two doubles a number, and the stretch of it doubles
       *  hold is that line rounded.
       */
      const Hits hits =
          Finer(Met(piece, across.Drawn), Met(Widened(work, moved.second), across.Wide));
      Point  best  = {};
      double least = std::numeric_limits<double>::infinity();
      for (const Point point : std::span{hits.At}.first(hits.Count))
      {
        const double run = Run(piece, point, work.Reach);
        if (run >= 0.0 && run < least)
        {
          least = run;
          best  = point;
        }
      }
      return {best, least < std::numeric_limits<double>::infinity()};
    }

    /**
     *  A wedge put among the wedges, and its outer edge among the curves.
     *
     *  @param work   The shape; the wedge and the curves are added.
     *  @param corner The corner, or the end, it is at.
     *  @param edge   Its outer edge from the first moved piece's end to the second's start, as the
     *                shape runs; each stretch along a moved piece's line or circle with that
     *                piece's number as its source.
     */
    void Wedged(Work &work, const Corner &corner, std::span<const Piece> edge)
    {
      const Point about = work.Shape[corner.After].Start;
      Wedge       wedge;
      wedge.Rim.push_back({.Start = about, .End = edge.front().Start});
      wedge.Rim.insert(wedge.Rim.end(), edge.begin(), edge.end());
      wedge.Rim.push_back({.Start = edge.back().End, .End = about});
      for (const Piece &piece : wedge.Rim)
      {
        wedge.Hull = Joined(wedge.Hull, Hulled(piece));
      }
      const Vector leaving  = edge.front().Start - about;
      const Vector entering = edge.back().End - about;
      wedge.Coned           = true;
      wedge.Left            = corner.Right;
      wedge.About           = about;
      wedge.Span = corner.Right ? std::array{leaving, entering} : std::array{entering, leaving};
      wedge.Edges.assign(wedge.Rim.size(), None);
      const double way = corner.Right ? -work.Asked.Distance : work.Asked.Distance;
      for (std::size_t index = 0; index < edge.size(); ++index)
      {
        /**
         *  With the band on its left: as the shape runs on its right, the other way and in the
         *  other order on its left.
         */
        Piece curve  = corner.Right ? edge[index] : Reversed(edge[edge.size() - 1 - index]);
        curve.Bounds = work.Wedges.size();
        curve.Flank  = corner.Right ? 1 : -1;
        /**
         *  A stretch along a moved piece's line or circle lies the distance from that piece, and
         *  the line across from the corner's two, where they are lines: an arc may curl nearer.
         */
        const auto lined = [&work](std::size_t beside, bool line) -> std::size_t
        {
          return work.Shape[beside].Round != line ? beside : None;
        };
        curve.Spared = {curve.Source, None};
        if (curve.Source == None)
        {
          curve.Spared = {lined(corner.Before, true), lined(corner.After, true)};
          curve.Beside = {lined(corner.Before, false), lined(corner.After, false)};
          curve.Pivot  = about;
        }
        curve.Joined = {None, None};
        wedge.Edges[1 + (corner.Right ? index : edge.size() - 1 - index)] =
            Put(work, curve, curve.Source == None ? 0.0 : way);
      }
      work.Wedges.push_back(std::move(wedge));
    }

    /**
     *  Whether a point lies to the left of a piece of the shape as it runs, in two doubles a
     *  number, and whether that is certain: beside the piece itself, not past an end of it, and not
     *  on it as near as two doubles tell.
     *
     *  @param  work  The polyline's vertices and the shape.
     *  @param  piece The piece.
     *  @param  point The point.
     *  @return       Whether it lies to the left, and whether that is certain.
     */
    auto Sided(const Work &work, std::size_t piece, Point point) -> std::pair<bool, bool>
    {
      const Piece          &shape = work.Shape[piece];
      const Carrier         wide  = Carried(shape, work.Vertices[piece].Bulge);
      const Twofold::Double offX  = Twofold::Difference(Twofold::Lifted(point.X), wide.X);
      const Twofold::Double offY  = Twofold::Difference(Twofold::Lifted(point.Y), wide.Y);
      if (shape.Round)
      {
        const Twofold::Double gap = Twofold::Difference(
            Twofold::Add(Twofold::Multiply(offX, offX), Twofold::Multiply(offY, offY)),
            Twofold::Multiply(wide.Radius, wide.Radius));
        const bool inside = gap.High < 0.0;
        return {inside == shape.Left,
                Facing(shape, point) && std::abs(gap.High) > 1e-26 * shape.Radius * shape.Radius};
      }
      const Twofold::Double cross = Twofold::Difference(Twofold::Multiply(wide.WayX, offY),
                                                        Twofold::Multiply(wide.WayY, offX));
      const double          along = Along(shape, point);
      const double size = (wide.WayX.High * wide.WayX.High) + (wide.WayY.High * wide.WayY.High);
      return {cross.High > 0.0, along > 0.0 && along < 1.0 && std::abs(cross.High) > 1e-26 * size};
    }

    /**
     *  How the winding about a corner itself, a hair out along a way, differs from the winding
     *  about a point farther out along it: by each other piece of the shape that passes between the
     *  two, the corner on one side of it and the point on the other.
     *
     *  A piece may pass a corner nearer than doubles can place a point beside it, and what the
     *  shape fills just outside the corner is what it fills on the corner's side of that piece.
     *
     *  @param  work   The shape.
     *  @param  corner The corner.
     *  @param  probe  The point farther out.
     *  @return        The winding about the corner's side less the winding about the point.
     */
    auto Slipped(const Work &work, const Corner &corner, Point probe) -> int
    {
      const Point  about = work.Shape[corner.After].Start;
      const double reach = 2.0 * Long(probe - about);
      int          slip  = 0;
      for (std::size_t piece = 0; piece < work.Shape.size(); ++piece)
      {
        if (piece                              == corner.Before
            || piece                           == corner.After
            || Apart(work.Shape[piece], about) > reach)
        {
          continue;
        }
        const std::pair<bool, bool> near = Sided(work, piece, about);
        if (near.second)
        {
          slip += static_cast<int>(near.first) - static_cast<int>(Sided(work, piece, probe).first);
        }
      }
      return slip;
    }

    /**
     *  Whether a corner's mitre is a part of the band asked: both ways, every one; outward, one
     *  that points out of what the shape fills, and inward into it, read a little out from the
     *  corner along the middle of its outer side, or straight ahead where it turns right back.
     *
     *  @param  work   The shape.
     *  @param  corner The corner.
     *  @param  sides  The two moved pieces.
     *  @return        Whether it is.
     */
    auto Wanted(const Work &work, const Corner &corner, const std::array<Piece, 2> &sides) -> bool
    {
      if (!work.Closed || work.Asked.Sides == Side::Both)
      {
        return true;
      }
      const Point  about    = work.Shape[corner.After].Start;
      const Vector middle   = (sides[0].End - about) + (sides[1].Start - about);
      const Vector out      = Long(middle) > work.Reach ? Aimed(middle) : corner.Into;
      const Point  probe    = about + (out * (work.Scale * 1e-9));
      const int    windings = Wound(work.Shape, probe) + Slipped(work, corner, probe);
      const bool   filled   = work.Asked.Fill == Fills::EvenOdd ? windings % 2 != 0 : windings != 0;
      return filled == (work.Asked.Sides == Side::Inward);
    }

    /**
     *  A corner's mitre as it is found: where its outer edge turns, and how far each of the two
     *  moved pieces is run to get there.
     */
    struct Mitre
    {
      std::size_t Count = 0;           /**< How many turns: none for no mitre, one where the
                                            pieces meet, two for the ends of the line across. */
      std::array<Point, 2>  At   = {}; /**< The turns, the first Count. */
      std::array<double, 2> Runs = {}; /**< How far the piece before is run on past its end,
                                            and the piece after back before its start.        */
      Carrier Wide = {};               /**< The line across, in two doubles a number. */
    };

    /**
     *  A corner's mitre: its two moved pieces run on to where they meet, or to the line across,
     *  tangent to the corner's circle halfway round it. None where a moved piece is an arc moved
     *  through its center, where the two never meet or do not reach the line, or where the corner
     *  points the other way from what is asked: the corner is then round.
     *
     *  @param  work   The shape.
     *  @param  corner The corner.
     *  @return        The mitre, of no turn where there is none.
     */
    auto Mitring(Work &work, const Corner &corner) -> Mitre
    {
      const std::array<Piece, 2> sides = Outer(work, corner);
      const Piece                back  = Reversed(sides[1]);
      const auto                 gone  = [](const Piece &piece) -> bool
      {
        return piece.Round && !(piece.Radius > 0.0);
      };
      /**
       *  Two moved pieces that end within the reach of each other, where points are one point,
       *  leave no room for a mitre: the corner is round.
       */
      const Vector gap = sides[1].Start - sides[0].End;
      if (gone(sides[0])
          || gone(sides[1])
          || !(Long(gap) > work.Reach)
          || !Wanted(work, corner, sides))
      {
        return {};
      }
      if (work.Asked.Gaps == Gap::Extend)
      {
        const std::pair<Point, bool> met = Meeting(work, corner, sides);
        return {.Count = met.second ? 1U : 0U,
                .At    = {met.first},
                .Runs  = {Run(sides[0], met.first, work.Reach), Run(back, met.first, work.Reach)}};
      }
      /**
       *  The way out of the corner halves its turn, in two doubles a number: square to the sum of
       *  the headings in and out, toward the corner's outer side; but where the shape turns back on
       *  itself that sum comes to nothing, and the way out is the heading in less the heading out,
       *  which does not. The line across passes the corner moved that way by the distance, square
       *  to it.
       */
      const std::array<Twofold::Double, 2> into = Unit(Headway(work, corner.Before, true));
      const std::array<Twofold::Double, 2> away = Unit(Headway(work, corner.After, false));
      const bool   turned = (into[0].High * away[0].High) + (into[1].High * away[1].High) < 0.0;
      const double side   = corner.Right ? 1.0 : -1.0;
      const std::array<Twofold::Double, 2> sum = {Twofold::Add(into[0], away[0]),
                                                  Twofold::Add(into[1], away[1])};
      const std::array<Twofold::Double, 2> out =
          turned
              ? Unit({Twofold::Difference(into[0], away[0]), Twofold::Difference(into[1], away[1])})
              : Unit({Twofold::Multiply(sum[1], Twofold::Lifted(side)),
                      Twofold::Multiply(sum[0], Twofold::Lifted(-side))});
      const Point           about    = work.Shape[corner.After].Start;
      const Twofold::Double distance = Twofold::Lifted(work.Asked.Distance);
      const Carrier         wide     = {
          .Round = false,
          .X     = Twofold::Add(Twofold::Lifted(about.X), Twofold::Multiply(out[0], distance)),
          .Y     = Twofold::Add(Twofold::Lifted(about.Y), Twofold::Multiply(out[1], distance)),
          .WayX  = Twofold::Negated(out[1]),
          .WayY  = out[0]};
      const Point foot = {.X = wide.X.High, .Y = wide.Y.High};
      /**
       *  The stretch of it doubles meet the pieces with is as long as the shape is large: one of a
       *  length of 1, about a foot thousands of millions out, is two points a few steps of a double
       *  apart, which tell no way.
       */
      const Vector along  = Vector{.X = wide.WayX.High, .Y = wide.WayY.High} * work.Scale;
      const Across across = {.Drawn = {.Start = foot - along, .End = foot + along}, .Wide = wide};
      const List<std::size_t>     &curves = corner.Right ? work.Rights : work.Lefts;
      const std::pair<Point, bool> first = Reached(work, {sides[0], curves[corner.Before]}, across);
      const std::pair<Point, bool> second = Reached(work, {back, curves[corner.After]}, across);
      return {.Count = first.second && second.second ? 2U : 0U,
              .At    = {first.first, second.first},
              .Runs = {Run(sides[0], first.first, work.Reach), Run(back, second.first, work.Reach)},
              .Wide = wide};
    }

    /**
     *  How the ends of the last mitre's edges are made, told: where the edge leaves and rejoins the
     *  moved pieces, the corner's foot on the first and the last of them; between, where each meets
     *  the next.
     *
     *  @param work The curves and the wedges.
     *  @param ends Where the edge leaves the first moved piece and rejoins the second.
     */
    void Recited(Work &work, const std::array<Point, 2> &ends)
    {
      const Wedge       wedge = work.Wedges.back();
      const std::size_t count = wedge.Edges.size() - 2;
      for (std::size_t index = 1; index <= count; ++index)
      {
        const std::size_t curve = wedge.Edges[index];
        for (const bool end : {false, true})
        {
          const Point where = end ? work.Curves[curve].End : work.Curves[curve].Start;
          Recipe      by    = {};
          if (where == ends[0])
          {
            by = {.Kind   = Making::Foot,
                  .First  = work.Curves[wedge.Edges[1]].Source,
                  .Second = 1,
                  .About  = wedge.About};
          }
          else if (where == ends[1])
          {
            by = {.Kind   = Making::Foot,
                  .First  = work.Curves[wedge.Edges[count]].Source,
                  .Second = 0,
                  .About  = wedge.About};
          }
          else
          {
            /**
             *  Where it meets the edge beside it on that side: the one before where the point is
             *  the rim piece's start, the one after where its end.
             */
            const bool before = where == wedge.Rim[index].Start;
            by                = {.Kind   = Making::Hit,
                                 .First  = curve,
                                 .Second = wedge.Edges[before ? index - 1 : index + 1]};
          }
          Tell(work, {curve, end}, by);
        }
      }
    }

    /**
     *  A corner's mitre put among the wedges, and its outer edge among the curves.
     *
     *  @param work   The shape; the wedge and its outer edge are added.
     *  @param corner The corner.
     *  @param mitre  Its mitre, of one turn or two.
     */
    void Mitred(Work &work, const Corner &corner, const Mitre &mitre)
    {
      const std::array<Piece, 2> sides    = Outer(work, corner);
      const Point                leaving  = sides[0].End;
      const Point                entering = sides[1].Start;
      if (mitre.Count == 1)
      {
        Wedged(work, corner,
               std::array{Stretched(sides[0], leaving, mitre.At[0]),
                          Stretched(sides[1], mitre.At[0], entering)});
        Recited(work, {leaving, entering});
        return;
      }
      Wedged(work, corner,
             std::array{Stretched(sides[0], leaving, mitre.At[0]),
                        Piece{.Start = mitre.At[0], .End = mitre.At[1]},
                        Stretched(sides[1], mitre.At[1], entering)});
      /**
       *  The line across keeps the carrier it was drawn by: lifted from its two rounded ends it
       *  would be another line, a step of a double off.
       */
      Borne(work, work.Wedges.back().Edges[2], mitre.Wide);
      Recited(work, {leaving, entering});
    }

    /**
     *  The signed curvature of a piece.
     *
     *  @param  piece The piece.
     *  @return       0 for a line; for an arc 1 over its radius, positive counterclockwise.
     */
    auto Bent(const Piece &piece) -> double
    {
      return piece.Round ? (piece.Left ? 1.0 : -1.0) / piece.Radius : 0.0;
    }

    /**
     *  Which way the shape turns at a corner, read in two doubles a number from its vertices and
     *  bulges: the cross product of the headings into it and out of it.
     *
     *  @param  work   The shape.
     *  @param  before The piece that ends at the corner.
     *  @param  after  The piece that starts at it.
     *  @return        Above 0 a turn to the left, below 0 one to the right; 0 where two doubles a
     *                 number tell no turn.
     */
    auto Veer(const Work &work, std::size_t before, std::size_t after) -> double
    {
      const std::array<Twofold::Double, 2> into = Headway(work, before, true);
      const std::array<Twofold::Double, 2> out  = Headway(work, after, false);
      const double cross = Twofold::Difference(Twofold::Multiply(into[0], out[1]),
                                               Twofold::Multiply(into[1], out[0]))
                               .High;
      const double size  = (std::abs(into[0].High) + std::abs(into[1].High))
                         * (std::abs(out[0].High) + std::abs(out[1].High));
      return std::abs(cross) > 1e-26 * size ? cross : 0.0;
    }

    /**
     *  Which side of a corner that turns right back, as near as doubles tell, is its outer one.
     *
     *  The turn is read again in two doubles a number: where the headings are not just opposite, it
     *  is the turn they make, however slight. Where they are, the corner turns the way its two
     *  pieces bend apart, its outer side the left where their curvatures sum counterclockwise; and
     *  where they bend neither way, a piece run back over the last, it has no outer side.
     *
     *  @param  work   The shape.
     *  @param  before The piece that ends at the corner.
     *  @param  after  The piece that starts at it.
     *  @return        Above 0 the right, below 0 the left, 0 neither.
     */
    auto Hairpin(const Work &work, std::size_t before, std::size_t after) -> int
    {
      const double cross = Veer(work, before, after);
      if (std::abs(cross) > 0.0)
      {
        return cross > 0.0 ? 1 : -1;
      }
      const double bend = Bent(work.Shape[before]) + Bent(work.Shape[after]);
      return bend > 0.0 ? -1 : static_cast<int>(bend < 0.0);
    }

    /**
     *  A vertex the shape turns at, as its circle or its mitre is made.
     */
    struct Turn
    {
      Corner                Bend   = {};    /**< The corner and its outer side.                */
      std::array<Vector, 2> Rights = {};    /**< The right normals of the headings in and out. */
      bool                  Round  = false; /**< Whether it has a circle on its right: it turns
                                                 left, or right back.                            */
      bool        Other = false;            /**< Whether it has one on its left: it turns right. */
      std::int8_t Flank = 0;                /**< The side its circle on the right is on: the right,
                                                 or where it turns right back its outer side, or
                                                 neither.                                        */
      Mitre Tip = {};                       /**< Its mitre, of no turn where it has none. */
    };

    /**
     *  The vertices the shape turns at, each with its mitre where corners are extended or
     *  chamfered. Of a shape offset on one side alone, those whose outer side is the other are left
     *  out.
     *
     *  @param  work The shape.
     *  @return      The turns, in order.
     */
    auto Turns(Work &work) -> List<Turn>
    {
      const std::size_t count = work.Shape.size();
      List<Turn>        turns;
      turns.reserve(count);
      work.Even.assign(count, 0);
      for (std::size_t index = work.Closed ? 0 : 1; index < count; ++index)
      {
        const std::size_t before = (index + count - 1) % count;
        const Vector      into   = Aimed(Heading(work.Shape[before], work.Shape[before].End));
        const Vector      out    = Aimed(Heading(work.Shape[index], work.Shape[index].Start));
        /**
         *  The turn as the same doubles tell it that place the moved pieces' ends: any at all is a
         *  turn, and none is none.
         */
        const double turn    = Cross(into, out);
        const bool   back    = !(turn > 0.0) && !(turn < 0.0) && Dot(into, out) < 0.0;
        const int    hairpin = back ? Hairpin(work, before, index) : 0;
        const bool   turned  = turn > 0.0 || turn < 0.0 || hairpin != 0;
        work.Even[index]     = static_cast<std::uint8_t>(!turned && !back);
        Turn made = {.Bend   = {.Before = before,
                                .After  = index,
                                .Right  = back ? hairpin > 0 : turn > 0.0,
                                .Into   = into},
                     .Rights = {Vector{.X = into.Y, .Y = -into.X}, Vector{.X = out.Y, .Y = -out.X}},
                     .Round  = turn > 0.0 || back,
                     .Other  = turn < 0.0,
                     .Flank  = static_cast<std::int8_t>(back ? std::clamp(hairpin, -1, 1) : 1)};
        if (work.Keeps != 0 && (work.Keeps > 0) != made.Bend.Right)
        {
          continue;
        }
        /**
         *  One side of an open shape asked, only the corners whose outer side it is have mitres.
         */
        const bool sided = work.Flank == 0 || (work.Flank > 0) == made.Bend.Right;
        if (work.Asked.Gaps != Gap::Fillet && turned && sided)
        {
          made.Tip = Mitring(work, made.Bend);
        }
        turns.push_back(made);
      }
      return turns;
    }

    /**
     *  An arc may not come round over itself: where an arc of the shape has a mitre at each end on
     *  one side, and its run on past its end, its run back before its start and the arc itself
     *  together pass a whole turn of its moved circle, the mitre with the longer run along it is
     *  taken away, and that corner is round.
     *
     *  @param work  The shape.
     *  @param turns The turns, whose mitres are taken away.
     */
    void Unwound(const Work &work, List<Turn> &turns)
    {
      for (std::size_t index = 0; index < turns.size(); ++index)
      {
        Turn        &start = turns[index];
        Turn        &end   = turns[(index + 1) % turns.size()];
        const Piece &arc   = work.Shape[start.Bend.After];
        if (!arc.Round
            || arc.Shallow
            || end.Bend.Before != start.Bend.After
            || start.Tip.Count == 0
            || end.Tip.Count == 0
            || start.Bend.Right != end.Bend.Right)
        {
          continue;
        }
        const double radius = Outer(work, start.Bend)[1].Radius;
        /**
         *  Past a whole turn by more than a billionth of it: three that make just one whole turn,
         *  as near as doubles tell, leave both mitres.
         */
        const double sweep = Swung(arc.Start - arc.Center, arc.End - arc.Center, arc.Left);
        if (end.Tip.Runs[0] + start.Tip.Runs[1] + (radius * sweep) > Whole * radius * (1.0 + 1e-13))
        {
          (end.Tip.Runs[0] >= start.Tip.Runs[1] ? end : start).Tip = {};
        }
      }
    }
  } // namespace

  auto Headway(const Work &work, std::size_t piece, bool ends) -> std::array<Twofold::Double, 2>
  {
    const Vertex &from  = work.Vertices[piece];
    const Point   onto  = work.Vertices[(piece + 1) % work.Vertices.size()].Position;
    const auto    apart = [](double end, double start) -> Twofold::Double
    {
      const Twofold::Two sum = Twofold::TwoSum(end, -start);
      return {.High = sum.High, .Low = sum.Low};
    };
    const Twofold::Double chordX = apart(onto.X, from.Position.X);
    const Twofold::Double chordY = apart(onto.Y, from.Position.Y);
    /**
     *  A line heads along itself.
     */
    if (!(from.Bulge < 0.0) && !(from.Bulge > 0.0))
    {
      return {chordX, chordY};
    }
    const Twofold::Two    square = Twofold::TwoProduct(from.Bulge, from.Bulge);
    const Twofold::Double along =
        Twofold::Difference(Twofold::Lifted(1.0), {.High = square.High, .Low = square.Low});
    const Twofold::Double side = Twofold::Lifted((ends ? 2.0 : -2.0) * from.Bulge);
    return {Twofold::Difference(Twofold::Multiply(along, chordX), Twofold::Multiply(side, chordY)),
            Twofold::Add(Twofold::Multiply(side, chordX), Twofold::Multiply(along, chordY))};
  }

  auto Unit(const std::array<Twofold::Double, 2> &way) -> std::array<Twofold::Double, 2>
  {
    const Twofold::Double length = Twofold::Root(
        Twofold::Add(Twofold::Multiply(way[0], way[0]), Twofold::Multiply(way[1], way[1])));
    return {Twofold::Divide(way[0], length), Twofold::Divide(way[1], length)};
  }

  void Cornered(Work &work)
  {
    List<Turn> turns = Turns(work);
    Unwound(work, turns);
    for (const Turn &turn : turns)
    {
      const std::size_t before = turn.Bend.Before;
      const std::size_t after  = turn.Bend.After;
      const Point       about  = work.Shape[after].Start;
      if (turn.Tip.Count > 0)
      {
        Mitred(work, turn.Bend, turn.Tip);
      }
      else if (turn.Round)
      {
        Hinged(work, {.About  = about,
                      .Sides  = turn.Rights,
                      .Beside = {before, after},
                      .Joined = {work.Rights[before], work.Rights[after]},
                      .Ended  = {before, after},
                      .Tails  = {1, 0},
                      .Flank  = turn.Flank});
      }
      else if (turn.Other)
      {
        Hinged(work, {.About  = about,
                      .Sides  = {-turn.Rights[1], -turn.Rights[0]},
                      .Beside = {before, after},
                      .Joined = {work.Lefts[after], work.Lefts[before]},
                      .Ended  = {after, before},
                      .Tails  = {0, 1},
                      .Flank  = -1});
      }
    }
  }
} // namespace Geoxy::Belt
