/***************************************************************************************************
 *  @file      source/belt/end.cpp
 *  @brief     The ends of an open polyline's band: round, flat or square.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "belt/end.hpp"

#include "belt/carrier.hpp"
#include "belt/corner.hpp"
#include "belt/make.hpp"
#include "belt/piece.hpp"
#include "belt/work.hpp"
#include "scratch.hpp"
#include "twofold.hpp"

#include "geoxy/box.hpp"
#include "geoxy/offset.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace Geoxy::Belt
{
  namespace
  {
    /**
     *  An end of an open shape.
     */
    struct Tip
    {
      Point       About = {}; /**< The end.                             */
      std::size_t Piece = 0;  /**< The piece that starts or ends there. */
      Vector      From  = {}; /**< The unit direction, square to the shape there, the band's
                                   edge round the end starts from: it ends the other way.    */
      Vector                         Out  = {}; /**< The unit heading out of the shape there. */
      std::array<Twofold::Double, 2> Wide = {}; /**< That heading in two doubles a number.    */
    };

    /**
     *  A line of an end's, by where it passes and which way it runs.
     */
    struct Sideline
    {
      double From  = 0.0;   /**< How many distances it passes from the end, sideways. */
      double Out   = 0.0;   /**< How many distances it passes from the end, outward.  */
      bool   Along = false; /**< Whether it runs outward, not sideways.               */
    };

    /**
     *  The carrier of a line of an end's, in two doubles a number: drawn from the end and the
     *  shape's heading there, not from the line's rounded ends.
     *
     *  @param  work The distance.
     *  @param  tip  The end.
     *  @param  edge The line.
     *  @return      Its carrier.
     */
    auto Laid(const Work &work, const Tip &tip, const Sideline &edge) -> Carrier
    {
      const std::array<Twofold::Double, 2> side = {tip.Wide[1], Twofold::Negated(tip.Wide[0])};
      const Twofold::Double                from = Twofold::Lifted(edge.From * work.Asked.Distance);
      const Twofold::Double                out  = Twofold::Lifted(edge.Out * work.Asked.Distance);
      const auto past = [&](double about, std::size_t axis) -> Twofold::Double
      {
        return Twofold::Add(Twofold::Lifted(about),
                            Twofold::Add(Twofold::Multiply(side.at(axis), from),
                                         Twofold::Multiply(tip.Wide.at(axis), out)));
      };
      const std::array<Twofold::Double, 2> &way = edge.Along ? tip.Wide : side;
      return {.Round = false,
              .X     = past(tip.About.X, 0),
              .Y     = past(tip.About.Y, 1),
              .WayX  = way[0],
              .WayY  = way[1]};
    }

    /**
     *  A point by an end, worked out in two doubles a number and kept: the end moved by so many
     *  distances sideways and outward.
     *
     *  @param  work The distance.
     *  @param  tip  The end.
     *  @param  by   How many distances sideways and outward.
     *  @return      The point's making.
     */
    auto Tipped(const Work &work, const Tip &tip, const Sideline &by) -> Recipe
    {
      const Carrier through = Laid(work, tip, by);
      return {.Kind = Making::Held, .About = {.X = through.X.High, .Y = through.Y.High}};
    }

    /**
     *  The line across an end, from one side of the shape to the other. What an arc sweeps ends at
     *  its center, so across an arc tighter than the distance the line stops there on that side.
     *
     *  @param work The shape; the curve is added.
     *  @param tip  The end.
     */
    void Flattened(Work &work, const Tip &tip)
    {
      const double distance = work.Asked.Distance;
      const Piece &piece    = work.Shape[tip.Piece];
      const bool   tight    = piece.Round && piece.Radius < distance;
      const bool   inner    = Dot(piece.Center - tip.About, tip.From) > 0.0;
      Piece across  = {.Start = tight && inner ? piece.Center : tip.About + (tip.From * distance),
                       .End   = tight && !inner ? piece.Center : tip.About - (tip.From * distance)};
      across.Spared = {tip.Piece, None};
      const std::size_t line = Put(work, across, 0.0);
      Borne(work, line, Laid(work, tip, {}));
      /**
       *  Its ends, where no arc's center, are the end moved the distance to each side.
       */
      if (!tight || !inner)
      {
        Tell(work, {line, false}, Tipped(work, tip, {.From = 1.0}));
      }
      if (!tight || inner)
      {
        Tell(work, {line, true}, Tipped(work, tip, {.From = -1.0}));
      }
    }

    /**
     *  The square past an end, as long as the distance: a part of the band, and its three far sides
     *  among the curves.
     *
     *  @param work The shape; the wedge and the curves are added.
     *  @param tip  The end.
     */
    void Squared(Work &work, const Tip &tip)
    {
      const double                  distance = work.Asked.Distance;
      const Vector                  out      = tip.Out * distance;
      const Point                   first    = tip.About + (tip.From * distance);
      const Point                   last     = tip.About - (tip.From * distance);
      const std::array<Point, 5>    corners  = {first, first + out, last + out, last, first};
      const std::array<Sideline, 3> edges    = {Sideline{.From = 1.0, .Along = true},
                                                Sideline{.Out = 1.0},
                                                Sideline{.From = -1.0, .Along = true}};
      const std::array<Sideline, 4> turns    = {
          Sideline{.From = 1.0}, Sideline{.From = 1.0, .Out = 1.0},
          Sideline{.From = -1.0, .Out = 1.0}, Sideline{.From = -1.0}};
      Wedge wedge;
      for (std::size_t side = 0; side < 4; ++side)
      {
        Piece edge = {.Start = corners.at(side), .End = corners.at(side + 1)};
        wedge.Rim.push_back(edge);
        wedge.Hull = Joined(wedge.Hull, Hulled(edge));
        /**
         *  The curve along each edge but the last, which runs through the end itself and is no edge
         *  of the region's.
         */
        wedge.Edges.push_back(None);
        if (side < 3)
        {
          edge.Bounds        = work.Wedges.size();
          wedge.Edges.back() = Put(work, edge, 0.0);
          Borne(work, wedge.Edges.back(), Laid(work, tip, edges.at(side)));
          Tell(work, {wedge.Edges.back(), false}, Tipped(work, tip, turns.at(side)));
          Tell(work, {wedge.Edges.back(), true}, Tipped(work, tip, turns.at(side + 1)));
        }
      }
      work.Wedges.push_back(std::move(wedge));
    }

    /**
     *  The whole circles about the vertices near an open shape's flat ends: what is within the
     *  distance of such a vertex reaches past the line across an end, where the end's own piece no
     *  longer holds it, so its edge there runs along the circle past the arc round the corner.
     *
     *  @param work The shape; the curves are added.
     */
    void Ringed(Work &work)
    {
      const double distance = work.Asked.Distance;
      const Vector reach    = {.X = distance};
      for (std::size_t index = 1; index < work.Shape.size(); ++index)
      {
        const Point about = work.Shape[index].Start;
        if (std::min(Long(about - work.Shape.front().Start), Long(about - work.Shape.back().End))
            < 2.0 * distance)
        {
          for (const double side : {1.0, -1.0})
          {
            Put(work,
                Arced({.Start  = about + (reach * side),
                       .End    = about - (reach * side),
                       .Center = about,
                       .Radius = distance,
                       .Left   = true}),
                0.0);
          }
        }
      }
    }
  } // namespace

  void Ended(Work &work)
  {
    if (work.Asked.Ends == Cap::Round)
    {
      Capped(work);
      return;
    }
    const Vector start = Aimed(Heading(work.Shape.front(), work.Shape.front().Start));
    const Vector end   = Aimed(Heading(work.Shape.back(), work.Shape.back().End));
    const std::array<Twofold::Double, 2> first = Unit(Headway(work, 0, false));
    const std::array<Tip, 2>             tips  = {
        Tip{.About = work.Shape.front().Start,
            .Piece = 0,
            .From  = {.X = -start.Y, .Y = start.X},
            .Out   = -start,
            .Wide  = {Twofold::Negated(first[0]), Twofold::Negated(first[1])}},
        Tip{.About = work.Shape.back().End,
            .Piece = work.Shape.size() - 1,
            .From  = {.X = end.Y, .Y = -end.X},
            .Out   = end,
            .Wide  = Unit(Headway(work, work.Shape.size() - 1, true))}};
    for (const Tip &tip : tips)
    {
      if (work.Asked.Ends == Cap::Flat)
      {
        Flattened(work, tip);
      }
      else
      {
        Squared(work, tip);
      }
    }
    Ringed(work);
  }
} // namespace Geoxy::Belt
