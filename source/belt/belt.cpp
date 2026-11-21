/***************************************************************************************************
 *  @file      source/belt/belt.cpp
 *  @brief     A polyline's offset as a region: the curves its edge can run along made, cut where
 *             they meet, kept where they lie on the edge, and joined.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "belt/belt.hpp"

#include "belt/carrier.hpp"
#include "belt/corner.hpp"
#include "belt/cut.hpp"
#include "belt/end.hpp"
#include "belt/join.hpp"
#include "belt/keep.hpp"
#include "belt/make.hpp"
#include "belt/piece.hpp"
#include "belt/work.hpp"
#include "scratch.hpp"

#include "geoxy/box.hpp"
#include "geoxy/offset.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

namespace Geoxy::Belt
{
  namespace
  {
    /**
     *  How many pieces a closed shape has at the most to be looked at for whether it crosses or
     *  touches itself, each piece's box tried against every other's.
     */
    constexpr std::size_t Looked = 128;

    /**
     *  Whether a point of two pieces' lines or circles falls within both.
     *
     *  @param  work  The reach.
     *  @param  pair  The pieces.
     *  @param  point The point.
     *  @return       Whether it does, or within the reach of an end of each.
     */
    auto Shared(const Work &work, const std::array<const Piece *, 2> &pair, Point point) -> bool
    {
      const auto within = [&work, point](const Piece &piece) -> bool
      {
        const double along = Along(piece, point);
        const double hair  = work.Reach / std::max(Long(piece), work.Reach);
        return along >= -hair && along <= 1.0 + hair;
      };
      return within(*pair[0]) && within(*pair[1]);
    }

    /**
     *  Whether two pieces of a shape along one line or one circle share more than a point: the end
     *  of the later lies within the earlier, or the middle of the earlier within the later. A later
     *  one that only starts within the earlier is come to by a piece that ends there, which meets
     *  the earlier first.
     *
     *  @param  work   The reach.
     *  @param  first  The earlier piece.
     *  @param  second The later, along the same line or circle.
     *  @return        Whether they do.
     */
    auto Lapped(const Work &work, const Piece &first, const Piece &second) -> bool
    {
      const auto inside = [&work](const Piece &piece, Point point) -> bool
      {
        const double along = Along(piece, point);
        const double hair  = work.Reach / std::max(Long(piece), work.Reach);
        return along > hair && along < 1.0 - hair;
      };
      return inside(first, second.End) || inside(second, Mid(first, {first.Start, first.End}, 1.0));
    }

    /**
     *  Whether two pieces of a closed shape meet as those of one that neither crosses nor touches
     *  itself may not: along one line or circle, at a point other than the vertex two in a row
     *  share, or so near touching that doubles cannot tell, but two in a row about that vertex.
     *
     *  @param  work   The shape.
     *  @param  first  One piece.
     *  @param  second A later one.
     *  @return        Whether they do.
     */
    auto Tangled(const Work &work, std::size_t first, std::size_t second) -> bool
    {
      const Piece &one   = work.Shape[first];
      const Piece &other = work.Shape[second];
      if (one.Round == other.Round && Together(work, one, other))
      {
        return Lapped(work, one, other);
      }
      const Hits  hits   = Met(one, other);
      const bool  next   = second == first + 1;
      const bool  joined = next || (first == 0 && second == work.Shape.size() - 1);
      const Point vertex = next ? one.End : one.Start;
      const auto  apart  = [&](Point point) -> bool
      {
        return (!joined || Long(point - vertex) > work.Reach)
            && Shared(work, {&one, &other}, point);
      };
      if (one.Shallow || other.Shallow)
      {
        /**
         *  A nearly straight arc is met in two doubles a number, wherever the other lies.
         */
        const Hits fine = Met(Carried(one, work.Vertices[first].Bulge),
                              Carried(other, work.Vertices[second].Bulge));
        return std::ranges::any_of(std::span{fine.At}.first(fine.Count), apart);
      }
      if (hits.Close)
      {
        /**
         *  Two in a row that run on smoothly, one tangent to the other or nearly, come near
         *  touching only about the vertex they share: past it the second runs ahead of where the
         *  first ends. Any other two so near touching may touch.
         */
        const Vector into = next ? Heading(one, one.End) : Heading(other, other.End);
        const Vector out  = next ? Heading(other, other.Start) : Heading(one, one.Start);
        const auto   size = [](const Piece &piece) -> double
        {
          return piece.Round ? piece.Radius : Long(piece);
        };
        const double near = 2e-3 * std::max(size(one), size(other));
        if (!joined || !(Dot(into, out) > 0.0))
        {
          return true;
        }
        /**
         *  Where they come nearest touching lies midway between the two points they share, the
         *  vertex one of them: near the vertex, so is the other.
         */
        if (Long(hits.Touch - vertex) <= near)
        {
          return false;
        }
        /**
         *  Two arcs in a row round nearly one circle meet again far from the vertex, where doubles
         *  place the point poorly: it is found in two doubles a number.
         */
        const Hits fine = Met(Carried(one, work.Vertices[first].Bulge),
                              Carried(other, work.Vertices[second].Bulge));
        return std::ranges::any_of(std::span{fine.At}.first(fine.Count), apart);
      }
      return std::ranges::any_of(std::span{hits.At}.first(hits.Count), apart);
    }

    /**
     *  Whether a closed shape neither crosses nor touches itself: no two of its pieces whose boxes
     *  come together meet but two in a row at the vertex they share.
     *
     *  @param  work The shape, closed.
     *  @return      Whether it is; false for one of more pieces than are looked at, or of fewer
     *               than three.
     */
    auto Plain(const Work &work) -> bool
    {
      const std::size_t count = work.Shape.size();
      if (count > Looked || count < 3)
      {
        return false;
      }
      const List<Box> &boxes = work.Hulls;
      for (std::size_t first = 0; first < count; ++first)
      {
        for (std::size_t second = first + 1; second < count; ++second)
        {
          if (!Parted(boxes[first], boxes[second], work.Reach) && Tangled(work, first, second))
          {
            return false;
          }
        }
      }
      return true;
    }

    /**
     *  The area a closed shape encloses, positive counterclockwise.
     *
     *  @param  work The shape.
     *  @return      The area: the polygon of its vertices by the shoelace formula, with each arc's
     *               segment.
     */
    auto Enclosed(const Work &work) -> double
    {
      double area = 0.0;
      for (std::size_t index = 0; index < work.Shape.size(); ++index)
      {
        const Piece &piece  = work.Shape[index];
        area               += Cross(piece.Start - Point{}, piece.End - Point{}) / 2.0;
        if (piece.Round)
        {
          /**
           *  The arc's sweep is four times the angle its bulge is the tangent of, and the sine of
           *  the sweep comes from the bulge with no more angles.
           */
          const double bulge = std::abs(work.Vertices[index].Bulge);
          const double once  = 1.0 + (bulge * bulge);
          const double sine  = 4.0 * bulge * (1.0 - (bulge * bulge)) / (once * once);
          /**
           *  A nearly straight arc's sweep less its sine is lost to rounding: its segment is a
           *  third of its bulge times its chord squared, to the bulge's fourth power.
           */
          const Vector chord = piece.End - piece.Start;
          const double segment =
              piece.Shallow ? Dot(chord, chord) * bulge * ((1.0 / 3.0) + (bulge * bulge / 15.0))
                            : piece.Radius * piece.Radius * ((4.0 * std::atan(bulge)) - sine) / 2.0;
          area += piece.Left ? segment : -segment;
        }
      }
      return area;
    }

    /**
     *  Where along one axis a shape is worked about, so that what counts as one point goes by its
     *  size and not by how far from 0 it is drawn, and so that moving it there rounds nothing: a
     *  whole number of a power of two no less than its extent, just below it. A shape within twice
     *  that power of 0 is worked where it is: moved, its numbers would no longer fit a double.
     *
     *  @param  least The least of its vertices along the axis.
     *  @param  most  The greatest.
     *  @return       The place: every vertex less it is a double, exactly.
     */
    auto Based(double least, double most) -> double
    {
      const double extent = most - least;
      if (!(extent > 0.0))
      {
        return least;
      }
      const double grain = std::exp2(std::ceil(std::log2(extent)) + 1.0);
      if (least <= 2.0 * grain && most >= -2.0 * grain)
      {
        return 0.0;
      }
      return (least > 0.0 ? std::floor(least / grain) : std::ceil(most / grain)) * grain;
    }

    /**
     *  Where a polyline is worked about.
     *
     *  @param  polyline The polyline, of one vertex or more.
     *  @return          The way from 0 to the place.
     */
    auto Origin(const Polyline &polyline) -> Vector
    {
      Point least = polyline.Vertices.front().Position;
      Point most  = least;
      for (const Vertex &vertex : polyline.Vertices)
      {
        least = {.X = std::min(least.X, vertex.Position.X),
                 .Y = std::min(least.Y, vertex.Position.Y)};
        most = {.X = std::max(most.X, vertex.Position.X), .Y = std::max(most.Y, vertex.Position.Y)};
      }
      return {.X = Based(least.X, most.X), .Y = Based(least.Y, most.Y)};
    }

    /**
     *  A polyline's offset, as a region or on one side.
     *
     *  @param  polyline  The polyline.
     *  @param  buffering The distance, and for a region the side, the filling, the corners and
     *                    ends.
     *  @param  flank     One side of an open polyline asked: 1 its right, -1 its left; 0 for a
     *                    region.
     *  @return           The region's loops, or the side's polylines.
     */
    auto Made(const Polyline &polyline, const Buffering &buffering, int flank)
        -> std::vector<Polyline>
    {
      const Vector origin = Origin(polyline);
      List<Vertex> shifted(polyline.Vertices.begin(), polyline.Vertices.end());
      for (Vertex &vertex : shifted)
      {
        vertex.Position = vertex.Position - origin;
      }
      Work work = {
          .Vertices = shifted, .Closed = polyline.Closed, .Asked = buffering, .Flank = flank};
      Shaped(work);
      Sized(work);
      /**
       *  A closed shape that neither crosses nor touches itself fills one side of itself: asked
       *  inward or outward, only that side's curves can lie on the region's edge.
       */
      if (work.Asked.Sides != Side::Both && Plain(work))
      {
        const bool outside = Enclosed(work) > 0.0;
        work.Keeps         = outside == (work.Asked.Sides == Side::Outward) ? 1 : -1;
      }
      Moved(work);
      Cornered(work);
      if (!work.Closed)
      {
        Ended(work);
      }
      /**
       *  About as many points as curves' ends and half as many again where they meet.
       */
      Nodes points(work.Reach, 3 * work.Curves.size());
      work.By.reserve(3 * work.Curves.size());
      const Marks marks = {.Curves = &work, .Points = &points};
      for (std::size_t curve = 0; curve < work.Curves.size(); ++curve)
      {
        const std::size_t start = Numbered(marks, work.Curves[curve].Start, Sourced(curve, false));
        const std::size_t end   = Numbered(marks, work.Curves[curve].End, Sourced(curve, true));
        work.Cuts.push_back({.Curve = curve, .Along = 0.0, .Node = start});
        work.Cuts.push_back({.Curve = curve, .Along = 1.0, .Node = end});
      }
      Met(marks);
      Kept(work, points);
      /**
       *  The points the edge turns at, each worked out again in two doubles a number by how it was
       *  made: the doubles nearest the point, where those that found it lay a few steps off.
       */
      work.Fined.assign(points.Count(), 0);
      for (const Stretch &stretch : work.Stretches)
      {
        for (const std::size_t node : {stretch.Start, stretch.End})
        {
          if (work.Fined[node] != 0)
          {
            continue;
          }
          work.Fined[node]  = 1;
          const Birth  from = work.By[node];
          const Recipe by =
              from.Other >= Starts
                  ? Ending(work, from.Curve, from.Other == Ends)
                  : Recipe{.Kind = Making::Hit, .First = from.Curve, .Second = from.Other};
          points.Move(node, Remade(work, by, points.At(node)));
        }
      }
      std::vector<Polyline> loops;
      if (flank != 0)
      {
        Flanking(work);
        loops = Sidelong(work, points);
      }
      else
      {
        Once(work);
        loops = Looped(work, points);
      }
      for (Polyline &loop : loops)
      {
        for (Vertex &vertex : loop.Vertices)
        {
          vertex.Position = vertex.Position + origin;
        }
      }
      return loops;
    }
  } // namespace

  auto Edged(const Polyline &polyline, const Buffering &buffering) -> std::vector<Polyline>
  {
    return Made(polyline, buffering, 0);
  }

  auto Flanked(const Polyline &polyline, double distance, Gap gaps) -> std::vector<Polyline>
  {
    return Made(polyline, {.Distance = std::abs(distance), .Sides = Side::Both, .Gaps = gaps},
                distance > 0.0 ? 1 : -1);
  }
} // namespace Geoxy::Belt
