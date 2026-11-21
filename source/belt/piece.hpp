/***************************************************************************************************
 *  @file      source/belt/piece.hpp
 *  @brief     A line or an arc as a region's offset asks of it many times over: its points and
 *             headings along it, where a point of its line or circle lies along it, its box and how
 *             far a point is from it, each in a few operations of doubles.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BELT_PIECE_HPP
#define GEOXY_BELT_PIECE_HPP

#include "geoxy/box.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>

namespace Geoxy::Belt
{
  /**
   *  A whole turn, in radians.
   */
  inline constexpr double Whole = 2.0 * std::numbers::pi;

  /**
   *  No piece: what a number of one holds where there is none.
   */
  inline constexpr std::size_t None = std::numeric_limits<std::size_t>::max();

  /**
   *  A line from its start to its end, or an arc between them about a center.
   */
  struct Piece
  {
    bool        Round  = false; /**< Whether it is an arc.                 */
    Point       Start  = {};    /**< Where it starts.                      */
    Point       End    = {};    /**< Where it ends.                        */
    Point       Center = {};    /**< An arc's center.                      */
    double      Radius = 0.0;   /**< An arc's radius.                      */
    bool        Left   = true;  /**< Whether an arc runs counterclockwise. */
    std::size_t Source = None;  /**< The shape's piece it is moved from; none for a circle
                                     about a vertex.                                     */
    std::array<std::size_t, 2> Spared = {None, None}; /**< The shape's pieces it lies the
                                                           distance from by its making, which
                                                           it is not measured against.     */
    std::array<std::size_t, 2> Joined = {None, None}; /**< A circle about a vertex: the two
                                                           curves it joins, which it meets at
                                                           its ends alone.                 */
    std::array<std::size_t, 2> Beside = {None, None}; /**< A curve about a vertex: the arcs of
                                                           the shape that meet there, measured
                                                           against only where they curl nearer
                                                           than the vertex itself.           */
    Point       Pivot  = {};                          /**< That vertex. */
    std::size_t Bounds = None; /**< The mitre or the square it is an edge of, which it is
                                    not asked whether it lies inside; none for any other. */
    bool Shallow = false;      /**< Whether it is an arc so nearly straight that its center
                                    lies too far off to measure from: it is measured from
                                    its ends, its radius and the way it bends instead.    */
    std::int8_t Flank = 0;     /**< Which side of the shape it lies on by its making: 1 its
                                    right, -1 its left, 0 neither.                        */
  };

  /**
   *  The length of a displacement, by a plain root: its squares are far from overflow at the sizes
   *  a drawing has.
   *
   *  @param  way The displacement.
   *  @return     Its length.
   */
  inline auto Long(Vector way) -> double
  {
    return std::sqrt(Dot(way, way));
  }

  /**
   *  An arc, marked one.
   *
   *  @param  arc The arc, its ends, center, radius and way set.
   *  @return     The arc.
   */
  inline auto Arced(Piece arc) -> Piece
  {
    arc.Round = true;
    return arc;
  }

  /**
   *  A piece run the other way.
   *
   *  @param  piece The piece.
   *  @return       The same line or arc from its end to its start.
   */
  inline auto Reversed(Piece piece) -> Piece
  {
    std::swap(piece.Start, piece.End);
    piece.Left = !piece.Left;
    return piece;
  }

  /**
   *  How far round one way lies from another, in radians: worked out only where a length along a
   *  circle is asked.
   *
   *  @param  from The way started from.
   *  @param  onto The way come to.
   *  @param  left Whether counterclockwise.
   *  @return      0 up to a whole turn.
   */
  inline auto Swung(Vector from, Vector onto, bool left) -> double
  {
    const double turn = std::atan2(left ? Cross(from, onto) : -Cross(from, onto), Dot(from, onto));
    return turn < 0.0 ? turn + Whole : turn;
  }

  /**
   *  How far round one way lies from another, as a number that grows with the angle and is worked
   *  out with no angle: a quarter turn a unit, from their dot and cross products alone.
   *
   *  @param  from The way started from, not of no length.
   *  @param  onto The way come to, not of no length.
   *  @param  left Whether counterclockwise.
   *  @return      0 up to 4, a whole turn.
   */
  inline auto Quartered(Vector from, Vector onto, bool left) -> double
  {
    const double along = Dot(from, onto);
    const double aside = left ? Cross(from, onto) : -Cross(from, onto);
    if (aside >= 0.0)
    {
      return along >= 0.0 ? aside / (along + aside) : 1.0 - (along / (aside - along));
    }
    return along < 0.0 ? 2.0 + (aside / (along + aside)) : 3.0 + (along / (along - aside));
  }

  /**
   *  A nearly straight arc's bulge as it runs, from its chord and its radius: the tangent of a
   *  quarter of its sweep, positive counterclockwise.
   *
   *  @param  arc The arc, of less than half a turn.
   *  @return     Its bulge.
   */
  inline auto Bowed(const Piece &arc) -> double
  {
    const double sine  = Long(arc.End - arc.Start) / (2.0 * arc.Radius);
    const double bulge = sine / (1.0 + std::sqrt((1.0 - sine) * (1.0 + sine)));
    return arc.Left ? bulge : -bulge;
  }

  /**
   *  The ways a nearly straight arc heads at its start and at its end, from its chord and its
   *  bulge: the chord turned back and on by half the arc's sweep.
   *
   *  @param  arc The arc.
   *  @return     The two headings, each as long as the chord.
   */
  inline auto Tangents(const Piece &arc) -> std::array<Vector, 2>
  {
    const double bulge  = Bowed(arc);
    const double whole  = 1.0 + (bulge * bulge);
    const double cosine = (1.0 - (bulge * bulge)) / whole;
    const double sine   = 2.0 * bulge / whole;
    const Vector chord  = arc.End - arc.Start;
    return {Vector{.X = (chord.X * cosine) + (chord.Y * sine),
                   .Y = (chord.Y * cosine) - (chord.X * sine)},
            Vector{.X = (chord.X * cosine) - (chord.Y * sine),
                   .Y = (chord.Y * cosine) + (chord.X * sine)}};
  }

  /**
   *  How far outside a nearly straight arc's circle a point lies, below 0 inside it: from the
   *  point's power about the circle, which the arc's ends and bulge give with nothing subtracted
   *  that is as large as the radius.
   *
   *  @param  arc   The arc.
   *  @param  point The point.
   *  @return       The point's distance from the center less the radius.
   */
  inline auto Outside(const Piece &arc, Point point) -> double
  {
    const double bulge = Bowed(arc);
    const Vector from  = point - arc.Start;
    const double power =
        Dot(from, point - arc.End)
        - (((1.0 - (bulge * bulge)) / (2.0 * bulge)) * Cross(arc.End - arc.Start, from));
    return power / (arc.Radius + std::sqrt(std::max(0.0, (arc.Radius * arc.Radius) + power)));
  }

  /**
   *  Where a point of a piece's line or circle lies along it: on a line its place; on an arc a
   *  number that grows with the angle from its start, worked out with no angle.
   *
   *  @param  piece The piece, an arc's two ends not one way from its center.
   *  @param  point The point.
   *  @return       0 at its start, 1 at its end, between them within it; a hair before an arc's
   *                start a hair below 0, not nearly a whole turn on.
   */
  inline auto Along(const Piece &piece, Point point) -> double
  {
    if (!piece.Round)
    {
      const Vector way = piece.End - piece.Start;
      return Dot(point - piece.Start, way) / Dot(way, way);
    }
    if (piece.Shallow)
    {
      const Vector chord = piece.End - piece.Start;
      return Dot(point - piece.Start, chord) / Dot(chord, chord);
    }
    const Vector start = piece.Start - piece.Center;
    const double span  = Quartered(start, piece.End - piece.Center, piece.Left);
    const double turn  = Quartered(start, point - piece.Center, piece.Left);
    return (turn > (span + 4.0) / 2.0 ? turn - 4.0 : turn) / span;
  }

  /**
   *  The way a piece heads at a point of it, with no angle worked out.
   *
   *  @param  piece The piece.
   *  @param  point The point.
   *  @return       The heading, of any length.
   */
  inline auto Heading(const Piece &piece, Point point) -> Vector
  {
    if (!piece.Round)
    {
      return piece.End - piece.Start;
    }
    const Vector out = point - piece.Center;
    return piece.Left ? Vector{.X = -out.Y, .Y = out.X} : Vector{.X = out.Y, .Y = -out.X};
  }

  /**
   *  Whether a point other than an arc's center lies in the directions the arc sweeps from it.
   *
   *  @param  arc   The arc.
   *  @param  point The point.
   *  @return       Whether it does; never for the center itself.
   */
  inline auto Facing(const Piece &arc, Point point) -> bool
  {
    if (arc.Shallow)
    {
      /**
       *  Ahead of the line square to it at its start and behind the one at its end.
       */
      const std::array<Vector, 2> ways = Tangents(arc);
      return Dot(point - arc.Start, ways[0]) >= 0.0 && Dot(point - arc.End, ways[1]) <= 0.0;
    }
    const Vector off = point - arc.Center;
    if (!(Dot(off, off) > 0.0))
    {
      return false;
    }
    /**
     *  By which side of the ways to the arc's two ends the point lies, with no angle worked out: an
     *  arc of half a turn or less, its end no way clockwise of its start, is swept between them, a
     *  longer one all round but between.
     */
    const Vector from   = (arc.Left ? arc.Start : arc.End) - arc.Center;
    const Vector onto   = (arc.Left ? arc.End : arc.Start) - arc.Center;
    const bool   after  = Cross(from, off) >= 0.0;
    const bool   before = Cross(off, onto) >= 0.0;
    return Cross(from, onto) >= 0.0 ? after && before : after || before;
  }

  /**
   *  The middle of a stretch of a piece between two points of it, with no angle worked out.
   *
   *  @param  piece The piece.
   *  @param  ends  Where the stretch starts and where it ends, as the piece runs.
   *  @param  share How much of the piece the stretch is, by where its ends lie along it: read only
   *                for an arc's stretch whose ends lie one way from the center, which is then all
   *                but none of the circle or all but the whole of it.
   *  @return       The point of the piece midway between the two.
   */
  inline auto Mid(const Piece &piece, const std::array<Point, 2> &ends, double share) -> Point
  {
    if (!piece.Round)
    {
      return ends[0] + ((ends[1] - ends[0]) / 2.0);
    }
    if (piece.Shallow)
    {
      /**
       *  The middle of the chord between the two, moved out by how far the arc rises from it.
       */
      const Vector chord = ends[1] - ends[0];
      const double half  = Long(chord) / 2.0;
      const double rise =
          (half * half) / (piece.Radius + std::sqrt((piece.Radius - half) * (piece.Radius + half)));
      const Vector right = Vector{.X = chord.Y, .Y = -chord.X} / (2.0 * half);
      return ends[0] + (chord / 2.0) + (right * (piece.Left ? rise : -rise));
    }
    const Vector from  = ends[0] - piece.Center;
    const Vector onto  = ends[1] - piece.Center;
    const Vector first = from / Long(from);
    const Vector last  = onto / Long(onto);
    const double aside = piece.Left ? Cross(first, last) : -Cross(first, last);
    if (std::abs(aside) <= 1e-9 && Dot(first, last) < 0.0)
    {
      /**
       *  Half a turn: a quarter turn on from its start.
       */
      const Vector square =
          piece.Left ? Vector{.X = -first.Y, .Y = first.X} : Vector{.X = first.Y, .Y = -first.X};
      return piece.Center + (square * piece.Radius);
    }
    const double span = Quartered(piece.Start - piece.Center, piece.End - piece.Center, piece.Left);
    const bool   wide = std::abs(aside) <= 1e-9 ? share * span > 2.0 : aside < 0.0;
    const Vector sum  = wide ? -(first + last) : first + last;
    return piece.Center + (sum * (piece.Radius / Long(sum)));
  }

  /**
   *  A piece's box.
   *
   *  @param  piece The piece.
   *  @return       The least box that holds it.
   */
  inline auto Hulled(const Piece &piece) -> Box
  {
    Box box = {
        .Min = {.X = std::min(piece.Start.X, piece.End.X),
                .Y = std::min(piece.Start.Y, piece.End.Y)},
        .Max = {.X = std::max(piece.Start.X, piece.End.X),
                .Y = std::max(piece.Start.Y, piece.End.Y)},
    };
    /**
     *  An arc reaches past its ends at each of its circle's rightmost, topmost, leftmost and lowest
     *  points that it passes.
     */
    static constexpr std::array<Vector, 4> Ways = {
        Vector{.X =  1.0, .Y =  0.0},
        Vector{.X =  0.0, .Y =  1.0},
        Vector{.X = -1.0, .Y =  0.0},
        Vector{.X =  0.0, .Y = -1.0},
    };
    if (piece.Shallow)
    {
      /**
       *  Its ends' box, grown by how far it rises from its chord.
       */
      const Vector rise = Vector{.X = 1.0, .Y = 1.0}
                        * (Long(piece.End - piece.Start) * std::abs(Bowed(piece)) / 2.0);
      return {.Min = box.Min - rise, .Max = box.Max + rise};
    }
    for (std::size_t quarter = 0; piece.Round && quarter < Ways.size(); ++quarter)
    {
      const Point reach = piece.Center + (Ways.at(quarter) * piece.Radius);
      if (Facing(piece, reach))
      {
        box.Min = {.X = std::min(box.Min.X, reach.X), .Y = std::min(box.Min.Y, reach.Y)};
        box.Max = {.X = std::max(box.Max.X, reach.X), .Y = std::max(box.Max.Y, reach.Y)};
      }
    }
    return box;
  }

  /**
   *  How far a point is from a piece, its ends counted.
   *
   *  @param  piece The piece.
   *  @param  point The point.
   *  @return       The distance to the piece's nearest point.
   */
  inline auto Apart(const Piece &piece, Point point) -> double
  {
    if (!piece.Round)
    {
      const Vector way   = piece.End - piece.Start;
      const double along = std::clamp(Dot(point - piece.Start, way) / Dot(way, way), 0.0, 1.0);
      return Long(point - (piece.Start + (way * along)));
    }
    if (Facing(piece, point))
    {
      return piece.Shallow ? std::abs(Outside(piece, point))
                           : std::abs(Long(point - piece.Center) - piece.Radius);
    }
    return std::min(Long(point - piece.Start), Long(point - piece.End));
  }

  /**
   *  The way to a point from a piece's nearest point, its ends counted.
   *
   *  @param  piece The piece.
   *  @param  point The point.
   *  @return       The way, as long as the point is far from the piece.
   */
  inline auto Away(const Piece &piece, Point point) -> Vector
  {
    if (!piece.Round)
    {
      const Vector way   = piece.End - piece.Start;
      const double along = std::clamp(Dot(point - piece.Start, way) / Dot(way, way), 0.0, 1.0);
      return point - (piece.Start + (way * along));
    }
    if (Facing(piece, point))
    {
      const Vector out = point - piece.Center;
      const double far = Long(out);
      return out * ((piece.Shallow ? Outside(piece, point) : far - piece.Radius) / far);
    }
    return Long(point - piece.Start) < Long(point - piece.End) ? point - piece.Start
                                                               : point - piece.End;
  }

  /**
   *  About how long a piece is, with no angle worked out: a line's length; an arc's within a fifth
   *  of it, by the number that grows with its angle. For what counts as a hair of a piece.
   *
   *  @param  piece The piece.
   *  @return       Its length, or about it.
   */
  inline auto Long(const Piece &piece) -> double
  {
    if (!piece.Round || piece.Shallow)
    {
      return Long(piece.End - piece.Start);
    }
    return piece.Radius
         * (Whole / 4.0)
         * Quartered(piece.Start - piece.Center, piece.End - piece.Center, piece.Left);
  }
} // namespace Geoxy::Belt

#endif /* GEOXY_BELT_PIECE_HPP */
