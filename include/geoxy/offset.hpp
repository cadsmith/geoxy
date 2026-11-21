/***************************************************************************************************
 *  @file      geoxy/offset.hpp
 *  @brief     The curve at a distance from another, as OFFSET draws it: its outer corners closed by
 *             a gap kind, the loops it makes crossing itself cleaned away.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_OFFSET_HPP
#define GEOXY_OFFSET_HPP

#include "geoxy/accuracy.hpp"
#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/error.hpp"
#include "geoxy/line.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/xline.hpp"

#include <cstdint>
#include <expected>
#include <vector>

namespace Geoxy
{
  /**
   *  How the gap an offset opens at an outer corner is closed: AutoCAD's OFFSETGAPTYPE.
   *
   *  Where a curve turns away from the side it is offset to, its two pieces' offsets part, and the
   *  gap between their ends is closed one of three ways. Inner corners, where the offsets cross,
   *  are cut where they cross whatever the kind.
   */
  enum class Gap : std::uint8_t
  {
    Extend,  /**< Each piece run on, a line along its line and an arc round its circle, to where
                  they meet; OFFSETGAPTYPE 0, AutoCAD's default.                                  */
    Fillet,  /**< An arc about the corner, its radius the distance; OFFSETGAPTYPE 1. */
    Chamfer, /**< A line across, as far from the corner as the distance, each piece run on to it;
                  OFFSETGAPTYPE 2.                                                                */
  };

  /**
   *  How far to offset, and how: what an offset takes beside the curve and the tolerance.
   */
  struct Offsetting
  {
    double Distance = 0.0;                /**< How far: to the right of the way the curve runs when
                                               positive, so outside a counterclockwise loop.      */
    Gap      Gaps          = Gap::Extend; /**< How the gaps at outer corners are closed. */
    Accuracy Approximation = {};          /**< How near an ellipse's or a spline's offset, which
                                               is neither, must come; nothing else reads it.      */
  };

  /**
   *  Which way a closed curve is offset as a region: against what it fills.
   */
  enum class Side : std::uint8_t
  {
    Inward,  /**< What it fills, less everything within the distance of it. */
    Outward, /**< What it fills, with everything within the distance of it. */
    Both,    /**< Everything within the distance of it.                     */
  };

  /**
   *  What a closed curve that crosses itself fills.
   */
  enum class Fills : std::uint8_t
  {
    NonZero, /**< Every point it winds about, however many times.    */
    EvenOdd, /**< Every point it winds about an odd number of times. */
  };

  /**
   *  How an open curve's offset as a region closes round its ends.
   */
  enum class Cap : std::uint8_t
  {
    Round,  /**< A half circle about the end, its radius the distance. */
    Flat,   /**< A line across the end, square to the curve there.     */
    Square, /**< A square past the end, as long as the distance.       */
  };

  /**
   *  How far to offset a closed curve as a region, and how.
   */
  struct Buffering
  {
    double Distance = 0.0;            /**< How far: more than 0.                      */
    Side   Sides    = Side::Outward;  /**< Which way against what it fills.           */
    Fills  Fill     = Fills::NonZero; /**< What it fills where it crosses itself.     */
    Gap    Gaps     = Gap::Fillet;    /**< How the buffer's outer corners are closed. */
    Cap    Ends     = Cap::Round;     /**< How it closes round an open curve's ends.  */
  };

  /**
   *  A line moved sideways.
   *
   *  @param  line     The line.
   *  @param  distance How far, to the right of the way it runs when positive.
   *  @return          The line as long, its ends moved square to it; Degenerate for a line of no
   *                   length, which has no side, and NotFinite where a number is not finite.
   */
  [[nodiscard]] GEOXY_API auto Offset(const Line &line, double distance)
      -> std::expected<Line, Error>;

  /**
   *  A ray moved sideways.
   *
   *  @param  ray      The ray.
   *  @param  distance How far, to the right of its way when positive.
   *  @return          The ray, its start moved square to its way; Degenerate for a ray of no way,
   *                   and NotFinite where a number is not finite.
   */
  [[nodiscard]] GEOXY_API auto Offset(const Ray &ray, double distance) -> std::expected<Ray, Error>;

  /**
   *  An endless line moved sideways.
   *
   *  @param  xline    The endless line.
   *  @param  distance How far, to the right of its way when positive.
   *  @return          The endless line through the point moved square to its way; Degenerate for
   *                   one of no way, and NotFinite where a number is not finite.
   */
  [[nodiscard]] GEOXY_API auto Offset(const Xline &xline, double distance)
      -> std::expected<Xline, Error>;

  /**
   *  An arc at a distance: the arc of the same center and sweep, its radius grown or shrunk.
   *
   *  To the right of a counterclockwise arc is away from its center, so a positive distance grows
   *  it and a negative one shrinks it; the other way round for a clockwise one. Worked out exactly
   *  but for the radius's one rounding.
   *
   *  @param  arc      The arc, or a circle as an arc of a whole turn.
   *  @param  distance How far, to the right of the way it runs when positive.
   *  @return          The arc; TooLarge where the distance takes the radius to 0 or past it, where
   *                   the arc collapses to its center and nothing is left; Degenerate for an arc of
   *                   no radius or no sweep; NotFinite where a number is not finite.
   */
  [[nodiscard]] GEOXY_API auto Offset(const Arc &arc, double distance) -> std::expected<Arc, Error>;

  /**
   *  A polyline at a distance, on one side of it.
   *
   *  An open polyline's offset is the part of the edge of everything within the distance of it that
   *  lies on the side asked: along a segment moved to that side, an arc to the radius r ± d about
   *  its own center; round a corner whose outer side it is, or along that corner's mitre as
   *  Offsetting::Gaps says, a fillet where an extension or a chamfer finds nothing to meet; and
   *  never nearer the polyline than the distance. A stretch the polyline runs over twice, there and
   *  back, has neither side, nor has a corner that turns right back, nor the round about an end:
   *  nothing is left of a polyline run out and straight back. The pieces run the way the polyline
   *  runs, each open, or closed where it comes back to its start.
   *
   *  A closed polyline's offset is its offset as a region, as Buffer gives it: to its right is
   *  outside what it fills where it runs counterclockwise and inside where clockwise, by the sign
   *  of the area it encloses; what it fills where it crosses itself is every point it winds about.
   *  One that encloses no area, a stretch run there and back or a figure of eight, has no such
   *  side, and its offset either way is the edge of all that lies within the distance of it. The
   *  loops run with what is left on their left, whichever way the polyline runs.
   *
   *  A segment is offset as the arc its bulge makes it, however slight its bulge: one so nearly
   *  straight that its center lies thousands of chords off is measured from its ends, its radius
   *  and the way it bends, not from that center. Only a bulge under 2^-60, which lifts the arc from
   *  its chord by less than a step of a double, is a line's.
   *
   *  @param  polyline   The polyline, open or closed, its vertices as given, those within the
   *                     tolerance of the one before taken as one.
   *  @param  offsetting The distance and the gap kind.
   *  @param  tolerance  Within which two vertices of the polyline are one. The offset itself is
   *                     worked out to the polyline's own numbers, whatever the tolerance.
   *  @return            The offset's polylines; none where nothing is left. Degenerate for a
   *                     polyline of fewer than two distinct vertices; NotFinite where a number is
   *                     not finite. A distance of 0 gives the polyline.
   */
  [[nodiscard]] GEOXY_API auto Offset(const Polyline &polyline, const Offsetting &offsetting,
                                      const Tolerance &tolerance)
      -> std::expected<std::vector<Polyline>, Error>;

  /**
   *  A polyline's offset as a region: the edge of what a closed one fills less, or with, everything
   *  within the distance of it, or of that buffer alone; an open one's buffer.
   *
   *  The buffer is every point within the distance of the polyline, so its corners are round, and
   *  an open polyline's ends; with corners extended or chamfered, the mitre each outer corner then
   *  gets as well, where it points the way asked: outward, out of what the polyline fills, inward,
   *  into it, both ways, every one. A corner that turns right back is outer on the side its two
   *  pieces bend away from, and round where they bend neither way. An open polyline's ends are
   *  round; or flat, the buffer then what its segments sweep square to themselves with the disks
   *  about the vertices between them; or square, that with a square as long as the distance past
   *  each end. With F what a closed polyline fills, inward is F less the buffer, outward F with the
   *  buffer, both the buffer; an open polyline fills nothing and has its buffer alone. The offset
   *  is the region's edge.
   *
   *  The buffer's edge can only run along a segment moved to one side by the distance, a circle of
   *  the distance about a vertex or an end, the outer edge of a corner's mitre, or the line across
   *  or the square past an end: those curves are made, cut wherever two of them meet, and a stretch
   *  between two cuts is kept where its middle lies outside the buffer, no nearer the polyline than
   *  the distance and in no mitre or square, and on the side of the filling asked. What is kept is
   *  joined end to start into loops, each closed; a loop that comes to a point several leave takes
   *  the first met turning clockwise, so that loops touch and never cross.
   *
   *  @param  polyline  The polyline, its vertices within the tolerance of the one before taken as
   *                    one.
   *  @param  buffering The distance, the side, the filling and how outer corners are closed; for an
   *                    open polyline both sides, the filling not read.
   *  @param  tolerance Within which two points are one, as the one-sided offset takes it.
   *  @return           The region's loops, each run with the region on its left: counterclockwise
   *                    about it, clockwise about a hole; none where nothing is left. Degenerate for
   *                    an open polyline asked inward or outward, a polyline of fewer than two
   *                    distinct vertices, or a distance not more than 0; NotFinite where a number
   *                    is not finite.
   */
  [[nodiscard]] GEOXY_API auto Buffer(const Polyline &polyline, const Buffering &buffering,
                                      const Tolerance &tolerance)
      -> std::expected<std::vector<Polyline>, Error>;

  /**
   *  An ellipse's offset, approximated by a polyline of arcs, since the curve at a distance from an
   *  ellipse is no ellipse.
   *
   *  The ellipse is laid with arcs, two to a stretch, tangent to it at the stretch's ends and to
   *  each other, halved until they lie within the accuracy of it; and that polyline is offset as a
   *  polyline is, each point of whose offset then lies within the accuracy of the ellipse's. Where
   *  the distance is past the radius of curvature, nothing is left round the cusps.
   *
   *  @param  ellipse    The ellipse, whole or an arc of it.
   *  @param  offsetting The distance; the accuracy, no finer than 64 steps of a double at the
   *                     ellipse's size; the gap kind is not read, an ellipse having no corner.
   *  @param  tolerance  Within which two points are one, as for a polyline.
   *  @return            The polylines, each point of them within the accuracy of the true offset;
   *                     closed for a whole ellipse. Spent where more than 65,536 arcs are needed;
   *                     Degenerate for an ellipse of no extent; NotFinite where a number is not
   *                     finite.
   */
  [[nodiscard]] GEOXY_API auto Offset(const Ellipse &ellipse, const Offsetting &offsetting,
                                      const Tolerance &tolerance)
      -> std::expected<std::vector<Polyline>, Error>;

  /**
   *  A spline's offset, approximated by a polyline of arcs, since the curve at a distance from a
   *  spline is no spline.
   *
   *  Each span, a rational Bézier piece, is laid with arcs as an ellipse is, and the polyline they
   *  make is offset as a polyline is. Where two spans meet at an angle, at a knot of full
   *  multiplicity, the gap is closed as Offsetting::Gaps says, as at a polyline's corner.
   *
   *  @param  spline     The spline, closed where its ends are within the tolerance of each other.
   *  @param  offsetting The distance, the gap kind at its corners, and the accuracy, as for an
   *                     ellipse.
   *  @param  tolerance  Within which two points are one, as for a polyline.
   *  @return            The polylines, as for an ellipse. Degenerate for a spline that is not valid
   *                     or has no extent; Spent and NotFinite as for an ellipse.
   */
  [[nodiscard]] GEOXY_API auto Offset(const Spline &spline, const Offsetting &offsetting,
                                      const Tolerance &tolerance)
      -> std::expected<std::vector<Polyline>, Error>;

  /**
   *  Any curve's offset.
   *
   *  A line, a ray, an endless line and an arc give one curve of their own kind, an arc that
   *  collapses none; a polyline, an ellipse and a spline their polylines; a Bézier curve the
   *  polylines of the spline it is.
   *
   *  @param  curve      The curve.
   *  @param  offsetting The distance, the gap kind and the accuracy, as its kind reads them.
   *  @param  tolerance  Within which two points are one, as for a polyline.
   *  @return            The curves; the errors of its kind's offset, but for an arc's collapse.
   */
  [[nodiscard]] GEOXY_API auto Offset(const Curve &curve, const Offsetting &offsetting,
                                      const Tolerance &tolerance)
      -> std::expected<std::vector<Curve>, Error>;
} // namespace Geoxy

#endif /* GEOXY_OFFSET_HPP */
