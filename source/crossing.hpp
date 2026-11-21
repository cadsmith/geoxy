/***************************************************************************************************
 *  @file      source/crossing.hpp
 *  @brief     Where two lines or conics cross, worked out exactly from the curves' own numbers and
 *             rounded once: the double nearest the true point.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CROSSING_HPP
#define GEOXY_CROSSING_HPP

#include "eased.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/line.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <optional>

namespace Geoxy::Crossings
{
  /**
   *  A crossing worked out exactly: the double nearest the point, and its place on each curve.
   */
  struct Crossed
  {
    Point                 At    = {}; /**< The point, rounded once.                */
    std::array<double, 2> Along = {}; /**< Its place on each curve, nearest.       */
    std::array<Point, 2>  On    = {}; /**< Each curve's point at its place, as At. */
  };

  /**
   *  A crossing found roughly: its places on the two curves as far as they matter, a polyline's
   *  segment, and its point, which the exact crossing starts from and chooses a root by.
   */
  struct Approach
  {
    std::array<double, 2> Along = {};    /**< Its places, the first curve's first.            */
    Point                 At    = {};    /**< Its point, roughly.                             */
    bool                  Touch = false; /**< Whether only where the carriers touch, exactly. */
  };

  /**
   *  Where two curves' carriers cross, found once in double-double for every hit near them: a
   *  crossing of two lines, or the two a square root's signs give for a circle, each the double
   *  nearest it where its bound settles it.
   */
  struct Roots
  {
    std::array<std::optional<Point>, 2> At    = {};    /**< Each crossing, where certain. */
    bool                                Twice = false; /**< Whether there are two.        */
    bool                                Apart = false; /**< Whether the carriers certainly
                                                            do not touch: a square root's
                                                            radicand clear of zero.         */
  };

  /**
   *  Where two curves' carriers cross, in double-double with error bounds, where both are lines or
   *  arcs: a hit near them takes the crossing it is at from these, with no more working.
   *
   *  @param  curves The two curves.
   *  @return        Their crossings, each where certain; none where they are not lines or arcs, or
   *                 a root does not settle.
   */
  auto Quickly(const std::array<const Curve *, 2> &curves) -> Roots;

  /**
   *  One of the two crossings of a line's carrier with an arc's circle, in double-double with error
   *  bounds: the same point Quickly gives for it among its two, the other left unworked.
   *
   *  @param  line    The line.
   *  @param  arc     The arc whose circle.
   *  @param  further Whether the crossing further along the line, of the two, or the nearer.
   *  @return         The double nearest the crossing, where certain; nothing where its bound does
   *                  not settle it, or the carrier does not cut clearly into the circle.
   */
  auto Quickly(const Line &line, const Arc &arc, bool further) -> std::optional<Point>;

  /**
   *  Where on a curve a point on its carrier lies, if on the curve: within a line's ends, past a
   *  ray's start, anywhere on an endless line or a whole conic; within an arc's span, or a few
   *  steps of a double from its end. The point is the true one rounded, so a few steps either way
   *  are its. A polyline's place is taken on the segment the hit's place names, an arc's and an
   *  ellipse's by its angle; any other curve's is its nearest, sought once.
   *
   *  @param  curve The curve.
   *  @param  along The hit's place on it, to choose a polyline's segment or a turn.
   *  @param  point The point.
   *  @return       The place; nothing where the point is not on the curve.
   */
  auto Located(const Curve &curve, double along, Point point) -> std::optional<Place>;

  /**
   *  A crossing's point placed on two curves, as Crossing places the point it works out: on the
   *  line, the ray, the endless line or the polyline's segment the hit's place names by the foot,
   *  on an arc or an ellipse by its angle, and on any other curve by its nearest place; each kept
   *  only where the point lies on its curve, to within a few steps of a double.
   *
   *  @param  curves The two curves.
   *  @param  places The hit's place on each, roughly, to choose a polyline's segment or a turn.
   *  @param  point  The point, the double nearest the crossing.
   *  @return        The crossing; nothing where the point is not on both.
   */
  auto Placed(const std::array<const Curve *, 2> &curves, const std::array<double, 2> &places,
              Point point) -> std::optional<Crossed>;

  /**
   *  Where two curves cross near a hit, exactly, as Crossing does, its crossing taken from roots
   *  found already where they hold it; for an approach that asks for a touch, only where the
   *  carriers touch there exactly, their quadratic's two roots one, and at once nothing where the
   *  roots show they do not.
   *
   *  @param  curves   The two curves.
   *  @param  approach The hit's place on each, and where it is.
   *  @param  roots    The curves' crossings found in double-double, as Quickly finds them.
   *  @return          The crossing; nothing where it has none.
   */
  auto Crossing(const std::array<const Curve *, 2> &curves, const Approach &approach,
                const Roots &roots) -> std::optional<Crossed>;

  /**
   *  A crossing worked out exactly, as Crossing works it out, with no point on the curves: its
   *  point, its place on each, and on an arc placed by its angle the direction there, which At
   *  turns its radius to and Derivative turns a quarter.
   */
  struct Bearings
  {
    Point                                At    = {}; /**< The point, rounded once.          */
    std::array<double, 2>                Along = {}; /**< Its place on each curve, nearest. */
    std::array<std::optional<Vector>, 2> Ways  = {}; /**< On an arc placed by its angle, the
                                                          direction there.                  */
  };

  /**
   *  Where two curves cross near a hit, exactly, as Crossing with roots found already works it out
   *  and places it, each curve's point there left out: on an arc placed by its angle, the direction
   *  At would turn to there kept instead, so that its tangent costs no second turn.
   *
   *  @param  curves   The two curves.
   *  @param  approach The hit's place on each, and where it is.
   *  @param  roots    The curves' crossings found in double-double, as Quickly finds them.
   *  @return          The crossing, as Crossing's but its points; nothing where Crossing gives
   *                   none.
   */
  auto Aiming(const std::array<const Curve *, 2> &curves, const Approach &approach,
              const Roots &roots) -> std::optional<Bearings>;

  /**
   *  Where two curves cross near a hit, exactly: where the part of each at the hit lies on a line
   *  or an ellipse, as a line, a ray, an endless line, a polyline's segment, an arc or an ellipse
   *  does, the point their carriers cross at nearest the hit, worked out in expansions from the
   *  curves' own numbers and rounded once to the nearest double; kept only where it lies on both
   *  curves, to within a few steps of a double at the size of their numbers. A line and a line by
   *  their determinant, a line and a conic by its quadratic, a line and a bulge by its exact
   *  rational circle's, two circles by their radical line; nothing for two ellipses, a bulge and
   *  any curve but a line, or a Bézier piece.
   *
   *  @param  curves The two curves.
   *  @param  places The hit's place on each.
   *  @param  near   Where the hit is.
   *  @return        The crossing; nothing where it has none.
   */
  auto Crossing(const std::array<const Curve *, 2> &curves, const std::array<double, 2> &places,
                Point near) -> std::optional<Crossed>;

  /**
   *  Where two polylines' bulges cross near a hit: of their exact circles' two crossings the one
   *  nearer the hit, the double nearest it. One Newton step from a guess settles most, as Hopped
   *  takes it; else their radical line met with the smaller circle in double-double with error
   *  bounds, as Crossbows works it out, the point settled where its bounds settle it and the two
   *  crossings lie clearly apart in distance from the hit; else the same worked out in expansions
   *  and rounded once, the nearer of the two rounded points taken, the root added where they are as
   *  near, a coordinate of nothing +0.
   *
   *  @pre              Each polyline's segment the hit's place names is a bulge.
   *  @param  polylines The two polylines.
   *  @param  places    The hit's place on each, to choose a segment.
   *  @param  near      Where the hit is.
   *  @return           The point; nothing where the circles do not cross, or are of one center.
   */
  auto Bowed(const std::array<const Polyline *, 2> &polylines, const std::array<double, 2> &places,
             Point near) -> std::optional<Point>;

  /**
   *  Where two curves' carriers cross near a hit, exactly, as Crossing works it out, not placed on
   *  the curves: so at a size where no place along them can be told, as along a line too short to
   *  square. Lifted, the crossing of curves scaled up by a power of two from their own size, taken
   *  back down: worked out from them in expansions alone and rounded once at their own size, so
   *  that curves scaled up from below the normal doubles cross where they cross as given, to the
   *  bit.
   *
   *  @param  curves   The two curves, their numbers below 2 where lifted.
   *  @param  approach The hit's place on each curve given, and where it is among them.
   *  @param  lift     The power of two the curves were scaled up by: 0 for the curves as given.
   *  @return          The double nearest the crossing; nothing where it has none.
   */
  auto Spotted(const std::array<const Curve *, 2> &curves, const Approach &approach, int lift)
      -> std::optional<Point>;

  /**
   *  Where two curves cross near a hit, polished in double-double: each curve worked out to about
   *  106 bits, as a polynomial or a ratio of two by de Casteljau's or de Boor's algorithm where it
   *  is one (a line, a ray, an endless line, a polyline's straight segment, a Bézier curve, a
   *  spline), else as a conic's quadratic (an arc, an ellipse, a polyline's bulge by its exact
   *  rational circle); Newton's method run on the gap between them until its steps are lost in that
   *  width, on both places, on the polynomial's place alone, or on the point for two conics; the
   *  point rounded once. Nothing where the curves run within 2^-48 radians of parallel, or where
   *  Newton's method does not settle, or settles far off or past an end.
   *
   *  @param  curves The two curves.
   *  @param  places The hit's place on each.
   *  @param  near   Where the hit is.
   *  @return        The crossing; nothing where it is not found so.
   */
  auto Polished(const std::array<const Curve *, 2> &curves, const std::array<double, 2> &places,
                Point near) -> std::optional<Crossed>;

  /**
   *  A curve as a crossing is polished on it: the curve itself, or the curve through a spline that
   *  stands for it, the stand-in's places the ones given and told, and the curve's own points the
   *  ones the crossing is worked out from.
   */
  struct Read
  {
    const Curve  *Shape = nullptr; /**< The curve the crossing is on.                        */
    const Curve  *Stand = nullptr; /**< The curve whose places are given: a stand-in, or it. */
    const Easing *Way   = nullptr; /**< How a place of the stand-in is one of the curve, or
                                        nothing where the curve stands for itself.             */
  };

  /**
   *  Where two curves cross near a hit, polished in double-double as Polished does, each read as
   *  given: a spline too stiff to be told by its own places is told by its stand-in's, its points
   *  worked out from its own numbers at the place the stand-in's stands for, in two doubles, and
   *  its slope taken from the stand-in.
   *
   *  @param  curves The two curves, as they are read.
   *  @param  places The hit's place on each, on its stand-in where it has one.
   *  @param  near   Where the hit is.
   *  @return        The crossing, its places the stand-ins'; nothing where it is not found so.
   */
  auto Retold(const std::array<Read, 2> &curves, const std::array<double, 2> &places, Point near)
      -> std::optional<Crossed>;
} // namespace Geoxy::Crossings

#endif /* GEOXY_CROSSING_HPP */
