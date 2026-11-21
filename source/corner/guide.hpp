/***************************************************************************************************
 *  @file      source/corner/guide.hpp
 *  @brief     A curve as a fillet or a chamfer reads it: what it lies on, its parameter there, the
 *             way it runs, where it was picked, and the curve cut or run on to a point.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CORNER_GUIDE_HPP
#define GEOXY_CORNER_GUIDE_HPP

#include "corner/circle.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/fillet.hpp"
#include "geoxy/line.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <cstddef>
#include <cstdint>
#include <expected>
#include <functional>

namespace Geoxy::Corner
{
  /**
   *  What a curve lies on, as a fillet touches it.
   */
  enum class Lies : std::uint8_t
  {
    Straight, /**< A line: a line, a ray, an endless line, or a polyline's straight segment. */
    Round,    /**< A circle: an arc, or a polyline's bulged segment.                         */
    Smooth,   /**< Itself: an ellipse, a spline, a Bézier curve.                             */
  };

  /**
   *  The curve a guide holds until it is given one: a line of no length.
   */
  inline const Curve Unheld = Line{};

  /**
   *  A curve as a fillet reads it.
   *
   *  Its parameter, the along of its places, is its own: a line's, a ray's and an endless line's
   *  from their point in their ways, an arc's or an ellipse's the fraction of its sweep, a
   *  polyline's segment's from 0 at its start to 1 at its end, a spline's and a Bézier curve's
   *  their own. An arc's, a bulged segment's and an ellipse's repeat every turn of the circle.
   */
  struct Guide
  {
    Lies    Kind   = Lies::Straight; /**< What it lies on.                                  */
    Carrier Holder = {};             /**< The line or the circle, held exactly.             */
    std::reference_wrapper<const Curve> Shape = std::cref(Unheld); /**< The curve given. */
    std::size_t Segment = 0;     /**< A polyline's segment picked.                      */
    Point       Origin  = {};    /**< A line's point at 0.                              */
    Vector      Way     = {};    /**< A line's way, rounded: its parameter's unit.      */
    Arc         Bow     = {};    /**< An arc, or a polyline's bulged segment as one.    */
    double      Period  = 0.0;   /**< How far its parameter runs in a turn; 0 for none. */
    double      Pick    = 0.0;   /**< Where it was picked, as its parameter.            */
    bool        Whole   = false; /**< Whether it is kept whole: closed.                 */
  };

  /**
   *  A curve as a fillet reads it.
   *
   *  @param  curve The curve, of finite numbers.
   *  @param  pick  Where it was picked, finite.
   *  @return       The guide; Degenerate for a line of no length, a ray or an endless line of no
   *                way, an arc of no radius or no sweep, a polyline of no segment or picked on a
   *                straight segment of no length, an ellipse flat to a line, or a spline that is
   *                not valid.
   */
  auto Guided(const Curve &curve, Point pick) -> std::expected<Guide, Misfit>;

  /**
   *  The parameter of a point on what a straight or round guide lies on.
   *
   *  @param  guide The guide, straight or round.
   *  @param  point The point.
   *  @return       Its parameter: along the line, or round the circle from its start the way it
   *                runs, below one period.
   */
  auto Along(const Guide &guide, Point point) -> double;

  /**
   *  Of the parameters a point has on a guide, the one nearest another.
   *
   *  @param  guide  The guide.
   *  @param  along  One of the point's parameters.
   *  @param  target The parameter to come near.
   *  @return        along, or along a whole number of periods on, within half a period of target.
   */
  auto Near(const Guide &guide, double along, double target) -> double;

  /**
   *  Whether a place of what a guide lies on is on the curve itself, its ends included.
   *
   *  @param  guide The guide.
   *  @param  along The place's parameter, below one period where it repeats.
   *  @return       Whether it is within the curve's domain, or a polyline's segment picked, to a
   *                few steps of a double; always on a curve kept whole.
   */
  auto Within(const Guide &guide, double along) -> bool;

  /**
   *  The way a guide runs at a place.
   *
   *  @param  guide The guide.
   *  @param  along The place's parameter.
   *  @param  point The place's point.
   *  @return       Its derivative's direction: a line's way, a circle's tangent the way its arc
   *                sweeps, a smooth curve's derivative.
   */
  auto Ahead(const Guide &guide, double along, Point point) -> Vector;

  /**
   *  A guide's curve cut or run on to a point, one side kept.
   *
   *  @param  guide The guide.
   *  @param  place The point, on what it lies on, and its parameter: for an arc or an ellipse, the
   *                one the part kept runs from or to.
   *  @param  keep  1 to keep what lies past it, the way its parameter grows, -1 what lies before.
   *  @return       The curve, of its own kind where it can be: a ray cut keeping its start a line,
   *                an endless line a ray; the curve as it is where it is kept whole; TooLarge where
   *                the point lies past the far end of the part kept, or where an arc or an ellipse
   *                would run on past a turn.
   */
  auto Cut(const Guide &guide, const Place &place, int keep) -> std::expected<Curve, Misfit>;
} // namespace Geoxy::Corner

#endif /* GEOXY_CORNER_GUIDE_HPP */
