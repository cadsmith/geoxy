/***************************************************************************************************
 *  @file      geoxy/fillet.hpp
 *  @brief     FILLET: the arc of a radius tangent to two curves, each cut or run on to where the
 *             arc touches it; and every corner of a polyline rounded where the radius fits.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_FILLET_HPP
#define GEOXY_FILLET_HPP

#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/tolerance.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <vector>

namespace Geoxy
{
  /**
   *  Why a fillet or a chamfer of two curves, or of a polyline's corner, has no answer: AutoCAD's
   *  own reasons, as values.
   */
  enum class Misfit : std::uint8_t
  {
    NotFinite,  /**< A number given is not finite. */
    Degenerate, /**< A curve with no way to go by: a line of no length, an arc of no radius, a
                     spline that is not valid; a radius or a distance below nothing; a chamfer
                     of a polyline's corner where an arc meets it.                             */
    Parallel,   /**< Lines that never meet and cannot be joined: on one line, an endless line
                     first, or a chamfer of parallel lines.                                    */
    TooLarge,   /**< The radius or a distance is too large: the arc or the chamfer would touch
                     a curve past the far end of the part of it kept, or no circle of the
                     radius fits the corner between the parts kept.                            */
    Apart,      /**< The curves never meet, and no circle of the radius touches both on the
                     parts kept.                                                               */
  };

  /**
   *  What FILLET asks of two curves: the radius, and where each was picked.
   */
  struct Filleting
  {
    double               Radius = 0.0; /**< The radius; 0 cuts or runs both on to their corner. */
    std::array<Point, 2> Picks  = {};  /**< Where each curve was picked, on the part of it to
                                            keep, as FILLET's picks say which.                */
  };

  /**
   *  Two curves filleted: each cut or run on to where the arc touches it, and the arc.
   */
  struct Filleted
  {
    Curve              First  = {};    /**< The first curve, cut or run on to its touch. */
    Curve              Second = {};    /**< The second, likewise.                        */
    std::optional<Arc> Round  = {};    /**< The arc, from the first's touch to the second's,
                                            either way round; none for a radius of 0.          */
    std::array<Point, 2> Touches = {}; /**< Where the arc touches each, or their corner for a
                                            radius of 0: the double nearest the exact point
                                            where both are lines or arcs.                      */
  };

  /**
   *  The arc of a radius tangent to two curves, and each cut or run on to where it touches it, as
   *  FILLET draws them.
   *
   *  Each curve is taken as what it lies on: a line, a ray and an endless line as their line, an
   *  arc as its circle, an ellipse as its whole ellipse, a polyline as the segment picked, a spline
   *  or a Bézier curve as itself. Where the two meet, the corner nearest the picks says which part
   *  of each is kept, the side of it the pick lies on, and the arc is the circle of the radius
   *  tangent to both in the angle between the parts kept nearest the corner; where they never meet,
   *  the arc is the circle tangent to both, nearest the picks, that leaves each curve smoothly onto
   *  the other. Each curve is then cut there, or run on to it: a line, a ray, an arc and an ellipse
   *  along themselves, a polyline along its picked segment, its segments past the touch on the
   *  corner's side dropped. A circle, a whole ellipse, and a closed polyline or spline are kept
   *  whole, as AutoCAD keeps a circle. Two parallel lines are joined by the half circle between
   *  them at the first's end nearest its pick, the second cut or run on to it, whatever the radius:
   *  AutoCAD's rule.
   *
   *  Between lines and arcs, every point is exact: each touch and the center is the double nearest
   *  the exact point, certified in double-double or found by exact comparisons, ties to even.
   *  Against an ellipse, a spline or a Bézier curve the circle is found by Newton's method on the
   *  curves' offsets, to a few steps of a double at their size.
   *
   *  @param  first     One curve.
   *  @param  second    The other.
   *  @param  filleting The radius and the picks.
   *  @return           The curves and the arc; why not, as a misfit.
   */
  [[nodiscard]] GEOXY_API auto Fillet(const Curve &first, const Curve &second,
                                      const Filleting &filleting)
      -> std::expected<Filleted, Misfit>;

  /**
   *  A corner of a polyline left as it was, and why.
   */
  struct Untouched
  {
    std::size_t Vertex = 0;                /**< The vertex at the corner.  */
    Misfit      Why    = Misfit::TooLarge; /**< Why it was left as it was. */
  };

  /**
   *  A polyline with its corners rounded or cut, and the corners left as they were.
   */
  struct Reshaped
  {
    Polyline               Shape = {}; /**< The polyline.                             */
    std::vector<Untouched> Left  = {}; /**< The corners left, in the vertices' order. */
  };

  /**
   *  A polyline with an arc of a radius at each corner where it fits, as FILLET's Polyline option
   *  rounds one.
   *
   *  Each corner where two segments meet at an angle, lines or arcs, is rounded by the arc of the
   *  radius tangent to both on the inside of the turn, its touches the doubles nearest the exact
   *  ones, as two curves are filleted, inserted as a bulge; the corners are rounded in order, each
   *  within what the one before left of their segment. A corner where the arc does not fit, its
   *  touch past what is left of a segment by more than the tolerance, is left, and reported; a
   *  segment left shorter than the tolerance is dropped. A corner where the two segments run on
   *  smoothly is no corner, and is passed by.
   *
   *  @param  polyline  The polyline.
   *  @param  radius    The radius, above nothing.
   *  @param  tolerance Within which a touch still fits and a segment is of no length.
   *  @return           The polyline and the corners left; NotFinite for a number that is not
   *                    finite, Degenerate for a radius not above nothing.
   */
  [[nodiscard]] GEOXY_API auto Fillet(const Polyline &polyline, double radius,
                                      const Tolerance &tolerance)
      -> std::expected<Reshaped, Misfit>;
} // namespace Geoxy

#endif /* GEOXY_FILLET_HPP */
