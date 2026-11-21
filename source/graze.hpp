/***************************************************************************************************
 *  @file      source/graze.hpp
 *  @brief     Two ready curves that cross shallowly, nearly touch or nearly miss, a line or an arc
 *             against an arc, told from what they keep as Intersect tells them, to the bit.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_GRAZE_HPP
#define GEOXY_GRAZE_HPP

#include "abut.hpp"

#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/point.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <functional>
#include <optional>

namespace Geoxy
{
  /**
   *  Two ready curves met where they come within a few tolerances of touching, or cross at a sine
   *  under Steep, as Intersect meets them, from what they keep: the pairs no pass that tells
   *  crossings at a steep angle can tell, which Intersect's direct way tells by its bands, in
   *  closed form, after a long way round.
   *
   *  Every test Intersect's direct way makes on its own rough numbers is made here on the same
   *  numbers, worked out by the same steps in the same order, so it comes out the same; every test
   *  Intersect makes on an exact point, or on an angle, is made here where the point or the angle
   *  lies clear of the test's threshold by more than its roundings can close. The points are worked
   *  out in lean double-double and certified the doubles nearest the true crossings, as Intersect's
   *  exact way rounds them. A pair any test is not clear for is left to Intersect.
   */
  struct Glancing
  {
    /**
     *  A line or an arc as the bands read it: what a ready curve keeps of it, or for an arc the
     *  same worked out by the same steps, with none of the rest a ready curve makes; and an arc's
     *  ends as At gives them, where they are known. A plain arc's ends are by the library's sine
     *  and cosine at its angles in radians, and its square is its radius squared, rounded.
     */
    struct Held
    {
      const Curve         *Shape   = nullptr; /**< The curve.                                 */
      std::array<Point, 2> Ends    = {};      /**< A line's ends, a plain arc's by sine.      */
      std::array<Point, 2> Tips    = {};      /**< An arc's ends as At gives them, if Tipped. */
      Vector               Slip    = {};      /**< What rounding left of a line's way.        */
      double               Sweep   = 0.0;     /**< A plain arc's sweep, in radians.           */
      double               Square  = 0.0;     /**< A line's way squared, or an arc's radius.  */
      double               Under   = 0.0;     /**< What that square leaves: the two are it.   */
      double               Inverse = 0.0;     /**< 1 / a line's square, rounded.              */
      Box                  Bounds  = {};      /**< Its box, as Confines gives it.             */
      double               Reach   = 0.0;     /**< The largest magnitude of its box's edges.  */
      bool                 Plain   = false;   /**< Whether a line or an arc with a box.       */
      bool                 Round   = false;   /**< Whether a plain arc of up to a turn.       */
      bool                 Tipped  = false;   /**< Whether Tips holds an arc's ends.          */
    };

    /**
     *  Two curves as the bands read them, the first first.
     */
    using Holds = std::array<std::reference_wrapper<const Held>, 2>;

    /**
     *  A line or an arc as the bands read it, worked out as a ready curve keeps it, by the same
     *  steps, and an arc's ends as At gives them: none of the rest a ready curve makes, which
     *  meeting it with a line or an arc by the bands needs nothing of.
     *
     *  @pre          The curve is a line or an arc.
     *  @param  curve The line or arc, which must outlive what is made of it.
     *  @param  tips  An arc's start and end, as At gives them; a line's are not asked.
     *  @return       It, held.
     */
    static auto Hold(const Curve &curve, const std::array<Point, 2> &tips) -> Held;

    /**
     *  A line and an arc, either first, or two arcs, as Intersect meets them, as Steeply tells
     *  them, from what the bands read of them: the length Intersect uses for them, nothing where
     *  their boxes lie too far apart, and their crossings, touches or nothing as Steeply tells
     *  them; where a crossing of two arcs lies at an end of each, within a few steps of a double of
     *  the end At gives, as Intersect places it on both without placing it.
     *
     *  @param  curves The two curves, held, the first first.
     *  @param  sane   The tolerance, as Sane makes it.
     *  @return        What they meet in; nothing where Intersect must decide.
     */
    static auto Hinged(const Holds &curves, const Tolerance &sane) -> std::optional<Abutment>;

    /**
     *  Two ready curves where they come near each other, each a line or an arc and one at least an
     *  arc, as Intersect meets them: nothing where they keep more than the tolerance apart; one
     *  touch where they come within it, or cross no deeper than twice it, at the point halfway
     *  between where the two come nearest, as Intersect's bands place it; or the crossings, each
     *  the double nearest its own, in order along the first curve.
     *
     *  @param  curves The two ready curves.
     *  @param  usable The tolerance they are met within: the length as Reached gives it, within
     *                 Coarse and Fine, and the angle as Sane makes it.
     *  @return        What they meet in; nothing for any other pair, or where Intersect must
     *                 decide.
     */
    static auto Met(const Prepared::Pair &curves, const Tolerance &usable)
        -> std::optional<Abutment>;

    /**
     *  A ready line and a ready arc, either first, as Intersect meets them, as Met tells them, and
     *  where they cross at a sine past Steep with an end near the other curve, or a whole circle's
     *  seam near the line: Roundly's crossings, hits at an end where an end lies at one, a hit at a
     *  whole circle's seam first; or, where Roundly gives them up, a line stopping short of the
     *  circle within the tolerance, as Shy tells it, a hit at its end. Asked after the passes that
     *  tell such pairs faster, as most are told.
     *
     *  @param  curves The two ready curves.
     *  @param  usable The tolerance they are met within, as Met takes it.
     *  @return        What they meet in; nothing for any other pair, or where Intersect must
     *                 decide.
     */
    static auto Steeply(const Prepared::Pair &curves, const Tolerance &usable)
        -> std::optional<Abutment>;

  private:
    /**
     *  A ready curve as the bands read it.
     *
     *  @param  curve The ready curve.
     *  @return       It, held.
     */
    static auto Hold(const Prepared &curve) -> Held;

    /**
     *  What is told of two curves near each other, from what the bands read of them.
     */
    struct Bands;
  };
} // namespace Geoxy

#endif /* GEOXY_GRAZE_HPP */
