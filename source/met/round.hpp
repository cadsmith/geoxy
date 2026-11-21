/***************************************************************************************************
 *  @file      source/met/round.hpp
 *  @brief     A straight curve and an arc, or two arcs, met directly: their crossings found
 *             roughly, and arcs of one center or one circle.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_MET_ROUND_HPP
#define GEOXY_MET_ROUND_HPP

#include "met/linear.hpp"
#include "thorough.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <span>

namespace Geoxy::Met
{
  /**
   *  An arc held as points for the direct ways, its angles turned to radians and its cuts worked
   *  out once: what tells whether a point lies on its sweep by turns, no angle taken. The cuts are
   *  worked by the library's sine and cosine, within a few steps of a double of the exact ends: the
   *  direct ways ask of them only whether they lie four tolerances, sixty-four steps at the least,
   *  from the other curve, and which side of them points that far off lie.
   */
  struct Round
  {
    Point                Center = {};  /**< Its center.                         */
    double               Radius = 0.0; /**< Its radius.                         */
    double               Sweep  = 0.0; /**< Its sweep in radians.               */
    std::array<Point, 2> Cuts   = {};  /**< Its start and end, or a whole seam. */
  };

  /**
   *  A plain arc held as points, its cuts worked out.
   *
   *  @param  arc   The arc, plain.
   *  @param  sweep Its sweep in radians, as Plain gives it.
   *  @return       It held as points.
   */
  auto Rounded(const Arc &arc, double sweep) -> Round;

  /**
   *  An arc as the direct ways take it: plain; nothing for any other.
   *
   *  @param  arc The arc.
   *  @return     It held as points.
   */
  auto Rounded(const Arc &arc) -> std::optional<Round>;

  /**
   *  How far a point is from a circle: from its center, less its radius.
   *
   *  @param  point The point.
   *  @param  round The arc whose circle.
   *  @return       The distance, no more than to any part of the arc.
   */
  auto Off(Point point, const Round &round) -> double;

  /**
   *  Whether a point on an arc's circle lies within its sweep, by the turns from its cuts seen from
   *  the center, no angle taken: within a sweep up to a half turn where it lies left of the start
   *  and right of the end, the sweep taken counterclockwise; within a longer one where it lies
   *  outside the rest of the turn so; anywhere on a whole turn. Asked only of points far from the
   *  cuts, where no turn is near zero.
   *
   *  @param  round The arc.
   *  @param  point The point.
   *  @return       Whether it does.
   */
  auto Swept(const Round &round, Point point) -> bool;

  /**
   *  Hits put in the thorough way's order along a whole circle first: those within the tolerance of
   *  its seam, found there at its first pass and polished to a hair short of its last, before the
   *  rest, each kept in the order it came. Each is rotated in behind those at the seam before it:
   *  the order a stable partition gives, in place, with no room taken.
   *
   *  @param circle    The circle, first of the two curves.
   *  @param tolerance The tolerance.
   *  @param into      The meeting, its hits in order along the circle.
   */
  void Seamwise(const Arc &circle, const Tolerance &tolerance, Intersections &into);

  /**
   *  The crossings a direct way finds, roughly: no more than two, held in place, so that finding
   *  them asks for no memory.
   */
  struct Nearby
  {
    std::array<Approach, 2> Items = {}; /**< The crossings.                          */
    std::size_t             Count = 0;  /**< How many there are.                     */
    std::array<Point, 4>    Ends  = {}; /**< The ends of either curve at a crossing. */
    std::size_t             Ended = 0;  /**< How many there are.                     */
    std::array<Hit, 2>      Ready = {}; /**< Hits found whole: a band's touch, or the
                                          ends arcs of one circle share.          */
    std::size_t Made = 0;               /**< How many there are. */

    /**
     *  One more crossing.
     *
     *  @param approach The crossing.
     */
    void Add(const Approach &approach)
    {
      Items.at(Count) = approach;
      ++Count;
    }

    /**
     *  One more end at a crossing.
     *
     *  @param end The end.
     */
    void End(Point end)
    {
      Ends.at(Ended) = end;
      ++Ended;
    }

    /**
     *  The crossings.
     *
     *  @return They.
     */
    [[nodiscard]] auto Held() const -> std::span<const Approach>
    {
      return std::span{Items}.first(Count);
    }

    /**
     *  Whether an end of either curve lies within a length of a point.
     *
     *  @param  point  The point.
     *  @param  length The length.
     *  @return        Whether one does.
     */
    [[nodiscard]] auto Ending(Point point, double length) const -> bool
    {
      return std::ranges::any_of(std::span{Ends}.first(Ended), [point, length](Point end) -> bool
                                 { return Distance(end, point) <= length; });
    }
  };

  /**
   *  Where the carriers of two curves cross, roughly: up to two points, each with its place on the
   *  first curve where that is a straight one.
   */
  struct Carried
  {
    std::array<Point, 2>  At    = {}; /**< The points.                       */
    std::array<double, 2> Along = {}; /**< Their places on a straight first. */
    std::size_t           Count = 0;  /**< How many there are.               */
  };

  /**
   *  A straight curve and an arc met directly, where nothing but crossings at a steep angle can
   *  come of them: the carrier missing the circle by more than four tolerances, or cutting into it
   *  deeper than that, at a sine above Steep; each end of either farther than four tolerances from
   *  the other's carrier or within half the tolerance of a crossing, which is then a hit at an end;
   *  a whole circle's seam far. A crossing is a hit where it lies on both, within the straight's
   *  stretch or at its end, within the arc's sweep or at its end. Two crossings lie a quarter of
   *  the radius apart at the least, more than 16 tolerances for a radius above 64.
   *
   *  @param  straight  The straight curve.
   *  @param  arc       The arc.
   *  @param  tolerance The tolerance, no finer than their grain.
   *  @return           Each crossing's places, the straight one first, and the ends at them;
   *                    nothing where the thorough way must decide.
   */
  auto Roundly(const Straight &straight, const Arc &arc, const Tolerance &tolerance)
      -> std::optional<Nearby>;

  /**
   *  An arc's ends: none for a whole circle, whose seam is no end, the place a band about it lies
   *  across and is met in closed form all the same.
   *
   *  @param  round The arc, held as points.
   *  @return       Its ends.
   */
  auto Cuts(const Round &round) -> std::span<const Point>;

  /**
   *  The ends two arcs share: pairs of ends, one of each, within half the tolerance of each other,
   *  each end in one pair at most, so two pairs at most.
   */
  struct Abutted
  {
    std::array<std::array<std::size_t, 2>, 2> Pairs = {}; /**< Which end of each.  */
    std::size_t                               Count = 0;  /**< How many there are. */
  };

  /**
   *  The ends two arcs share, a whole circle's seam not an end: each pair of ends, one of each,
   *  within half the tolerance of each other or farther than four tolerances apart; nothing where a
   *  pair lies between, nor where an end lies within half the tolerance of both ends of the other
   *  arc. That arc's ends lie that near each other, no longer than the tolerance or all but a whole
   *  circle: as many as four pairs do, their hits one or along a stretch, the thorough way's to
   *  tell.
   *
   *  @param  rounds    The two arcs, held as points.
   *  @param  tolerance The tolerance.
   *  @return           The pairs within half the tolerance.
   */
  auto Abutting(const std::array<const Round *, 2> &rounds, const Tolerance &tolerance)
      -> std::optional<Abutted>;

  /**
   *  Two arcs of one circle, its center and radius the same to the bit, not both whole, a whole
   *  one's seam off the other or at its end, met directly where they share stretches and nothing
   *  else: every end of each farther than four tolerances from the other's ends, so either within
   *  the other's sweep, and an end of a stretch they share, or outside it and far from it; so the
   *  two meet at no point alone. There the thorough way's answer is the stretches the sweeps share,
   *  worked out by its own steps and joined.
   *
   *  @param  arcs      The two arcs.
   *  @param  tolerance The tolerance, no finer than their grain.
   *  @param  into      Where what they meet in goes, empty on entry.
   *  @return           Whether it was decided so.
   */
  auto Concentrically(const std::array<const Arc *, 2> &arcs, const Tolerance &tolerance,
                      Intersections &into) -> bool;

  /**
   *  Two arcs of one circle, or of two within a quarter of the tolerance of each other, centers and
   *  radii together, met as Intersect's direct way meets them: by Concentrically where they are of
   *  one circle to the bit; where it leaves them, their centers so near that no rough crossing is
   *  found, as Banded meets arcs of one circle, at ends they share alone or along stretches. Its
   *  answer is the direct way's, by the same steps.
   *
   *  @param  arcs      The two arcs.
   *  @param  tolerance The tolerance, no finer than their grain.
   *  @param  into      Where what they meet in goes, empty on entry.
   *  @return           Whether it was decided so; never for circles farther apart.
   */
  auto Coincident(const std::array<const Arc *, 2> &arcs, const Tolerance &tolerance,
                  Intersections &into) -> bool;

  /**
   *  Two arcs met directly, where nothing but crossings at a steep angle can come of them: their
   *  circles apart, crossing or one inside the other by more than four tolerances, crossing at a
   *  sine above Steep; each end of either farther than four tolerances from the other's circle or
   *  within half the tolerance of a crossing, which is then a hit at an end; a whole circle's seam
   *  far. A crossing is a hit where it lies on both, within each sweep or at its end. Two crossings
   *  lie an eighth of the smaller radius apart at the least, by the law of sines, more than eight
   *  tolerances for radii above 64.
   *
   *  @param  first     One arc.
   *  @param  second    The other.
   *  @param  tolerance The tolerance, no finer than their grain.
   *  @return           Each crossing's places, and the ends at them; nothing where the thorough way
   *                    must decide.
   */
  auto Circularly(const Arc &first, const Arc &second, const Tolerance &tolerance)
      -> std::optional<Nearby>;
} // namespace Geoxy::Met

#endif /* GEOXY_MET_ROUND_HPP */
