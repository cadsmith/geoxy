/***************************************************************************************************
 *  @file      source/met/judge.hpp
 *  @brief     The candidates every pair of parts found, gathered with the ends of parts, each
 *             polished, judged by how the curves pass there, and merged into hits.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_MET_JUDGE_HPP
#define GEOXY_MET_JUDGE_HPP

#include "met/floor.hpp"
#include "met/parts.hpp"
#include "scratch.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <cstddef>
#include <optional>

namespace Geoxy::Met
{
  /**
   *  Two curves being met, and what is worked out once about them.
   */
  struct Scene
  {
    const Curve      *First        = nullptr; /**< The first curve.            */
    const Curve      *Second       = nullptr; /**< The second.                 */
    Geoxy::Tolerance  Tolerance    = {};      /**< The tolerance.              */
    List<Point>       FirstEnds    = {};      /**< The first curve's ends.     */
    List<Point>       SecondEnds   = {};      /**< The second's.               */
    bool              Splined      = false;   /**< Whether either is a spline. */
    List<Floor>       FirstFloors  = {};      /**< The first's parts, bounded. */
    List<Floor>       SecondFloors = {};      /**< The second's.               */
    const List<Part> *FirstParts   = nullptr; /**< The first curve's parts.    */
    const List<Part> *SecondParts  = nullptr; /**< The second's.               */
  };

  /**
   *  A candidate put in the curves' domains and polished by Newton's method on the gap between
   *  them, kept in their domains; a step that does not narrow the gap is not taken, so a tangency
   *  stays where it was found.
   *
   *  @param  first  The first curve.
   *  @param  second The second.
   *  @param  local  The candidate, in the curves' parameters.
   *  @return        The candidate polished.
   */
  auto Neared(const Curve &first, const Curve &second, Local local) -> Local;

  /**
   *  Whether a shared stretch runs on into the next: the same way round, and on both curves the
   *  next entering within it or at its end.
   *
   *  @param  before The stretch before, along the first curve.
   *  @param  after  The stretch after.
   *  @return        Whether they are one.
   */
  auto Continues(const Overlap &before, const Overlap &after) -> bool;

  /**
   *  What every pair of parts found, in the curves' parameters.
   *
   *  @param  parts  The first curve's parts.
   *  @param  others The second's.
   *  @param  scene  The curves being met.
   *  @return        The candidates, the shared stretches in order along the first, and whether a
   *                 budget ran out.
   */
  auto Gathered(const List<Part> &parts, const List<Part> &others, const Scene &scene) -> Found;

  /**
   *  The ends of one curve's parts as candidates, and whether one was left out.
   */
  struct Tips
  {
    List<Local> Locals  = {};    /**< The candidates, in the curves' parameters.  */
    bool        Skipped = false; /**< Whether an end far from the other was left. */
  };

  /**
   *  The ends of one curve's parts, each with its nearest place on the other curve: every one, or
   *  only those its parts' bounds do not put farther than Far tolerances from the other, whose
   *  nearest places, the dearest thing a meeting asks for, are then not sought. A far end's
   *  candidate is polished from where the other curve is nearest it, a long way from any hit:
   *  Newton's method takes it to a hit another candidate finds, or to none. An end of the curve
   *  itself is also met with every other segment of a polyline that passes within the tolerance of
   *  it, at that segment's nearest place.
   *
   *  @param  scene The curves being met.
   *  @param  side  Whose parts' ends: 0 for the first curve's, 1 for the second's.
   *  @param  all   Whether every end, or only those that may lie near the other curve.
   *  @return       The candidates, in the curves' parameters, the first curve's first; and whether
   *                a far end was left out.
   */
  auto Ended(const Scene &scene, std::size_t side, bool all) -> Tips;

  /**
   *  Which side of another curve's tangent a curve passes to: 1 or -1 where it comes in from one
   *  side and goes out to the other, each beyond the tolerance's angle; 0 where it turns back at a
   *  corner or runs along the tangent.
   *
   *  @param  ways    The curve's ways in and out.
   *  @param  other   The other's tangent.
   *  @param  radians The tolerance's angle.
   *  @return         The side, or 0.
   */
  auto Passing(const std::array<Vector, 2> &ways, Vector other, double radians) -> int;

  /**
   *  A hit as judged, with the side the first curve passes to, 0 where it does not cross.
   */
  struct Seen
  {
    Hit    Found  = {};  /**< The hit.                                               */
    int    Side   = 0;   /**< The side the first passes to, or 0.                    */
    double Within = 0.0; /**< The length its two points were held to, as Judged did. */
  };

  /**
   *  A candidate judged: polished, unless it is at an end of a curve already within the tolerance
   *  of the other, which Newton's method would slide off the end along a touch's band; kept only
   *  where both curves are within the tolerance and not in a shared stretch on both, and told how
   *  the curves meet there.
   *
   *  @param  scene     The curves being met.
   *  @param  overlaps  The stretches they share.
   *  @param  candidate The candidate, in the curves' parameters.
   *  @return           The hit, if it holds.
   */
  auto Judged(const Scene &scene, const List<Overlap> &overlaps, Local candidate)
      -> std::optional<Seen>;

  /**
   *  The hits of two curves from their candidates as judged: in order along the first, each put in
   *  the run before it, or in an earlier run it is one with that ends within the tolerance of it,
   *  as two places of a loop that pass through one point interleave, or that is a closed first
   *  curve's first, as its last hit may be its first again; the runs merged.
   *
   *  @param  scene    The curves being met.
   *  @param  overlaps The stretches they share.
   *  @param  seen     The candidates that hold, judged, in the order the candidates came.
   *  @return          The hits.
   */
  auto Hits(const Scene &scene, const List<Overlap> &overlaps, const List<Seen> &seen) -> List<Hit>;

  /**
   *  The hits just past a shared stretch's ends that the curves stay within the tolerance up to,
   *  dropped: a tangency the stretch runs on into, or the end of a curve it already reaches, the
   *  stretch's ends being where the curves' parts end.
   *
   *  @param scene   The curves being met.
   *  @param meeting What meeting them found.
   */
  void Absorb(const Scene &scene, Intersections &meeting);
} // namespace Geoxy::Met

#endif /* GEOXY_MET_JUDGE_HPP */
