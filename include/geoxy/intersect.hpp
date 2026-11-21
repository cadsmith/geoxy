/***************************************************************************************************
 *  @file      geoxy/intersect.hpp
 *  @brief     Where two curves meet: the places they cross or touch, and the stretches they share.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_INTERSECT_HPP
#define GEOXY_INTERSECT_HPP

#include "geoxy/api.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"

#include <array>
#include <cstdint>
#include <functional>
#include <vector>

namespace Geoxy
{
  /**
   *  How two curves meet at a place.
   */
  enum class Contact : std::uint8_t
  {
    Cross, /**< They cross, each passing from one side of the other to its other side.  */
    Touch, /**< They touch: tangent, within the tolerance, or turning back at a corner. */
    End,   /**< The place is an end of one of them, within the tolerance.               */
  };

  /**
   *  A place where two curves meet.
   */
  struct Hit
  {
    Point   At     = {};             /**< Where, halfway between the two curves' points. */
    Place   First  = {};             /**< The place on the first curve.                  */
    Place   Second = {};             /**< The place on the second curve.                 */
    Contact Kind   = Contact::Cross; /**< How they meet there.                           */
  };

  /**
   *  A stretch two curves share.
   */
  struct Overlap
  {
    Range First  = {};    /**< The stretch of the first curve's parameter, forward.  */
    Range Second = {};    /**< The stretch of the second curve's parameter, forward. */
    bool  Along  = false; /**< Whether the second runs the way the first does.       */
  };

  /**
   *  All that intersecting two curves found.
   *
   *  A hit may be missing, and Spent says so, where the work budget ran out; where a curve is swept
   *  round more than the 1024 turns an intersection looks at, or two such curves pass a point on
   *  more than 2048 pairs of turns; and where a spline's weights lie farther apart than doubles
   *  hold beside each other, some 2^1080 times. Every hit that is given is one.
   *
   *  A spline with a weight far heavier than its neighbours is met through one that stands for it,
   *  cut until its weights are even: its crossings are found at any weight, their points the
   *  doubles nearest the true ones. A hit's place on such a spline is the double nearest the true
   *  place, which beside a weight 2^52 times its neighbours' or more may be the end of a span the
   *  hit's point is a long way from: the place cannot say more, the point does.
   */
  struct Intersections
  {
    std::vector<Hit>     Hits     = {};    /**< The hits, in order along the first curve.    */
    std::vector<Overlap> Overlaps = {};    /**< The stretches shared, along the first curve. */
    bool                 Spent    = false; /**< Whether a hit may be missing: never silent.  */
  };

  /**
   *  Where two curves meet.
   *
   *  One geometry gives one answer, however its curves are spelled: a line or the same line as a
   *  Bézier curve, an arc or the same arc as a rational spline. Two curves share a stretch where
   *  they stay within the tolerance of each other along it and each of its ends is where a part of
   *  either curve ends, a curve's end, a vertex, a knot or a closed curve's seam, or where the
   *  second turns back on itself; a stretch that ends where the curves part instead is the band
   *  about a touch or a shallow crossing, and is one hit, a touch where the curves are tangent
   *  somewhere in it. A touch is where the curves come nearest: for a line, a ray or an endless
   *  line near tangent to an arc, or two arcs, each curve's point there worked out in closed form;
   *  for any other two, where their tangents are parallel and the way between them square to both,
   *  polished in double-double from the curves' own numbers and rounded once. Its point is the same
   *  to the bit whichever curve is first, and swapped, turned a quarter, mirrored or scaled by a
   *  power of two with the plane, save at a cusp, where a curve's slope comes to nothing, and along
   *  curves that bend alike, which come nearest all along. Two hits are one where the curves stay
   *  within twice the tolerance all the way between them along both, so a curve that passes through
   *  a point twice, as a loop does at its double point, meets a curve there twice; one run of hits
   *  is a touch where the first curve passes to both sides of the second, going in and coming back
   *  out. An end of a curve within the tolerance of the other is a hit where it stands.
   *
   *  Each curve is cut into straight runs, elliptic arcs and rational Bézier pieces, and each pair
   *  of parts solved by its own rule: runs together in closed form, a run and an ellipse in the
   *  unit circle's frame, two circles in closed form, an ellipse with an ellipse or a piece by the
   *  ellipse's implicit equation on the other, a run with a piece by the run's, two pieces by
   *  subdivision within a budget and Newton's method. Runs, and pieces straight to within the
   *  tolerance's angle, lie on one line where each end of what they share is within the tolerance
   *  of the other's line; circles and ellipses on one carrier are decided exactly, and two arcs of
   *  circles that part share a stretch from an end to an end where each stays within the tolerance
   *  of the other's circle all the way, in closed form; other pieces lie along each other where
   *  nine places of one have the other within the tolerance.
   *
   *  @param  first     The first curve.
   *  @param  second    The second.
   *  @param  tolerance Within which two points are one: hits closer merge, a curve that comes as
   *                    near another without crossing it touches it, curves as near all along share
   *                    a stretch; and the angle within which two tangents are one, and a piece is
   *                    straight. A length finer than 16 steps of a double at the size of the
   *                    curves' numbers is taken as that: finer, where they lie is not known. A
   *                    length or an angle that is not a number, or below 0, is taken as 0, and one
   *                    past every double as the largest double.
   *  @return           The hits, merged and in order along the first curve, none within a shared
   *                    stretch; the shared stretches, joined where they run on, in order; and
   *                    whether the budget ran out.
   */
  [[nodiscard]] GEOXY_API auto Intersect(const Curve &first, const Curve &second,
                                         const Tolerance &tolerance) -> Intersections;

  /**
   *  Two curves to meet, the first first.
   */
  using Curves = std::array<std::reference_wrapper<const Curve>, 2>;

  /**
   *  Where two curves meet, as Intersect answers, written into intersections kept by the caller:
   *  their hits and stretches replaced, the room they had kept, so that intersecting many pairs
   *  into one asks the system for no memory once it has room enough.
   *
   *  @param curves    The two curves, the first first.
   *  @param tolerance The tolerance, as Intersect takes it.
   *  @param into      The intersections, whatever they held replaced.
   */
  GEOXY_API void Intersect(const Curves &curves, const Tolerance &tolerance, Intersections &into);

  /**
   *  Where a curve crosses or touches itself, as PEDIT and a region's check ask.
   *
   *  The curve is cut into pieces that cannot cross themselves: a polyline into its segments; a
   *  Bézier curve, and a spline span by span, halved until each piece's control points run one way
   *  in x or in y; an arc or an ellipse of more than a turn into turns. Every two pieces are met
   *  through Intersect's exact kernels, the pairs found by a Grid's cells, and where two pieces
   *  follow one another the end they share is no hit. Hits of several pairs at one place, as where
   *  the curve passes through one of its vertices again, are one; where Intersect tells only that
   *  pieces end there, the hit is a crossing where one pass's ways in and out lie either side of
   *  the other's, decided by Orient, a touch where not, and an end at an open curve's own ends.
   *
   *  @param  curve     The curve.
   *  @param  tolerance Within which two places of it are one.
   *  @return           The hits, the first place of each before its second, both in the curve's own
   *                    parameter, in order along it; the stretches it runs along itself, as pairs
   *                    of stretches; and whether a part of the answer may be missing.
   */
  [[nodiscard]] GEOXY_API auto SelfIntersect(const Curve &curve, const Tolerance &tolerance)
      -> Intersections;
} // namespace Geoxy

#endif /* GEOXY_INTERSECT_HPP */
