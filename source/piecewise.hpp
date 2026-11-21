/***************************************************************************************************
 *  @file      source/piecewise.hpp
 *  @brief     A polyline met segment by segment: each segment a line or an arc of its own, met
 *             directly with each piece of the other curve whose box comes near it, where the
 *             thorough way, every pair of parts weighed at once, would give the same.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_PIECEWISE_HPP
#define GEOXY_PIECEWISE_HPP

#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/tolerance.hpp"

#include <array>
#include <memory>

namespace Geoxy
{
  /**
   *  What a polyline or an arc made ready keeps to be met piece by piece: its legs, each segment's
   *  box, ends and, for a bulge, its circle and its arc, and the grid they are filed in; an arc's
   *  one leg, its box and ends; worked out once as Piecewise would work them out at each meeting.
   */
  struct Legwork;

  /**
   *  A curve's legwork, worked out once: for a polyline, its legs as Piecewise makes them, every
   *  bulge's arc made, and their grid where they are many; for an arc, its one leg, whose box and
   *  ends take the most working out of any curve met with a polyline.
   *
   *  @param  curve The curve.
   *  @return       Its legwork; nothing for any other curve.
   */
  auto Legworked(const Curve &curve) -> std::shared_ptr<const Legwork>;

  /**
   *  Two curves to meet piece by piece, and the legwork kept for each made ready.
   */
  struct Stance
  {
    std::array<const Curve *, 2>   Shapes = {}; /**< The two curves, the first first. */
    std::array<const Legwork *, 2> Legs   = {}; /**< Each one's legwork, where kept.  */
  };

  /**
   *  Two curves met piece by piece, where one or both is a polyline and the other a line, a ray, an
   *  endless line or an arc: each polyline's segment met as a line or an arc of its own, with each
   *  piece of the other curve whose box comes within a few hundred tolerances of its, by the direct
   *  ways; each hit then put on the curves themselves as the thorough way ends: a crossing of two
   *  straight pieces or of one and an arc as the direct way rounded it once from the curves' own
   *  numbers, one on a polyline's arc worked out again from the bulge's exact circle, or polished
   *  in double-double where that finds none.
   *
   *  The thorough way weighs every pair of parts together, and answers so only where nothing ties
   *  one pair's hits to another's: where every hit is a plain crossing at a sine of an eighth or
   *  more, far from every vertex, every end and every other hit; where the curves come near nowhere
   *  else, no pair of straight pieces lies along one line and no arc nearly touches a piece's
   *  carrier, so that no band and no shared stretch forms; and where the curves part between each
   *  hit and the next, as the thorough way looks to tell them apart. Anything less is left to it,
   *  and so is a crossing neither exact way places on the curves.
   *
   *  A polyline's legwork kept from when it was made ready is used as it is; any other's is worked
   *  out at each meeting, the same.
   *
   *  @param  curves    The two curves, and the legwork kept for each.
   *  @param  tolerance The tolerance, no finer than the curves' grain, as the thorough way holds
   *                    them to.
   *  @param  into      Where what they meet in goes, empty on entry, and left so where the thorough
   *                    way must decide.
   *  @return           Whether it was decided so.
   */
  auto Piecewise(const Stance &curves, const Tolerance &tolerance, Intersections &into) -> bool;
} // namespace Geoxy

#endif /* GEOXY_PIECEWISE_HPP */
