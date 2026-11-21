/***************************************************************************************************
 *  @file      source/thorough.hpp
 *  @brief     Intersections two curves the thorough way, every rule applied, and directly, where
 *             nothing is left to decide: the two ways Intersect takes, for tests that hold one to
 *             the other.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_THOROUGH_HPP
#define GEOXY_THOROUGH_HPP

#include "crossing.hpp"
#include "piecewise.hpp"

#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/tolerance.hpp"

#include <array>
#include <optional>
#include <span>

namespace Geoxy
{
  /**
   *  Where two curves meet, found the thorough way, as Intersect finds it where no direct way
   *  applies.
   *
   *  @param curves    The two curves, the first first.
   *  @param tolerance The tolerance.
   *  @param into      Where what Intersect would find goes, whatever it held replaced.
   */
  void Thorough(const Curves &curves, const Tolerance &tolerance, Intersections &into);

  /**
   *  Two curves with remote arcs met as Intersect goes on with curves within reach of each other,
   *  which leaves them to the thorough way, each such arc a piece: what Slightly asks where the
   *  arcs' chords do not tell what the arcs meet in.
   *
   *  @param curves The two curves.
   *  @param sane   The tolerance, as Sane makes it.
   *  @param into   Where what they meet in goes, whatever it held replaced.
   */
  void Deeply(const Curves &curves, const Tolerance &sane, Intersections &into);

  /**
   *  Two curves within reach of each other, as Intersect goes on with them once their boxes are
   *  weighed.
   */
  struct Entering
  {
    Stance Pair;         /**< The two curves, and the legwork kept for each. */
    double Length = 0.0; /**< The length they are met within, as Reachable or a ready curve's
                              boxes give it.                                  */
    const Prepared::Pair *Kept = nullptr; /**< The curves made ready, where they are: what they
                                            keep tells whether Shrunk leaves them as they
                                            are.                                          */
  };

  /**
   *  Two curves within reach of each other met as Intersect meets them: through others that stand
   *  for them where they are too large, too small or too stiff, directly where nothing is left to
   *  decide, the thorough way otherwise. Both Intersect and Intersect for curves made ready go on
   *  so.
   *
   *  @param entering The two curves, the legwork kept for each, and the length.
   *  @param sane     The tolerance, as Sane makes it.
   *  @param into     Where what they meet in goes, empty on entry.
   */
  void Entered(const Entering &entering, const Tolerance &sane, Intersections &into);

  /**
   *  Where two curves meet, found directly, where their kinds and places leave nothing to decide.
   *
   *  @param  curves    The two curves, the first first.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @param  into      Where what they meet in goes, empty on entry, and left so where the thorough
   *                    way must decide.
   *  @return           Whether it was decided directly.
   */
  auto Direct(const Curves &curves, const Tolerance &tolerance, Intersections &into) -> bool;

  /**
   *  A crossing found roughly, as the exact crossing takes it.
   */
  using Crossings::Approach;

  /**
   *  Crossings found approximately made hits by the exact crossing the thorough way ends in, in
   *  order along the first curve; nothing if any is not found so.
   *
   *  @param  curves The two curves.
   *  @param  near   Each crossing, roughly.
   *  @param  into   Where the hits go, empty on entry, and left so where one is not found.
   *  @return        Whether every exact crossing was found.
   */
  auto Exactly(const std::array<const Curve *, 2> &curves, std::span<const Approach> near,
               Intersections &into) -> bool;

  /**
   *  Two curves near tangent met in closed form, where a line, a ray or an endless line and an arc,
   *  or two arcs, come within four tolerances of touching with every end far from where they come
   *  near: nothing farther apart than the tolerance; one touch, halfway between each curve's point
   *  where they come nearest, where they stay within it all the way, apart or crossing so
   *  shallowly; the two exact crossings where they cross deeper. The thorough way finds the touch
   *  somewhere in the band they stay within the tolerance over, by sampling it; here it is where
   *  they come nearest.
   *
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @param  into      Where what they meet in goes, empty on entry, and left so where it was not
   *                    decided.
   *  @return           Whether it was decided so.
   */
  auto Banded(const Curves &curves, const Tolerance &tolerance, Intersections &into) -> bool;

  /**
   *  The tolerance a meeting uses: the one asked for, no finer than the curves' numbers can hold.
   *
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance asked for.
   *  @return           The tolerance used.
   */
  auto Usable(const Curves &curves, const Tolerance &tolerance) -> Tolerance;

  /**
   *  The box a curve cannot leave, as Reachable reads it: a line's between its ends, an arc's about
   *  its whole circle, any other bounded curve's about what holds it in.
   *
   *  @param  curve The curve.
   *  @return       The box; nothing for a curve that runs on without end, nor where a number of it
   *                is not finite.
   */
  auto Confines(const Curve &curve) -> std::optional<Box>;

  /**
   *  The size of a curve's numbers, as its grain is taken from: the largest magnitude of the points
   *  it is given by, its center and reach for an arc or an ellipse, and for a polyline each bulge's
   *  circle's too, finite ones only.
   *
   *  @param  curve The curve.
   *  @return       The size, 0 for none.
   */
  auto Sized(const Curve &curve) -> double;

  /**
   *  Whether two curves of one part each share a stretch as the thorough way's pair of their parts
   *  finds one: a stretch it keeps, or one it weighs as a band about a touch or a shallow crossing,
   *  which no hit of the direct ways tells.
   *
   *  @pre              Each curve is a line, a ray, an endless line or an arc.
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance the thorough way holds them to.
   *  @return           Whether they do.
   */
  auto Shares(const Curves &curves, const Tolerance &tolerance) -> bool;

  /**
   *  The length two curves are met within, where they may meet at all, as Intersect asks before
   *  anything else: nothing where neither runs on without end and their boxes lie more than four
   *  times that length apart, sixteen where either is not a line or an arc, so that they meet in
   *  nothing.
   *
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance asked for.
   *  @return           The length used, as Usable gives it; nothing where they lie too far apart.
   */
  auto Reachable(const Curves &curves, const Tolerance &tolerance) -> std::optional<double>;

  /**
   *  The lengths past which two curves may be too large, or too small, to meet as they are: where
   *  Reachable gives more than the one or less than the other, Shrunk is asked. Curves of numbers
   *  past 2^98 have a grain past the first, and curves of numbers under 2^-80 met within a
   *  tolerance fit for them a length under the second; nothing a drawing holds comes near either.
   */
  inline constexpr double Coarse = 0x1p50;
  inline constexpr double Fine   = 0x1p-110;

  /**
   *  Two curves whose numbers are so large or so small that their products run past every double or
   *  come to nothing, met at a size where they do not: both scaled by one power of two, which
   *  changes no bit of any number but its exponent, met, and each hit's points scaled back. A place
   *  along a curve is the same at either size. Past 2^100 or under 2^-100, measured on every pair
   *  of kinds, meetings came out wrong; here they come out as at a size near 1, to the bit.
   *
   *  Save where a number falls below the doubles, scaled down or back: a coordinate a sliver of the
   *  rest's size, or a point below the normal doubles. A crossing of two straight curves is then
   *  worked out again on the curves as given, the double nearest the true one, and where a number
   *  of the curves was lost each curve's point is taken on it as given, at its place.
   *
   *  A ray's or an endless line's way is of any length its maker chose, and one so long or so short
   *  that its square is past every double or nothing is scaled likewise, on its own, and the places
   *  along it scaled back.
   *
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance asked for.
   *  @param  into      The meeting, filled where they are so large or a way is so odd.
   *  @return           Whether they are, and were met.
   */
  auto Shrunk(const Curves &curves, const Tolerance &tolerance, Intersections &into) -> bool;

  /**
   *  Whether a size Shrunk weighs is one it meets at another: two curves' size, or an endless
   *  curve's way's largest coordinate, past 2^96 or under 2^-96, nothing and what is not finite
   *  left as they are.
   *
   *  @param  size The size.
   *  @return      Whether it is.
   */
  auto Outscaled(double size) -> bool;

  /**
   *  Whether two ready curves are of a size and of ways Shrunk leaves as they are, as what they
   *  keep tells it: their size not past 2^96 nor under 2^-96, and each ray or endless line one
   *  Prepared keeps a way for, whose largest coordinate lies within them.
   *
   *  @param  kept The two ready curves.
   *  @return      Whether they are; false where it is not known.
   */
  auto Shrinkless(const Prepared::Pair &kept) -> bool;

  /**
   *  Two curves of which one is a spline too stiff to be told by its own places, met through the
   *  splines that stand for such ones, as Eased makes them: every place of the curve is told there
   *  as finely as its coordinates. Each crossing found is polished again from the spline's own
   *  numbers, so that its point is the double nearest the true one, and each place is given back as
   *  the spline's, the double nearest it: beside a weight 2^52 times its neighbours' or more that
   *  may be the end of a span the crossing is a long way from.
   *
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance asked for.
   *  @param  into      The meeting, filled where either is so stiff.
   *  @return           Whether one is, and they were met.
   */
  auto Loosened(const Curves &curves, const Tolerance &tolerance, Intersections &into) -> bool;
} // namespace Geoxy

#endif /* GEOXY_THOROUGH_HPP */
