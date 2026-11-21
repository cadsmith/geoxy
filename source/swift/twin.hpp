/***************************************************************************************************
 *  @file      source/swift/twin.hpp
 *  @brief     Two quadratic or cubic Bézier curves met by a swift kernel where the thorough way
 *             would leave nothing to decide: each crossing the double nearest the true point.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_TWIN_HPP
#define GEOXY_SWIFT_TWIN_HPP

#include "scratch.hpp"
#include "swift/section.hpp"

#include "geoxy/domain.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"

#include <array>
#include <optional>

namespace Geoxy::Swift
{
  /**
   *  A crossing the twin kernel proved, and how near the curves come about it: what a meeting of
   *  pieces weighs it against the crossings of other pieces by.
   */
  struct Twain
  {
    Point                 At    = {};  /**< The point, the double nearest the true one. */
    std::array<double, 2> Along = {};  /**< Its places, each the curve's own parameter and
                                            the double nearest the true one.                */
    std::array<Range, 2> Near = {};    /**< Each piece's places, from 0 to 1, where the
                                            curves may come within the tolerance about it.  */
    std::array<Range, 2> Twice  = {};  /**< Likewise within twice the tolerance.           */
    double               Runs   = 0.0; /**< How far either curve may run on its near band. */
    double               Within = 0.0; /**< How far apart the curves may be where the
                                            thorough way finds them within the tolerance,
                                            and its own working's grain.                    */
  };

  /**
   *  The crossings the twin kernel proved.
   */
  using Twains = List<Twain>;

  /**
   *  Two quadratic or cubic pieces met, as the thorough way meets them as parts of their curves,
   *  where that leaves nothing to decide: wherever the pieces come within the reach of each other,
   *  a crossing there, proved by Krawczyk's operator alone in its box; every crossing steeper than
   *  twice the tolerance's angle, far from either piece's ends and from any other, and every
   *  stretch of the first piece within the tolerance of the second about it parting from it at more
   *  than the angle; each crossing the double nearest the true point, its places the doubles
   *  nearest the true ones. A piece's control points within its slack of the true ones are charged
   *  that slack wherever they are weighed, and a crossing is proved on the true pieces.
   *
   *  @param  sections  The two pieces, the first curve's first.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @return           The crossings, placed, the first curve's place first; nothing where the
   *                    thorough way must decide.
   */
  auto Twined(const std::array<const Section *, 2> &sections, const Tolerance &tolerance)
      -> std::optional<Twains>;
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_TWIN_HPP */
