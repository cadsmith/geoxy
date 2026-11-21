/***************************************************************************************************
 *  @file      source/weave/find.hpp
 *  @brief     Every crossing and shared stretch of an arrangement's curves, found from the pairs
 *             that come near.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_WEAVE_FIND_HPP
#define GEOXY_WEAVE_FIND_HPP

#include "weave/found.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/tolerance.hpp"

#include <cstddef>
#include <span>

namespace Geoxy::Weave
{
  /**
   *  How far Find reaches: the length the curves are met within, at its most for any pair of them;
   *  the length within which vertices merge; and whether runs of lines may be one edge.
   */
  struct Compass
  {
    double Reach   = 0.0;   /**< The length the curves are met within.   */
    double Length  = 0.0;   /**< The length within which vertices merge. */
    bool   Chained = false; /**< Whether runs of lines may be one edge.  */
  };

  /**
   *  Every crossing and shared stretch of the curves: each curve's ends its first two vertices, an
   *  arc past a half turn cut in half, so that no edge is a loop of its own; then what each pair
   *  whose boxes come near meets in, the curves made ready once, the pairs met as IntersectAll
   *  meets a drawing of them: the same pairs that meet, in the same order, each the same to the
   *  bit, with no drawing made. Two lines joined end to end that Jointed tells meet in their joint
   *  alone, which cuts nothing, and are dropped as they are found, before the rest are put in
   *  order. The vertices are laid once the pairs are known, in the curves' order, but where chained
   *  a lone joint's, as Lonely tells it.
   *
   *  @param  sources   The lines and arcs.
   *  @param  tolerance The tolerance.
   *  @param  compass   The length the curves are met within, at its most for any pair of them, the
   *                    length within which vertices merge, and whether chained.
   *  @return           The findings.
   */
  auto Find(std::span<const Curve> sources, const Tolerance &tolerance, const Compass &compass)
      -> Findings;
} // namespace Geoxy::Weave

#endif /* GEOXY_WEAVE_FIND_HPP */
