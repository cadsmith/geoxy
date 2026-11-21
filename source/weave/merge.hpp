/***************************************************************************************************
 *  @file      source/weave/merge.hpp
 *  @brief     Cuts carried across shared stretches, and vertices within the tolerance merged in a
 *             stated order.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_WEAVE_MERGE_HPP
#define GEOXY_WEAVE_MERGE_HPP

#include "scratch.hpp"
#include "unions.hpp"
#include "weave/found.hpp"

#include "geoxy/domain.hpp"
#include "geoxy/intersect.hpp"

#include <algorithm>
#include <cstddef>

namespace Geoxy::Weave
{
  /**
   *  Where a place on one curve of a shared stretch is on the other: the stretch's parameters run
   *  in step on both, as they do on lines and on arcs of one circle.
   *
   *  @param  common The stretch.
   *  @param  along  The place, on the first curve, or on the second where back.
   *  @param  back   Whether from the second curve to the first.
   *  @return        The place on the other, within its stretch.
   */
  inline auto Reflected(const Common &common, double along, bool back) -> double
  {
    const Range  span = back ? common.Shared.Second : common.Shared.First;
    const Range  onto = back ? common.Shared.First : common.Shared.Second;
    const double part = (along - span.Start) / (span.End - span.Start);
    const double wide = onto.End - onto.Start;
    return std::clamp(common.Shared.Along ? onto.Start + (part * wide) : onto.End - (part * wide),
                      onto.Start, onto.End);
  }

  /**
   *  Every cut inside a shared stretch carried to the other curve, back and forth, until none is
   *  added, as each adds a vertex to a curve that had none there: so the two curves are cut alike
   *  all along it. The cuts carried are kept with the rest; which are carried, and where to, does
   *  not hang on the order a curve's cuts are asked in, as no two of a curve's are at one vertex.
   *
   *  @param found The findings.
   */
  void Carry(Findings &found);

  /**
   *  What merging the vertices made: each vertex's number, and those merged into one at another
   *  point.
   */
  struct Merging
  {
    List<std::size_t> Ids   = {}; /**< Each vertex's number.                             */
    List<std::size_t> Moved = {}; /**< The vertices merged into one elsewhere, in order. */
  };

  /**
   *  Vertices merged in order, each into the first earlier one within a length that stands for
   *  itself, or standing for itself; so no two that stand for themselves are within the length.
   *  Those that stand for themselves are kept in a grid of cells four lengths wide, and each vertex
   *  looks in the cells the square of its length about it reaches, most often one.
   *
   *  Only those Crowds marks look, or are kept: one it leaves has no other within the length, so it
   *  stands for itself and no other would find it. The grid's table, a hash of cells scattered over
   *  memory, is asked of a few vertices, not of every one in the order the meetings made them.
   *
   *  A vertex at the very point of an earlier one, its twin, joins the twin's set without a look:
   *  were the twin to stand for itself, no earlier one standing for itself would lie within the
   *  length of it, and the twin would be the first; were it merged into one, that one would be the
   *  first for both, those standing for themselves after it all lying beyond it.
   *
   *  Each vertex is numbered as it is merged: the next number where it stands for itself, its set's
   *  otherwise, which is the first's of the set, as no merge here joins a set to one standing for
   *  itself after it. Each vertex merged into one at another point is kept.
   *
   *  @param  found  The findings: the vertices and their twins.
   *  @param  length The length, above 0.
   *  @param  sets   The sets the merged ones are joined in.
   *  @return        Each vertex's number among those that stand for themselves, and the vertices
   *                 merged into one at another point.
   */
  auto Merge(const Findings &found, double length, Unions &sets) -> Merging;
} // namespace Geoxy::Weave

#endif /* GEOXY_WEAVE_MERGE_HPP */
