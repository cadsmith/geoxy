/***************************************************************************************************
 *  @file      source/belt/join.hpp
 *  @brief     The stretches kept joined end to start into the offset's polylines.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BELT_JOIN_HPP
#define GEOXY_BELT_JOIN_HPP

#include "belt/work.hpp"

#include "geoxy/polyline.hpp"

#include <cstddef>
#include <vector>

namespace Geoxy::Belt
{
  /**
   *  The stretches joined into loops: each chain from a stretch not yet in one, as a polyline run
   *  with the region on its left, which inward is the other way round from the band.
   *
   *  @param  work   The curves and the stretches, in order of their starts.
   *  @param  points The points' numbers.
   *  @return        The loops, about the polyline's first vertex.
   */
  auto Looped(const Work &work, const Nodes &points) -> std::vector<Polyline>;

  /**
   *  The stretches on one side of an open shape, each once: those of the curves that lie on that
   *  side by their making. Where several run from one point to another, the shape running over that
   *  stretch more than once, they are one if every one of them is on the side, and none if any is
   *  on the other or on neither: the shape's nearest points then disagree.
   *
   *  @param work The stretches kept, which are sorted and thinned.
   */
  void Flanking(Work &work);

  /**
   *  One side of an open shape as polylines: the stretches kept on that side joined end to start, a
   *  chain starting wherever no stretch arrives, and what is left then in loops; each run the way
   *  the shape runs.
   *
   *  @param  work   The curves and the stretches kept on the side.
   *  @param  points The points' numbers.
   *  @return        The polylines.
   */
  auto Sidelong(const Work &work, const Nodes &points) -> std::vector<Polyline>;
} // namespace Geoxy::Belt

#endif /* GEOXY_BELT_JOIN_HPP */
