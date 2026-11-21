/***************************************************************************************************
 *  @file      source/belt/keep.hpp
 *  @brief     The stretches of a band's curves that lie on its edge, kept.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BELT_KEEP_HPP
#define GEOXY_BELT_KEEP_HPP

#include "belt/work.hpp"

#include <cstddef>

namespace Geoxy::Belt
{
  /**
   *  The stretches between each curve's cuts that lie on the region's edge.
   *
   *  @param work   The curves and their cuts, which are put in order; the stretches are set.
   *  @param points The points' numbers.
   */
  void Kept(Work &work, const Nodes &points);

  /**
   *  A stretch kept twice, by two curves along one line or one circle, kept once: the stretches are
   *  put in order of their ends, and one whose ends are another's is dropped. Two that are not one
   *  curve cannot both run from one point to another with the region on their left.
   *
   *  @param work The curves and the stretches.
   */
  void Once(Work &work);
} // namespace Geoxy::Belt

#endif /* GEOXY_BELT_KEEP_HPP */
