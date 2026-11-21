/***************************************************************************************************
 *  @file      source/belt/end.hpp
 *  @brief     The ends of an open polyline's band: round, flat or square.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BELT_END_HPP
#define GEOXY_BELT_END_HPP

#include "belt/work.hpp"

#include <cstddef>

namespace Geoxy::Belt
{
  /**
   *  The band's edge round an open shape's two ends: half circles, the lines across, or squares.
   *
   *  @param work The shape; the curves and the wedges are added.
   */
  void Ended(Work &work);
} // namespace Geoxy::Belt

#endif /* GEOXY_BELT_END_HPP */
