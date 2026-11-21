/***************************************************************************************************
 *  @file      source/sound.hpp
 *  @brief     Whether every number a curve is held by is finite.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SOUND_HPP
#define GEOXY_SOUND_HPP

#include "geoxy/curve.hpp"

namespace Geoxy
{
  /**
   *  Whether every number a curve is held by is finite: its points, ways, radius and angles, a
   *  polyline's bulges, a spline's knots and weights.
   *
   *  @param  curve The curve.
   *  @return       Whether they are; a ray or an endless line whose numbers are finite is sound.
   */
  auto Sound(const Curve &curve) -> bool;
} // namespace Geoxy

#endif /* GEOXY_SOUND_HPP */
