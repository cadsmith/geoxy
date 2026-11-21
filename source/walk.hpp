/***************************************************************************************************
 *  @file      source/walk.hpp
 *  @brief     The segment of a polyline a value of its parameter falls in.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_WALK_HPP
#define GEOXY_WALK_HPP

#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace Geoxy
{
  /**
   *  A polyline's segment, as its two vertices and its bulge, and how far along it a value falls.
   */
  struct Step
  {
    Point  Start    = {};  /**< The segment's start.                */
    Point  End      = {};  /**< The segment's end.                  */
    double Bulge    = 0.0; /**< The segment's bulge.                */
    double Fraction = 0.0; /**< How far along it, past 0 and 1 too. */
  };

  /**
   *  The segment a value falls in: its index the whole part, the first before it and the last past
   *  it.
   *
   *  @pre             The polyline has a segment.
   *  @param  polyline The polyline.
   *  @param  along    The value.
   *  @return          The segment and the fraction along it.
   */
  inline auto Locate(const Polyline &polyline, double along) -> Step
  {
    const auto    last  = static_cast<double>(Count(polyline) - 1);
    const double  index = along > 0.0 ? std::min(std::floor(along), last) : 0.0;
    const auto    first = static_cast<std::size_t>(index);
    const Vertex &start = polyline.Vertices[first];
    const Vertex &end   = polyline.Vertices[(first + 1) % polyline.Vertices.size()];
    return {.Start    = start.Position,
            .End      = end.Position,
            .Bulge    = start.Bulge,
            .Fraction = along - index};
  }

  /**
   *  Whether a step's segment is a line, as Bent decides it.
   *
   *  @param  step The step.
   *  @return      Whether its bulge is under Straight in magnitude or its ends are one point.
   */
  inline auto Straightness(const Step &step) -> bool
  {
    return !(std::abs(step.Bulge) >= Straight) || step.Start == step.End;
  }

} // namespace Geoxy

#endif /* GEOXY_WALK_HPP */
