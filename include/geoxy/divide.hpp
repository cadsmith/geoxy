/***************************************************************************************************
 *  @file      geoxy/divide.hpp
 *  @brief     Places at equal lengths along a curve, as DIVIDE, MEASURE and ARRAY path set them.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_DIVIDE_HPP
#define GEOXY_DIVIDE_HPP

#include "geoxy/accuracy.hpp"
#include "geoxy/api.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/place.hpp"

#include <cstddef>
#include <vector>

namespace Geoxy
{
  /**
   *  The most places Divide and Measure set along one curve.
   */
  inline constexpr std::size_t MostPlaces = 16'777'216;

  /**
   *  The places that cut a curve into parts of equal length, as DIVIDE sets them: the length times
   *  k over the count, for k from 1 to one less than the count.
   *
   *  The curve is walked once, a polyline segment by segment and a spline span by span, each place
   *  found by length within the piece it falls in, as PlaceAt finds it there.
   *
   *  @param  curve    The curve.
   *  @param  parts    How many parts.
   *  @param  accuracy How far each place's length from the start may be off.
   *  @return          The places, in order; none for fewer than two parts or more than MostPlaces,
   *                   a curve without end, or one of no segment or span.
   */
  [[nodiscard]] GEOXY_API auto Divide(const Curve &curve, std::size_t parts, Accuracy accuracy)
      -> std::vector<Place>;

  /**
   *  The places a length apart along a curve from its start, as MEASURE sets them: every multiple
   *  of the length up to the curve's.
   *
   *  @param  curve    The curve.
   *  @param  length   The length.
   *  @param  accuracy How far each place's length from the start may be off.
   *  @return          The places, in order; none for a length not above 0 or not finite, more than
   *                   MostPlaces of them, a curve without end, or one of no segment or span.
   */
  [[nodiscard]] GEOXY_API auto Measure(const Curve &curve, double length, Accuracy accuracy)
      -> std::vector<Place>;
} // namespace Geoxy

#endif /* GEOXY_DIVIDE_HPP */
