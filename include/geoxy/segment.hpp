/***************************************************************************************************
 *  @file      geoxy/segment.hpp
 *  @brief     A polyline's segments: a line, or an arc by its bulge.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SEGMENT_HPP
#define GEOXY_SEGMENT_HPP

#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"

#include <cstddef>
#include <variant>

namespace Geoxy
{
  /**
   *  One segment of a polyline: a line, or an arc.
   */
  using Segment = std::variant<Line, Arc>;

  /**
   *  The bulge under which a segment is its chord: 2^-60. Such a bulge lifts the arc from the chord
   *  by less than a thousandth of a step of a double at the chord's own size, so the arc is the
   *  chord to every double there is. Above it a segment is the arc its bulge makes it, however
   *  slight. (AutoCAD reads a bulge under 1e-6 as a line; that is its rule, not this library's.)
   */
  inline constexpr double Straight = 0x1p-60;

  /**
   *  The segment from one vertex to the next by the bulge between them.
   *
   *  @param  start Where it starts.
   *  @param  end   Where it ends.
   *  @param  bulge The tangent of a quarter of its sweep, positive counterclockwise.
   *  @return       A line when the bulge is under Straight in magnitude or the ends are one point,
   *                an arc otherwise.
   */
  [[nodiscard]] GEOXY_API auto Bent(Point start, Point end, double bulge) -> Segment;

  /**
   *  How many segments a polyline has.
   *
   *  @param  polyline The polyline.
   *  @return          One less than its vertices when open, as many when closed with two or more,
   *                   none for none.
   */
  [[nodiscard]] GEOXY_API auto Count(const Polyline &polyline) -> std::size_t;

  /**
   *  A polyline's segment.
   *
   *  @pre             The index is under Count(polyline).
   *  @param  polyline The polyline.
   *  @param  index    Which segment, from the first vertex.
   *  @return          The segment from that vertex to the next, the last joining the first when
   *                   closed.
   */
  [[nodiscard]] GEOXY_API auto Nth(const Polyline &polyline, std::size_t index) -> Segment;
} // namespace Geoxy

#endif /* GEOXY_SEGMENT_HPP */
