/***************************************************************************************************
 *  @file      geoxy/outline.hpp
 *  @brief     The outline of a wide polyline, as TRACE and a wide PLINE are drawn: a band per
 *             segment, its corners mitered.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_OUTLINE_HPP
#define GEOXY_OUTLINE_HPP

#include "geoxy/accuracy.hpp"
#include "geoxy/api.hpp"
#include "geoxy/polyline.hpp"

#include <span>
#include <vector>

namespace Geoxy
{
  /**
   *  A segment's width where it starts and where it ends, as LWPOLYLINE holds them per vertex.
   */
  struct Width
  {
    double Start = 0.0; /**< Its width at its start. */
    double End   = 0.0; /**< Its width at its end.   */
  };

  /**
   *  The outline of a wide polyline: one closed polyline per segment, counterclockwise, its band.
   *
   *  A line's band is a trapezoid, its edges each half the width from it either side. An arc's of
   *  one width is two concentric arcs, the polyline's bulge on both; one whose width changes along
   *  it has edges that are spirals, drawn as chords within the accuracy, as AutoCAD draws them.
   *
   *  Where two straight segments meet, each side's edges run on to where they cross, and both bands
   *  end there: mitered as AutoCAD and TRACE miter a corner, however sharp. A corner with an arc,
   *  or between parallel edges, is butt: each band ends square at the vertex, seamless where the
   *  arc is tangent and the widths there agree.
   *
   *  @param  polyline The polyline.
   *  @param  widths   One width for every segment, or one per segment, the last going on for those
   *                   past them; a negative width counts as its magnitude.
   *  @param  accuracy How far a chord of a changing arc's edge may stray from it; at or below 0, as
   *                   many chords as Flatten draws at most.
   *  @return          The bands, in the segments' order; none for a segment of no width at either
   *                   end or of no length, or for no widths.
   */
  [[nodiscard]] GEOXY_API auto Outline(const Polyline &polyline, std::span<const Width> widths,
                                       Accuracy accuracy) -> std::vector<Polyline>;
} // namespace Geoxy

#endif /* GEOXY_OUTLINE_HPP */
