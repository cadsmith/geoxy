/***************************************************************************************************
 *  @file      geoxy/dash.hpp
 *  @brief     A linetype's dashes, gaps and dots laid along a curve by length.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_DASH_HPP
#define GEOXY_DASH_HPP

#include "geoxy/accuracy.hpp"
#include "geoxy/api.hpp"
#include "geoxy/curve.hpp"

#include <cstddef>
#include <vector>

namespace Geoxy
{
  /**
   *  The most times a pattern is laid along one walk: past it the walk is drawn whole, as a pattern
   *  too fine to see is.
   */
  inline constexpr std::size_t MostPatterns = 1'048'576;

  /**
   *  A linetype's pattern, as DXF's LTYPE holds it, scaled as it is to be drawn.
   */
  struct Pattern
  {
    std::vector<double> Elements = {};   /**< A dash positive, a gap negative, a dot 0.        */
    double              Phase    = 0.0;  /**< How far into the pattern the curve's start lies. */
    bool                Across   = true; /**< Whether it runs on across a polyline's vertices,
                                              as PLINEGEN on; else it starts again at each.    */
  };

  /**
   *  The dashes of a linetype laid along a curve by length, from its start.
   *
   *  The pattern repeats from its phase, every dash cut from the curve between the places at the
   *  lengths where it starts and ends, the last one where the curve ends; each a curve of the
   *  curve's own kind, its ends the points At gives at those places, which are worked out by length
   *  as PlaceAt finds them. A polyline is walked segment by segment and a spline span by span, so
   *  that a long one costs a pass, not a search per dash. With Across off, each segment of a
   *  polyline is laid on its own, from the pattern's start.
   *
   *  A walk shorter than one whole pattern is drawn whole, as AutoCAD draws it, and so is one along
   *  which the pattern would be laid more than MostPatterns times.
   *
   *  @param  curve    The curve.
   *  @param  pattern  The pattern.
   *  @param  accuracy How far a dash's length may be off, for a curve whose length is approximated.
   *  @return          The dashes, in order, a dot as a line of no length; the curve whole for a
   *                   pattern of no length or one not finite, a ray or an xline, which never end,
   *                   and a curve of no segment or span to walk.
   */
  [[nodiscard]] GEOXY_API auto Dashes(const Curve &curve, const Pattern &pattern, Accuracy accuracy)
      -> std::vector<Curve>;
} // namespace Geoxy

#endif /* GEOXY_DASH_HPP */
