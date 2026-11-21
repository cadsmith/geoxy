/***************************************************************************************************
 *  @file      source/belt/corner.hpp
 *  @brief     The corners of a band's edge: the round about a vertex, its wedge and its mitre.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BELT_CORNER_HPP
#define GEOXY_BELT_CORNER_HPP

#include "belt/work.hpp"
#include "twofold.hpp"

#include <array>
#include <cstddef>

namespace Geoxy::Belt
{
  /**
   *  A piece's heading at one of its ends in two doubles a number: its chord turned by half its
   *  sweep, worked out from its bulge with no trigonometry.
   *
   *  @param  work  The polyline's vertices.
   *  @param  piece The piece.
   *  @param  ends  Whether at its end, not its start.
   *  @return       The heading's x and y, of any length.
   */
  auto Headway(const Work &work, std::size_t piece, bool ends) -> std::array<Twofold::Double, 2>;

  /**
   *  A way brought to a length of 1, in two doubles a number.
   *
   *  @param  way The way, not of no length.
   *  @return     The same way, as long as 1.
   */
  auto Unit(const std::array<Twofold::Double, 2> &way) -> std::array<Twofold::Double, 2>;

  /**
   *  The circles about the vertices the shape turns at: an arc round the outer side of each, from
   *  where one moved piece ends to where the next starts; or, with corners extended or chamfered,
   *  the corner's mitre, which holds that arc. A vertex it turns right back at has its circle on
   *  its right, and its mitre on its outer side as its headings or its pieces' bends tell it, none
   *  where they tell neither.
   *
   *  @param work The shape; the curves and the wedges are added.
   */
  void Cornered(Work &work);
} // namespace Geoxy::Belt

#endif /* GEOXY_BELT_CORNER_HPP */
