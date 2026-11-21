/***************************************************************************************************
 *  @file      source/belt/belt.hpp
 *  @brief     A polyline's offset as a region: the edge of its band, everything within a distance
 *             of it with its corners' mitres and its ends, of what it fills less that, or with it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BELT_BELT_HPP
#define GEOXY_BELT_BELT_HPP

#include "geoxy/offset.hpp"
#include "geoxy/polyline.hpp"

#include <vector>

namespace Geoxy::Belt
{
  /**
   *  A polyline's offset as a region.
   *
   *  The band is everything within the distance of the polyline; with corners extended or
   *  chamfered, each outer corner's mitre as well, where it points the way asked; with an open
   *  polyline's ends flat, only what its segments sweep square to themselves and the disks about
   *  the vertices between them, and with square ends a square past each end besides. Its edge can
   *  only run along a segment moved to one side by the distance, a circle of the distance about a
   *  vertex or an end, the outer edge of a mitre, or the line across or the square past an end. So
   *  those curves are made, each run with the band on its left; cut wherever two of them meet; and
   *  a stretch between two cuts is kept where its middle lies outside the band and, for a closed
   *  polyline asked inward or outward, on the side of what the polyline fills that is asked. What
   *  is kept is joined end to start, a loop that comes to a point several leave taking the first
   *  met turning clockwise from the way it came. Doubles decide all of it but what they are too
   *  coarse for: two curves that come near touching are met again in two doubles a number, and a
   *  stretch whose middle lies within a billionth of the shape's size of the distance from another
   *  part of the polyline is measured again likewise.
   *
   *  @param  polyline  The polyline, of finite numbers, no two vertices in a row one point.
   *  @param  buffering The distance, more than 0; the side, both for an open polyline; the filling;
   *                    how outer corners are closed, and an open polyline's ends.
   *  @return           The region's loops, each run with the region on its left.
   */
  auto Edged(const Polyline &polyline, const Buffering &buffering) -> std::vector<Polyline>;

  /**
   *  One side of an open polyline's offset.
   *
   *  The band is everything within the distance of the polyline, round about its ends, with the
   *  mitres of the corners whose outer side is the side asked. Its edge is found as a region's is;
   *  of it, the stretches along the curves that lie on the side asked by their making are the
   *  offset: a segment moved to that side, a circle about a corner whose outer side it is, the
   *  outer edge of such a corner's mitre. Not a circle about an end, nor about a corner that turns
   *  right back, nor a stretch the polyline runs over twice with a different side to it each time.
   *
   *  @param  polyline The polyline, open, of finite numbers, no two vertices in a row one point.
   *  @param  distance How far, not 0: to the right of the way it runs when positive.
   *  @param  gaps     How outer corners are closed.
   *  @return          The offset's polylines, each run the way the polyline runs: open, or closed
   *                   where one comes back to its start.
   */
  auto Flanked(const Polyline &polyline, double distance, Gap gaps) -> std::vector<Polyline>;
} // namespace Geoxy::Belt

#endif /* GEOXY_BELT_BELT_HPP */
