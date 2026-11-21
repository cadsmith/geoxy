/***************************************************************************************************
 *  @file      source/plain.hpp
 *  @brief     Two ready lines met by their exact pass alone, where it tells them, with no marking
 *             made: what meeting many pairs one after another asks for.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_PLAIN_HPP
#define GEOXY_PLAIN_HPP

#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/tolerance.hpp"

#include <cstddef>

namespace Geoxy
{
  /**
   *  What two ready lines meet in, where it is plain: nothing, or one hit, the double nearest the
   *  carriers' crossing, a crossing or an end of one on the other; else that it is not, and what
   *  they meet in is written into a marking.
   */
  struct Plain
  {
    Point       At    = {};             /**< The hit, where they meet.                       */
    std::size_t Count = 0;              /**< How many hits: none, or one.                    */
    Contact     Kind  = Contact::Cross; /**< How they meet there.                            */
    bool        Told  = false;          /**< Whether it is plain; else it is in the marking. */
  };

  /**
   *  Two ready lines met as Intersect meets them made ready, to the bit, with nothing made where
   *  their exact passes tell them: a hit then needs no more than its point and kind, and nothing no
   *  marking.
   *
   *  Most pairs of a drawing's lines that come near cross plainly, or plainly do not, and the pass
   *  tells them in a few products; a level line and an upright one at a corner or a joint, an end
   *  of one on the other, Abutting tells from their gaps. The rest go on as Intersect goes on with
   *  them, into the marking: an end near the other line but not on it, a narrow angle, a crossing
   *  exactly halfway between two doubles, or one on an axis. Nothing is worked out twice.
   *
   *  @pre          Both curves are lines.
   *  @param  lines The two ready lines, the first first.
   *  @param  sane  The tolerance, as Sane makes it.
   *  @param  into  The marking, whatever it held replaced where it is not told, and left where it
   *                is.
   *  @return       What they meet in, where it is told.
   */
  auto Plainly(const Prepared::Pair &lines, const Tolerance &sane, Marking &into) -> Plain;
} // namespace Geoxy

#endif /* GEOXY_PLAIN_HPP */
