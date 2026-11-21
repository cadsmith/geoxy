/***************************************************************************************************
 *  @file      source/met/ellipse.hpp
 *  @brief     Two parts met where one is an elliptic arc: with a piece, with a run, with another
 *             elliptic arc, of one ellipse or of two.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_MET_ELLIPSE_HPP
#define GEOXY_MET_ELLIPSE_HPP

#include "met/parts.hpp"
#include "straight.hpp"

#include "geoxy/ellipse.hpp"
#include "geoxy/tolerance.hpp"

#include <cstddef>

namespace Geoxy::Met
{
  /**
   *  Where two arcs of one carrier overlap: the stretches of angle both sweep, each turn of either
   *  tried against the other; none for an arc of no sweep, a point, which the ends find where it
   *  lies on the other.
   *
   *  @param  first  One arc.
   *  @param  second The other, its angles the first's plus a turn offset.
   *  @param  offset What to add to the second's angles to have the first's.
   *  @return        The stretches shared, and where they only touch at an end.
   */
  auto Overlapped(const Ellipse &first, const Ellipse &second, double offset) -> Found;

  /**
   *  Where a piece meets an ellipse: the zeros of the ellipse's implicit equation on it, or, where
   *  the piece lies on the ellipse, the stretches they share.
   *
   *  @param  piece     The piece.
   *  @param  ellipse   The ellipse.
   *  @param  tolerance How near the parts must come, for an ellipse flattened to a line.
   *  @return           What they found, the piece first.
   */
  auto Against(const Piece &piece, const Ellipse &ellipse, const Tolerance &tolerance) -> Found;

  /**
   *  Where a run meets an ellipse, in closed form: the run in the ellipse's unit circle's frame, a
   *  quadratic whose roots are the crossings; and, where they come near without crossing, the
   *  closest approach.
   *
   *  @param  run       The run.
   *  @param  ellipse   The ellipse.
   *  @param  tolerance How near the parts must come, for an ellipse flattened to a line.
   *  @return           What they found, the run first.
   */
  auto Against(const Run &run, const Ellipse &ellipse, const Tolerance &tolerance) -> Found;

  /**
   *  Where two ellipses meet: one carrier, decided exactly, shares stretches; two circles meet in
   *  closed form, and two arcs of circles that part share a stretch between their ends where each
   *  stays within the tolerance of the other's circle all the way; any other pair by one's implicit
   *  equation on the other's quarters.
   *
   *  @param  first     One ellipse.
   *  @param  second    The other.
   *  @param  tolerance How near two circles must stay to be one.
   *  @return           What they found; a budget spent for two ellipses both flattened to lines.
   */
  auto Against(const Ellipse &first, const Ellipse &second, const Tolerance &tolerance) -> Found;
} // namespace Geoxy::Met

#endif /* GEOXY_MET_ELLIPSE_HPP */
