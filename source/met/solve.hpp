/***************************************************************************************************
 *  @file      source/met/solve.hpp
 *  @brief     Two parts met by the rule for their kinds, two runs, a run and a piece, two pieces,
 *             and any two parts by the rule for theirs.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_MET_SOLVE_HPP
#define GEOXY_MET_SOLVE_HPP

#include "met/parts.hpp"
#include "straight.hpp"

#include "geoxy/domain.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"

#include <cstddef>
#include <optional>

namespace Geoxy::Met
{
  /**
   *  A point of a piece.
   *
   *  @param  piece The piece.
   *  @param  along Its parameter.
   *  @return       The point.
   */
  auto Evaluated(const Piece &piece, double along) -> Point;

  /**
   *  The angle on an ellipse's unit circle of a point.
   *
   *  @param  unit  The map that takes the ellipse to the unit circle.
   *  @param  point The point.
   *  @return       The angle, in radians.
   */
  auto Angled(const Matrix &unit, Point point) -> double;

  /**
   *  Whether two runs lie on one line within the tolerance over the stretch of the first that the
   *  second's projection covers: that stretch's ends each within its length of the second's line,
   *  and so all of it, and longer than the tolerance, a shorter one being where their ends meet; an
   *  endless stretch only where the lines are parallel exactly. Neither run is of no length.
   *
   *  @param  first     One run.
   *  @param  second    The other.
   *  @param  tolerance The tolerance.
   *  @return           Whether they are, and the stretch of the first, from and to.
   */
  auto Aligned(const Run &first, const Run &second, const Tolerance &tolerance)
      -> std::optional<Range>;

  /**
   *  Where two runs meet: one stretch where they lie on one line within the tolerance, the second's
   *  points at its ends within twice it of the first's line, else the crossing of their lines if it
   *  lies on both.
   *
   *  @param  first     One run.
   *  @param  second    The other.
   *  @param  tolerance How near two lines must be to be one, and how far past an end a crossing may
   *                    lie.
   *  @return           What they found.
   */
  auto Against(const Run &first, const Run &second, const Tolerance &tolerance) -> Found;

  /**
   *  Where a run meets a piece: the stretch they share where the piece is straight and lies on the
   *  run's line within the tolerance, as two runs would; where it lies on the line to rounding, the
   *  stretch too; else the zeros of the run's line on it.
   *
   *  @param  run       The run.
   *  @param  piece     The piece.
   *  @param  tolerance How near a straight piece must be to the run's line to share it.
   *  @return           What they found, the run first.
   */
  auto Against(const Run &run, const Piece &piece, const Tolerance &tolerance) -> Found;

  /**
   *  How near two lines come: nothing where each has the other's ends on both sides of it or on it,
   *  decided exactly; else, as where all four ends are on one line, the nearest of each one's ends
   *  to the other.
   *
   *  @param  one   One line.
   *  @param  other The other.
   *  @return       The distance.
   */
  auto Apart(const Line &one, const Line &other) -> double;

  /**
   *  The candidate meetings of two pieces, by subdivision: pairs whose hulls come within the
   *  tolerance halved, the larger first, until both are nearly straight, their chords then met for
   *  a start that Newton's method polishes later. A part of the first that lies along the second,
   *  sampled, is a stretch they share, and is not looked into further; a part found not to is not
   *  sampled again while only the second is halved, nor where the two part, as a stretch along this
   *  part of the second would run its way. A pair whose chords are farther apart than the tolerance
   *  and how far the pieces stray from them has nothing. A pair is halved no further, its chords
   *  met as the start for all of it, where every hit in it is one: where the pieces are within the
   *  tolerance of each other all along and too short to share a stretch, or where both are within
   *  the tolerance of their chords and part at more than its angle, so that they cross at most once
   *  and touch nowhere. Else two curves leaving one tangent, or crossing, stay within a loose
   *  tolerance a long way and would be halved until the budget ran out.
   *
   *  @param  first     One piece.
   *  @param  second    The other.
   *  @param  tolerance How near the hulls must come.
   *  @return           What they found, and whether the budget ran out.
   */
  auto Against(const Piece &first, const Piece &second, const Tolerance &tolerance) -> Found;

  /**
   *  What two parts found, whatever their shapes, each pair by its own rule.
   *
   *  @param  first     One part's shape.
   *  @param  second    The other's.
   *  @param  tolerance How near they must come.
   *  @return           What they found, in the parts' own parameters.
   */
  auto Matched(const Part &first, const Part &second, const Tolerance &tolerance) -> Found;
} // namespace Geoxy::Met

#endif /* GEOXY_MET_SOLVE_HPP */
