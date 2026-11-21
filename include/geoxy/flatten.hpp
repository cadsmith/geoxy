/***************************************************************************************************
 *  @file      geoxy/flatten.hpp
 *  @brief     Curves as chords: points along a curve, joined by lines that stay within a tolerance
 *             of it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_FLATTEN_HPP
#define GEOXY_FLATTEN_HPP

#include "geoxy/accuracy.hpp"
#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/spline.hpp"

#include <cstddef>
#include <vector>

namespace Geoxy
{
  /**
   *  The most chords a curve, or a polyline's segment, is flattened into: past it the accuracy is
   *  not met, and the call still ends.
   */
  inline constexpr std::size_t MostChords = 65'536;

  /**
   *  A line as chords: no point of it farther than a tolerance from them.
   *
   *  @param  line     The line.
   *  @param  accuracy How far the curve may stray from its chords; one not positive takes
   *                   MostChords.
   *  @return          Its two ends.
   */
  [[nodiscard]] GEOXY_API auto Flatten(const Line &line, Accuracy accuracy) -> std::vector<Point>;

  /**
   *  An arc as chords: no point of it farther than a tolerance from them.
   *
   *  @param  arc      The arc.
   *  @param  accuracy How far the curve may stray from its chords; one not positive takes
   *                   MostChords.
   *  @return          Its points at even steps of its sweep, as few as keep every chord within the
   *                   tolerance and a quarter turn, its ends exact.
   */
  [[nodiscard]] GEOXY_API auto Flatten(const Arc &arc, Accuracy accuracy) -> std::vector<Point>;

  /**
   *  An ellipse as chords: no point of it farther than a tolerance from them.
   *
   *  @param  ellipse  The ellipse.
   *  @param  accuracy How far the curve may stray from its chords; one not positive takes
   *                   MostChords.
   *  @return          Its points, its ends as At gives them: in each quarter between the directions
   *                   of its axes, runs of even steps of its parameter, each the longest whose
   *                   chords stray no more than the tolerance by the exact sagitta of an affine
   *                   circle, sized where the run is nearer a tip of the major axis; one chord a
   *                   quarter at least. MostChords even steps where those would be more.
   */
  [[nodiscard]] GEOXY_API auto Flatten(const Ellipse &ellipse, Accuracy accuracy)
      -> std::vector<Point>;

  /**
   *  A polyline as chords: no point of it farther than a tolerance from them.
   *
   *  @param  polyline The polyline.
   *  @param  accuracy How far the curve may stray from its chords; one not positive takes
   *                   MostChords.
   *  @return          Its vertices exactly, and between them each bulged segment's points at even
   *                   steps of its sweep, worked along its chord; none for no vertices.
   */
  [[nodiscard]] GEOXY_API auto Flatten(const Polyline &polyline, Accuracy accuracy)
      -> std::vector<Point>;

  /**
   *  A spline as chords: no point of it farther than a tolerance from them.
   *
   *  @param  spline   The spline.
   *  @param  accuracy How far the curve may stray from its chords; one not positive takes
   *                   MostChords.
   *  @return          A polynomial spline's of degree two or three span by span, each span's piece
   *                   as a Bézier curve's chords are taken, its ends as At gives them; any other's
   *                   pieces halved until each is within the tolerance of its chord, the ends of
   *                   the halves in order; none for a spline that is not valid.
   */
  [[nodiscard]] GEOXY_API auto Flatten(const Spline &spline, Accuracy accuracy)
      -> std::vector<Point>;

  /**
   *  A quadratic Bézier curve as chords: no point of it farther than a tolerance from them.
   *
   *  @param  quadratic The curve.
   *  @param  accuracy  How far the curve may stray from its chords; one not positive takes
   *                    MostChords.
   *  @return           Its points, its ends exactly: from each, the longest chord of a few tries
   *                    whose stretch of the curve strays no more than the tolerance, by its exact
   *                    distance; halved until each half is within the tolerance of its chord where
   *                    the tolerance is not past rounding, or those would be more than MostChords.
   */
  [[nodiscard]] GEOXY_API auto Flatten(const Bezier::Quadratic &quadratic, Accuracy accuracy)
      -> std::vector<Point>;

  /**
   *  A cubic Bézier curve as chords: no point of it farther than a tolerance from them.
   *
   *  @param  cubic    The curve.
   *  @param  accuracy How far the curve may stray from its chords; one not positive takes
   *                   MostChords.
   *  @return          Its points, its ends exactly: from each, the longest chord of a few tries
   *                   whose stretch of the curve strays no more than the tolerance, by its exact
   *                   distance, the curve cut where it turns hardest; halved until each half is
   *                   within the tolerance of its chord where the tolerance is not past rounding,
   *                   or those would be more than MostChords.
   */
  [[nodiscard]] GEOXY_API auto Flatten(const Bezier::Cubic &cubic, Accuracy accuracy)
      -> std::vector<Point>;

  /**
   *  Any curve as chords: no point of it farther than a tolerance from them.
   *
   *  @param  curve    The curve.
   *  @param  accuracy How far the curve may stray from its chords; one not positive takes
   *                   MostChords.
   *  @return          Its points as its own kind gives them; none for a ray or an xline, which
   *                   never end: clip one to a box first.
   */
  [[nodiscard]] GEOXY_API auto Flatten(const Curve &curve, Accuracy accuracy) -> std::vector<Point>;

} // namespace Geoxy

#endif /* GEOXY_FLATTEN_HPP */
