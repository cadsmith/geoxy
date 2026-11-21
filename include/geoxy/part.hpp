/***************************************************************************************************
 *  @file      geoxy/part.hpp
 *  @brief     A curve's stretch between two places, a curve cut at or run on to a place, and a
 *             curve split at places: TRIM, EXTEND, BREAK and LENGTHEN.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_PART_HPP
#define GEOXY_PART_HPP

#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace Geoxy
{
  /**
   *  Which end of a curve a trim keeps.
   */
  enum class Keep : std::uint8_t
  {
    Start, /**< The part from the curve's start to the place. */
    End,   /**< The part from the place to the curve's end.   */
  };

  /**
   *  The stretch of a line between two places.
   *
   *  @param  line    The line.
   *  @param  stretch Where it starts and ends, past 0 and 1 the line run on, the other way round
   *                  when its start is after its end.
   *  @return         The line between the points At gives there: At(part, u) is At(line, s + u (e -
   *                  s)) to a rounding, and exactly so at its ends.
   */
  [[nodiscard]] GEOXY_API auto Part(const Line &line, Range stretch) -> Line;

  /**
   *  The stretch of a ray between two places, on the line it runs along.
   *
   *  @param  ray     The ray.
   *  @param  stretch Where it starts and ends, before 0 the ray run back.
   *  @return         The line between the points At gives where both are finite; a ray from the
   *                  finite one toward the other's infinity where one is infinite, held from its
   *                  finite end, so the other way round when that is its end; the xline of its
   *                  line, the way from the start's infinity to the end's, where both are.
   */
  [[nodiscard]] GEOXY_API auto Part(const Ray &ray, Range stretch) -> Curve;

  /**
   *  The stretch of an xline between two places.
   *
   *  @param  xline   The xline.
   *  @param  stretch Where it starts and ends.
   *  @return         As a ray's part: a line, a ray or an xline.
   */
  [[nodiscard]] GEOXY_API auto Part(const Xline &xline, Range stretch) -> Curve;

  /**
   *  The stretch of an arc between two places.
   *
   *  Its start is the angle At sums at the stretch's start, so it starts at At(arc, s) to the bit,
   *  and its sweep the difference of the angles At sums at either end. Added back to the start, as
   *  At adds it, that comes to the end's angle, and the part ends at At(arc, e) to the bit,
   *  wherever the difference is exact, as it is for angles within a factor of two of each other and
   *  for angles on one grid of doubles, such as whole degrees and their binary fractions, and
   *  wherever it is no larger than the end's angle. Else, where the end's angle is a finer double
   *  than the start and the sweep can sum to, it comes a step of a double short of it or past it:
   *  about one stretch in four of every digit, in radians.
   *
   *  @param  arc     The arc.
   *  @param  stretch Where it starts and ends, past 0 and 1 its circle run on, past a turn too.
   *  @return         The arc: At(part, u) is At(arc, s + u (e - s)) to a rounding.
   */
  [[nodiscard]] GEOXY_API auto Part(const Arc &arc, Range stretch) -> Arc;

  /**
   *  The stretch of an ellipse between two places, as an arc's is cut, its map kept.
   *
   *  @param  ellipse The ellipse.
   *  @param  stretch Where it starts and ends, past 0 and 1 the ellipse run on.
   *  @return         The ellipse: At(part, u) is At(ellipse, s + u (e - s)) to a rounding, and at
   *                  its ends as exactly as an arc's part is.
   */
  [[nodiscard]] GEOXY_API auto Part(const Ellipse &ellipse, Range stretch) -> Ellipse;

  /**
   *  The stretch of a polyline between two places.
   *
   *  Its vertices are the point At gives at the stretch's start, the polyline's vertices between,
   *  and the point At gives at its end. A segment kept whole keeps its bulge; one cut keeps its
   *  circle, its bulge tan(atan(b) f) for f the fraction kept, and a segment Bent takes for a line
   *  stays one. The part is open, of a closed polyline too.
   *
   *  @param  polyline The polyline.
   *  @param  stretch  Where it starts and ends, before the first segment or past the last that
   *                   segment run on, the other way round when its start is after its end.
   *  @return          The polyline, its segment k the polyline's segment first + k between the
   *                   fractions kept; the polyline as it is where it has no segment; none where a
   *                   bulged segment would be run on through a whole turn of its circle, which no
   *                   bulge holds.
   */
  [[nodiscard]] GEOXY_API auto Part(const Polyline &polyline, Range stretch)
      -> std::optional<Polyline>;

  /**
   *  The stretch of a spline between two places, its parameter kept.
   *
   *  Its knots are the stretch's ends, each a degree and one times, and the spline's knots strictly
   *  between them within its domain; each pole the blossom of the spline's span it lies over, at
   *  the knots it weighs, so At(part, t) is At(spline, t) to a rounding. Its end poles are the
   *  points At gives at the stretch's ends, worked on the spans At reads; a rational spline's end
   *  weights are those At divides by there, so At reads its ends back as its quotient times its
   *  divisor, divided again: the same points, to the bit, in every case the tests try. The other
   *  way round, it is the part reversed, its knots negated.
   *
   *  @param  spline  The spline.
   *  @param  stretch Where it starts and ends, before its first span or past its last that span's
   *                  polynomial run on.
   *  @return         The spline; an empty spline, a curve of nothing, for a stretch of no length,
   *                  not finite, or of a spline that is not valid. Run on where its weights sum to
   *                  nothing or less, it is not valid either.
   */
  [[nodiscard]] GEOXY_API auto Part(const Spline &spline, Range stretch) -> Spline;

  /**
   *  The stretch of a quadratic Bézier curve between two places.
   *
   *  @param  quadratic The curve.
   *  @param  stretch   Where it starts and ends, past 0 and 1 the curve run on.
   *  @return           The curve from and to the points At gives there, its control point the
   *                    blossom at the two: At(part, u) is At(quadratic, s + u (e - s)) to a
   *                    rounding.
   */
  [[nodiscard]] GEOXY_API auto Part(const Bezier::Quadratic &quadratic, Range stretch)
      -> Bezier::Quadratic;

  /**
   *  The stretch of a cubic Bézier curve between two places.
   *
   *  @param  cubic   The curve.
   *  @param  stretch Where it starts and ends, past 0 and 1 the curve run on.
   *  @return         The curve from and to the points At gives there, its control points the
   *                  blossoms at s, s, e and at s, e, e: At(part, u) is At(cubic, s + u (e - s)) to
   *                  a rounding.
   */
  [[nodiscard]] GEOXY_API auto Part(const Bezier::Cubic &cubic, Range stretch) -> Bezier::Cubic;

  /**
   *  The stretch of any curve between two places.
   *
   *  @param  curve   The curve.
   *  @param  stretch Where it starts and ends, as its own kind reads them.
   *  @return         The part, as its own kind gives it; none where a polyline's is none.
   */
  [[nodiscard]] GEOXY_API auto Part(const Curve &curve, Range stretch) -> std::optional<Curve>;

  /**
   *  A curve cut at a place, one side of it kept, as TRIM and BREAK keep it.
   *
   *  @param  curve The curve.
   *  @param  along The place, as its own kind reads it; past an end, the curve is run on to it.
   *  @param  keep  Which end is kept.
   *  @return       The part from the start of its domain to the place, or from the place to the end
   *                of its domain: an xline cut keeping its start is a ray from the place back along
   *                it. None where Part gives none.
   */
  [[nodiscard]] GEOXY_API auto Trim(const Curve &curve, double along, Keep keep)
      -> std::optional<Curve>;

  /**
   *  A curve run on to a place past one of its ends, as EXTEND and LENGTHEN run it.
   *
   *  @param  curve The curve.
   *  @param  along The place, as its own kind reads it.
   *  @return       The curve with its end moved out to the place where the place is past it, its
   *                start where the place is before it; as it is where the place is within its
   *                domain or not a number. None where Part gives none.
   */
  [[nodiscard]] GEOXY_API auto Extend(const Curve &curve, double along) -> std::optional<Curve>;

  /**
   *  A curve split at places, as BREAK AT POINT splits it.
   *
   *  The places within its domain, its ends left out, are taken in order, each once, and the curve
   *  is cut at each. No seam is closed: a closed curve is cut at its start too.
   *
   *  @param  curve  The curve.
   *  @param  alongs The places, as its own kind reads them, in any order.
   *  @return        The parts, in order: one more than the places within its domain.
   */
  [[nodiscard]] GEOXY_API auto Split(const Curve &curve, std::span<const double> alongs)
      -> std::vector<Curve>;
} // namespace Geoxy

#endif /* GEOXY_PART_HPP */
