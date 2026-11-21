/***************************************************************************************************
 *  @file      geoxy/reverse.hpp
 *  @brief     Each curve run the other way.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_REVERSE_HPP
#define GEOXY_REVERSE_HPP

#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <optional>

namespace Geoxy
{
  /**
   *  A line the other way.
   *
   *  @param  line The line.
   *  @return      Its ends swapped: At(Reversed(line), u) is At(line, 1 - u). Twice, the same line.
   */
  [[nodiscard]] GEOXY_API auto Reversed(const Line &line) -> Line;

  /**
   *  An xline the other way.
   *
   *  @param  xline The xline.
   *  @return       Its way negated: At(Reversed(xline), u) is At(xline, -u). Twice, the same xline.
   */
  [[nodiscard]] GEOXY_API auto Reversed(const Xline &xline) -> Xline;

  /**
   *  An arc the other way.
   *
   *  Its start is the end angle At sums, start + sweep, so it starts at At(arc, 1) to the bit. An
   *  arc's end is not held, but summed, so twice reversed it is the same arc to the bit only where
   *  its start and sweep add exactly, as whole degrees and their binary fractions do; else its
   *  start is within a step of a double of its own.
   *
   *  @param  arc The arc.
   *  @return     The same circle, from its end angle, its sweep negated.
   */
  [[nodiscard]] GEOXY_API auto Reversed(const Arc &arc) -> Arc;

  /**
   *  An ellipse the other way, as an arc is reversed: it starts at At(ellipse, 1) to the bit.
   *
   *  @param  ellipse The ellipse.
   *  @return         The same map, from its end parameter, its sweep negated.
   */
  [[nodiscard]] GEOXY_API auto Reversed(const Ellipse &ellipse) -> Ellipse;

  /**
   *  A polyline the other way.
   *
   *  Its vertices in the other order, each segment's bulge negated and held by the vertex the
   *  segment now leaves; a closed polyline keeps its first vertex first. Twice, the same polyline:
   *  an open one's last bulge, which joins nothing, is negated and kept for that.
   *
   *  @param  polyline The polyline.
   *  @return          The polyline: At(Reversed(polyline), u) is At(polyline, n - u), n its count
   *                   of segments.
   */
  [[nodiscard]] GEOXY_API auto Reversed(const Polyline &polyline) -> Polyline;

  /**
   *  A spline the other way.
   *
   *  Its poles and weights in the other order, and its knots negated and in the other order, which
   *  rounds nothing: twice, the same spline to the bit.
   *
   *  @param  spline The spline.
   *  @return        The spline: At(Reversed(spline), -t) is At(spline, t), and its domain is the
   *                 spline's negated.
   */
  [[nodiscard]] GEOXY_API auto Reversed(const Spline &spline) -> Spline;

  /**
   *  A quadratic Bézier curve the other way.
   *
   *  @param  quadratic The curve.
   *  @return           Its control points in the other order: At(Reversed(q), u) is At(q, 1 - u).
   */
  [[nodiscard]] GEOXY_API auto Reversed(const Bezier::Quadratic &quadratic) -> Bezier::Quadratic;

  /**
   *  A cubic Bézier curve the other way.
   *
   *  @param  cubic The curve.
   *  @return       Its control points in the other order: At(Reversed(c), u) is At(c, 1 - u).
   */
  [[nodiscard]] GEOXY_API auto Reversed(const Bezier::Cubic &cubic) -> Bezier::Cubic;

  /**
   *  Any curve the other way.
   *
   *  @param  curve The curve.
   *  @return       It reversed as its own kind is; none for a ray, which reversed would start
   *                nowhere.
   */
  [[nodiscard]] GEOXY_API auto Reversed(const Curve &curve) -> std::optional<Curve>;
} // namespace Geoxy

#endif /* GEOXY_REVERSE_HPP */
