/***************************************************************************************************
 *  @file      geoxy/domain.hpp
 *  @brief     The range of a curve's parameter.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_DOMAIN_HPP
#define GEOXY_DOMAIN_HPP

#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <functional>
#include <limits>

namespace Geoxy
{
  /**
   *  A stretch of a parameter, from where it starts to where it ends.
   */
  struct Range
  {
    double Start = 0.0; /**< Where it starts. */
    double End   = 0.0; /**< Where it ends.   */

    /**
     *  Whether two ranges are equal exactly.
     *
     *  @param  lhs One range.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Range &lhs, const Range &rhs) -> bool
    {
      return std::equal_to<double>{}(lhs.Start, rhs.Start)
          && std::equal_to<double>{}(lhs.End, rhs.End);
    }
  };

  /**
   *  The parameter a line runs over: 0 at its start, 1 at its end.
   *
   *  @return The range.
   */
  [[nodiscard]] constexpr auto Domain(const Line & /*line*/) -> Range
  {
    return {.Start = 0.0, .End = 1.0};
  }

  /**
   *  The parameter a ray runs over, in its ways: from 0 at its start, without end.
   *
   *  @return The range.
   */
  [[nodiscard]] constexpr auto Domain(const Ray & /*ray*/) -> Range
  {
    return {.Start = 0.0, .End = std::numeric_limits<double>::infinity()};
  }

  /**
   *  The parameter an xline runs over, in its ways, 0 at its point: without end either way.
   *
   *  @return The range.
   */
  [[nodiscard]] constexpr auto Domain(const Xline & /*xline*/) -> Range
  {
    return {.Start = -std::numeric_limits<double>::infinity(),
            .End   = std::numeric_limits<double>::infinity()};
  }

  /**
   *  The parameter an arc runs over, a fraction of its sweep: 0 at its start, 1 at its end.
   *
   *  @return The range.
   */
  [[nodiscard]] constexpr auto Domain(const Arc & /*arc*/) -> Range
  {
    return {.Start = 0.0, .End = 1.0};
  }

  /**
   *  The parameter an ellipse runs over, a fraction of its sweep: 0 at its start, 1 at its end.
   *
   *  @return The range.
   */
  [[nodiscard]] constexpr auto Domain(const Ellipse & /*ellipse*/) -> Range
  {
    return {.Start = 0.0, .End = 1.0};
  }

  /**
   *  The parameter a quadratic Bézier curve runs over: 0 at its start, 1 at its end.
   *
   *  @return The range.
   */
  [[nodiscard]] constexpr auto Domain(const Bezier::Quadratic & /*quadratic*/) -> Range
  {
    return {.Start = 0.0, .End = 1.0};
  }

  /**
   *  The parameter a cubic Bézier curve runs over: 0 at its start, 1 at its end.
   *
   *  @return The range.
   */
  [[nodiscard]] constexpr auto Domain(const Bezier::Cubic & /*cubic*/) -> Range
  {
    return {.Start = 0.0, .End = 1.0};
  }

  /**
   *  The parameter a polyline runs over: a segment's index and the fraction along it, from 0 at its
   *  first vertex to its count of segments at its end.
   *
   *  @param  polyline The polyline.
   *  @return          The range.
   */
  [[nodiscard]] GEOXY_API auto Domain(const Polyline &polyline) -> Range;

  /**
   *  The parameter a spline runs over: its knots, from the one at its degree to the one past its
   *  last pole.
   *
   *  @param  spline The spline.
   *  @return        The range; from 0 to 0 for a spline whose counts do not agree.
   */
  [[nodiscard]] GEOXY_API auto Domain(const Spline &spline) -> Range;

  /**
   *  The parameter any curve runs over.
   *
   *  @param  curve The curve.
   *  @return       The range.
   */
  [[nodiscard]] GEOXY_API auto Domain(const Curve &curve) -> Range;
} // namespace Geoxy

#endif /* GEOXY_DOMAIN_HPP */
