/***************************************************************************************************
 *  @file      source/knots.hpp
 *  @brief     A spline's poles worked out again over other knots, by its blossoms: knots inserted,
 *             a part cut out, a degree raised.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_KNOTS_HPP
#define GEOXY_KNOTS_HPP

#include "geoxy/spline.hpp"

#include <span>

namespace Geoxy
{
  /**
   *  A valid spline over other knots, of its own degree or of one more: the same curve over them.
   *
   *  Every pole of the new spline is the blossom at its knots, of degree one more the mean of the
   *  blossoms with one of them left out, each taken over the span of the old spline that holds a
   *  span of the new one's support within its domain. So the knots must refine the old ones within
   *  the new domain, and the new domain lie within the old one: knots inserted, a part's knots
   *  clamped at its ends, or every knot repeated once more for a degree more.
   *
   *  @param  spline The spline, valid.
   *  @param  knots  The new knots.
   *  @param  raise  0 to keep the degree, 1 to raise it by one.
   *  @return        The spline over the new knots, weighted where the old one is.
   */
  auto Respanned(const Spline &spline, std::span<const double> knots, int raise) -> Spline;
} // namespace Geoxy

#endif /* GEOXY_KNOTS_HPP */
