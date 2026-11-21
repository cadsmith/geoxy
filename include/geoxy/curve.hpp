/***************************************************************************************************
 *  @file      geoxy/curve.hpp
 *  @brief     Any one curve: the closed set of the curves a drawing is made of.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CURVE_HPP
#define GEOXY_CURVE_HPP

#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <variant>

namespace Geoxy
{
  /**
   *  Any one curve, as a selection holds them.
   *
   *  A closed set: every operation is written for each kind, and std::visit chooses, so no kind is
   *  ever forgotten and nothing allocates but a polyline's vertices and a spline's poles.
   */
  using Curve = std::variant<Line, Ray, Xline, Arc, Ellipse, Polyline, Spline, Bezier::Quadratic,
                             Bezier::Cubic>;
} // namespace Geoxy

#endif /* GEOXY_CURVE_HPP */
