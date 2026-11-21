/***************************************************************************************************
 *  @file      source/met/linear.hpp
 *  @brief     Two straight curves met directly, lines, rays and endless lines: crossing, lying
 *             along each other, or meeting end to end.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_MET_LINEAR_HPP
#define GEOXY_MET_LINEAR_HPP

#include "straight.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"

#include <array>
#include <cstddef>
#include <optional>

namespace Geoxy::Met
{
  /**
   *  The sine below which two curves crossing are left to the thorough way: under it they stay near
   *  each other so long that a band about the crossing may form.
   */
  inline constexpr double Steep = 0.125;

  /**
   *  A straight curve as the direct way takes it: its run, and its ends, a line's two, a ray's
   *  start, none for an endless line.
   */
  struct Straight
  {
    Geoxy::Run            Run    = {};      /**< Its run.            */
    std::array<Point, 2>  Ends   = {};      /**< Its ends.           */
    std::size_t           Count  = 0;       /**< How many there are. */
    const Curve          *Shape  = nullptr; /**< The curve it is.    */
    std::array<double, 2> Places = {};      /**< Each end's place.   */
  };

  /**
   *  A line, a ray or an endless line as the direct way takes it; nothing for any other curve.
   *
   *  @param  curve The curve.
   *  @return       It, straight.
   */
  auto Straightened(const Curve &curve) -> std::optional<Straight>;

  /**
   *  How far a point is from a run: to the foot on it, or to its nearer end.
   *
   *  @param  point The point.
   *  @param  run   The run, of some length.
   *  @return       The distance to its nearest point.
   */
  auto Away(Point point, const Run &run) -> double;

  /**
   *  Two straight curves met directly, where nothing but a crossing or nothing at all can come of
   *  them: each of some length, every end farther than four tolerances from the other, so that no
   *  end is a hit and a crossing lies well within both. Two that part by more than twice the angle
   *  and lie along each other nowhere cross where their carriers do, if that lies on both, by the
   *  same exact crossing the thorough way finds; two lines nearer parallel meet nowhere where they
   *  do not cross, decided exactly, as the nearest they come is then at an end. Two nearer parallel
   *  that cross on both, every end farther than twice the tolerance from the other, meet once
   *  there, as Shallowly tells. All else is left to the thorough way.
   *
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @param  into      Where what they meet in goes, empty on entry.
   *  @return           Whether it was decided so; else the thorough way must decide.
   */
  auto Straightly(const Curves &curves, const Tolerance &tolerance, Intersections &into) -> bool;
} // namespace Geoxy::Met

#endif /* GEOXY_MET_LINEAR_HPP */
