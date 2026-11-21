/***************************************************************************************************
 *  @file      source/elliptic.hpp
 *  @brief     Arcs as the ellipses they are, and where an ellipse passes an angle of its circle, so
 *             that what works on ellipses works on arcs.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_ELLIPTIC_HPP
#define GEOXY_ELLIPTIC_HPP

#include "real.hpp"
#include "scratch.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/ellipse.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>

namespace Geoxy
{
  /**
   *  The most turns of an ellipse whose passes are sought.
   */
  inline constexpr double MostTurns = 65'536.0;

  /**
   *  A whole turn in radians, as Turn.InRadians() gives it, worked out once: the double nearest 2
   *  pi, twice the double nearest pi.
   */
  inline constexpr double Revolution = 2.0 * Pi;

  /**
   *  An arc as an ellipse: the unit circle scaled by its radius and moved to its center.
   *
   *  @param  arc The arc.
   *  @return     The ellipse, its parameter the arc's.
   */
  inline auto Elliptic(const Arc &arc) -> Ellipse
  {
    return {.Map   = {.X      = {.X = arc.Radius, .Y = 0.0},
                      .Y      = {.X = 0.0, .Y = arc.Radius},
                      .Origin = arc.Center},
            .Start = arc.Start,
            .Sweep = arc.Sweep};
  }

  /**
   *  Where an ellipse passes an angle of its unit circle, in every turn it sweeps.
   *
   *  @param  ellipse The ellipse.
   *  @param  angle   The angle, in radians, any turn of it.
   *  @return         Its parameters there, from 0 to 1, unordered; none for a sweep of nothing, and
   *                  those of its first 65,536 turns only.
   */
  inline auto Passes(const Ellipse &ellipse, double angle) -> List<double>
  {
    List<double> passes;
    const double sweep = ellipse.Sweep.InRadians();
    if (Zero(sweep))
    {
      return passes;
    }
    const double whole  = Revolution;
    const double turned = angle - ellipse.Start.InRadians();
    const double first  = std::ceil((std::min(0.0, sweep) - turned) / whole);
    const double last   = std::floor((std::max(0.0, sweep) - turned) / whole);
    const auto   count =
        static_cast<std::size_t>(std::max(0.0, std::min(last - first + 1.0, MostTurns)));
    for (std::size_t turn = 0; turn < count; ++turn)
    {
      const double times = first + static_cast<double>(turn);
      passes.push_back(std::clamp((turned + (times * whole)) / sweep, 0.0, 1.0));
    }
    return passes;
  }

  /**
   *  An ellipse's start and sweep in radians, as InRadians gives them, worked out once.
   */
  struct Radiance
  {
    double Start = 0.0; /**< Its start. */
    double Sweep = 0.0; /**< Its sweep. */
  };

  /**
   *  An ellipse's start and sweep in radians.
   *
   *  @param  ellipse The ellipse.
   *  @return         Its start and sweep.
   */
  inline auto Radiant(const Ellipse &ellipse) -> Radiance
  {
    return {.Start = ellipse.Start.InRadians(), .Sweep = ellipse.Sweep.InRadians()};
  }

  /**
   *  The first place an ellipse passes an angle of its unit circle: the least of Passes, worked as
   *  it is for the one turn that gives it, the first turn for a sweep counterclockwise and the last
   *  for one clockwise, as each step of the sum and quotient keeps the order of the turns; so with
   *  no list.
   *
   *  @param  radiance The ellipse's start and sweep in radians.
   *  @param  angle    The angle, in radians, any turn of it.
   *  @return          The least of its passes; nothing where there is none.
   */
  inline auto Earliest(Radiance radiance, double angle) -> std::optional<double>
  {
    const double sweep = radiance.Sweep;
    if (Zero(sweep))
    {
      return std::nullopt;
    }
    const double whole  = Revolution;
    const double turned = angle - radiance.Start;
    const double first  = std::ceil((std::min(0.0, sweep) - turned) / whole);
    const double last   = std::floor((std::max(0.0, sweep) - turned) / whole);
    const double count  = std::max(0.0, std::min(last - first + 1.0, MostTurns));
    if (!(count >= 1.0))
    {
      return std::nullopt;
    }
    const double times = sweep > 0.0 ? first : first + (count - 1.0);
    return std::clamp((turned + (times * whole)) / sweep, 0.0, 1.0);
  }

  /**
   *  The first place an ellipse passes an angle of its unit circle, as Earliest tells it from its
   *  start and sweep.
   *
   *  @param  ellipse The ellipse.
   *  @param  angle   The angle, in radians, any turn of it.
   *  @return         The least of its passes; nothing where there is none.
   */
  inline auto Earliest(const Ellipse &ellipse, double angle) -> std::optional<double>
  {
    return Earliest(Radiant(ellipse), angle);
  }

  /**
   *  Whether an ellipse sweeps round more than once, and so passes each angle on more than one
   *  turn.
   *
   *  @param  ellipse The ellipse.
   *  @return         Whether it does.
   */
  inline auto Coiled(const Ellipse &ellipse) -> bool
  {
    return std::abs(ellipse.Sweep.InRadians()) > Revolution;
  }

  /**
   *  Where an ellipse passes an angle of its unit circle nearest a parameter of its own: of Passes,
   *  the turn a place already near it is on, with no list. One swept round more than once passes
   *  each angle on every turn, and a hit found on one turn is to stay on it; one swept round no
   *  more than once is passed where Earliest tells. Its start and sweep are turned into radians
   *  once.
   *
   *  @param  ellipse The ellipse.
   *  @param  angle   The angle, in radians, any turn of it.
   *  @param  along   The parameter.
   *  @return         The pass nearest the parameter; nothing where there is none.
   */
  inline auto Abreast(const Ellipse &ellipse, double angle, double along) -> std::optional<double>
  {
    const Radiance radiance = Radiant(ellipse);
    if (!(std::abs(radiance.Sweep) > Revolution))
    {
      return Earliest(radiance, angle);
    }
    const double sweep  = radiance.Sweep;
    const double whole  = Revolution;
    const double turned = angle - radiance.Start;
    const double first  = std::ceil((std::min(0.0, sweep) - turned) / whole);
    const double last   = std::floor((std::max(0.0, sweep) - turned) / whole);
    const double times  = std::clamp(std::round(((along * sweep) - turned) / whole), first, last);
    return std::clamp((turned + (times * whole)) / sweep, 0.0, 1.0);
  }
} // namespace Geoxy

#endif /* GEOXY_ELLIPTIC_HPP */
