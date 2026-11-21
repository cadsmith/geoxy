/***************************************************************************************************
 *  @file      tests/geoxy/polylines.hpp
 *  @brief     Polylines drawn at random, near every edge a polyline has: vertices on a coarse grid
 *             or anywhere, given twice, segments straight, bulged slightly, bulged past a half turn
 *             or by a bulge just under where one turns straight, open and closed.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TESTS_POLYLINES_HPP
#define GEOXY_TESTS_POLYLINES_HPP

#include "draws.hpp"

#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"

#include <array>
#include <cstddef>

namespace Geoxy
{
  /**
   *  A bulge drawn at random: none half the time, else one of a few that matter or any.
   *
   *  @param  draws The draws.
   *  @return       The bulge.
   */
  inline auto Bulged(Draws &draws) -> double
  {
    constexpr std::array Kept = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, Straight, -Straight / 2};
    const std::size_t    pick = draws.Next() % (Kept.size() + 4);
    if (pick < Kept.size())
    {
      return Kept.at(pick);
    }
    if (pick == Kept.size())
    {
      return draws.Near(0, 1e-3);
    }
    if (pick == Kept.size() + 1)
    {
      return draws.Near(0, 40);
    }
    return draws.Near(0, 2);
  }

  /**
   *  A polyline drawn at random about the unit square: its vertices on a grid of quarters a third
   *  of the time, so that segments share lines and points lie on them exactly; one vertex in eight
   *  given twice; closed a third of the time.
   *
   *  @param  draws The draws.
   *  @param  count How many vertices, before any given twice.
   *  @return       The polyline.
   */
  inline auto Strung(Draws &draws, std::size_t count) -> Polyline
  {
    Polyline   polyline;
    const bool grid = draws.Next() % 3 == 0;
    for (std::size_t index = 0; index < count; ++index)
    {
      Point position = {.X = draws.Near(0.5, 0.75), .Y = draws.Near(0.5, 0.75)};
      if (grid)
      {
        position = {.X = (static_cast<double>(draws.Next() % 9) / 4) - 0.5,
                    .Y = (static_cast<double>(draws.Next() % 9) / 4) - 0.5};
      }
      polyline.Vertices.push_back({.Position = position, .Bulge = Bulged(draws)});
      if (draws.Next() % 8 == 0)
      {
        polyline.Vertices.push_back({.Position = position, .Bulge = Bulged(draws)});
      }
    }
    polyline.Closed = draws.Next() % 3 == 0;
    return polyline;
  }

  /**
   *  A polyline zigzagging left to right across the unit square, a third of its segments bulged,
   *  its heights a little off by draws.
   *
   *  @param  draws The draws.
   *  @param  count How many vertices, two at the least.
   *  @return       The polyline.
   */
  inline auto Zigzagged(Draws &draws, std::size_t count) -> Polyline
  {
    Polyline polyline;
    for (std::size_t index = 0; index < count; ++index)
    {
      const double along = static_cast<double>(index) / static_cast<double>(count - 1);
      const double high  = index % 2 == 0 ? 0.2 : 0.8;
      polyline.Vertices.push_back({.Position = {.X = along, .Y = high + draws.Near(0, 0.05)},
                                   .Bulge    = index % 3 == 0 ? draws.Near(0, 0.4) : 0.0});
    }
    return polyline;
  }
} // namespace Geoxy

#endif /* GEOXY_TESTS_POLYLINES_HPP */
