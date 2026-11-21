/***************************************************************************************************
 *  @file      tests/geoxy/waves.hpp
 *  @brief     Polynomial splines drawn at random: clamped or not, their knots even or not and some
 *             repeated, starting anywhere and spaced at any scale, their poles waving from left to
 *             right or looping anywhere, two of them one now and then.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TESTS_WAVES_HPP
#define GEOXY_TESTS_WAVES_HPP

#include "draws.hpp"

#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace Geoxy::Tests
{
  /**
   *  A spline's knots: clamped or not, even or not, some repeated within the domain as often as the
   *  degree, starting anywhere and spaced at any scale.
   *
   *  @param  draws  The draws.
   *  @param  degree The degree.
   *  @param  poles  How many poles.
   *  @return        The knots.
   */
  inline auto Knotted(Draws &draws, int degree, std::size_t poles) -> std::vector<double>
  {
    const auto          order   = static_cast<std::size_t>(degree);
    const std::size_t   count   = poles + order + 1;
    const bool          clamped = draws.Next() % 3 != 0;
    const bool          even    = draws.Next() % 2 == 0;
    const double        start   = draws.Next() % 3 == 0 ? draws.Near(0.0, 1e3) : 0.0;
    const double        scale = draws.Next() % 3 == 0 ? std::pow(10.0, draws.Near(0.0, 2.0)) : 1.0;
    std::vector<double> knots = {start};
    double              knot  = start;
    std::size_t         run   = 1;
    for (std::size_t index = 1; index < count; ++index)
    {
      const bool head   = clamped && index <= order;
      const bool tail   = clamped && index + order + 1 > count;
      const bool again  = !even && index > order && index < poles && run < order;
      const bool same   = head || tail || (again && draws.Next() % 5 == 0);
      knot             += same ? 0.0 : scale * (even ? 1.0 : 0.25 + std::abs(draws.Near(0.0, 1.0)));
      run               = same ? run + 1 : 1;
      knots.push_back(knot);
    }
    return knots;
  }

  /**
   *  A polynomial spline of a degree between two: its poles from left to right, up and down, as the
   *  benchmark draws them, or anywhere, looping, now and then two of them one; its weights none or
   *  all one.
   *
   *  @param  draws  The draws.
   *  @param  lowest The least degree.
   *  @param  most   How many degrees above it at most.
   *  @return        The spline, valid.
   */
  inline auto Splined(Draws &draws, int lowest, int most) -> Spline
  {
    const int degree =
        lowest + static_cast<int>(draws.Next() % static_cast<std::uint64_t>(most + 1));
    const std::size_t poles  = static_cast<std::size_t>(degree) + 1 + (draws.Next() % 7);
    const bool        wavy   = draws.Next() % 2 == 0;
    Spline            spline = {.Degree = degree, .Knots = {}, .Poles = {}, .Weights = {}};
    for (std::size_t index = 0; index < poles; ++index)
    {
      const double along = static_cast<double>(index) / static_cast<double>(poles - 1);
      spline.Poles.push_back(
          wavy ? Point{.X = along + draws.Near(0.0, 0.05), .Y = draws.Near(0.0, 0.4)}
               : Point{.X = draws.Near(0.0, 1.0), .Y = draws.Near(0.0, 1.0)});
      if (index > 0 && draws.Next() % 12 == 0)
      {
        spline.Poles.back() = spline.Poles[index - 1];
      }
    }
    spline.Knots = Knotted(draws, degree, poles);
    if (draws.Next() % 4 == 0)
    {
      spline.Weights.assign(poles, draws.Next() % 2 == 0 ? 0.75 : 3.0);
    }
    return spline;
  }
} // namespace Geoxy::Tests

#endif /* GEOXY_TESTS_WAVES_HPP */
