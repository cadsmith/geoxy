/***************************************************************************************************
 *  @file      source/spline.cpp
 *  @brief     Whether a spline is one.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/spline.hpp"

#include "basis.hpp"
#include "finite.hpp"

#include "geoxy/point.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether a weight is one a rational spline can hold.
     *
     *  @param  weight The weight.
     *  @return        Whether it is finite and positive.
     */
    auto Positive(double weight) -> bool
    {
      return std::isfinite(weight) && weight > 0.0;
    }

    /**
     *  Whether a spline's knots hold it in one piece: no knot within its domain as often as a
     *  degree and one, where the curve may jump from one pole to another, and none anywhere as
     *  often as a degree and two, which leaves a pole that nothing weighs.
     *
     *  @pre           The spline's counts agree, and its knots are finite and none falls.
     *  @param  spline The spline.
     *  @return        Whether they do.
     */
    auto Unbroken(const Spline &spline) -> bool
    {
      const auto                    degree = static_cast<std::size_t>(spline.Degree);
      const std::span<const double> knots  = spline.Knots;
      const double                  start  = knots[degree];
      const double                  end    = knots[spline.Poles.size()];
      for (std::size_t index = 0; index + degree < knots.size(); ++index)
      {
        const double knot     = knots[index];
        const bool   repeated = !(knot < knots[index + degree]);
        const bool   inside   = start < knot && knot < end;
        const bool   unheld =
            index + degree + 1 < knots.size() && !(knot < knots[index + degree + 1]);
        if ((repeated && inside) || unheld)
        {
          return false;
        }
      }
      return true;
    }

    /**
     *  Whether a power of two brings the weights about each span of a spline within the doubles:
     *  the one Leveling gives, which leaves the lightest a double of all its digits, also leaves
     *  every pole times its weight finite. Past that the spline's points are worked out from no
     *  number, whatever power of two its weights are taken at.
     *
     *  @pre           The spline's counts agree, its knots hold it in one piece, and its weights
     *                 are finite and positive.
     *  @param  spline The spline.
     *  @return        Whether they do.
     */
    auto Levelable(const Spline &spline) -> bool
    {
      if (spline.Weights.empty())
      {
        return true;
      }
      const auto degree = static_cast<std::size_t>(spline.Degree);
      for (std::size_t span = degree; span < spline.Poles.size(); ++span)
      {
        if (!(spline.Knots[span] < spline.Knots[span + 1]))
        {
          continue;
        }
        const int power = Leveling(spline, span);
        for (std::size_t index = span - degree; index <= span; ++index)
        {
          const double weight = std::ldexp(spline.Weights[index], power);
          const Point  pole   = spline.Poles[index];
          if (!std::isfinite(weight * std::max({1.0, std::abs(pole.X), std::abs(pole.Y)})))
          {
            return false;
          }
        }
      }
      return true;
    }
  } // namespace

  auto Valid(const Spline &spline) -> bool
  {
    if (spline.Degree < 1 || spline.Degree > MaxDegree)
    {
      return false;
    }
    const auto        degree = static_cast<std::size_t>(spline.Degree);
    const std::size_t poles  = spline.Poles.size();
    if (poles <= degree || spline.Knots.size() != poles + degree + 1)
    {
      return false;
    }
    const bool finite =
        std::ranges::all_of(spline.Knots, [](double knot) -> bool { return std::isfinite(knot); });
    const bool rising = std::ranges::is_sorted(spline.Knots);
    const bool spans  = spline.Knots[degree] < spline.Knots[poles];
    const bool placed =
        std::ranges::all_of(spline.Poles, [](Point point) -> bool { return Finite(point); });
    const bool weighted =
        spline.Weights.empty()
        || (spline.Weights.size() == poles && std::ranges::all_of(spline.Weights, Positive));
    return finite && rising && spans && placed && weighted && Unbroken(spline) && Levelable(spline);
  }
} // namespace Geoxy
