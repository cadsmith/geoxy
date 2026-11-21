/***************************************************************************************************
 *  @file      source/straight.hpp
 *  @brief     Lines, rays and xlines as one straight run: a start, a way and a stretch of its
 *             parameter.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_STRAIGHT_HPP
#define GEOXY_STRAIGHT_HPP

#include "real.hpp"
#include "scratch.hpp"

#include "geoxy/domain.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <optional>

namespace Geoxy
{
  /**
   *  A straight run: where it starts, its way, and the stretch of its parameter it runs over.
   */
  struct Run
  {
    Point  Start = {}; /**< Where it is at 0.           */
    Vector Way   = {}; /**< How far it goes for each 1. */
    Range  Along = {}; /**< Its parameter's stretch.    */
  };

  /**
   *  A line as a run.
   *
   *  @param  line The line.
   *  @return      The run from its start to its end, over 0 to 1.
   */
  inline auto Running(const Line &line) -> Run
  {
    return {.Start = line.Start, .Way = line.End - line.Start, .Along = Domain(line)};
  }

  /**
   *  A ray as a run.
   *
   *  @param  ray The ray.
   *  @return     The run from its start, over 0 to infinity.
   */
  inline auto Running(const Ray &ray) -> Run
  {
    return {.Start = ray.Start, .Way = ray.Way, .Along = Domain(ray)};
  }

  /**
   *  An xline as a run.
   *
   *  @param  xline The xline.
   *  @return       The run through its point, over all of its parameter.
   */
  inline auto Running(const Xline &xline) -> Run
  {
    return {.Start = xline.Through, .Way = xline.Way, .Along = Domain(xline)};
  }

  /**
   *  Where a run comes nearest a point: the foot of the perpendicular from it.
   *
   *  @param  run   The run.
   *  @param  point The point.
   *  @return       The foot's parameter, if the run has a way and the foot lies on its stretch.
   */
  inline auto Footing(const Run &run, Point point) -> std::optional<double>
  {
    const double square = Dot(run.Way, run.Way);
    if (Zero(square))
    {
      return std::nullopt;
    }
    const double along = Dot(point - run.Start, run.Way) / square;
    if (run.Along.Start <= along && along <= run.Along.End)
    {
      return along;
    }
    return std::nullopt;
  }

  /**
   *  Where a run comes nearest a point, as a list: the foot, or none.
   *
   *  @param  run   The run.
   *  @param  point The point.
   *  @return       The foot's parameter, if the run has a way and the foot lies on its stretch.
   */
  inline auto Foot(const Run &run, Point point) -> List<double>
  {
    if (const std::optional<double> along = Footing(run, point))
    {
      return {*along};
    }
    return {};
  }
} // namespace Geoxy

#endif /* GEOXY_STRAIGHT_HPP */
