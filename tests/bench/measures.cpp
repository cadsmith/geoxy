/***************************************************************************************************
 *  @file      tests/bench/measures.cpp
 *  @brief     Benchmarks of Nearest, Length, Area, Extents and Flatten, curve kind by curve kind.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "../geoxy/draws.hpp"
#include "shapes.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/area.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/flatten.hpp"
#include "geoxy/length.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"

#include <benchmark/benchmark.h>
#include <cstddef>
#include <vector>

/**
 *  Each benchmark asks one curve a call, going round 256 curves of a kind drawn from a fixed seed,
 *  so its time is that of one call. Length and Area are asked as near as can be, Accuracy{0}, and
 *  Flatten to a thousandth of a unit, a ten thousandth of the drawing; its Points counter is how
 *  many points a curve flattens into, on average.
 */

namespace Geoxy::Bench
{
  namespace
  {
    auto Lines() -> std::vector<Curve>
    {
      return Drawn(201, AnyLine);
    }

    auto Rays() -> std::vector<Curve>
    {
      return Drawn(202, AnyRay);
    }

    auto Xlines() -> std::vector<Curve>
    {
      return Drawn(203, AnyXline);
    }

    auto Arcs() -> std::vector<Curve>
    {
      return Drawn(204, AnyArc);
    }

    auto Ellipses() -> std::vector<Curve>
    {
      return Drawn(205, AnyEllipse);
    }

    auto Polylines() -> std::vector<Curve>
    {
      return Drawn(206, AnyPolyline);
    }

    auto Splines() -> std::vector<Curve>
    {
      return Drawn(207, AnySpline);
    }

    auto Quadratics() -> std::vector<Curve>
    {
      return Drawn(208, AnyQuadratic);
    }

    auto Cubics() -> std::vector<Curve>
    {
      return Drawn(209, AnyCubic);
    }

    constexpr Accuracy Finest = {.Error = 0};
    constexpr Accuracy Chords = {.Error = 1e-3};

    auto Points() -> std::vector<Point>
    {
      Draws              draws{299};
      std::vector<Point> points;
      points.reserve(Pool);
      for (std::size_t index = 0; index < Pool; ++index)
      {
        points.push_back(Anywhere(draws));
      }
      return points;
    }

    /**
     *  Calls a query on each curve in turn.
     */
    template <typename QueryT>
    void Each(benchmark::State &state, const std::vector<Curve> &curves, QueryT query)
    {
      std::size_t index = 0;
      for ([[maybe_unused]] auto iteration : state)
      {
        auto answer = query(curves[index], index);
        benchmark::DoNotOptimize(answer);
        index = (index + 1) % Pool;
      }
    }

    /**
     *  The place on a curve nearest a point anywhere in the drawing.
     */
    void ByNearest(benchmark::State &state, auto (*make)()->std::vector<Curve>)
    {
      const std::vector<Point> points = Points();
      Each(state, make(), [&points](const Curve &curve, std::size_t index) -> Place
           { return Nearest(curve, points[index]); });
    }

    void ByLength(benchmark::State &state, auto (*make)()->std::vector<Curve>)
    {
      Each(state, make(), [](const Curve &curve, [[maybe_unused]] std::size_t index) -> double
           { return Length(curve, Finest); });
    }

    void ByArea(benchmark::State &state, auto (*make)()->std::vector<Curve>)
    {
      Each(state, make(), [](const Curve &curve, [[maybe_unused]] std::size_t index) -> double
           { return Area(curve, Finest); });
    }

    void ByExtents(benchmark::State &state, auto (*make)()->std::vector<Curve>)
    {
      Each(state, make(), [](const Curve &curve, [[maybe_unused]] std::size_t index) -> Box
           { return Extents(curve); });
    }

    void ByFlatten(benchmark::State &state, auto (*make)()->std::vector<Curve>)
    {
      const std::vector<Curve> curves = make();

      std::size_t points = 0;
      for (const Curve &curve : curves)
      {
        points += Flatten(curve, Chords).size();
      }
      Each(state, curves,
           [](const Curve &curve, [[maybe_unused]] std::size_t index) -> std::vector<Point>
           { return Flatten(curve, Chords); });
      state.counters["Points"] = static_cast<double>(points) / static_cast<double>(Pool);
    }
  } // namespace

  BENCHMARK_CAPTURE(ByNearest, Line, Lines);
  BENCHMARK_CAPTURE(ByNearest, Ray, Rays);
  BENCHMARK_CAPTURE(ByNearest, Xline, Xlines);
  BENCHMARK_CAPTURE(ByNearest, Arc, Arcs);
  BENCHMARK_CAPTURE(ByNearest, Ellipse, Ellipses);
  BENCHMARK_CAPTURE(ByNearest, Polyline, Polylines);
  BENCHMARK_CAPTURE(ByNearest, Spline, Splines);
  BENCHMARK_CAPTURE(ByNearest, Quadratic, Quadratics);
  BENCHMARK_CAPTURE(ByNearest, Cubic, Cubics);

  BENCHMARK_CAPTURE(ByLength, Line, Lines);
  BENCHMARK_CAPTURE(ByLength, Arc, Arcs);
  BENCHMARK_CAPTURE(ByLength, Ellipse, Ellipses);
  BENCHMARK_CAPTURE(ByLength, Polyline, Polylines);
  BENCHMARK_CAPTURE(ByLength, Spline, Splines);
  BENCHMARK_CAPTURE(ByLength, Quadratic, Quadratics);
  BENCHMARK_CAPTURE(ByLength, Cubic, Cubics);

  BENCHMARK_CAPTURE(ByArea, Line, Lines);
  BENCHMARK_CAPTURE(ByArea, Arc, Arcs);
  BENCHMARK_CAPTURE(ByArea, Ellipse, Ellipses);
  BENCHMARK_CAPTURE(ByArea, Polyline, Polylines);
  BENCHMARK_CAPTURE(ByArea, Spline, Splines);
  BENCHMARK_CAPTURE(ByArea, Quadratic, Quadratics);
  BENCHMARK_CAPTURE(ByArea, Cubic, Cubics);

  BENCHMARK_CAPTURE(ByExtents, Line, Lines);
  BENCHMARK_CAPTURE(ByExtents, Arc, Arcs);
  BENCHMARK_CAPTURE(ByExtents, Ellipse, Ellipses);
  BENCHMARK_CAPTURE(ByExtents, Polyline, Polylines);
  BENCHMARK_CAPTURE(ByExtents, Spline, Splines);
  BENCHMARK_CAPTURE(ByExtents, Quadratic, Quadratics);
  BENCHMARK_CAPTURE(ByExtents, Cubic, Cubics);

  BENCHMARK_CAPTURE(ByFlatten, Line, Lines);
  BENCHMARK_CAPTURE(ByFlatten, Arc, Arcs);
  BENCHMARK_CAPTURE(ByFlatten, Ellipse, Ellipses);
  BENCHMARK_CAPTURE(ByFlatten, Polyline, Polylines);
  BENCHMARK_CAPTURE(ByFlatten, Spline, Splines);
  BENCHMARK_CAPTURE(ByFlatten, Quadratic, Quadratics);
  BENCHMARK_CAPTURE(ByFlatten, Cubic, Cubics);
} // namespace Geoxy::Bench
