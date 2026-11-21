/***************************************************************************************************
 *  @file      tests/bench/intersect.cpp
 *  @brief     Benchmarks of Intersect, pair kind by pair kind, of curves as given and made ready.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/intersect.hpp"

#include "../geoxy/draws.hpp"
#include "shapes.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/tolerance.hpp"

#include <benchmark/benchmark.h>
#include <cstddef>
#include <cstdint>
#include <vector>

/**
 *  Each benchmark meets one pair a call, going round 256 pairs drawn from a fixed seed, so its time
 *  is that of one Intersect. Its Hits counter is how many hits a call finds, on average.
 *  ByIntersect meets the curves as given, a meeting returned; ByPrepared meets them made ready
 *  beforehand, into a meeting kept, as an application meeting a drawing's curves pair after pair
 *  does.
 */

namespace Geoxy::Bench
{
  namespace
  {
    struct Pairs
    {
      std::vector<Curve> First;
      std::vector<Curve> Second;
    };

    template <typename FirstT, typename SecondT>
    auto Paired(std::uint64_t seed, FirstT first, SecondT second) -> Pairs
    {
      Draws draws{seed};
      Pairs pairs;
      for (std::size_t index = 0; index < Pool; ++index)
      {
        pairs.First.push_back(first(draws));
        pairs.Second.push_back(second(draws));
      }
      return pairs;
    }

    auto LineLine() -> Pairs
    {
      return Paired(101, AnyLine, AnyLine);
    }

    auto LineArc() -> Pairs
    {
      return Paired(102, AnyLine, AnyArc);
    }

    auto ArcArc() -> Pairs
    {
      return Paired(103, AnyArc, AnyArc);
    }

    auto LineEllipse() -> Pairs
    {
      return Paired(104, AnyLine, AnyEllipse);
    }

    auto LineCubic() -> Pairs
    {
      return Paired(105, AnyLine, AnyCubic);
    }

    auto QuadraticQuadratic() -> Pairs
    {
      return Paired(110, AnyQuadratic, AnyQuadratic);
    }

    auto QuadraticCubic() -> Pairs
    {
      return Paired(111, AnyQuadratic, AnyCubic);
    }

    auto CubicCubic() -> Pairs
    {
      return Paired(106, AnyCubic, AnyCubic);
    }

    auto EllipseEllipse() -> Pairs
    {
      return Paired(107, AnyEllipse, AnyEllipse);
    }

    auto ArcEllipse() -> Pairs
    {
      return Paired(112, AnyArc, AnyEllipse);
    }

    auto ArcCubic() -> Pairs
    {
      return Paired(113, AnyArc, AnyCubic);
    }

    auto EllipseQuadratic() -> Pairs
    {
      return Paired(114, AnyEllipse, AnyQuadratic);
    }

    auto EllipseCubic() -> Pairs
    {
      return Paired(115, AnyEllipse, AnyCubic);
    }

    auto SplineSpline() -> Pairs
    {
      return Paired(108, AnySpline, AnySpline);
    }

    /**
     *  Lines across a zigzag of so many vertices, each crossing every segment.
     */
    template <std::size_t VerticesT>
    auto LinePolyline() -> Pairs
    {
      return {.First = Drawn(109, Across), .Second = std::vector<Curve>(Pool, Zigzag(VerticesT))};
    }

    /**
     *  Hits counted over the calls, as an average a call.
     */
    auto Hits(std::size_t hits) -> benchmark::Counter
    {
      return {static_cast<double>(hits), benchmark::Counter::kAvgIterations};
    }

    /**
     *  Intersect of curves as given.
     */
    void ByIntersect(benchmark::State &state, auto (*make)()->Pairs)
    {
      const Pairs     pairs     = make();
      const Tolerance tolerance = Within();

      std::size_t index = 0;
      std::size_t hits  = 0;
      for ([[maybe_unused]] auto iteration : state)
      {
        Intersections meeting = Intersect(pairs.First[index], pairs.Second[index], tolerance);
        benchmark::DoNotOptimize(meeting);
        hits  += meeting.Hits.size();
        index  = (index + 1) % Pool;
      }
      state.counters["Hits"] = Hits(hits);
    }

    /**
     *  Intersect of curves made ready, into a meeting kept.
     */
    void ByPrepared(benchmark::State &state, auto (*make)()->Pairs)
    {
      const Pairs     pairs     = make();
      const Tolerance tolerance = Within();

      std::vector<Prepared> first;
      std::vector<Prepared> second;
      first.reserve(Pool);
      second.reserve(Pool);
      for (std::size_t index = 0; index < Pool; ++index)
      {
        first.emplace_back(pairs.First[index]);
        second.emplace_back(pairs.Second[index]);
      }

      Intersections meeting;
      std::size_t   index = 0;
      std::size_t   hits  = 0;
      for ([[maybe_unused]] auto iteration : state)
      {
        Intersect(Prepared::Pair{first[index], second[index]}, tolerance, meeting);
        benchmark::DoNotOptimize(meeting);
        hits  += meeting.Hits.size();
        index  = (index + 1) % Pool;
      }
      state.counters["Hits"] = Hits(hits);
    }
  } // namespace

  BENCHMARK_CAPTURE(ByIntersect, LineLine, LineLine);
  BENCHMARK_CAPTURE(ByPrepared, LineLine, LineLine);
  BENCHMARK_CAPTURE(ByIntersect, LineArc, LineArc);
  BENCHMARK_CAPTURE(ByPrepared, LineArc, LineArc);
  BENCHMARK_CAPTURE(ByIntersect, ArcArc, ArcArc);
  BENCHMARK_CAPTURE(ByPrepared, ArcArc, ArcArc);
  BENCHMARK_CAPTURE(ByIntersect, LineEllipse, LineEllipse);
  BENCHMARK_CAPTURE(ByPrepared, LineEllipse, LineEllipse);
  BENCHMARK_CAPTURE(ByIntersect, LineCubic, LineCubic);
  BENCHMARK_CAPTURE(ByPrepared, LineCubic, LineCubic);
  BENCHMARK_CAPTURE(ByIntersect, QuadraticQuadratic, QuadraticQuadratic);
  BENCHMARK_CAPTURE(ByPrepared, QuadraticQuadratic, QuadraticQuadratic);
  BENCHMARK_CAPTURE(ByIntersect, QuadraticCubic, QuadraticCubic);
  BENCHMARK_CAPTURE(ByPrepared, QuadraticCubic, QuadraticCubic);
  BENCHMARK_CAPTURE(ByIntersect, CubicCubic, CubicCubic);
  BENCHMARK_CAPTURE(ByPrepared, CubicCubic, CubicCubic);
  BENCHMARK_CAPTURE(ByIntersect, EllipseEllipse, EllipseEllipse);
  BENCHMARK_CAPTURE(ByPrepared, EllipseEllipse, EllipseEllipse);
  BENCHMARK_CAPTURE(ByIntersect, ArcEllipse, ArcEllipse);
  BENCHMARK_CAPTURE(ByPrepared, ArcEllipse, ArcEllipse);
  BENCHMARK_CAPTURE(ByIntersect, ArcCubic, ArcCubic);
  BENCHMARK_CAPTURE(ByPrepared, ArcCubic, ArcCubic);
  BENCHMARK_CAPTURE(ByIntersect, EllipseQuadratic, EllipseQuadratic);
  BENCHMARK_CAPTURE(ByPrepared, EllipseQuadratic, EllipseQuadratic);
  BENCHMARK_CAPTURE(ByIntersect, EllipseCubic, EllipseCubic);
  BENCHMARK_CAPTURE(ByPrepared, EllipseCubic, EllipseCubic);
  BENCHMARK_CAPTURE(ByIntersect, SplineSpline, SplineSpline);
  BENCHMARK_CAPTURE(ByPrepared, SplineSpline, SplineSpline);
  BENCHMARK_CAPTURE(ByIntersect, LinePolyline8, LinePolyline<8>);
  BENCHMARK_CAPTURE(ByPrepared, LinePolyline8, LinePolyline<8>);
  BENCHMARK_CAPTURE(ByIntersect, LinePolyline64, LinePolyline<64>);
  BENCHMARK_CAPTURE(ByPrepared, LinePolyline64, LinePolyline<64>);
  BENCHMARK_CAPTURE(ByIntersect, LinePolyline512, LinePolyline<512>);
  BENCHMARK_CAPTURE(ByPrepared, LinePolyline512, LinePolyline<512>);
} // namespace Geoxy::Bench
