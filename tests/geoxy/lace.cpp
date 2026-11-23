/***************************************************************************************************
 *  @file      tests/geoxy/lace.cpp
 *  @brief     Tests of the swift kernel of a spline, an ellipse or a Bézier curve and a polyline:
 *             every answer it is sure of the thorough way's to the bit, and the pairs it must
 *             leave.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "draws.hpp"
#include "swift/met.hpp"
#include "thorough.hpp"
#include "waves.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether two meetings are the same to the bit, hit for hit.
     *
     *  @param  one   One.
     *  @param  other The other.
     *  @return       Whether they are.
     */
    auto LaceSame(const Intersections &one, const Intersections &other) -> bool
    {
      const auto hit = [](const Hit &lhs, const Hit &rhs) -> bool
      {
        return lhs.At     == rhs.At
            && lhs.First  == rhs.First
            && lhs.Second == rhs.Second
            && lhs.Kind   == rhs.Kind;
      };
      return std::ranges::equal(one.Hits, other.Hits, hit)
          && one.Overlaps.size() == other.Overlaps.size()
          && one.Spent == other.Spent;
    }

    /**
     *  Two curves met by the swift kernels within the length Intersect would use, and where they
     *  are sure by the thorough way: whether they were sure, and whether the two agree to the bit.
     *
     *  @param  curves    The two curves.
     *  @param  tolerance The tolerance asked for.
     *  @return           Whether sure, and whether agreeing.
     */
    auto LaceBoth(const Curves &curves, const Tolerance &tolerance) -> std::array<bool, 2>
    {
      const std::optional<double> length = Reachable(curves, tolerance);
      if (!length)
      {
        return {false, true};
      }
      Intersections swiftly;
      const bool    sure =
          Swift::Met(curves, {.Length = *length, .Radians = tolerance.Radians}, swiftly);
      if (!sure)
      {
        return {false, true};
      }
      Intersections thorough;
      Thorough(curves, tolerance, thorough);
      return {true, LaceSame(swiftly, thorough)};
    }

    /**
     *  A polyline of a few vertices about the unit square, half its segments bulged, closed now and
     *  then.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The polyline.
     */
    auto Chained(Draws &draws) -> Polyline
    {
      Polyline          chain = {};
      const std::size_t count = 2 + (draws.Next() % 5);
      chain.Closed            = count >= 3 && draws.Next() % 3 == 0;
      for (std::size_t index = 0; index < count; ++index)
      {
        const double bulge = draws.Next() % 2 == 0 ? 0.0 : draws.Near(0.0, 1.5);
        chain.Vertices.push_back(
            {.Position = {.X = draws.Near(0.0, 1.0), .Y = draws.Near(0.0, 1.0)}, .Bulge = bulge});
      }
      return chain;
    }

    /**
     *  A curve scaled and moved, every point p to s p + m.
     *
     *  @param  curve The curve, a spline, an ellipse, a Bézier curve or a polyline.
     *  @param  scale s.
     *  @param  by    m.
     *  @return       The curve moved.
     */
    auto LaceMoved(const Curve &curve, double scale, Vector by) -> Curve
    {
      const auto moved = [scale, by](Point point) -> Point
      {
        return {.X = (point.X * scale) + by.X, .Y = (point.Y * scale) + by.Y};
      };
      if (const auto *spline = std::get_if<Spline>(&curve))
      {
        Spline out = *spline;
        std::ranges::transform(out.Poles, out.Poles.begin(), moved);
        return out;
      }
      if (const auto *ellipse = std::get_if<Ellipse>(&curve))
      {
        Ellipse out    = *ellipse;
        out.Map.Origin = moved(out.Map.Origin);
        out.Map.X      = out.Map.X * scale;
        out.Map.Y      = out.Map.Y * scale;
        return out;
      }
      if (const auto *quadratic = std::get_if<Bezier::Quadratic>(&curve))
      {
        return Bezier::Quadratic{.Start   = moved(quadratic->Start),
                                 .Control = moved(quadratic->Control),
                                 .End     = moved(quadratic->End)};
      }
      if (const auto *cubic = std::get_if<Bezier::Cubic>(&curve))
      {
        return Bezier::Cubic{.Start = moved(cubic->Start),
                             .Out   = moved(cubic->Out),
                             .In    = moved(cubic->In),
                             .End   = moved(cubic->End)};
      }
      Polyline out = std::get<Polyline>(curve);
      for (Vertex &vertex : out.Vertices)
      {
        vertex.Position = moved(vertex.Position);
      }
      return out;
    }

    /**
     *  An ellipse about the unit square, its axes at any angle and of any lengths up to one and a
     *  half, whole now and then, or a quadratic or cubic Bézier curve of control points about it.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The curve.
     */
    auto Laceworked(Draws &draws) -> Curve
    {
      const auto point = [&draws]() -> Point
      {
        const double across = draws.Near(0.0, 1.0);
        return {.X = across, .Y = draws.Near(0.0, 1.0)};
      };
      const auto kind = draws.Next() % 3;
      if (kind == 0)
      {
        const Point  center = {.X = draws.Near(0.0, 0.5), .Y = draws.Near(0.0, 0.5)};
        const Vector major  = {.X = draws.Near(0.0, 1.0), .Y = draws.Near(0.0, 1.0)};
        const Vector minor  = {.X = draws.Near(0.0, 1.0), .Y = draws.Near(0.0, 1.0)};
        const bool   whole  = draws.Next() % 4 == 0;
        const double start  = draws.Near(0.0, 180.0);
        return Ellipse{.Map   = {.X = major, .Y = minor, .Origin = center},
                       .Start = Degrees(start),
                       .Sweep = Degrees(whole ? 360.0 : draws.Near(0.0, 350.0))};
      }
      if (kind == 1)
      {
        const Point start   = point();
        const Point control = point();
        return Bezier::Quadratic{.Start = start, .Control = control, .End = point()};
      }
      const Point start  = point();
      const Point out    = point();
      const Point inward = point();
      return Bezier::Cubic{.Start = start, .Out = out, .In = inward, .End = point()};
    }
  } // namespace

  TEST_CASE("lace: ellipses and Bézier curves with polylines as met thoroughly, to the bit")
  {
    const int rounds = Rounds(20000);
    Draws     draws{0x0BE1ACEU};
    int       sure   = 0;
    int       differ = 0;
    for (int round = 0; round < rounds; ++round)
    {
      /**
       *  Scaled by 2^-40 to 2^40, or moved 1e8 off, the tolerance then no finer than 1e-7, or 1e6
       *  off, or as made.
       */
      constexpr std::array<Vector, 3> Moves = {
          Vector{.X = 0.0, .Y = 0.0}, Vector{.X = 1e8, .Y = -3e7}, Vector{.X = 1e6, .Y = -3e5}};
      constexpr std::array<double, 3> Fines  = {1e-10, 1e-7, 1e-10};
      const auto                      moving = draws.Next() % 6;
      const std::size_t               move   = moving < 3 ? moving : 0;
      const double scale  = std::ldexp(moving == 0 ? 1.37 : 1.0,
                                       moving == 0 ? static_cast<int>(draws.Next() % 80) - 40 : 0);
      const Curve  one    = LaceMoved(Laceworked(draws), scale, Moves.at(move));
      const Curve  other  = LaceMoved(Chained(draws), scale, Moves.at(move));
      const double near   = draws.Next() % 4 == 0 ? 1e-4 : Fines.at(move);
      const bool   swap   = draws.Next() % 2 == 0;
      const Curves curves = {swap ? other : one, swap ? one : other};
      const auto [certain, agrees]  = LaceBoth(curves, {.Length = near * scale, .Radians = 1e-10});
      sure                         += static_cast<int>(certain);
      differ                       += static_cast<int>(!agrees);
      if (!agrees)
      {
        INFO("round " << round);
        CHECK(agrees);
      }
    }
    MESSAGE("sure of " << sure << " of " << rounds);
    CHECK(differ == 0);
    CHECK(sure > rounds / 2);
  }

  TEST_CASE("lace: an ellipse, a Bézier curve and a spline with bulges far from nothing")
  {
    /**
     *  The comparison's quadratic across a polyline at 1e8 (case 3935) and a cubic (4324), each
     *  crossing two bulges: met about the bulges' rounded centers, a step of the doubles there from
     *  the true ones, the kernel's equation was off by that much and no crossing's place settled;
     *  met about the true centers in double-double, both are sure.
     */
    const Tolerance tolerance = {.Length = 0x1p-21, .Radians = 0x1p-22};
    const Curve     quadratic =
        Bezier::Quadratic{.Start   = {.X = -0x1.7d783fe145f2bp+26, .Y = 0x1.56bb588d50c14p+24},
                          .Control = {.X = -0x1.7d78403418560p+26, .Y = 0x1.56bb599887187p+24},
                          .End     = {.X = -0x1.7d784021f711dp+26, .Y = 0x1.56bb5929163a4p+24}};
    const Curve across = Polyline{
        .Vertices = {{.Position = {.X = -0x1.7d78401a75402p+26, .Y = 0x1.56bb587c23b29p+24},
                      .Bulge    = 0x1.7b82fe9ba7c98p-2},
                     {.Position = {.X = -0x1.7d783fe179d0fp+26, .Y = 0x1.56bb5974907a0p+24},
                      .Bulge    = 0.0},
                     {.Position = {.X = -0x1.7d783fe197d54p+26, .Y = 0x1.56bb5952305b4p+24},
                      .Bulge    = -0x1.226f07df65807p-1},
                     {.Position = {.X = -0x1.7d784035f03c4p+26, .Y = 0x1.56bb57c65d262p+24},
                      .Bulge    = 0.0}},
        .Closed   = false};
    CHECK(LaceBoth({quadratic, across}, tolerance) == std::array{true, true});
    CHECK(LaceBoth({across, quadratic}, tolerance) == std::array{true, true});
    const Curve cubic =
        Bezier::Cubic{.Start = {.X = 0x1.7d783fd647ae3p+26, .Y = -0x1.957f0a8aee9d8p+21},
                      .Out   = {.X = 0x1.7d783fd0deea2p+26, .Y = -0x1.957f080f19ed7p+21},
                      .In    = {.X = 0x1.7d78403079720p+26, .Y = -0x1.957f0fe764eddp+21},
                      .End   = {.X = 0x1.7d784025fb174p+26, .Y = -0x1.957f06a3ac454p+21}};
    const Curve bowed = Polyline{
        .Vertices = {{.Position = {.X = 0x1.7d783fddf68bap+26, .Y = -0x1.957f0ee529e7bp+21},
                      .Bulge    = 0.0},
                     {.Position = {.X = 0x1.7d78401bcb686p+26, .Y = -0x1.957f081cba2d7p+21},
                      .Bulge    = -0x1.c6567050d9e9cp-2},
                     {.Position = {.X = 0x1.7d783fd8dc1e9p+26, .Y = -0x1.957f07778fa2ap+21},
                      .Bulge    = -0x1.9366294099722p-1},
                     {.Position = {.X = 0x1.7d783fef6d6acp+26, .Y = -0x1.957f0e32a0fa8p+21},
                      .Bulge    = 0x1.a4b8be4680de0p-3}},
        .Closed   = false};
    CHECK(LaceBoth({cubic, bowed}, tolerance) == std::array{true, true});
    CHECK(LaceBoth({bowed, cubic}, tolerance) == std::array{true, true});
    /**
     *  An ellipse whose quarters lie but one near a bulge, and a cubic far from a bulge it does not
     *  reach, across a straight segment: each sure, the thorough way's.
     */
    const Tolerance fine = {.Length = 1e-9, .Radians = 1e-10};
    const Curve     oval = Ellipse{
        .Map   = {.X = {.X = 2.0, .Y = 0.0}, .Y = {.X = 0.0, .Y = 1.0}, .Origin = {.X = 0, .Y = 0}},
        .Start = Degrees(0.0),
        .Sweep = Degrees(360.0)};
    const Curve arch = Polyline{.Vertices = {{.Position = {.X = 1.5, .Y = 0.3}, .Bulge = 0.4},
                                             {.Position = {.X = 2.5, .Y = 0.3}, .Bulge = 0.0}},
                                .Closed   = false};
    CHECK(LaceBoth({oval, arch}, fine) == std::array{true, true});
    CHECK(LaceBoth({arch, oval}, fine) == std::array{true, true});
    const Curve wave  = Bezier::Cubic{.Start = {.X = -1, .Y = -1},
                                      .Out   = {.X = -0.5, .Y = 2},
                                      .In    = {.X = 0.5, .Y = -2},
                                      .End   = {.X = 1, .Y = 1}};
    const Curve apart = Polyline{.Vertices = {{.Position = {.X = -2, .Y = 0.1}, .Bulge = 0.0},
                                              {.Position = {.X = 2, .Y = 0.1}, .Bulge = 0.0},
                                              {.Position = {.X = 4, .Y = 0.1}, .Bulge = 1.0},
                                              {.Position = {.X = 6, .Y = 0.1}, .Bulge = 0.0}},
                                 .Closed   = false};
    CHECK(LaceBoth({wave, apart}, fine) == std::array{true, true});
    /**
     *  An ellipse and a Bézier curve of numbers past the kernels' range are left.
     */
    const Curve wild = Ellipse{
        .Map = {.X = {.X = 1e250, .Y = 0.0}, .Y = {.X = 0.0, .Y = 1.0}, .Origin = {.X = 0, .Y = 0}},
        .Start = Degrees(0.0),
        .Sweep = Degrees(360.0)};
    const Curve tall = Bezier::Quadratic{
        .Start = {.X = -1, .Y = -1}, .Control = {.X = 0, .Y = 1e250}, .End = {.X = 1, .Y = -1}};
    Intersections met;
    CHECK_FALSE(Swift::Met({wild, apart}, fine, met));
    CHECK_FALSE(Swift::Met({tall, apart}, fine, met));
  }

  TEST_CASE("lace: splines with polylines as met thoroughly, to the bit")
  {
    const int rounds = Rounds(20000);
    Draws     draws{0x1ACEDU};
    int       sure   = 0;
    int       differ = 0;
    for (int round = 0; round < rounds; ++round)
    {
      Spline spline = Tests::Splined(draws, 2, 5);
      if (draws.Next() % 3 == 0)
      {
        spline.Weights.resize(spline.Poles.size());
        for (double &weight : spline.Weights)
        {
          const double fraction = 1.0 + (0.5 * draws.Near(0.0, 1.0));
          weight                = std::ldexp(fraction, static_cast<int>(draws.Next() % 3) - 1);
        }
      }
      const auto   moving = draws.Next() % 6;
      const double scale  = std::ldexp(moving == 0 ? 1.37 : 1.0,
                                       moving == 0 ? static_cast<int>(draws.Next() % 80) - 40 : 0);
      const Vector by     = {.X = moving == 1 ? 1e6 : 0.0, .Y = moving == 1 ? -3e5 : 0.0};
      const Curve  one    = LaceMoved(spline, scale, by);
      const Curve  other  = LaceMoved(Chained(draws), scale, by);
      const double near   = draws.Next() % 4 == 0 ? 1e-4 : 1e-10;
      const bool   swap   = draws.Next() % 2 == 0;
      const Curves curves = {swap ? other : one, swap ? one : other};
      const auto [certain, agrees]  = LaceBoth(curves, {.Length = near * scale, .Radians = 1e-10});
      sure                         += static_cast<int>(certain);
      differ                       += static_cast<int>(!agrees);
      if (!agrees)
      {
        INFO("round " << round);
        CHECK(agrees);
      }
    }
    MESSAGE("sure of " << sure << " of " << rounds);
    CHECK(differ == 0);
    CHECK(sure > rounds / 4);
  }

  TEST_CASE("lace: a vertex on the spline, a segment of no length, and pairs it leaves")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-10};
    const Curve     spline = Spline{.Degree = 2,
                                    .Knots  = {0, 0, 0, 1, 1, 1},
                                    .Poles = {{.X = 0, .Y = 0}, {.X = 1, .Y = 2}, {.X = 2, .Y = 0}},
                                    .Weights = {}};
    const Point     top    = At(spline, 0.5);
    const auto      sure   = [&spline](const Polyline &chain, const Tolerance &within) -> bool
    {
      const Curve   other = chain;
      Intersections met;
      return Swift::Met({spline, other}, within, met);
    };
    /**
     *  A straight segment across, a bulge across, both crossing once, far from each other.
     */
    const Curve across = Polyline{.Vertices = {{.Position = {.X = -1, .Y = 0.5}, .Bulge = 0.0},
                                               {.Position = {.X = 0.6, .Y = 0.5}, .Bulge = 0.4},
                                               {.Position = {.X = 1.4, .Y = 0.5}, .Bulge = 0.0},
                                               {.Position = {.X = 1.4, .Y = 0.5}, .Bulge = 0.0},
                                               {.Position = {.X = 3, .Y = 0.5}, .Bulge = 0.0}},
                                  .Closed   = false};
    CHECK(LaceBoth({spline, across}, tolerance) == std::array{true, true});
    CHECK(LaceBoth({across, spline}, tolerance) == std::array{true, true});
    /**
     *  A vertex on the spline, its segments a line and a bulge: neither is sure.
     */
    const Polyline touching = {.Vertices = {{.Position = {.X = top.X, .Y = 3}, .Bulge = 0.0},
                                            {.Position = top, .Bulge = 0.0},
                                            {.Position = {.X = 3, .Y = 3}, .Bulge = 0.0}},
                               .Closed   = false};
    CHECK_FALSE(sure(touching, tolerance));
    const Polyline bowing = {
        .Vertices = {{.Position = top, .Bulge = 0.5}, {.Position = {.X = 3, .Y = 3}, .Bulge = 0.0}},
        .Closed   = false};
    CHECK_FALSE(sure(bowing, tolerance));
    /**
     *  A polyline far off, all its segments passed; one whose segments have no length; one of a
     *  segment too short to be a line; a bulge past the range the kernels take; a lone vertex.
     */
    const Polyline far      = {.Vertices = {{.Position = {.X = 50, .Y = 50}, .Bulge = 0.0},
                                            {.Position = {.X = 60, .Y = 50}, .Bulge = 0.0}},
                               .Closed   = false};
    const Curve    outlying = far;
    Intersections  none;
    CHECK(Swift::Met({spline, outlying}, tolerance, none));
    CHECK(none.Hits.empty());
    const Polyline point = {
        .Vertices = {{.Position = top, .Bulge = 0.0}, {.Position = top, .Bulge = 0.0}},
        .Closed   = false};
    CHECK_FALSE(sure(point, tolerance));
    const Polyline sliver = {
        .Vertices = {{.Position = {.X = 0x1p-300, .Y = 0x1p-300}, .Bulge = 0.0},
                     {.Position = {.X = 0x1p-299, .Y = 0x1p-300}, .Bulge = 0.0}},
        .Closed   = false};
    CHECK_FALSE(sure(sliver, {.Length = 0x1p-1070, .Radians = 1e-10}));
    const Polyline huge = {.Vertices = {{.Position = {.X = 0, .Y = 0.5}, .Bulge = 1e250},
                                        {.Position = {.X = 2, .Y = 0.5}, .Bulge = 0.0}},
                           .Closed   = false};
    CHECK_FALSE(sure(huge, tolerance));
    const Polyline lone = {.Vertices = {{.Position = top, .Bulge = 0.0}}, .Closed = false};
    CHECK_FALSE(sure(lone, tolerance));
    /**
     *  A bulge on a chord so short its radius rounds to nothing; three crossings, the spline
     *  halfway between the outer two on the polyline at the middle one.
     */
    const Polyline tiny = {.Vertices = {{.Position = {.X = 0x1p-700, .Y = 0.5}, .Bulge = 0.5},
                                        {.Position = {.X = 0x1p-699, .Y = 0.5}, .Bulge = 0.0}},
                           .Closed   = false};
    CHECK_FALSE(sure(tiny, tolerance));
    const Polyline flat = {.Vertices = {{.Position = {.X = -1e60, .Y = 0.5}, .Bulge = 1.5e-6},
                                        {.Position = {.X = 1e60, .Y = 0.5}, .Bulge = 0.0}},
                           .Closed   = false};
    CHECK_FALSE(sure(flat, tolerance));
    const Polyline thrice = {
        .Vertices = {{.Position = {.X = -1, .Y = 0.5}, .Bulge = 0.0},
                     {.Position = {.X = 1, .Y = 0.5}, .Bulge = 0.0},
                     {.Position = {.X = 1, .Y = 3}, .Bulge = 0.0},
                     {.Position = {.X = 1.0 + std::numbers::sqrt2, .Y = -2}, .Bulge = 0.0}},
        .Closed   = false};
    CHECK_FALSE(sure(thrice, tolerance));
    /**
     *  A spline of no poles, and one with a pole out of the kernels' range.
     */
    const Curve   empty   = Spline{};
    const Curve   across2 = thrice;
    Intersections met;
    CHECK_FALSE(Swift::Met({empty, across2}, tolerance, met));
    const Curve wild = Spline{.Degree  = 2,
                              .Knots   = {0, 0, 0, 1, 1, 1},
                              .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1e250}, {.X = 2, .Y = 0}},
                              .Weights = {}};
    CHECK_FALSE(Swift::Met({wild, across2}, tolerance, met));
  }
} // namespace Geoxy
