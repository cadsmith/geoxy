/***************************************************************************************************
 *  @file      tests/geoxy/crossings/self.cpp
 *  @brief     Tests of where a curve crosses itself: polylines worked by hand, a cubic's loop
 *             against its double point in closed form, splines, and curves that never cross.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "crossings/self.hpp"

#include "../curves.hpp"
#include "real.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/fit.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

namespace Geoxy::Crossings
{
  namespace
  {
    constexpr Tolerance Fine = {.Length = 1e-9, .Radians = 1e-9};

    /**
     *  A polyline through points, every bulge 0.
     *
     *  @param  points The points.
     *  @param  closed Whether it closes.
     *  @return        The polyline.
     */
    auto Through(const std::vector<Point> &points, bool closed) -> Polyline
    {
      Polyline polyline = {.Vertices = {}, .Closed = closed};
      for (const Point point : points)
      {
        polyline.Vertices.push_back({.Position = point, .Bulge = 0.0});
      }
      return polyline;
    }
  } // namespace

  TEST_CASE("self: a figure eight and a star")
  {
    for (const bool closed : {false, true})
    {
      const Polyline eight =
          Through({{.X = 0, .Y = 0}, {.X = 2, .Y = 2}, {.X = 2, .Y = 0}, {.X = 0, .Y = 2}}, closed);
      const Intersections met = Self(eight, Fine);
      CHECK(SelfIntersect(eight, Fine).Hits.size() == met.Hits.size());
      REQUIRE(met.Hits.size()        == 1);
      CHECK(met.Hits[0].Kind         == Contact::Cross);
      CHECK(met.Hits[0].First.Along  == doctest::Approx(0.5).epsilon(1e-15));
      CHECK(met.Hits[0].Second.Along == doctest::Approx(2.5).epsilon(1e-15));
      CHECK(met.Overlaps.empty());
      CHECK_FALSE(met.Spent);

      Polyline back = eight;
      std::ranges::reverse(back.Vertices);
      const Intersections reversed = Self(back, Fine);
      REQUIRE(reversed.Hits.size()       == 1);
      CHECK(reversed.Hits[0].First.Along == doctest::Approx(0.5).epsilon(1e-15));
    }
    std::vector<Point> star;
    for (int index = 0; index < 5; ++index)
    {
      const double angle = (std::numbers::pi / 2) + (index * 4 * std::numbers::pi / 5);
      star.push_back({.X = 5e6 + std::cos(angle), .Y = 5e6 + std::sin(angle)});
    }
    const Intersections pentagram = Self(Through(star, true), {.Length = 1e-8, .Radians = 1e-9});
    CHECK(pentagram.Hits.size() == 5);
    CHECK(std::ranges::all_of(
        pentagram.Hits, [](const Hit &hit) -> bool
        { return hit.Kind == Contact::Cross && hit.First.Along < hit.Second.Along; }));
  }

  TEST_CASE("self: a vertex passed through again")
  {
    const Polyline      crossing = Through({{.X = -1, .Y = 0},
                                            {.X = 0, .Y = 0},
                                            {.X = 1, .Y = 0},
                                            {.X = 1, .Y = 1},
                                            {.X = 0, .Y = 0},
                                            {.X = 0, .Y = -1}},
                                           false);
    const Intersections crossed  = Self(crossing, Fine);
    REQUIRE(crossed.Hits.size() == 1);
    CHECK(crossed.Hits[0].Kind  == Contact::Cross);
    CHECK(Exactly(crossed.Hits[0].First.Along, 1));
    CHECK(Exactly(crossed.Hits[0].Second.Along, 4));

    const Polyline      touching = Through({{.X = -1, .Y = 0},
                                            {.X = 0, .Y = 0},
                                            {.X = 1, .Y = 0},
                                            {.X = 1, .Y = 1},
                                            {.X = 0, .Y = 0},
                                            {.X = -1, .Y = 1}},
                                           false);
    const Intersections touched  = Self(touching, Fine);
    REQUIRE(touched.Hits.size() == 1);
    CHECK(touched.Hits[0].Kind  == Contact::Touch);

    const Polyline ended =
        Through({{.X = 0, .Y = 0}, {.X = 2, .Y = 0}, {.X = 2, .Y = 2}, {.X = 1, .Y = 0}}, false);
    const Intersections end = Self(ended, Fine);
    REQUIRE(end.Hits.size() == 1);
    CHECK(end.Hits[0].Kind  == Contact::End);

    const Polyline      twice = Through({{.X = 0, .Y = 0},
                                         {.X =  1, .Y =  1},
                                         {.X =  1, .Y = -1},
                                         {.X =  0, .Y =  0},
                                         {.X = -1, .Y =  1},
                                         {.X = -1, .Y = -1}},
                                        true);
    const Intersections again = Self(twice, Fine);
    REQUIRE(again.Hits.size() == 1);
    CHECK(again.Hits[0].Kind  == Contact::Cross);
    CHECK(Exactly(again.Hits[0].First.Along, 0));
    CHECK(Exactly(again.Hits[0].Second.Along, 3));

    const Polyline      hooked = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                               {.Position = {.X = 4, .Y = 0}, .Bulge = -1},
                                               {.Position = {.X = 2, .Y = 0}, .Bulge = 0}},
                                  .Closed   = false};
    const Intersections hook   = Self(hooked, Fine);
    REQUIRE(hook.Hits.size() == 1);
    CHECK(hook.Hits[0].Kind  == Contact::End);

    const Polyline started =
        Through({{.X = 2, .Y = 0}, {.X = 2, .Y = 2}, {.X = 0, .Y = 0}, {.X = 4, .Y = 0}}, false);
    const Intersections start = Self(started, Fine);
    REQUIRE(start.Hits.size() == 1);
    CHECK(start.Hits[0].Kind  == Contact::End);

    const Polyline      bent  = Through({{.X = -1, .Y = -1},
                                         {.X = 0, .Y =  0},
                                         {.X = 1, .Y = -1},
                                         {.X = 1, .Y =  1},
                                         {.X = 0, .Y =  0},
                                         {.X = 0, .Y = -2}},
                                        false);
    const Intersections under = Self(bent, Fine);
    REQUIRE(under.Hits.size() == 1);
    CHECK(under.Hits[0].Kind  == Contact::Cross);
    const Polyline      aside = Through({{.X = -1, .Y = -1},
                                         {.X = 0, .Y =  0},
                                         {.X = 1, .Y = -1},
                                         {.X = 1, .Y =  1},
                                         {.X = 0, .Y =  0},
                                         {.X = -1, .Y = 1}},
                                        false);
    const Intersections over  = Self(aside, Fine);
    REQUIRE(over.Hits.size() == 1);
    CHECK(over.Hits[0].Kind  == Contact::Touch);

    const Polyline      along  = Through({{.X = -1, .Y = 0},
                                          {.X = 0, .Y = 0},
                                          {.X = 1, .Y = 0},
                                          {.X = 1, .Y = 1},
                                          {.X = 0, .Y = 0},
                                          {.X = 1, .Y = 0}},
                                         false);
    const Intersections shared = Self(along, Fine);
    CHECK_FALSE(shared.Overlaps.empty());
    CHECK(
        std::ranges::any_of(shared.Hits, [](const Hit &hit) -> bool
                            { return Exactly(hit.First.Along, 1) && hit.Kind == Contact::Touch; }));

    const Polyline back = Through({{.X = 0, .Y = 0}, {.X = 4, .Y = 0}, {.X = 2, .Y = 0}}, false);
    CHECK(Self(back, Fine).Overlaps.size() == 1);
  }

  TEST_CASE("self: a cubic's loop and a spline's")
  {
    const Bezier::Cubic loop  = {.Start = {.X = 0, .Y = 0},
                                 .Out   = {.X = 2, .Y = 1},
                                 .In    = {.X = -1, .Y = 1},
                                 .End   = {.X = 1, .Y = 0}};
    const Intersections met   = Self(loop, Fine);
    const double        first = (1 - std::sqrt(0.6)) / 2;
    REQUIRE(met.Hits.size() == 1);
    CHECK(met.Hits[0].Kind  == Contact::Cross);
    CHECK(std::abs(met.Hits[0].First.Along - first) <= 1e-12);
    CHECK(std::abs(met.Hits[0].Second.Along - (1 - first)) <= 1e-12);

    const std::array<Point, 8> eight  = {Point{.X = 0, .Y = 0},   Point{.X = 1, .Y = 1},
                                         Point{.X = 2, .Y =  0},
                                         Point{.X = 3, .Y = -1},
                                         Point{.X = 4, .Y =  0},
                                         Point{.X = 3, .Y =  1},
                                         Point{.X = 2, .Y = 0.5}, Point{.X = 1, .Y = -1}};
    const auto                 spline = Fit(eight, {.Closed = true});
    REQUIRE(spline.has_value());
    const Intersections twisted = Self(*spline, Fine);
    CHECK(twisted.Hits.size() == 1);
    CHECK_FALSE(twisted.Spent);
    CHECK(Self(Tests::Uneven(), Fine).Hits.empty());
  }

  TEST_CASE("self: curves that never cross themselves")
  {
    const Polyline square =
        Through({{.X = 0, .Y = 0}, {.X = 1, .Y = 0}, {.X = 1, .Y = 1}, {.X = 0, .Y = 1}}, true);
    CHECK(Self(square, Fine).Hits.empty());
    CHECK(Self(Circle({}, 2), Fine).Hits.empty());
    CHECK(Self(Line{.Start = {}, .End = {.X = 1, .Y = 0}}, Fine).Hits.empty());
    CHECK(Self(Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}, Fine).Hits.empty());
    CHECK(Self(Bezier::Quadratic{.Start = {}, .Control = {.X = 1, .Y = 2}, .End = {.X = 2, .Y = 0}},
               Fine)
              .Hits.empty());
    CHECK(Self(Ellipse{.Map = Scaling(3, 1), .Start = {}, .Sweep = Turn}, Fine).Hits.empty());
    CHECK(Self(Spline{}, Fine).Hits.empty());
    const Bezier::Cubic point = {};
    CHECK(Self(point, Fine).Hits.empty());

    const Intersections wound =
        Self(Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = Degrees(540)}, Fine);
    CHECK(wound.Overlaps.size() == 1);
    const Intersections folded = Self(Bezier::Quadratic{.Start   = {.X = 0, .Y = 0},
                                                        .Control = {.X = 2, .Y = 2},
                                                        .End     = {.X = 1, .Y = 1}},
                                      Fine);
    CHECK_FALSE(folded.Overlaps.empty());
    const Intersections deep = Self(Bezier::Quadratic{.Start   = {.X = 0, .Y = 0},
                                                      .Control = {.X = 2e6, .Y = 0},
                                                      .End     = {.X = 1e6, .Y = 0}},
                                    {.Length = 1e-3, .Radians = 1e-9});
    CHECK(deep.Spent);
  }
} // namespace Geoxy::Crossings
