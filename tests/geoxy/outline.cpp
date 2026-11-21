/***************************************************************************************************
 *  @file      tests/geoxy/outline.cpp
 *  @brief     Tests of wide polylines' outlines: corners worked by hand, and bands' areas against
 *             the trapezoid's, the annulus's and the spiral band's closed forms.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/outline.hpp"

#include "real.hpp"

#include "geoxy/area.hpp"
#include "geoxy/flatten.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A polyline through points, every bulge 0 but those given.
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

    /**
     *  A band's corners.
     *
     *  @param  band The band.
     *  @return      Its vertices' positions.
     */
    auto Corners(const Polyline &band) -> std::vector<Point>
    {
      std::vector<Point> corners;
      corners.reserve(band.Vertices.size());
      for (const Vertex &vertex : band.Vertices)
      {
        corners.push_back(vertex.Position);
      }
      return corners;
    }

    constexpr std::array<Width, 1> WidthTwo = {Width{.Start = 2, .End = 2}};
  } // namespace

  TEST_CASE("outline: a straight segment's trapezoid")
  {
    const Polyline              line  = Through({{.X = 0, .Y = 0}, {.X = 10, .Y = 0}}, false);
    const std::array<Width, 1>  wide  = {Width{.Start = 2, .End = -4}};
    const std::vector<Polyline> bands = Outline(line, wide, {});
    REQUIRE(bands.size() == 1);
    CHECK(bands[0].Closed);
    CHECK(Corners(bands[0])
          == std::vector<Point>{
              {.X = 0, .Y = -1}, {.X = 10, .Y = -2}, {.X = 10, .Y = 2}, {.X = 0, .Y = 1}});
    CHECK(Area(bands[0]) == doctest::Approx(30).epsilon(1e-15));
  }

  TEST_CASE("outline: straight corners mitered, both bands ending at one point")
  {
    const Polyline ell = Through({{.X = 0, .Y = 0}, {.X = 10, .Y = 0}, {.X = 10, .Y = 10}}, false);
    const std::vector<Polyline> bands = Outline(ell, WidthTwo, {});
    REQUIRE(bands.size() == 2);
    CHECK(Corners(bands[0])
          == std::vector<Point>{
              {.X = 0, .Y = -1}, {.X = 11, .Y = -1}, {.X = 9, .Y = 1}, {.X = 0, .Y = 1}});
    CHECK(Corners(bands[1])
          == std::vector<Point>{
              {.X = 11, .Y = -1}, {.X = 11, .Y = 10}, {.X = 9, .Y = 10}, {.X = 9, .Y = 1}});

    const Polyline             square = Through({{.X = 5e6, .Y = 5e6},
                                                 {.X = 5e6 + 3, .Y = 5e6 + 1},
                                                 {.X = 5e6 + 2, .Y = 5e6 + 4},
                                                 {.X = 5e6 - 1, .Y = 5e6 + 3}},
                                                true);
    const std::array<Width, 4> widths = {
        Width{.Start =  0.5, .End = 0.25},
        Width{.Start = 0.25, .End =    1},
        Width{.Start = 1, .End = 1}, Width{.Start = 1, .End = 0.5}};
    const std::vector<Polyline> loop = Outline(square, widths, {});
    REQUIRE(loop.size() == 4);
    for (std::size_t index = 0; index < 4; ++index)
    {
      const Polyline &band = loop[index];
      const Polyline &next = loop[(index + 1) % 4];
      CHECK(band.Vertices[1].Position == next.Vertices[0].Position);
      CHECK(band.Vertices[2].Position == next.Vertices[3].Position);
      CHECK(Area(band) > 0.0);
    }

    const std::array<Width, 4> gapped = {Width{}, Width{.Start = 1, .End = 1},
                                         Width{.Start = 1, .End = 1}, Width{.Start = 1, .End = 1}};
    CHECK(Outline(square, gapped, {}).size() == 3);
    const std::array<Width, 4> last = {Width{.Start = 1, .End = 1}, Width{.Start = 1, .End = 1},
                                       Width{.Start = 1, .End = 1}, Width{}};
    CHECK(Outline(square, last, {}).size() == 3);

    const Polyline sharp = Through({{.X = 0, .Y = 0}, {.X = 10, .Y = 0}, {.X = 0, .Y = 1}}, false);
    const std::vector<Polyline> spike = Outline(sharp, WidthTwo, {});
    REQUIRE(spike.size() == 2);
    CHECK(spike[0].Vertices[1].Position.X > 20);

    const Polyline straight =
        Through({{.X = 0, .Y = 0}, {.X = 5, .Y = 0}, {.X = 10, .Y = 0}}, false);
    const std::vector<Polyline> butt = Outline(straight, WidthTwo, {});
    REQUIRE(butt.size() == 2);
    CHECK(butt[0].Vertices[1].Position == Point{.X = 5, .Y = -1});
    CHECK(butt[1].Vertices[0].Position == Point{.X = 5, .Y = -1});
  }

  TEST_CASE("outline: an arc's band, of one width or of a changing one")
  {
    for (const double bulge : {1.0, -1.0})
    {
      const Polyline half = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = bulge},
                                          {.Position = {.X = 2, .Y = 0}}},
                             .Closed   = false};
      const std::array<Width, 1>  wide  = {Width{.Start = 0.5, .End = 0.5}};
      const std::vector<Polyline> bands = Outline(half, wide, {});
      REQUIRE(bands.size()             == 1);
      REQUIRE(bands[0].Vertices.size() == 4);
      CHECK(Exactly(bands[0].Vertices[0].Bulge, bulge));
      CHECK(Exactly(bands[0].Vertices[2].Bulge, -bulge));
      CHECK(Area(bands[0]) == doctest::Approx(std::numbers::pi / 2).epsilon(1e-15));

      const std::array<Width, 1>  taper    = {Width{.Start = 0.2, .End = 0.6}};
      const double                accuracy = 1e-6;
      const std::vector<Polyline> spirals  = Outline(half, taper, {.Error = accuracy});
      REQUIRE(spirals.size() == 1);
      const std::size_t count = spirals[0].Vertices.size();
      CHECK(count > 100);
      CHECK(count < 4000);
      const double exact = 1.0 * std::numbers::pi * (0.1 + 0.3);
      CHECK(std::abs(Area(spirals[0]) - exact) <= accuracy * 8);
    }
    const Polyline quarter = {
        .Vertices = {{.Position = {.X = 1, .Y = 0}, .Bulge = 0.4142135623730951},
                     {.Position = {.X = 0, .Y = 1}}},
        .Closed   = false};
    const std::array<Width, 1> taper = {Width{.Start = 0, .End = 0.5}};
    CHECK(Outline(quarter, taper, {}).at(0).Vertices.size() == 2 * (MostChords + 1));
  }

  TEST_CASE("outline: arcs and lines meet butt")
  {
    const Polyline              turn  = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                                      {.Position = {.X = 4, .Y = 0}, .Bulge = 1},
                                                      {.Position = {.X = 4, .Y = 4}, .Bulge = 0},
                                                      {.Position = {.X = 0, .Y = 4}, .Bulge = 0}},
                                         .Closed   = false};
    const std::vector<Polyline> bands = Outline(turn, WidthTwo, {});
    REQUIRE(bands.size() == 3);
    CHECK(Distance(bands[0].Vertices[1].Position, bands[1].Vertices[0].Position) <= 1e-15);
    CHECK(Distance(bands[1].Vertices[1].Position, bands[2].Vertices[0].Position) <= 1e-14);
    CHECK(bands[0].Vertices[1].Position == Point{.X = 4, .Y = -1});
  }

  TEST_CASE("outline: what has no band")
  {
    const Polyline path = Through(
        {{.X = 0, .Y = 0}, {.X = 4, .Y = 0}, {.X = 4, .Y = 0}, {.X = 4, .Y = 3}, {.X = 0, .Y = 3}},
        false);
    CHECK(Outline(path, {}, {}).empty());
    const std::array<Width, 1> none = {Width{}};
    CHECK(Outline(path, none, {}).empty());
    CHECK(Outline(path, WidthTwo, {}).size() == 3);
    const std::array<Width, 3>  some  = {Width{.Start = 1, .End = 1}, Width{.Start = 1, .End = 1},
                                         Width{}};
    const std::vector<Polyline> first = Outline(path, some, {});
    REQUIRE(first.size() == 1);
    CHECK(first[0].Vertices[1].Position == Point{.X = 4, .Y = -0.5});
    const std::array<Width, 2> growing = {Width{.Start = 0, .End = 1}, Width{}};
    CHECK(Outline(path, growing, {}).size() == 1);
  }
} // namespace Geoxy
