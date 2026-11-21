/***************************************************************************************************
 *  @file      tests/geoxy/chamfer.cpp
 *  @brief     Tests of CHAMFER: two lines cut by two distances and by a distance and an angle, to
 *             the exact cuts' doubles; lines and curves; a polyline's corners; every misfit.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/chamfer.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/fillet.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/tolerance.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <expected>
#include <limits>
#include <numbers>
#include <optional>

namespace Geoxy
{
  namespace
  {
    /**
     *  A chamfer's misfit, or none.
     *
     *  @param  bevelled What Chamfer gave.
     *  @return          Its misfit.
     */
    auto Failing(const std::expected<Bevelled, Misfit> &bevelled) -> std::optional<Misfit>
    {
      if (bevelled)
      {
        return std::nullopt;
      }
      return bevelled.error();
    }

    constexpr Tolerance Close = {.Length = 1e-9, .Radians = 1e-9};
  } // namespace

  TEST_CASE("chamfer: two lines cut back by two distances, and by a distance and an angle")
  {
    const Line                 across = {.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const Line                 upward = {.Start = {.X = 12, .Y = -2}, .End = {.X = 12, .Y = 10}};
    const std::array<Point, 2> picks  = {Point{.X = 5, .Y = 0}, Point{.X = 12, .Y = 5}};
    const std::expected<Bevelled, Misfit> two =
        Chamfer(across, upward, {.Cut = {.First = 1, .Second = 2}, .Picks = picks});
    REQUIRE(two);
    CHECK(std::get<Line>(two->First) == Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 11, .Y = 0}});
    CHECK(std::get<Line>(two->Second)
          == Line{.Start = {.X = 12, .Y = 2}, .End = {.X = 12, .Y = 10}});
    CHECK(*two->Edge == Line{.Start = {.X = 11, .Y = 0}, .End = {.X = 12, .Y = 2}});

    const std::expected<Bevelled, Misfit> slant =
        Chamfer(across, upward, {.Cut = {.First = 1, .Slant = Degrees(45)}, .Picks = picks});
    REQUIRE(slant);
    CHECK(*slant->Edge == Line{.Start = {.X = 11, .Y = 0}, .End = {.X = 12, .Y = 1}});
    const std::expected<Bevelled, Misfit> thirty =
        Chamfer(across, upward, {.Cut = {.First = 3, .Slant = Degrees(30)}, .Picks = picks});
    REQUIRE(thirty);
    CHECK(thirty->Edge->End.Y == doctest::Approx(std::sqrt(3.0)));

    const std::expected<Bevelled, Misfit> sharp =
        Chamfer(across, Line{.Start = {.X = 4, .Y = -2}, .End = {.X = 12, .Y = 10}},
                {.Cut   = {.First = std::numbers::sqrt2, .Second = std::sqrt(13.0)},
                 .Picks = {Point{.X = 9, .Y = 0}, Point{.X = 10, .Y = 7}}});
    REQUIRE(sharp);
    CHECK(sharp->Edge->Start == Point{.X = 0x1.afd7ceef52449p+2, .Y = 0});
    CHECK(sharp->Edge->End   == Point{.X = 0x1.d555555555555p+2, .Y = 3});

    const std::expected<Bevelled, Misfit> none =
        Chamfer(across, upward, {.Cut = {}, .Picks = picks});
    REQUIRE(none);
    CHECK(!none->Edge);
    CHECK(std::get<Line>(none->First).End == Point{.X = 12, .Y = 0});

    CHECK(Failing(Chamfer(across, upward, {.Cut = {.First = 20, .Second = 1}, .Picks = picks}))
          == Misfit::TooLarge);
    CHECK(Failing(
              Chamfer(across, upward, {.Cut = {.First = 1, .Slant = Degrees(91)}, .Picks = picks}))
          == Misfit::TooLarge);
    CHECK(Failing(
              Chamfer(across, upward, {.Cut = {.First = 1, .Slant = Degrees(90)}, .Picks = picks}))
          == Misfit::TooLarge);
    CHECK(Failing(Chamfer(across, Line{.Start = {.X = 0, .Y = 1}, .End = {.X = 9, .Y = 1}},
                          {.Cut = {.First = 1, .Second = 1}, .Picks = picks}))
          == Misfit::Parallel);
  }

  TEST_CASE("chamfer: curves cut where they are the distance from the corner")
  {
    const Line across = {.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const Arc  ring   = Circle({.X = 12, .Y = 0}, 3);
    const std::expected<Bevelled, Misfit> round =
        Chamfer(across, ring,
                {.Cut   = {.First = 1, .Second = 1},
                 .Picks = {Point{.X = 5, .Y = 0}, Point{.X = 12, .Y = 3}}});
    REQUIRE(round);
    CHECK(round->Edge->Start == Point{.X = 8, .Y = 0});
    CHECK(std::abs(Distance(round->Edge->End, Point{.X = 9, .Y = 0}) - 1) < 1e-15);
    CHECK(round->Edge->End.Y > 0.0);
    CHECK(std::get<Arc>(round->Second) == ring);
    const std::expected<Bevelled, Misfit> slant =
        Chamfer(ring, across,
                {.Cut   = {.First = 1, .Slant = Degrees(60)},
                 .Picks = {Point{.X = 12, .Y = 3}, Point{.X = 5, .Y = 0}}});
    REQUIRE(slant);
    CHECK(std::abs(slant->Edge->End.Y) < 1e-15);
    CHECK(slant->Edge->End.X < 9.0);
    const Arc bow = {
        .Center = {.X = 12, .Y = 0}, .Radius = 3, .Start = Degrees(170), .Sweep = Degrees(5)};
    CHECK(Failing(Chamfer(
              across, bow,
              {.Cut = {.First = 1, .Second = 1}, .Picks = {Point{.X = 5, .Y = 0}, At(bow, 0.5)}}))
          == Misfit::TooLarge);
    CHECK(Failing(Chamfer(across, ring,
                          {.Cut   = {.First = 1, .Second = 7},
                           .Picks = {Point{.X = 5, .Y = 0}, Point{.X = 12, .Y = 3}}}))
          == Misfit::TooLarge);
    CHECK(Failing(Chamfer(across, ring,
                          {.Cut   = {.First = 1, .Slant = Degrees(170)},
                           .Picks = {Point{.X = 5, .Y = 0}, Point{.X = 12, .Y = 3}}}))
          == Misfit::TooLarge);
    const Ellipse oval = {.Map   = {.X = {.X = 4, .Y = 0}, .Y = {.X = 0, .Y = 2}},
                          .Start = Degrees(0),
                          .Sweep = Degrees(360)};
    const std::expected<Bevelled, Misfit> curved =
        Chamfer(Line{.Start = {.X = -6, .Y = 1}, .End = {.X = 6, .Y = 1}}, oval,
                {.Cut   = {.First = 0.5, .Second = 0.5},
                 .Picks = {Point{.X = 5, .Y = 1}, Point{.X = 0, .Y = 2}}});
    REQUIRE(curved);
    CHECK(std::abs(Distance(curved->Edge->End, curved->Edge->Start) - 0.5) < 0.6);
    CHECK(Failing(Chamfer(Line{.Start = {.X = -6, .Y = 5}, .End = {.X = 6, .Y = 5}}, oval,
                          {.Cut   = {.First = 1, .Second = 1},
                           .Picks = {Point{.X = 5, .Y = 5}, Point{.X = 0, .Y = 2}}}))
          == Misfit::Apart);
  }

  TEST_CASE("chamfer: misfits for numbers not finite, distances below nothing and wild angles")
  {
    const Line                 across = {.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const Line                 upward = {.Start = {.X = 12, .Y = -2}, .End = {.X = 12, .Y = 10}};
    const double               nan    = std::numeric_limits<double>::quiet_NaN();
    const std::array<Point, 2> picks  = {Point{.X = 5, .Y = 0}, Point{.X = 12, .Y = 5}};
    CHECK(Failing(Chamfer(across, upward, {.Cut = {.First = nan}, .Picks = picks}))
          == Misfit::NotFinite);
    CHECK(Failing(Chamfer(across, upward, {.Cut = {.Second = nan}, .Picks = picks}))
          == Misfit::NotFinite);
    CHECK(Failing(Chamfer(across, upward, {.Cut = {.Slant = Radians(nan)}, .Picks = picks}))
          == Misfit::NotFinite);
    CHECK(Failing(Chamfer(across, upward, {.Cut = {}, .Picks = {Point{.X = nan}, picks[1]}}))
          == Misfit::NotFinite);
    CHECK(Failing(Chamfer(across, upward, {.Cut = {}, .Picks = {picks[0], Point{.Y = nan}}}))
          == Misfit::NotFinite);
    CHECK(Failing(Chamfer(Line{.End = {.X = nan}}, upward, {.Cut = {}, .Picks = picks}))
          == Misfit::NotFinite);
    CHECK(Failing(Chamfer(across, Line{.End = {.X = nan}}, {.Cut = {}, .Picks = picks}))
          == Misfit::NotFinite);
    CHECK(Failing(Chamfer(across, upward, {.Cut = {.First = -1}, .Picks = picks}))
          == Misfit::Degenerate);
    CHECK(Failing(Chamfer(across, upward, {.Cut = {.Second = -1}, .Picks = picks}))
          == Misfit::Degenerate);
    CHECK(Failing(Chamfer(across, upward, {.Cut = {.Slant = Degrees(0)}, .Picks = picks}))
          == Misfit::Degenerate);
    CHECK(Failing(Chamfer(across, upward, {.Cut = {.Slant = Degrees(180)}, .Picks = picks}))
          == Misfit::Degenerate);
    CHECK(Failing(Chamfer(Line{}, upward, {.Cut = {}, .Picks = picks})) == Misfit::Degenerate);
    CHECK(Failing(Chamfer(across, Line{}, {.Cut = {}, .Picks = picks})) == Misfit::Degenerate);
  }

  TEST_CASE("chamfer: a polyline's corners cut where the distances fit, the rest reported")
  {
    const Polyline                        box = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                              {.Position = {.X = 10, .Y = 0}},
                                                              {.Position = {.X = 10, .Y = 4}},
                                                              {.Position = {.X = 0, .Y = 4}}},
                                                 .Closed   = true};
    const std::expected<Reshaped, Misfit> two = Chamfer(box, {.First = 1, .Second = 2}, Close);
    REQUIRE(two);
    CHECK(two->Left.empty());
    REQUIRE(two->Shape.Vertices.size()    == 8);
    CHECK(two->Shape.Vertices[0].Position == Point{.X = 2, .Y = 0});
    CHECK(two->Shape.Vertices[1].Position == Point{.X = 9, .Y = 0});
    CHECK(two->Shape.Vertices[2].Position == Point{.X = 10, .Y = 2});

    const std::expected<Reshaped, Misfit> slant =
        Chamfer(box, {.First = 1, .Slant = Degrees(45)}, Close);
    REQUIRE(slant);
    CHECK(slant->Shape.Vertices[2].Position == Point{.X = 10, .Y = 1});

    const std::expected<Reshaped, Misfit> halves = Chamfer(box, {.First = 2, .Second = 2}, Close);
    REQUIRE(halves);
    CHECK(halves->Shape.Vertices.size() == 6);

    const std::expected<Reshaped, Misfit> wide = Chamfer(box, {.First = 3, .Second = 3}, Close);
    REQUIRE(wide);
    CHECK(wide->Left.size() == 2);

    const Polyline bent = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                        {.Position = {.X = 10, .Y = 0}, .Bulge = 0.5},
                                        {.Position = {.X = 10, .Y = 10}}}};
    const std::expected<Reshaped, Misfit> bowed = Chamfer(bent, {.First = 1, .Second = 1}, Close);
    REQUIRE(bowed);
    CHECK(bowed->Left.empty());
    CHECK(bowed->Shape.Vertices.size() == 4);

    CHECK(Chamfer(box, {}, Close).error() == Misfit::Degenerate);
    CHECK(Chamfer(box, {.Slant = Degrees(30)}, Close).error()     == Misfit::Degenerate);
    CHECK(Chamfer(box, {.First = -1, .Second = 1}, Close).error() == Misfit::Degenerate);
    CHECK(Chamfer(
              Polyline{.Vertices = {{.Position = {.X = std::numeric_limits<double>::infinity()}}}},
              {.First = 1}, Close)
              .error()
          == Misfit::NotFinite);
    const Polyline                        back = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                               {.Position = {.X = 10, .Y = 0}},
                                                               {.Position = {.X = 5, .Y = 0}}}};
    const std::expected<Reshaped, Misfit> reversed =
        Chamfer(back, {.First = 1, .Second = 1}, Close);
    REQUIRE(reversed);
    CHECK(reversed->Left[0].Why == Misfit::Parallel);
  }

  TEST_CASE("chamfer: rays crossing twice, cuts of nothing, picks on the corner, arcs first")
  {
    const Line                 across = {.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const Arc                  ring   = Circle({.X = 12, .Y = 0}, 3);
    const std::array<Point, 2> picks  = {Point{.X = 5, .Y = 0}, Point{.X = 12, .Y = 3}};
    const std::expected<Bevelled, Misfit> twice =
        Chamfer(across, ring, {.Cut = {.First = 1, .Slant = Degrees(30)}, .Picks = picks});
    REQUIRE(twice);
    CHECK(twice->Edge->End.X < 10.0);
    const std::expected<Bevelled, Misfit> nought =
        Chamfer(across, ring, {.Cut = {.First = 0, .Second = 1}, .Picks = picks});
    REQUIRE(nought);
    CHECK(nought->Edge->Start == Point{.X = 9, .Y = 0});
    const std::expected<Bevelled, Misfit> slant = Chamfer(
        ring, across, {.Cut = {.First = 0, .Slant = Degrees(30)}, .Picks = {picks[1], picks[0]}});
    REQUIRE(slant);
    CHECK(!slant->Edge);
    CHECK(Failing(Chamfer(ring, across,
                          {.Cut = {.First = 7, .Second = 1}, .Picks = {picks[1], picks[0]}}))
          == Misfit::TooLarge);
    const std::expected<Bevelled, Misfit> picked =
        Chamfer(across, Line{.Start = {.X = 12, .Y = -2}, .End = {.X = 12, .Y = 10}},
                {.Cut   = {.First = 1, .Second = 1},
                 .Picks = {Point{.X = 5, .Y = 0}, Point{.X = 12, .Y = 0}}});
    REQUIRE(picked);
    CHECK(picked->Edge->End == Point{.X = 12, .Y = 1});

    const Polyline                        box    = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                                 {.Position = {.X = 10, .Y = 0}},
                                                                 {.Position = {.X = 10, .Y = 4}},
                                                                 {.Position = {.X = 0, .Y = 4}}},
                                                    .Closed   = true};
    const std::expected<Reshaped, Misfit> second = Chamfer(box, {.Second = 1}, Close);
    REQUIRE(second);
    CHECK(second->Shape.Vertices[0].Position == Point{.X = 1, .Y = 0});
    const Polyline arched = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0.5},
                                          {.Position = {.X = 10, .Y = 0}},
                                          {.Position = {.X = 10, .Y = 10}}}};
    const std::expected<Reshaped, Misfit> arcFirst =
        Chamfer(arched, {.First = 1, .Second = 1}, Close);
    REQUIRE(arcFirst);
    CHECK(arcFirst->Left.empty());
  }
} // namespace Geoxy
