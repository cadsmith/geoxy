/***************************************************************************************************
 *  @file      tests/geoxy/loop.cpp
 *  @brief     Tests of loops: winding numbers worked by hand for lines, arcs, polylines, ellipses,
 *             Bézier curves and splines, a ray through every vertex, points a tolerance either side
 *             of an edge, at the origin and far from it; areas and perimeters in closed form.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/loop.hpp"

#include "curves.hpp"
#include "loops.hpp"
#include "winding.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/area.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/region.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A tolerance fine enough for every case drawn about unit size.
     */
    constexpr Tolerance LoopFine = {.Length = 1e-9, .Radians = 1e-9};

    /**
     *  The places each case is drawn at: the origin, survey coordinates and far past them.
     */
    constexpr std::array<double, 3> LoopOffsets = {0.0, 5e6, 1e8};

    /**
     *  A tolerance for a case drawn about a place.
     *
     *  @param  offset The place's coordinates.
     *  @return        The tolerance.
     */
    auto Near(double offset) -> Tolerance
    {
      return Tolerance::Of(
          {.Min = {.X = offset, .Y = offset}, .Max = {.X = offset + 10, .Y = offset + 10}});
    }

    /**
     *  A point moved from the origin by a place in both coordinates.
     *
     *  @param  offset The place.
     *  @param  x      Its x from there.
     *  @param  y      Its y from there.
     *  @return        The point.
     */
    auto Spot(double offset, double x, double y) -> Point
    {
      return {.X = offset + x, .Y = offset + y};
    }

    /**
     *  A loop moved by a place in both coordinates.
     *
     *  @param  loop   A loop of lines and arcs.
     *  @param  offset The place.
     *  @return        It moved.
     */
    auto Placed(const Loop &loop, double offset) -> Loop
    {
      return Tests::Moved(Region{.Loops = {loop}}, {.X = offset, .Y = offset}).Loops.front();
    }

    /**
     *  The winding number of a loop about a point, 99 where it has none.
     *
     *  @param  loop      The loop.
     *  @param  point     The point.
     *  @param  tolerance The tolerance.
     *  @return           The winding number.
     */
    auto Wound(const Loop &loop, Point point, const Tolerance &tolerance) -> int
    {
      return Winding(loop, point, tolerance).value_or(99);
    }

    /**
     *  A closed polyline: a square of side 2 from the origin, its corners rounded by quarter
     *  circles of radius 1/2 as bulges.
     *
     *  @return The polyline.
     */
    auto RoundedSquare() -> Polyline
    {
      const double bulge = std::tan(std::numbers::pi / 8);
      return {.Vertices = {{.Position = {.X = 0.5, .Y = 0}},
                           {.Position = {.X = 1.5, .Y = 0}, .Bulge = bulge},
                           {.Position = {.X = 2, .Y = 0.5}},
                           {.Position = {.X = 2, .Y = 1.5}, .Bulge = bulge},
                           {.Position = {.X = 1.5, .Y = 2}},
                           {.Position = {.X = 0.5, .Y = 2}, .Bulge = bulge},
                           {.Position = {.X = 0, .Y = 1.5}},
                           {.Position = {.X = 0, .Y = 0.5}, .Bulge = bulge}},
              .Closed   = true};
    }
  } // namespace

  TEST_CASE("loop: a square winds once about a point inside, at every scale")
  {
    for (const double offset : LoopOffsets)
    {
      const Tolerance tolerance = Near(offset);
      const Loop      square    = Placed(Tests::Square({}, 2), offset);
      const Loop      back      = Tests::Reversed(square);
      CHECK(Wound(square, Spot(offset, 1, 1), tolerance)  == 1);
      CHECK(Wound(back, Spot(offset, 1, 1), tolerance)    == -1);
      CHECK(Wound(square, Spot(offset, 3, 1), tolerance)  == 0);
      CHECK(Wound(square, Spot(offset, -1, 1), tolerance) == 0);
      CHECK(Wound(square, Spot(offset, 1, 3), tolerance)  == 0);
      CHECK(Wound(square, Spot(offset, 1, -1), tolerance) == 0);
      CHECK(Wound(square, Spot(offset, 1, 0), tolerance)  == 99);
      CHECK(Wound(square, Spot(offset, 2, 2), tolerance)  == 99);
      CHECK(square       == Placed(Tests::Square({}, 2), offset));
      CHECK_FALSE(square == back);
    }
  }

  TEST_CASE("loop: a point a tolerance either side of an edge is on it, or not")
  {
    const Loop square = Tests::Square({}, 2);
    CHECK(Wound(square, {.X = 1, .Y = 0.5e-9}, LoopFine)  == 99);
    CHECK(Wound(square, {.X = 1, .Y = -0.5e-9}, LoopFine) == 99);
    CHECK(Wound(square, {.X = 1, .Y = 2e-9}, LoopFine)    == 1);
    CHECK(Wound(square, {.X = 1, .Y = -2e-9}, LoopFine)   == 0);
    CHECK(Wound(square, {.X = 1, .Y = 1e-300}, {})        == 1);
    CHECK(Wound(square, {.X = 1, .Y = -1e-300}, {.Length = -1})    == 0);
    CHECK(Wound(square, {.X = 1, .Y = 1e-300}, {.Length = 1e-299}) == 99);
    CHECK(Wound(square, {.X = 1, .Y = 0}, {}) == 99);
  }

  TEST_CASE("loop: a ray through a vertex counts it once")
  {
    const Loop diamond =
        Tests::Polygon({{.X = 0, .Y = -1}, {.X = 1, .Y = 0}, {.X = 0, .Y = 1}, {.X = -1, .Y = 0}});
    CHECK(Wound(diamond, {.X = 0, .Y = 0}, LoopFine)    == 1);
    CHECK(Wound(diamond, {.X = -2, .Y = 0}, LoopFine)   == 0);
    CHECK(Wound(diamond, {.X = 2, .Y = 0}, LoopFine)    == 0);
    CHECK(Wound(diamond, {.X = -0.5, .Y = 0}, LoopFine) == 1);
    const Loop teeth = Tests::Polygon({{.X = 0, .Y = 0},
                                       {.X = 4, .Y = 0},
                                       {.X = 4, .Y = 2},
                                       {.X = 3, .Y = 1},
                                       {.X = 2, .Y = 2},
                                       {.X = 1, .Y = 1},
                                       {.X = 0, .Y = 2}});
    CHECK(Wound(teeth, {.X = 0.5, .Y = 1}, LoopFine)    == 1);
    CHECK(Wound(teeth, {.X = 1.5, .Y = 1}, LoopFine)    == 1);
    CHECK(Wound(teeth, {.X = -1, .Y = 1}, LoopFine)     == 0);
    CHECK(Wound(teeth, {.X = 1.5, .Y = 1.75}, LoopFine) == 0);
    CHECK(Wound(teeth, {.X = -1, .Y = 2}, LoopFine)     == 0);
    CHECK(Wound(teeth, {.X = 2, .Y = 1.9}, LoopFine)    == 1);
  }

  TEST_CASE("loop: a bow-tie winds once each way, and a loop twice about its middle")
  {
    const Loop bow =
        Tests::Polygon({{.X = 0, .Y = 0}, {.X = 2, .Y = 2}, {.X = 2, .Y = 0}, {.X = 0, .Y = 2}});
    CHECK(Wound(bow, {.X = 1.5, .Y = 1}, LoopFine) == -1);
    CHECK(Wound(bow, {.X = 0.5, .Y = 1}, LoopFine) == 1);
    CHECK(Wound(bow, {.X = 1, .Y = 1}, LoopFine)   == 99);
    CHECK(Area(bow, Accuracy{}) == 0.0);
    const Loop twice = {
        .Edges = {Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = Degrees(720)}}};
    CHECK(Wound(twice, {.X = 0.5, .Y = 0}, LoopFine)  == 2);
    CHECK(Wound(twice, {.X = 1.5, .Y = 0}, LoopFine)  == 0);
    CHECK(Wound(Tests::Reversed(twice), {}, LoopFine) == -2);
    CHECK(Area(twice, Accuracy{}) == doctest::Approx(2 * std::numbers::pi).epsilon(1e-15));
  }

  TEST_CASE("loop: arcs wind by their quarters, inside and outside their circle, at every scale")
  {
    for (const double offset : LoopOffsets)
    {
      const Tolerance tolerance = Near(offset);
      const Point     center    = Spot(offset, 0, 0);
      const Loop      half      = {
          .Edges = {
              Arc{.Center = center, .Radius = 1, .Start = Degrees(-90), .Sweep = Degrees(180)},
              Line{.Start = Spot(offset, 0, 1), .End = Spot(offset, 0, -1)}}};
      CHECK(Wound(half, Spot(offset, 0.5, 0), tolerance)    == 1);
      CHECK(Wound(half, Spot(offset, 0.5, 0.5), tolerance)  == 1);
      CHECK(Wound(half, Spot(offset, 0.9, 0.5), tolerance)  == 0);
      CHECK(Wound(half, Spot(offset, -0.5, 0), tolerance)   == 0);
      CHECK(Wound(half, Spot(offset, -0.5, 0.9), tolerance) == 0);
      CHECK(Wound(half, Spot(offset, 0.1, 0.1), tolerance)  == 1);
      CHECK(Wound(half, Spot(offset, 0.71, 0.7), tolerance) == 1);
      CHECK(Wound(half, Spot(offset, 0.72, 0.7), tolerance) == 0);
      CHECK(Wound(half, Spot(offset, -1, 0), tolerance)     == 0);
      CHECK(Wound(half, Spot(offset, 2, 0), tolerance)      == 0);
      CHECK(Wound(half, Spot(offset, 0.5, 1), tolerance)    == 0);
      CHECK(Wound(half, Spot(offset, 0.5, -1), tolerance)   == 0);
      CHECK(Wound(Tests::Reversed(half), Spot(offset, 0.5, 0), tolerance) == -1);
      const Loop west = {
          .Edges = {Arc{.Center = center, .Radius = 1, .Start = Degrees(90), .Sweep = Degrees(180)},
                    Line{.Start = Spot(offset, 0, -1), .End = Spot(offset, 0, 1)}}};
      CHECK(Wound(west, Spot(offset, -0.5, 0), tolerance)    == 1);
      CHECK(Wound(west, Spot(offset, -0.5, 0.5), tolerance)  == 1);
      CHECK(Wound(west, Spot(offset, -0.72, 0.7), tolerance) == 0);
      CHECK(Wound(west, Spot(offset, -0.71, 0.7), tolerance) == 1);
      CHECK(Wound(west, Spot(offset, -2, 0.5), tolerance)    == 0);
      CHECK(Wound(west, Spot(offset, 0.5, 0.5), tolerance)   == 0);
      CHECK(Wound(Tests::Reversed(west), Spot(offset, -0.5, 0.1), tolerance) == -1);
    }
  }

  TEST_CASE("loop: a circle winds about its inside, a point on it beyond an arc's sweep not on it")
  {
    for (const double offset : LoopOffsets)
    {
      const Tolerance tolerance = Near(offset);
      const Loop      disc      = Placed(Tests::Disc({}, 1), offset);
      CHECK(Wound(disc, Spot(offset, 0, 0), tolerance)     == 1);
      CHECK(Wound(disc, Spot(offset, 0.99, 0), tolerance)  == 1);
      CHECK(Wound(disc, Spot(offset, 1.01, 0), tolerance)  == 0);
      CHECK(Wound(disc, Spot(offset, 0, 0.99), tolerance)  == 1);
      CHECK(Wound(disc, Spot(offset, 0, -0.99), tolerance) == 1);
      CHECK(Wound(disc, Spot(offset, 1, 0), tolerance)     == 99);
      CHECK(Wound(Tests::Reversed(disc), Spot(offset, 0, 0.5), tolerance) == -1);
      const Loop cap = {.Edges = {Arc{.Center = Spot(offset, 0, 0),
                                      .Radius = 1,
                                      .Start  = Degrees(30),
                                      .Sweep  = Degrees(120)},
                                  Line{.Start = Spot(offset, -std::sqrt(0.75), 0.5),
                                       .End   = Spot(offset, std::sqrt(0.75), 0.5)}}};
      CHECK(Wound(cap, Spot(offset, 0, 0.75), tolerance) == 1);
      CHECK(Wound(cap, Spot(offset, 0, -1), tolerance)   == 0);
      CHECK(Wound(cap, Spot(offset, std::sqrt(0.75), -0.5), tolerance)  == 0);
      CHECK(Wound(cap, Spot(offset, -std::sqrt(0.75), -0.5), tolerance) == 0);
      CHECK(Wound(cap, Spot(offset, -0.6, 0.6), tolerance)  == 1);
      CHECK(Wound(cap, Spot(offset, -0.85, 0.6), tolerance) == 0);
      CHECK(Wound(cap, Spot(offset, -0.8, 0.6), tolerance)  == 99);
    }
  }

  TEST_CASE("loop: an arc of a negative radius, held in radians, winds as the arc it draws")
  {
    const Loop flipped = {.Edges = {Arc{.Center = {},
                                        .Radius = -1,
                                        .Start  = Radians(std::numbers::pi / 2),
                                        .Sweep  = Radians(std::numbers::pi)},
                                    Line{.Start = {.X = 0, .Y = 1}, .End = {.X = 0, .Y = -1}}}};
    CHECK(Wound(flipped, {.X = 0.5, .Y = 0}, LoopFine)   == 1);
    CHECK(Wound(flipped, {.X = -0.5, .Y = 0}, LoopFine)  == 0);
    CHECK(Wound(flipped, {.X = 0.5, .Y = 0.5}, LoopFine) == 1);
    const Loop point = {
        .Edges = {Arc{.Center = {.X = 1, .Y = 1}, .Radius = 1, .Start = {}, .Sweep = {}},
                  Line{.Start = {.X = 2, .Y = 1}, .End = {.X = 0, .Y = 0}},
                  Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 2, .Y = 1}}}};
    CHECK(Wound(point, {.X = 5, .Y = 5}, LoopFine) == 0);
  }

  TEST_CASE("loop: an arc of many turns counts each, and a count past an int has none")
  {
    const Loop many = {
        .Edges = {Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = Degrees(360.0 * 0x1p30)}}};
    CHECK(Wound(many, {}, LoopFine) == 1 << 24);
    Loop coils = {};
    for (int index = 0; index < 200; ++index)
    {
      coils.Edges.push_back(many.Edges.front());
    }
    CHECK(Wound(coils, {}, LoopFine) == 99);
    CHECK(Wound(Tests::Reversed(coils), {}, LoopFine) == 99);
    CHECK(Wound(coils, {.X = 2, .Y = 0}, LoopFine)    == 0);
  }

  TEST_CASE("loop: a gap is closed by its line, which a point may be on")
  {
    const Loop open = {.Edges = {Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 2, .Y = 0}},
                                 Line{.Start = {.X = 2, .Y = 0}, .End = {.X = 2, .Y = 2}}}};
    CHECK(Wound(open, {.X = 1.5, .Y = 0.5}, LoopFine) == 1);
    CHECK(Wound(open, {.X = 0.5, .Y = 1.5}, LoopFine) == 0);
    CHECK(Wound(open, {.X = 1, .Y = 1}, LoopFine)     == 99);
    CHECK(Area(open, Accuracy{})   == 2.0);
    CHECK(Length(open, Accuracy{}) == doctest::Approx(4 + std::sqrt(8.0)).epsilon(1e-15));
    CHECK_FALSE(Closed(open, LoopFine));
    CHECK(Closed(open, {.Length = 3}));
  }

  TEST_CASE("loop: polylines wind segment by segment, alone or among other edges")
  {
    const Polyline rounded = RoundedSquare();
    const Loop     alone   = {.Edges = {rounded}};
    CHECK(Wound(alone, {.X = 1, .Y = 1}, LoopFine)       == 1);
    CHECK(Wound(alone, {.X = 0.05, .Y = 0.05}, LoopFine) == 0);
    CHECK(Wound(alone, {.X = 0.2, .Y = 0.2}, LoopFine)   == 1);
    CHECK(Area(alone, Accuracy{})
          == doctest::Approx(4 - 1 + (std::numbers::pi / 4)).epsilon(1e-15));
    CHECK(Area(alone, Accuracy{})   == Area(Curve{rounded}, Accuracy{}));
    CHECK(Length(alone, Accuracy{}) == doctest::Approx(4 + std::numbers::pi).epsilon(1e-15));
    CHECK(Closed(alone, LoopFine));
    const Polyline open = {
        .Vertices = {{.Position = {.X = 2, .Y = 0}, .Bulge = 1}, {.Position = {.X = 0, .Y = 0}}}};
    const Loop mixed = {.Edges = {open, Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 2, .Y = 0}}}};
    CHECK(Wound(mixed, {.X = 1, .Y = 0.5}, LoopFine)  == 1);
    CHECK(Wound(mixed, {.X = 1, .Y = -0.5}, LoopFine) == 0);
    CHECK(Area(mixed, Accuracy{}) == doctest::Approx(std::numbers::pi / 2).epsilon(1e-15));
    CHECK(Closed(mixed, LoopFine));
  }

  TEST_CASE("loop: ellipses, Bézier curves and splines wind by halving, near and far")
  {
    const Ellipse ellipse = {.Map   = {.X = {.X = 2, .Y = 0}, .Y = {.X = 0, .Y = 1}, .Origin = {}},
                             .Start = {},
                             .Sweep = Turn};
    const Loop    oval    = {.Edges = {ellipse}};
    CHECK(Wound(oval, {}, LoopFine) == 1);
    CHECK(Wound(oval, {.X = 1.9, .Y = 0.3}, LoopFine)    == 1);
    CHECK(Wound(oval, {.X = 1.99, .Y = 0.2}, LoopFine)   == 0);
    CHECK(Wound(oval, {.X = 2 - 1e-8, .Y = 0}, LoopFine) == 1);
    CHECK(Wound(oval, {.X = 2 + 1e-8, .Y = 0}, LoopFine) == 0);
    CHECK(Wound(oval, {.X = 5, .Y = 0}, LoopFine)        == 0);
    CHECK(Wound(oval, {.X = 2, .Y = 0}, LoopFine)        == 99);
    CHECK(Area(oval, Accuracy{}) == doctest::Approx(2 * std::numbers::pi).epsilon(1e-15));
    const Loop wave = {
        .Edges = {Tests::Wave, Line{.Start = Tests::Wave.End, .End = Tests::Wave.Start}}};
    CHECK(Wound(wave, {.X = 1, .Y = 0.5}, LoopFine) == -1);
    CHECK(Wound(wave, {.X = 3, .Y = 0.4}, LoopFine) == 1);
    CHECK(Wound(wave, {.X = 2, .Y = 3}, LoopFine)   == 0);
    const Loop clamped = {
        .Edges = {Tests::Clamped(), Line{.Start = Tests::Wave.End, .End = Tests::Wave.Start}}};
    CHECK(Wound(clamped, {.X = 1, .Y = 0.5}, LoopFine) == -1);
    CHECK(Wound(clamped, {.X = 3, .Y = 0.4}, LoopFine) == 1);
    const Bezier::Quadratic cap = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 2, .Y = 0}};
    const Loop dome = {.Edges = {cap, Line{.Start = cap.End, .End = cap.Start}}};
    CHECK(Wound(dome, {.X = 1, .Y = 0.5}, LoopFine) == -1);
    CHECK(Wound(dome, {.X = 1, .Y = 1.5}, LoopFine) == 0);
    CHECK(Area(dome, Accuracy{}) == doctest::Approx(-4.0 / 3.0).epsilon(1e-15));
    const Loop quarter = {.Edges = {Tests::Quarter(), Line{.Start = {.X = 0, .Y = 1}, .End = {}},
                                    Line{.Start = {}, .End = {.X = 1, .Y = 0}}}};
    CHECK(Wound(quarter, {.X = 0.7, .Y = 0.7}, LoopFine)   == 1);
    CHECK(Wound(quarter, {.X = 0.71, .Y = 0.71}, LoopFine) == 0);
  }

  TEST_CASE("loop: halving stops where a point is on the curve")
  {
    const Curve ellipse =
        Ellipse{.Map   = {.X = {.X = 2, .Y = 0}, .Y = {.X = 0, .Y = 1}, .Origin = {}},
                .Start = {},
                .Sweep = Degrees(180)};
    const Traverse traverse = {
        .Shape = &ellipse, .From = {.X = 2, .Y = 0}, .To = {.X = -2, .Y = 0}};
    CHECK(std::abs(Ascents(traverse, {.X = 0, .Y = 1})) <= 1);
  }

  TEST_CASE("loop: a point or an edge that is not finite, or an edge without end, has no winding")
  {
    constexpr double Nan    = std::numeric_limits<double>::quiet_NaN();
    constexpr double Huge   = std::numeric_limits<double>::infinity();
    const Loop       square = Tests::Square({}, 2);
    CHECK(Wound(square, {.X = Nan, .Y = 1}, LoopFine)  == 99);
    CHECK(Wound(square, {.X = 1, .Y = Huge}, LoopFine) == 99);
    Loop broken          = square;
    broken.Edges.front() = Line{.Start = {.X = Nan, .Y = 0}, .End = {.X = 2, .Y = 0}};
    CHECK(Wound(broken, {.X = 1, .Y = 1}, LoopFine) == 99);
    CHECK_FALSE(Closed(broken, LoopFine));
    const Loop rayed = {.Edges = {Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}}};
    CHECK(Wound(rayed, {.X = 1, .Y = 1}, LoopFine) == 99);
    CHECK(Area(rayed, Accuracy{}) == 0.0);
    CHECK(std::isinf(Length(rayed, Accuracy{})));
    CHECK_FALSE(Closed(rayed, LoopFine));
    const Loop endless = {.Edges = {Line{.Start = {}, .End = {.X = 1, .Y = 0}},
                                    Xline{.Through = {}, .Way = {.X = 0, .Y = 1}}}};
    CHECK(Area(endless, Accuracy{}) == 0.0);
    CHECK(Wound(endless, {.X = 1, .Y = 1}, LoopFine) == 99);
    const Loop empty = {};
    CHECK(Wound(empty, {}, LoopFine) == 0);
    CHECK(Area(empty, Accuracy{})    == 0.0);
    CHECK(Length(empty, Accuracy{})  == 0.0);
    CHECK_FALSE(Closed(empty, LoopFine));
    const Loop none = {.Edges = {Polyline{}}};
    CHECK(Wound(none, {}, LoopFine) == 0);
    const Loop nan = {.Edges = {Arc{.Center = {}, .Radius = Nan, .Start = {}, .Sweep = Turn}}};
    CHECK(Wound(nan, {}, LoopFine) == 99);
  }

  TEST_CASE("loop: areas and perimeters of squares, circles and slots, at every scale")
  {
    for (const double offset : LoopOffsets)
    {
      CHECK(Area(Placed(Tests::Square({}, 2), offset), Accuracy{}) == 4.0);
      CHECK(Area(Tests::Reversed(Placed(Tests::Square({}, 2), offset)), Accuracy{}) == -4.0);
      CHECK(Length(Placed(Tests::Square({}, 2), offset), Accuracy{}) == 8.0);
      CHECK(Area(Placed(Tests::Disc({}, 1), offset), Accuracy{})
            == doctest::Approx(std::numbers::pi).epsilon(1e-15));
      CHECK(Length(Placed(Tests::Disc({}, 1), offset), Accuracy{})
            == doctest::Approx(2 * std::numbers::pi).epsilon(1e-15));
      const Loop slot = Placed(Tests::Slot({}, 3, 1), offset);
      CHECK(Area(slot, Accuracy{})   == doctest::Approx(6 + std::numbers::pi).epsilon(1e-15));
      CHECK(Length(slot, Accuracy{}) == doctest::Approx(6 + (2 * std::numbers::pi)).epsilon(1e-15));
      CHECK(Closed(slot, Near(offset)));
      CHECK(Wound(slot, Spot(offset, -0.5, 0.5), Near(offset)) == 1);
      CHECK(Wound(slot, Spot(offset, 3.5, -0.5), Near(offset)) == 1);
      CHECK(Wound(slot, Spot(offset, 3.9, -0.9), Near(offset)) == 0);
    }
  }
} // namespace Geoxy
