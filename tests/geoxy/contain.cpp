/***************************************************************************************************
 *  @file      tests/geoxy/contain.cpp
 *  @brief     Tests of what a window and a crossing window pick: lines against an exact reference
 *             of their ends and the box's edges, and curves of other kinds inside, across, near and
 *             far from a box.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/contain.hpp"

#include "draws.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/line.hpp"
#include "geoxy/loop.hpp"
#include "geoxy/orient.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/region.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether a point collinear with two others lies between them, ends included.
     *
     *  @param  start One end.
     *  @param  end   The other.
     *  @param  point The point.
     *  @return       Whether it does.
     */
    auto Between(Point start, Point end, Point point) -> bool
    {
      return std::min(start.X, end.X) <= point.X
          && point.X                  <= std::max(start.X, end.X)
          && std::min(start.Y, end.Y) <= point.Y
          && point.Y                  <= std::max(start.Y, end.Y);
    }

    /**
     *  Whether two closed segments share a point, decided exactly from the sides each end lies on.
     *
     *  @param  one   One segment.
     *  @param  other The other.
     *  @return       Whether they do.
     */
    auto Touch(const Line &one, const Line &other) -> bool
    {
      const Orientation first  = Orient(one.Start, one.End, other.Start);
      const Orientation second = Orient(one.Start, one.End, other.End);
      const Orientation third  = Orient(other.Start, other.End, one.Start);
      const Orientation fourth = Orient(other.Start, other.End, one.End);
      if (first != second && third != fourth)
      {
        return true;
      }
      return (first == Orientation::Collinear && Between(one.Start, one.End, other.Start))
          || (second == Orientation::Collinear && Between(one.Start, one.End, other.End))
          || (third == Orientation::Collinear && Between(other.Start, other.End, one.Start))
          || (fourth == Orientation::Collinear && Between(other.Start, other.End, one.End));
    }

    /**
     *  Whether a line meets a box that is not empty, by an exact reference: an end within it, or a
     *  point shared with one of its edges.
     *
     *  @param  line The line.
     *  @param  box  The box, of finite edges.
     *  @return      Whether they meet.
     */
    auto Reference(const Line &line, const Box &box) -> bool
    {
      const Point               low   = box.Min;
      const Point               right = {.X = box.Max.X, .Y = box.Min.Y};
      const Point               high  = box.Max;
      const Point               left  = {.X = box.Min.X, .Y = box.Max.Y};
      const std::array<Line, 4> edges = {
          Line{.Start = low  , .End = right},
          Line{.Start = right, .End = high },
          Line{.Start = high, .End = left}, Line{.Start = left, .End = low}};
      return Contains(box, line.Start)
          || Contains(box, line.End)
          || std::ranges::any_of(edges,
                                 [&line](const Line &edge) -> bool { return Touch(line, edge); });
    }
  } // namespace

  /**
   *  Lines and boxes of small whole numbers, which touch, run along edges and pass through corners
   *  often, and of any doubles: every line meets a box as the exact reference says.
   */
  TEST_CASE("contain: a line meets a box as its ends and the box's edges say, exactly")
  {
    Draws     draws{0x91C4};
    int       differ = 0;
    int       met    = 0;
    const int rounds = Rounds(20'000);
    for (int round = 0; round < rounds; ++round)
    {
      const bool whole = round % 2 == 0;
      const auto pick  = [&draws, whole]() -> double
      {
        return whole ? std::round(draws.Near(0, 4.5)) : draws.Near(0, 4);
      };
      const Line   line  = {.Start = {.X = pick(), .Y = pick()}, .End = {.X = pick(), .Y = pick()}};
      const double one   = pick();
      const double two   = pick();
      const double three = pick();
      const double four  = pick();
      const Box    box   = {.Min = {.X = std::min(one, two), .Y = std::min(three, four)},
                            .Max = {.X = std::max(one, two), .Y = std::max(three, four)}};
      const bool   meets = Intersects(line, box);
      differ += meets == Reference(line, box) ? 0 : 1;
      met    += meets ? 1 : 0;
    }
    CHECK(differ == 0);
    CHECK(met > rounds / 8);
    CHECK(met < rounds - (rounds / 8));
  }

  TEST_CASE("contain: lines at a box's corners and edges, and boxes empty, lost or endless")
  {
    const Box    unit    = {.Min = {.X = 0, .Y = 0}, .Max = {.X = 1, .Y = 1}};
    const double none    = std::numeric_limits<double>::quiet_NaN();
    const double endless = std::numeric_limits<double>::infinity();
    const double past    = std::nextafter(1.0, 2.0);
    const double over    = 2 + 0x1p-50;
    CHECK(Intersects(Line{.Start = {.X = 0, .Y = 2}, .End = {.X = 2, .Y = 0}}, unit));
    CHECK(!Intersects(Line{.Start = {.X = 0, .Y = over}, .End = {.X = over, .Y = 0}}, unit));
    CHECK(Intersects(Line{.Start = {.X = 1, .Y = -1}, .End = {.X = 1, .Y = 2}}, unit));
    CHECK(!Intersects(Line{.Start = {.X = past, .Y = -1}, .End = {.X = past, .Y = 2}}, unit));
    CHECK(Intersects(Line{.Start = {.X = 0.5, .Y = 0.5}, .End = {.X = 0.5, .Y = 0.5}}, unit));
    CHECK(!Intersects(Line{.Start = {.X = 2, .Y = 2}, .End = {.X = 2, .Y = 2}}, unit));
    CHECK(Intersects(Line{.Start = {.X = -1, .Y = 1.5}, .End = {.X = 0.5, .Y = -0.5}}, unit));
    CHECK(!Intersects(Line{.Start = {.X = -1, .Y = 0.5}, .End = {.X = 0.5, .Y = -1.0}}, unit));
    CHECK(!Intersects(Line{.Start = {.X = 0.5, .Y = 3}, .End = {.X = 3, .Y = 0.5}}, unit));

    /**
     *  An empty box, a box of a number that is not a number, and one running without end.
     */
    const Line across = {.Start = {.X = -1, .Y = -2}, .End = {.X = 2, .Y = 1}};
    CHECK(!Intersects(across, Box{}));
    CHECK(!Intersects(across, {.Min = {.X = 1, .Y = 0}, .Max = {.X = 0, .Y = 1}}));
    CHECK(!Intersects(across, {.Min = {.X = none, .Y = 0}, .Max = {.X = 1, .Y = 1}}));
    CHECK(!Intersects(across, {.Min = {.X = -endless, .Y = 0}, .Max = {.X = 0, .Y = endless}}));
    CHECK(Intersects(across, {.Min = {.X = 0, .Y = -endless}, .Max = {.X = endless, .Y = 0}}));
    CHECK(Intersects(across,
                     {.Min = {.X = -endless, .Y = -endless}, .Max = {.X = endless, .Y = endless}}));
  }

  TEST_CASE("contain: curves within a box, across its edges, near it and far from it")
  {
    const Box    unit  = {.Min = {.X = 0, .Y = 0}, .Max = {.X = 1, .Y = 1}};
    const Box    whole = {.Min = {.X = -std::numeric_limits<double>::infinity(),
                                  .Y = -std::numeric_limits<double>::infinity()},
                          .Max = {.X = std::numeric_limits<double>::infinity(),
                                  .Y = std::numeric_limits<double>::infinity()}};
    const double past  = std::nextafter(1.0, 2.0);
    const Curve  inner = Line{.Start = {.X = 0, .Y = 0.5}, .End = {.X = 1, .Y = 1}};
    const Curve  poked = Line{.Start = {.X = 0, .Y = 0.5}, .End = {.X = past, .Y = 1}};
    CHECK(Inside(inner, unit));
    CHECK(!Inside(poked, unit));
    CHECK(Intersects(poked, unit));

    const Curve small  = Circle({.X = 0.5, .Y = 0.5}, 0.25);
    const Curve large  = Circle({.X = 0.5, .Y = 0.5}, 0.6);
    const Curve around = Circle({.X = 0.5, .Y = 0.5}, 0.75);
    const Curve corner =
        Arc{.Center = {.X = 0, .Y = 0}, .Radius = 1, .Start = Degrees(0), .Sweep = Degrees(90)};
    const Curve away = Circle({.X = 5, .Y = 5}, 1);
    const Box   nook = {.Min = {.X = 0.8, .Y = 0.8}, .Max = {.X = 1, .Y = 1}};
    CHECK(Inside(small, unit));
    CHECK(Intersects(small, unit));
    CHECK(!Inside(large, unit));
    CHECK(Intersects(large, unit));
    CHECK(!Intersects(around, unit));
    CHECK(Intersects(corner, unit));
    CHECK(!Intersects(corner, nook));
    CHECK(!Intersects(away, unit));
    CHECK(!Inside(away, unit));

    const Curve bent = Bezier::Cubic{.Start = {.X = -1, .Y = 0.5},
                                     .Out   = {.X = 0, .Y = 2},
                                     .In    = {.X = 1, .Y = -1},
                                     .End   = {.X = 2, .Y = 0.5}};
    const Curve none = Polyline{};
    const Curve ray  = Ray{.Start = {.X = -1, .Y = 0.5}, .Way = {.X = 1, .Y = 0}};
    const Curve wide = Xline{.Through = {.X = 0.5, .Y = 5}, .Way = {.X = 0, .Y = 1}};
    CHECK(Intersects(bent, unit));
    CHECK(!Inside(bent, unit));
    CHECK(!Inside(none, unit));
    CHECK(!Intersects(none, unit));
    CHECK(!Inside(none, whole));
    CHECK(Intersects(ray, unit));
    CHECK(!Inside(ray, unit));
    CHECK(Inside(ray, whole));
    CHECK(Intersects(wide, unit));
    CHECK(!Intersects(wide, nook));
    CHECK(!Intersects(small, Box{}));
    CHECK(!Inside(small, Box{}));
  }

  TEST_CASE("contain: within or across a polygon, a polygon with a bulge, and a region with a hole")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Polyline  triangle  = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                              {.Position = {.X = 10, .Y = 0}},
                                              {.Position = {.X = 5, .Y = 8}}},
                                 .Closed   = true};
    const Curve     inner     = Line{.Start = {.X = 3, .Y = 1}, .End = {.X = 7, .Y = 1}};
    const Curve     across    = Line{.Start = {.X = 3, .Y = 1}, .End = {.X = 12, .Y = 1}};
    const Curve     outer     = Line{.Start = {.X = 11, .Y = 5}, .End = {.X = 12, .Y = 6}};
    const Curve     edge      = Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const Curve     round     = Circle({.X = 5, .Y = 3}, 1);
    const Curve     around    = Circle({.X = 5, .Y = 3}, 20);
    const Curve     outward   = Ray{.Start = {.X = 5, .Y = 3}, .Way = {.X = 1, .Y = 0}};
    const Curve     endless   = Xline{.Through = {.X = -5, .Y = -5}, .Way = {.X = 1, .Y = 1}};
    const Curve     none      = Polyline{};
    CHECK(Inside(inner, triangle, tolerance));
    CHECK(Intersects(inner, triangle, tolerance));
    CHECK(!Inside(across, triangle, tolerance));
    CHECK(Intersects(across, triangle, tolerance));
    CHECK(!Inside(outer, triangle, tolerance));
    CHECK(!Intersects(outer, triangle, tolerance));
    CHECK(Inside(edge, triangle, tolerance));
    CHECK(Inside(round, triangle, tolerance));
    CHECK(!Inside(around, triangle, tolerance));
    CHECK(!Intersects(around, triangle, tolerance));
    CHECK(!Inside(outward, triangle, tolerance));
    CHECK(Intersects(outward, triangle, tolerance));
    CHECK(Intersects(endless, triangle, tolerance));
    CHECK(!Inside(none, triangle, tolerance));
    CHECK(!Intersects(none, triangle, tolerance));
    CHECK(!Intersects(Curve{Line{.End = {.X = std::numeric_limits<double>::quiet_NaN()}}}, triangle,
                      tolerance));
    CHECK(!Inside(Curve{Line{.End = {.X = std::numeric_limits<double>::quiet_NaN()}}}, triangle,
                  tolerance));
    CHECK(!Intersects(Curve{Xline{.Through = {.X = 50, .Y = -5}, .Way = {.X = 0, .Y = 1}}},
                      triangle, tolerance));
    CHECK(!Intersects(Curve{Xline{.Through = {.X = 9.9, .Y = 7.9}, .Way = {.X = 1, .Y = -1}}},
                      triangle, tolerance));
    CHECK(!Inside(inner, Region{}, tolerance));
    CHECK(!Intersects(inner, Region{}, tolerance));

    const Polyline bulged = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                                          {.Position = {.X = 10, .Y = 0}},
                                          {.Position = {.X = 10, .Y = 4}},
                                          {.Position = {.X = 0, .Y = 4}}},
                             .Closed   = true};
    CHECK(Inside(Line{.Start = {.X = 5, .Y = -4}, .End = {.X = 5, .Y = 3}}, bulged, tolerance));
    CHECK(!Inside(Line{.Start = {.X = 5, .Y = -6}, .End = {.X = 5, .Y = 3}}, bulged, tolerance));

    const Region holed = {
        .Loops = {Loop{.Edges = {Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                       {.Position = {.X = 10, .Y = 0}},
                                                       {.Position = {.X = 10, .Y = 10}},
                                                       {.Position = {.X = 0, .Y = 10}}},
                                          .Closed   = true}}},
                  Loop{.Edges = {Arc{.Center = {.X = 5, .Y = 5}, .Radius = 2, .Sweep = -Turn}}}}};
    CHECK(Inside(Line{.Start = {.X = 1, .Y = 1}, .End = {.X = 2, .Y = 1}}, holed, tolerance));
    CHECK(!Inside(Line{.Start = {.X = 1, .Y = 5}, .End = {.X = 9, .Y = 5}}, holed, tolerance));
    CHECK(Intersects(Line{.Start = {.X = 1, .Y = 5}, .End = {.X = 9, .Y = 5}}, holed, tolerance));
    CHECK(!Intersects(Circle({.X = 5, .Y = 5}, 1), holed, tolerance));
    CHECK(Intersects(Circle({.X = 5, .Y = 5}, 3), holed, tolerance));
  }
} // namespace Geoxy
