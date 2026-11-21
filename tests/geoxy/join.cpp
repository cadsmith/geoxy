/***************************************************************************************************
 *  @file      tests/geoxy/join.cpp
 *  @brief     Tests of curves linked end to end: into a polyline where they are lines and arcs,
 *             else into a spline, each turned as it must be, gaps refused.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/join.hpp"

#include "curves.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The tolerance the tests link within.
     */
    constexpr Tolerance LinkFuzz = {.Length = 1e-9, .Radians = 1e-9};

    /**
     *  A curve linked into a polyline.
     *
     *  @param  curves The curves.
     *  @return        The polyline, or an empty one where they do not link into one.
     */
    auto Polygonal(const std::vector<Curve> &curves) -> Polyline
    {
      const std::optional<Curve> linked = Join(curves, LinkFuzz);
      return linked && std::holds_alternative<Polyline>(*linked) ? std::get<Polyline>(*linked)
                                                                 : Polyline{};
    }

    /**
     *  A curve linked into a spline.
     *
     *  @param  curves The curves.
     *  @return        The spline, or an empty one where they do not link into one.
     */
    auto Splined(const std::vector<Curve> &curves) -> Spline
    {
      const std::optional<Curve> linked = Join(curves, LinkFuzz);
      return linked && std::holds_alternative<Spline>(*linked) ? std::get<Spline>(*linked)
                                                               : Spline{};
    }

    /**
     *  Whether a point lies on a circle, to a few steps of a double.
     *
     *  @param  point  The point.
     *  @param  center The circle's center.
     *  @param  radius Its radius.
     *  @return        Whether it does.
     */
    auto Round(Point point, Point center, double radius) -> bool
    {
      return std::abs(Distance(point, center) - radius) <= 1e-14 * radius;
    }
  } // namespace

  TEST_CASE("join: lines end to end, each turned as it must be, into one polyline")
  {
    const Line first  = {.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}};
    const Line second = {.Start = {.X = 4, .Y = 3}, .End = {.X = 4, .Y = 1e-10}};
    const Line third  = {.Start = {.X = 4, .Y = 3}, .End = {.X = 0, .Y = 3}};

    const Polyline chain = Polygonal({first, second, third});
    REQUIRE(chain.Vertices.size()    == 4);
    CHECK(chain.Vertices[1].Position == Point{.X = 4, .Y = 0});
    CHECK(chain.Vertices[2].Position == Point{.X = 4, .Y = 3});
    CHECK(chain.Vertices[3].Position == Point{.X = 0, .Y = 3});
    CHECK_FALSE(chain.Closed);

    const Polyline turned = Polygonal({Line{.Start = first.End, .End = first.Start}, second});
    REQUIRE(turned.Vertices.size()    == 3);
    CHECK(turned.Vertices[0].Position == Point{.X = 0, .Y = 0});

    const Line     back   = {.Start = {.X = 0, .Y = 3}, .End = {.X = 0, .Y = 1e-10}};
    const Polyline closed = Polygonal({first, second, third, back});
    REQUIRE(closed.Vertices.size() == 4);
    CHECK(closed.Closed);
    CHECK(Polygonal({first}).Vertices.size() == 2);

    CHECK_FALSE(Join(std::vector<Curve>{first, third}, LinkFuzz));
    CHECK_FALSE(Join(std::vector<Curve>{third, first}, LinkFuzz));
    CHECK_FALSE(Join(std::vector<Curve>{first, second, first}, LinkFuzz));
    CHECK_FALSE(Join(std::vector<Curve>{}, LinkFuzz));
  }

  TEST_CASE("join: arcs as bulges, a whole circle as two halves, polylines as they are")
  {
    const Line first = {.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}};
    const Arc  bend  = {
        .Center = {.X = 4, .Y = 2}, .Radius = 2, .Start = Degrees(-90), .Sweep = Degrees(180)};
    const Polyline chain = Polygonal({first, bend});
    REQUIRE(chain.Vertices.size()    == 3);
    CHECK(chain.Vertices[1].Bulge    == doctest::Approx(1));
    CHECK(chain.Vertices[2].Position == At(bend, 1));

    const Polyline circle = Polygonal({Circle({.X = 1, .Y = 1}, 2)});
    REQUIRE(circle.Vertices.size() == 2);
    CHECK(circle.Closed);
    CHECK(circle.Vertices[0].Bulge == doctest::Approx(1));

    const Polyline loop   = {.Vertices = {{.Position = {.X = 4, .Y = 4}, .Bulge = 0},
                                          {.Position = {.X = 0, .Y = 4}, .Bulge = 0.3},
                                          {.Position = {.X = 0, .Y = 0}, .Bulge = 0}},
                             .Closed   = false};
    const Polyline joined = Polygonal({first, bend, loop});
    REQUIRE(joined.Vertices.size() == 4);
    CHECK(joined.Closed);
    CHECK(joined.Vertices[3].Bulge == doctest::Approx(0.3));

    const Polyline square = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                          {.Position = {.X = 1, .Y = 0}, .Bulge = 0},
                                          {.Position = {.X = 1, .Y = 1}, .Bulge = -0.2}},
                             .Closed   = true};
    CHECK(Polygonal({square}) == square);

    CHECK_FALSE(Join(std::vector<Curve>{Arc{
                         .Center = {}, .Radius = 1, .Start = Degrees(0), .Sweep = Degrees(400)}},
                     LinkFuzz));
    CHECK_FALSE(
        Join(std::vector<Curve>{Arc{.Center = {},
                                    .Radius = 1,
                                    .Start  = Degrees(0),
                                    .Sweep  = Radians(std::numeric_limits<double>::quiet_NaN())}},
             LinkFuzz));
    CHECK_FALSE(Join(std::vector<Curve>{Polyline{}}, LinkFuzz));
  }

  TEST_CASE("join: a half circle held in degrees has a bulge of 1 exactly, a quarter in radians "
            "its tangent")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Curve     half =
        Arc{.Center = {.X = 4, .Y = 1}, .Radius = 1, .Start = Degrees(-90), .Sweep = Degrees(180)};
    const auto linked = Join(std::array{half}, tolerance);
    REQUIRE(linked);
    CHECK(std::equal_to<double>{}(std::get<Polyline>(*linked).Vertices.front().Bulge, 1.0));

    const Curve whole = Circle({.X = 0, .Y = 0}, 3);
    const auto  round = Join(std::array{whole}, tolerance);
    REQUIRE(round);
    CHECK(std::equal_to<double>{}(std::get<Polyline>(*round).Vertices.front().Bulge, 1.0));

    const Curve quarter =
        Arc{.Center = {}, .Radius = 2, .Start = Radians(0), .Sweep = Radians(Pi / 2)};
    const auto turned = Join(std::array{quarter}, tolerance);
    REQUIRE(turned);
    CHECK(std::get<Polyline>(*turned).Vertices.front().Bulge == doctest::Approx(std::tan(Pi / 8)));
  }

  TEST_CASE("join: any other curve among them makes a spline, its pieces over units")
  {
    const Bezier::Cubic cubic = Tests::Wave;
    const Line          line  = {.Start = cubic.End, .End = {.X = 6, .Y = 1}};
    const Arc           bend  = {
        .Center = {.X = 6, .Y = 3}, .Radius = 2, .Start = Degrees(-90), .Sweep = Degrees(135)};
    const Spline spline = Splined({cubic, line, bend});
    REQUIRE(Valid(spline));
    CHECK(spline.Degree  == 3);
    CHECK(Domain(spline) == Range{.Start = 0, .End = 4});
    for (const double along : {0.0, 0.25, 0.5, 0.75, 1.0})
    {
      CHECK(Distance(At(spline, along), At(cubic, along)) <= 1e-14);
      CHECK(Distance(At(spline, 1 + along), At(line, along)) <= 1e-14);
      CHECK(Round(At(spline, 2 + (2 * along)), bend.Center, bend.Radius));
    }
    CHECK(At(spline, 4) == At(bend, 1));

    const Bezier::Quadratic quadratic = {
        .Start = At(bend, 1), .Control = {.X = 9, .Y = 9}, .End = {.X = 0, .Y = 0}};
    const Spline closed = Splined({quadratic, cubic, line, bend});
    REQUIRE(Valid(closed));
    CHECK(closed.Poles.back() == closed.Poles.front());

    const Spline turned = Splined(
        {Bezier::Cubic{.Start = cubic.End, .Out = cubic.In, .In = cubic.Out, .End = cubic.Start},
         line});
    REQUIRE(Valid(turned));
    CHECK(Distance(At(turned, 0.25), At(cubic, 0.25)) <= 1e-14);
  }

  TEST_CASE("join: ellipses, polylines and rational splines, their weights carried on")
  {
    const Ellipse ellipse = {.Map = Scaling(3, 1), .Start = Degrees(0), .Sweep = Degrees(360)};
    const Spline  whole   = Splined({ellipse});
    REQUIRE(Valid(whole));
    CHECK(whole.Degree      == 2);
    CHECK(Domain(whole).End == doctest::Approx(4));
    for (int step = 0; step <= 32; ++step)
    {
      const Point point = At(whole, step / 8.0);
      CHECK(std::abs((point.X * point.X / 9) + (point.Y * point.Y) - 1) <= 1e-14);
    }

    const Spline   heavy    = {.Degree  = 2,
                               .Knots   = {0, 0, 0, 1, 1, 1},
                               .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 2}, {.X = 2, .Y = 0}},
                               .Weights = {1, 1, 2}};
    const Polyline polyline = {.Vertices = {{.Position = {.X = 2, .Y = 0}, .Bulge = 0},
                                            {.Position = {.X = 3, .Y = 0}, .Bulge = 1},
                                            {.Position = {.X = 5, .Y = 0}, .Bulge = 0}},
                               .Closed   = false};
    const Spline   chain    = Splined({heavy, polyline});
    REQUIRE(Valid(chain));
    CHECK_FALSE(chain.Weights.empty());
    CHECK(Distance(At(chain, 0.5), At(heavy, 0.5)) <= 1e-15);
    CHECK(Distance(At(chain, 1.5), Point{.X = 2.5, .Y = 0}) <= 1e-15);
    CHECK(Round(At(chain, 2.5), {.X = 4, .Y = 0}, 1));

    CHECK(Valid(Splined({Tests::Uneven(), Line{.Start = {.X = 8, .Y = 1}, .End = {}}})));
  }

  TEST_CASE("join: no spline from curves with no end, swept past a turn, or apart")
  {
    const Bezier::Cubic cubic = Tests::Wave;
    CHECK_FALSE(Join(std::vector<Curve>{cubic, Ray{.Start = cubic.End, .Way = {.X = 1, .Y = 0}}},
                     LinkFuzz));
    CHECK_FALSE(Join(
        std::vector<Curve>{cubic, Xline{.Through = cubic.End, .Way = {.X = 1, .Y = 0}}}, LinkFuzz));
    CHECK_FALSE(Join(std::vector<Curve>{cubic, Spline{}}, LinkFuzz));
    CHECK_FALSE(Join(
        std::vector<Curve>{
            cubic, Arc{.Center = {}, .Radius = 1, .Start = Degrees(0), .Sweep = Degrees(400)}},
        LinkFuzz));
    CHECK_FALSE(
        Join(std::vector<Curve>{cubic,
                                Ellipse{.Map = Scaling(3, 1), .Start = {}, .Sweep = Degrees(-361)}},
             LinkFuzz));
    CHECK_FALSE(
        Join(std::vector<Curve>{cubic, Line{.Start = {.X = 9, .Y = 9}, .End = {.X = 8, .Y = 8}}},
             LinkFuzz));
  }
} // namespace Geoxy
