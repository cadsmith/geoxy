/***************************************************************************************************
 *  @file      tests/geoxy/construct.cpp
 *  @brief     Tests of arcs, circles and lines built from what they pass through or touch: against
 *             exact rational circumcenters, mpmath's Apollonius circles, and every circle's
 *             tangency to what it was built on.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/construct.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether two points are within a distance.
     *
     *  @param  lhs      One point.
     *  @param  rhs      The other.
     *  @param  distance The distance.
     *  @return          Whether they are.
     */
    auto Closeby(Point lhs, Point rhs, double distance) -> bool
    {
      return Distance(lhs, rhs) <= distance;
    }

    /**
     *  Whether a circle touches a line or a circle where it says, as it says, within a part of the
     *  scale.
     *
     *  @param  segment The line, taken whole, or the arc, taken as its circle.
     *  @param  circle  The circle built.
     *  @param  touch   Where it touches and how it lies.
     *  @return         Whether it does.
     */
    auto Touches(const Segment &segment, const Arc &circle, const Grazing &touch) -> bool
    {
      const double scale = 1e-12 * std::max({1.0, circle.Radius, Length(circle.Center - Point{})});
      const bool   onto  = std::abs(Distance(touch.At, circle.Center) - circle.Radius) <= scale;
      if (const Line *line = std::get_if<Line>(&segment))
      {
        const Vector way    = *Unit(line->End - line->Start);
        const double offset = Cross(way, circle.Center - line->Start);
        const bool   left   = touch.Kind == Qualifier::Left;
        return onto
            && std::abs(Cross(way, touch.At - line->Start)) <= scale
            && std::abs(std::abs(offset) - circle.Radius) <= scale
            && (left ? offset > 0 : (touch.Kind == Qualifier::Right && offset < 0));
      }
      const Arc   &guide  = std::get<Arc>(segment);
      const double radius = std::abs(guide.Radius);
      const double apart  = Distance(circle.Center, guide.Center);
      const bool   upon   = std::abs(Distance(touch.At, guide.Center) - radius) <= scale;
      switch (touch.Kind)
      {
      case Qualifier::Outside:
        return onto && upon && std::abs(apart - (radius + circle.Radius)) <= scale;
      case Qualifier::Enclosed:
        return onto
            && upon
            && circle.Radius < radius
            && std::abs(apart - (radius - circle.Radius)) <= scale;
      case Qualifier::Enclosing:
        return onto
            && upon
            && circle.Radius > radius
            && std::abs(apart - (circle.Radius - radius)) <= scale;
      default:
        return false;
      }
    }

    /**
     *  Whether every circle built touches every curve it was built on, as it says.
     *
     *  @tparam CountT  How many curves.
     *  @param  curves  The curves.
     *  @param  circles The circles.
     *  @return         Whether they do.
     */
    template <std::size_t CountT>
    auto TouchesAll(const std::array<Segment, CountT>     &curves,
                    const std::vector<Tangential<CountT>> &circles) -> bool
    {
      bool                           tangent = true;
      const std::span<const Segment> guides  = curves;
      for (const Tangential<CountT> &circle : circles)
      {
        const std::span<const Grazing> touches = circle.Touches;
        for (std::size_t index = 0; index < CountT; ++index)
        {
          tangent = tangent && Touches(guides[index], circle.Circle, touches[index]);
        }
      }
      return tangent;
    }

    /**
     *  Whether a line is tangent to a circle: the circle's center its radius from the line's
     *  carrier, the line touching it at its start or its end.
     *
     *  @param  line   The line.
     *  @param  circle The circle.
     *  @return        Whether it is.
     */
    auto Touches(const Line &line, const Arc &circle) -> bool
    {
      const double scale = 4e-15 * std::max({1.0, circle.Radius, Length(circle.Center - Point{})});
      const double start = std::abs(Distance(line.Start, circle.Center) - circle.Radius);
      const double end   = std::abs(Distance(line.End, circle.Center) - circle.Radius);
      const Point  touch = start < end ? line.Start : line.End;
      const Vector way   = Perpendicular(touch - circle.Center);
      const Vector along = line.End - line.Start;
      return std::min(start, end)               <= scale
          && std::abs(Cross(*Unit(way), along)) <= scale * std::max(1.0, Length(along));
    }
  } // namespace

  TEST_CASE("construct: an arc through three points, its center the exact circumcenter's double")
  {
    const std::optional<Arc> arc =
        Arc::Through({.X = 0, .Y = 0}, {.X = 4, .Y = 2}, {.X = 1, .Y = 5});
    REQUIRE(arc);
    CHECK(Closeby(arc->Center, {.X = 1.3333333333333333, .Y = 2.3333333333333335}, 1e-15));
    CHECK(arc->Sweep.InRadians() > 0);
    CHECK(Closeby(At(*arc, 0), {.X = 0, .Y = 0}, 1e-15));
    CHECK(Closeby(At(*arc, 1), {.X = 1, .Y = 5}, 1e-14));

    const std::optional<Arc> far =
        Arc::Through({.X = 5e6 + 0.25, .Y = -3e6}, {.X = 5e6 + 4.5, .Y = -3e6 + 2},
                     {.X = 5e6 + 1, .Y = -3e6 + 5.75});
    REQUIRE(far);
    CHECK(Closeby(far->Center, {.X = 5000001.549386921, .Y = -2999997.245572207}, 4e-9));

    const std::optional<Arc> down =
        Arc::Through({.X = 2, .Y = -3}, {.X = 5, .Y = 2}, {.X = 8, .Y = -3});
    REQUIRE(down);
    CHECK(down->Center == Point{.X = 5, .Y = -1.4});
    CHECK(down->Sweep.InRadians() < 0);
    CHECK(Closeby(At(*down, 0.5), {.X = 5, .Y = 2}, 1e-14));

    const std::optional<Arc> clockwise = Arc::Through(
        {.X = 0, .Y = 1}, {.X = std::sqrt(0.5), .Y = std::sqrt(0.5)}, {.X = 1, .Y = 0});
    REQUIRE(clockwise);
    CHECK(clockwise->Sweep.InRadians() == doctest::Approx(-Pi / 2));

    CHECK_FALSE(Arc::Through({.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 3, .Y = 3}));
    CHECK_FALSE(Arc::Through({.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 0, .Y = 0}));
    CHECK_FALSE(Arc::Through({.X = 0, .Y = 0}, {.X = 1e-200, .Y = 0}, {.X = 0, .Y = 1e-200}));
    CHECK_FALSE(Arc::Through({.X = 0, .Y = 0}, {.X = 1e300, .Y = 1}, {.X = 2e300, .Y = 0}));
  }

  TEST_CASE("construct: an arc by its start, center and end, counterclockwise")
  {
    const std::optional<Arc> quarter = Arc::Around({.X = 2, .Y = 0}, {}, {.X = 0, .Y = 5});
    REQUIRE(quarter);
    CHECK(quarter->Radius            == 2);
    CHECK(quarter->Sweep.InRadians() == doctest::Approx(Pi / 2));
    CHECK(Arc::Around({.X = 2, .Y = 0}, {}, {.X = 0, .Y = -5})->Sweep.InRadians()
          == doctest::Approx(3 * Pi / 2));
    CHECK(Arc::Around({.X = 2, .Y = 0}, {}, {.X = 5, .Y = 0})->Sweep.InRadians()
          == doctest::Approx(2 * Pi));
    CHECK_FALSE(Arc::Around({}, {}, {.X = 1, .Y = 0}));
    CHECK_FALSE(Arc::Around({.X = 1, .Y = 0}, {}, {}));
  }

  TEST_CASE("construct: an arc by its start, end and radius, minor or major")
  {
    const std::optional<Arc> half = Arc::Spanning({.X = 0, .Y = 0}, {.X = 2, .Y = 0}, 1);
    REQUIRE(half);
    CHECK(half->Center            == Point{.X = 1, .Y = 0});
    CHECK(half->Sweep.InRadians() == doctest::Approx(Pi));
    CHECK(Closeby(At(*half, 0.5), {.X = 1, .Y = -1}, 1e-15));

    const std::optional<Arc> minor = Arc::Spanning({.X = 0, .Y = 0}, {.X = 2, .Y = 0}, 2);
    REQUIRE(minor);
    CHECK(Closeby(minor->Center, {.X = 1, .Y = std::sqrt(3.0)}, 1e-15));
    CHECK(minor->Sweep.InRadians() == doctest::Approx(Pi / 3));
    const std::optional<Arc> major = Arc::Spanning({.X = 0, .Y = 0}, {.X = 2, .Y = 0}, -2);
    REQUIRE(major);
    CHECK(major->Radius == 2);
    CHECK(Closeby(major->Center, {.X = 1, .Y = -std::sqrt(3.0)}, 1e-15));
    CHECK(major->Sweep.InRadians() == doctest::Approx(5 * Pi / 3));
    CHECK(Closeby(At(*major, 1), {.X = 2, .Y = 0}, 1e-15));

    CHECK_FALSE(Arc::Spanning({.X = 0, .Y = 0}, {.X = 2, .Y = 0}, 0.5));
    CHECK_FALSE(Arc::Spanning({.X = 1, .Y = 1}, {.X = 1, .Y = 1}, 3));
  }

  TEST_CASE("construct: an arc by its start, end and the way it leaves its start")
  {
    const std::optional<Arc> right =
        Arc::Leaving({.X = 0, .Y = 0}, {.X = 2, .Y = 0}, {.X = 0, .Y = 3});
    REQUIRE(right);
    CHECK(right->Center            == Point{.X = 1, .Y = 0});
    CHECK(right->Radius            == 1);
    CHECK(right->Sweep.InRadians() == doctest::Approx(-Pi));
    CHECK(Closeby(Point{} + *Unit(Derivative(*right, 0)), {.X = 0, .Y = 1}, 1e-15));

    const std::optional<Arc> left =
        Arc::Leaving({.X = 5e6, .Y = 1}, {.X = 5e6 + 3, .Y = 4}, {.X = 1, .Y = 0});
    REQUIRE(left);
    CHECK(left->Sweep.InRadians() > 0);
    CHECK(Closeby(At(*left, 1), {.X = 5e6 + 3, .Y = 4}, 1e-8));
    CHECK(Closeby(Point{} + *Unit(Derivative(*left, 0)), {.X = 1, .Y = 0}, 1e-15));

    CHECK_FALSE(Arc::Leaving({.X = 0, .Y = 0}, {.X = 2, .Y = 0}, {.X = 0, .Y = 0}));
    CHECK_FALSE(Arc::Leaving({.X = 0, .Y = 0}, {.X = 2, .Y = 0}, {.X = -1, .Y = 0}));
    CHECK_FALSE(Arc::Leaving({.X = 0, .Y = 0}, {.X = 1e300, .Y = 1e300}, {.X = 1, .Y = 0}));
  }

  TEST_CASE("construct: circles on a diameter and through three points")
  {
    CHECK(Circle(Point{.X = 1, .Y = 1}, Point{.X = 5, .Y = 4})
          == Arc{.Center = {.X = 3, .Y = 2.5}, .Radius = 2.5, .Start = {}, .Sweep = Turn});
    const std::optional<Arc> three = Circle({.X = 2, .Y = -3}, {.X = 5, .Y = 2}, {.X = 8, .Y = -3});
    REQUIRE(three);
    CHECK(three->Center == Point{.X = 5, .Y = -1.4});
    CHECK(three->Sweep  == Turn);
    CHECK_FALSE(Circle({.X = 0, .Y = 0}, {.X = 1, .Y = 0}, {.X = 2, .Y = 0}));
  }

  TEST_CASE("construct: circles of a radius tangent to two lines, a line and a circle, two circles")
  {
    const std::array<Segment, 2> lines = {
        Line{.Start = {.X = -5, .Y = 0}, .End = {.X = 5, .Y = 0}},
        Line{.Start = {.X = 5e-1, .Y = -4}, .End = {.X = 3, .Y = 6}}};
    const std::vector<Tangential<2>> four = TangentCircles(lines[0], lines[1], 1.5);
    CHECK(four.size() == 4);
    CHECK(TouchesAll(lines, four));

    const std::array<Segment, 2> crossing = {
        Line{.Start = {.X = -10, .Y = 0}, .End = {.X = 10, .Y = 0}}, Circle({.X = 0, .Y = 0.5}, 2)};
    const std::vector<Tangential<2>> eight = TangentCircles(crossing[0], crossing[1], 0.5);
    CHECK(eight.size() == 8);
    CHECK(TouchesAll(crossing, eight));
    const std::array<Segment, 2> swapped = {crossing[1], crossing[0]};
    CHECK(TouchesAll(swapped, TangentCircles(swapped[0], swapped[1], 0.5)));

    const std::array<Segment, 2>     apart = {Circle({.X = 5e6, .Y = 0}, 1),
                                              Circle({.X = 5e6 + 5, .Y = 0}, -1.5)};
    const std::vector<Tangential<2>> both  = TangentCircles(apart[0], apart[1], 4);
    CHECK(both.size() == 8);
    CHECK(TouchesAll(apart, both));
    CHECK(both[0].Circle.Radius == 4);

    const std::array<Segment, 2> kissing = {
        Line{.Start = {.X = -5, .Y = 0}, .End = {.X = 5, .Y = 0}}, Circle({.X = 0, .Y = 3}, 1)};
    const std::vector<Tangential<2>> one = TangentCircles(kissing[0], kissing[1], 1);
    REQUIRE(one.size()           == 1);
    CHECK(one[0].Circle.Center   == Point{.X = 0, .Y = 1});
    CHECK(one[0].Touches[0].Kind == Qualifier::Left);
    CHECK(one[0].Touches[1].Kind == Qualifier::Outside);

    const std::array<Segment, 2> itself = {
        Circle({.X = 0, .Y = 0}, 2), Line{.Start = {.X = -5, .Y = -2}, .End = {.X = 5, .Y = -2}}};
    const std::vector<Tangential<2>> three = TangentCircles(itself[0], itself[1], 2);
    CHECK(three.size() == 3);
    CHECK(TouchesAll(itself, three));

    const std::array<Segment, 2> point = {Line{.Start = {.X = -5, .Y = 0}, .End = {.X = 5, .Y = 0}},
                                          Circle({.X = 0, .Y = 2}, 0)};
    const std::vector<Tangential<2>> through = TangentCircles(point[0], point[1], 1.5);
    REQUIRE(through.size() == 2);
    CHECK(Closeby(through[0].Touches[1].At, {.X = 0, .Y = 2}, 1e-15));

    CHECK(TangentCircles(lines[0], Line{.Start = {.X = 0, .Y = 3}, .End = {.X = 1, .Y = 3}}, 1.5)
              .empty());
    CHECK(TangentCircles(lines[0], lines[1], 0).empty());
    CHECK(TangentCircles(lines[0], Line{.Start = {}, .End = {.X = 1, .Y = 1e-300}}, 1e10).size()
          == 2);
    CHECK(TangentCircles(Circle({}, 1), Circle({}, 3), 1).empty());
    CHECK(TangentCircles(lines[0], Line{.Start = {.X = 1, .Y = 1}, .End = {.X = 1, .Y = 1}}, 1)
              .empty());
  }

  TEST_CASE("construct: circles tangent to three curves, Apollonius's eight against mpmath")
  {
    const std::array<Segment, 3> circles = {
        Circle({.X = 0, .Y = 0}, 1), Circle({.X = 4, .Y = 0}, 1.5), Circle({.X = 1, .Y = 5}, 0.5)};
    std::vector<Tangential<3>> found = TangentCircles(circles[0], circles[1], circles[2]);
    REQUIRE(found.size() == 8);
    CHECK(TouchesAll(circles, found));
    std::ranges::sort(found, {},
                      [](const Tangential<3> &circle) -> double { return circle.Circle.Radius; });
    constexpr std::array<std::array<double, 3>, 8> Mpmath   = {{
        {1.5921375869193368855, 2.5578624130806631145, 2.0128993046453049163},
        {1.5354432727998168274, 3.1078474907204762488, 2.466453817601465381},
        {0.084070433626922018482, 1.8135397214155381652, 2.8154873061969247704},
        {-0.28661128487822326048, 2.3914644513951289304, 3.4085780558051572168},
        {2.288881051524507497, 1.1489092660362805079, 3.5610484121960599756},
        {4.1685414696831513919, 2.2132583412126739443, 3.719666351493042227},
        {2.3385075743709856952, 1.8114924256290143048, 3.9580605949678855614},
        {4.8041585080577432083, 3.1351643822561680983, 4.7366536128923891333},
    }};
    const std::span<const std::array<double, 3>>   expected = Mpmath;
    for (std::size_t index = 0; index < 8; ++index)
    {
      const std::array<double, 3> &circle = expected[index];
      CHECK(Closeby(found[index].Circle.Center, {.X = circle[0], .Y = circle[1]}, 1e-14));
      CHECK(found[index].Circle.Radius == doctest::Approx(circle[2]).epsilon(1e-14));
    }

    const std::array<Segment, 3> triangle = {
        Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}},
        Line{.Start = {.X = 4, .Y = 0}, .End = {.X = 0, .Y = 3}},
        Line{.Start = {.X = 0, .Y = 3}, .End = {.X = 0, .Y = 0}}};
    std::vector<Tangential<3>> inscribed = TangentCircles(triangle[0], triangle[1], triangle[2]);
    REQUIRE(inscribed.size() == 4);
    CHECK(TouchesAll(triangle, inscribed));
    std::ranges::sort(inscribed, {},
                      [](const Tangential<3> &circle) -> double { return circle.Circle.Radius; });
    CHECK(Closeby(inscribed[0].Circle.Center, {.X = 1, .Y = 1}, 1e-15));
    CHECK(inscribed[0].Circle.Radius == doctest::Approx(1));
    CHECK(inscribed[1].Circle.Radius == doctest::Approx(2));
    CHECK(inscribed[2].Circle.Radius == doctest::Approx(3));
    CHECK(inscribed[3].Circle.Radius == doctest::Approx(6));

    const std::array<Segment, 3> parallel = {
        Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}},
        Line{.Start = {.X = 0, .Y = 2}, .End = {.X = 4, .Y = 2}},
        Line{.Start = {.X = 1, .Y = -1}, .End = {.X = 2, .Y = 3}}};
    const std::vector<Tangential<3>> two = TangentCircles(parallel[0], parallel[1], parallel[2]);
    CHECK(two.size() == 2);
    CHECK(TouchesAll(parallel, two));

    const std::array<Segment, 3> mixed = {
        Line{.Start = {.X = -5, .Y = 0}, .End = {.X = 5, .Y = 0}}, Circle({.X = 1, .Y = 3}, 1),
        Line{.Start = {.X = -1, .Y = -5}, .End = {.X = 2, .Y = 5}}};
    std::vector<Tangential<3>> some = TangentCircles(mixed[0], mixed[1], mixed[2]);
    REQUIRE(some.size() == 4);
    CHECK(TouchesAll(mixed, some));
    std::ranges::sort(some, {},
                      [](const Tangential<3> &circle) -> double { return circle.Circle.Radius; });
    constexpr std::array<std::array<double, 3>, 4> Mixed  = {{
        {2.01820831662260528, 1.12959352200492995, 1.12959352200492995},
        {-0.435745509301551904, 1.25767064593494658, 1.25767064593494658},
        {5.93403689050583486, 4.04309000460906098, 4.04309000460906098},
        {-8.31649969782688824, 11.8496458274510625, 11.8496458274510625},
    }};
    const std::span<const std::array<double, 3>>   mpmath = Mixed;
    for (std::size_t index = 0; index < 4; ++index)
    {
      const std::array<double, 3> &circle = mpmath[index];
      CHECK(Closeby(some[index].Circle.Center, {.X = circle[0], .Y = circle[1]}, 1e-13));
      CHECK(some[index].Circle.Radius == doctest::Approx(circle[2]).epsilon(1e-14));
    }
    const std::array<Segment, 3> rounder = {
        Circle({.X = 5e6, .Y = 5e6}, 2),
        Line{.Start = {.X = 5e6 - 9, .Y = 5e6 - 3}, .End = {.X = 5e6 + 9, .Y = 5e6 - 3}},
        Circle({.X = 5e6 + 6, .Y = 5e6 + 1}, 1.5)};
    const std::vector<Tangential<3>> more = TangentCircles(rounder[0], rounder[1], rounder[2]);
    CHECK(more.size() == 8);
    CHECK(TouchesAll(rounder, more));

    const std::array<Segment, 3> cone = {Circle({.X = 0, .Y = 0}, 1), Circle({.X = 2, .Y = 0}, 1),
                                         Circle({.X = 0, .Y = 2}, 3)};
    const std::vector<Tangential<3>> held = TangentCircles(cone[0], cone[1], cone[2]);
    CHECK(held.size() == 3);
    CHECK(TouchesAll(cone, held));
    CHECK(TangentCircles(triangle[0], triangle[1],
                         Line{.Start = {.X = 1, .Y = 1}, .End = {.X = 1, .Y = 1}})
              .empty());
  }

  TEST_CASE("construct: lines from a point tangent to a circle, and tangent to two circles")
  {
    const Arc               circle = Circle({}, 3);
    const std::vector<Line> lines  = TangentLines({.X = 5, .Y = 0}, circle);
    REQUIRE(lines.size() == 2);
    CHECK(Closeby(lines[0].End, {.X = 1.8, .Y = 2.4}, 1e-15));
    CHECK(Closeby(lines[1].End, {.X = 1.8, .Y = -2.4}, 1e-15));
    CHECK(TangentLines({.X = 1, .Y = 0}, circle).empty());
    CHECK(TangentLines(Point{}, Circle({}, 0)).empty());
    CHECK(TangentLines({.X = 3, .Y = 0},
                       Arc{.Center = {}, .Radius = -3, .Start = {}, .Sweep = Degrees(10)})
          == std::vector<Line>{{.Start = {.X = 3, .Y = 0}, .End = {.X = 3, .Y = 0}}});
    for (const Line &line : TangentLines({.X = 5e6 + 2, .Y = -7}, Circle({.X = 5e6, .Y = 1}, 3)))
    {
      CHECK(Touches(line, Circle({.X = 5e6, .Y = 1}, 3)));
    }

    const Arc               first  = Circle({.X = 5e6, .Y = 1}, 1);
    const Arc               second = Circle({.X = 5e6 + 6, .Y = 3}, 2);
    const std::vector<Line> four   = TangentLines(first, second);
    REQUIRE(four.size() == 4);
    for (const Line &line : four)
    {
      CHECK(Touches(line, first));
      CHECK(Touches(line, second));
    }
    const double outer = Cross(four[0].End - four[0].Start, first.Center - four[0].Start);
    CHECK(outer * Cross(four[0].End - four[0].Start, second.Center - four[0].Start) > 0);
    const double inner = Cross(four[2].End - four[2].Start, first.Center - four[2].Start);
    CHECK(inner * Cross(four[2].End - four[2].Start, second.Center - four[2].Start) < 0);

    CHECK(TangentLines(Circle({}, 2), Circle({.X = 3, .Y = 0}, 2)).size() == 2);
    CHECK(TangentLines(Circle({}, 5), Circle({.X = 1, .Y = 0}, 1)).empty());
    CHECK(TangentLines(Circle({}, 5), Circle({}, 1)).empty());
    const std::vector<Line> kissing = TangentLines(Circle({}, 1), Circle({.X = 3, .Y = 0}, 2));
    REQUIRE(kissing.size() == 4);
    CHECK(Closeby(kissing[2].Start, {.X = 1, .Y = 0}, 1e-15));
    CHECK(Closeby(kissing[2].End, {.X = 1, .Y = 0}, 1e-15));
    CHECK(kissing[2] == kissing[3]);
  }
} // namespace Geoxy
