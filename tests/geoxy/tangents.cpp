/***************************************************************************************************
 *  @file      tests/geoxy/tangents.cpp
 *  @brief     Tests of the places of a curve whose tangent passes through a point: each tangent
 *             through it, or a vertex the curve turns about it at.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/tangents.hpp"

#include "curves.hpp"
#include "stiff.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether a curve turns about a point the other way at a parameter: its tangent's side of the
     *  way from the point differs a step either side.
     *
     *  @param  curve The curve.
     *  @param  point The point.
     *  @param  along The parameter.
     *  @return       Whether it does.
     */
    auto TurnsAbout(const Curve &curve, Point point, double along) -> bool
    {
      const Range  domain = Domain(curve);
      const double back   = along - 1e-7 < domain.Start ? domain.End - 1e-7 : along - 1e-7;
      const double before = Cross(At(curve, back) - point, Derivative(curve, back));
      const double after  = Cross(At(curve, along + 1e-7) - point, Derivative(curve, along + 1e-7));
      return (before > 0.0) != (after > 0.0);
    }

    /**
     *  Whether every tangent place of a curve has its tangent through the point, or turns about it
     *  at a vertex, the places in order and none of them the point's own.
     *
     *  @param  curve The curve.
     *  @param  point The point.
     *  @param  slack How far a place's point may lie from the curve's at its parameter: nothing,
     *                but for a spline whose weights are uneven, whose place is told as found on a
     *                piece of it and not through the parameter's rounding.
     *  @return       Whether they do.
     */
    auto Touch(const Curve &curve, Point point, double slack = 0.0) -> bool
    {
      const std::vector<Place> places = Tangents(curve, point);
      for (std::size_t index = 0; index < places.size(); ++index)
      {
        const Place &place   = places[index];
        const Vector tangent = Derivative(curve, place.Along);
        const Vector away    = point - place.At;
        const double scale   = std::max(1.0, Length(away) * Length(tangent));
        const bool   through = std::abs(Cross(away, tangent)) <= 1e-9 * scale;
        const bool   ordered = index == 0 || places[index - 1].Along < place.Along;
        if (!(Distance(place.At, At(curve, place.Along)) <= slack)
            || !(through || TurnsAbout(curve, point, place.Along))
            || !ordered)
        {
          return false;
        }
      }
      return true;
    }

    /**
     *  Whether a place is a point to within a few bits.
     *
     *  @param  place The place.
     *  @param  point The point expected.
     *  @return       Whether it is.
     */
    auto Near(const Place &place, Point point) -> bool
    {
      return Distance(place.At, point) < 1e-12;
    }
  } // namespace

  TEST_CASE("tangents: a line, a ray and an xline touch nothing")
  {
    CHECK(Tangents(Line{.Start = {}, .End = {.X = 1, .Y = 0}}, {.X = 0, .Y = 1}).empty());
    CHECK(Tangents(Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}, {.X = 0, .Y = 1}).empty());
    CHECK(Tangents(Xline{.Through = {}, .Way = {.X = 1, .Y = 0}}, {.X = 0, .Y = 1}).empty());
    CHECK(Tangents(Curve{Line{.Start = {}, .End = {.X = 1, .Y = 0}}}, {.X = 0, .Y = 1}).empty());
  }

  TEST_CASE("tangents: an arc in closed form, two from outside, the point on it, none within")
  {
    const Arc                circle = Circle({}, 1);
    const std::vector<Place> two    = Tangents(circle, {.X = 2, .Y = 0});
    REQUIRE(two.size() == 2);
    CHECK(Near(two[0], {.X = 0.5, .Y = std::sqrt(0.75)}));
    CHECK(Near(two[1], {.X = 0.5, .Y = -std::sqrt(0.75)}));
    CHECK(two[0].Along == doctest::Approx(1.0 / 6.0));

    CHECK(Tangents(circle, {.X = 0.5, .Y = 0}).empty());
    CHECK(Tangents(circle, {.X = 0, .Y = 1}).size() == 1);
    CHECK(Tangents(Circle({}, 0), {.X = 2, .Y = 0}).empty());

    const Arc backward = {
        .Center = {.X = 1, .Y = 0.5}, .Radius = -0.8, .Start = Degrees(30), .Sweep = Degrees(-300)};
    for (const Point point :
         {Point{.X = 3, .Y = 3}, Point{.X = 1.2, .Y = 0.4}, Point{.X = -2, .Y = 0}})
    {
      CHECK(Touch(backward, point));
    }
    CHECK(Tangents(backward, {.X = -2, .Y = 0}).size() == 1);
    CHECK(Tangents(backward, {.X = 1, .Y = 3}).size()  == 2);
  }

  TEST_CASE("tangents: an ellipse by its quarters")
  {
    const Ellipse            wide = {.Map = Scaling(3, 1), .Start = {}, .Sweep = Turn};
    const std::vector<Place> two  = Tangents(wide, {.X = 0, .Y = 2});
    REQUIRE(two.size() == 2);
    CHECK(Near(two[0], {.X = 3 * std::sqrt(0.75), .Y = 0.5}));
    CHECK(Near(two[1], {.X = -3 * std::sqrt(0.75), .Y = 0.5}));
    CHECK(Tangents(wide, {}).empty());

    const Ellipse turned = {
        .Map   = Translation({.X = 1, .Y = 0.5}) * Rotation(Degrees(30)) * Scaling(1.5, 0.4),
        .Start = Degrees(20),
        .Sweep = Degrees(300)};
    for (const Point point :
         {Point{.X = 3, .Y = 3}, Point{.X = -2, .Y = 0}, Point{.X = 1, .Y = -3}})
    {
      CHECK(Touch(turned, point));
    }
  }

  TEST_CASE("tangents: a polyline's bulges and the vertices it turns about the point at")
  {
    const Polyline bulged = {.Vertices = {{.Position = {.X = -0.5, .Y = 0.2}, .Bulge = 0.8},
                                          {.Position = {.X = 1.5, .Y = 0.2}, .Bulge = -1.5},
                                          {.Position = {.X = 2.5, .Y = 0.9}, .Bulge = 0},
                                          {.Position = {.X = 0.5, .Y = 1.5}, .Bulge = 0}},
                             .Closed   = true};
    for (const Point point :
         {Point{.X = 3, .Y = 3}, Point{.X = 1, .Y = 0.8}, Point{.X = -2, .Y = 0}})
    {
      CHECK(Touch(bulged, point));
    }
    CHECK(Tangents(bulged, {.X = -2, .Y = 0}).size() >= 2);

    const Polyline open = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                        {.Position = {.X = 1, .Y = 1}, .Bulge = 0},
                                        {.Position = {.X = 2, .Y = 0}, .Bulge = 0}},
                           .Closed   = false};
    CHECK(Tangents(open, {.X = 1, .Y = 3}).empty());
    const std::vector<Place> peak = Tangents(open, {.X = 3, .Y = 1});
    REQUIRE(peak.size() == 1);
    CHECK(peak[0]       == Place{.At = {.X = 1, .Y = 1}, .Along = 1});
  }

  TEST_CASE("tangents: a spline of any weights, its touches as exact arithmetic gives them")
  {
    /**
     *  Before, touches beside a weight of 2^200 or more were missed, and under a weight of 2^-52
     *  each end was a touch from every point.
     */
    for (const Tests::Asked &asked : Tests::Asks(Tests::Touches))
    {
      CAPTURE(asked.Curve.Weights[1]);
      CAPTURE(asked.From);
      const std::vector<Place> places = Tangents(asked.Curve, asked.From);
      REQUIRE(places.size() == asked.Places.size());
      for (std::size_t index = 0; index < places.size(); ++index)
      {
        CHECK(Tests::Lands(places[index].At, asked.Places[index]));
        CHECK(Tests::Tells(places[index].Along, asked.Places[index]));
      }
    }
    CHECK(Tangents(Tests::Stiff("conic", -100), {.X = 1, .Y = -1}).empty());
  }

  TEST_CASE("tangents: a polyline's vertex given twice is still the corner it turns at")
  {
    const Polyline once  = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                         {.Position = {.X = 1, .Y = 0}},
                                         {.Position = {.X = 1, .Y = 1}}}};
    const Polyline twice = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                         {.Position = {.X = 1, .Y = 0}},
                                         {.Position = {.X = 1, .Y = 0}},
                                         {.Position = {.X = 1, .Y = 1}}}};
    const Point    point = {.X = 2, .Y = 2};
    REQUIRE(Tangents(once, point).size()  == 1);
    REQUIRE(Tangents(twice, point).size() == 1);
    CHECK(Tangents(twice, point)[0].At    == Point{.X = 1, .Y = 0});
  }

  TEST_CASE("tangents: splines and Bézier curves by the touches on their pieces")
  {
    const Bezier::Quadratic arch = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 2, .Y = 0}};
    const std::vector<Place> over = Tangents(arch, {.X = 1, .Y = 5});
    CHECK(over.empty());
    CHECK(Tangents(arch, {.X = 1, .Y = -1}).empty());
    const std::vector<Place> low = Tangents(arch, {.X = 1, .Y = 1.5});
    REQUIRE(low.size() == 2);
    CHECK(low[0].Along == doctest::Approx((1 - std::sqrt(0.5)) / 2));

    for (const Point point : {Point{.X = 3, .Y = 3}, Point{.X = 1, .Y = -1}, Point{.X = -2, .Y = 0},
                              Point{.X = 6, .Y = -1}})
    {
      CHECK(Touch(Tests::Uneven(), point, 1e-13));
      CHECK(Touch(Tests::Quarter(), point));
      CHECK(Touch(Tests::Wave, point));
      CHECK(Touch(arch, point));
    }
    CHECK(!Tangents(Tests::Uneven(), {.X = 3, .Y = 3}).empty());
    CHECK(!Tangents(Tests::Wave, {.X = 2, .Y = -1}).empty());
  }

  /**
   *  The places of an arc or an ellipse swept round a million times whose tangent passes through a
   *  point: two a turn, of its first 65,536 turns. And an ellipse thin past what doubles hold, its
   *  axes some 2^12 and 2^116 long, whose quarters give five spots in a turn, a place found twice:
   *  four places a turn at most.
   */
  TEST_CASE("tangents: an arc or ellipse swept past 65,536 turns, and one too thin for doubles")
  {
    constexpr std::size_t Turns = 65'536;
    const Point           point = {.X = 3, .Y = 0.5};
    const Arc arc = {.Center = {}, .Radius = 1, .Start = Degrees(0), .Sweep = Degrees(360e6)};
    const std::vector<Place> places = Tangents(arc, point);
    REQUIRE(places.size() == 2 * Turns);
    CHECK(places.back().Along <= static_cast<double>(Turns) / 1e6);
    const Ellipse ellipse = {.Map   = {.X = {.X = 2, .Y = 0}, .Y = {.X = 0, .Y = 1}, .Origin = {}},
                             .Start = Degrees(0),
                             .Sweep = Degrees(360e6)};
    CHECK(Tangents(ellipse, point).size() == 2 * Turns);

    Ellipse     thin = {.Map = {.X      = {.X = -0x1.07371dc9f28b3p+12, .Y = 0},
                                .Y      = {.X = 0x1.aa166e0c4b24p+116, .Y = 0x1.04713bbc39374p+11},
                                .Origin = {.X = -0x1.dceae57b97c51p-1, .Y = -0x1.f28fb970a8212p-2}},
                        .Start = Degrees(0),
                        .Sweep = Degrees(360)};
    const Point asked = {.X = -0x1.15985b2ff109ap+3, .Y = 0x1.5fe1bb0d22556p+9};
    CHECK(Tangents(thin, asked).size() <= 4);
    thin.Sweep = Degrees(360e4);
    CHECK(Tangents(thin, asked).size() <= 40'001);

    /**
     *  Axes 1e279 and 1e255 long nearly along one line, whose feet are a thousand spots of rounding
     *  noise a turn: the places whose tangent passes through a point are few even there.
     */
    const Ellipse noise = {.Map   = {.X      = {.X = -7.92994e+07, .Y = 7.18053e+279},
                                     .Y      = {.X = -6.82878e+21, .Y = -1.56768e+255},
                                     .Origin = {.X = 0.754822, .Y = 0.724193}},
                           .Start = Degrees(0),
                           .Sweep = Degrees(360)};
    CHECK(Tangents(noise, {.X = -0.344748, .Y = -43.4504}).size() <= 4);
  }
} // namespace Geoxy
