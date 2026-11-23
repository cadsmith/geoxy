/***************************************************************************************************
 *  @file      tests/geoxy/nearest.cpp
 *  @brief     Tests of a curve's places nearest a point: each foot square to the curve, and the
 *             nearest place no farther than any dense sample.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/nearest.hpp"

#include "curves.hpp"
#include "draws.hpp"
#include "polylines.hpp"
#include "stiff.hpp"
#include "waves.hpp"

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
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether the distance from a point turns at a parameter: no less, or no more, a step either
     *  side, the step before the start of a closed curve taken from its end.
     *
     *  @param  curve The curve.
     *  @param  point The point.
     *  @param  along The parameter.
     *  @return       Whether it does.
     */
    auto Turns(const Curve &curve, Point point, double along) -> bool
    {
      const Range  domain = Domain(curve);
      const double back   = along - 1e-7 < domain.Start ? domain.End - 1e-7 : along - 1e-7;
      const double here   = Distance(point, At(curve, along));
      const double before = Distance(point, At(curve, back));
      const double after  = Distance(point, At(curve, along + 1e-7));
      return (before >= here && after >= here) || (before <= here && after <= here);
    }

    /**
     *  How far the point of a place found on a piece of a spline of uneven weights may lie from the
     *  spline's point at the place's parameter, the parameter a double and the spline fast.
     */
    constexpr double Told = 1e-13;

    /**
     *  Whether every foot of a curve is a place on it whose tangent is square to the way to the
     *  point, or where the distance turns at a corner, the feet in order.
     *
     *  @param  curve The curve.
     *  @param  point The point.
     *  @param  slack How far a foot's point may lie from the curve's at its parameter: nothing, but
     *                for an arc or a polyline, whose foot is the exact one rounded, and a spline
     *                whose weights are uneven, whose foot is told as found on a piece of it: not
     *                through the parameter's rounding.
     *  @return       Whether they are.
     */
    auto Square(const Curve &curve, Point point, double slack = 0.0) -> bool
    {
      const std::vector<Place> feet = Perpendiculars(curve, point);
      for (std::size_t index = 0; index < feet.size(); ++index)
      {
        const Place &foot    = feet[index];
        const Vector tangent = Derivative(curve, foot.Along);
        const Vector away    = point - foot.At;
        const double scale   = std::max(1.0, Length(away) * Length(tangent));
        const bool   ordered = index == 0 || feet[index - 1].Along < foot.Along;
        const bool   square  = std::abs(Dot(away, tangent)) <= 1e-9 * scale;
        if (!(Distance(foot.At, At(curve, foot.Along)) <= slack)
            || !(square || Turns(curve, point, foot.Along))
            || !ordered)
        {
          return false;
        }
      }
      return true;
    }

    /**
     *  Whether a curve's nearest place is no farther from a point than any dense sample of it.
     *
     *  @param  curve The curve.
     *  @param  point The point.
     *  @param  slack How far its point may lie from the curve's at its parameter, as Square takes.
     *  @return       Whether it is.
     */
    auto Settles(const Curve &curve, Point point, double slack = 0.0) -> bool
    {
      const Place   nearest = Nearest(curve, point);
      const Range   domain  = Domain(curve);
      const double  least   = Distance(point, nearest.At);
      constexpr int Samples = 20'000;
      for (int index = 0; index <= Samples; ++index)
      {
        const double along = domain.Start + ((domain.End - domain.Start) * index / Samples);
        if (Distance(point, At(curve, along)) < least - 1e-12)
        {
          return false;
        }
      }
      return Distance(nearest.At, At(curve, nearest.Along)) <= slack;
    }

    /**
     *  Whether a place is a point and a parameter to within a few bits.
     *
     *  @param  place The place.
     *  @param  point The point expected.
     *  @param  along The parameter expected.
     *  @return       Whether it is.
     */
    auto Near(const Place &place, Point point, double along) -> bool
    {
      return Distance(place.At, point) < 1e-12 && std::abs(place.Along - along) < 1e-12;
    }

    /**
     *  The place of a polyline nearest a point as every foot and vertex weighed gives it: its feet,
     *  chained and kept once, and its vertices, put in order by parameter, the first strictly
     *  nearest of them.
     *
     *  @param  polyline The polyline.
     *  @param  point    The point.
     *  @return          The place.
     */
    auto Weighed(const Polyline &polyline, Point point) -> Place
    {
      std::vector<Place> places = Perpendiculars(polyline, point);
      for (std::size_t index = 0; index <= Count(polyline); ++index)
      {
        const auto along = static_cast<double>(index);
        places.push_back({.At = At(polyline, along), .Along = along});
      }
      std::ranges::stable_sort(places, {}, &Place::Along);
      Place  best  = places.front();
      double least = Distance(point, best.At);
      for (const Place &place : places)
      {
        if (Distance(point, place.At) < least)
        {
          best  = place;
          least = Distance(point, place.At);
        }
      }
      return best;
    }

    /**
     *  Points to ask a polyline for its nearest place: anywhere about it, on vertices and a step of
     *  a double off them, on segments, and at a bulge's center.
     *
     *  @param  draws    The draws.
     *  @param  polyline The polyline.
     *  @return          The points.
     */
    auto Asked(Draws &draws, const Polyline &polyline) -> std::vector<Point>
    {
      std::vector<Point> points;
      points.reserve(28 + Count(polyline));
      for (int index = 0; index < 4; ++index)
      {
        points.push_back({.X = draws.Near(0.5, 1.5), .Y = draws.Near(0.5, 1.5)});
      }
      for (int index = 0; index < 6 && !polyline.Vertices.empty(); ++index)
      {
        const Point corner = polyline.Vertices[draws.Next() % polyline.Vertices.size()].Position;
        points.push_back(corner);
        points.push_back({.X = std::nextafter(corner.X, 9.0), .Y = corner.Y});
      }
      const auto count = static_cast<double>(Count(polyline));
      for (int index = 0; index < 4 && count > 0; ++index)
      {
        const double along = draws.Near(count / 2, count / 2);
        points.push_back(At(polyline, along));
        points.push_back(At(polyline, std::floor(along) + 0.5));
      }
      for (std::size_t index = 0; index < Count(polyline); ++index)
      {
        const Segment segment = Nth(polyline, index);
        if (const auto *arc = std::get_if<Arc>(&segment))
        {
          points.push_back(arc->Center);
        }
      }
      return points;
    }
  } // namespace

  TEST_CASE("nearest: a line, a ray and an xline by the foot of the perpendicular")
  {
    const Line line = {.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}};
    CHECK(Perpendiculars(line, {.X = 1, .Y = 2})
          == std::vector<Place>{{.At = {.X = 1, .Y = 0}, .Along = 0.25}});
    CHECK(Perpendiculars(line, {.X = 5, .Y = 1}).empty());
    CHECK(Nearest(line, {.X = 5, .Y = 1})  == Place{.At = {.X = 4, .Y = 0}, .Along = 1});
    CHECK(Nearest(line, {.X = 1, .Y = -2}) == Place{.At = {.X = 1, .Y = 0}, .Along = 0.25});

    const Line point = {.Start = {.X = 1, .Y = 1}, .End = {.X = 1, .Y = 1}};
    CHECK(Perpendiculars(point, {.X = 3, .Y = 3}).empty());
    CHECK(Nearest(point, {.X = 3, .Y = 3}) == Place{.At = {.X = 1, .Y = 1}, .Along = 0});

    const Ray ray = {.Start = {.X = 0, .Y = 0}, .Way = {.X = 2, .Y = 0}};
    CHECK(Perpendiculars(ray, {.X = 4, .Y = 3})
          == std::vector<Place>{{.At = {.X = 4, .Y = 0}, .Along = 2}});
    CHECK(Perpendiculars(ray, {.X = -3, .Y = 1}).empty());
    CHECK(Nearest(ray, {.X = -3, .Y = 1}) == Place{.At = {.X = 0, .Y = 0}, .Along = 0});
    CHECK(Nearest(ray, {.X = 4, .Y = 3})  == Place{.At = {.X = 4, .Y = 0}, .Along = 2});

    const Xline xline = {.Through = {.X = 0, .Y = 1}, .Way = {.X = 1, .Y = 0}};
    CHECK(Nearest(xline, {.X = -7, .Y = 5}) == Place{.At = {.X = -7, .Y = 1}, .Along = -7});
    const Xline still = {.Through = {.X = 0, .Y = 1}, .Way = {}};
    CHECK(Nearest(still, {.X = -7, .Y = 5}) == Place{.At = {.X = 0, .Y = 1}, .Along = 0});
  }

  TEST_CASE("nearest: an arc in closed form, every turn of it")
  {
    const Arc                circle = Circle({}, 2);
    const std::vector<Place> feet   = Perpendiculars(circle, {.X = 3, .Y = 4});
    const double             toward = std::atan2(4.0, 3.0) / (2 * std::numbers::pi);
    REQUIRE(feet.size() == 2);
    CHECK(Near(feet[0], {.X = 1.2, .Y = 1.6}, toward));
    CHECK(Near(feet[1], {.X = -1.2, .Y = -1.6}, toward + 0.5));
    CHECK(Near(Nearest(circle, {.X = 3, .Y = 4}), {.X = 1.2, .Y = 1.6}, toward));

    CHECK(Perpendiculars(circle, {}).empty());
    CHECK(Nearest(circle, {}) == Place{.At = {.X = 2, .Y = 0}, .Along = 0});

    const Arc twice = {.Center = {}, .Radius = 1, .Start = {}, .Sweep = Turn * 2};
    CHECK(Perpendiculars(twice, {.X = 3, .Y = 4}).size() == 4);
    const Arc none = {.Center = {}, .Radius = 1, .Start = Degrees(10), .Sweep = {}};
    CHECK(Perpendiculars(none, {.X = 3, .Y = 4}).empty());

    /**
     *  Swept round many times it passes each foot on every turn, all as near: its nearest place is
     *  on the first pass, and is found without weighing the rest, which made a meeting with such an
     *  arc take the square of its turns.
     */
    const double start = std::numbers::pi / 18;
    for (const double sweep : {1e5, -1e5, 7.0, -7.0})
    {
      CAPTURE(sweep);
      const Arc coil = {
          .Center = {}, .Radius = 2, .Start = Radians(start), .Sweep = Radians(sweep)};
      const double ahead = std::atan2(4.0, 3.0) - start;
      const double first = (sweep > 0 ? ahead : (2 * std::numbers::pi) - ahead) / std::abs(sweep);
      const Place  place = Nearest(coil, {.X = 3, .Y = 4});
      CHECK(place.Along == doctest::Approx(first).epsilon(1e-12));
      CHECK(Distance(place.At, {.X = 1.2, .Y = 1.6}) < 1e-9);
      CHECK(Nearest(coil, {}) == Place{.At = At(coil, 0), .Along = 0});
    }

    const Arc bent = {
        .Center = {.X = 1, .Y = 0.5}, .Radius = 0.8, .Start = Degrees(30), .Sweep = Degrees(-300)};
    for (const Point point :
         {Point{.X = 3, .Y = 3}, Point{.X = 1.2, .Y = 0.4}, Point{.X = -2, .Y = 0}})
    {
      CHECK(Square(bent, point, Told));
      CHECK(Settles(bent, point, Told));
    }
  }

  /**
   *  An arc's foot was At's point at the foot's parameter, a few steps of a double off the exact
   *  foot: of a circle of radius 13 about the origin, the foot toward (10, 24) came out at
   *  (5.0000000000000036, 11.999999999999998). It is now the circle's point toward the point,
   *  worked in double-double and rounded once: the double nearest the exact foot.
   */
  TEST_CASE("nearest: an arc's foot is the double nearest the exact one")
  {
    struct Triple
    {
      double Across = 0.0;
      double Up     = 0.0;
      double Long   = 0.0;
    };

    for (const Triple triple :
         {Triple{.Across = 3, .Up = 4, .Long = 5}, Triple{.Across = 5, .Up = 12, .Long = 13},
          Triple{.Across = 7, .Up = 24, .Long = 25}, Triple{.Across = 33, .Up = 56, .Long = 65}})
    {
      for (const Point center : {Point{.X = 0, .Y = 0}, Point{.X = -7, .Y = -3.5}})
      {
        CAPTURE(triple.Long);
        CAPTURE(center.X);
        const Point foot  = {.X = center.X + triple.Across, .Y = center.Y + triple.Up};
        const Point point = {.X = center.X + (2 * triple.Across), .Y = center.Y + (2 * triple.Up)};
        for (const double radius : {triple.Long, -triple.Long})
        {
          const Arc arc = {.Center = center,
                           .Radius = radius,
                           .Start  = Degrees(radius > 0 ? 1 : 181),
                           .Sweep  = Degrees(170)};
          CHECK(Nearest(arc, point).At == foot);
          CHECK(Perpendiculars(arc, point).front().At == foot);
          const Arc whole = {.Center = center, .Radius = radius, .Start = {}, .Sweep = Turn * 3};
          CHECK(Nearest(whole, point).At == foot);
        }
      }
    }

    /**
     *  A point past 2^450 from the center, or within 2^-450 of it, whose square would leave the
     *  doubles: the foot in plain doubles.
     */
    const Arc unit = {.Center = {}, .Radius = 1, .Start = Degrees(-90), .Sweep = Degrees(180)};
    for (const Point point : {Point{.X = 0x1p460, .Y = 0}, Point{.X = 0x1p-460, .Y = 0}})
    {
      CHECK(std::abs(Length(Perpendiculars(unit, point).front().At - Point{}) - 1.0) <= 1e-15);
    }
  }

  /**
   *  An arc swept round once at most is answered by its angles where they tell its nearest place
   *  clearly, with one point worked out: the same place to the bit as its start, its feet in order
   *  and its end weighed, of arcs held in degrees and radians, of either sign of radius, of every
   *  sweep to a whole turn, and of points on the arc, at its center, in line with its ends and
   *  anywhere about it.
   */
  TEST_CASE("nearest: an arc by its angles, as its feet and ends weighed give it, to the bit")
  {
    CHECK(Turn.InRadians() == 2 * std::numbers::pi);
    const auto weighed = [](const Arc &arc, Point point) -> Place
    {
      std::vector<Place>       places = {{.At = At(arc, 0.0), .Along = 0.0}};
      const std::vector<Place> feet   = Perpendiculars(arc, point);
      places.insert(places.end(), feet.begin(), feet.end());
      places.push_back({.At = At(arc, 1.0), .Along = 1.0});
      Place  best  = places.front();
      double least = Distance(point, best.At);
      for (const Place &place : places)
      {
        if (Distance(point, place.At) < least)
        {
          best  = place;
          least = Distance(point, place.At);
        }
      }
      return best;
    };
    Draws       draws{0x5EED'0A2C};
    std::size_t asked = 0;
    std::size_t wrong = 0;
    for (int round = 0; round < Rounds(4000); ++round)
    {
      const bool   inDegrees = round % 3 != 0;
      const double scale     = round % 11 == 10 ? 1e8 : 1.0;
      const Point  center    = {.X = draws.Near(0, 5) * scale, .Y = draws.Near(0, 5) * scale};
      const double radius    = draws.Near(2, 1.9) * (round % 4 == 3 ? -1.0 : 1.0);
      double       start     = inDegrees ? draws.Near(0, 720) : draws.Near(0, 12);
      double       sweep     = inDegrees ? draws.Near(0, 360) : draws.Near(0, 6.28);
      if (round % 7 == 6)
      {
        sweep = inDegrees ? 360.0 : 2 * std::numbers::pi;
      }
      if (round % 13 == 12)
      {
        sweep = inDegrees ? 1e-9 : -1e-11;
      }
      if (round % 17 == 16)
      {
        start = inDegrees ? 1e9 : 1e7;
      }
      const Arc   arc    = {.Center = center,
                            .Radius = radius,
                            .Start  = inDegrees ? Degrees(start) : Radians(start),
                            .Sweep  = inDegrees ? Degrees(sweep) : Radians(sweep)};
      const Point end    = At(arc, 1.0);
      const Point ahead  = {.X = center.X + (3 * (end.X - center.X)),
                            .Y = center.Y + (3 * (end.Y - center.Y))};
      const Point inside = {.X = center.X + ((end.X - center.X) / 4),
                            .Y = center.Y + ((end.Y - center.Y) / 4)};
      for (const Point point :
           {Point{.X = center.X + draws.Near(0, 8), .Y = center.Y + draws.Near(0, 8)},
            At(arc, draws.Near(0.5, 0.5)), center, ahead, inside, end})
      {
        ++asked;
        wrong += Nearest(arc, point) == weighed(arc, point) ? 0U : 1U;
      }
    }
    CHECK(asked > 0);
    CHECK(wrong == 0);
    const Arc nothing = {.Center = {}, .Radius = 0, .Start = {}, .Sweep = Degrees(90)};
    CHECK(Nearest(nothing, {.X = 1, .Y = 1}) == weighed(nothing, {.X = 1, .Y = 1}));
  }

  TEST_CASE("nearest: an ellipse by its quarters, up to four feet a turn")
  {
    const Ellipse wide = {.Map = Scaling(3, 1), .Start = {}, .Sweep = Turn};
    CHECK(Perpendiculars(wide, {.X = 0.1, .Y = 0.1}).size() == 4);
    CHECK(Perpendiculars(wide, {.X = 0, .Y = 5}).size()     == 4);
    CHECK(Perpendiculars(wide, {.X = 0, .Y = 9}).size()     == 2);
    CHECK(Near(Nearest(wide, {.X = 0, .Y = 5}), {.X = 0, .Y = 1}, 0.25));
    CHECK(Perpendiculars(Ellipse{.Map = Scaling(2), .Start = {}, .Sweep = Turn}, {}).empty());
    CHECK(Perpendiculars(Ellipse{.Map = Scaling(2), .Start = {}, .Sweep = Turn}, {.X = 3, .Y = 0})
              .size()
          == 3);
    CHECK(Perpendiculars(Ellipse{.Map = Scaling(0), .Start = {}, .Sweep = Turn}, {.X = 1, .Y = 1})
              .empty());

    const Ellipse turned = {
        .Map   = Translation({.X = 1, .Y = 0.5}) * Rotation(Degrees(30)) * Scaling(1.5, 0.4),
        .Start = Degrees(20),
        .Sweep = Degrees(300)};
    const Ellipse folded = {.Map = Scaling(-1.5, 0.8), .Start = Degrees(-40), .Sweep = Turn * 2};
    for (const Point point :
         {Point{.X = 3, .Y = 3}, Point{.X = 1.2, .Y = 0.4}, Point{.X = -2, .Y = 0}})
    {
      CHECK(Square(turned, point));
      CHECK(Settles(turned, point));
      CHECK(Square(folded, point));
      CHECK(Settles(folded, point));
    }
  }

  /**
   *  An ellipse swept round many times had every pass of every foot listed and weighed for its
   *  nearest place, the cost of a meeting with one the square of its turns: its first turn holds
   *  the first of the nearest, and is all that is looked at.
   */
  /**
   *  An ellipse's feet by the roots of the secular equation: four within its evolute and two
   *  outside it, every one square to the way to the point, in order; and where the roots are not
   *  sure, an ellipse thin past 2^-20, a point on an axis or by the evolute, or a map that is not
   *  finite, as its quarters find them, as before.
   */
  TEST_CASE("nearest: an ellipse's feet by the secular equation, four within the evolute")
  {
    Draws       draws{0x5EED'0E11};
    std::size_t wrong = 0;
    for (int round = 0; round < Rounds(3000); ++round)
    {
      const double  major  = draws.Near(2, 1.5);
      const double  minor  = major * draws.Near(0.5, 0.45);
      const Matrix  moved  = Translation({.X = draws.Near(0, 5), .Y = draws.Near(0, 5)});
      const Matrix  map    = moved * Rotation(Degrees(draws.Near(0, 180))) * Scaling(major, minor);
      const Ellipse whole  = {.Map = map, .Start = Degrees(draws.Near(0, 180)), .Sweep = Turn};
      const double  along  = draws.Near(0, 1.5 * major);
      const double  across = draws.Near(0, 1.5 * major);
      const Point   point  = map * Point{.X = along / major, .Y = across / minor};
      const double  spread = (major * major) - (minor * minor);
      const double  reach =
          std::cbrt(major * major * along * along) + std::cbrt(minor * minor * across * across);
      const std::size_t expected = reach * reach * reach < spread * spread ? 4 : 2;
      wrong += Perpendiculars(whole, point).size() == expected && Square(whole, point) ? 0U : 1U;
    }
    CHECK(wrong == 0);

    const Ellipse wide = {.Map = Scaling(3, 1), .Start = Degrees(10), .Sweep = Turn};
    const double  cusp = std::pow(4 - std::cbrt(9.0), 1.5);
    CHECK(Perpendiculars(wide, {.X = 1, .Y = 0}).size()    == 4);
    CHECK(Perpendiculars(wide, {.X = 0, .Y = 0.5}).size()  == 4);
    CHECK(Perpendiculars(wide, {.X = 1, .Y = cusp}).size() == 4);
    for (const Point point :
         {Point{.X = 1, .Y = 0}, Point{.X = 0, .Y = 0.5}, Point{.X = 1, .Y = cusp}})
    {
      CHECK(Square(wide, point));
    }
    const Ellipse thin = {.Map = Scaling(1, 0x1p-25), .Start = Degrees(10), .Sweep = Turn};
    CHECK(Perpendiculars(thin, {.X = 0.3, .Y = 0.2}).size() == 4);
    CHECK(Square(thin, {.X = 0.3, .Y = 0.2}));
    const double  hollow = std::numeric_limits<double>::quiet_NaN();
    const Ellipse lost = {.Map = {.X = {.X = hollow, .Y = 0}, .Y = {.X = 0, .Y = 1}, .Origin = {}},
                          .Start = {},
                          .Sweep = Turn};
    CHECK(std::ranges::all_of(Perpendiculars(lost, {.X = 0.3, .Y = 0.2}),
                              [](const Place &foot) -> bool
                              { return foot.Along >= 0 && foot.Along <= 1; }));
  }

  TEST_CASE("nearest: an ellipse swept round many times, by its first turn")
  {
    const Matrix map = {.X = {.X = 3, .Y = 0}, .Y = {.X = 0, .Y = 1}, .Origin = {.X = 1, .Y = 2}};
    for (const double sweep : {1e5, -1e5, 9.0, -9.0})
    {
      CAPTURE(sweep);
      const Ellipse coil = {.Map = map, .Start = Degrees(20), .Sweep = Radians(sweep)};
      const Ellipse once = {.Map = map, .Start = Degrees(20), .Sweep = sweep > 0 ? Turn : -Turn};
      for (const Point point :
           {Point{.X = 5, .Y = 4}, Point{.X = 1, .Y = 2.5}, Point{.X = -3, .Y = 0}})
      {
        const Place owed  = Nearest(once, point);
        const Place place = Nearest(coil, point);
        CHECK(
            place.Along
            == doctest::Approx(owed.Along * 2 * std::numbers::pi / std::abs(sweep)).epsilon(1e-12));
        CHECK(Distance(place.At, owed.At) < 1e-9);
      }
    }
  }

  TEST_CASE("nearest: a polyline, segment by segment")
  {
    const Polyline bulged = {.Vertices = {{.Position = {.X = -0.5, .Y = 0.2}, .Bulge = 0.8},
                                          {.Position = {.X = 1.5, .Y = 0.2}, .Bulge = -1.5},
                                          {.Position = {.X = 2.5, .Y = 0.9}, .Bulge = 0},
                                          {.Position = {.X = 0.5, .Y = 1.5}, .Bulge = 0}},
                             .Closed   = true};
    for (const Point point :
         {Point{.X = 3, .Y = 3}, Point{.X = 1, .Y = 0.8}, Point{.X = -2, .Y = 0}})
    {
      CHECK(Square(bulged, point, Told));
      CHECK(Settles(bulged, point, Told));
    }
    CHECK(!Perpendiculars(bulged, {.X = 1, .Y = 0.8}).empty());

    /**
     *  A bulge past a half turn has its far foot on it, the circle's point away from the point,
     *  where its near one is off it.
     */
    const Polyline           wide = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 3},
                                                  {.Position = {.X = 1, .Y = 0}, .Bulge = 0}},
                                     .Closed   = false};
    const std::vector<Place> both = Perpendiculars(wide, {.X = 0.5, .Y = -0.6});
    REQUIRE(both.size() == 1);
    CHECK(Distance(both.front().At, {.X = 0.5, .Y = -1.5}) <= 1e-15);
    CHECK(Square(wide, {.X = 0.5, .Y = -0.6}, Told));

    const Polyline lone = {.Vertices = {{.Position = {.X = 1, .Y = 1}, .Bulge = 0}},
                           .Closed   = false};
    CHECK(Perpendiculars(lone, {}).empty());
    CHECK(Nearest(lone, {}) == Place{.At = {.X = 1, .Y = 1}, .Along = 0});
    CHECK(Nearest(Polyline{}, {.X = 3, .Y = 3}) == Place{});
  }

  /**
   *  A polyline's foot was At's point at the foot's parameter: on this bulge of -1.14, 16 steps of
   *  a double off the exact foot in x. Each foot is now the exact one worked in double-double from
   *  the segment's start, a bulge's circle with no center rounded, and rounded once.
   */
  TEST_CASE("nearest: a polyline's foot is the double nearest the exact one")
  {
    const Polyline polyline = {
        .Vertices =
            {{.Position = {.X = 4.4326217109732955, .Y = -1.0083910230751867}, .Bulge = 0.0},
             {.Position = {.X = -1.0296071670230507, .Y = 3.826161358428565},
              .Bulge    = -1.1356368862624984},
             {.Position = {.X = -2.2304894499093804, .Y = -0.5948343584680036}, .Bulge = 0.0},
             {.Position = {.X = 0.5185594838881735, .Y = -2.3189161217751066}, .Bulge = 0.0}},
        .Closed = false};
    const Place place = Nearest(polyline, {.X = -1.103748684122844, .Y = 0.2273960669458619});
    CHECK(place.At == Point{.X = -0x1.d9b28b5e799dcp-1, .Y = -0x1.7645eef1172e0p-1});
    CHECK(Perpendiculars(polyline, {.X = -1.103748684122844, .Y = 0.2273960669458619}).size() >= 1);

    /**
     *  A straight segment's foot is its exact projection: (3, 4) on the line from the origin to (6,
     *  8), asked from (7, 1), where the share 0.5 times the way is exact too.
     */
    const Polyline straight = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                            {.Position = {.X = 6, .Y = 8}, .Bulge = 0}},
                               .Closed   = false};
    CHECK(Nearest(straight, {.X = 7, .Y = 1}).At == Point{.X = 3, .Y = 4});
    const Polyline slanted = {.Vertices = {{.Position = {.X = 0.1, .Y = 0.7}, .Bulge = 0},
                                           {.Position = {.X = 3.3, .Y = -1.9}, .Bulge = 0}},
                              .Closed   = false};
    const Point    foot    = Nearest(slanted, {.X = 2.0, .Y = 1.5}).At;
    const Vector   along   = Point{.X = 3.3, .Y = -1.9} - Point{.X = 0.1, .Y = 0.7};
    CHECK(std::abs(Cross(along, foot - Point{.X = 0.1, .Y = 0.7})) <= 1e-15);

    /**
     *  A bulge far from the origin, past 2^450, takes At's point there.
     */
    const Polyline far = {.Vertices = {{.Position = {.X = 0x1p460, .Y = 0}, .Bulge = 0.5},
                                       {.Position = {.X = 0x1p461, .Y = 0}, .Bulge = 0}},
                          .Closed   = false};
    const Place    out = Nearest(far, {.X = 0x1.8p460, .Y = -0x1p462});
    CHECK(out.At.Y < 0.0);
    const Polyline line = {.Vertices = {{.Position = {.X = 0x1p460, .Y = 0}, .Bulge = 0},
                                        {.Position = {.X = 0x1p461, .Y = 0}, .Bulge = 0}},
                           .Closed   = false};
    CHECK(Nearest(line, {.X = 0x1.8p460, .Y = 0x1p460}).At == Point{.X = 0x1.8p460, .Y = 0});

    /**
     *  A point a hair from a bulge's center, its center at the origin: At's point there.
     */
    const Polyline round = {.Vertices = {{.Position = {.X = -1, .Y = -0.75}, .Bulge = 0.5},
                                         {.Position = {.X = 1, .Y = -0.75}, .Bulge = 0}},
                            .Closed   = false};
    CHECK(std::abs(Distance(Nearest(round, {.X = 0, .Y = -0x1p-460}).At, Point{}) - 1.25) <= 1e-15);
  }

  TEST_CASE("nearest: a polyline's nearest place as every foot and vertex weighed gives it")
  {
    /**
     *  Only the segments whose boxes come as near as the nearest place so far are solved: the
     *  answer is the one every segment solved gives, to the bit, on polylines of every kind.
     */
    Draws       draws{0x5EED'0015};
    std::size_t asked = 0;
    std::size_t wrong = 0;
    for (int round = 0; round < Rounds(2000); ++round)
    {
      const std::size_t count = 1 + (draws.Next() % 12);
      const Polyline    polyline =
          round % 5 == 0 ? Zigzagged(draws, 2 + (draws.Next() % 200)) : Strung(draws, count);
      for (const Point point : Asked(draws, polyline))
      {
        ++asked;
        wrong += Nearest(polyline, point) == Weighed(polyline, point) ? 0U : 1U;
      }
    }
    CHECK(asked > 0);
    CHECK(wrong == 0);

    /**
     *  A vertex past every double, or two whose sum is, leaves every segment solved, as before: the
     *  same place to the bit, numbers that are not numbers alike.
     */
    const auto same = [](const Place &lhs, const Place &rhs) -> bool
    {
      const auto bits = [](const Place &place) -> std::array<std::uint64_t, 3>
      {
        return {std::bit_cast<std::uint64_t>(place.At.X), std::bit_cast<std::uint64_t>(place.At.Y),
                std::bit_cast<std::uint64_t>(place.Along)};
      };
      return bits(lhs) == bits(rhs);
    };
    const double huge = std::numeric_limits<double>::max();
    for (const Point far :
         {Point{.X = std::numeric_limits<double>::infinity(), .Y = 1}, Point{.X = huge, .Y = huge}})
    {
      const Polyline beyond = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0.5},
                                            {.Position = {.X = 1, .Y = 0}},
                                            {.Position = far}}};
      for (const Point point : {Point{.X = 0.5, .Y = 1}, Point{.X = 2, .Y = 2}})
      {
        CHECK(same(Nearest(beyond, point), Weighed(beyond, point)));
      }
    }
  }

  TEST_CASE("nearest: a spline's nearest place as every span solved gives it, rational or not")
  {
    /**
     *  Only the spans whose boxes come as near as the nearest place so far are solved: the answer
     *  is the one every span solved gives, every foot and every knot weighed, to the bit, on
     *  splines of every degree, clamped or not, of knots repeated, scaled and moved.
     */
    const auto weighed = [](const Spline &spline, Point point) -> Place
    {
      std::vector<Place> places = Perpendiculars(spline, point);
      const auto         degree = static_cast<std::size_t>(spline.Degree);
      for (std::size_t span = degree; span < spline.Poles.size(); ++span)
      {
        if (spline.Knots[span] < spline.Knots[span + 1])
        {
          for (const double along : {spline.Knots[span], spline.Knots[span + 1]})
          {
            places.push_back({.At = At(spline, along), .Along = along});
          }
        }
      }
      std::ranges::stable_sort(places, {}, &Place::Along);
      Place  best  = places.front();
      double least = Distance(point, best.At);
      for (const Place &place : places)
      {
        if (Distance(point, place.At) < least)
        {
          best  = place;
          least = Distance(point, place.At);
        }
      }
      return best;
    };
    const auto same = [](const Place &lhs, const Place &rhs) -> bool
    {
      const auto bits = [](const Place &place) -> std::array<std::uint64_t, 3>
      {
        return {std::bit_cast<std::uint64_t>(place.At.X), std::bit_cast<std::uint64_t>(place.At.Y),
                std::bit_cast<std::uint64_t>(place.Along)};
      };
      return bits(lhs) == bits(rhs);
    };
    Draws       draws{0x5EED'5B1E};
    std::size_t asked = 0;
    std::size_t wrong = 0;
    for (int round = 0; round < Rounds(4000); ++round)
    {
      Spline       spline = Tests::Splined(draws, 1, 4);
      const double scale =
          draws.Next() % 5 == 0 ? std::ldexp(1.37, static_cast<int>(draws.Next() % 80) - 40) : 1.0;
      const double by = draws.Next() % 5 == 0 ? 1e6 : 0.0;
      for (Point &pole : spline.Poles)
      {
        pole = {.X = (pole.X * scale) + by, .Y = (pole.Y * scale) - by};
      }
      REQUIRE(Valid(spline));
      const Range        domain = Domain(spline);
      std::vector<Point> points;
      for (int index = 0; index < 6; ++index)
      {
        points.push_back(
            {.X = (draws.Near(0.5, 1.5) * scale) + by, .Y = (draws.Near(0.0, 1.5) * scale) - by});
        points.push_back(At(
            spline, draws.Near((domain.Start + domain.End) / 2, (domain.End - domain.Start) / 2)));
      }
      for (const double knot : spline.Knots)
      {
        const Point onto = At(spline, std::clamp(knot, domain.Start, domain.End));
        points.push_back(onto);
        points.push_back({.X = std::nextafter(onto.X, 9e9), .Y = onto.Y + (0.25 * scale)});
      }
      points.push_back(spline.Poles[draws.Next() % spline.Poles.size()]);
      /**
       *  And weighted, its weights within a few tenths of one, so that its spans are even and are
       *  weighed by branch and bound, or now and then far apart, its feet found on a paving.
       */
      Spline       rational = spline;
      const double spread   = round % 4 == 3 ? 3.0 : 0.3;
      for (std::size_t index = 0; index < spline.Poles.size(); ++index)
      {
        rational.Weights.push_back(1.0 + (spread * draws.Near(0.0, 1.0) / (1.0 + spread)));
      }
      const bool valid = Valid(rational);
      for (const Point point : points)
      {
        asked += 2;
        wrong += same(Nearest(spline, point), weighed(spline, point)) ? 0U : 1U;
        wrong += !valid || same(Nearest(rational, point), weighed(rational, point)) ? 0U : 1U;
      }
    }
    CHECK(asked > 0);
    CHECK(wrong == 0);

    /**
     *  A point or a pole past 2^500, where a distance may run past the doubles, leaves every span
     *  solved, as before.
     */
    const Spline wave = {.Degree  = 3,
                         .Knots   = {0, 0, 0, 0, 1, 2, 2, 2, 2},
                         .Poles   = {{.X = 0, .Y = 0},
                                     {.X = 1, .Y =  2},
                                     {.X = 2, .Y = -1},
                                     {.X = 3, .Y =  1},
                                     {.X = 4, .Y = 0}},
                         .Weights = {}};
    Spline       far  = wave;
    far.Poles[2].Y    = 1e160;
    Spline heavy      = wave;
    heavy.Weights     = {1, 0x1p41, 1, 1, 1};
    Spline lifted     = far;
    lifted.Weights    = {1, 1.25, 1, 1.5, 1};
    for (const Point point : {Point{.X = 1e200, .Y = 0}, Point{.X = 2, .Y = 0.5}})
    {
      CHECK(same(Nearest(wave, point), weighed(wave, point)));
      CHECK(same(Nearest(far, point), weighed(far, point)));
      CHECK(same(Nearest(heavy, point), weighed(heavy, point)));
      CHECK(same(Nearest(lifted, point), weighed(lifted, point)));
    }
  }

  TEST_CASE("nearest: splines and Bézier curves by the feet on their pieces")
  {
    const Spline             quarter = Tests::Quarter();
    const std::vector<Place> middle  = Perpendiculars(quarter, {.X = 2, .Y = 2});
    REQUIRE(middle.size() == 1);
    CHECK(Near(middle[0], {.X = std::sqrt(0.5), .Y = std::sqrt(0.5)}, 0.5));

    const Bezier::Quadratic arch = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 2, .Y = 0}};
    for (const Point point : {Point{.X = 3, .Y = 3}, Point{.X = 1, .Y = 0.8},
                              Point{.X = -2, .Y = 0}, Point{.X = 6, .Y = -1}})
    {
      CHECK(Square(Tests::Uneven(), point, Told));
      CHECK(Settles(Tests::Uneven(), point, Told));
      CHECK(Square(quarter, point));
      CHECK(Settles(quarter, point));
      CHECK(Square(Tests::Wave, point));
      CHECK(Settles(Tests::Wave, point));
      CHECK(Square(arch, point));
      CHECK(Settles(arch, point));
    }
    CHECK(Near(Nearest(arch, {.X = 1, .Y = 5}), {.X = 1, .Y = 1}, 0.5));
  }

  TEST_CASE("nearest: a spline of any weights, its feet and its nearest place as exact arithmetic")
  {
    /**
     *  Before, a foot beside a weight of 256 lay a hundred steps of a double off, one beside 2^52
     *  on the wrong leg of the curve, and under a weight of 2^-52 each end was a foot of every
     *  point.
     */
    for (const Tests::Asked &asked : Tests::Asks(Tests::Perpendiculars))
    {
      CAPTURE(asked.Curve.Weights[1]);
      CAPTURE(asked.From);
      const std::vector<Place> feet = Perpendiculars(asked.Curve, asked.From);
      REQUIRE(feet.size() == asked.Places.size());
      for (std::size_t index = 0; index < feet.size(); ++index)
      {
        CHECK(Tests::Lands(feet[index].At, asked.Places[index]));
        CHECK(Tests::Tells(feet[index].Along, asked.Places[index]));
      }
    }
    for (const Tests::Asked &asked : Tests::Asks(Tests::Nearests))
    {
      CAPTURE(asked.Curve.Weights[1]);
      CAPTURE(asked.From);
      const Place nearest = Nearest(asked.Curve, asked.From);
      CHECK(Tests::Lands(nearest.At, asked.Places.front()));
      CHECK(Tests::Tells(nearest.Along, asked.Places.front()));
    }
  }

  TEST_CASE("nearest: a light weight's end is no foot, and a weight past doubles is answered")
  {
    /**
     *  A quadratic whose middle weight is 2^-100 runs along its chord: one foot from above it,
     *  where both its ends used to be feet too.
     */
    const Spline             light = Tests::Stiff("conic", -100);
    const std::vector<Place> feet  = Perpendiculars(light, {.X = 0.7, .Y = 0.9});
    REQUIRE(feet.size() == 1);
    CHECK(feet[0].At.X  == doctest::Approx(0.7));
    CHECK(std::abs(feet[0].At.Y) < 1e-25);

    /**
     *  Weights 2^1400 apart in one span are more than doubles hold beside each other: the spline is
     *  answered as its own parameter tells it, with a place of the curve.
     */
    const Spline beyond = {
        .Degree  = 3,
        .Knots   = {0, 0, 0, 0, 1, 1, 1, 1},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 0.5, .Y = 1.5}, {.X = 1.5, .Y = 1.5}, {.X = 2, .Y = 0}},
        .Weights = {0x1p700, 1, 1, 0x1p-700}};
    const Place nearest = Nearest(beyond, {.X = 1, .Y = 3});
    CHECK(std::isfinite(nearest.At.X));
    CHECK(std::isfinite(nearest.At.Y));
    CHECK(nearest.Along >= 0.0);
    CHECK(nearest.Along <= 1.0);
  }

  TEST_CASE("nearest: a polyline's vertex given twice is still the corner the distance turns at")
  {
    /**
     *  Before, the segment of no length between the two hid the turn from both its neighbours.
     */
    const Polyline once  = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                         {.Position = {.X = 1, .Y = 0}},
                                         {.Position = {.X = 1, .Y = 1}}}};
    const Polyline twice = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                         {.Position = {.X = 1, .Y = 0}},
                                         {.Position = {.X = 1, .Y = 0}},
                                         {.Position = {.X = 1, .Y = 1}}}};
    const Point    point = {.X = 2, .Y = -1};
    REQUIRE(Perpendiculars(once, point).size()  == 1);
    REQUIRE(Perpendiculars(twice, point).size() == 1);
    CHECK(Perpendiculars(twice, point)[0].At    == Point{.X = 1, .Y = 0});
    CHECK(Perpendiculars(twice, point)[0].Along == 2);
    CHECK(Perpendiculars(twice, {.X = 0.5, .Y = 3}).size()
          == Perpendiculars(once, {.X = 0.5, .Y = 3}).size());
  }

  TEST_CASE("nearest: any curve as its own kind")
  {
    const Curve line = Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}};
    CHECK(Nearest(line, {.X = 1, .Y = 2}) == Place{.At = {.X = 1, .Y = 0}, .Along = 0.25});
    CHECK(Perpendiculars(line, {.X = 1, .Y = 2}).size() == 1);
    CHECK(Perpendiculars(Curve{Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}}, {.X = 1, .Y = 1}).size()
          == 1);
    CHECK(Perpendiculars(Curve{Xline{.Through = {}, .Way = {.X = 1, .Y = 0}}}, {.X = 1, .Y = 1})
              .size()
          == 1);
  }

  /**
   *  An arc or an ellipse swept round a million times has two feet a turn, and the feet of its
   *  first 65,536 turns are given, as Clip seeks its crossings, the cap documented where the
   *  documentation said every turn; swept 1e300 times, its turns are closer in its parameter than a
   *  few steps of a double, and its feet come once.
   */
  TEST_CASE("nearest: the feet of an arc or ellipse swept past 65,536 turns, of its first 65,536")
  {
    constexpr std::size_t Turns = 65'536;
    const Point           point = {.X = 3, .Y = 0.5};
    const Arc arc = {.Center = {}, .Radius = 1, .Start = Degrees(0), .Sweep = Degrees(360e6)};
    const std::vector<Place> feet = Perpendiculars(arc, point);
    REQUIRE(feet.size() == 2 * Turns);
    CHECK(feet.back().Along <= static_cast<double>(Turns) / 1e6);
    const Ellipse ellipse = {.Map   = {.X = {.X = 2, .Y = 0}, .Y = {.X = 0, .Y = 1}, .Origin = {}},
                             .Start = Degrees(0),
                             .Sweep = Degrees(360e6)};
    CHECK(Perpendiculars(ellipse, point).size() == 2 * Turns);
    const Arc endless = {.Center = {}, .Radius = 1, .Start = Degrees(0), .Sweep = Degrees(360e300)};
    CHECK(Perpendiculars(endless, point).size() == 1);
  }

  /**
   *  An ellipse thin past what doubles hold, its axes 1e279 and 1e255 long nearly along one line:
   *  its quarters' roots were 1,126 spots of rounding noise in a turn, each passed on every turn,
   *  89 feet in one turn, 620,001 in 10,000 turns and past 2 GB in a million. At most four spots a
   *  turn are kept now, those nearest square to the way to the point. And the robustness audit's
   *  ellipse with an axis of 1e300: four feet a turn, its first 65,536 turns.
   */
  TEST_CASE("nearest: an ellipse thin past what doubles hold has four feet a turn at most")
  {
    Ellipse     thin  = {.Map   = {.X      = {.X = -7.92994e+07, .Y = 7.18053e+279},
                                   .Y      = {.X = -6.82878e+21, .Y = -1.56768e+255},
                                   .Origin = {.X = 0.754822, .Y = 0.724193}},
                         .Start = Degrees(0),
                         .Sweep = Degrees(360)};
    const Point asked = {.X = -0.344748, .Y = -43.4504};
    CHECK(Perpendiculars(thin, asked).size() <= 4);
    thin.Sweep = Degrees(360e4);
    CHECK(Perpendiculars(thin, asked).size() <= 40'001);
    const Ellipse wide = {.Map   = {.X = {.X = 1, .Y = 0}, .Y = {.X = 1e300, .Y = 1}, .Origin = {}},
                          .Start = Degrees(0),
                          .Sweep = Degrees(360e6)};
    CHECK(Perpendiculars(wide, {.X = 3, .Y = 0.5}).size() <= (4 * 65'536) + 1);
  }
} // namespace Geoxy
