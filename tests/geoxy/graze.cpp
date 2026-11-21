/***************************************************************************************************
 *  @file      tests/geoxy/graze.cpp
 *  @brief     Tests that ready lines and arcs that cross shallowly, nearly touch or nearly miss
 *             meet as the curves themselves do, to the bit, and that the grazing pass tells them.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "graze.hpp"

#include "abut.hpp"
#include "draws.hpp"
#include "real.hpp"
#include "sane.hpp"
#include "thorough.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  How a pair met: whether made ready it met as the curves themselves, either first, and
     *  whether the grazing pass told it, either first.
     */
    struct Glanced
    {
      bool Same = true; /**< Whether ready and plain agree to the bit.               */
      int  Told = 0;    /**< How many of the two orders either way of the pass told. */
    };

    /**
     *  Two curves met one way round, made ready and not: whether the marks are the meeting's hits
     *  to the bit, with its stretches and its Spent; and whether the grazing pass, asked directly,
     *  told it, and told the same.
     *
     *  @param  curves    The two curves, the first first.
     *  @param  tolerance The tolerance.
     *  @return           How they met.
     */
    auto GlanceOnce(const Curves &curves, const Tolerance &tolerance) -> Glanced
    {
      const Prepared one{curves[0].get()};
      const Prepared other{curves[1].get()};
      Intersections  meeting;
      Marking        marking;
      Intersect(curves, tolerance, meeting);
      Intersect({one, other}, tolerance, marking);
      const auto mark = [](const Hit &lhs, const Mark &rhs) -> bool
      {
        return rhs == Mark{.At = lhs.At, .Kind = lhs.Kind};
      };
      Glanced         met  = {.Same = std::ranges::equal(meeting.Hits, marking.Marks, mark)
                                   && marking.Overlaps.size() == meeting.Overlaps.size()
                                   && meeting.Spent == marking.Spent,
                              .Told = 0};
      const Tolerance sane = Sane(tolerance);
      const auto      same = [&meeting, &mark](const Abutment &told) -> bool
      {
        return meeting.Overlaps.size() == told.Stretches
            && std::ranges::equal(meeting.Hits, std::span{told.Marks}.first(told.Count), mark);
      };
      if (std::holds_alternative<Arc>(curves[0].get())
          || std::holds_alternative<Arc>(curves[1].get()))
      {
        const auto tips = [](const Curve &arc) -> std::array<Point, 2>
        {
          return {At(arc, 0.0), At(arc, 1.0)};
        };
        const Glancing::Held          first   = Glancing::Hold(curves[0], tips(curves[0]));
        const Glancing::Held          second  = Glancing::Hold(curves[1], tips(curves[1]));
        const std::optional<Abutment> hinged  = Glancing::Hinged({first, second}, sane);
        met.Told                             += static_cast<int>(hinged.has_value());
        met.Same = met.Same && (!hinged || (!meeting.Spent && same(*hinged)));
      }
      const std::optional<double> length = Reachable(curves, sane);
      if (!length)
      {
        return met;
      }
      const Tolerance usable = {.Length = *length, .Radians = sane.Radians};
      for (const std::optional<Abutment> &told :
           {Glancing::Met({one, other}, usable), Glancing::Steeply({one, other}, usable)})
      {
        if (told)
        {
          met.Told = 1;
          met.Same =
              met.Same
              && meeting.Overlaps.empty()
              && std::ranges::equal(meeting.Hits, std::span{told->Marks}.first(told->Count), mark);
        }
      }
      return met;
    }

    /**
     *  Two curves met either first, made ready and not.
     *
     *  @param  first     One curve.
     *  @param  second    The other.
     *  @param  tolerance The tolerance.
     *  @return           How they met.
     */
    auto GlanceBoth(const Curve &first, const Curve &second, const Tolerance &tolerance) -> Glanced
    {
      const Glanced one   = GlanceOnce({first, second}, tolerance);
      const Glanced other = GlanceOnce({second, first}, tolerance);
      return {.Same = one.Same && other.Same, .Told = one.Told + other.Told};
    }

    /**
     *  A random arc of a circle: whole, or of a sweep either way, its angles in radians or degrees.
     *
     *  @param  draws  The draws.
     *  @param  center Its center.
     *  @param  radius Its radius.
     *  @return        The arc.
     */
    auto Penned(Draws &draws, Point center, double radius) -> Arc
    {
      const auto pick = [&draws](std::uint64_t count) -> std::uint64_t
      {
        return draws.Next() % count;
      };
      const double start = draws.Near(0.0, Pi);
      const double sweep = pick(4) == 0 ? 2.0 * Pi : draws.Near(0.0, 1.9 * Pi);
      const bool   whole = Exactly(sweep, 2.0 * Pi);
      return {.Center = center,
              .Radius = radius,
              .Start  = pick(2) == 0 ? Radians(start) : Degrees(start * 180.0 / Pi),
              .Sweep  = pick(2) == 0 || whole ? Radians(sweep) : Degrees(sweep * 180.0 / Pi)};
    }

    /**
     *  The gaps, in tolerances, a line or a circle is put at from a circle where they come nearest:
     *  through every band of Intersect's, either side of each bound and on it.
     */
    constexpr std::array Gaps = {-1e6,  -4096.0, -64.0, -8.0, -4.5,  -4.0, -3.0, -2.25,
                                 -2.0,  -1.75,   -1.0,  -0.5, -0.25, 0.0,  0.25, 0.5,
                                 0.875, 1.0,     1.125, 2.0,  3.0,   4.0,  6.0,  8.0};

    /**
     *  A scale the plane is taken to, as a power of two.
     *
     *  @param  draws The draws.
     *  @return       The scale.
     */
    auto Scaled(Draws &draws) -> double
    {
      constexpr std::array Powers = {-40, -20, 0, 0, 0, 20, 40};
      return std::ldexp(1.0, Powers.at(draws.Next() % Powers.size()));
    }
  } // namespace

  TEST_CASE("graze: lines near a circle at every band's gap meet as Intersect answers, to the bit")
  {
    Draws draws{0x5EED1};
    int   told  = 0;
    int   pairs = 0;
    for (int round = 0; round < Rounds(3000); ++round)
    {
      const double scale  = Scaled(draws);
      const double radius = draws.Next() % 8 == 0 ? 1e6 : draws.Near(20.0, 19.0);
      const Point  center = {.X = draws.Near(0.0, 1e3) * scale, .Y = draws.Near(0.0, 1e3) * scale};
      const double length =
          1e-9 * scale * (1.0 + (((std::abs(center.X) + std::abs(center.Y)) / scale) * 1e-3));
      const double  gap     = Gaps.at(draws.Next() % Gaps.size());
      const double  bearing = draws.Near(0.0, Pi);
      const Vector  normal  = {.X = std::cos(bearing), .Y = std::sin(bearing)};
      const Vector  along   = Perpendicular(normal);
      const double  reach   = (radius * scale) + (gap * length);
      const Point   foot    = center + (normal * reach);
      const double  before  = radius * scale * draws.Near(1.0, 0.9);
      const double  after   = radius * scale * draws.Near(1.0, 0.9);
      const Line    line    = {.Start = foot - (along * before), .End = foot + (along * after)};
      const Arc     arc     = Penned(draws, center, radius * scale);
      const Glanced met     = GlanceBoth(line, arc, {.Length = length, .Radians = 1e-9});
      CHECK(met.Same);
      told  += met.Told;
      pairs += 2;
    }
    CHECK(told > pairs / 3);
  }

  TEST_CASE("graze: arcs of bulged polygons joined end to end, held, meet as Intersect answers")
  {
    Draws draws{0x5EED'4848};
    int   told   = 0;
    int   joints = 0;
    for (int round = 0; round < Rounds(600); ++round)
    {
      /**
       *  A wavy ring of 3 to 40 corners at a scale from 2^-40 to 2^40, about the origin or far from
       *  it, every side bent alike: by 0.15, gently, both ways, or to near a half turn.
       */
      const double      scale = Scaled(draws) * std::ldexp(1.0, static_cast<int>(draws.Next() % 9));
      const double      far   = draws.Next() % 2 == 0 ? 0.0 : draws.Near(0.0, 1e6) * scale;
      const std::size_t corners = 3 + (draws.Next() % 38);
      const double      phase   = draws.Near(0.0, Pi);
      const double      bulge =
          std::array{0.15, draws.Near(0.0, 0.3), draws.Near(0.0, 0.99), -0.15, 1e-6}.at(draws.Next()
                                                                                        % 5);
      std::vector<Curve> sides;
      for (std::size_t index = 0; index < corners; ++index)
      {
        const auto corner = [&](std::size_t place) -> Point
        {
          const double angle =
              2.0 * Pi * static_cast<double>(place % corners) / static_cast<double>(corners);
          const double radius = scale * (1.0 + (0.3 * std::sin((3.0 * angle) + phase)));
          return {.X = far + (radius * std::cos(angle)), .Y = far + (radius * std::sin(angle))};
        };
        sides.emplace_back(std::visit([](const auto &shape) -> Curve { return shape; },
                                      Bent(corner(index), corner(index + 1), bulge)));
      }
      const double length = scale * std::pow(10.0, -static_cast<double>(2 + (draws.Next() % 11)));
      for (std::size_t index = 0; index < corners; ++index)
      {
        const Glanced met = GlanceOnce({sides[index], sides[(index + 1) % corners]},
                                       {.Length = length, .Radians = 1e-9});
        CHECK(met.Same);
        told   += static_cast<int>(met.Told > 0);
        joints += 1;
      }
    }
    CHECK(told > joints / 2);
  }

  TEST_CASE("graze: curves held that lie apart, or that the bands leave to Intersect")
  {
    const Tolerance sane   = Sane({.Length = 1e-9, .Radians = 1e-9});
    const auto      hinged = [](const Curve &first, const Curve &second,
                                const Tolerance &tolerance) -> std::optional<Abutment>
    {
      const Glancing::Held one   = Glancing::Hold(first, {At(first, 0.0), At(first, 1.0)});
      const Glancing::Held other = Glancing::Hold(second, {At(second, 0.0), At(second, 1.0)});
      return Glancing::Hinged({one, other}, tolerance);
    };
    const Curve near  = Arc{.Center = {.X = 0.0, .Y = 0.0},
                            .Radius = 1.0,
                            .Start  = Degrees(0.0),
                            .Sweep  = Degrees(90.0)};
    const Curve far   = Arc{.Center = {.X = 9.0, .Y = 0.0},
                            .Radius = 1.0,
                            .Start  = Degrees(0.0),
                            .Sweep  = Degrees(90.0)};
    const Curve line  = Line{.Start = {.X = 1.0, .Y = -1.0}, .End = {.X = 1.0, .Y = 1.0}};
    const Curve other = Line{.Start = {.X = 0.0, .Y = -1.0}, .End = {.X = 2.0, .Y = 1.0}};
    /**
     *  Arcs whose circles' boxes lie far apart meet in nothing, as Intersect finds; two lines are
     *  left to Intersect, as are an arc whose box runs past every double, or whose center and
     *  radius sum past it, an arc of no sweep, a tolerance past Coarse and one, at a grain finer
     *  still, under Fine.
     */
    const Curve above = Arc{.Center = {.X = 0.0, .Y = 9.0},
                            .Radius = 1.0,
                            .Start  = Degrees(0.0),
                            .Sweep  = Degrees(90.0)};
    for (const auto &[first, second] : {std::pair{&near, &far}, std::pair{&far, &near},
                                        std::pair{&near, &above}, std::pair{&above, &near}})
    {
      const std::optional<Abutment> apart = hinged(*first, *second, sane);
      CHECK(apart.has_value());
      CHECK(apart.value_or(Abutment{.Count = 1}).Count == 0);
      CHECK(Intersect(*first, *second, sane).Hits.empty());
    }
    CHECK_FALSE(hinged(line, other, sane).has_value());
    const Curve huge   = Arc{.Center = {.X = 1e308, .Y = -1e308},
                             .Radius = 1e308,
                             .Start  = Degrees(0.0),
                             .Sweep  = Degrees(90.0)};
    const Curve beyond = Arc{.Center = {.X = 1e308, .Y = 1e308},
                             .Radius = 1e308,
                             .Start  = Degrees(0.0),
                             .Sweep  = Degrees(90.0)};
    const Curve none   = Arc{.Center = {.X = 0.5, .Y = 0.0},
                             .Radius = 1.0,
                             .Start  = Degrees(0.0),
                             .Sweep  = Degrees(0.0)};
    CHECK_FALSE(hinged(huge, near, sane).has_value());
    CHECK_FALSE(hinged(near, beyond, sane).has_value());
    CHECK_FALSE(hinged(none, near, sane).has_value());
    CHECK_FALSE(hinged(line, none, sane).has_value());
    CHECK_FALSE(hinged(near, line, {.Length = 0x1p60, .Radians = 1e-9}).has_value());
    const Curve tiny  = Arc{.Center = {.X = 0.0, .Y = 0.0},
                            .Radius = 1e-20,
                            .Start  = Degrees(0.0),
                            .Sweep  = Degrees(90.0)};
    const Curve small = Arc{.Center = {.X = 1e-20, .Y = 0.0},
                            .Radius = 1e-20,
                            .Start  = Degrees(90.0),
                            .Sweep  = Degrees(90.0)};
    CHECK_FALSE(hinged(tiny, small, {.Length = 1e-40, .Radians = 1e-9}).has_value());
  }

  TEST_CASE("graze: lines grazing a circle to 1e-12 of its radius cross where Intersect answers")
  {
    Draws draws{0x5EED2};
    int   told  = 0;
    int   pairs = 0;
    for (int round = 0; round < Rounds(3000); ++round)
    {
      const double  scale  = Scaled(draws);
      const double  radius = (draws.Next() % 4 == 0 ? 1e6 : draws.Near(50.0, 49.0)) * scale;
      const Point   center = {.X = draws.Near(0.0, 100.0) * scale,
                              .Y = draws.Near(0.0, 100.0) * scale};
      const double  depth  = radius * std::pow(10.0, -static_cast<double>(2 + (draws.Next() % 11)));
      const double  bearing = draws.Near(0.0, Pi);
      const Vector  normal  = {.X = std::cos(bearing), .Y = std::sin(bearing)};
      const Point   foot    = center + (normal * (radius - depth));
      const double  half    = std::sqrt(2.0 * radius * depth);
      const Vector  along   = Perpendicular(normal) * (half * draws.Near(2.0, 0.9));
      const Line    line    = {.Start = foot - along, .End = foot + along};
      const Arc     arc     = Penned(draws, center, radius);
      const double  size    = std::max({std::abs(center.X), std::abs(center.Y), radius});
      const Glanced met     = GlanceBoth(line, arc, {.Length = 1e-12 * size, .Radians = 1e-9});
      CHECK(met.Same);
      told  += met.Told;
      pairs += 2;
    }
    CHECK(told > pairs / 2);
  }

  TEST_CASE(
      "graze: circles near each other at every band's gap meet as Intersect answers, to the bit")
  {
    Draws draws{0x5EED3};
    int   told  = 0;
    int   pairs = 0;
    for (int round = 0; round < Rounds(3000); ++round)
    {
      const double scale  = Scaled(draws);
      const double inner  = draws.Near(20.0, 19.0) * scale;
      const double outer  = (draws.Next() % 8 == 0 ? 1e6 : draws.Near(20.0, 19.0)) * scale;
      const Point  center = {.X = draws.Near(0.0, 1e3) * scale, .Y = draws.Near(0.0, 1e3) * scale};
      const double length =
          1e-9 * scale * (1.0 + (((std::abs(center.X) + std::abs(center.Y)) / scale) * 1e-3));
      const double gap    = Gaps.at(draws.Next() % Gaps.size());
      const bool   within = draws.Next() % 2 == 0;
      const double apart =
          within ? std::abs(outer - inner) - (gap * length) : inner + outer + (gap * length);
      const double bearing = draws.Near(0.0, Pi);
      const Point other = center + (Vector{.X = std::cos(bearing), .Y = std::sin(bearing)} * apart);
      const Arc   first = Penned(draws, center, inner);
      const Arc   second = Penned(draws, other, outer);
      const Glanced met  = GlanceBoth(first, second, {.Length = length, .Radians = 1e-9});
      CHECK(met.Same);
      told  += met.Told;
      pairs += 2;
    }
    CHECK(told > pairs / 3);
  }

  TEST_CASE(
      "graze: circles crossing shallowly to 1e-12 of their radius cross where Intersect answers")
  {
    Draws draws{0x5EED4};
    int   told  = 0;
    int   pairs = 0;
    for (int round = 0; round < Rounds(3000); ++round)
    {
      const double scale  = Scaled(draws);
      const double inner  = draws.Near(50.0, 49.0) * scale;
      const double outer  = (draws.Next() % 4 == 0 ? 1e6 : draws.Near(50.0, 49.0)) * scale;
      const Point  center = {.X = draws.Near(0.0, 100.0) * scale,
                             .Y = draws.Near(0.0, 100.0) * scale};
      const double depth =
          std::min(inner, outer) * std::pow(10.0, -static_cast<double>(2 + (draws.Next() % 11)));
      const bool   within  = draws.Next() % 2 == 0 && std::abs(outer - inner) > 4.0 * depth;
      const double apart   = within ? std::abs(outer - inner) + depth : inner + outer - depth;
      const double bearing = draws.Near(0.0, Pi);
      const Point other = center + (Vector{.X = std::cos(bearing), .Y = std::sin(bearing)} * apart);
      const Arc   first = Penned(draws, center, inner);
      const Arc   second = Penned(draws, other, outer);
      const double  size = std::max({std::abs(center.X), std::abs(center.Y), inner, outer});
      const Glanced met  = GlanceBoth(first, second, {.Length = 1e-12 * size, .Radians = 1e-9});
      CHECK(met.Same);
      told  += met.Told;
      pairs += 2;
    }
    CHECK(told > pairs / 3);
  }

  TEST_CASE("graze: touches, misses and shallow crossings told, and edges left to Intersect")
  {
    const Tolerance tolerance = {.Length = 0x1p-30, .Radians = 1e-9};
    const Arc       circle    = {.Center = {.X = 3.0, .Y = 2.0},
                                 .Radius = 5.0,
                                 .Start  = Radians(0.0),
                                 .Sweep  = Radians(2.0 * Pi)};
    const auto      level     = [](double y) -> Line
    {
      return {.Start = {.X = -9.0, .Y = y}, .End = {.X = 15.0, .Y = y}};
    };
    const auto told = [&tolerance](const Curve &first, const Curve &second) -> bool
    {
      const Glanced met = GlanceBoth(first, second, tolerance);
      CHECK(met.Same);
      return met.Told == 2;
    };
    /**
     *  A level line half a tolerance off the top, two tolerances off, and a
     *  thousand-and-twenty-fourth of radius into the circle: a touch, nothing, and two shallow
     *  crossings.
     */
    CHECK(told(level(7.0 + 0x1p-31), circle));
    CHECK(told(level(7.0 + 0x1p-29), circle));
    CHECK(told(level(7.0 - (5.0 / 1024.0)), circle));
    /**
     *  On a bound of a band, and a line whose ends lie on the circle: left to Intersect.
     */
    CHECK(!told(level(7.0 + (0x1p-30 * 7.0 / 8.0)), circle));
    CHECK(!told(Line{.Start = {.X = 3.0, .Y = 7.0}, .End = {.X = 13.0, .Y = 7.0}}, circle));
    /**
     *  A circle half a tolerance outside another, and one sixteen tolerances inside crossing it
     *  shallowly.
     */
    const Arc beside = {.Center = {.X = 3.0, .Y = 2.0 + 8.0 + 0x1p-31},
                        .Radius = 3.0,
                        .Start  = Degrees(-120.0),
                        .Sweep  = Degrees(60.0)};
    CHECK(told(circle, beside));
    const Arc shallow = {.Center = {.X = 3.0, .Y = 2.0 + 2.0 + 0x1p-26},
                         .Radius = 3.0,
                         .Start  = Degrees(0.0),
                         .Sweep  = Degrees(360.0)};
    CHECK(GlanceBoth(circle, shallow, tolerance).Same);
  }

  TEST_CASE("graze: lines and circles touching exactly, or a step of a double off, told")
  {
    const Tolerance tolerance = {.Length = 0x1p-30, .Radians = 1e-9};
    for (const double scale : {1.0, 0x1p-20, 0x1p20})
    {
      const Arc  circle = {.Center = {.X = 3.0 * scale, .Y = 2.0 * scale},
                           .Radius = 5.0 * scale,
                           .Start  = Radians(0.0),
                           .Sweep  = Radians(2.0 * Pi)};
      const Arc  upper  = {.Center = circle.Center,
                           .Radius = circle.Radius,
                           .Start  = Degrees(20.0),
                           .Sweep  = Degrees(140.0)};
      const auto told   = [&tolerance, scale](const Curve &first, const Curve &second) -> bool
      {
        const Glanced met =
            GlanceBoth(first, second, {.Length = tolerance.Length * scale, .Radians = 1e-9});
        CHECK(met.Same);
        return met.Told == 2;
      };
      const double top   = 7.0 * scale;
      const auto   level = [scale](double y) -> Line
      {
        return {.Start = {.X = -9.0 * scale, .Y = y}, .End = {.X = 15.0 * scale, .Y = y}};
      };
      /**
       *  A level line on the top, a step of a double above and below it; a slanted one on the
       *  circle's point (6, 6), its way (4, -3).
       */
      for (const Arc &arc : {circle, upper})
      {
        CHECK(told(level(top), arc));
        CHECK(told(level(std::nextafter(top, 2.0 * top)), arc));
        CHECK(told(level(std::nextafter(top, 0.0)), arc));
        const Point  touch = {.X = 6.0 * scale, .Y = 6.0 * scale};
        const Vector way   = {.X = 4.0 * scale, .Y = -3.0 * scale};
        CHECK(told(Line{.Start = touch - way, .End = touch + (way * 0.75)}, arc));
        /**
         *  Circles touching outside and inside on the line of centers through (6, 6).
         */
        const Arc outside = {.Center = {.X = 9.0 * scale, .Y = 10.0 * scale},
                             .Radius = 5.0 * scale,
                             .Start  = Degrees(0.0),
                             .Sweep  = Degrees(360.0)};
        const Arc inside  = {.Center = {.X = 4.5 * scale, .Y = 4.0 * scale},
                             .Radius = 2.5 * scale,
                             .Start  = Degrees(-90.0),
                             .Sweep  = Degrees(270.0)};
        CHECK(told(arc, outside));
        CHECK(told(inside, arc));
      }
    }
  }

  TEST_CASE("graze: arcs coiled, ends at a crossing and crossings on an axis, left to Intersect")
  {
    const Tolerance tolerance = {.Length = 0x1p-30, .Radians = 1e-9};
    const Arc       circle    = {.Center = {.X = 3.0, .Y = 2.0},
                                 .Radius = 5.0,
                                 .Start  = Radians(0.0),
                                 .Sweep  = Radians(2.0 * Pi)};
    const auto      left      = [&tolerance](const Curve &first, const Curve &second) -> bool
    {
      const Glanced met = GlanceBoth(first, second, tolerance);
      CHECK(met.Same);
      return met.Told == 0;
    };
    /**
     *  An arc swept round one and a half turns, against a line and a circle grazing it.
     */
    const Arc coiled = {.Center = {.X = 3.0, .Y = 2.0},
                        .Radius = 5.0,
                        .Start  = Radians(0.0),
                        .Sweep  = Radians(3.0 * Pi)};
    CHECK(left(coiled, Line{.Start = {.X = -9.0, .Y = 7.0 + 0x1p-31},
                            .End   = {.X = 15.0, .Y = 7.0 + 0x1p-31}}));
    CHECK(left(coiled, Arc{.Center = {.X = 3.0, .Y = 10.0 + 0x1p-31},
                           .Radius = 3.0,
                           .Start  = Radians(0.0),
                           .Sweep  = Radians(2.0 * Pi)}));
    /**
     *  A circle within another, crossing it sixty-four tolerances deep, an arc of it ending at one
     *  of the crossings.
     */
    const Point  centre = {.X = 3.0, .Y = 4.0 + 0x1p-24};
    const double apart  = centre.Y - circle.Center.Y;
    const double along  = ((apart * apart) + 9.0 - 25.0) / (2.0 * apart);
    const double rise   = std::sqrt(9.0 - (along * along));
    CHECK(left(circle, Arc{.Center = centre,
                           .Radius = 3.0,
                           .Start  = Radians(0.0),
                           .Sweep  = Radians(std::atan2(-along, rise))}));
    /**
     *  Two circles of one radius crossing shallowly at x = 0, which no bound tells.
     */
    const double offset = 5.0 - 0x1p-24;
    CHECK(left(Arc{.Center = {.X = -offset, .Y = 0.0},
                   .Radius = 5.0,
                   .Start  = Radians(0.0),
                   .Sweep  = Radians(2.0 * Pi)},
               Arc{.Center = {.X = offset, .Y = 0.0},
                   .Radius = 5.0,
                   .Start  = Radians(0.0),
                   .Sweep  = Radians(2.0 * Pi)}));
  }

  TEST_CASE(
      "graze: lines at a steep angle ending at, short of and through an arc, as Intersect answers")
  {
    const Tolerance tolerance = {.Length = 0x1p-30, .Radians = 1e-9};
    const double    length    = tolerance.Length;
    const Point     centre    = {.X = 3.0, .Y = 2.0};
    const Arc       circle    = {
        .Center = centre, .Radius = 5.0, .Start = Radians(0.0), .Sweep = Radians(2.0 * Pi)};
    const Arc part = {
        .Center = centre, .Radius = 5.0, .Start = Radians(0.3), .Sweep = Radians(1.0)};
    const Vector way  = {.X = 0.6, .Y = 0.8};
    const auto   same = [&tolerance](const Curve &first, const Curve &second) -> int
    {
      const Glanced met = GlanceBoth(first, second, tolerance);
      CHECK(met.Same);
      return met.Told;
    };
    int told = 0;
    /**
     *  From the center out, stopping short of the circle by a quarter, two and twenty tolerances,
     *  and ending on it; and the other way round.
     */
    for (const double shy : {0.25, 2.0, 20.0, 0.0})
    {
      const Point end  = centre + (way * (5.0 - (shy * length)));
      told            += same(Line{.Start = centre, .End = end}, circle);
      told            += same(Line{.Start = end, .End = centre}, circle);
      told            += same(Line{.Start = centre, .End = end}, part);
    }
    /**
     *  Stopping short of an arc's circle where the arc does not reach, and past the circle by two
     *  tolerances.
     */
    const Vector aside = {.X = -0.8, .Y = 0.6};
    told += same(Line{.Start = centre, .End = centre + (aside * (5.0 - (2.0 * length)))}, part);
    told += same(Line{.Start = centre, .End = centre + (aside * (5.0 - (0.25 * length)))}, part);
    told += same(Line{.Start = centre, .End = centre + (way * (5.0 + (2.0 * length)))}, circle);
    /**
     *  Across the circle, stopping a quarter of a tolerance short of its far side.
     */
    const Point enter = centre - (way * 6.0);
    told += same(Line{.Start = enter, .End = centre + (way * (5.0 - (0.25 * length)))}, circle);
    told += same(Line{.Start = enter, .End = centre + (way * (5.0 - (0.25 * length)))}, part);
    /**
     *  Through a whole circle's seam, and two tolerances beside it; through an arc's start.
     */
    const Point seam   = {.X = 8.0, .Y = 2.0};
    told              += same(Line{.Start = seam + (way * 0.5), .End = seam - (way * 6.0)}, circle);
    told              += same(Line{.Start = seam + (aside * (2.0 * length)) + (way * 0.5),
                                   .End   = seam + (aside * (2.0 * length)) - (way * 6.0)},
                              circle);
    const Point start  = centre + (Vector{.X = std::cos(0.3), .Y = std::sin(0.3)} * 5.0);
    told += same(Line{.Start = start + (aside * 0.5), .End = start - (aside * 6.0)}, part);
    /**
     *  Cutting the circle at a sine a hair past Steep's, both ends far.
     */
    const double off   = 5.0 * std::sqrt((63.0 / 64.0) * (1.0 - 0x1p-21));
    const Point  steep = centre + (aside * off);
    told += same(Line{.Start = steep - (way * 7.0), .End = steep + (way * 7.0)}, circle);
    /**
     *  Through the seam of a circle at the origin, its crossings on the axes.
     */
    told += same(Line{.Start = {.X = -1.0, .Y = -6.0}, .End = {.X = 6.0, .Y = 1.0}},
                 Arc{.Center = {.X = 0.0, .Y = 0.0},
                     .Radius = 5.0,
                     .Start  = Radians(0.0),
                     .Sweep  = Radians(2.0 * Pi)});
    CHECK(told > 0);
  }

  TEST_CASE("graze: arcs crossing steeply, an end at, near or off a crossing, as Intersect answers")
  {
    const Tolerance tolerance = {.Length = 0x1p-30, .Radians = 1e-9};
    const double    length    = tolerance.Length;
    int             told      = 0;
    const auto      same      = [&tolerance, &told](const Curve &first, const Curve &second) -> void
    {
      const Glanced met  = GlanceBoth(first, second, tolerance);
      told              += met.Told;
      CHECK(met.Same);
    };
    /**
     *  Arcs of radius 5 about the origin and (-5, 5), crossing at (0, 5) and (-5, 0): the first
     *  ending at the second crossing, off the second arc, as given and swept to the crossing.
     */
    const Arc right = {.Center = {.X = 0.0, .Y = 0.0},
                       .Radius = 5.0,
                       .Start  = Radians(Pi),
                       .Sweep  = Radians(-2.498091544796509)};
    const Arc left  = {.Center = {.X = -5.0, .Y = 5.0},
                       .Radius = 5.0,
                       .Start  = Radians(Pi / 2.0),
                       .Sweep  = Radians(-2.498091544796509)};
    same(right, left);
    same(Arc{.Center = {.X = 0.0, .Y = 0.0},
             .Radius = 5.0,
             .Start  = Degrees(180.0),
             .Sweep  = Degrees(-180.0)},
         Arc{.Center = {.X = -5.0, .Y = 5.0},
             .Radius = 5.0,
             .Start  = Degrees(90.0),
             .Sweep  = Degrees(-90.0)});
    /**
     *  Circles of radius 5 and 3 crossing steeply; an arc of the second ending at a crossing, a
     *  quarter, two and twenty tolerances short of it or past it, and starting there; against the
     *  whole first circle, its seam far or at the other crossing, and an arc of it.
     */
    const Point  centre = {.X = 4.0, .Y = 1.0};
    const double apart  = std::sqrt(17.0);
    const double along  = ((apart * apart) + 9.0 - 25.0) / (2.0 * apart);
    const double angle  = std::atan2(-1.0, -4.0) + std::acos(along / 3.0);
    const Arc    circle = {.Center = {.X = 0.0, .Y = 0.0},
                           .Radius = 5.0,
                           .Start  = Radians(0.0),
                           .Sweep  = Radians(2.0 * Pi)};
    const Arc    part   = {.Center = {.X = 0.0, .Y = 0.0},
                           .Radius = 5.0,
                           .Start  = Radians(-1.0),
                           .Sweep  = Radians(2.5)};
    for (const double shy : {0.0, 0.25, 2.0, 20.0, -0.25, -2.0})
    {
      const double sweep = angle - 0.5 - (shy * length / 3.0);
      for (const Arc *round : {&circle, &part})
      {
        same(*round,
             Arc{.Center = centre, .Radius = 3.0, .Start = Radians(0.5), .Sweep = Radians(sweep)});
        same(*round, Arc{.Center = centre,
                         .Radius = 3.0,
                         .Start  = Radians(angle + (shy * length / 3.0)),
                         .Sweep  = Radians(-1.5)});
      }
    }
    const Point crossing = centre + (Vector{.X = std::cos(angle), .Y = std::sin(angle)} * 3.0);
    same(Arc{.Center = {.X = 0.0, .Y = 0.0},
             .Radius = 5.0,
             .Start  = Radians(std::atan2(crossing.Y, crossing.X)),
             .Sweep  = Radians(2.0 * Pi)},
         Arc{.Center = centre, .Radius = 3.0, .Start = Radians(0.0), .Sweep = Radians(2.0 * Pi)});
    same(Arc{.Center = {.X = 0.0, .Y = 0.0},
             .Radius = 5.0,
             .Start  = Radians(std::atan2(crossing.Y, crossing.X) + (3.0 * length / 5.0)),
             .Sweep  = Radians(2.0 * Pi)},
         Arc{.Center = centre, .Radius = 3.0, .Start = Radians(1.0), .Sweep = Radians(1.0)});
    /**
     *  At random, about circles crossing at every angle past Steep's, arcs ending at a crossing, a
     *  hair from it, or anywhere.
     */
    Draws draws{0x5EED'2047};
    for (int round = 0; round < Rounds(4000); ++round)
    {
      const double scale  = Scaled(draws);
      const double radius = draws.Near(1.0, 0.5) * scale;
      const double other  = draws.Near(1.0, 0.5) * scale;
      const double reach  = draws.Near(0.5, 0.45) * (radius + other);
      const double gap    = std::max(std::abs(radius - other) * 1.05, reach);
      const double bend   = draws.Near(0.0, Pi);
      const Point  near   = {.X = gap * std::cos(bend), .Y = gap * std::sin(bend)};
      const double cosine =
          ((gap * gap) + (radius * radius) - (other * other)) / (2.0 * gap * radius);
      if (!(std::abs(cosine) < 1.0))
      {
        continue;
      }
      const double meet  = bend + ((draws.Next() % 2 == 0 ? 1.0 : -1.0) * std::acos(cosine));
      const double hair  = std::array{0.0, 0.0, 0.25, 2.0, 1e3}.at(draws.Next() % 5) * length;
      const double start = meet + (hair / radius);
      const Arc    first = Penned(draws, {.X = 0.0, .Y = 0.0}, radius);
      const Arc    ended = {.Center = {.X = 0.0, .Y = 0.0},
                            .Radius = radius,
                            .Start  = Radians(start),
                            .Sweep  = Radians(draws.Near(0.0, 3.0))};
      same(draws.Next() % 2 == 0 ? first : ended, Penned(draws, near, other));
    }
    CHECK(told > 0);
  }

} // namespace Geoxy
