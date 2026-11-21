/***************************************************************************************************
 *  @file      tests/geoxy/strand.cpp
 *  @brief     Tests for a line, a ray or an endless line met with a spline of any degree, rational
 *             or not, by the kernel that holds the line's equation as a spline of its own: every
 *             answer the thorough way's to the bit, or left to it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "basis.hpp"
#include "draws.hpp"
#include "swift/lane.hpp"
#include "swift/met.hpp"
#include "swift/pieced.hpp"
#include "thorough.hpp"
#include "waves.hpp"

#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether two meetings are the same to the bit, hit for hit.
     *
     *  @param  one   One.
     *  @param  other The other.
     *  @return       Whether they are.
     */
    auto StrandSame(const Intersections &one, const Intersections &other) -> bool
    {
      const auto hit = [](const Hit &lhs, const Hit &rhs) -> bool
      {
        return lhs.At     == rhs.At
            && lhs.First  == rhs.First
            && lhs.Second == rhs.Second
            && lhs.Kind   == rhs.Kind;
      };
      return std::ranges::equal(one.Hits, other.Hits, hit)
          && one.Overlaps.size() == other.Overlaps.size()
          && one.Spent == other.Spent;
    }

    /**
     *  Two curves met by the swift kernels and by the thorough way.
     */
    struct StrandMet
    {
      bool Sure   = false;         /**< Whether the kernels were sure. */
      bool Agrees = true;          /**< Whether they found what the thorough
                                        way finds, to the bit.               */
      Intersections Swiftly  = {}; /**< What the kernels found.      */
      Intersections Thorough = {}; /**< What the thorough way found. */
    };

    /**
     *  Two curves met by the swift kernels, within the length Intersect would use, and by the
     *  thorough way where the kernels are sure.
     *
     *  @param  curves    The two curves.
     *  @param  tolerance The tolerance asked for.
     *  @return           What each found.
     */
    auto StrandBoth(const Curves &curves, const Tolerance &tolerance) -> StrandMet
    {
      StrandMet                   met;
      const std::optional<double> length = Reachable(curves, tolerance);
      if (!length)
      {
        return met;
      }
      met.Sure = Swift::Met(curves, {.Length = *length, .Radians = tolerance.Radians}, met.Swiftly);
      if (met.Sure)
      {
        Thorough(curves, tolerance, met.Thorough);
        met.Agrees = StrandSame(met.Swiftly, met.Thorough);
      }
      return met;
    }

    /**
     *  A spline of a degree from two to nine, now and then to fifteen: its poles and knots as the
     *  waves draw them, its weights none, or mild, within four times one another.
     *
     *  @param  draws The draws.
     *  @return       The spline, valid.
     */
    auto StrandSpline(Draws &draws) -> Spline
    {
      const int lowest = draws.Next() % 8 == 0 ? 9 : 2;
      Spline    spline = Tests::Splined(draws, lowest, lowest == 9 ? 6 : 7);
      if (draws.Next() % 2 == 0)
      {
        spline.Weights.clear();
        for (std::size_t index = 0; index < spline.Poles.size(); ++index)
        {
          spline.Weights.push_back(std::exp2(2.0 * draws.Near(0.0, 1.0)));
        }
      }
      return spline;
    }

    /**
     *  A straight curve: a line, a ray or an endless line, anywhere about the unit square, or now
     *  and then through a place of the spline, a knot or between them.
     *
     *  @param  draws  The draws.
     *  @param  spline The spline.
     *  @return        The curve.
     */
    auto StrandLine(Draws &draws, const Spline &spline) -> Curve
    {
      Point        start = {.X = draws.Near(0.0, 1.0), .Y = draws.Near(0.0, 1.0)};
      const Vector way   = {.X = draws.Near(0.0, 1.0), .Y = draws.Near(0.0, 1.0)};
      if (draws.Next() % 4 == 0)
      {
        const auto   degree = static_cast<std::size_t>(spline.Degree);
        const double knot  = spline.Knots[degree + (draws.Next() % (spline.Poles.size() - degree))];
        const double along = draws.Next() % 2 == 0 ? knot : knot + draws.Near(0.0, 0.1);
        start              = At(spline, along) - (way * 0.5);
      }
      switch (draws.Next() % 3)
      {
      case 0:
        return Line{.Start = start, .End = start + way};
      case 1:
        return Ray{.Start = start, .Way = way};
      default:
        return Xline{.Through = start, .Way = way};
      }
    }

    /**
     *  A spline of degree four through five poles, one span.
     *
     *  @param  poles The poles.
     *  @return       The spline.
     */
    auto StrandQuartic(const std::array<Point, 5> &poles) -> Spline
    {
      return {.Degree  = 4,
              .Knots   = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1},
              .Poles   = {poles.begin(), poles.end()},
              .Weights = {}};
    }
  } // namespace

  TEST_CASE("strand: lines and splines of any degree, rational or not, as met thoroughly")
  {
    const int          rounds = Rounds(40000);
    Draws              draws{0x57A4D5U};
    std::array<int, 4> count = {};
    for (int round = 0; round < rounds; ++round)
    {
      Spline       spline = StrandSpline(draws);
      const Curve  line   = StrandLine(draws, spline);
      const double scale =
          draws.Next() % 5 == 0 ? std::ldexp(1.0, static_cast<int>(draws.Next() % 81) - 40) : 1.0;
      const double near = draws.Next() % 4 == 0
                            ? std::pow(10.0, -3.0 - (9.0 * std::abs(draws.Near(0.0, 1.0))))
                            : 2e-10;
      for (Point &pole : spline.Poles)
      {
        pole = {.X = pole.X * scale, .Y = pole.Y * scale};
      }
      Curve shifted = line;
      if (auto *segment = std::get_if<Line>(&shifted))
      {
        *segment = {.Start = {.X = segment->Start.X * scale, .Y = segment->Start.Y * scale},
                    .End   = {.X = segment->End.X * scale, .Y = segment->End.Y * scale}};
      }
      else if (auto *ray = std::get_if<Ray>(&shifted))
      {
        ray->Start = {.X = ray->Start.X * scale, .Y = ray->Start.Y * scale};
      }
      else if (auto *xline = std::get_if<Xline>(&shifted))
      {
        xline->Through = {.X = xline->Through.X * scale, .Y = xline->Through.Y * scale};
      }
      REQUIRE(Valid(spline));
      const bool      swapped  = draws.Next() % 2 == 0;
      const Curve     curve    = spline;
      const Curves    curves   = swapped ? Curves{curve, shifted} : Curves{shifted, curve};
      const StrandMet met      = StrandBoth(curves, {.Length = near * scale, .Radians = 1e-10});
      const bool      skeined  = spline.Degree > 3 || Rational(spline);
      count[0]                += static_cast<int>(met.Sure);
      count[1]                += static_cast<int>(met.Sure && skeined);
      count[2]                += static_cast<int>(met.Sure && skeined && !met.Swiftly.Hits.empty());
      count[3]                += static_cast<int>(!met.Agrees);
      INFO("round " << round);
      CHECK(met.Agrees);
    }
    MESSAGE("sure of "
            << count[0]
            << " of "
            << rounds
            << ", past a cubic or rational "
            << count[1]
            << ", crossing in "
            << count[2]
            << ", unlike the thorough way in "
            << count[3]);
    CHECK(count[1] > rounds / 4);
    CHECK(count[2] > rounds / 10);
  }

  TEST_CASE("strand: what the kernel leaves to the thorough way")
  {
    const Tolerance                   tolerance = {.Length = 1e-9, .Radians = 1e-10};
    std::vector<std::array<Curve, 2>> pairs;
    /**
     *  A degree past fifteen.
     */
    Spline high = {.Degree = 16, .Knots = {}, .Poles = {}, .Weights = {}};
    for (int index = 0; index <= 16; ++index)
    {
      high.Poles.push_back({.X = index / 16.0, .Y = index % 2 == 0 ? -1.0 : 1.0});
    }
    high.Knots.assign(17, 0.0);
    high.Knots.insert(high.Knots.end(), 17, 1.0);
    pairs.push_back({Line{.Start = {.X = -1, .Y = 0.01}, .End = {.X = 2, .Y = 0.02}}, high});
    /**
     *  A fraction of knot differences below the normal doubles.
     */
    Spline faint = StrandQuartic({{{.X = 0, .Y = -1},
                                   {.X = 1, .Y =  1},
                                   {.X = 2, .Y = -1},
                                   {.X = 3, .Y =  1},
                                   {.X = 4, .Y = -1}}});
    faint.Poles.push_back({.X = 5, .Y = 1});
    faint.Knots = {0, 0, 0, 0, 0, 0x1p-1010, 1, 1, 1, 1, 1};
    REQUIRE(Valid(faint));
    pairs.push_back({Xline{.Through = {.X = 0, .Y = 0.1}, .Way = {.X = 1, .Y = 0.01}}, faint});
    /**
     *  An S whose chord runs along a short line at an angle that keeps the chord's ends clear of
     *  it.
     */
    const double angle  = 1.5e-9;
    const auto   turned = [angle](Point point) -> Point
    {
      return {.X = (point.X * std::cos(angle)) - (point.Y * std::sin(angle)),
              .Y = (point.X * std::sin(angle)) + (point.Y * std::cos(angle))};
    };
    pairs.push_back({Line{.Start = {.X = -0.5, .Y = 0}, .End = {.X = 0.5, .Y = 0}},
                     StrandQuartic({turned({.X = -20, .Y = 0}), turned({.X = -10, .Y = 6}),
                                    turned({.X = 0, .Y = 0}), turned({.X = 10, .Y = -6}),
                                    turned({.X = 20, .Y = 0})})});
    /**
     *  A crossing at the line's end.
     */
    pairs.push_back({Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 1, .Y = 0}},
                     StrandQuartic({{{.X = 0, .Y = -1},
                                     {.X =  0.1, .Y = -0.5},
                                     {.X = -0.1, .Y =    0},
                                     {.X =  0.1, .Y =  0.5},
                                     {.X = 0, .Y = 1}}})});
    /**
     *  A crossing whose y lies halfway between two doubles: settled by no bound.
     */
    pairs.push_back({Line{.Start = {.X = -1, .Y = 1}, .End = {.X = 1, .Y = 1 + 0x1p-52}},
                     StrandQuartic({{{.X = 0, .Y = -1},
                                     {.X = 0, .Y = 0},
                                     {.X = 0, .Y = 1},
                                     {.X = 0, .Y = 2},
                                     {.X = 0, .Y = 3}}})});
    /**
     *  Knots far from nothing, where a step of the parameter moves the curve past the tolerance.
     */
    Spline distant = StrandQuartic({{{.X = 0, .Y = -1},
                                     {.X = 1, .Y =  1},
                                     {.X = 2, .Y = -1},
                                     {.X = 3, .Y =  1},
                                     {.X = 4, .Y = -1}}});
    for (double &knot : distant.Knots)
    {
      knot += 0x1p40;
    }
    pairs.push_back({Line{.Start = {.X = 0.5, .Y = -5}, .End = {.X = 0.6, .Y = 5}}, distant});
    for (std::size_t index = 0; index < pairs.size(); ++index)
    {
      CAPTURE(index);
      const Curves  curves = {pairs[index][0], pairs[index][1]};
      Intersections meeting;
      CHECK(!Swift::Met(curves, tolerance, meeting));
      CHECK(meeting.Hits.empty());
      const Spline                    &spline = std::get<Spline>(pairs[index][1]);
      const std::optional<Swift::Lane> lane   = Swift::Laned(pairs[index][0]);
      REQUIRE(lane);
      CHECK(Swift::Skeined(*lane, spline, tolerance) == std::nullopt);
      const Intersections met = Intersect(pairs[index][0], pairs[index][1], tolerance);
      CHECK(met.Hits.size() == Intersect(pairs[index][1], pairs[index][0], tolerance).Hits.size());
    }
    /**
     *  A crossing too shallow for the tolerance's angle, the curve clear of the line elsewhere.
     */
    const Tolerance wide = {.Length = 1e-9, .Radians = 1e-6};
    const Curve     line = Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 1, .Y = 0}};
    const Curve     flat = StrandQuartic({{{.X = -1, .Y = -1e-6},
                                           {.X = -0.5, .Y = -0.5e-6},
                                           {.X =    0, .Y =       0},
                                           {.X =  0.5, .Y =  0.5e-6},
                                           {.X = 1, .Y = 1e-6}}});
    Intersections   shallow;
    CHECK(!Swift::Met({line, flat}, wide, shallow));
    /**
     *  A ray whose carrier crosses the spline behind its start: met, and nothing found.
     */
    const Curve     ray    = Ray{.Start = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 0}};
    const Curve     behind = StrandQuartic({{{.X = -0.5, .Y = -1},
                                             {.X = -0.4, .Y = -0.5},
                                             {.X = -0.6, .Y =    0},
                                             {.X = -0.4, .Y =  0.5},
                                             {.X = -0.5, .Y = 1}}});
    const StrandMet met    = StrandBoth({ray, behind}, tolerance);
    CHECK(met.Sure);
    CHECK(met.Agrees);
    CHECK(met.Swiftly.Hits.empty());
    /**
     *  A span that closes on itself, its chord of no length: crossed twice, and met.
     */
    const Curve     loop   = StrandQuartic({{{.X = 0, .Y = 0},
                                             {.X =  2, .Y = 1},
                                             {.X =  0, .Y = 2},
                                             {.X = -2, .Y = 1},
                                             {.X = 0, .Y = 0}}});
    const Curve     across = Line{.Start = {.X = -3, .Y = 0.9}, .End = {.X = 3, .Y = 1.1}};
    const StrandMet closed = StrandBoth({across, loop}, tolerance);
    CHECK(closed.Sure);
    CHECK(closed.Agrees);
    CHECK(closed.Swiftly.Hits.size() == 2);
  }
} // namespace Geoxy
