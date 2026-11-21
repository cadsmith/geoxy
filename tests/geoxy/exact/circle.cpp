/***************************************************************************************************
 *  @file      tests/geoxy/exact/circle.cpp
 *  @brief     Tests of the circle test's sign, against Shewchuk's own.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/circle.hpp"

#include "../draws.hpp"
#include "exact/bounds.hpp"
#include "exact/triangle.hpp"
#include "exact/two.hpp"
#include "original.hpp"
#include "same.hpp"

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>

namespace Geoxy::Exact
{
  namespace
  {
    auto CircleOriginal(const Triangle &circle, Point point) -> double
    {
      std::array<double, 2> a = {circle.First.X, circle.First.Y};
      std::array<double, 2> b = {circle.Second.X, circle.Second.Y};
      std::array<double, 2> c = {circle.Third.X, circle.Third.Y};
      std::array<double, 2> d = {point.X, point.Y};
      return incircle(a.data(), b.data(), c.data(), d.data());
    }

    /**
     *  Every point of a triangle and the point asked about, moved alike.
     */
    template <typename MoveT>
    auto Moved(const Triangle &circle, MoveT move) -> Triangle
    {
      return {
          .First = move(circle.First), .Second = move(circle.Second), .Third = move(circle.Third)};
    }

    /**
     *  InCircle gives the original's sign, its adaptive stages alone too, even where the first
     *  decides, and keeps the laws of the determinant: the same when the corners turn, the other
     *  when two swap or the plane is mirrored, and the same when the plane is scaled by a power of
     *  two.
     */
    auto CircleAgrees(const Triangle &circle, Point point) -> bool
    {
      const int sign = Sign(CircleOriginal(circle, point));

      const auto mirror = [](Point p) -> Point
      {
        return Point{.X = -p.X, .Y = p.Y};
      };
      const auto scale = [](Point p) -> Point
      {
        return Point{.X = std::ldexp(p.X, 32), .Y = std::ldexp(p.Y, 32)};
      };

      const Triangle turnedCircle = {
          .First = circle.Second, .Second = circle.Third, .Third = circle.First};
      const Triangle swappedCircle = {
          .First = circle.Second, .Second = circle.First, .Third = circle.Third};

      const bool same     = InCircle(circle, point) == sign;
      const bool adaptive = CircleAdapt(circle, point, Lift(circle, point).Permanent) == sign;
      const bool turned   = InCircle(turnedCircle, point) == sign;
      const bool swapped  = InCircle(swappedCircle, point) == -sign;
      const bool mirrored = InCircle(Moved(circle, mirror), mirror(point)) == -sign;
      const bool scaled   = InCircle(Moved(circle, scale), scale(point)) == sign;
      return same && adaptive && turned && swapped && mirrored && scaled;
    }

    /**
     *  The double the given number of steps away.
     */
    auto CircleNudged(double value, int steps) -> double
    {
      for (; steps > 0; --steps)
      {
        value = std::nextafter(value, HUGE_VAL);
      }
      for (; steps < 0; ++steps)
      {
        value = std::nextafter(value, -HUGE_VAL);
      }
      return value;
    }

    struct CircleFixture
    {
      CircleFixture() { exactinit(); }
    };
  } // namespace

  TEST_CASE_FIXTURE(CircleFixture,
                    "circle: its first and third stages' constants are the original's")
  {
    CHECK(Bits(CircleBoundA) == Bits(iccerrboundA));
    CHECK(Bits(CircleBoundC) == Bits(iccerrboundC));
  }

  /**
   *  The compensated second stage's bounds are the proof's (exact/circle.cpp): the original's
   *  second bound plus at least the estimate's own error, 47u^2, and the exact case's above it.
   */
  TEST_CASE_FIXTURE(CircleFixture,
                    "circle: its second stage's bounds cover the compensated estimate")
  {
    constexpr double SquaredUnit = Epsilon * Epsilon;
    CHECK(CircleBoundB >= iccerrboundB + (47 * SquaredUnit));
    CHECK(CircleBoundCompensated >= 47 * SquaredUnit);
  }

  TEST_CASE_FIXTURE(CircleFixture,
                    "circle: inside is positive, outside negative, on the circle zero")
  {
    const Point    south            = {.X = 0, .Y = -1};
    const Point    east             = {.X = 1, .Y = 0};
    const Point    north            = {.X = 0, .Y = 1};
    const Point    centre           = {.X = 0, .Y = 0};
    const Triangle counterclockwise = {.First = south, .Second = east, .Third = north};
    const Triangle clockwise        = {.First = north, .Second = east, .Third = south};
    CHECK(InCircle(counterclockwise, centre)            == 1);
    CHECK(InCircle(counterclockwise, {.X = 2, .Y = 0})  == -1);
    CHECK(InCircle(counterclockwise, {.X = -1, .Y = 0}) == 0);
    CHECK(InCircle(clockwise, centre) == -1);
    CHECK(CircleAgrees(counterclockwise, centre));
  }

  /**
   *  robust-predicates' cases: the circle through (0, x), (-x, -x) and (x, -x), for x from 10^-64
   *  up by tens to 10^63, with its centre inside, (0, 2x) outside and (0, x) on it.
   */
  TEST_CASE_FIXTURE(CircleFixture, "circle: a circle at every scale from 10^-64 to 10^63")
  {
    double x       = 1e-64;
    int    correct = 0;
    for (int index = 0; index < 128; ++index)
    {
      const Triangle circle = {
          .First = {.X = 0, .Y = x}, .Second = {.X = -x, .Y = -x}, .Third = {.X = x, .Y = -x}};
      const Point outside  = {.X = 0, .Y = 2 * x};
      correct             += static_cast<int>(InCircle(circle, {.X = 0, .Y = 0}) == 1
                                              && InCircle(circle, outside) == -1
                                              && InCircle(circle, {.X = 0, .Y = x}) == 0
                                              && CircleAgrees(circle, outside));
      x                   *= 10;
    }
    CHECK(correct == 128);
  }

  TEST_CASE_FIXTURE(CircleFixture, "circle: random places agree with the original at every scale")
  {
    Draws draws{0x6369'7263U};

    int agreed = 0;
    for (const double spread : {1e-20, 1.0, 1e20})
    {
      for (int index = 0; index < 20'000; ++index)
      {
        const auto next = [&draws, spread] -> Point
        {
          return Point{.X = draws.Near(0, spread), .Y = draws.Near(0, spread)};
        };
        const Triangle circle  = {.First = next(), .Second = next(), .Third = next()};
        agreed                += static_cast<int>(CircleAgrees(circle, next()));
      }
    }
    CHECK(agreed == 60'000);
  }

  /**
   *  Places on the unit circle, rounded, and a point on it moved a few doubles each way: inside,
   *  outside or on it as the rounding falls, and decided only by the later stages.
   */
  TEST_CASE_FIXTURE(CircleFixture,
                    "circle: places steps off the unit circle agree with the original")
  {
    Draws draws{0x756E'6974U};

    const auto onCircle = [&draws] -> Point
    {
      const double angle = draws.Near(0, std::numbers::pi);
      return Point{.X = std::cos(angle), .Y = std::sin(angle)};
    };

    int agreed = 0;
    int count  = 0;
    for (int index = 0; index < 1'000; ++index)
    {
      const Triangle circle = {.First = onCircle(), .Second = onCircle(), .Third = onCircle()};
      const Point    point  = onCircle();
      for (int dx = -2; dx <= 2; ++dx)
      {
        for (int dy = -2; dy <= 2; ++dy)
        {
          agreed += static_cast<int>(CircleAgrees(
              circle, {.X = CircleNudged(point.X, dx), .Y = CircleNudged(point.Y, dy)}));
          ++count;
        }
      }
    }
    CHECK(agreed == count);
  }

  /**
   *  Integer places on circles of integer radius, from Pythagorean triples about an integer centre:
   *  every difference exact, the determinant zero, and one step off it, not.
   */
  TEST_CASE_FIXTURE(CircleFixture, "circle: integer places on a circle are on it")
  {
    constexpr int Count = 20'000;

    Draws draws{0x7079'7468U};

    int onCircle = 0;
    int agreed   = 0;
    for (int index = 0; index < Count; ++index)
    {
      const auto larger = static_cast<double>(2 + (draws.Next() % 2000U));
      const auto smaller =
          static_cast<double>(1 + (draws.Next() % static_cast<unsigned>(larger - 1)));
      const double a      = (larger * larger) - (smaller * smaller);
      const double b      = 2 * larger * smaller;
      const Point  centre = {.X = static_cast<double>(draws.Next() >> 24U),
                             .Y = -static_cast<double>(draws.Next() >> 24U)};

      /**
       *  One of the eight places (+-a, +-b) and (+-b, +-a) at distance a^2 + b^2 from the centre.
       */
      const auto around = [&draws, a, b] -> Vector
      {
        const std::uint64_t pick    = draws.Next();
        const bool          swapped = (pick & 4U) != 0;
        const double        x       = swapped ? b : a;
        const double        y       = swapped ? a : b;
        return {.X = (pick & 1U) != 0 ? -x : x, .Y = (pick & 2U) != 0 ? -y : y};
      };
      std::array<Point, 4> picked = {};
      for (Point &place : picked)
      {
        place = centre + around();
      }
      const Triangle circle  = {.First = picked[0], .Second = picked[1], .Third = picked[2]};
      onCircle              += static_cast<int>(InCircle(circle, picked[3]) == 0);
      agreed                += static_cast<int>(CircleAgrees(circle, picked[3]));
      agreed += static_cast<int>(CircleAgrees(circle, {.X = picked[3].X + 1, .Y = picked[3].Y}));
    }
    CHECK(onCircle == Count);
    CHECK(agreed   == 2 * Count);
  }

  /**
   *  A circle of radius 5 2^50 about the origin through integer places, and points on it and a unit
   *  off it. Every difference is exact, and the determinant, at most 2^-100 of the permanent, is
   *  below what the compensated estimate can tell from zero: the exact sum decides.
   */
  TEST_CASE_FIXTURE(CircleFixture,
                    "circle: a unit off a vast circle is past the compensated estimate")
  {
    constexpr double Scale = 0x1.0p50;

    const Triangle circle   = {.First  = {.X = 3 * Scale, .Y = 4 * Scale},
                               .Second = {.X = -4 * Scale, .Y = 3 * Scale},
                               .Third  = {.X = -3 * Scale, .Y = -4 * Scale}};
    const Point    onCircle = {.X = 5 * Scale, .Y = 0};
    const Point    inside   = {.X = (5 * Scale) - 1, .Y = 0};
    const Point    beyond   = {.X = 5 * Scale, .Y = 1};

    const double permanent = Lift(circle, beyond).Permanent;
    CHECK(std::abs(CircleOriginal(circle, beyond)) < CircleBoundCompensated * permanent);
    CHECK(InCircle(circle, onCircle) == 0);
    CHECK(InCircle(circle, inside)   == 1);
    CHECK(InCircle(circle, beyond)   == -1);
    CHECK(CircleAgrees(circle, onCircle));
    CHECK(CircleAgrees(circle, inside));
    CHECK(CircleAgrees(circle, beyond));
  }

  /**
   *  The point on a grid of neighbouring doubles around (2, 2), on the circle through (0, 0), (2,
   *  0) and (0, 2), all moved by an offset. An offset that rounds in x, in y or in both makes those
   *  differences inexact, and so reaches each of the last stage's branches.
   */
  TEST_CASE_FIXTURE(CircleFixture,
                    "circle: grids around a place on a circle agree with the original")
  {
    constexpr int Side = 48;

    int agreed = 0;
    int count  = 0;
    for (const Point offset : {Point{.X = 0, .Y = 0}, Point{.X = 0.3, .Y = 0},
                               Point{.X = 0, .Y = 0.3}, Point{.X = 0.3, .Y = 0.7}})
    {
      const Triangle circle = {.First  = {.X = offset.X, .Y = offset.Y},
                               .Second = {.X = 2 + offset.X, .Y = offset.Y},
                               .Third  = {.X = offset.X, .Y = 2 + offset.Y}};
      const Point    corner = {.X = 2 + offset.X, .Y = 2 + offset.Y};
      for (int column = -Side; column < Side; ++column)
      {
        for (int row = -Side; row < Side; ++row)
        {
          agreed += static_cast<int>(CircleAgrees(
              circle, {.X = CircleNudged(corner.X, column), .Y = CircleNudged(corner.Y, row)}));
          ++count;
        }
      }
    }
    CHECK(agreed == count);
  }

  /**
   *  Circles and points over a few values, signed zeros among them: repeated places, corners on a
   *  line, products of zero.
   */
  TEST_CASE_FIXTURE(CircleFixture, "circle: places over a few values agree with the original")
  {
    static constexpr std::array Values{0.0, -0.0, 1.0, -1.0, 0.1, 0x1.0p26};
    constexpr int               Count = 50'000;

    Draws draws{0x7661'6CU};

    const auto pick = [&draws] -> double
    {
      return Values.at(draws.Next() % Values.size());
    };

    int agreed = 0;
    for (int index = 0; index < Count; ++index)
    {
      const Triangle circle  = {.First  = {.X = pick(), .Y = pick()},
                                .Second = {.X = pick(), .Y = pick()},
                                .Third  = {.X = pick(), .Y = pick()}};
      agreed                += static_cast<int>(CircleAgrees(circle, {.X = pick(), .Y = pick()}));
    }
    CHECK(agreed == Count);
  }
} // namespace Geoxy::Exact
