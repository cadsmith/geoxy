/***************************************************************************************************
 *  @file      tests/geoxy/exact/whole.cpp
 *  @brief     Tests of both predicates past Shewchuk's limits, of their whole number sums, and of
 *             where two lines cross in whole numbers.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/whole.hpp"

#include "../draws.hpp"
#include "exact/batch.hpp"
#include "exact/circle.hpp"
#include "exact/orient.hpp"
#include "exact/triangle.hpp"
#include "exact/two.hpp"
#include "expansion.hpp"
#include "original.hpp"

#include "geoxy/point.hpp"

#include <doctest/doctest.h>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>

namespace Geoxy::Exact
{
  namespace
  {
    auto OriginalOrient(const Triangle &corners) -> int
    {
      std::array<double, 2> a = {corners.First.X, corners.First.Y};
      std::array<double, 2> b = {corners.Second.X, corners.Second.Y};
      std::array<double, 2> c = {corners.Third.X, corners.Third.Y};
      return Sign(orient2d(a.data(), b.data(), c.data()));
    }

    auto OriginalCircle(const Triangle &circle, Point point) -> int
    {
      std::array<double, 2> a = {circle.First.X, circle.First.Y};
      std::array<double, 2> b = {circle.Second.X, circle.Second.Y};
      std::array<double, 2> c = {circle.Third.X, circle.Third.Y};
      std::array<double, 2> d = {point.X, point.Y};
      return Sign(incircle(a.data(), b.data(), c.data(), d.data()));
    }

    /**
     *  A point scaled by a power of two, or none where a coordinate would round, underflow or
     *  overflow.
     */
    auto Scaled(Point point, int exponent, bool &exact) -> Point
    {
      const auto scale = [&exact, exponent](double value) -> double
      {
        const double scaled = std::ldexp(value, exponent);
        exact               = exact
                           && std::bit_cast<std::uint64_t>(std::ldexp(scaled, -exponent))
                                  == std::bit_cast<std::uint64_t>(value);
        return scaled;
      };
      const double x = scale(point.X);
      return {.X = x, .Y = scale(point.Y)};
    }

    /**
     *  Three places near a line, the third a few steps off it.
     */
    auto NearLine(Draws &draws) -> Triangle
    {
      const Point  first  = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)};
      const Point  second = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)};
      const double along  = draws.Near(0.5, 1);
      return {.First  = first,
              .Second = second,
              .Third  = {.X = first.X + (along * (second.X - first.X)),
                         .Y = first.Y + (along * (second.Y - first.Y))}};
    }

    /**
     *  Four places on the unit circle, rounded.
     */
    auto NearCircle(Draws &draws, Point &asked) -> Triangle
    {
      const auto onCircle = [&draws] -> Point
      {
        const double angle = draws.Near(0, std::numbers::pi);
        return {.X = std::cos(angle), .Y = std::sin(angle)};
      };
      asked = onCircle();
      return {.First = onCircle(), .Second = onCircle(), .Third = onCircle()};
    }

    /**
     *  Where two straight carriers cross, worked out in expansions and rounded once, as the exact
     *  crossing works it out where sums of doubles hold its products: the oracle the whole numbers
     *  are held to there.
     *
     *  @param  ways The two carriers.
     *  @return      The double nearest the crossing.
     */
    auto Summed(const std::array<Straightway, 2> &ways) -> Point
    {
      using Expansion::Difference;
      using Expansion::Number;
      using Expansion::Of;
      using Expansion::Rounded;
      using Expansion::Sum;
      const auto way = [](const Straightway &line) -> std::array<Number, 2>
      {
        if (line.Way)
        {
          return {Of(line.Toward.X), Of(line.Toward.Y)};
        }
        return {Difference(Of(line.Toward.X), Of(line.Through.X)),
                Difference(Of(line.Toward.Y), Of(line.Through.Y))};
      };
      const std::array<Number, 2> one   = way(ways[0]);
      const std::array<Number, 2> other = way(ways[1]);
      const std::array<Number, 2> gap = {Difference(Of(ways[1].Through.X), Of(ways[0].Through.X)),
                                         Difference(Of(ways[1].Through.Y), Of(ways[0].Through.Y))};
      const Number                turn =
          Difference(Expansion::Product(one[0], other[1]), Expansion::Product(one[1], other[0]));
      const Number reach =
          Difference(Expansion::Product(gap[0], other[1]), Expansion::Product(gap[1], other[0]));
      const auto axis = [&turn, &reach](double start, const Number &step) -> double
      {
        return Rounded(
            {.Top    = Sum(Expansion::Product(Of(start), turn), Expansion::Product(step, reach)),
             .Bottom = {},
             .Root   = {},
             .Below  = turn},
            1);
      };
      return {.X = axis(ways[0].Through.X, one[0]), .Y = axis(ways[0].Through.Y, one[1])};
    }

    /**
     *  A carrier scaled by a power of two, and whether every number of it scaled exactly.
     *
     *  @param  line     The carrier.
     *  @param  exponent The power.
     *  @param  exact    Cleared where a number rounds, underflows or overflows.
     *  @return          The carrier scaled.
     */
    auto Scaled(const Straightway &line, int exponent, bool &exact) -> Straightway
    {
      return {.Through = Scaled(line.Through, exponent, exact),
              .Toward  = Scaled(line.Toward, exponent, exact),
              .Way     = line.Way};
    }

    /**
     *  Whether two points are the same to the bit, the signs of zeros too.
     *
     *  @param  lhs One point.
     *  @param  rhs The other.
     *  @return     Whether they are.
     */
    auto Bitwise(Point lhs, Point rhs) -> bool
    {
      return std::bit_cast<std::uint64_t>(lhs.X) == std::bit_cast<std::uint64_t>(rhs.X)
          && std::bit_cast<std::uint64_t>(lhs.Y) == std::bit_cast<std::uint64_t>(rhs.Y);
    }

    struct WholeFixture
    {
      WholeFixture() { exactinit(); }
    };
  } // namespace

  TEST_CASE_FIXTURE(WholeFixture, "whole: the sums agree with the original inside its limits")
  {
    Draws draws{0x7768'6F6CU};

    int agreed = 0;
    for (int index = 0; index < 20'000; ++index)
    {
      const Triangle spread = {.First  = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)},
                               .Second = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)},
                               .Third  = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)}};
      const Triangle line   = NearLine(draws);
      Point          asked;
      const Triangle circle = NearCircle(draws, asked);
      agreed += static_cast<int>(OrientWhole(spread) == OriginalOrient(spread)
                                 && OrientWhole(line) == OriginalOrient(line)
                                 && CircleWhole(circle, asked) == OriginalCircle(circle, asked)
                                 && CircleWhole(spread, line.Third)
                                        == OriginalCircle(spread, line.Third));
    }
    CHECK(agreed == 20'000);
  }

  /**
   *  Scaling every coordinate by 2^k scales the orientation determinant by 2^2k and the circle
   *  determinant by 2^4k: neither sign changes. So places near a line or a circle, scaled from
   *  2^-1070 to 2^1000, far past Shewchuk's limits, keep the sign the original gives unscaled,
   *  wherever the scaling is exact.
   */
  TEST_CASE_FIXTURE(WholeFixture, "whole: places scaled past the limits keep the original's sign")
  {
    Draws draws{0x7363'616CU};

    int checked = 0;
    int agreed  = 0;
    for (int index = 0; index < 4'000; ++index)
    {
      const Triangle line = NearLine(draws);
      Point          asked;
      const Triangle circle      = NearCircle(draws, asked);
      const int      orientation = OriginalOrient(line);
      const int      inside      = OriginalCircle(circle, asked);
      for (const int exponent : {-1070, -1030, -900, -500, -160, 160, 500, 900, 1000})
      {
        bool           exact       = true;
        const Triangle scaledLine  = {.First  = Scaled(line.First, exponent, exact),
                                      .Second = Scaled(line.Second, exponent, exact),
                                      .Third  = Scaled(line.Third, exponent, exact)};
        const Triangle scaledShape = {.First  = Scaled(circle.First, exponent, exact),
                                      .Second = Scaled(circle.Second, exponent, exact),
                                      .Third  = Scaled(circle.Third, exponent, exact)};
        const Point    scaledAsked = Scaled(asked, exponent, exact);
        if (!exact)
        {
          continue;
        }
        ++checked;
        agreed += static_cast<int>(
            Exact::Orient(scaledLine.First, scaledLine.Second, scaledLine.Third) == orientation
            && OrientWhole(scaledLine)                                           == orientation
            && InCircle(scaledShape, scaledAsked)                                == inside
            && CircleWhole(scaledShape, scaledAsked)                             == inside);
      }
    }
    CHECK(checked > 20'000);
    CHECK(agreed == checked);
  }

  /**
   *  Ozaki et al.'s case, where the product (ay - cy)(bx - cx) underflows: a first stage without
   *  u_N decides it zero.
   */
  TEST_CASE_FIXTURE(WholeFixture, "whole: an underflowing product is not taken for zero")
  {
    const Triangle corners = {.First  = {.X = 0x1.0p-702, .Y = 0x1.0p-701},
                              .Second = {.X = 0x1.0p-700, .Y = 0x1.0p-700},
                              .Third  = {.X = 0x1.0p-699, .Y = 0x1.0p-700}};
    bool           exact   = true;
    const Triangle scaled  = {.First  = Scaled(corners.First, 1000, exact),
                              .Second = Scaled(corners.Second, 1000, exact),
                              .Third  = Scaled(corners.Third, 1000, exact)};
    REQUIRE(exact);
    REQUIRE(OriginalOrient(scaled) != 0);
    CHECK(Exact::Orient(corners.First, corners.Second, corners.Third) == OriginalOrient(scaled));
    CHECK(OrientWhole(corners) == OriginalOrient(scaled));
  }

  /**
   *  Places at either end of the doubles: differences and products that overflow, a product of
   *  2^1000 and 2^-1000, the smallest doubles of all.
   */
  TEST_CASE_FIXTURE(WholeFixture, "whole: places at the ends of the doubles")
  {
    constexpr double Large = std::numeric_limits<double>::max();
    constexpr double Tiny  = std::numeric_limits<double>::denorm_min();

    CHECK(Exact::Orient({.X = 0x1.0p1023, .Y = 0}, {.X = 0, .Y = 0x1.0p1023},
                        {.X = -0x1.0p1023, .Y = 0})
          == 1);
    CHECK(Exact::Orient({.X = 0x1.0p1000, .Y = 0}, {.X = 0, .Y = 0x1.0p-1000}, {.X = 0, .Y = 0})
          == 1);
    CHECK(
        Exact::Orient({.X = Large, .Y = Large}, {.X = -Large, .Y = -Large}, {.X = Tiny, .Y = Tiny})
        == 0);
    CHECK(Exact::Orient({.X = Large, .Y = Large}, {.X = -Large, .Y = -Large}, {.X = Tiny, .Y = 0})
          == 1);
    CHECK(Exact::Orient({.X = 0, .Y = 0}, {.X = Tiny, .Y = 0}, {.X = 0, .Y = Tiny}) == 1);

    const Triangle unit = {
        .First = {.X = Tiny, .Y = 0}, .Second = {.X = 0, .Y = Tiny}, .Third = {.X = -Tiny, .Y = 0}};
    CHECK(InCircle(unit, {.X = 0, .Y = -Tiny})    == 0);
    CHECK(InCircle(unit, {.X = 0, .Y = 0})        == 1);
    CHECK(InCircle(unit, {.X = 0, .Y = 2 * Tiny}) == -1);

    const Triangle huge = {.First  = {.X = 0x1.0p1022, .Y = 0},
                           .Second = {.X = 0, .Y = 0x1.0p1022},
                           .Third  = {.X = -0x1.0p1022, .Y = 0}};
    CHECK(InCircle(huge, {.X = 0, .Y = -0x1.0p1022}) == 0);
    CHECK(InCircle(huge, {.X = 0x1.0p-1074, .Y = 0}) == 1);
    CHECK(InCircle(huge, {.X = 0, .Y = Large})       == -1);
  }

  /**
   *  The fuzzer's find: near the top of the doubles, a difference's roundoff overflows on its way
   *  and comes out NaN, which must not pass for an exact difference.
   */
  TEST_CASE_FIXTURE(WholeFixture, "whole: a roundoff that overflows is not taken for exact")
  {
    const Point a = {.X = 0x1.fffffffffffffp+1023, .Y = -0x1.e9420010624ddp+18};
    const Point b = {.X = 0x1.ffffffffffff7p+1022, .Y = -0x1.e94205cac0831p+18};
    const Point c = {.X = 0x1.7fffffffffffdp+1023, .Y = -0x1.e94202ed91687p+18};
    REQUIRE(OrientWhole({.First = a, .Second = b, .Third = c}) == -1);
    CHECK(Exact::Orient(a, b, c) == -1);
    CHECK(Exact::Orient(b, c, a) == -1);
    CHECK(Exact::Orient(c, a, b) == -1);
    CHECK(Exact::Orient(b, a, c) == 1);
    CHECK(OrientWhole({.First = b, .Second = c, .Third = a}) == -1);
  }

  TEST_CASE_FIXTURE(WholeFixture, "whole: a coordinate NaN or infinite gives zero, everywhere")
  {
    constexpr double Infinity = std::numeric_limits<double>::infinity();
    constexpr double Nothing  = std::numeric_limits<double>::quiet_NaN();

    const Triangle circle = {
        .First = {.X = 0, .Y = 1}, .Second = {.X = -1, .Y = -1}, .Third = {.X = 1, .Y = -1}};
    for (const double bad : {Infinity, -Infinity, Nothing})
    {
      CHECK(Exact::Orient({.X = bad, .Y = 0}, {.X = 1, .Y = 0}, {.X = 0, .Y = 1})   == 0);
      CHECK(Exact::Orient({.X = 0, .Y = 0}, {.X = 1, .Y = bad}, {.X = 0, .Y = 1})   == 0);
      CHECK(Exact::Orient({.X = 0, .Y = 0}, {.X = 1, .Y = 0}, {.X = bad, .Y = bad}) == 0);
      CHECK(InCircle(circle, {.X = bad, .Y = 0}) == 0);
      CHECK(InCircle({.First = circle.First, .Second = {.X = 0, .Y = bad}, .Third = circle.Third},
                     {.X = 0, .Y = 0})
            == 0);

      const std::array<Point, 3> points = {Point{.X = bad, .Y = 0}, Point{.X = 2, .Y = 2},
                                           Point{.X = 0, .Y = bad}};
      std::array<std::int8_t, 3> signs  = {7, 7, 7};
      OrientSigns({.First = {.X = 0, .Y = 0}, .Second = {.X = 1, .Y = 0}}, points, signs);
      CHECK(signs[0] == 0);
      CHECK(signs[1] == 1);
      CHECK(signs[2] == 0);
    }
  }

  /**
   *  The integer lattice of (3, 4) about a centre, scaled down to subnormals or up to 2^1000: on
   *  its circle, and a step off it, past the grid and the limits both.
   */
  TEST_CASE_FIXTURE(WholeFixture, "whole: a lattice circle place either end of the doubles")
  {
    for (const int exponent : {-1072, -1000, 900, 1015})
    {
      const auto place = [exponent](double x, double y) -> Point
      {
        return {.X = std::ldexp(x, exponent), .Y = std::ldexp(y, exponent)};
      };
      const Triangle circle = {
          .First = place(3, 4), .Second = place(-4, 3), .Third = place(-3, -4)};
      CHECK(InCircle(circle, place(5, 0))    == 0);
      CHECK(InCircle(circle, place(4, 0))    == 1);
      CHECK(InCircle(circle, place(5, 0.5))  == -1);
      CHECK(CircleWhole(circle, place(5, 0)) == 0);
    }
  }

  /**
   *  Lines and rays drawn at random cross, in whole numbers, at the double sums of doubles round
   *  their crossing to, where those hold its products: and, every number scaled by a power of two
   *  from 2^-1070 to 2^1000, far past where they hold them, at that double scaled alike, wherever
   *  the numbers and the crossing scale exactly.
   */
  TEST_CASE("whole: lines cross where sums of doubles put them, at every scale")
  {
    Draws draws{0x6C69'6E65U};
    int   agreed = 0;
    int   scaled = 0;
    int   held   = 0;
    for (int round = 0; round < Rounds(4000); ++round)
    {
      const auto line = [&draws](bool way) -> Straightway
      {
        return {.Through = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)},
                .Toward  = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)},
                .Way     = way};
      };
      const std::array<Straightway, 2> ways   = {line((draws.Next() & 1U) != 0), line(false)};
      const Point                      summed = Summed(ways);
      agreed += static_cast<int>(Bitwise(CrossWhole(ways).value_or(Point{}), summed));
      for (const int exponent : {-1070, -1000, -600, 600, 1000})
      {
        bool                             exact = true;
        const std::array<Straightway, 2> moved = {Scaled(ways[0], exponent, exact),
                                                  Scaled(ways[1], exponent, exact)};
        const Point                      there = Scaled(summed, exponent, exact);
        if (exact)
        {
          ++scaled;
          held += static_cast<int>(Bitwise(CrossWhole(moved).value_or(Point{}), there));
        }
      }
    }
    CHECK(agreed == Rounds(4000));
    CHECK(held   == scaled);
    CHECK(scaled > Rounds(4000) * 2);
  }

  /**
   *  Crossings whose doubles the whole numbers alone tell: a halfway point between the least
   *  doubles, to the even one; an end three steps of the least double off the axis, shared;
   *  crossings at the largest double and past it, an infinity; a crossing at zero, of the ways'
   *  cross product's sign; and none for parallel carriers, numbers that are not finite, or carriers
   *  of no numbers at all.
   */
  TEST_CASE("whole: lines crossing between the least doubles, past the largest, and at zero")
  {
    constexpr double Infinity = std::numeric_limits<double>::infinity();
    constexpr double Nothing  = std::numeric_limits<double>::quiet_NaN();

    const Straightway upright = {.Through = {.X = 1, .Y = -1}, .Toward = {.X = 1, .Y = 1}};
    const auto        down    = [](double height) -> Straightway
    {
      return {.Through = {.X = 0, .Y = height}, .Toward = {.X = 2, .Y = 0}};
    };
    const Point half  = CrossWhole({down(0x1p-1074), upright}).value_or(Point{.X = 0, .Y = 1});
    const Point three = CrossWhole({down(0x3p-1074), upright}).value_or(Point{});
    CHECK(Bitwise(half, Point{.X = 1, .Y = 0}));
    CHECK(Bitwise(three, Point{.X = 1, .Y = 0x1p-1073}));

    const Point shared = {.X = 1, .Y = 0x3p-1074};
    const auto  joint  = CrossWhole({Straightway{.Through = {.X = 1, .Y = -1}, .Toward = shared},
                                     Straightway{.Through = shared, .Toward = {.X = 1.5, .Y = 1}}});
    CHECK(Bitwise(joint.value_or(Point{}), shared));

    const Straightway level = {.Through = {}, .Toward = {.X = 1, .Y = 0}, .Way = true};
    const auto        tilt  = [](double slope) -> Straightway
    {
      return {.Through = {.X = 0, .Y = 1}, .Toward = {.X = 1, .Y = -slope}, .Way = true};
    };
    const Point largest = CrossWhole({level, tilt(0x1p-1023)}).value_or(Point{});
    const Point past    = CrossWhole({level, tilt(0x1p-1074)}).value_or(Point{});
    CHECK(Bitwise(largest, Point{.X = 0x1p1023, .Y = -0.0}));
    CHECK(Bitwise(past, Point{.X = Infinity, .Y = -0.0}));

    const Straightway lattice = {.Through = {.X = 0, .Y = 0}, .Toward = {.X = 4, .Y = 2}};
    const Straightway across  = {.Through = {.X = 0, .Y = 2}, .Toward = {.X = 4, .Y = 0}};
    CHECK(Bitwise(CrossWhole({lattice, across}).value_or(Point{}), Point{.X = 2, .Y = 1}));

    const Straightway rising  = {.Through = {.X = -1, .Y = -1}, .Toward = {.X = 1, .Y = 1}};
    const Straightway falling = {.Through = {.X = -1, .Y = 1}, .Toward = {.X = 1, .Y = -1}};
    CHECK(Bitwise(CrossWhole({rising, falling}).value_or(Point{}), Point{.X = -0.0, .Y = -0.0}));
    CHECK(Bitwise(CrossWhole({falling, rising}).value_or(Point{}), Point{.X = 0.0, .Y = 0.0}));

    CHECK(!CrossWhole({rising, rising}).has_value());
    CHECK(!CrossWhole({rising, Straightway{.Through = {.X = Nothing, .Y = 0}}}).has_value());
    CHECK(!CrossWhole({rising, Straightway{.Through = {.X = 0, .Y = Infinity}}}).has_value());
    CHECK(!CrossWhole({Straightway{.Way = true}, Straightway{.Way = true}}).has_value());
  }
} // namespace Geoxy::Exact
